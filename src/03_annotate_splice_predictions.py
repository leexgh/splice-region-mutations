"""Step 3 — Splice prediction annotations.

Three signals merged into one TSV:

1. **SpliceAI** (local install). Build a minimal VCF from the unique variants,
   run `spliceai -I in.vcf -O out.vcf -R hg19.fa -A grch37`, parse the
   `SpliceAI=...` INFO field for `DS_AG/AL/DG/DL` (delta scores).
2. **SpliceVault** (REST API). Per variant, query for observed mis-splicing
   events at the same junction in GTEx RNA-seq. Cached per-variant.
3. **MaxEntScan** (local `maxentpy`). For variants in or adjacent to a splice
   site, score the donor (9nt) and acceptor (23nt) windows on ref and alt;
   report `delta_mes`.

Output `out/annotated/splice_predictions.tsv` with 5-column join key + scores.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _apiclient import (  # noqa: E402
    cached_read,
    cached_write,
    polite_get,
    variant_key,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
UNIQUE_VARIANTS = PROJECT_ROOT / "out" / "unique_variants.tsv"
OUTPUT_TSV = PROJECT_ROOT / "out" / "annotated" / "splice_predictions.tsv"
HG19_FASTA = PROJECT_ROOT / "data" / "refs" / "hg19.fa"
SPLICEAI_VCF_IN = PROJECT_ROOT / "out" / "cache" / "spliceai" / "input.vcf"
SPLICEAI_VCF_OUT = PROJECT_ROOT / "out" / "cache" / "spliceai" / "output.vcf"

log = logging.getLogger("step3")
logging.basicConfig(
    format="[%(name)s] %(message)s", level=logging.INFO, stream=sys.stderr
)


# --------------- SpliceAI ---------------


def write_spliceai_vcf(variants: pd.DataFrame, vcf_path: Path, fasta_path: Path) -> None:
    """Emit a minimal VCF with chr-prefixed contig names + contig length
    headers (pysam, used by SpliceAI to write output, requires them)."""
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    # Collect chromosomes we'll actually use + their lengths from the .fai.
    chrom_lengths: dict[str, int] = {}
    fai_path = Path(str(fasta_path) + ".fai")
    if fai_path.exists():
        for line in fai_path.read_text().splitlines():
            cols = line.split("\t")
            if len(cols) >= 2:
                chrom_lengths[cols[0]] = int(cols[1])

    used_chroms = sorted(
        {
            (c if c.startswith("chr") else f"chr{c}")
            for c in variants["chromosome"].astype(str).unique()
        },
        key=lambda x: (len(x), x),
    )

    with vcf_path.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in used_chroms:
            length = chrom_lengths.get(c, 250_000_000)
            f.write(f"##contig=<ID={c},length={length}>\n")
        f.write("##INFO=<ID=KEY,Number=1,Type=String,Description=\"variant_key\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, r in variants.iterrows():
            chrom = r["chromosome"]
            chrom_vcf = chrom if chrom.startswith("chr") else f"chr{chrom}"
            key = variant_key(
                chrom,
                int(r["start_position"]),
                int(r["end_position"]),
                r["reference_allele"],
                r["variant_allele"],
            )
            f.write(
                f"{chrom_vcf}\t{r['start_position']}\t.\t{r['reference_allele']}\t"
                f"{r['variant_allele']}\t.\t.\tKEY={key}\n"
            )


def run_spliceai(in_vcf: Path, out_vcf: Path, fasta: Path, spliceai_bin: str = "spliceai") -> None:
    """Invoke the spliceai CLI. Pre-built annotation file `grch37` ships with
    the package."""
    if out_vcf.exists() and out_vcf.stat().st_size > 1024:
        log.info("  spliceai output exists, skipping run (delete %s to force)", out_vcf)
        return
    cmd = [
        spliceai_bin,
        "-I", str(in_vcf),
        "-O", str(out_vcf),
        "-R", str(fasta),
        "-A", "grch37",
    ]
    log.info("  running: %s", " ".join(cmd))
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        log.error("  spliceai stderr (tail):\n%s", proc.stderr[-2000:])
        raise RuntimeError(f"spliceai exited {proc.returncode}")
    log.info("  spliceai completed in %.1fs", elapsed)


_INFO_RE = re.compile(r"SpliceAI=([^;\t\n]+)")


def parse_spliceai_vcf(out_vcf: Path) -> dict[str, dict]:
    """Map variant_key -> {ds_ag, ds_al, ds_dg, ds_dl, dp_ag, dp_al, dp_dg, dp_dl,
    spliceai_max, spliceai_symbol}. If multiple transcript annotations exist,
    take the maximum delta-score across them."""
    out: dict[str, dict] = {}
    if not out_vcf.exists():
        return out
    with out_vcf.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            info = cols[7]
            # KEY=... — preserve the variant_key written in input
            kmatch = re.search(r"KEY=([^;\t]+)", info)
            if not kmatch:
                continue
            key = kmatch.group(1)
            m = _INFO_RE.search(info)
            entry = out.setdefault(
                key,
                {
                    "spliceai_symbol": "",
                    "spliceai_ds_ag": float("nan"),
                    "spliceai_ds_al": float("nan"),
                    "spliceai_ds_dg": float("nan"),
                    "spliceai_ds_dl": float("nan"),
                    "spliceai_dp_ag": "",
                    "spliceai_dp_al": "",
                    "spliceai_dp_dg": "",
                    "spliceai_dp_dl": "",
                    "spliceai_max": float("nan"),
                    "spliceai_status": "no_annotation",
                },
            )
            if not m:
                continue
            # SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
            # comma-separated for each gene annotated
            for ann in m.group(1).split(","):
                parts = ann.split("|")
                if len(parts) < 10:
                    continue
                _allele, symbol, ds_ag, ds_al, ds_dg, ds_dl, dp_ag, dp_al, dp_dg, dp_dl = parts[:10]
                try:
                    dsv = [float(x) for x in (ds_ag, ds_al, ds_dg, ds_dl)]
                except ValueError:
                    continue
                cur_max = entry["spliceai_max"]
                if (entry["spliceai_status"] == "no_annotation") or (
                    pd.notna(cur_max) and max(dsv) > cur_max
                ) or pd.isna(cur_max):
                    entry["spliceai_symbol"] = symbol
                    entry["spliceai_ds_ag"] = dsv[0]
                    entry["spliceai_ds_al"] = dsv[1]
                    entry["spliceai_ds_dg"] = dsv[2]
                    entry["spliceai_ds_dl"] = dsv[3]
                    entry["spliceai_dp_ag"] = dp_ag
                    entry["spliceai_dp_al"] = dp_al
                    entry["spliceai_dp_dg"] = dp_dg
                    entry["spliceai_dp_dl"] = dp_dl
                    entry["spliceai_max"] = max(dsv)
                    entry["spliceai_status"] = "scored"
    return out


# --------------- SpliceVault ---------------

# SpliceVault REST API base. See https://kidsneuro.shinyapps.io/splicevault/.
# The community-known endpoint (subject to availability):
SPLICEVAULT_API = "https://splicevault-api.org/api/variant"


def query_splicevault(variant_key_str: str, hgvsc: str) -> dict:
    """Best-effort SpliceVault lookup by HGVSc (their query unit). Returns a
    flat dict; on failure returns a stub with `splicevault_status=error`.

    Note: SpliceVault's REST endpoint is research-grade and may be down. We
    cache `__not_found__` so a 404 doesn't keep retrying."""
    if not hgvsc:
        return {
            "splicevault_observed": False,
            "splicevault_top_event": "",
            "splicevault_n_events": 0,
            "splicevault_status": "no_hgvsc",
        }
    try:
        payload = polite_get(
            SPLICEVAULT_API,
            api="splicevault",
            key=variant_key_str,
            params={"query": hgvsc},
            delay=0.5,
            timeout=20.0,
        )
    except Exception as e:  # noqa: BLE001
        log.debug("  splicevault fetch error for %s: %s", hgvsc, e)
        return {
            "splicevault_observed": False,
            "splicevault_top_event": "",
            "splicevault_n_events": 0,
            "splicevault_status": "error",
        }

    if payload is None:
        return {
            "splicevault_observed": False,
            "splicevault_top_event": "",
            "splicevault_n_events": 0,
            "splicevault_status": "not_found",
        }

    events = (payload or {}).get("events") or []
    return {
        "splicevault_observed": bool(events),
        "splicevault_top_event": (events[0] if events else {}).get("event_type", ""),
        "splicevault_n_events": len(events),
        "splicevault_status": "ok",
    }


# --------------- MaxEntScan ---------------


def score_maxentscan(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    fasta,
    score3=None,
    score5=None,
) -> dict:
    """Score 9nt donor and 23nt acceptor windows around the variant on ref and
    alt. Report the larger absolute ΔMES across the two windows.

    Skips indels (maxentpy expects fixed-window SNVs).
    """
    out = {
        "maxent_donor_ref": float("nan"),
        "maxent_donor_alt": float("nan"),
        "maxent_acceptor_ref": float("nan"),
        "maxent_acceptor_alt": float("nan"),
        "maxent_delta": float("nan"),
        "maxent_status": "skipped_indel",
    }
    if len(ref) != 1 or len(alt) != 1:
        return out
    if score3 is None or score5 is None:
        out["maxent_status"] = "maxent_unavailable"
        return out

    chrom_key = chrom if chrom.startswith("chr") else f"chr{chrom}"
    try:
        # Donor: 9nt window centred near the splice donor (3 exon + 6 intron)
        # We approximate by scoring all 9nt windows containing pos, and take
        # the most-changed one.
        seq_around = fasta[chrom_key][pos - 25 : pos + 25].seq.upper()
        offset_in_chunk = 25  # pos-25 is index 0
        # Build alt sequence (substitution).
        if seq_around[offset_in_chunk - 1] != ref:
            out["maxent_status"] = "ref_mismatch"
            return out
        seq_alt = (
            seq_around[: offset_in_chunk - 1]
            + alt
            + seq_around[offset_in_chunk:]
        )

        best_donor = (float("nan"), float("nan"), 0.0)
        for start_off in range(-8, 1):
            if offset_in_chunk - 1 + start_off < 0 or offset_in_chunk + start_off + 9 > len(seq_around):
                continue
            ref_window = seq_around[
                offset_in_chunk - 1 + start_off : offset_in_chunk - 1 + start_off + 9
            ]
            alt_window = seq_alt[
                offset_in_chunk - 1 + start_off : offset_in_chunk - 1 + start_off + 9
            ]
            try:
                rs = score5(ref_window)
                als = score5(alt_window)
            except Exception:
                continue
            d = abs(als - rs)
            if d > best_donor[2]:
                best_donor = (rs, als, d)

        best_accept = (float("nan"), float("nan"), 0.0)
        for start_off in range(-22, 1):
            if offset_in_chunk - 1 + start_off < 0 or offset_in_chunk + start_off + 23 > len(seq_around):
                continue
            ref_window = seq_around[
                offset_in_chunk - 1 + start_off : offset_in_chunk - 1 + start_off + 23
            ]
            alt_window = seq_alt[
                offset_in_chunk - 1 + start_off : offset_in_chunk - 1 + start_off + 23
            ]
            try:
                rs = score3(ref_window)
                als = score3(alt_window)
            except Exception:
                continue
            d = abs(als - rs)
            if d > best_accept[2]:
                best_accept = (rs, als, d)

        out["maxent_donor_ref"] = best_donor[0]
        out["maxent_donor_alt"] = best_donor[1]
        out["maxent_acceptor_ref"] = best_accept[0]
        out["maxent_acceptor_alt"] = best_accept[1]
        out["maxent_delta"] = max(best_donor[2], best_accept[2])
        out["maxent_status"] = "ok"
    except Exception as e:  # noqa: BLE001
        log.debug("  maxent error: %s", e)
        out["maxent_status"] = f"error: {type(e).__name__}"
    return out


# --------------- main ---------------


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(UNIQUE_VARIANTS))
    parser.add_argument("--output", default=str(OUTPUT_TSV))
    parser.add_argument("--fasta", default=str(HG19_FASTA))
    parser.add_argument("--skip-spliceai", action="store_true")
    parser.add_argument("--skip-splicevault", action="store_true")
    parser.add_argument("--skip-maxent", action="store_true")
    parser.add_argument(
        "--spliceai-bin",
        default=str(PROJECT_ROOT / ".venv" / "bin" / "spliceai"),
        help="Path to spliceai CLI (default: project venv to avoid AVX issues).",
    )
    args = parser.parse_args()

    variants = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    log.info("loaded %d unique variants", len(variants))

    # -------- SpliceAI --------
    spliceai_by_key: dict[str, dict] = {}
    if not args.skip_spliceai:
        write_spliceai_vcf(variants, SPLICEAI_VCF_IN, Path(args.fasta))
        log.info("wrote SpliceAI input VCF: %s", SPLICEAI_VCF_IN)
        try:
            run_spliceai(SPLICEAI_VCF_IN, SPLICEAI_VCF_OUT, Path(args.fasta), args.spliceai_bin)
        except Exception as e:  # noqa: BLE001
            log.error("SpliceAI failed: %s — leaving SpliceAI columns blank", e)
        spliceai_by_key = parse_spliceai_vcf(SPLICEAI_VCF_OUT)
        log.info("SpliceAI parsed for %d variants", len(spliceai_by_key))
    else:
        log.info("SpliceAI skipped")

    # -------- SpliceVault --------
    # Need HGVSc from Genome Nexus step to query SpliceVault by transcript.
    gn_path = PROJECT_ROOT / "out" / "annotated" / "genome_nexus.tsv"
    hgvsc_by_key: dict[str, str] = {}
    if gn_path.exists():
        gn_df = pd.read_csv(gn_path, sep="\t", dtype=str, keep_default_na=False)
        for _, r in gn_df.iterrows():
            k = variant_key(
                r["chromosome"],
                int(r["start_position"]),
                int(r["end_position"]),
                r["reference_allele"],
                r["variant_allele"],
            )
            hgvsc_by_key[k] = r.get("gn_hgvsc", "")
    else:
        log.warning("Genome Nexus output missing — SpliceVault queries will be HGVSc-less")

    # -------- MaxEntScan --------
    score3 = score5 = None
    if not args.skip_maxent:
        try:
            from maxentpy import maxent
            score3 = maxent.score3
            score5 = maxent.score5
        except ImportError as e:
            log.error("maxentpy import failed: %s — leaving MES columns blank", e)

    fasta = None
    if Path(args.fasta).exists():
        from pyfaidx import Fasta
        fasta = Fasta(args.fasta, as_raw=False, sequence_always_upper=True)

    # -------- Merge --------
    rows = []
    n_sv_calls = 0
    for _, r in variants.iterrows():
        chrom = r["chromosome"]
        start = int(r["start_position"])
        end = int(r["end_position"])
        ref = r["reference_allele"]
        alt = r["variant_allele"]
        k = variant_key(chrom, start, end, ref, alt)

        out_row: dict = {
            "chromosome": chrom,
            "start_position": start,
            "end_position": end,
            "reference_allele": ref,
            "variant_allele": alt,
        }

        sa = spliceai_by_key.get(k) or {
            "spliceai_symbol": "",
            "spliceai_ds_ag": float("nan"),
            "spliceai_ds_al": float("nan"),
            "spliceai_ds_dg": float("nan"),
            "spliceai_ds_dl": float("nan"),
            "spliceai_dp_ag": "",
            "spliceai_dp_al": "",
            "spliceai_dp_dg": "",
            "spliceai_dp_dl": "",
            "spliceai_max": float("nan"),
            "spliceai_status": "not_run" if args.skip_spliceai else "no_annotation",
        }
        out_row.update(sa)

        if args.skip_splicevault:
            sv = {
                "splicevault_observed": False,
                "splicevault_top_event": "",
                "splicevault_n_events": 0,
                "splicevault_status": "skipped",
            }
        else:
            sv = query_splicevault(k, hgvsc_by_key.get(k, ""))
            if sv["splicevault_status"] == "ok":
                n_sv_calls += 1
        out_row.update(sv)

        if score3 is not None and score5 is not None and fasta is not None:
            mes = score_maxentscan(chrom, start, ref, alt, fasta, score3, score5)
        else:
            mes = {
                "maxent_donor_ref": float("nan"),
                "maxent_donor_alt": float("nan"),
                "maxent_acceptor_ref": float("nan"),
                "maxent_acceptor_alt": float("nan"),
                "maxent_delta": float("nan"),
                "maxent_status": "skipped",
            }
        out_row.update(mes)
        rows.append(out_row)

    df = pd.DataFrame(rows)
    df["_generated_at"] = pd.Timestamp.utcnow().isoformat()
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)

    n_sa = (df["spliceai_status"] == "scored").sum() if "spliceai_status" in df else 0
    n_sa_strong = (df.get("spliceai_max", pd.Series(dtype=float)) >= 0.5).sum()
    n_mes = (df["maxent_status"] == "ok").sum() if "maxent_status" in df else 0
    log.info(
        "wrote %d rows; SpliceAI scored: %d (≥0.5: %d), SpliceVault hits: %d, MaxEnt scored: %d",
        len(df), n_sa, n_sa_strong, n_sv_calls, n_mes,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
