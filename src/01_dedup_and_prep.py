"""Step 1 — Deduplicate splice_region variants and emit unique_variants.tsv.

Reads `data/New_Query_2026_06_02_16_16_54.csv`, deduplicates by
(chromosome, start_position, end_position, reference_allele, variant_allele),
preserves per-variant occurrence counts and cancer-type tallies, normalises
the MSK `-` placeholder for indels into the form downstream tools expect
(using the FASTA anchor base), and writes `out/unique_variants.tsv`.

The 5 join-key columns at the start match the variant_key contract used by
every other annotation step.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_CSV = PROJECT_ROOT / "data" / "New_Query_2026_06_02_16_16_54.csv"
REVUE_TSV = PROJECT_ROOT / "data" / "VUEs.txt"
HG19_FASTA = PROJECT_ROOT / "data" / "refs" / "hg19.fa"
OUT_TSV = PROJECT_ROOT / "out" / "unique_variants.tsv"


def load_revue_keys() -> set[str]:
    """Return a set of `chr,start,end,ref,alt` keys present in reVUE."""
    keys = set()
    if not REVUE_TSV.exists():
        return keys
    with REVUE_TSV.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gl = row.get("genomicLocation", "").strip()
            if gl:
                keys.add(gl)
    return keys


def _anchor_base(chrom: str, pos: int, fasta) -> str:
    """Return the base immediately before `pos` (1-based) on `chrom`."""
    chrom_key = chrom if chrom.startswith("chr") else f"chr{chrom}"
    seq = fasta[chrom_key][pos - 1 : pos].seq
    return seq.upper()


def normalise_indel(
    chrom: str,
    start: int,
    end: int,
    ref: str,
    alt: str,
    fasta,
) -> tuple[int, str, str]:
    """Convert MSK '-' placeholder indel notation to anchor-base form.

    MSK MAF uses ref='-' for insertions and alt='-' for deletions, with
    start/end pointing at the affected base(s). Downstream tools (VCF, SpliceAI,
    Genome Nexus genomic endpoint) expect VCF-style notation: ref and alt both
    begin with the anchor base immediately preceding the indel.

    Returns (new_start, new_ref, new_alt).
    """
    if ref == "-" and alt != "-":
        # Insertion: anchor at base `start`, insertion happens after it.
        # MSK convention has start = base before insertion (no shift needed).
        anchor = _anchor_base(chrom, start, fasta)
        return start, anchor, anchor + alt
    if alt == "-" and ref != "-":
        # Deletion: shift start to the anchor base immediately preceding the
        # deletion (start - 1), prepend it to both ref and alt.
        anchor = _anchor_base(chrom, start - 1, fasta)
        return start - 1, anchor + ref, anchor
    return start, ref, alt


def verify_hg19(fasta) -> None:
    """Spot-check a known hg19 coordinate to confirm the build."""
    # TP53 chr17:7,579,472 in hg19 is the start codon region, ref 'C'.
    # Use a less risky check: chr1:1 should not raise, and length of chr1 should be ~249M.
    chr1_len = len(fasta["chr1"])
    assert 249_000_000 < chr1_len < 250_000_000, f"chr1 length {chr1_len} doesn't match hg19"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default=str(INPUT_CSV),
        help="Path to the input MSK MAF-style CSV.",
    )
    parser.add_argument(
        "--output",
        default=str(OUT_TSV),
        help="Path to write the deduplicated TSV.",
    )
    parser.add_argument(
        "--fasta",
        default=str(HG19_FASTA),
        help="Path to hg19 reference FASTA (for indel anchor-base normalisation).",
    )
    parser.add_argument(
        "--skip-fasta",
        action="store_true",
        help="Skip FASTA-based normalisation (will write raw MSK '-' notation; downstream steps may fail on indels).",
    )
    args = parser.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)
    fasta_path = Path(args.fasta)

    if not in_path.exists():
        print(f"ERROR: input not found: {in_path}", file=sys.stderr)
        return 1

    fasta = None
    if not args.skip_fasta:
        if not fasta_path.exists():
            print(
                f"WARNING: {fasta_path} not found. Falling back to --skip-fasta; "
                "indel rows will keep MSK '-' notation and downstream tools will need it normalised later.",
                file=sys.stderr,
            )
        else:
            from pyfaidx import Fasta

            fasta = Fasta(str(fasta_path), as_raw=False, sequence_always_upper=True)
            verify_hg19(fasta)
            print(f"[step1] hg19 reference verified at {fasta_path}", file=sys.stderr)

    revue_keys = load_revue_keys()
    print(f"[step1] loaded {len(revue_keys)} reVUE genomic-location keys", file=sys.stderr)

    rows = pd.read_csv(in_path, dtype=str, keep_default_na=False)
    print(f"[step1] read {len(rows)} input rows from {in_path.name}", file=sys.stderr)

    # Group by genomic identity.
    grouped: dict[tuple, dict] = defaultdict(
        lambda: {
            "cancer_types": Counter(),
            "samples": [],
            "oncokb_calls": Counter(),
            "highest_levels": Counter(),
            "hugo_symbols": Counter(),
            "hgvsp_short": Counter(),
            "consequences": Counter(),
            "variant_types": Counter(),
            "in_revue": False,
        }
    )

    for _, r in rows.iterrows():
        chrom = str(r["Chromosome"]).strip()
        start = int(r["Start_Position"])
        end = int(r["End_Position"])
        ref = str(r["Reference_Allele"]).strip()
        alt = str(r["Tumor_Seq_Allele2"]).strip()
        key = (chrom, start, end, ref, alt)
        bucket = grouped[key]
        bucket["cancer_types"][r.get("CANCER_TYPE", "")] += 1
        bucket["samples"].append(r.get("Tumor_Sample_Barcode", ""))
        bucket["oncokb_calls"][r.get("ONCOGENIC", "")] += 1
        bucket["highest_levels"][r.get("HIGHEST_LEVEL", "")] += 1
        bucket["hugo_symbols"][r.get("Hugo_Symbol", "")] += 1
        bucket["hgvsp_short"][r.get("HGVSp_Short", "")] += 1
        bucket["consequences"][r.get("Consequence", "")] += 1
        bucket["variant_types"][r.get("Variant_Type", "")] += 1
        if str(r.get("inReVUE", "")).upper() == "TRUE":
            bucket["in_revue"] = True

    print(f"[step1] {len(grouped)} unique variants by chr/start/end/ref/alt", file=sys.stderr)

    output_rows = []
    n_indel_normalised = 0
    for (chrom, start, end, ref, alt), bucket in grouped.items():
        norm_start, norm_ref, norm_alt = start, ref, alt
        if fasta is not None and (ref == "-" or alt == "-"):
            norm_start, norm_ref, norm_alt = normalise_indel(chrom, start, end, ref, alt, fasta)
            n_indel_normalised += 1

        # Also build a reVUE-style genomic_location key: "{chr},{start},{end},{ref},{alt}"
        # using the *original* MSK notation (matches VUEs.txt format).
        msk_genomic_location = f"{chrom},{start},{end},{ref},{alt}"
        in_revue = bucket["in_revue"] or (msk_genomic_location in load_revue_keys())

        most_common = lambda c: c.most_common(1)[0][0] if c else ""

        output_rows.append(
            {
                # 5-column join key (VCF-style after normalisation):
                "chromosome": chrom,
                "start_position": norm_start,
                "end_position": end if (ref != "-" and alt != "-") else norm_start,
                "reference_allele": norm_ref,
                "variant_allele": norm_alt,
                # Original MSK notation kept for traceability and reVUE matching:
                "msk_start_position": start,
                "msk_end_position": end,
                "msk_reference_allele": ref,
                "msk_variant_allele": alt,
                "msk_genomic_location": msk_genomic_location,
                # Variant-level metadata:
                "hugo_symbol": most_common(bucket["hugo_symbols"]),
                "hgvsp_short_input": most_common(bucket["hgvsp_short"]),
                "variant_type": most_common(bucket["variant_types"]),
                "consequence_input": most_common(bucket["consequences"]),
                # Aggregated cohort context:
                "occurrence_count": sum(bucket["cancer_types"].values()),
                "n_distinct_samples": len(set(bucket["samples"])),
                "input_sample_ids": ",".join(sorted(set(s for s in bucket["samples"] if s))),
                "cancer_types_json": json.dumps(dict(bucket["cancer_types"])),
                "oncokb_call_input": most_common(bucket["oncokb_calls"]),
                "oncokb_highest_level_input": most_common(bucket["highest_levels"]),
                # reVUE membership (gold-standard positive anchor):
                "in_revue": bool(in_revue),
            }
        )

    df = pd.DataFrame(output_rows)
    df.sort_values(["hugo_symbol", "chromosome", "start_position"], inplace=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)

    print(f"[step1] wrote {len(df)} unique variants to {out_path}", file=sys.stderr)
    print(f"[step1] {n_indel_normalised} indels normalised from MSK '-' notation", file=sys.stderr)
    print(f"[step1] {df['in_revue'].sum()} variants flagged inReVUE", file=sys.stderr)
    print(f"[step1] top genes: {df['hugo_symbol'].value_counts().head(10).to_dict()}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
