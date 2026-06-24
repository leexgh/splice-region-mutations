"""Step 6 — Population frequency via gnomAD GraphQL API.

For each variant, query gnomAD v4 (joint exomes + genomes) for AF_popmax,
homozygote counts, and FILTER status. Falls back to gnomAD v2 (hg19) if v4
can't lift the coordinate.

Results cached per-variant at `out/cache/gnomad/{variant_key}.json`.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _apiclient import cached_read, cached_write, variant_key  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
UNIQUE_VARIANTS = PROJECT_ROOT / "out" / "unique_variants.tsv"
OUTPUT_TSV = PROJECT_ROOT / "out" / "annotated" / "population.tsv"

GNOMAD_API = "https://gnomad.broadinstitute.org/api"

log = logging.getLogger("step6")
logging.basicConfig(
    format="[%(name)s] %(message)s", level=logging.INFO, stream=sys.stderr
)

GRAPHQL_V4 = """
query VariantSearch($variantId: String!) {
  variant(variantId: $variantId, dataset: gnomad_r4) {
    variantId
    rsid
    joint { ac an filters populations { id ac an } }
  }
}
"""

# Fallback for variants not present in v4 (hg19 coords still useful).
GRAPHQL_V2 = """
query VariantSearchV2($variantId: String!) {
  variant(variantId: $variantId, dataset: gnomad_r2_1) {
    variantId
    rsid
    exome { ac an filters populations { id ac an } }
  }
}
"""


def _af(group: dict) -> float:
    """Compute allele frequency from ac/an, or NaN if missing."""
    if not group:
        return float("nan")
    ac = group.get("ac")
    an = group.get("an")
    if ac is None or not an:
        return float("nan")
    return ac / an


def gnomad_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """gnomAD GraphQL expects 'CHROM-POS-REF-ALT' with no 'chr' prefix."""
    chrom_clean = chrom[3:] if chrom.startswith("chr") else chrom
    return f"{chrom_clean}-{pos}-{ref}-{alt}"


def _post_graphql(query: str, variables: dict, delay: float, max_retries: int = 8) -> dict:
    sleep = 4.0
    last = None
    for attempt in range(max_retries):
        time.sleep(delay)
        try:
            resp = requests.post(
                GNOMAD_API,
                json={"query": query, "variables": variables},
                headers={"Content-Type": "application/json"},
                timeout=60,
            )
        except requests.exceptions.Timeout:
            log.warning(
                "  gnomAD timeout, retry in %.1fs (attempt %d)",
                sleep, attempt + 1,
            )
            time.sleep(sleep)
            sleep = min(sleep * 2, 120)
            continue
        if resp.status_code == 200:
            return resp.json()
        if resp.status_code in (429, 502, 503, 504):
            log.warning(
                "  gnomAD HTTP %s, retry in %.1fs (attempt %d)",
                resp.status_code, sleep, attempt + 1,
            )
            time.sleep(sleep)
            sleep = min(sleep * 2, 120)
            last = resp
            continue
        # Other error — return raw
        return {"errors": [{"message": f"HTTP {resp.status_code}: {resp.text[:200]}"}]}
    return {"errors": [{"message": f"max retries exceeded (last={last.status_code if last else 'timeout'})"}]}


def af_popmax_from_populations(pops: list) -> float:
    """Compute AF_popmax over non-bottlenecked populations.

    gnomAD's official popmax excludes Ashkenazi Jewish and 'oth' (Other).
    For simplicity we compute it client-side here from the per-population AC/AN.
    """
    excluded = {"asj", "oth", "remaining"}
    best = float("nan")
    for pop in pops or []:
        pop_id = (pop.get("id") or "").lower()
        if pop_id in excluded:
            continue
        # exclude sex- and PCR-split groups (gnomAD returns subpopulation breakdowns too)
        if pop_id.count("_") > 0:
            continue
        an = pop.get("an") or 0
        ac = pop.get("ac") or 0
        if an >= 200:
            af = ac / an
            if pd.isna(best) or af > best:
                best = af
    return best


def parse_payload(payload: dict) -> dict:
    """Flatten gnomAD GraphQL response into our column schema."""
    out = {
        "gnomad_dataset": "",
        "gnomad_variant_id": "",
        "gnomad_rsid": "",
        "gnomad_exome_af": "",
        "gnomad_genome_af": "",
        "gnomad_joint_af": "",
        "gnomad_af_popmax": "",
        "gnomad_filters": "",
        "gnomad_status": "not_found",
    }
    if not payload:
        out["gnomad_status"] = "error"
        return out
    # gnomAD GraphQL returns errors=[{"message":"Variant not found"}] when the
    # variant is genuinely absent from the dataset — that's signal, not failure.
    errs = payload.get("errors") or []
    is_notfound = any(
        "not found" in (
            (e.get("message") or "").lower() if isinstance(e, dict) else str(e).lower()
        )
        for e in errs
    )
    if errs and not is_notfound:
        out["gnomad_status"] = "error"
        return out
    data = (payload.get("data") or {}).get("variant")
    if data is None:
        # Explicit "Variant not found" or null variant — legitimate absence.
        return out

    out["gnomad_variant_id"] = data.get("variantId", "") or ""
    out["gnomad_rsid"] = data.get("rsid", "") or ""

    exome = data.get("exome") or {}
    genome = data.get("genome") or {}
    joint = data.get("joint") or {}

    af_e = _af(exome)
    af_g = _af(genome)
    af_j = _af(joint)
    out["gnomad_exome_af"] = af_e if pd.notna(af_e) else ""
    out["gnomad_genome_af"] = af_g if pd.notna(af_g) else ""
    out["gnomad_joint_af"] = af_j if pd.notna(af_j) else ""

    # In v4 we only fetched joint; in v2 fallback we only fetched exome.
    # Make sure the primary AF column is populated.

    # Popmax from joint when present (v4), else exome (v2).
    pops_source = (joint or exome or {}).get("populations") or []
    popmax = af_popmax_from_populations(pops_source)
    out["gnomad_af_popmax"] = popmax if pd.notna(popmax) else ""

    filters = []
    for src in (exome, genome, joint):
        if src and src.get("filters"):
            filters.extend(src["filters"])
    out["gnomad_filters"] = ",".join(sorted(set(filters))) if filters else ""

    out["gnomad_status"] = "ok"
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(UNIQUE_VARIANTS))
    parser.add_argument("--output", default=str(OUTPUT_TSV))
    parser.add_argument("--delay", type=float, default=1.0, help="delay between API calls")
    parser.add_argument("--no-v2-fallback", action="store_true")
    args = parser.parse_args()

    variants = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    log.info("loaded %d unique variants", len(variants))

    rows = []
    n_cache_hits = 0
    n_v4_hits = 0
    n_v2_hits = 0
    n_not_found = 0
    n_errors = 0

    for i, (_, r) in enumerate(variants.iterrows()):
        chrom = r["chromosome"]
        pos = int(r["start_position"])
        end = int(r["end_position"])
        ref = r["reference_allele"]
        alt = r["variant_allele"]
        k = variant_key(chrom, pos, end, ref, alt)
        var_id = gnomad_variant_id(chrom, pos, ref, alt)

        cached = cached_read("gnomad", k)
        payload = cached
        dataset = ""
        if payload is None:
            # First try v4 (default).
            payload = _post_graphql(GRAPHQL_V4, {"variantId": var_id}, args.delay)
            dataset = "gnomad_r4"
            if (
                not (payload.get("data") or {}).get("variant")
                and not args.no_v2_fallback
            ):
                payload_v2 = _post_graphql(GRAPHQL_V2, {"variantId": var_id}, args.delay)
                if (payload_v2.get("data") or {}).get("variant"):
                    payload = payload_v2
                    dataset = "gnomad_r2_1"
            payload["__dataset__"] = dataset
            # Only cache responses we trust: real data, or an explicit "Variant
            # not found". Transient 429s / max-retry errors should NOT be cached,
            # so the next run can retry them.
            errs = payload.get("errors") or []
            should_cache = False
            if (payload.get("data") or {}).get("variant") is not None:
                should_cache = True
            else:
                for e in errs:
                    msg = (e.get("message") or "").lower() if isinstance(e, dict) else str(e).lower()
                    if "variant not found" in msg or "not found" in msg:
                        should_cache = True
                        break
            if should_cache:
                cached_write("gnomad", k, payload)
            else:
                log.warning("  not caching transient failure for %s", var_id)
        else:
            n_cache_hits += 1
            dataset = payload.get("__dataset__", "")

        parsed = parse_payload(payload)
        parsed["gnomad_dataset"] = dataset
        row = {
            "chromosome": chrom,
            "start_position": pos,
            "end_position": end,
            "reference_allele": ref,
            "variant_allele": alt,
        }
        row.update(parsed)
        rows.append(row)

        if parsed["gnomad_status"] == "ok":
            if dataset == "gnomad_r4":
                n_v4_hits += 1
            else:
                n_v2_hits += 1
        elif parsed["gnomad_status"] == "not_found":
            n_not_found += 1
        else:
            n_errors += 1

        if (i + 1) % 100 == 0:
            log.info(
                "  progress %d/%d (cache=%d v4=%d v2=%d nf=%d err=%d)",
                i + 1, len(variants), n_cache_hits, n_v4_hits, n_v2_hits, n_not_found, n_errors,
            )

    df = pd.DataFrame(rows)
    df["_generated_at"] = pd.Timestamp.utcnow().isoformat()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    log.info(
        "wrote %d rows; cache hits %d, v4 hits %d, v2 fallback %d, not found %d, errors %d",
        len(df), n_cache_hits, n_v4_hits, n_v2_hits, n_not_found, n_errors,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
