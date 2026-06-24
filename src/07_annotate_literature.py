"""Step 7 — Literature aggregation.

Collects PMIDs from:
- reVUE (`data/VUEs.txt` — direct `pubmedId` column for inReVUE variants)
- LitVar2 NCBI API (best-effort lookup by gene + HGVSc; cached per variant)

The scoring axis in Step 8 dedupes across reVUE + ClinVar + CIViC + LitVar2.
Output `out/annotated/literature.tsv`.
"""

from __future__ import annotations

import argparse
import csv
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
GN_TSV = PROJECT_ROOT / "out" / "annotated" / "genome_nexus.tsv"
REVUE_TSV = PROJECT_ROOT / "data" / "VUEs.txt"
OUTPUT_TSV = PROJECT_ROOT / "out" / "annotated" / "literature.tsv"

LITVAR2_API = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/"

log = logging.getLogger("step7")
logging.basicConfig(
    format="[%(name)s] %(message)s", level=logging.INFO, stream=sys.stderr
)


def load_revue_index() -> dict[str, dict]:
    """Map MSK-style genomic_location -> reVUE row."""
    out: dict[str, dict] = {}
    if not REVUE_TSV.exists():
        return out
    with REVUE_TSV.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gl = (row.get("genomicLocation") or "").strip()
            if gl:
                out[gl] = row
    return out


def query_litvar2(gene: str, hgvsc: str, key: str, delay: float = 0.4) -> list[str]:
    """Best-effort LitVar2 autocomplete lookup. Returns PMID list."""
    if not gene:
        return []
    q = f"{gene} {hgvsc}".strip() if hgvsc else gene
    cached = cached_read("litvar2", key)
    if cached is not None:
        return cached.get("pmids", []) if isinstance(cached, dict) else []
    time.sleep(delay)
    try:
        r = requests.get(LITVAR2_API, params={"query": q}, timeout=15)
    except requests.RequestException as e:
        log.debug("  litvar2 fetch error: %s", e)
        return []
    if r.status_code != 200:
        return []
    data = r.json() if r.text.strip() else []
    pmids: set[str] = set()
    for entry in data if isinstance(data, list) else []:
        for p in entry.get("pmids", []) or []:
            pmids.add(str(p))
    payload = {"pmids": sorted(pmids), "query": q}
    cached_write("litvar2", key, payload)
    return payload["pmids"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(UNIQUE_VARIANTS))
    parser.add_argument("--gn", default=str(GN_TSV))
    parser.add_argument("--output", default=str(OUTPUT_TSV))
    parser.add_argument("--skip-litvar2", action="store_true")
    args = parser.parse_args()

    variants = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    revue_idx = load_revue_index()
    log.info("loaded %d variants, %d reVUE entries", len(variants), len(revue_idx))

    gn_hgvsc: dict[str, str] = {}
    if Path(args.gn).exists():
        gn = pd.read_csv(args.gn, sep="\t", dtype=str, keep_default_na=False)
        for _, r in gn.iterrows():
            k = variant_key(
                r["chromosome"], int(r["start_position"]), int(r["end_position"]),
                r["reference_allele"], r["variant_allele"],
            )
            gn_hgvsc[k] = r.get("gn_hgvsc", "")

    rows = []
    for i, (_, v) in enumerate(variants.iterrows()):
        chrom = v["chromosome"]
        start = int(v["start_position"])
        end = int(v["end_position"])
        ref = v["reference_allele"]
        alt = v["variant_allele"]
        k = variant_key(chrom, start, end, ref, alt)

        # reVUE — match by MSK genomic_location.
        msk_gl = v.get("msk_genomic_location", "")
        revue_row = revue_idx.get(msk_gl) or {}
        revue_pmid = revue_row.get("pubmedId", "").strip() if revue_row else ""
        revue_pmids = [p.strip() for p in revue_pmid.replace(";", ",").split(",") if p.strip().isdigit()]

        # LitVar2.
        litvar_pmids: list[str] = []
        if not args.skip_litvar2:
            hgvsc = gn_hgvsc.get(k, "") or v.get("hgvsp_short_input", "")
            try:
                litvar_pmids = query_litvar2(v.get("hugo_symbol", ""), hgvsc, k)
            except Exception as e:  # noqa: BLE001
                log.debug("litvar2 error for %s: %s", k, e)

        rows.append({
            "chromosome": chrom,
            "start_position": start,
            "end_position": end,
            "reference_allele": ref,
            "variant_allele": alt,
            "revue_pmids": ",".join(revue_pmids),
            "litvar_pmids": ",".join(litvar_pmids),
            "literature_status": "ok",
        })
        if (i + 1) % 200 == 0:
            log.info("  progress %d/%d", i + 1, len(variants))

    df = pd.DataFrame(rows)
    df["_generated_at"] = pd.Timestamp.utcnow().isoformat()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    n_revue = (df["revue_pmids"] != "").sum()
    n_litvar = (df["litvar_pmids"] != "").sum()
    log.info("wrote %d rows; %d with reVUE PMIDs, %d with LitVar2 PMIDs", len(df), n_revue, n_litvar)
    return 0


if __name__ == "__main__":
    sys.exit(main())
