"""Step 4 — Clinical evidence (CIViC + ClinVar PMIDs + OncoKB gene flags + COSMIC stub).

Reads `out/annotated/genome_nexus.tsv` for ClinVar variation IDs and HGVSc.
Then for each variant:

- **ClinVar PMIDs:** look up the submitter PMID list for `clinvar_id` in
  `data/refs/clinvar_submission_summary.txt.gz` (one-pass gzip scan).
- **CIViC:** query the GraphQL/REST endpoint by gene_symbol + HGVSc; cache per
  variant; capture top evidence items with PMIDs and evidence levels.
- **OncoKB Cancer Gene List:** static JSON download → join on hugo_symbol to
  set `gene_is_tsg`, `gene_is_oncogene`, `gene_is_oncokb_annotated`. **OncoKB
  per-variant oncogenicity calls are NOT pulled.**
COSMIC removed — user has no license. No axis, no stub.

Output: `out/annotated/clinical.tsv`.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _apiclient import cached_read, cached_write, variant_key  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GN_TSV = PROJECT_ROOT / "out" / "annotated" / "genome_nexus.tsv"
UNIQUE_VARIANTS = PROJECT_ROOT / "out" / "unique_variants.tsv"
OUTPUT_TSV = PROJECT_ROOT / "out" / "annotated" / "clinical.tsv"

CLINVAR_SUBM = PROJECT_ROOT / "data" / "refs" / "clinvar_submission_summary.txt.gz"
ONCOKB_GENES = PROJECT_ROOT / "data" / "refs" / "oncokb_cancer_gene_list.json"

CIVIC_API = "https://civicdb.org/api/graphql"

log = logging.getLogger("step4")
logging.basicConfig(
    format="[%(name)s] %(message)s", level=logging.INFO, stream=sys.stderr
)


# --------------- ClinVar PMIDs ---------------

# Submitters write PMIDs as "PMID:12345", "PMID: 12345", "PubMed:12345",
# "PMID 12345", and comma-separated lists like "PMID: 17576681, 9536098".
# We catch the leading prefix once, then collect every digit run inside the
# same parenthetical or sentence.
_PMID_PATTERN = re.compile(r"(?:PMID|PubMed)\s*[:\s]\s*([\d,\s]+)", re.IGNORECASE)
_DIGITS_PATTERN = re.compile(r"\d{4,9}")


_REVIEW_RANK = {
    # Higher number = more authoritative (more stars).
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no assertion provided": 0,
}


def load_clinvar_pmids_and_status(needed_ids: set[str]) -> dict[str, dict]:
    """One-pass scan of submission_summary.txt.gz, collecting both PMIDs and the
    highest review-status seen per ClinVar VariationID.

    Returns {varid: {"pmids": [...], "review_status": "...", "review_rank": int}}.
    """
    out: dict[str, dict] = {
        k: {"pmids": set(), "review_status": "", "review_rank": -1}
        for k in needed_ids
    }
    if not CLINVAR_SUBM.exists() or not needed_ids:
        log.warning("ClinVar submission_summary not found or no IDs needed")
        return out
    log.info("scanning ClinVar submission_summary for %d variation IDs", len(needed_ids))
    n_rows = 0
    with gzip.open(CLINVAR_SUBM, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            varid = cols[0]
            if varid not in needed_ids:
                continue
            n_rows += 1
            bucket = out[varid]
            # PMIDs anywhere in the row.
            for col in cols[1:]:
                for m in _PMID_PATTERN.finditer(col):
                    chunk = m.group(1)
                    for d in _DIGITS_PATTERN.findall(chunk):
                        bucket["pmids"].add(d)
            # Review status is in column 6 (0-indexed).
            review_status = cols[6].strip().lower()
            rank = _REVIEW_RANK.get(review_status, 0)
            if rank > bucket["review_rank"]:
                bucket["review_status"] = review_status
                bucket["review_rank"] = rank
    n_pmids = sum(1 for v in out.values() if v["pmids"])
    n_status = sum(1 for v in out.values() if v["review_rank"] >= 0)
    log.info(
        "ClinVar scan: %d submission rows; %d IDs with ≥1 PMID; %d with review-status",
        n_rows, n_pmids, n_status,
    )
    for v in out.values():
        v["pmids"] = sorted(v["pmids"])
    return out


# --------------- OncoKB Cancer Gene List ---------------


def load_oncokb_genes() -> dict[str, dict]:
    """Load the JSON gene list → flag dict keyed by hugoSymbol."""
    out: dict[str, dict] = {}
    if not ONCOKB_GENES.exists():
        log.warning("OncoKB cancer gene list missing — gene_is_tsg/oncogene will be blank")
        return out
    with ONCOKB_GENES.open() as f:
        data = json.load(f)
    for g in data:
        symbol = g.get("hugoSymbol")
        if not symbol:
            continue
        gtype = (g.get("geneType") or "").upper()
        out[symbol] = {
            "gene_is_tsg": gtype in ("TSG", "ONCOGENE_AND_TSG"),
            "gene_is_oncogene": gtype in ("ONCOGENE", "ONCOGENE_AND_TSG"),
            "gene_is_oncokb_annotated": bool(g.get("oncokbAnnotated")),
            "gene_type_oncokb": gtype,
        }
    log.info("loaded OncoKB gene-type flags for %d genes", len(out))
    return out


# --------------- CIViC ---------------
# CIViC updated its GraphQL schema (2025+): search() no longer takes a `type`
# argument and returns a flat list of SearchResult{id,name,resultType}.
# Evidence is now on MolecularProfile, accessed via variant.singleVariantMolecularProfile.

CIVIC_SEARCH_QUERY = """
query VariantSearch($name: String!) {
  search(query: $name) {
    id
    name
    resultType
  }
}
"""

CIVIC_VARIANT_QUERY = """
query VariantEvidence($id: Int!) {
  variant(id: $id) {
    id
    name
    variantAliases
    singleVariantMolecularProfile {
      id
      evidenceItems {
        nodes {
          evidenceLevel
          evidenceDirection
          significance
          status
          source { citationId sourceType }
        }
      }
    }
  }
}
"""


def query_civic(gene: str, hgvsp_short: str, hgvsc: str, key: str, delay: float = 0.25) -> dict:
    """Look up a variant in CIViC using the updated two-step GraphQL API.

    Step 1: free-text search → get variant IDs.
    Step 2: fetch evidence for each matching variant.

    Search terms tried (in order):
      - "{gene} {hgvsc_bare}"  (bare c. notation, most likely to match)
      - "{gene} {hgvsp_short}" (MSK placeholder, rarely matches)
    """
    cached = cached_read("civic", key)
    if cached is not None:
        return _summarise_civic(cached)

    import re as _re
    hgvsc_bare = _re.sub(r'^ENS[A-Z]\d+[^:]*:', '', hgvsc or '').strip()
    queries: list[str] = []
    if gene:
        if hgvsc_bare:
            queries.append(f"{gene} {hgvsc_bare}")
        queries.append(f"{gene} {hgvsp_short}")

    ev_nodes: list[dict] = []
    matched_variant: dict = {}

    for q in queries:
        if not q.strip():
            continue
        time.sleep(delay)
        try:
            r = requests.post(
                CIVIC_API,
                json={"query": CIVIC_SEARCH_QUERY, "variables": {"name": q}},
                timeout=20,
            )
            if r.status_code != 200:
                continue
            results = (r.json().get("data") or {}).get("search") or []
            variant_ids = [res["id"] for res in results if res.get("resultType") == "VARIANT"]
            for vid in variant_ids[:3]:
                time.sleep(delay)
                vr = requests.post(
                    CIVIC_API,
                    json={"query": CIVIC_VARIANT_QUERY, "variables": {"id": vid}},
                    timeout=20,
                )
                if vr.status_code != 200:
                    continue
                vdata = (vr.json().get("data") or {}).get("variant") or {}
                nodes = ((vdata.get("singleVariantMolecularProfile") or {}).get("evidenceItems") or {}).get("nodes") or []
                if nodes:
                    ev_nodes = nodes
                    matched_variant = vdata
                    break
            if ev_nodes:
                break
        except requests.RequestException:
            continue

    cache_payload = {
        "data": {
            "ev_nodes": ev_nodes,
            "variant": matched_variant,
        }
    }
    cached_write("civic", key, cache_payload)
    return _summarise_civic(cache_payload)


def _summarise_civic(payload: dict) -> dict:
    data = payload.get("data") or {}
    ev_nodes = data.get("ev_nodes") or []
    variant = data.get("variant") or {}

    no_match = {
        "civic_variant_id": "",
        "civic_variant_name": "",
        "civic_evidence_count": 0,
        "civic_best_level": "",
        "civic_significance_summary": "",
        "civic_pmids": "",
        "civic_status": "no_match" if "data" in payload else "error",
    }
    if not ev_nodes:
        return no_match

    pmids: list[str] = []
    levels: list[str] = []
    sigs: list[str] = []
    for ei in ev_nodes:
        if ei.get("status") != "ACCEPTED":
            continue
        if ei.get("evidenceLevel"):
            levels.append(ei["evidenceLevel"])
        if ei.get("significance"):
            sigs.append(ei["significance"])
        src = ei.get("source") or {}
        if src.get("sourceType") == "PUBMED" and src.get("citationId"):
            pmids.append(str(src["citationId"]))

    best_level = next((L for L in ("A", "B", "C", "D", "E") if L in levels), "")
    from collections import Counter
    return {
        "civic_variant_id": variant.get("id", ""),
        "civic_variant_name": variant.get("name", ""),
        "civic_evidence_count": len(ev_nodes),
        "civic_best_level": best_level,
        "civic_significance_summary": ",".join(f"{s}:{n}" for s, n in Counter(sigs).items()),
        "civic_pmids": ",".join(sorted(set(pmids))),
        "civic_status": "ok",
    }


# --------------- main ---------------


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(UNIQUE_VARIANTS))
    parser.add_argument("--gn", default=str(GN_TSV))
    parser.add_argument("--output", default=str(OUTPUT_TSV))
    parser.add_argument("--skip-civic", action="store_true")
    parser.add_argument("--skip-clinvar-pmids", action="store_true")
    args = parser.parse_args()

    variants = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    log.info("loaded %d unique variants", len(variants))

    if not Path(args.gn).exists():
        log.error("Genome Nexus output missing at %s — run Step 2 first", args.gn)
        return 1
    gn = pd.read_csv(args.gn, sep="\t", dtype=str, keep_default_na=False)
    log.info("loaded Genome Nexus annotation for %d variants", len(gn))

    # Build per-variant index keyed by variant_key for joining.
    gn_idx: dict[str, dict] = {}
    for _, r in gn.iterrows():
        k = variant_key(
            r["chromosome"], int(r["start_position"]), int(r["end_position"]),
            r["reference_allele"], r["variant_allele"],
        )
        gn_idx[k] = r.to_dict()

    # ClinVar PMID prefetch.
    needed_ids: set[str] = set()
    for r in gn_idx.values():
        cv = (r.get("clinvar_id") or "").strip()
        if cv and cv != "nan":
            needed_ids.add(str(cv))
    if args.skip_clinvar_pmids:
        clinvar_data: dict[str, dict] = {
            k: {"pmids": [], "review_status": "", "review_rank": -1} for k in needed_ids
        }
        log.info("ClinVar PMID/review-status scan skipped")
    else:
        clinvar_data = load_clinvar_pmids_and_status(needed_ids)

    oncokb_genes = load_oncokb_genes()

    rows = []
    n_civic_hits = 0
    for i, (_, v) in enumerate(variants.iterrows()):
        chrom = v["chromosome"]
        start = int(v["start_position"])
        end = int(v["end_position"])
        ref = v["reference_allele"]
        alt = v["variant_allele"]
        k = variant_key(chrom, start, end, ref, alt)
        gn_row = gn_idx.get(k, {})

        gene = v.get("hugo_symbol") or gn_row.get("gn_hugo_symbol", "") or ""
        hgvsp_short = gn_row.get("gn_hgvsp_short", "") or v.get("hgvsp_short_input", "")
        hgvsc = gn_row.get("gn_hgvsc", "") or ""
        clinvar_id = gn_row.get("clinvar_id", "") or ""

        out_row: dict = {
            "chromosome": chrom,
            "start_position": start,
            "end_position": end,
            "reference_allele": ref,
            "variant_allele": alt,
            "clinvar_id": clinvar_id,
            "clinvar_significance": gn_row.get("clinvar_significance", "") or "",
            "clinvar_conflicting": gn_row.get("clinvar_conflicting", "") or "",
            "clinvar_pmids": ",".join(
                clinvar_data.get(str(clinvar_id), {}).get("pmids", []) if clinvar_id else []
            ),
            "clinvar_review_status": (
                clinvar_data.get(str(clinvar_id), {}).get("review_status", "") if clinvar_id else ""
            ),
            "clinvar_review_rank": (
                clinvar_data.get(str(clinvar_id), {}).get("review_rank", -1) if clinvar_id else -1
            ),
        }

        # OncoKB Cancer Gene flags (TSG/oncogene only, no oncogenicity call).
        gene_flags = oncokb_genes.get(gene)
        if gene_flags:
            out_row.update({
                "gene_is_tsg": gene_flags["gene_is_tsg"],
                "gene_is_oncogene": gene_flags["gene_is_oncogene"],
                "gene_is_oncokb_annotated": gene_flags["gene_is_oncokb_annotated"],
                "gene_type_oncokb": gene_flags["gene_type_oncokb"],
            })
        else:
            out_row.update({
                "gene_is_tsg": False,
                "gene_is_oncogene": False,
                "gene_is_oncokb_annotated": False,
                "gene_type_oncokb": "",
            })

        # CIViC.
        if args.skip_civic:
            out_row.update({
                "civic_variant_id": "",
                "civic_variant_name": "",
                "civic_evidence_count": 0,
                "civic_best_level": "",
                "civic_significance_summary": "",
                "civic_pmids": "",
                "civic_status": "skipped",
            })
        else:
            civic = query_civic(gene, hgvsp_short, hgvsc, k)
            out_row.update(civic)
            if civic["civic_status"] == "ok":
                n_civic_hits += 1

        rows.append(out_row)

        if (i + 1) % 100 == 0:
            log.info("  progress %d/%d (civic hits=%d)", i + 1, len(variants), n_civic_hits)

    df = pd.DataFrame(rows)
    df["_generated_at"] = pd.Timestamp.utcnow().isoformat()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    n_tsg = df["gene_is_tsg"].sum()
    n_onc = df["gene_is_oncogene"].sum()
    n_pmids = (df["clinvar_pmids"] != "").sum()
    log.info(
        "wrote %d rows; %d TSG / %d oncogene by OncoKB CGL, %d variants with ClinVar PMIDs, %d CIViC hits",
        len(df), n_tsg, n_onc, n_pmids, n_civic_hits,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
