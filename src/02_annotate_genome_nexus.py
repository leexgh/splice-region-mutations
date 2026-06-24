"""Step 2 — Annotate via Genome Nexus genomic endpoint.

POST batches of variants to `/annotation/genomic` with
`fields=annotation_summary,clinvar,my_variant_info`. Extract canonical
Ensembl transcript, HGVSc/HGVSp_Short, exon, variant classification,
ClinVar variation ID + significance, dbSNP rsid.

Per-variant JSON cached at `out/cache/genome_nexus/{variant_key}.json` so
re-runs are skip-only.
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
from _apiclient import CACHE_ROOT, cached_read, cached_write, variant_key  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
UNIQUE_VARIANTS = PROJECT_ROOT / "out" / "unique_variants.tsv"
OUTPUT_TSV = PROJECT_ROOT / "out" / "annotated" / "genome_nexus.tsv"
GN_API = "https://www.genomenexus.org/annotation/genomic"
FIELDS = "annotation_summary,clinvar,my_variant_info"

log = logging.getLogger("step2")
logging.basicConfig(
    format="[%(name)s] %(message)s", level=logging.INFO, stream=sys.stderr
)


def request_batch(batch: list[dict]) -> list[dict]:
    """One POST against /annotation/genomic with retries on 429/5xx."""
    sleep = 2.0
    for attempt in range(5):
        resp = requests.post(
            GN_API,
            params={"fields": FIELDS},
            headers={"accept": "application/json", "Content-Type": "application/json"},
            json=batch,
            timeout=120,
        )
        if resp.status_code < 500 and resp.status_code != 429:
            if not resp.ok:
                raise RuntimeError(
                    f"GN POST {resp.status_code}: {resp.text[:300]}"
                )
            return resp.json()
        log.warning(
            "  HTTP %s on batch of %d, retry in %.1fs (attempt %d)",
            resp.status_code,
            len(batch),
            sleep,
            attempt + 1,
        )
        time.sleep(sleep)
        sleep *= 2
    raise RuntimeError(f"GN gave up after 5 retries (last={resp.status_code})")


def extract_canonical_summary(payload: dict) -> dict:
    """Pull the fields we care about out of an annotation_summary payload."""
    out = {
        "gn_canonical_transcript_id": "",
        "gn_hugo_symbol": "",
        "gn_hgvsc": "",
        "gn_hgvsp_short": "",
        "gn_variant_classification": "",
        "gn_exon": "",
        "gn_consequence_terms": "",
        "gn_entrez_gene_id": "",
        "gn_refseq": "",
        "gn_protein_position": "",
        "gn_most_severe_consequence": payload.get("most_severe_consequence", ""),
        "gn_successfully_annotated": payload.get("successfully_annotated", False),
        "clinvar_id": "",
        "clinvar_significance": "",
        "clinvar_conflicting": "",
        "dbsnp_rsid": "",
        "gnomad_af_popmax": "",
        "gnomad_af_overall": "",
    }
    summary = payload.get("annotation_summary") or {}
    out["gn_canonical_transcript_id"] = summary.get("canonicalTranscriptId", "")
    canonical = summary.get("transcriptConsequenceSummary") or {}
    if canonical:
        out["gn_hugo_symbol"] = canonical.get("hugoGeneSymbol", "")
        out["gn_hgvsc"] = canonical.get("hgvsc", "")
        out["gn_hgvsp_short"] = canonical.get("hgvspShort", "")
        out["gn_variant_classification"] = canonical.get("variantClassification", "")
        out["gn_exon"] = canonical.get("exon", "")
        out["gn_consequence_terms"] = canonical.get("consequenceTerms", "")
        out["gn_entrez_gene_id"] = canonical.get("entrezGeneId", "")
        out["gn_refseq"] = canonical.get("refSeq", "")
        pp = canonical.get("proteinPosition") or {}
        if pp:
            out["gn_protein_position"] = f"{pp.get('start','')}-{pp.get('end','')}"

    clinvar_ann = ((payload.get("clinvar") or {}).get("annotation")) or {}
    if clinvar_ann:
        out["clinvar_id"] = clinvar_ann.get("clinvarId", "") or ""
        out["clinvar_significance"] = clinvar_ann.get("clinicalSignificance", "") or ""
        out["clinvar_conflicting"] = (
            clinvar_ann.get("conflictingClinicalSignificance", "") or ""
        )

    # my_variant_info contains both dbSNP rsid and gnomAD AFs (subject to mv.info coverage).
    mvi = (payload.get("my_variant_info") or {}).get("annotation") or {}
    if mvi:
        dbsnp = mvi.get("dbsnp") or {}
        out["dbsnp_rsid"] = dbsnp.get("rsid", "") or ""
        gnomad = mvi.get("gnomad_exome") or {}
        af = gnomad.get("af") or {}
        if isinstance(af, dict):
            out["gnomad_af_overall"] = af.get("af", "") or ""
            out["gnomad_af_popmax"] = af.get("af_popmax", "") or ""

    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(UNIQUE_VARIANTS))
    parser.add_argument("--output", default=str(OUTPUT_TSV))
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--delay", type=float, default=0.5, help="delay between batches")
    args = parser.parse_args()

    variants = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    log.info("loaded %d unique variants from %s", len(variants), args.input)

    # Identify variants needing API call (cache miss).
    todo = []
    keys = []
    for _, r in variants.iterrows():
        k = variant_key(
            r["chromosome"],
            int(r["start_position"]),
            int(r["end_position"]),
            r["reference_allele"],
            r["variant_allele"],
        )
        keys.append(k)
        if cached_read("genome_nexus", k) is None:
            todo.append(
                (
                    k,
                    {
                        "chromosome": r["chromosome"],
                        "start": int(r["start_position"]),
                        "end": int(r["end_position"]),
                        "referenceAllele": r["reference_allele"],
                        "variantAllele": r["variant_allele"],
                    },
                )
            )

    log.info("cache: %d hits / %d new", len(variants) - len(todo), len(todo))

    for i in range(0, len(todo), args.batch_size):
        batch = todo[i : i + args.batch_size]
        batch_keys = [k for k, _ in batch]
        batch_payload = [p for _, p in batch]
        log.info(
            "  batch %d-%d / %d (%d items)",
            i,
            i + len(batch),
            len(todo),
            len(batch),
        )
        try:
            results = request_batch(batch_payload)
        except Exception as e:
            log.error("  batch failed: %s", e)
            continue

        # GN returns results in input order (per their docs) — but be defensive
        # and re-key off the `originalVariantQuery` field when available.
        by_query = {}
        for r in results:
            q = r.get("originalVariantQuery") or r.get("variant") or ""
            by_query[q] = r

        for k, payload in zip(batch_keys, batch_payload):
            query_str = (
                f"{payload['chromosome']},{payload['start']},"
                f"{payload['end']},{payload['referenceAllele']},{payload['variantAllele']}"
            )
            matched = by_query.get(query_str)
            if matched is None and len(results) == len(batch_payload):
                matched = results[batch_keys.index(k)]
            if matched is None:
                log.warning("  no GN result for %s", query_str)
                cached_write("genome_nexus", k, {"__not_found__": True})
                continue
            cached_write("genome_nexus", k, matched)

        time.sleep(args.delay)

    # Build the annotated TSV.
    rows = []
    n_missing = 0
    for (_, r), k in zip(variants.iterrows(), keys):
        payload = cached_read("genome_nexus", k)
        out_row = {
            "chromosome": r["chromosome"],
            "start_position": r["start_position"],
            "end_position": r["end_position"],
            "reference_allele": r["reference_allele"],
            "variant_allele": r["variant_allele"],
        }
        if payload is None or payload.get("__not_found__"):
            n_missing += 1
            out_row.update({k: "" for k in extract_canonical_summary({})})
        else:
            out_row.update(extract_canonical_summary(payload))
        rows.append(out_row)

    df = pd.DataFrame(rows)
    df["_generated_at"] = pd.Timestamp.utcnow().isoformat()
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)

    log.info(
        "wrote %d rows to %s (%d missing GN annotation)", len(df), args.output, n_missing
    )
    if not df.empty:
        n_clinvar = (df["clinvar_id"] != "").sum()
        n_rsid = (df["dbsnp_rsid"] != "").sum()
        log.info("  %d have ClinVar IDs, %d have dbSNP rsids", n_clinvar, n_rsid)
    return 0


if __name__ == "__main__":
    sys.exit(main())
