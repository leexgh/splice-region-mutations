"""Step 4b — PubMed paper context search.

For every PMID cited by any source (ClinVar, CIViC, reVUE, LitVar2), fetch the
abstract from NCBI E-utilities and check whether the paper mentions:
  (a) the variant directly (by HGVSc / HGVSp / chr:pos / dbSNP rsid), and
  (b) a cancer context (keyword list).

Per-variant outputs the lists of validated / unvalidated PMIDs and the
manual-review payload (matched pattern + snippet) for the report card.

Used downstream by Step 8 to gate the HIGH `flag_literature_validated`,
MEDIUM `flag_clinvar_plp_validated`, and MEDIUM `flag_civic_validated` flags.

Cache layout (idempotent — never re-fetched):
  out/cache/pubmed_abstract/{pmid}.txt — raw abstract
  out/cache/pubmed_context/{pmid}__{variant_key}.json — per-(pmid,variant)
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _apiclient import CACHE_ROOT, variant_key  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
UNIQUE_VARIANTS = PROJECT_ROOT / "out" / "unique_variants.tsv"
GN_TSV = PROJECT_ROOT / "out" / "annotated" / "genome_nexus.tsv"
CLINICAL_TSV = PROJECT_ROOT / "out" / "annotated" / "clinical.tsv"
LITERATURE_TSV = PROJECT_ROOT / "out" / "annotated" / "literature.tsv"
OUTPUT_TSV = PROJECT_ROOT / "out" / "annotated" / "paper_context.tsv"
ABSTRACT_CACHE_DIR = CACHE_ROOT / "pubmed_abstract"
CONTEXT_CACHE_DIR = CACHE_ROOT / "pubmed_context"

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

log = logging.getLogger("step4b")
logging.basicConfig(
    format="[%(name)s] %(message)s", level=logging.INFO, stream=sys.stderr
)


CANCER_KEYWORDS = re.compile(
    r"\b("
    r"cancer|cancers|tumor|tumors|tumour|tumours|carcinom\w*|"
    r"leukem\w*|leukaem\w*|lymphom\w*|sarcom\w*|melanom\w*|glioma\w*|"
    r"neoplas\w*|oncogen\w*|metasta\w*|malignan\w*|adenocarcinom\w*|"
    r"myelodysplas\w*|myeloprolif\w*|chemo\w*|tumorigen\w*"
    r")\b",
    re.IGNORECASE,
)


def fetch_abstract(pmid: str, delay: float) -> str | None:
    """Fetch and cache a single PMID's abstract via NCBI efetch."""
    ABSTRACT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = ABSTRACT_CACHE_DIR / f"{pmid}.txt"
    if cache_path.exists():
        return cache_path.read_text(errors="replace")
    time.sleep(delay)
    try:
        r = requests.get(
            EFETCH,
            params={
                "db": "pubmed",
                "id": pmid,
                "rettype": "abstract",
                "retmode": "text",
            },
            timeout=30,
        )
    except requests.RequestException as e:
        log.debug("efetch error for PMID %s: %s", pmid, e)
        return None
    if r.status_code != 200 or not r.text.strip():
        return None
    cache_path.write_text(r.text)
    return r.text


def build_variant_patterns(row: dict) -> list[tuple[str, re.Pattern]]:
    """Return a list of (label, pattern) for matching this variant in text.

    Patterns are anchored to specific notation forms. Cheap ones first so the
    first match wins for the manual-review snippet.
    """
    chrom = (row.get("chromosome") or "").lstrip("chr")
    start = row.get("start_position", "")
    hgvsc = row.get("gn_hgvsc", "") or ""
    hgvsp_short = row.get("gn_hgvsp_short", "") or ""
    rsid = row.get("dbsnp_rsid", "") or ""

    patterns: list[tuple[str, re.Pattern]] = []

    # 1. dbSNP rsid is most specific. Avoid false positives by requiring word boundary.
    if rsid and rsid.startswith("rs"):
        patterns.append((f"rsid {rsid}", re.compile(rf"\b{re.escape(rsid)}\b")))

    # 2. Full HGVSc (e.g. ENST00000394351.3:c.345G>A) — match just the c.xxx tail
    if hgvsc and ":" in hgvsc:
        c_tail = hgvsc.split(":", 1)[1]  # "c.345G>A"
        patterns.append((f"HGVSc {c_tail}", re.compile(re.escape(c_tail))))

    # 3. HGVSp_short — the splice placeholder p.X{pos}_splice gives us the protein
    #    position; rebuild it as p.X{pos}, c.{pos}, codon {pos}, residue {pos}
    if hgvsp_short:
        m = re.search(r"p\.X?(\d+)_?splice", hgvsp_short)
        if m:
            pos = m.group(1)
            patterns.append((
                f"protein pos {pos}",
                re.compile(rf"\bp\.\w?{pos}\b|\bcodon\s*{pos}\b|\bresidue\s*{pos}\b"),
            ))

    # 4. chr:pos
    if chrom and start:
        patterns.append((
            f"chr{chrom}:{start}",
            re.compile(rf"chr0*{chrom}[\s:_-]0*{start}\b|\b{chrom}[:.]0*{start}\b"),
        ))

    return patterns


def evaluate_pmid_against_variant(
    pmid: str,
    abstract: str,
    variant_row: dict,
) -> dict:
    """Run the regex match + cancer-keyword check. Returns a dict with
    `mentions_variant`, `mentions_cancer`, `evidence_quality`,
    `matched_pattern`, `snippet`."""
    out = {
        "pmid": pmid,
        "mentions_variant": False,
        "mentions_cancer": False,
        "evidence_quality": "unvalidated",
        "matched_pattern": "",
        "snippet": "",
    }
    if not abstract:
        out["evidence_quality"] = "no_abstract"
        return out

    out["mentions_cancer"] = bool(CANCER_KEYWORDS.search(abstract))

    patterns = build_variant_patterns(variant_row)
    for label, pat in patterns:
        m = pat.search(abstract)
        if not m:
            continue
        out["mentions_variant"] = True
        out["matched_pattern"] = label
        s, e = m.span()
        ctx_start = max(0, s - 100)
        ctx_end = min(len(abstract), e + 100)
        snippet = abstract[ctx_start:ctx_end].replace("\n", " ").strip()
        out["snippet"] = snippet
        break

    if out["mentions_variant"] and out["mentions_cancer"]:
        out["evidence_quality"] = "validated_with_cancer"
    elif out["mentions_variant"]:
        out["evidence_quality"] = "validated_variant_only"
    elif out["mentions_cancer"]:
        out["evidence_quality"] = "cancer_context_only"
    else:
        out["evidence_quality"] = "unvalidated"
    return out


def cached_context(pmid: str, vkey: str) -> dict | None:
    p = CONTEXT_CACHE_DIR / f"{pmid}__{vkey}.json"
    if p.exists():
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            return None
    return None


def write_context(pmid: str, vkey: str, payload: dict) -> None:
    CONTEXT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    p = CONTEXT_CACHE_DIR / f"{pmid}__{vkey}.json"
    p.write_text(json.dumps(payload))


def _pmid_list(s) -> list[str]:
    if not s or pd.isna(s):
        return []
    return [p.strip() for p in str(s).split(",") if p.strip().isdigit()]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(UNIQUE_VARIANTS))
    parser.add_argument("--gn", default=str(GN_TSV))
    parser.add_argument("--clinical", default=str(CLINICAL_TSV))
    parser.add_argument("--literature", default=str(LITERATURE_TSV))
    parser.add_argument("--output", default=str(OUTPUT_TSV))
    parser.add_argument("--delay", type=float, default=0.35, help="delay between efetch calls (NCBI = 3/sec without key)")
    args = parser.parse_args()

    variants = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    log.info("loaded %d unique variants", len(variants))

    gn = pd.read_csv(args.gn, sep="\t", dtype=str, keep_default_na=False)
    gn_idx = {
        variant_key(
            r["chromosome"], int(r["start_position"]), int(r["end_position"]),
            r["reference_allele"], r["variant_allele"],
        ): r.to_dict()
        for _, r in gn.iterrows()
    }

    clinical = pd.read_csv(args.clinical, sep="\t", dtype=str, keep_default_na=False)
    clinical_idx = {
        variant_key(
            r["chromosome"], int(r["start_position"]), int(r["end_position"]),
            r["reference_allele"], r["variant_allele"],
        ): r.to_dict()
        for _, r in clinical.iterrows()
    }

    if Path(args.literature).exists():
        lit = pd.read_csv(args.literature, sep="\t", dtype=str, keep_default_na=False)
        lit_idx = {
            variant_key(
                r["chromosome"], int(r["start_position"]), int(r["end_position"]),
                r["reference_allele"], r["variant_allele"],
            ): r.to_dict()
            for _, r in lit.iterrows()
        }
    else:
        lit_idx = {}

    # Build per-variant PMID source map: (variant_key) -> dict of {pmid: [sources...]}
    variant_pmids: dict[str, dict[str, list[str]]] = {}
    all_pmids: set[str] = set()

    for _, v in variants.iterrows():
        k = variant_key(
            v["chromosome"], int(v["start_position"]), int(v["end_position"]),
            v["reference_allele"], v["variant_allele"],
        )
        per_pmid: dict[str, list[str]] = {}
        cl = clinical_idx.get(k, {})
        for src, col in (("clinvar", "clinvar_pmids"), ("civic", "civic_pmids")):
            for p in _pmid_list(cl.get(col, "")):
                per_pmid.setdefault(p, []).append(src)
        litr = lit_idx.get(k, {})
        for src, col in (("revue", "revue_pmids"), ("litvar2", "litvar_pmids")):
            for p in _pmid_list(litr.get(col, "")):
                per_pmid.setdefault(p, []).append(src)
        variant_pmids[k] = per_pmid
        all_pmids.update(per_pmid.keys())

    log.info("aggregated %d unique PMIDs to validate across %d variants", len(all_pmids), len(variant_pmids))

    # Fetch abstracts (cached; never re-fetched).
    n_fetched = 0
    n_cached = 0
    abstracts: dict[str, str | None] = {}
    for pmid in sorted(all_pmids):
        p = ABSTRACT_CACHE_DIR / f"{pmid}.txt"
        if p.exists():
            abstracts[pmid] = p.read_text(errors="replace")
            n_cached += 1
        else:
            txt = fetch_abstract(pmid, args.delay)
            abstracts[pmid] = txt
            if txt is not None:
                n_fetched += 1
    log.info("abstracts: %d cache hits, %d new fetches", n_cached, n_fetched)

    # Run context check per (variant, pmid).
    rows = []
    for _, v in variants.iterrows():
        k = variant_key(
            v["chromosome"], int(v["start_position"]), int(v["end_position"]),
            v["reference_allele"], v["variant_allele"],
        )
        variant_full = {**v.to_dict(), **gn_idx.get(k, {})}
        per_pmid = variant_pmids.get(k, {})

        validated_with_cancer: list[str] = []
        validated_variant_only: list[str] = []
        unvalidated: list[str] = []
        manual_review: list[dict] = []

        clinvar_pmids = set(_pmid_list(clinical_idx.get(k, {}).get("clinvar_pmids", "")))
        civic_pmids = set(_pmid_list(clinical_idx.get(k, {}).get("civic_pmids", "")))

        for pmid in per_pmid:
            res = cached_context(pmid, k)
            if res is None:
                abstract = abstracts.get(pmid)
                res = evaluate_pmid_against_variant(pmid, abstract or "", variant_full)
                write_context(pmid, k, res)

            quality = res.get("evidence_quality", "unvalidated")
            if quality == "validated_with_cancer":
                validated_with_cancer.append(pmid)
            elif quality == "validated_variant_only":
                validated_variant_only.append(pmid)
            else:
                unvalidated.append(pmid)

            # Manual-review payload for ClinVar / CIViC validated hits.
            if quality == "validated_with_cancer" and (pmid in clinvar_pmids or pmid in civic_pmids):
                manual_review.append({
                    "pmid": pmid,
                    "source": "ClinVar" if pmid in clinvar_pmids else "CIViC",
                    "pattern": res.get("matched_pattern", ""),
                    "snippet": res.get("snippet", ""),
                })

        clinvar_paper_validation = "n/a"
        if clinvar_pmids:
            clinvar_paper_validation = "yes" if any(p in validated_with_cancer for p in clinvar_pmids) else "no"

        rows.append({
            "chromosome": v["chromosome"],
            "start_position": v["start_position"],
            "end_position": v["end_position"],
            "reference_allele": v["reference_allele"],
            "variant_allele": v["variant_allele"],
            "n_pmids_total": len(per_pmid),
            "validated_with_cancer_pmids": ",".join(validated_with_cancer),
            "validated_variant_only_pmids": ",".join(validated_variant_only),
            "unvalidated_pmids": ",".join(unvalidated),
            "clinvar_paper_validation": clinvar_paper_validation,
            "manual_review_papers_json": json.dumps(manual_review) if manual_review else "",
        })

    df = pd.DataFrame(rows)
    df["_generated_at"] = pd.Timestamp.utcnow().isoformat()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)

    n_validated = (df["validated_with_cancer_pmids"] != "").sum()
    n_clinvar_yes = (df["clinvar_paper_validation"] == "yes").sum()
    n_with_manual = (df["manual_review_papers_json"] != "").sum()
    log.info(
        "wrote %d rows; %d variants have ≥1 validated_with_cancer PMID; "
        "%d ClinVar-paper-validated; %d variants surfaced for manual review",
        len(df), n_validated, n_clinvar_yes, n_with_manual,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
