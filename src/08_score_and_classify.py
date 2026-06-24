"""Step 8 — Evidence flags + 3-tier classification (no continuous score).

Joins all `out/annotated/*.tsv` on the 5-column variant_key, applies the
evidence-flag rules from the plan, derives a 3-tier label (Likely / Maybe /
Unlikely), and emits `out/variants_scored.tsv`. Every boolean flag is stored
as its own column so the rule can be retuned without re-running annotation.

Weight order (user-supplied):
  reVUE ≥ literature-validated >
  cBio (MSK+GENIE) high-recurrence >
  {high-class ClinVar P/LP, CIViC Level A/B} (paper-validated) >
  {low-class ClinVar P/LP, lower-level CIViC} (unvalidated) >
  {SpliceAI, SpliceVault, gnomAD-rare}.

OncoKB calls are explicitly NOT consumed anywhere.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
UNIQUE_VARIANTS = PROJECT_ROOT / "out" / "unique_variants.tsv"
ANNOTATED_DIR = PROJECT_ROOT / "out" / "annotated"
OUTPUT_TSV = PROJECT_ROOT / "out" / "variants_scored.tsv"

log = logging.getLogger("step8")
logging.basicConfig(
    format="[%(name)s] %(message)s", level=logging.INFO, stream=sys.stderr
)

JOIN_KEY = ["chromosome", "start_position", "end_position", "reference_allele", "variant_allele"]


def _bool(x) -> bool:
    return str(x).strip().lower() in ("true", "1")


def _to_float(x) -> float:
    try:
        if x in ("", None):
            return float("nan")
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _to_int(x) -> int | None:
    try:
        if x in ("", None):
            return None
        return int(float(x))
    except (TypeError, ValueError):
        return None


# Tunable thresholds — change here, re-run Step 8 only.
THRESHOLDS = {
    "cbio_high_recurrence": 5,    # cbio_total_n_samples ≥ this → MEDIUM flag
    "spliceai_strong": 0.5,        # spliceai_max ≥ this → LOW flag
    "gnomad_rare_max": 1e-4,       # AF_popmax < this (or absent) → LOW rare flag
    "gnomad_common_min": 1e-3,     # AF_popmax > this → VETO common flag
}

CLINVAR_HIGH_STARS = {
    # Genome Nexus / ClinVar review-status strings that mean ≥3★.
    "reviewed_by_expert_panel",
    "practice_guideline",
}


def _normalise_review(s: str) -> str:
    return (s or "").strip().lower().replace(" ", "_").replace(",", "")


def compute_flags(row: dict) -> dict:
    """Apply the evidence-flag rule table to a single merged row."""
    flags = {
        # HIGH
        "flag_revue": False,
        "flag_literature_validated": False,
        # MEDIUM
        "flag_cbio_high_recurrence": False,
        "flag_clinvar_plp_validated": False,
        "flag_civic_validated": False,
        # LOW
        "flag_clinvar_plp_unvalidated": False,
        "flag_civic_unvalidated": False,
        "flag_spliceai": False,
        "flag_splicevault": False,
        "flag_gnomad_rare": False,
        # VETO
        "flag_gnomad_common": False,
        "flag_clinvar_benign": False,
    }
    reasons: list[str] = []

    # ----- reVUE (HIGH) -----
    if _bool(row.get("in_revue")):
        flags["flag_revue"] = True
        reasons.append("reVUE-confirmed")

    # ----- Literature validated (HIGH) -----
    validated_with_cancer = (row.get("validated_with_cancer_pmids") or "").strip()
    if validated_with_cancer:
        flags["flag_literature_validated"] = True
        reasons.append(f"Literature: validated_with_cancer PMIDs [{validated_with_cancer}]")

    # ----- cBio recurrence (MEDIUM) -----
    n_total = _to_int(row.get("cbio_total_n_samples"))
    if n_total is not None and n_total >= THRESHOLDS["cbio_high_recurrence"]:
        flags["flag_cbio_high_recurrence"] = True
        reasons.append(f"cBio (MSK+GENIE) n={n_total} ≥ {THRESHOLDS['cbio_high_recurrence']}")

    # ----- ClinVar significance + review status -----
    clinvar_sig = (row.get("clinvar_significance") or "").lower()
    clinvar_conflict = (row.get("clinvar_conflicting") or "").lower()
    is_plp = ("pathogenic" in clinvar_sig) and ("conflict" not in clinvar_conflict and clinvar_conflict == "")
    is_benign = (
        ("benign" in clinvar_sig)
        and "conflict" not in clinvar_conflict
        and clinvar_conflict == ""
    )
    # Review status comes from Genome Nexus's clinvar field. Older runs may not have
    # populated it — gracefully degrade. For now: paper-validation flag drives the
    # high-class gate even when review-status text isn't present.
    review_status = _normalise_review(row.get("clinvar_review_status", ""))
    is_clinvar_high_class = review_status in CLINVAR_HIGH_STARS

    clinvar_paper_validation = (row.get("clinvar_paper_validation") or "").lower()
    clinvar_validated = clinvar_paper_validation == "yes"

    if is_plp:
        # MEDIUM only if high-class AND paper-validated; otherwise LOW.
        if is_clinvar_high_class and clinvar_validated:
            flags["flag_clinvar_plp_validated"] = True
            reasons.append(f"ClinVar P/LP, {review_status}, paper-validated")
        else:
            flags["flag_clinvar_plp_unvalidated"] = True
            extra = []
            if not is_clinvar_high_class:
                extra.append(f"review_status={review_status or '?'}")
            if not clinvar_validated:
                extra.append("no paper validation")
            reasons.append("ClinVar P/LP, " + ", ".join(extra))
    if is_benign:
        flags["flag_clinvar_benign"] = True
        reasons.append("VETO: ClinVar B/LB")

    # ----- CIViC -----
    civic_level = (row.get("civic_best_level") or "").upper()
    civic_sig = (row.get("civic_significance_summary") or "").lower()
    civic_pmids = (row.get("civic_pmids") or "").strip()
    has_civic_oncogenic = bool(civic_sig) and any(
        kw in civic_sig for kw in ("oncogenic", "pathogenic", "positive", "supports oncogenic")
    )
    civic_high_level = civic_level in ("A", "B")
    # CIViC paper validation = at least one CIViC PMID is in validated_with_cancer_pmids.
    if civic_pmids and validated_with_cancer:
        civic_pmid_set = {p.strip() for p in civic_pmids.split(",") if p.strip()}
        validated_set = {p.strip() for p in validated_with_cancer.split(",") if p.strip()}
        civic_paper_validated = bool(civic_pmid_set & validated_set)
    else:
        civic_paper_validated = False

    if has_civic_oncogenic:
        if civic_high_level and civic_paper_validated:
            flags["flag_civic_validated"] = True
            reasons.append(f"CIViC Level {civic_level} oncogenic, paper-validated")
        else:
            flags["flag_civic_unvalidated"] = True
            extras = []
            if not civic_high_level:
                extras.append(f"level={civic_level or '?'}")
            if not civic_paper_validated:
                extras.append("no paper validation")
            reasons.append("CIViC oncogenic, " + ", ".join(extras))

    # ----- SpliceAI (LOW) -----
    sa_max = _to_float(row.get("spliceai_max"))
    if pd.notna(sa_max) and sa_max >= THRESHOLDS["spliceai_strong"]:
        flags["flag_spliceai"] = True
        reasons.append(f"SpliceAI Δ_max = {sa_max:.2f} ≥ {THRESHOLDS['spliceai_strong']}")

    # ----- SpliceVault (LOW) -----
    if _bool(row.get("splicevault_observed")):
        flags["flag_splicevault"] = True
        reasons.append("SpliceVault observed mis-splicing")

    # ----- gnomAD (LOW positive / VETO negative) -----
    af = _to_float(row.get("gnomad_af_popmax"))
    if pd.isna(af):
        # Fall back through other AF columns if popmax isn't computed.
        for k in ("gnomad_joint_af", "gnomad_exome_af", "gnomad_genome_af"):
            af = _to_float(row.get(k))
            if pd.notna(af):
                break

    if pd.isna(af):
        # Variant not in gnomAD → counts as rare.
        flags["flag_gnomad_rare"] = True
        reasons.append("gnomAD: absent → rare")
    elif af > THRESHOLDS["gnomad_common_min"]:
        flags["flag_gnomad_common"] = True
        reasons.append(f"VETO: gnomAD AF {af:.4g} > {THRESHOLDS['gnomad_common_min']}")
    elif af < THRESHOLDS["gnomad_rare_max"]:
        flags["flag_gnomad_rare"] = True
        reasons.append(f"gnomAD rare (AF {af:.2g} < {THRESHOLDS['gnomad_rare_max']})")
    # else: in-between gnomAD frequency — no flag

    return flags, reasons


def assign_tier(flags: dict) -> str:
    """Apply the 3-tier rule from the plan."""
    veto = flags["flag_gnomad_common"] or flags["flag_clinvar_benign"]
    if veto:
        return "Unlikely"
    if flags["flag_revue"] or flags["flag_literature_validated"]:
        return "Likely"
    medium_count = (
        int(flags["flag_cbio_high_recurrence"])
        + int(flags["flag_clinvar_plp_validated"])
        + int(flags["flag_civic_validated"])
    )
    low_count = (
        int(flags["flag_clinvar_plp_unvalidated"])
        + int(flags["flag_civic_unvalidated"])
        + int(flags["flag_spliceai"])
        + int(flags["flag_splicevault"])
        + int(flags["flag_gnomad_rare"])
    )
    if medium_count >= 1 or low_count >= 2:
        return "Maybe"
    return "Unlikely"


def _load_annotated(path: Path) -> pd.DataFrame:
    if not path.exists():
        log.warning("  %s missing — that step's columns will be blank", path.name)
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(UNIQUE_VARIANTS))
    parser.add_argument("--output", default=str(OUTPUT_TSV))
    args = parser.parse_args()

    variants = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    log.info("loaded %d unique variants", len(variants))

    annotated_sources = {
        "genome_nexus": ANNOTATED_DIR / "genome_nexus.tsv",
        "splice_predictions": ANNOTATED_DIR / "splice_predictions.tsv",
        "clinical": ANNOTATED_DIR / "clinical.tsv",
        "paper_context": ANNOTATED_DIR / "paper_context.tsv",
        "recurrence": ANNOTATED_DIR / "recurrence.tsv",
        "population": ANNOTATED_DIR / "population.tsv",
        "literature": ANNOTATED_DIR / "literature.tsv",
    }

    merged = variants.copy()
    merged["start_position"] = merged["start_position"].astype(str)
    merged["end_position"] = merged["end_position"].astype(str)
    for name, path in annotated_sources.items():
        df = _load_annotated(path)
        if df.empty:
            continue
        if "_generated_at" in df.columns:
            df = df.rename(columns={"_generated_at": f"{name}_generated_at"})
        df["start_position"] = df["start_position"].astype(str)
        df["end_position"] = df["end_position"].astype(str)
        merged = merged.merge(df, on=JOIN_KEY, how="left", suffixes=("", f"_{name}"))
    log.info("merged annotations: %d cols", len(merged.columns))

    flag_rows = []
    tiers = []
    reason_strs = []
    for _, r in merged.iterrows():
        flags, reasons = compute_flags(r.to_dict())
        tier = assign_tier(flags)
        flag_rows.append(flags)
        tiers.append(tier)
        reason_strs.append(" | ".join(reasons))

    flags_df = pd.DataFrame(flag_rows)
    out = pd.concat([merged.reset_index(drop=True), flags_df], axis=1)
    out["tier"] = tiers
    out["tier_reasons"] = reason_strs
    out["_generated_at"] = pd.Timestamp.utcnow().isoformat()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    log.info("wrote %d scored variants to %s", len(out), args.output)

    tier_counts = out["tier"].value_counts().to_dict()
    log.info("tier distribution: %s", tier_counts)
    flag_counts = {c: int(flags_df[c].sum()) for c in flags_df.columns}
    log.info("flag counts: %s", flag_counts)
    return 0


if __name__ == "__main__":
    sys.exit(main())
