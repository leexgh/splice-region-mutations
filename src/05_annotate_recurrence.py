"""Step 5 — cBioPortal recurrence annotation.

The user produces the cBioPortal queries externally and drops two TSV files in
`data/cbio/`. This script consumes them and joins on the 5-column variant_key.

If neither file exists, every variant gets blank/null recurrence columns and a
`cbio_status="not_provided"` marker — downstream scoring will skip the
recurrence axis for those.

Expected schemas — see plan Step 5b:
  - `data/cbio/cbio_variant_summary.tsv` (one row per unique variant)
  - `data/cbio/cbio_variant_cancer_types.tsv` (long: variant × cancer_type)
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _apiclient import variant_key  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
UNIQUE_VARIANTS = PROJECT_ROOT / "out" / "unique_variants.tsv"
SUMMARY_TSV = PROJECT_ROOT / "data" / "cbio" / "cbio_variant_summary.tsv"
CANCER_TYPES_TSV = PROJECT_ROOT / "data" / "cbio" / "cbio_variant_cancer_types.tsv"
OUTPUT_TSV = PROJECT_ROOT / "out" / "annotated" / "recurrence.tsv"

log = logging.getLogger("step5")
logging.basicConfig(
    format="[%(name)s] %(message)s", level=logging.INFO, stream=sys.stderr
)

JOIN_KEY_COLS = [
    "chromosome",
    "start_position",
    "end_position",
    "reference_allele",
    "variant_allele",
]

SUMMARY_FIELDS = [
    "cbio_msk_n_samples",
    "cbio_genie_n_samples",
    "cbio_tcga_n_samples",
    "cbio_total_n_samples",
    "splice_rank_in_gene",
    "splice_rank_total",
    "gene_truncating_median_n",
    "cbio_studies",
    "splice_top5_in_gene",
    "avg_mrna_zscore",
    "sample_ids",
]


def empty_recurrence_row() -> dict:
    return {f: "" for f in SUMMARY_FIELDS} | {
        "cbio_cancer_types_json": "",
        "cbio_status": "not_provided",
    }


def _normalise_key_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Cast join columns to string + int matching the variants TSV's types."""
    df = df.copy()
    df["chromosome"] = df["chromosome"].astype(str).str.replace(r"^chr", "", regex=True)
    for col in ("start_position", "end_position"):
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    for col in ("reference_allele", "variant_allele"):
        df[col] = df[col].astype(str).str.upper().str.strip()
    return df


def _build_variant_index(
    variants: pd.DataFrame,
) -> dict[tuple, tuple]:
    """For each input variant, return a dict mapping each candidate join key
    (both VCF-style and MSK-style notations) to its canonical row id.

    Returns: {alt_key_tuple: canonical_5tuple} so that whichever notation the
    user picked in their cBio TSV, we can resolve it back to the canonical row.
    """
    alt_to_canonical: dict[tuple, tuple] = {}
    for _, r in variants.iterrows():
        canonical = (
            r["chromosome"],
            int(r["start_position"]),
            int(r["end_position"]),
            r["reference_allele"].upper(),
            r["variant_allele"].upper(),
        )
        alt_to_canonical[canonical] = canonical
        # MSK-style fallback (only differs for indels).
        msk_ref = (r.get("msk_reference_allele") or "").upper()
        msk_alt = (r.get("msk_variant_allele") or "").upper()
        if msk_ref and msk_alt and (msk_ref == "-" or msk_alt == "-"):
            try:
                msk_key = (
                    r["chromosome"],
                    int(r["msk_start_position"]),
                    int(r["msk_end_position"]),
                    msk_ref,
                    msk_alt,
                )
                alt_to_canonical[msk_key] = canonical
            except (KeyError, ValueError):
                pass
    return alt_to_canonical


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(UNIQUE_VARIANTS))
    parser.add_argument("--summary", default=str(SUMMARY_TSV))
    parser.add_argument("--cancer-types", default=str(CANCER_TYPES_TSV))
    parser.add_argument("--output", default=str(OUTPUT_TSV))
    args = parser.parse_args()

    variants = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    variants["start_position"] = pd.to_numeric(variants["start_position"])
    variants["end_position"] = pd.to_numeric(variants["end_position"])
    log.info("loaded %d unique variants from %s", len(variants), args.input)

    summary_present = Path(args.summary).exists()
    cancer_present = Path(args.cancer_types).exists()
    log.info("summary TSV present: %s", summary_present)
    log.info("cancer-type TSV present: %s", cancer_present)

    if not summary_present and not cancer_present:
        log.warning(
            "No cBioPortal data provided. Emitting empty recurrence annotation. "
            "When user drops the TSVs into %s, re-run this step.",
            Path(args.summary).parent,
        )
        rows = []
        for _, r in variants.iterrows():
            row = {col: r[col] for col in JOIN_KEY_COLS}
            row.update(empty_recurrence_row())
            rows.append(row)
        df = pd.DataFrame(rows)
        df["_generated_at"] = pd.Timestamp.utcnow().isoformat()
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, sep="\t", index=False)
        log.info("wrote %d empty-recurrence rows to %s", len(df), args.output)
        return 0

    # Build alt-key index so we can accept either VCF-style or MSK-style notation
    # in the user-provided TSVs (cBioPortal returns MSK-style for indels).
    alt_to_canonical = _build_variant_index(variants)

    def _resolve(k_in: tuple) -> tuple | None:
        return alt_to_canonical.get(k_in)

    # --- Summary table ---
    summary_idx: dict[tuple, dict] = {}
    if summary_present:
        s = pd.read_csv(args.summary, sep="\t", dtype=str, keep_default_na=False)
        s = _normalise_key_columns(s)
        n_unresolved = 0
        for _, r in s.iterrows():
            raw_key = (
                r["chromosome"],
                int(r["start_position"]),
                int(r["end_position"]),
                r["reference_allele"],
                r["variant_allele"],
            )
            canonical = _resolve(raw_key)
            if canonical is None:
                n_unresolved += 1
                continue
            summary_idx[canonical] = {f: r.get(f, "") for f in SUMMARY_FIELDS if f in r}
        log.info(
            "indexed %d rows from summary TSV (%d rows could not be matched to any input variant)",
            len(summary_idx), n_unresolved,
        )

    # --- Cancer-types long table -> per-variant JSON-encoded Counter ---
    ct_by_key: dict[tuple, dict[str, int]] = defaultdict(dict)
    if cancer_present:
        c = pd.read_csv(args.cancer_types, sep="\t", dtype=str, keep_default_na=False)
        c = _normalise_key_columns(c)
        n_unresolved = 0
        for _, r in c.iterrows():
            raw_key = (
                r["chromosome"],
                int(r["start_position"]),
                int(r["end_position"]),
                r["reference_allele"],
                r["variant_allele"],
            )
            canonical = _resolve(raw_key)
            if canonical is None:
                n_unresolved += 1
                continue
            try:
                n = int(r.get("n_samples", "0") or 0)
            except ValueError:
                n = 0
            ct_by_key[canonical][r.get("cancer_type", "")] = n
        log.info(
            "indexed cancer-type breakdown for %d variants (%d rows unmatched)",
            len(ct_by_key), n_unresolved,
        )

    rows = []
    n_matched = 0
    for _, r in variants.iterrows():
        k = (
            r["chromosome"],
            int(r["start_position"]),
            int(r["end_position"]),
            r["reference_allele"],
            r["variant_allele"],
        )
        row = {col: r[col] for col in JOIN_KEY_COLS}
        summary = summary_idx.get(k)
        if summary is not None:
            n_matched += 1
            for f in SUMMARY_FIELDS:
                row[f] = summary.get(f, "")
        else:
            for f in SUMMARY_FIELDS:
                row[f] = ""
        row["cbio_cancer_types_json"] = (
            json.dumps(ct_by_key[k]) if k in ct_by_key else ""
        )
        row["cbio_status"] = (
            "matched" if summary is not None or k in ct_by_key else "not_in_cbio_tsv"
        )
        rows.append(row)

    df = pd.DataFrame(rows)
    df["_generated_at"] = pd.Timestamp.utcnow().isoformat()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    log.info(
        "wrote %d rows to %s (%d matched cBio summary, %d cancer-type breakdowns)",
        len(df), args.output, n_matched, len(ct_by_key),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
