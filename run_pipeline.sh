#!/usr/bin/env bash
# Pipeline entrypoint. Re-runs every step in order; each step is idempotent
# (skip-on-cache-hit), so it's cheap to re-invoke after new external data
# (cBioPortal TSVs, COSMIC license) lands.
#
# Steps that hit external APIs (2, 3 SpliceVault, 4 CIViC, 6 gnomAD, 7
# LitVar2) read from a per-variant JSON cache in out/cache/. To force a
# fresh API call delete the relevant subdirectory.

set -euo pipefail

PY="${PY:-/Users/lix2/Documents/GitHub/splice-region-mutations/.venv/bin/python}"
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "==> Step 1 — dedup"
$PY src/01_dedup_and_prep.py

echo "==> Step 2 — Genome Nexus"
$PY src/02_annotate_genome_nexus.py

echo "==> Step 3 — Splice predictions (SpliceAI + SpliceVault)"
$PY src/03_annotate_splice_predictions.py --skip-maxent

echo "==> Step 4 — Clinical (CIViC + ClinVar PMIDs + OncoKB CGL)"
$PY src/04_annotate_clinical.py

echo "==> Step 4b — PubMed paper context search"
$PY src/04b_annotate_paper_context.py

echo "==> Step 5 — Recurrence (cBioPortal user-provided TSVs)"
$PY src/05_annotate_recurrence.py

echo "==> Step 6 — Population (gnomAD)"
$PY src/06_annotate_population.py

echo "==> Step 7 — Literature (LitVar2)"
$PY src/07_annotate_literature.py

echo "==> Step 8 — Score & classify"
$PY src/08_score_and_classify.py

echo "==> Step 9 — Build HTML report"
$PY src/09_build_report.py

echo
echo "Done. Open out/report.html in a browser."
