# Splice-Region Mutations — Oncogenicity Re-classification

Independent evidence-based re-classification of MSK-IMPACT splice-region variants into three tiers (**Likely** / **Maybe** / **Unlikely** oncogenic), without relying on OncoKB per-variant calls.

**Live report → https://leexgh.github.io/splice-region-mutations/**

---

## Overview

MSK-IMPACT frequently flags splice-region variants as "Likely Oncogenic" via OncoKB, but this class is mechanistically heterogeneous and those calls are often unvalidated. This pipeline independently scores each variant using converging evidence from multiple sources and applies a transparent rule-based tier assignment.

### Evidence sources

| Source | What it contributes |
|--------|---------------------|
| **SpliceAI** (local) | Δ_max splice-disruption score |
| **SpliceVault** | Observed mis-splicing events at the junction |
| **ClinVar** | P/LP classification + review star rating + PMIDs |
| **CIViC** | Clinical oncogenic evidence (Level A/B) |
| **reVUE** | Gold-standard curated splice-altering variants |
| **cBioPortal** | Sample recurrence across MSK-IMPACT + GENIE + 32 TCGA studies |
| **gnomAD** | Population allele frequency (rare/common filter) |
| **PubMed** | Abstract-level validation of ClinVar/CIViC PMIDs |

### Tier rule

```
if any VETO flag (gnomAD common OR ClinVar Benign)  → Unlikely
elif any HIGH flag (reVUE OR literature-validated)   → Likely
elif any MEDIUM flag OR ≥2 LOW flags                 → Maybe
else                                                  → Unlikely
```

### Current results (927 unique variants)

| Tier | Count |
|------|------:|
| Likely | 26 |
| Maybe | 362 |
| Unlikely | 539 |

---

## Installation

### Requirements

- Python 3.11 (arm64 for Apple Silicon; required for TensorFlow/SpliceAI)
- ~4 GB disk for the hg19 reference FASTA (optional — only needed to re-run Step 1 for indel normalisation)

### 1. Clone the repo

```bash
git clone https://github.com/leexgh/splice-region-mutations.git
cd splice-region-mutations
```

### 2. Create a virtual environment

```bash
python3.11 -m venv .venv
source .venv/bin/activate
```

### 3. Install dependencies

```bash
pip install --upgrade pip
pip install pandas requests numpy==1.26.4
pip install tensorflow==2.21.0        # for SpliceAI (arm64 only)
pip install spliceai pyfaidx
```

> **Note:** `tensorflow` and `spliceai` are only needed for Step 3 (SpliceAI scoring). All other steps run on standard `pandas` + `requests`.

### 4. Provide input data

Place these files in `data/`:

| File | Description |
|------|-------------|
| `data/New_Query_*.csv` | MSK-IMPACT MAF-style export (occurrence-level, with `ONCOGENIC` and `HIGHEST_LEVEL` columns) |
| `data/VUEs.txt` | reVUE export (tab-separated, `genomicLocation` column) |
| `data/refs/hg19.fa` + `.fai` | GRCh37 reference FASTA (for indel anchor-base normalisation) |
| `data/cbio/cbio_variant_summary.tsv` | cBioPortal recurrence summary (one row per variant; see schema below) |
| `data/cbio/cbio_variant_cancer_types.tsv` | cBioPortal cancer-type breakdown (long format) |

#### cBioPortal TSV schema

The `cbio_variant_summary.tsv` must have these columns:

```
chromosome  start_position  end_position  reference_allele  variant_allele
hugo_symbol  hgvsp_short
cbio_msk_n_samples  cbio_genie_n_samples  cbio_tcga_n_samples  cbio_total_n_samples
splice_rank_in_gene  splice_rank_total  gene_truncating_median_n
cbio_studies  splice_top5_in_gene  avg_mrna_zscore  sample_ids
```

`sample_ids` format: MSK samples as `P-XXXXXXX-T01-IM*`, GENIE as `GENIE-CENTER-*`, TCGA as `TCGA-XX-YYYY-01(study_id)`.

---

## Running the pipeline

```bash
bash run_pipeline.sh
```

All annotation steps are **idempotent** — API responses are cached in `out/cache/{api}/`. Re-running is cheap; only changed steps need recomputing.

To force a fresh API call, delete the relevant cache directory:

```bash
rm -rf out/cache/gnomad/   # re-fetch gnomAD
```

### Running individual steps

```bash
source .venv/bin/activate

python src/01_dedup_and_prep.py          # dedup input → out/unique_variants.tsv
python src/05_annotate_recurrence.py     # join cBio TSVs → out/annotated/recurrence.tsv
python src/08_score_and_classify.py      # apply tier rules → out/variants_scored.tsv
python src/09_build_report.py            # build HTML report → out/report.html
```

Steps 2–7 (Genome Nexus, SpliceAI, ClinVar, CIViC, gnomAD, LitVar2) hit external APIs and are cached. After new cBioPortal data arrives, only steps 5 → 8 → 9 need to re-run.

### Pipeline architecture

```
data/
  New_Query_*.csv          ← input MAF
  VUEs.txt                 ← reVUE gold standard
  cbio/                    ← user-provided cBioPortal TSVs
  refs/hg19.fa             ← GRCh37 reference

src/
  01_dedup_and_prep.py     → out/unique_variants.tsv
  02_annotate_genome_nexus.py
  03_annotate_splice_predictions.py   (SpliceAI + SpliceVault)
  04_annotate_clinical.py             (CIViC + ClinVar + OncoKB gene list)
  04b_annotate_paper_context.py       (PubMed abstract search)
  05_annotate_recurrence.py           (cBioPortal TSVs)
  06_annotate_population.py           (gnomAD)
  07_annotate_literature.py           (reVUE + LitVar2)
  08_score_and_classify.py            → out/variants_scored.tsv
  09_build_report.py                  → out/report.html

out/
  variants_scored.tsv      ← full audit TSV (every flag + tier)
  report.html              ← self-contained HTML deliverable
  cache/{api}/*.json       ← per-variant API response cache
  annotated/*.tsv          ← per-step intermediate files
```

---

## Viewing the report locally

The report uses [IGV.js](https://github.com/igvteam/igv.js/) to render a genomic view per variant card. IGV requires HTTP (not `file://`) to load the hg19 reference from Broad's CDN.

```bash
python serve.py          # starts http://localhost:8080/report.html
```

Or specify a different port:

```bash
python serve.py 9090
```

---

## Report features

- **Master table** — sortable and searchable across all 927 variants. Columns: Tier, Gene, HGVSp, HGVSc, Position, SpliceAI Δmax, gnomAD popmax, ClinVar, cBio N, reVUE, Occ count, Occ/Total ratio, OncoKB Oncogenic (from input), OncoKB Level (from input), mRNA Z-score (TCGA avg).
- **Per-variant evidence cards** — evidence checkboxes grouped by weight (HIGH / MEDIUM / LOW / VETO), Recurrence section with MSK / GENIE / TCGA breakdown, cBioPortal sample hyperlinks with show-more, IGV.js track zoomed to the variant locus.
- **Manual-review panel** — surfaces ClinVar P/LP and CIViC variants where a PubMed abstract matched both a variant-specific notation and a cancer keyword.

---

## Updating the live report

After re-running the pipeline:

```bash
cp out/report.html docs/index.html
git add docs/index.html
git commit -m "Update report"
git push
```

GitHub Pages rebuilds automatically (~1–2 min).
