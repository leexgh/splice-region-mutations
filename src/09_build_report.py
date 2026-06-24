"""Step 9 — Build the HTML report (3-tier evidence-checkbox layout).

Reads `out/variants_scored.tsv` and emits a single self-contained
`out/report.html` with:
- Header coverage stats (no OncoKB comparison).
- Master sortable/searchable table.
- Three tier sections: Likely / Maybe / Unlikely (oncogenic).
- Per-variant cards with evidence checkboxes grouped by weight (HIGH/MEDIUM/LOW + VETO).
- cBioPortal sample linkouts (MSK, GENIE, TCGA) with show-more.
- IGV.js track per card zoomed to variant location.
- Manual-review panel for variants where ClinVar P/LP or CIViC was paper-validated.
"""

from __future__ import annotations

import argparse
import html
import json
import logging
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCORED_TSV = PROJECT_ROOT / "out" / "variants_scored.tsv"
REPORT_HTML = PROJECT_ROOT / "out" / "report.html"

log = logging.getLogger("step9")
logging.basicConfig(
    format="[%(name)s] %(message)s", level=logging.INFO, stream=sys.stderr
)

TIER_ORDER = ["Likely", "Maybe", "Unlikely"]
TIER_COLORS = {
    "Likely": "#c62828",
    "Maybe": "#ef6c00",
    "Unlikely": "#546e7a",
}

FLAG_LABELS = [
    ("flag_revue", "HIGH", "reVUE confirmed (curated splice-altering variant)"),
    ("flag_literature_validated", "HIGH", "Literature support"),
    ("flag_cbio_high_recurrence", "MEDIUM", "cBioPortal (MSK+GENIE+TCGA) high recurrence (≥5 samples)"),
    ("flag_clinvar_plp_validated", "MEDIUM", "ClinVar P/LP at ≥3★, paper-validated"),
    ("flag_civic_validated", "MEDIUM", "CIViC Level A/B oncogenic, paper-validated"),
    ("flag_clinvar_plp_unvalidated", "LOW", "ClinVar P/LP (low review status or no paper validation)"),
    ("flag_civic_unvalidated", "LOW", "CIViC oncogenic, no paper validation"),
    ("flag_spliceai", "LOW", "SpliceAI Δ_max ≥ 0.5"),
    ("flag_splicevault", "LOW", "SpliceVault observed mis-splicing at this junction"),
    ("flag_gnomad_rare", "LOW", "gnomAD rare (popmax < 1e-4 or absent)"),
    ("flag_gnomad_common", "VETO", "gnomAD common (popmax > 1e-3) — VETO"),
    ("flag_clinvar_benign", "VETO", "ClinVar Benign / Likely Benign — VETO"),
]

WEIGHT_COLORS = {
    "HIGH": "#c62828",
    "MEDIUM": "#ef6c00",
    "LOW": "#9e9e9e",
    "VETO": "#1565c0",
}


def _h(s) -> str:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return ""
    return html.escape(str(s))


def _bool(x) -> bool:
    return str(x).strip().lower() in ("true", "1")


def _format_num(x, prec=2) -> str:
    try:
        if x in (None, "", "nan"):
            return ""
        f = float(x)
        return f"{f:.{prec}f}"
    except (TypeError, ValueError):
        return _h(x)


def _pmid_links(pmid_str: str) -> str:
    if not pmid_str or pmid_str == "nan":
        return "<span class='muted'>—</span>"
    out = []
    for p in str(pmid_str).split(","):
        p = p.strip()
        if not p.isdigit():
            continue
        out.append(f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}/" target="_blank">{p}</a>')
    return ", ".join(out) if out else "<span class='muted'>—</span>"


def _clinvar_link(cvid) -> str:
    if not cvid or str(cvid) in ("", "nan"):
        return "<span class='muted'>—</span>"
    return f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{cvid}/" target="_blank">{cvid}</a>'


def _format_cancer_types(json_str: str) -> str:
    if not json_str or json_str == "nan":
        return "<span class='muted'>—</span>"
    try:
        d = json.loads(json_str)
    except (json.JSONDecodeError, TypeError):
        return _h(json_str)
    items = sorted(d.items(), key=lambda kv: -kv[1])
    return ", ".join(f"{_h(k)} ({v})" for k, v in items)


def _card_id(r: dict) -> str:
    chrom = str(r.get('chromosome', '')).replace(' ', '')
    pos   = str(r.get('start_position', '')).replace(' ', '')
    ref   = re.sub(r'[^A-Za-z0-9]', 'x', str(r.get('reference_allele', '')))
    alt   = re.sub(r'[^A-Za-z0-9]', 'x', str(r.get('variant_allele', '')))
    return f"card-{chrom}-{pos}-{ref}-{alt}"


def _cbio_sample_links(sample_ids_str: str, uid: str, max_shown: int = 10) -> str:
    """Generate cBioPortal hyperlinks for each sample_id, with show-more if >max_shown."""
    if not sample_ids_str or str(sample_ids_str) in ("", "nan"):
        return "<span class='muted'>—</span>"

    samples = [s.strip() for s in str(sample_ids_str).split(",") if s.strip()]
    links = []

    for s in samples:
        # TCGA format: TCGA-XX-YYYY-01(study_id)
        tcga_m = re.match(r'^(TCGA-[^(]+)\(([^)]+)\)$', s)
        if tcga_m:
            sid = tcga_m.group(1).strip()
            study = tcga_m.group(2).strip()
            url = f"https://cbioportal.mskcc.org/patient?sampleId={sid}&studyId={study}"
            links.append(f'<a href="{url}" target="_blank" class="lnk-tcga">{_h(sid)}</a>')
        elif s.upper().startswith("GENIE-"):
            url = f"https://genie.cbioportal.org/patient?sampleId={s}&studyId=genie_public"
            links.append(f'<a href="{url}" target="_blank" class="lnk-genie">{_h(s)}</a>')
        else:
            # MSK sample
            url = f"https://cbioportal.mskcc.org/patient?sampleId={s}&studyId=mskimpact"
            links.append(f'<a href="{url}" target="_blank" class="lnk-msk">{_h(s)}</a>')

    if not links:
        return "<span class='muted'>—</span>"

    if len(links) <= max_shown:
        return " ".join(links)

    shown = links[:max_shown]
    hidden = links[max_shown:]
    more_id = f"sm-{uid}"
    return (
        " ".join(shown) +
        f'<span id="{more_id}" style="display:none;"> ' + " ".join(hidden) + "</span>"
        f' <button class="more-btn" '
        f'onclick="var e=document.getElementById(\'{more_id}\');'
        f'e.style.display=e.style.display===\'none\'?\'inline\':\'\';'
        f'this.textContent=this.textContent===\'+{len(hidden)} more\'?\'show less\':\'+{len(hidden)} more\';">'
        f'+{len(hidden)} more</button>'
    )


def build_evidence_checklist(r: dict) -> str:
    sections: dict[str, list[str]] = {"HIGH": [], "MEDIUM": [], "LOW": [], "VETO": []}
    for col, weight, label in FLAG_LABELS:
        checked = _bool(r.get(col))
        mark = "checkbox-checked" if checked else "checkbox-empty"
        symbol = "✓" if checked else "&nbsp;"
        sections[weight].append(
            f'<li class="{mark}"><span class="check">{symbol}</span> {_h(label)}</li>'
        )

    out = []
    for weight in ("HIGH", "MEDIUM", "LOW", "VETO"):
        color = WEIGHT_COLORS[weight]
        items = "".join(sections[weight])
        out.append(
            f'<div class="flag-group" style="border-left:4px solid {color};">'
            f'<h5 style="color:{color};">{weight} weight</h5>'
            f'<ul>{items}</ul></div>'
        )
    return "\n".join(out)


def build_manual_review(r: dict) -> str:
    raw = r.get("manual_review_papers_json") or ""
    if not raw or raw == "nan":
        return ""
    try:
        papers = json.loads(raw)
    except (json.JSONDecodeError, TypeError):
        return ""
    if not papers:
        return ""
    cards = []
    for p in papers:
        pmid = _h(p.get("pmid", ""))
        source = _h(p.get("source", ""))
        pattern = _h(p.get("pattern", ""))
        snippet = _h(p.get("snippet", ""))
        cards.append(
            f'<div class="manual-review-card">'
            f'<div class="mr-head"><a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">PMID {pmid}</a> · {source} · matched <code>{pattern}</code></div>'
            f'<div class="mr-snippet">… {snippet} …</div>'
            f'</div>'
        )
    return (
        '<div class="manual-review">'
        '<h5>📚 Manual-review panel — paper-validated ClinVar / CIViC evidence</h5>'
        + "".join(cards)
        + "</div>"
    )


def build_card(r: dict) -> str:
    tier = r.get("tier", "")
    color = TIER_COLORS.get(tier, "#666")
    title = (
        f"{_h(r.get('hugo_symbol'))} · "
        f"{_h(r.get('gn_hgvsc'))} · "
        f"{_h(r.get('gn_hgvsp_short'))} · "
        f"chr{_h(r.get('chromosome'))}:{_h(r.get('start_position'))} "
        f"{_h(r.get('reference_allele'))}&gt;{_h(r.get('variant_allele'))}"
    )

    cancers_input = _format_cancer_types(r.get("cancer_types_json", ""))
    cancers_cbio = _format_cancer_types(r.get("cbio_cancer_types_json", ""))

    splice_html = f"""
    <tr><td><b>SpliceAI Δ_max</b></td><td>{_format_num(r.get('spliceai_max'), 3)}
        <span class='muted'>(AG {_format_num(r.get('spliceai_ds_ag'), 2)} /
        AL {_format_num(r.get('spliceai_ds_al'), 2)} /
        DG {_format_num(r.get('spliceai_ds_dg'), 2)} /
        DL {_format_num(r.get('spliceai_ds_dl'), 2)})</span></td></tr>
    <tr><td>SpliceVault</td><td>{'observed mis-splicing' if _bool(r.get('splicevault_observed')) else _h(r.get('splicevault_status') or '—')}</td></tr>
    <tr><td>MaxEntScan Δ</td><td>{_format_num(r.get('maxent_delta'), 2)} ({_h(r.get('maxent_status', ''))})</td></tr>
    """

    clinical_html = f"""
    <tr><td><b>ClinVar</b></td><td>{_clinvar_link(r.get('clinvar_id'))}: {_h(r.get('clinvar_significance')) or '—'}</td></tr>
    <tr><td>ClinVar review status</td><td>{_h(r.get('clinvar_review_status') or '—')}</td></tr>
    <tr><td>ClinVar PMIDs</td><td>{_pmid_links(r.get('clinvar_pmids', ''))}</td></tr>
    <tr><td>ClinVar paper validation</td><td>{_h(r.get('clinvar_paper_validation') or '—')}</td></tr>
    <tr><td>CIViC</td><td>level {_h(r.get('civic_best_level') or '—')}, {_h(r.get('civic_evidence_count', 0))} evidence items</td></tr>
    <tr><td>CIViC PMIDs</td><td>{_pmid_links(r.get('civic_pmids', ''))}</td></tr>
    """

    n_total = r.get("cbio_total_n_samples") or ""
    n_msk   = r.get("cbio_msk_n_samples") or "—"
    n_genie = r.get("cbio_genie_n_samples") or "—"
    n_tcga  = r.get("cbio_tcga_n_samples") or "—"
    zscore  = _format_num(r.get("avg_mrna_zscore"), 3)
    cid     = _card_id(r)
    sample_links = _cbio_sample_links(r.get("sample_ids", ""), cid)

    rec_html = f"""
    <tr><td><b>cBio status</b></td><td>{_h(r.get('cbio_status') or 'not_provided')}</td></tr>
    <tr><td>Total samples (MSK+GENIE+TCGA deduped)</td><td>{_h(n_total)} (MSK {_h(n_msk)} / GENIE {_h(n_genie)} / TCGA {_h(n_tcga)})</td></tr>
    <tr><td>Splice rank in gene</td><td>{_h(r.get('splice_rank_in_gene') or '—')} of {_h(r.get('splice_rank_total') or '—')}</td></tr>
    <tr><td>TSG truncating median</td><td>{_h(r.get('gene_truncating_median_n') or '—')}</td></tr>
    <tr><td>Top 5 in gene</td><td>{_h(r.get('splice_top5_in_gene') or '—')}</td></tr>
    <tr><td>Cancer types (cBio)</td><td>{cancers_cbio}</td></tr>
    <tr><td>Studies</td><td>{_h(r.get('cbio_studies') or '—')}</td></tr>
    <tr><td>mRNA z-score (TCGA avg)</td><td>{zscore if zscore else '<span class="muted">—</span>'}</td></tr>
    <tr><td>Sample links</td><td class="sample-links">{sample_links}</td></tr>
    """

    pop_html = f"""
    <tr><td><b>gnomAD AF popmax</b></td><td>{_format_num(r.get('gnomad_af_popmax'), 5) or '<span class="muted">—</span>'}</td></tr>
    <tr><td>gnomAD joint AF</td><td>{_format_num(r.get('gnomad_joint_af'), 5) or '<span class="muted">—</span>'}</td></tr>
    <tr><td>gnomAD dataset</td><td>{_h(r.get('gnomad_dataset') or '—')}</td></tr>
    <tr><td>Filters</td><td>{_h(r.get('gnomad_filters') or '—')}</td></tr>
    """

    lit_html = f"""
    <tr><td><b>PMIDs (validated w/ cancer)</b></td><td>{_pmid_links(r.get('validated_with_cancer_pmids', ''))}</td></tr>
    <tr><td>PMIDs (variant only)</td><td>{_pmid_links(r.get('validated_variant_only_pmids', ''))}</td></tr>
    <tr><td>PMIDs (unvalidated)</td><td>{_pmid_links(r.get('unvalidated_pmids', ''))}</td></tr>
    <tr><td>reVUE PMIDs</td><td>{_pmid_links(r.get('revue_pmids', ''))}</td></tr>
    """

    cohort_html = f"""
    <tr><td><b>Input occurrences</b></td><td>{_h(r.get('occurrence_count'))} ({_h(r.get('n_distinct_samples'))} samples)</td></tr>
    <tr><td>Cancer types (input MAF)</td><td>{cancers_input}</td></tr>
    <tr><td>Gene</td><td>{_h(r.get('hugo_symbol'))} · OncoKB type {_h(r.get('gene_type_oncokb') or '—')} <span class='muted'>(TSG/oncogene flag only)</span></td></tr>
    <tr><td>reVUE</td><td>{'<b>confirmed</b>' if _bool(r.get('in_revue')) else '—'}</td></tr>
    """

    reasons = _h(r.get("tier_reasons", "")).replace(" | ", "<br>")
    checklist_html = build_evidence_checklist(r)
    manual_review_html = build_manual_review(r)

    chrom = str(r.get('chromosome', '')).replace('chr', '')
    pos   = str(r.get('start_position', '0'))
    end   = str(r.get('end_position', pos))
    gene  = _h(r.get('hugo_symbol', ''))
    hgvsc = _h(r.get('gn_hgvsc', ''))

    return f"""
<details class="card" id="{cid}" data-tier="{tier}"
  data-chrom="{chrom}" data-pos="{pos}" data-end="{end}"
  data-gene="{gene}" data-hgvsc="{hgvsc}">
<summary>
  <span class="tier-badge" style="background:{color};">{tier}</span>
  <span class="title">{title}</span>
</summary>
<div class="card-body">
  <div class="checklist">{checklist_html}</div>
  {manual_review_html}
  <div class="grid">
    <section><h4>Splice prediction</h4><table>{splice_html}</table></section>
    <section><h4>Clinical</h4><table>{clinical_html}</table></section>
    <section><h4>Recurrence (cBio MSK+GENIE+TCGA)</h4><table>{rec_html}</table></section>
    <section><h4>Population (gnomAD)</h4><table>{pop_html}</table></section>
    <section><h4>Literature</h4><table>{lit_html}</table></section>
    <section><h4>Cohort / input</h4><table>{cohort_html}</table></section>
  </div>
  <div class="reasons"><b>Reasons:</b><br>{reasons}</div>
  <div class="igv-wrap">
    <h4>Genomic view (IGV) — <span class="igv-locus">chr{chrom}:{pos}</span></h4>
    <div class="igv-container" id="igv-{cid}"></div>
  </div>
</div>
</details>
"""


def _ratio(occ, total) -> str:
    """Return formatted ratio or 'NA' if denominator is 0 or missing."""
    try:
        o = float(occ)
        t = float(total)
        if t == 0:
            return "NA"
        return f"{o / t:.3f}"
    except (TypeError, ValueError):
        return "NA"


def build_master_row(r: dict) -> str:
    tier = r.get("tier", "")
    color = TIER_COLORS.get(tier, "#666")
    cid = _card_id(r)
    ratio = _ratio(r.get('occurrence_count', ''), r.get('cbio_total_n_samples', ''))
    oncogenic = _h(r.get('oncokb_call_input', ''))
    level     = _h(r.get('oncokb_highest_level_input', ''))
    zscore    = _format_num(r.get('avg_mrna_zscore'), 2)
    return f"""
<tr data-card="{cid}" title="Click to open evidence card">
  <td><span class="tier-badge" style="background:{color};">{_h(tier)}</span></td>
  <td>{_h(r.get('hugo_symbol'))}</td>
  <td>{_h(r.get('gn_hgvsp_short'))}</td>
  <td>{_h(r.get('gn_hgvsc'))}</td>
  <td>chr{_h(r.get('chromosome'))}:{_h(r.get('start_position'))} {_h(r.get('reference_allele'))}&gt;{_h(r.get('variant_allele'))}</td>
  <td>{_format_num(r.get('spliceai_max'), 2)}</td>
  <td>{_format_num(r.get('gnomad_af_popmax'), 5)}</td>
  <td>{_h(r.get('clinvar_significance'))}</td>
  <td>{_h(r.get('cbio_total_n_samples') or '—')}</td>
  <td>{'✓' if _bool(r.get('in_revue')) else ''}</td>
  <td>{_h(r.get('occurrence_count'))}</td>
  <td data-order="{ratio if ratio != 'NA' else -1}">{ratio}</td>
  <td>{oncogenic}</td>
  <td>{level}</td>
  <td data-order="{zscore if zscore else ''}">{zscore if zscore else '—'}</td>
</tr>
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(SCORED_TSV))
    parser.add_argument("--output", default=str(REPORT_HTML))
    args = parser.parse_args()

    if not Path(args.input).exists():
        log.error("Scored TSV not found at %s — run Step 8 first.", args.input)
        return 1

    df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    log.info("loaded %d scored variants", len(df))

    df["_high_count"] = sum(df.get(c, "").map(_bool) for c, w, _ in FLAG_LABELS if w == "HIGH")
    df["_medium_count"] = sum(df.get(c, "").map(_bool) for c, w, _ in FLAG_LABELS if w == "MEDIUM")
    df["_low_count"] = sum(df.get(c, "").map(_bool) for c, w, _ in FLAG_LABELS if w == "LOW")
    df = df.sort_values(
        ["tier", "_high_count", "_medium_count", "_low_count"],
        ascending=[True, False, False, False],
    ).reset_index(drop=True)

    tier_counts = Counter(df["tier"])

    n_total = len(df)
    n_sa = (df.get("spliceai_max", pd.Series(dtype=str)) != "").sum()
    n_sa_strong = (pd.to_numeric(df.get("spliceai_max"), errors="coerce") >= 0.5).sum()
    n_clinvar = (df.get("clinvar_id", pd.Series(dtype=str)) != "").sum()
    n_clinvar_pmids = (df.get("clinvar_pmids", pd.Series(dtype=str)) != "").sum()
    n_validated = (df.get("validated_with_cancer_pmids", pd.Series(dtype=str)) != "").sum()
    n_manual_review = (df.get("manual_review_papers_json", pd.Series(dtype=str)) != "").sum()
    n_cbio = (df.get("cbio_status", pd.Series(dtype=str)) == "matched").sum()
    n_gnomad = (df.get("gnomad_status", pd.Series(dtype=str)) == "ok").sum()
    n_revue = df["in_revue"].apply(_bool).sum()

    tier_summary_rows = "".join(
        f'<tr><td><span class="tier-badge" style="background:{TIER_COLORS.get(t,"#666")};">{t}</span></td>'
        f'<td>{tier_counts.get(t,0)}</td></tr>'
        for t in TIER_ORDER
    )

    coverage_rows = f"""
      <tr><td>SpliceAI scored</td><td>{n_sa}/{n_total} ({n_sa_strong} with Δ≥0.5)</td></tr>
      <tr><td>ClinVar IDs</td><td>{n_clinvar}/{n_total}</td></tr>
      <tr><td>ClinVar PMIDs</td><td>{n_clinvar_pmids}/{n_total}</td></tr>
      <tr><td>Literature-validated</td><td>{n_validated}/{n_total} variants have ≥1 PMID with direct genomic + cancer mention</td></tr>
      <tr><td>Manual-review panels</td><td>{n_manual_review}/{n_total} variants surfaced for ClinVar/CIViC paper review</td></tr>
      <tr><td>cBioPortal matched</td><td>{n_cbio}/{n_total}</td></tr>
      <tr><td>gnomAD ok</td><td>{n_gnomad}/{n_total}</td></tr>
      <tr><td>reVUE confirmed</td><td>{n_revue}/{n_total}</td></tr>
    """

    sections_html = []
    for tier in TIER_ORDER:
        sub = df[df["tier"] == tier]
        if sub.empty:
            continue
        cards = "\n".join(build_card(r.to_dict()) for _, r in sub.iterrows())
        sections_html.append(
            f'<section class="tier-section">'
            f'<h2 style="color:{TIER_COLORS.get(tier,"#666")};">{tier} (oncogenic): {len(sub)} variants</h2>'
            f'{cards}</section>'
        )
    sections = "\n".join(sections_html)

    master_rows = "\n".join(build_master_row(r.to_dict()) for _, r in df.iterrows())

    html_out = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Splice-region oncogenicity report</title>
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 0; padding: 0; color: #222; }}
header {{ background: #f7f7f7; border-bottom: 1px solid #ddd; padding: 16px 24px; }}
h1 {{ margin: 0; font-size: 22px; }}
.summary-grid {{ display: grid; grid-template-columns: 1fr 2fr; gap: 24px; margin-top: 16px; }}
.summary-grid > div {{ background: white; padding: 12px; border: 1px solid #e0e0e0; border-radius: 4px; }}
.summary-grid h3 {{ margin-top: 0; font-size: 14px; color: #555; }}
table {{ border-collapse: collapse; }}
.summary-grid table td {{ padding: 3px 8px; font-size: 14px; }}
.tier-badge {{ display: inline-block; padding: 2px 8px; border-radius: 3px; color: white; font-size: 11px; font-weight: 600; }}
main {{ padding: 16px 24px 64px; }}
.tier-section {{ margin-top: 24px; }}
.tier-section h2 {{ font-size: 18px; margin-bottom: 8px; border-bottom: 2px solid currentColor; padding-bottom: 4px; }}
.card {{ background: #fafafa; border: 1px solid #e0e0e0; border-radius: 4px; margin-bottom: 6px; }}
.card summary {{ cursor: pointer; padding: 8px 12px; display: flex; gap: 12px; align-items: center; font-family: monospace; font-size: 13px; }}
.card .title {{ flex: 1; }}
.card-body {{ padding: 12px; border-top: 1px solid #e0e0e0; background: white; }}
.checklist {{ display: grid; grid-template-columns: 1fr 1fr 1fr 1fr; gap: 8px; margin-bottom: 12px; }}
.flag-group {{ padding: 8px 12px; background: #fafafa; border-radius: 4px; }}
.flag-group h5 {{ margin: 4px 0 6px; font-size: 12px; font-weight: 600; }}
.flag-group ul {{ list-style: none; margin: 0; padding: 0; font-size: 12px; }}
.flag-group li {{ padding: 2px 0; }}
.flag-group li.checkbox-checked {{ font-weight: 600; }}
.flag-group li.checkbox-empty .check {{ color: #ccc; }}
.flag-group li .check {{ font-family: monospace; display: inline-block; width: 14px; text-align: center; border: 1px solid #aaa; border-radius: 2px; margin-right: 4px; }}
.flag-group li.checkbox-checked .check {{ background: #e8f5e9; color: #2e7d32; border-color: #2e7d32; }}
.manual-review {{ margin: 12px 0; padding: 10px 14px; background: #fff8e1; border-left: 4px solid #ffa000; border-radius: 4px; }}
.manual-review h5 {{ margin: 0 0 6px; color: #e65100; font-size: 13px; }}
.manual-review-card {{ margin: 6px 0; padding: 6px 10px; background: white; border: 1px solid #ffe082; border-radius: 4px; font-size: 12px; }}
.mr-head {{ margin-bottom: 4px; font-weight: 600; color: #5d4037; }}
.mr-head code {{ background: #ffecb3; padding: 1px 4px; border-radius: 2px; }}
.mr-snippet {{ color: #555; font-style: italic; line-height: 1.4; }}
.grid {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 12px; }}
.grid section {{ background: #fafafa; padding: 8px 12px; border-radius: 4px; }}
.grid h4 {{ margin: 4px 0 6px; font-size: 13px; color: #555; }}
.grid table {{ width: 100%; font-size: 12px; }}
.grid td {{ padding: 2px 6px; vertical-align: top; border-bottom: 1px dashed #eee; }}
.grid td:first-child {{ color: #555; width: 40%; }}
.reasons {{ margin-top: 12px; padding: 8px 12px; background: #fffaef; border-left: 3px solid #ef6c00; font-size: 12px; }}
.muted {{ color: #999; }}
#master-table {{ width: 100%; font-size: 12px; }}
#master-table th {{ background: #f7f7f7; padding: 6px 8px; text-align: left; border-bottom: 1px solid #ccc; }}
#master-table td {{ padding: 4px 8px; border-bottom: 1px solid #eee; }}
#master-section {{ margin-bottom: 24px; }}
.note {{ font-size: 12px; color: #777; margin-top: 4px; }}
th[title] {{ cursor: help; border-bottom: 1px dashed #aaa; }}
#master-table tbody tr {{ cursor: pointer; }}
#master-table tbody tr:hover {{ background: #fffde7 !important; }}
.card-flash {{ outline: 3px solid #ef6c00; outline-offset: 2px; background: #fff8e1 !important; transition: outline 0.3s; }}
.sample-links a {{ font-size: 11px; padding: 1px 4px; border-radius: 2px; margin: 1px; display: inline-block; text-decoration: none; }}
.lnk-msk  {{ background: #e3f2fd; color: #1565c0; }}
.lnk-genie {{ background: #e8f5e9; color: #2e7d32; }}
.lnk-tcga  {{ background: #fce4ec; color: #880e4f; }}
.more-btn {{ font-size: 11px; padding: 1px 6px; cursor: pointer; background: #eee; border: 1px solid #ccc; border-radius: 3px; margin-left: 4px; }}
.igv-wrap {{ margin-top: 16px; border-top: 1px solid #eee; padding-top: 12px; }}
.igv-wrap h4 {{ margin: 0 0 8px; font-size: 13px; color: #555; }}
.igv-locus {{ font-family: monospace; color: #1565c0; }}
.igv-container {{ width: 100%; height: 300px; border: 1px solid #ddd; border-radius: 4px; overflow: hidden; background: #fff; }}
</style>
</head>
<body>

<header>
<h1>Splice-region oncogenicity report — independent re-classification</h1>
<div class="note">
{n_total} unique variants. Sources: SpliceAI (local) + SpliceVault + ClinVar + CIViC + reVUE + cBioPortal (MSK + GENIE + TCGA) + gnomAD + PubMed paper context search.
OncoKB oncogenicity calls and therapy levels shown from input data only (not independently verified in this pipeline).
</div>

<div class="summary-grid">
  <div>
    <h3>Tier distribution</h3>
    <table>{tier_summary_rows}</table>
  </div>
  <div>
    <h3>Pipeline coverage</h3>
    <table>{coverage_rows}</table>
  </div>
</div>
</header>

<main>

<section id="master-section">
<h2>Master table (sortable, searchable)</h2>
<table id="master-table" class="display">
<thead>
<tr>
  <th>Tier</th><th>Gene</th><th>HGVSp_short</th><th>HGVSc</th><th>Position</th>
  <th>SpliceAI Δmax</th><th>gnomAD popmax</th><th>ClinVar</th>
  <th title="Unique samples carrying this variant across MSK-IMPACT + GENIE + TCGA. MSK↔GENIE deduped by patient ID; TCGA patients are added separately. ≥5 triggers the MEDIUM recurrence flag.">cBio N (MSK+GENIE+TCGA)</th>
  <th>reVUE</th>
  <th title="Occurrences in the input MSK-IMPACT MAF before dedup. One patient can contribute multiple rows from repeat biopsies or timepoints.">Occ count</th>
  <th title="Occ count divided by cBio total unique samples (MSK+GENIE+TCGA). NA if cBio N = 0.">Occ/Total ratio</th>
  <th title="OncoKB oncogenicity annotation from input data.">OncoKB Oncogenic</th>
  <th title="OncoKB highest therapeutic level from input data.">OncoKB Level</th>
  <th title="Average mRNA expression z-score across TCGA samples carrying this variant (relative to diploid; negative = reduced expression).">mRNA Z-score (TCGA)</th>
</tr>
</thead>
<tbody>
{master_rows}
</tbody>
</table>
</section>

<h2>Per-tier evidence cards</h2>
{sections}

</main>

<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
<script>
$(document).ready(function() {{
  $('#master-table').DataTable({{
    pageLength: 25,
    order: [[0, 'asc'], [8, 'desc']],
    columnDefs: [
      {{ type: 'num', targets: [5, 6, 8, 10, 14] }},
      {{ type: 'num', targets: [11] }},
      {{ orderable: true, targets: '_all' }}
    ]
  }});

  $('#master-table tbody').on('click', 'tr', function() {{
    var cardId = $(this).data('card');
    if (!cardId) return;
    var card = document.getElementById(cardId);
    if (!card) return;
    card.open = true;
    card.scrollIntoView({{ behavior: 'smooth', block: 'start' }});
    card.classList.add('card-flash');
    setTimeout(function() {{ card.classList.remove('card-flash'); }}, 1800);
  }});
}});

// IGV.js lazy initialization — load browser only when card is first expanded.
var igvBrowsers = {{}};

document.querySelectorAll('details.card').forEach(function(card) {{
  card.addEventListener('toggle', function() {{
    if (!this.open) return;
    var cid = this.id;
    if (igvBrowsers[cid]) return; // already initialised

    var chrom = this.dataset.chrom;
    var pos   = parseInt(this.dataset.pos, 10);
    var end   = parseInt(this.dataset.end || pos, 10);
    var gene  = this.dataset.gene || '';
    var hgvsc = this.dataset.hgvsc || '';

    var container = document.getElementById('igv-' + cid);
    if (!container) return;

    var locus = 'chr' + chrom + ':' + Math.max(1, pos - 300) + '-' + (end + 300);
    var options = {{
      genome: 'hg19',
      locus: locus,
      tracks: [
        {{
          type: 'annotation',
          name: gene + (hgvsc ? ' ' + hgvsc : ''),
          features: [{{
            chr: 'chr' + chrom,
            start: pos - 1,
            end: end,
            name: gene + (hgvsc ? ' ' + hgvsc : '')
          }}],
          color: '#c62828',
          displayMode: 'EXPANDED',
          order: 1
        }}
      ]
    }};

    igv.createBrowser(container, options).then(function(browser) {{
      igvBrowsers[cid] = browser;
    }}).catch(function(err) {{
      container.innerHTML = '<div style="padding:8px;color:#999;font-size:12px;">IGV failed to load: ' + err + '</div>';
    }});
  }});
}});
</script>
</body>
</html>
"""

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(html_out)
    log.info("wrote report to %s (%.1f KB)", args.output, Path(args.output).stat().st_size / 1024)
    log.info("tier distribution: %s", dict(tier_counts))
    return 0


if __name__ == "__main__":
    sys.exit(main())
