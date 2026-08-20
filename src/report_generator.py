"""
report_generator.py — Generates HTML (and optionally PDF) summary reports.

Output structure:
  output/reports/<construct_name>/
    ├── report.html       — main human-readable report
    ├── report.pdf        — PDF version (if weasyprint available)
    ├── mrna_gc.png       — GC content profile plot
    ├── cleavage.png      — proteasomal cleavage plot
    ├── kd_scan.png       — hydrophobicity scan plot
    ├── mhc_heatmap.png   — MHC binding heatmap
    └── mfe_windows.png   — mRNA MFE sliding window plot
"""

from __future__ import annotations

import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

log = logging.getLogger(__name__)

plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 9,
    "axes.titlesize": 10,
})


# =============================================================================
# Plotting functions
# =============================================================================

def plot_gc_profile(gc_data: Dict, output_path: str) -> None:
    """Plot GC content sliding window across the mRNA."""
    windows = gc_data.get("windows", [])
    if not windows:
        return

    centers = [(w["start"] + w["end"]) / 2 for w in windows]
    gc_vals = [w["gc"] for w in windows]
    lo, hi  = gc_data.get("target_range", (0.45, 0.65))

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(centers, gc_vals, color="steelblue", lw=1.2, label="GC content")
    ax.axhline(lo, ls="--", color="orange", lw=1, label=f"Target min ({lo:.0%})")
    ax.axhline(hi, ls="--", color="red",    lw=1, label=f"Target max ({hi:.0%})")
    ax.fill_between(centers, lo, hi, alpha=0.08, color="green", label="Target zone")
    ax.set_xlabel("Position (nt)")
    ax.set_ylabel("GC fraction")
    ax.set_title(f"GC Content Profile — global {gc_data.get('global_gc_pct', 0):.1f}%")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    log.info("Saved GC profile plot: %s", output_path)


def plot_cleavage(cleavage_data: Dict, aa_sequence: str,
                  feature_map: Optional[Dict], output_path: str) -> None:
    """Plot per-residue proteasomal cleavage scores."""
    scores = cleavage_data.get("per_position", [])
    if not scores:
        return

    fig, ax = plt.subplots(figsize=(12, 3))
    positions = list(range(len(scores)))
    ax.plot(positions, scores, color="darkgreen", lw=0.8, alpha=0.8)
    ax.axhline(0.5, ls="--", color="red", lw=1, label="Threshold 0.5")

    # Shade neoantigen regions
    if feature_map:
        colors = plt.cm.Set2.colors
        ci = 0
        for feat, (start, end) in feature_map.items():
            if "neoantigen" in feat.lower():
                ax.axvspan(start, end, alpha=0.15, color=colors[ci % len(colors)],
                           label=feat.split("_", 2)[-1][:20])
                ci += 1

    ax.set_xlabel("Residue position")
    ax.set_ylabel("Cleavage probability")
    ax.set_title(f"Proteasomal Cleavage ({cleavage_data.get('method', 'N/A')})")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    log.info("Saved cleavage plot: %s", output_path)


def plot_kd_hydrophobicity(solubility_data: Dict, output_path: str) -> None:
    """Plot Kyte-Doolittle hydrophobicity sliding window."""
    windows = solubility_data.get("windows", [])
    if not windows:
        return

    centers = [w["center"] for w in windows]
    kd_vals = [w["kd_mean"] for w in windows]
    risks   = [w["risk"]    for w in windows]
    thresh  = solubility_data.get("aggregation_threshold", 1.6)

    fig, ax = plt.subplots(figsize=(10, 3))
    colors = ["tomato" if r else "cornflowerblue" for r in risks]
    ax.bar(centers, kd_vals, width=1, color=colors, alpha=0.8)
    ax.axhline(thresh, ls="--", color="red", lw=1, label=f"Aggregation threshold ({thresh})")
    ax.axhline(0, ls="-", color="gray", lw=0.5)
    ax.set_xlabel("Residue position (window center)")
    ax.set_ylabel("Mean KD hydrophobicity")
    ax.set_title("Kyte-Doolittle Hydrophobicity Scan")
    red_patch  = mpatches.Patch(color="tomato",         label="Aggregation risk")
    blue_patch = mpatches.Patch(color="cornflowerblue", label="Safe")
    ax.legend(handles=[red_patch, blue_patch], fontsize=7)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    log.info("Saved KD plot: %s", output_path)


def plot_mhc_heatmap(epitope_df: pd.DataFrame, output_path: str) -> None:
    """Heatmap of best MHC affinity (nM) per neoantigen × allele."""
    if epitope_df.empty:
        return

    pivot = (
        epitope_df
        .groupby(["neoantigen", "allele"])["affinity_nM"]
        .min()
        .unstack(fill_value=10000)
    )

    fig, ax = plt.subplots(figsize=(max(6, pivot.shape[1] * 1.2),
                                    max(4, pivot.shape[0] * 0.6)))
    sns.heatmap(
        np.log10(pivot + 1),
        ax=ax,
        cmap="RdYlGn_r",
        linewidths=0.5,
        annot=True,
        fmt=".1f",
        cbar_kws={"label": "log10(IC50 nM + 1)"},
    )
    ax.set_title("MHC-I Binding Affinity (log10 IC50, lower = stronger)")
    ax.set_xlabel("HLA Allele")
    ax.set_ylabel("Neoantigen")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    log.info("Saved MHC heatmap: %s", output_path)


def plot_mfe_windows(mfe_data: Dict, output_path: str) -> None:
    """Plot sliding-window MFE across the CDS."""
    data = mfe_data.get("data", [])
    if not data:
        return

    starts = [w["start"] for w in data]
    mfes   = [w["mfe_kcal"] for w in data]

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(starts, mfes, color="purple", lw=1.2)
    ax.fill_between(starts, mfes, 0, alpha=0.15, color="purple")
    ax.axhline(0, ls="-", color="gray", lw=0.5)
    ax.set_xlabel("CDS position (nt)")
    ax.set_ylabel("MFE (kcal/mol)")
    ax.set_title("mRNA CDS Secondary Structure — Sliding Window MFE")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    log.info("Saved MFE windows plot: %s", output_path)


# =============================================================================
# HTML report rendering
# =============================================================================

REPORT_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>mRNA Vaccine Design Report — {name}</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #333; }}
  h1 {{ color: #1a5276; border-bottom: 2px solid #1a5276; }}
  h2 {{ color: #1f618d; margin-top: 30px; }}
  h3 {{ color: #2874a6; }}
  table {{ border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13px; }}
  th {{ background: #1a5276; color: white; padding: 6px 10px; text-align: left; }}
  td {{ padding: 5px 10px; border-bottom: 1px solid #ddd; }}
  tr:nth-child(even) {{ background: #f2f6fc; }}
  .ok   {{ color: #1e8449; font-weight: bold; }}
  .warn {{ color: #e67e22; font-weight: bold; }}
  .fail {{ color: #c0392b; font-weight: bold; }}
  .metric {{ display: inline-block; background: #eaf2ff; border-radius: 6px;
             padding: 8px 16px; margin: 6px; font-size: 14px; }}
  .metric b {{ display: block; font-size: 11px; color: #666; }}
  img {{ max-width: 100%; border: 1px solid #ccc; border-radius: 4px; margin: 8px 0; }}
  pre {{ background: #f4f4f4; padding: 12px; border-radius: 4px; font-size: 11px;
         overflow-x: auto; white-space: pre-wrap; word-break: break-all; }}
  .section {{ margin: 20px 0; padding: 15px; border-left: 4px solid #1a5276;
              background: #f8fbff; }}
</style>
</head>
<body>
<h1>mRNA Neoantigen Vaccine Design Report</h1>
<p><b>Construct:</b> {name} &nbsp;|&nbsp; <b>Generated:</b> {timestamp} &nbsp;|&nbsp;
   <b>Pipeline version:</b> 1.0.0</p>

<div class="section">
  <h2>Design Summary</h2>
  <div class="metric">{total_length} nt<b>Total mRNA length</b></div>
  <div class="metric">{gc_pct:.1f}%<b>GC Content</b></div>
  <div class="metric">{cai:.3f}<b>CAI (codon adaptation)</b></div>
  <div class="metric">{n_epitopes}<b>Neoantigen epitopes</b></div>
  <div class="metric">{fusion_length} aa<b>Fusion protein</b></div>
</div>

<h2>Architecture</h2>
<pre>
5'─[m7G Cap]─[5'UTR {utr5_len} nt]─[Kozak]─ATG─[Neoantigens+Linkers {neo_len} aa]─[HBB {hbb_len} aa]─[3'UTR {utr3_len} nt]─[PolyA {polya_len} nt]─3'
</pre>

{protein_section}

{mrna_section}

{epitope_section}

{stability_section}

{cleavage_section}

{plots_section}

<h2>Full mRNA Sequence (RNA)</h2>
<pre>{mrna_rna}</pre>

<h2>Full Fusion Protein Sequence</h2>
<pre>{fusion_aa}</pre>

<footer>
<hr>
<p style="font-size:11px;color:#888">
Generated by mRNA Neoantigen Vaccine Pipeline v1.0.0 |
References: Kreiter et al. 2015 Nature; Orlandini von Niessen et al. 2019 Nat Commun;
O'Donnell et al. 2020 Cell Systems (mhcflurry); Lorenz et al. 2011 (ViennaRNA)
</p>
</footer>
</body>
</html>"""


def _flag_html(flag_str: str) -> str:
    if "OK" in flag_str:
        return f'<span class="ok">{flag_str}</span>'
    if "WARNING" in flag_str:
        return f'<span class="warn">{flag_str}</span>'
    return f'<span class="fail">{flag_str}</span>'


def _table_from_dict(d: Dict, title: str) -> str:
    rows = "".join(
        f"<tr><td><b>{k}</b></td><td>{v}</td></tr>"
        for k, v in d.items()
        if not isinstance(v, (dict, list))
    )
    return f"<h3>{title}</h3><table><tr><th>Property</th><th>Value</th></tr>{rows}</table>"


def generate_html_report(
    construct_name:  str,
    mrna_construct,          # MRNAConstruct
    fusion_construct,        # FusionConstruct
    epitope_df:      pd.DataFrame,
    protein_report:  Dict,
    stability_report: Dict,
    codon_report:    Dict,
    output_dir:      str,
) -> str:
    """
    Generate a full HTML report and all supporting plots.

    Returns the path to the generated HTML file.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ── Plots ─────────────────────────────────────────────────────────────────
    plots_html = "<h2>Analysis Plots</h2>"

    gc_path = str(out / "mrna_gc.png")
    if stability_report.get("gc_profile"):
        plot_gc_profile(stability_report["gc_profile"], gc_path)
        plots_html += f'<h3>GC Content Profile</h3><img src="mrna_gc.png">'

    cl_path = str(out / "cleavage.png")
    if protein_report.get("proteasomal_cleavage"):
        plot_cleavage(
            protein_report["proteasomal_cleavage"],
            fusion_construct.full_aa,
            fusion_construct.feature_map,
            cl_path,
        )
        plots_html += f'<h3>Proteasomal Cleavage</h3><img src="cleavage.png">'

    kd_path = str(out / "kd_scan.png")
    if protein_report.get("solubility"):
        plot_kd_hydrophobicity(protein_report["solubility"], kd_path)
        plots_html += f'<h3>Hydrophobicity Scan (Kyte-Doolittle)</h3><img src="kd_scan.png">'

    mhc_path = str(out / "mhc_heatmap.png")
    if not epitope_df.empty:
        plot_mhc_heatmap(epitope_df, mhc_path)
        plots_html += f'<h3>MHC-I Binding Heatmap</h3><img src="mhc_heatmap.png">'

    mfe_path = str(out / "mfe_windows.png")
    if stability_report.get("cds_mfe_windows"):
        plot_mfe_windows(stability_report["cds_mfe_windows"], mfe_path)
        plots_html += f'<h3>mRNA MFE Sliding Window (CDS)</h3><img src="mfe_windows.png">'

    # ── Protein section ────────────────────────────────────────────────────────
    pc = protein_report.get("physicochemical", {})
    flags_html = "<br>".join(_flag_html(f) for f in pc.get("flags", []))
    protein_section = f"""
    <div class="section">
    <h2>Protein Quality</h2>
    {_table_from_dict({
        "Length": f"{pc.get('length_aa', '?')} aa",
        "Molecular weight": f"{pc.get('molecular_weight_Da', '?'):.0f} Da",
        "Isoelectric point": f"{pc.get('isoelectric_point', '?'):.2f}",
        "Instability index": f"{pc.get('instability_index', '?'):.1f} ({'stable' if pc.get('is_stable') else 'unstable'})",
        "GRAVY": f"{pc.get('gravy', '?'):.4f} ({'OK' if pc.get('gravy_ok') else 'out of range'})",
        "Aromaticity": f"{pc.get('aromaticity', '?'):.4f}",
    }, "Physicochemical Properties")}
    <p><b>Flags:</b> {flags_html}</p>
    {_table_from_dict({
        "Mean KD hydrophobicity": protein_report.get('solubility', {}).get('mean_kd', '?'),
        "Aggregation risk windows": protein_report.get('solubility', {}).get('n_risk_windows', 0),
        "Net charge (pH 7.4)": protein_report.get('charge', {}).get('net_charge', '?'),
    }, "Solubility & Aggregation")}
    </div>"""

    # ── mRNA section ───────────────────────────────────────────────────────────
    gc_p = stability_report.get("gc_profile", {})
    polya_p = stability_report.get("polya_signal", {})
    mfe_g = stability_report.get("global_mfe", {})
    mrna_section = f"""
    <div class="section">
    <h2>mRNA Properties</h2>
    {_table_from_dict({
        "Total length": f"{mrna_construct.length()} nt",
        "GC content": f"{mrna_construct.gc_content()*100:.1f}%",
        "CAI": codon_report.get('cai', 'N/A'),
        "Global MFE": f"{mfe_g.get('mfe_kcal', 'N/A')} kcal/mol",
        "MFE flag": mfe_g.get('flag', ''),
        "PolyA signal": f"{polya_p.get('signal', 'N/A')} at pos {polya_p.get('position', 'N/A')} — {polya_p.get('status', '')}",
        "Uridine fraction": stability_report.get('uridine', {}).get('U_fraction', 'N/A'),
        "m1Ψ recommendation": stability_report.get('uridine', {}).get('m1psi_recommendation', ''),
    }, "mRNA Metrics")}
    </div>"""

    # ── Epitope section ────────────────────────────────────────────────────────
    if not epitope_df.empty:
        top_rows = (
            epitope_df
            .sort_values("presentation_score", ascending=False)
            .head(20)
        )
        epi_table = top_rows.to_html(index=False, border=0, classes="")
        epitope_section = f"""
        <div class="section">
        <h2>Top MHC-I Epitopes</h2>
        {epi_table}
        </div>"""
    else:
        epitope_section = ""

    # ── Stability section ──────────────────────────────────────────────────────
    utr5_acc = stability_report.get("utr5_accessibility", {})
    stability_section = f"""
    <div class="section">
    <h2>mRNA Stability & Accessibility</h2>
    {_table_from_dict({
        "5'UTR accessibility": utr5_acc.get('accessibility', 'N/A'),
        "5'UTR AUG region flag": utr5_acc.get('flag', ''),
        "GC out-of-range windows": gc_p.get('out_of_range_windows', 0),
        "GC flag": gc_p.get('flag', ''),
    }, "Structural Metrics")}
    </div>"""

    # ── Cleavage section ───────────────────────────────────────────────────────
    ep_scores = protein_report.get("proteasomal_cleavage", {}).get("epitope_scores", {})
    if ep_scores:
        epi_rows = "".join(
            f"<tr><td>{feat.split('_', 2)[-1]}</td>"
            f"<td>{d['start']}–{d['end']}</td>"
            f"<td>{d['mean_cleavage']:.3f}</td>"
            f"<td>{d['c_term_score']:.3f}</td>"
            f"<td>{_flag_html(d['flag'])}</td></tr>"
            for feat, d in ep_scores.items()
        )
        cleavage_section = f"""
        <div class="section">
        <h2>Proteasomal Cleavage per Epitope</h2>
        <table>
        <tr><th>Epitope</th><th>Position</th><th>Mean Score</th><th>C-term Score</th><th>Flag</th></tr>
        {epi_rows}
        </table>
        </div>"""
    else:
        cleavage_section = ""

    # ── Feature map ───────────────────────────────────────────────────────────
    feat_rows = "".join(
        f"<tr><td>{feat}</td><td>{s}</td><td>{e}</td>"
        f"<td>{fusion_construct.full_aa[s:min(s+15,e)]}…</td></tr>"
        for feat, (s, e) in fusion_construct.feature_map.items()
    )

    # ── Assemble full HTML ─────────────────────────────────────────────────────
    html = REPORT_HTML_TEMPLATE.format(
        name             = construct_name,
        timestamp        = datetime.now().strftime("%Y-%m-%d %H:%M"),
        total_length     = mrna_construct.length(),
        gc_pct           = mrna_construct.gc_content() * 100,
        cai              = codon_report.get("cai", 0),
        n_epitopes       = len(fusion_construct.epitope_entries),
        fusion_length    = len(fusion_construct.full_aa),
        utr5_len         = len(mrna_construct.five_utr),
        utr3_len         = len(mrna_construct.three_utr),
        polya_len        = len(mrna_construct.poly_a),
        neo_len          = sum(len(e.sequence) for e in fusion_construct.epitope_entries),
        hbb_len          = len(fusion_construct.hbb_aa),
        protein_section  = protein_section,
        mrna_section     = mrna_section,
        epitope_section  = epitope_section,
        stability_section = stability_section,
        cleavage_section = cleavage_section,
        plots_section    = plots_html,
        mrna_rna         = "\n".join(
            mrna_construct.mrna_rna[i: i+80]
            for i in range(0, len(mrna_construct.mrna_rna), 80)
        ),
        fusion_aa        = "\n".join(
            fusion_construct.full_aa[i: i+60]
            for i in range(0, len(fusion_construct.full_aa), 60)
        ),
    )

    html_path = str(out / "report.html")
    with open(html_path, "w") as fh:
        fh.write(html)
    log.info("HTML report saved: %s", html_path)

    # ── PDF (optional) ─────────────────────────────────────────────────────────
    pdf_path = str(out / "report.pdf")
    try:
        from weasyprint import HTML
        HTML(filename=html_path).write_pdf(pdf_path)
        log.info("PDF report saved: %s", pdf_path)
    except Exception as exc:
        log.info("PDF generation skipped (%s). HTML report available.", exc)

    return html_path
