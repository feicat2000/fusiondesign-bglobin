"""
protein_analyzer.py — Fusion protein quality analysis.

Metrics computed:
  1. Physicochemical properties (Biopython ProteinAnalysis)
       - Molecular weight, isoelectric point, GRAVY, instability index
       - Amino acid composition, aromaticity
  2. Solubility prediction
       - GRAVY score (Kyte-Doolittle) — proxy for solubility
       - CamSol-like window-based hydrophobicity scan
       - Charge at physiological pH
  3. Aggregation risk (sliding-window Kyte-Doolittle)
  4. Proteasomal cleavage site scoring (per-residue)
       - Uses predict_cleavage() from linker_design (pepsickle → NetChop → heuristic)
  5. Epitope accessibility
       - Score each epitope region's average cleavage probability
       - Flag under-scored regions
  6. Secondary structure propensity (Chou-Fasman via Bio.SeqUtils.ProtParam)

References:
  Kyte & Doolittle (1982) J Mol Biol — hydrophobicity scale
  Boman (2003) J Intern Med — charge/hydrophobicity solubility correlation
  Gasteiger et al. (2005) in Proteomics Protocols — ProtParam
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from config import (
    INSTABILITY_INDEX_CUTOFF,
    GRAVY_TARGET_RANGE,
    AGGREGATION_WINDOW,
    AGGREGATION_THRESHOLD,
    NETCHOP_THRESHOLD,
)
from src.linker_design import predict_cleavage

log = logging.getLogger(__name__)


# =============================================================================
# Physicochemical properties
# =============================================================================

def physicochemical_properties(aa_sequence: str) -> Dict:
    """
    Compute physicochemical properties using Biopython ProtParam.

    Returns dict with all key protein quality metrics.
    """
    seq = aa_sequence.upper()
    # ProteinAnalysis requires standard 20 AAs
    clean_seq = "".join(c for c in seq if c.isalpha() and c not in "BJOUXZ")

    pa = ProteinAnalysis(clean_seq)

    mw        = pa.molecular_weight()
    pi        = pa.isoelectric_point()
    instab    = pa.instability_index()
    gravy     = pa.gravy()
    aromatics = pa.aromaticity()

    secondary = pa.secondary_structure_fraction()   # (helix, turn, sheet)
    aa_comp   = pa.get_amino_acids_percent()

    gravy_lo, gravy_hi = GRAVY_TARGET_RANGE
    return {
        "length_aa":          len(clean_seq),
        "molecular_weight_Da": round(mw, 1),
        "isoelectric_point":   round(pi, 2),
        "instability_index":   round(instab, 2),
        "is_stable":           instab < INSTABILITY_INDEX_CUTOFF,
        "gravy":               round(gravy, 4),
        "gravy_ok":            gravy_lo <= gravy <= gravy_hi,
        "aromaticity":         round(aromatics, 4),
        "ss_helix_fraction":   round(secondary[0], 3),
        "ss_turn_fraction":    round(secondary[1], 3),
        "ss_sheet_fraction":   round(secondary[2], 3),
        "aa_composition":      {k: round(v, 4) for k, v in aa_comp.items()},
        "flags": _physicochemical_flags(instab, gravy, pi),
    }


def _physicochemical_flags(instab: float, gravy: float, pi: float) -> List[str]:
    flags = []
    if instab >= INSTABILITY_INDEX_CUTOFF:
        flags.append(f"INSTABILITY: index {instab:.1f} ≥ {INSTABILITY_INDEX_CUTOFF} (unstable)")
    lo, hi = GRAVY_TARGET_RANGE
    if gravy > hi:
        flags.append(f"HYDROPHOBIC: GRAVY {gravy:.3f} > {hi} (aggregation risk)")
    if gravy < lo:
        flags.append(f"HYDROPHILIC: GRAVY {gravy:.3f} < {lo} (may be too charged)")
    if pi < 5.0 or pi > 9.0:
        flags.append(f"CHARGE: pI {pi:.1f} outside 5–9 (may affect solubility)")
    if not flags:
        flags.append("OK")
    return flags


# =============================================================================
# Solubility prediction
# =============================================================================

# Kyte-Doolittle scale
_KD_SCALE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8,  "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def solubility_scan(
    aa_sequence: str,
    window: int = 9,
    threshold: float = AGGREGATION_THRESHOLD,
) -> Dict:
    """
    Sliding-window Kyte-Doolittle hydrophobicity scan.

    Regions with mean KD > threshold are flagged as aggregation-prone.

    Returns dict with per-window scores and aggregation risk regions.
    """
    seq = aa_sequence.upper()
    half = window // 2
    kd_vals = [_KD_SCALE.get(aa, 0.0) for aa in seq]

    windows = []
    for center in range(half, len(seq) - half):
        chunk = kd_vals[center - half: center + half + 1]
        mean_kd = float(np.mean(chunk))
        windows.append({
            "center": center,
            "aa":     seq[center],
            "kd_mean": round(mean_kd, 3),
            "risk":   mean_kd > threshold,
        })

    risk_regions = [w for w in windows if w["risk"]]
    overall_kd   = float(np.mean(kd_vals))

    return {
        "mean_kd":              round(overall_kd, 3),
        "aggregation_threshold": threshold,
        "n_risk_windows":       len(risk_regions),
        "risk_centers":         [w["center"] for w in risk_regions],
        "windows":              windows,
        "flag": (
            "OK" if not risk_regions else
            f"WARNING: {len(risk_regions)} aggregation-prone windows detected"
        ),
    }


# =============================================================================
# Charge profile
# =============================================================================

def charge_profile(aa_sequence: str, pH: float = 7.4) -> Dict:
    """
    Compute net charge at a given pH and identify charged clusters.

    Uses Biopython ProteinAnalysis for charge calculation.
    """
    pa = ProteinAnalysis(aa_sequence.upper())
    charge = pa.charge_at_pH(pH)
    basic  = sum(1 for aa in aa_sequence if aa in "KRH")
    acidic = sum(1 for aa in aa_sequence if aa in "DE")

    return {
        "pH":          pH,
        "net_charge":  round(charge, 2),
        "basic_aa":    basic,
        "acidic_aa":   acidic,
        "charge_bias":  "basic" if charge > 2 else ("acidic" if charge < -2 else "neutral"),
    }


# =============================================================================
# Proteasomal cleavage accessibility
# =============================================================================

def protease_accessibility(
    aa_sequence: str,
    feature_map: Optional[Dict[str, Tuple[int, int]]] = None,
) -> Dict:
    """
    Compute proteasomal cleavage scores and assess epitope accessibility.

    Parameters
    ----------
    aa_sequence : full fusion protein AA sequence
    feature_map : positions of each neoantigen in the full sequence
                  (from FusionConstruct.feature_map); if provided, each
                  epitope's mean cleavage score is reported.

    Returns
    -------
    Dict with per-position cleavage scores and per-epitope summaries.
    """
    scores, method = predict_cleavage(aa_sequence)
    log.info("Proteasomal cleavage scored via %s.", method)

    epitope_scores = {}
    if feature_map:
        for feature, (start, end) in feature_map.items():
            if "neoantigen" in feature.lower():
                region_scores = scores[start:end]
                epitope_scores[feature] = {
                    "start":       start,
                    "end":         end,
                    "mean_cleavage": round(float(np.mean(region_scores)), 3),
                    "max_cleavage":  round(float(np.max(region_scores)),  3),
                    "c_term_score":  round(float(scores[end - 1]),         3),
                    "flag": (
                        "OK" if scores[end - 1] >= NETCHOP_THRESHOLD else
                        "WARNING: low C-terminal cleavage score"
                    ),
                }

    # Find top-cleavage positions
    top_positions = np.argsort(scores)[-10:][::-1]
    top_sites = [
        {"position": int(i), "aa": aa_sequence[i], "score": round(float(scores[i]), 3)}
        for i in top_positions
    ]

    return {
        "method":          method,
        "per_position":    [round(float(s), 3) for s in scores],
        "mean_score":      round(float(np.mean(scores)), 3),
        "top_sites":       top_sites,
        "epitope_scores":  epitope_scores,
    }


# =============================================================================
# Cross-epitope junction peptide check
# =============================================================================

def junction_peptide_scan(
    full_aa:     str,
    feature_map: Dict[str, Tuple[int, int]],
    window:      int = 9,
    alleles:     Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Scan peptides spanning epitope-linker and linker-epitope junctions.

    Cross-junction peptides can create spurious MHC-I epitopes that may
    cause unwanted immune responses.

    Returns DataFrame of junction-spanning peptides for manual review.
    """
    junction_rows = []
    neoantigen_features = [
        (name, start, end)
        for name, (start, end) in feature_map.items()
        if "neoantigen" in name.lower()
    ]

    for i in range(len(neoantigen_features) - 1):
        name1, _, end1   = neoantigen_features[i]
        name2, start2, _ = neoantigen_features[i + 1]
        junction_center  = end1   # last AA of epitope1; junction starts here

        # Scan window around the junction
        scan_start = max(0, junction_center - window // 2)
        scan_end   = min(len(full_aa), junction_center + window)

        for j in range(scan_start, scan_end - window + 1):
            pep = full_aa[j: j + window]
            is_junction = j < junction_center < j + window
            junction_rows.append({
                "junction":          f"{name1}|{name2}",
                "peptide":           pep,
                "start":             j,
                "end":               j + window,
                "spans_junction":    is_junction,
            })

    df = pd.DataFrame(junction_rows)
    if df.empty:
        return df

    # Optionally score with mhcflurry
    if alleles:
        try:
            from src.epitope_analysis import predict_mhcflurry
            mhc = predict_mhcflurry(df["peptide"].tolist(), alleles)
            df = df.merge(
                mhc.groupby("peptide")["affinity_nM"].min().reset_index()
                   .rename(columns={"affinity_nM": "best_affinity_nM"}),
                on="peptide", how="left",
            )
        except Exception as exc:
            log.warning("MHC scoring of junction peptides failed: %s", exc)

    return df


# =============================================================================
# Comprehensive protein report
# =============================================================================

def analyze_protein(
    aa_sequence:  str,
    feature_map:  Optional[Dict] = None,
    hla_alleles:  Optional[List[str]] = None,
) -> Dict:
    """
    Run full protein quality analysis suite.

    Returns nested dict with all metrics.
    """
    report = {}

    log.info("Computing physicochemical properties …")
    report["physicochemical"] = physicochemical_properties(aa_sequence)

    log.info("Running solubility scan …")
    report["solubility"] = solubility_scan(aa_sequence)

    log.info("Computing charge profile …")
    report["charge"] = charge_profile(aa_sequence)

    log.info("Predicting proteasomal cleavage …")
    report["proteasomal_cleavage"] = protease_accessibility(aa_sequence, feature_map)

    if feature_map:
        log.info("Scanning junction peptides …")
        junction_df = junction_peptide_scan(aa_sequence, feature_map, alleles=hla_alleles)
        if not junction_df.empty:
            spanning = junction_df[junction_df["spans_junction"]]
            report["junction_peptides"] = {
                "total":               len(junction_df),
                "spanning_junctions":  len(spanning),
                "dataframe":           junction_df,   # stored for report_generator
            }

    log.info("Protein analysis complete.")
    return report
