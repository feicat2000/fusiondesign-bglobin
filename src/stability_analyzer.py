"""
stability_analyzer.py — mRNA secondary structure and stability analysis.

Methods (in priority order):
  1. ViennaRNA (RNA package) — gold standard; full MFE + ensemble analysis
  2. Fallback: GC-content, CAI, and sliding-window estimates

Key metrics computed:
  - MFE (minimum free energy, kcal/mol) — more negative = more structured
  - Ensemble diversity — spread of structures in the Boltzmann ensemble
  - GC content (overall + windowed)
  - 5' UTR accessibility — ideally unstructured near AUG for ribosome binding
  - PolyA signal (AAUAAA) integrity check
  - Uridine content (relevant for immunogenicity of unmodified mRNA)

References:
  Lorenz et al. (2011) Algorithms Mol Biol — ViennaRNA 2.0
  Kozak (1989) J Cell Biol — 5' UTR secondary structure affects translation
  Sample et al. (2019) Nat Biotechnol — 5'UTR optimization for translation
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

from config import MFE_THRESHOLD_KCAL, GC_CONTENT_TARGET, CAI_TARGET

log = logging.getLogger(__name__)


# =============================================================================
# ViennaRNA interface
# =============================================================================

def _rna_available() -> bool:
    try:
        import RNA  # noqa: F401
        return True
    except ImportError:
        return False


def fold_mfe(rna_sequence: str) -> Tuple[Optional[str], Optional[float]]:
    """
    Compute the MFE secondary structure using ViennaRNA RNAfold.

    Returns (dot-bracket structure, MFE in kcal/mol) or (None, None).
    """
    if not _rna_available():
        log.info("ViennaRNA not available — MFE folding skipped.")
        return None, None
    import RNA
    seq = rna_sequence.upper().replace("T", "U")
    structure, mfe = RNA.fold(seq)
    return structure, mfe


def ensemble_analysis(rna_sequence: str) -> Dict:
    """
    Full ensemble analysis: MFE, partition function, ensemble diversity.

    Returns dict with:
      mfe_structure, mfe_kcal, ensemble_energy_kcal,
      centroid_structure, ensemble_diversity, frequency_mfe
    """
    if not _rna_available():
        return {"error": "ViennaRNA not available"}

    import RNA
    seq = rna_sequence.upper().replace("T", "U")
    md = RNA.md()
    fc = RNA.fold_compound(seq, md)

    mfe_structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    pp = fc.pf()   # returns [structure_str, ensemble_energy_kcal] in ViennaRNA 2.x
    centroid, dist = fc.centroid()
    diversity = fc.mean_bp_distance()

    # pf() returns a list [structure, energy] in newer ViennaRNA versions
    ensemble_energy = float(pp[1]) if isinstance(pp, (list, tuple)) else float(pp)

    return {
        "mfe_structure":        mfe_structure,
        "mfe_kcal":             round(mfe, 2),
        "ensemble_energy_kcal": round(ensemble_energy, 2),
        "centroid_structure":   centroid,
        "ensemble_diversity":   round(diversity, 3),
    }


def sliding_window_mfe(
    rna_sequence: str,
    window: int = 120,
    step: int = 10,
) -> List[Dict]:
    """
    Compute MFE in sliding windows across the mRNA.

    Useful for identifying highly structured regions that may impede
    ribosome scanning or affect mRNA half-life.

    Returns list of {start, end, mfe_kcal, gc_fraction}.
    """
    if not _rna_available():
        log.info("ViennaRNA not available — sliding window skipped.")
        return []

    import RNA
    seq = rna_sequence.upper().replace("T", "U")
    results = []
    for start in range(0, len(seq) - window + 1, step):
        chunk = seq[start: start + window]
        _, mfe = RNA.fold(chunk)
        gc = sum(1 for c in chunk if c in "GC") / len(chunk)
        results.append({
            "start":      start,
            "end":        start + window,
            "mfe_kcal":   round(mfe, 2),
            "gc_fraction": round(gc, 3),
        })
    return results


# =============================================================================
# 5' UTR accessibility
# =============================================================================

def check_utr5_accessibility(
    utr5_rna: str,
    aug_flank: int = 30,
) -> Dict:
    """
    Assess 5' UTR accessibility around the AUG start codon.

    A highly structured 5'UTR (low MFE, many paired bases near AUG)
    impedes 43S ribosome scanning and reduces translation efficiency.

    Parameters
    ----------
    utr5_rna  : 5'UTR sequence (RNA, ending just before AUG)
    aug_flank : number of nucleotides on each side of AUG to analyze
    """
    region = utr5_rna[-aug_flank:] if len(utr5_rna) > aug_flank else utr5_rna
    structure, mfe = fold_mfe(region)
    if mfe is None:
        gc = sum(1 for c in region if c in "GCgc") / len(region)
        # Heuristic: higher GC near AUG → less accessible
        accessibility = max(0.0, 1.0 - gc)
        return {
            "region_length":  len(region),
            "gc_fraction":    round(gc, 3),
            "accessibility":  round(accessibility, 3),
            "method":         "gc_heuristic",
        }

    # Count paired bases near AUG
    paired = structure.count("(") + structure.count(")")
    unpaired = structure.count(".")
    total = len(structure)
    accessibility = unpaired / total if total > 0 else 0.0

    return {
        "region_length":  len(region),
        "structure":      structure,
        "mfe_kcal":       round(mfe, 2),
        "paired_fraction": round(paired / total, 3),
        "accessibility":  round(accessibility, 3),
        "method":         "ViennaRNA",
        "flag":           "OK" if accessibility > 0.5 else "WARNING: AUG region may be occluded",
    }


# =============================================================================
# GC content analysis
# =============================================================================

def gc_content_profile(
    sequence: str,
    window: int = 50,
    step: int = 10,
) -> Dict:
    """
    Compute GC content globally and in a sliding window.

    Returns dict with overall GC and list of windowed GC values.
    """
    seq = sequence.upper()
    global_gc = sum(1 for c in seq if c in "GC") / len(seq)

    windows = []
    for start in range(0, len(seq) - window + 1, step):
        chunk = seq[start: start + window]
        gc = sum(1 for c in chunk if c in "GC") / len(chunk)
        windows.append({"start": start, "end": start + window, "gc": round(gc, 3)})

    lo, hi = GC_CONTENT_TARGET
    out_of_range = [w for w in windows if w["gc"] < lo or w["gc"] > hi]

    return {
        "global_gc":        round(global_gc, 4),
        "global_gc_pct":    round(global_gc * 100, 1),
        "target_range":     GC_CONTENT_TARGET,
        "windows":          windows,
        "out_of_range_windows": len(out_of_range),
        "flag":             "OK" if not out_of_range else
                            f"WARNING: {len(out_of_range)} windows outside GC target.",
    }


# =============================================================================
# Uridine / immunogenicity
# =============================================================================

def uridine_analysis(rna_sequence: str) -> Dict:
    """
    Compute U content and UU/AU run statistics.

    Unmodified uridine is recognized by TLR7/8; high U content raises
    innate immune activation. N1-methylpseudouridine (m1Ψ) substitution
    mitigates this.
    """
    seq = rna_sequence.upper().replace("T", "U")
    n = len(seq)
    u_count  = seq.count("U")
    uu_count = sum(1 for i in range(n - 1) if seq[i] == "U" and seq[i+1] == "U")
    au_count = sum(1 for i in range(n - 1) if seq[i] in "AU" and seq[i+1] in "AU")

    return {
        "length_nt":    n,
        "U_count":      u_count,
        "U_fraction":   round(u_count / n, 3),
        "UU_dinucleotides": uu_count,
        "AU_dinucleotides": au_count,
        "m1psi_recommendation": (
            "Apply N1-methylpseudouridine substitution during IVT"
            if u_count / n > 0.20 else
            "Uridine fraction acceptable; m1Ψ still recommended for clinical use"
        ),
    }


# =============================================================================
# polyadenylation signal check
# =============================================================================

def check_polya_signal(three_utr: str) -> Dict:
    """Check for canonical AATAAA polyadenylation signal in 3'UTR."""
    utr = three_utr.upper().replace("U", "T")
    signal = "AATAAA"
    alt_signals = ["ATTAAA", "AGTAAA", "TATAAA"]

    pos = utr.find(signal)
    if pos >= 0:
        return {"signal": signal, "position": pos, "status": "CANONICAL"}

    for alt in alt_signals:
        pos = utr.find(alt)
        if pos >= 0:
            return {"signal": alt, "position": pos, "status": "ALTERNATIVE — verify efficiency"}

    return {"signal": None, "position": None, "status": "WARNING: No polyadenylation signal found"}


# =============================================================================
# Public API
# =============================================================================

def analyze_mrna_stability(
    mrna_rna:  str,
    utr5_rna:  str,
    utr3_seq:  str,
    cds_start: int,
    cds_end:   int,
    run_full_fold: bool = True,
) -> Dict:
    """
    Run the complete mRNA stability analysis suite.

    Parameters
    ----------
    mrna_rna   : full mRNA RNA sequence
    utr5_rna   : 5' UTR sequence (RNA)
    utr3_seq   : 3' UTR sequence
    cds_start  : 0-based start of CDS in mrna_rna
    cds_end    : 0-based end of CDS in mrna_rna
    run_full_fold : if True, run ensemble_analysis (slow for long sequences)

    Returns
    -------
    Nested dict with all stability metrics.
    """
    report = {}

    # Global MFE
    log.info("Computing MFE for full mRNA (%d nt) …", len(mrna_rna))
    mfe_struct, mfe_val = fold_mfe(mrna_rna)
    report["global_mfe"] = {
        "mfe_kcal":  mfe_val,
        "structure": mfe_struct[:80] + "…" if mfe_struct and len(mfe_struct) > 80 else mfe_struct,
        "flag":      (
            "OK" if mfe_val is None or mfe_val > MFE_THRESHOLD_KCAL
            else f"WARNING: MFE {mfe_val:.1f} kcal/mol — highly structured"
        ),
    }

    # Ensemble analysis (optional, slow)
    if run_full_fold and _rna_available():
        log.info("Running ensemble analysis …")
        report["ensemble"] = ensemble_analysis(mrna_rna)

    # 5' UTR accessibility
    report["utr5_accessibility"] = check_utr5_accessibility(utr5_rna)

    # GC content
    report["gc_profile"] = gc_content_profile(mrna_rna)

    # Uridine / immunogenicity
    report["uridine"] = uridine_analysis(mrna_rna)

    # PolyA signal
    report["polya_signal"] = check_polya_signal(utr3_seq)

    # Sliding window MFE (CDS only)
    cds_rna = mrna_rna[cds_start: cds_end]
    log.info("Sliding window MFE over CDS (%d nt) …", len(cds_rna))
    windows = sliding_window_mfe(cds_rna)
    if windows:
        mfes = [w["mfe_kcal"] for w in windows]
        report["cds_mfe_windows"] = {
            "n_windows": len(windows),
            "mean_mfe":  round(float(np.mean(mfes)), 2),
            "max_mfe":   round(float(np.max(mfes)),  2),
            "min_mfe":   round(float(np.min(mfes)),  2),
            "data":      windows,
        }

    log.info("Stability analysis complete.")
    return report
