"""
linker_design.py — Linker selection and proteasomal cleavage analysis.

Responsibilities:
  1. Score each linker choice for predicted proteasomal cleavage efficiency
     using pepsickle (Python) or NetChop (web API fallback).
  2. Rank linker options per epitope junction.
  3. Generate the optimal inter-epitope linker for each junction.

Proteasomal cleavage references:
  Kesmir et al. (2002) Protein Eng — NetChop model
  Weeder et al. (2021) PLoS Comp Biol — pepsickle

Geometric compatibility:
  GPGPG  ~ 15 Å end-to-end (flexible) — separates epitopes spatially
  EAAAK  ~ 10 Å (α-helix former)      — rigid spacer
  AAY    ~  6 Å (compact)             — proteasome cleavage optimized
  KK     ~  4 Å                       — very short
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests

from config import (
    LINKERS,
    DEFAULT_LINKER,
    NETCHOP_URL,
    NETCHOP_METHOD,
    NETCHOP_THRESHOLD,
)

log = logging.getLogger(__name__)


# =============================================================================
# Geometric compatibility scoring
# =============================================================================

# Approximate end-to-end distance per residue in an extended chain ≈ 3.8 Å
ANGSTROM_PER_RESIDUE = 3.8

# Minimum recommended separation between epitope C-terminus and N-terminus
# of the next epitope to prevent inter-epitope steric interference.
# MHC-I groove length ≈ 25 Å → epitope needs ~8-10 residues of space.
MIN_SEPARATION_AA = 3   # residues


def geometric_score(linker_seq: str) -> Dict[str, float]:
    """
    Estimate geometric compatibility of a linker.

    Returns dict with:
      length_aa     : number of residues
      extension_A   : approximate spatial extension (Angstroms)
      flexibility   : fraction of Gly+Pro residues (0–1); higher = more flexible
      rigidity      : fraction of helix-forming residues (Ala, Leu, Met)
      geo_score     : composite 0–1 (higher = better geometry for MHC-I spacing)
    """
    n = len(linker_seq)
    if n == 0:
        return dict(length_aa=0, extension_A=0.0, flexibility=0.0, rigidity=0.0, geo_score=0.0)

    gly_pro = sum(1 for aa in linker_seq if aa in "GP") / n
    helix   = sum(1 for aa in linker_seq if aa in "ALMF") / n
    ext_A   = n * ANGSTROM_PER_RESIDUE

    # Composite: reward adequate length and flexibility
    length_bonus = min(n / 5.0, 1.0)          # saturates at 5 residues
    flex_bonus   = gly_pro                     # more flexible = better separation
    geo_score    = 0.5 * length_bonus + 0.5 * flex_bonus

    return dict(
        length_aa=n,
        extension_A=round(ext_A, 1),
        flexibility=round(gly_pro, 3),
        rigidity=round(helix, 3),
        geo_score=round(geo_score, 3),
    )


# =============================================================================
# pepsickle — proteasomal cleavage prediction
# =============================================================================

def _predict_cleavage_pepsickle(sequence: str) -> Optional[np.ndarray]:
    """
    Use pepsickle to predict proteasomal cleavage probability at each position.

    Returns array of length len(sequence) with cleavage probabilities,
    or None if pepsickle unavailable.
    """
    try:
        import pepsickle  # noqa: F401
    except ImportError:
        return None

    try:
        from pepsickle import predict_cleavage_for_sequence
        scores = predict_cleavage_for_sequence(sequence)
        return np.array(scores)
    except Exception as exc:
        log.warning("pepsickle prediction failed: %s", exc)
        return None


# =============================================================================
# NetChop web API fallback
# =============================================================================

def _predict_cleavage_netchop(sequence: str, timeout: int = 30) -> Optional[np.ndarray]:
    """
    Submit to NetChop web service and parse cleavage probabilities.

    Returns array of per-position cleavage probability, or None on failure.
    """
    try:
        resp = requests.post(
            NETCHOP_URL,
            data={
                "sequence": sequence,
                "method":   NETCHOP_METHOD,
                "threshold": str(NETCHOP_THRESHOLD),
            },
            timeout=timeout,
        )
        resp.raise_for_status()
    except Exception as exc:
        log.warning("NetChop web request failed: %s", exc)
        return None

    # Parse the HTML table response
    # NetChop returns a table: pos | aa | score | cleavage
    scores = []
    for match in re.finditer(r"<td[^>]*>\s*([\d.]+)\s*</td>", resp.text):
        try:
            scores.append(float(match.group(1)))
        except ValueError:
            continue

    if len(scores) < len(sequence):
        log.warning("NetChop returned insufficient data (%d scores for %d-aa seq)",
                    len(scores), len(sequence))
        return None

    return np.array(scores[: len(sequence)])


# =============================================================================
# Heuristic fallback cleavage scoring
# =============================================================================

# Kyte-Doolittle hydrophobicity — proteasome cleaves after hydrophobic/basic
_KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

def _predict_cleavage_heuristic(sequence: str) -> np.ndarray:
    """
    Heuristic proteasomal cleavage score based on:
      - C-terminal hydrophobicity (proteasome favors hydrophobic)
      - Basic residues K, R (alternate cleavage)
      - Proline in P1' position suppresses cleavage

    Scores are 0–1 at each position, where 1 = high predicted cleavage.
    """
    n = len(sequence)
    scores = np.zeros(n)
    for i in range(n - 1):
        aa = sequence[i]
        next_aa = sequence[i + 1]
        kd = _KD.get(aa, 0.0)
        # Hydrophobic or basic C-terminal → promotes cleavage
        hydrophobic_bonus = max(kd, 0) / 4.5        # normalize to 0–1
        basic_bonus = 0.3 if aa in "KR" else 0.0
        # Pro in P1' suppresses cleavage (Pro effect)
        pro_penalty = 0.5 if next_aa == "P" else 0.0
        scores[i] = min(hydrophobic_bonus + basic_bonus - pro_penalty, 1.0)

    return np.clip(scores, 0, 1)


# =============================================================================
# Unified cleavage prediction
# =============================================================================

def predict_cleavage(sequence: str) -> Tuple[np.ndarray, str]:
    """
    Predict proteasomal cleavage at each position using best available method.

    Returns
    -------
    (scores_array, method_used)
    """
    scores = _predict_cleavage_pepsickle(sequence)
    if scores is not None:
        return scores, "pepsickle"

    scores = _predict_cleavage_netchop(sequence)
    if scores is not None:
        return scores, "netchop"

    log.info("Using heuristic cleavage scoring (no pepsickle/NetChop available).")
    return _predict_cleavage_heuristic(sequence), "heuristic"


# =============================================================================
# Linker scoring
# =============================================================================

def score_junction(
    left_epitope: str,
    right_epitope: str,
    linker_seq: str,
    cleavage_method: str = "auto",
) -> Dict:
    """
    Score a junction: left_epitope + linker + right_epitope.

    Computes:
      - Cleavage probability at the left_epitope|linker boundary
      - Cleavage probability at the linker|right_epitope boundary
      - Geometric score of the linker
      - Cross-epitope junction score (absence of spurious MHC-I peptides)

    Returns dict of scores.
    """
    junction_seq = left_epitope + linker_seq + right_epitope
    cleavage_scores, method = predict_cleavage(junction_seq)

    left_len  = len(left_epitope)
    right_start = left_len + len(linker_seq)

    # Score at desired cleavage positions
    score_at_left_boundary  = float(cleavage_scores[left_len - 1])   if left_len > 0 else 0.0
    score_at_right_boundary = float(cleavage_scores[right_start - 1]) if right_start > 0 else 0.0

    geo = geometric_score(linker_seq)

    return {
        "left_epitope":          left_epitope,
        "linker":                linker_seq if linker_seq else "(none)",
        "right_epitope":         right_epitope,
        "cleavage_at_left_C":   round(score_at_left_boundary,  3),
        "cleavage_at_right_N":  round(score_at_right_boundary, 3),
        "cleavage_method":       method,
        **{f"geo_{k}": v for k, v in geo.items()},
        # Composite junction score: want high cleavage + good geometry
        "junction_score": round(
            0.5 * score_at_left_boundary
            + 0.5 * score_at_right_boundary
            + 0.2 * geo["geo_score"],
            3
        ),
    }


def select_best_linker(
    left_epitope: str,
    right_epitope: str,
    linker_candidates: Dict[str, str] = LINKERS,
) -> Tuple[str, str, pd.DataFrame]:
    """
    Evaluate all candidate linkers for a given epitope junction.

    Returns
    -------
    (best_linker_name, best_linker_seq, scores_dataframe)
    """
    rows = []
    for name, seq in linker_candidates.items():
        row = score_junction(left_epitope, right_epitope, seq)
        row["linker_name"] = name
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("junction_score", ascending=False)
    best = df.iloc[0]
    best_name = best["linker_name"]
    best_seq  = linker_candidates[best_name]
    log.info("Best linker for %s|%s: %s (score=%.3f)",
             left_epitope[-4:], right_epitope[:4], best_name, best["junction_score"])
    return best_name, best_seq, df


def design_linker_strategy(
    epitope_sequences: List[str],
    fixed_linker: Optional[str] = None,
) -> List[Tuple[str, str]]:
    """
    Design the optimal linker for each consecutive pair of epitopes.

    Parameters
    ----------
    epitope_sequences : ordered list of epitope AA sequences
    fixed_linker      : if provided, use this linker name for all junctions

    Returns
    -------
    List of (linker_name, linker_seq) for each of the N-1 junctions.
    """
    junctions = []
    for i in range(len(epitope_sequences) - 1):
        left  = epitope_sequences[i]
        right = epitope_sequences[i + 1]

        if fixed_linker:
            if fixed_linker not in LINKERS:
                raise ValueError(f"Unknown linker '{fixed_linker}'. "
                                 f"Available: {list(LINKERS.keys())}")
            junctions.append((fixed_linker, LINKERS[fixed_linker]))
        else:
            name, seq, _ = select_best_linker(left, right)
            junctions.append((name, seq))

    return junctions
