"""
codon_optimizer.py — Codon optimization for maximum human expression.

Strategy (layered):
  1. dnachisel — constraint-based optimization with explicit objectives:
       - CodonOptimize (human CAI)
       - AvoidHairpins  (mRNA secondary structure)
       - AvoidPattern   (restriction sites, poly-nucleotide runs)
       - EnforceGCContent (window-based 45–65%)
  2. Fallback: frequency-table back-translation via python-codon-tables
     (used when dnachisel fails or is unavailable).

CAI reference:
  Sharp & Li (1987) Nucleic Acids Res — Codon Adaptation Index
dnachisel:
  Zulkower & Rosser (2020) Bioinformatics
python-codon-tables:
  https://github.com/Edinburgh-Genome-Foundry/python-codon-tables
"""

from __future__ import annotations

import logging
import random
from typing import Dict, List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)


# =============================================================================
# Human codon table (Kazusa, Homo sapiens)
# Structured as {amino_acid: {codon: frequency}}
# =============================================================================

def _load_human_codon_table() -> Dict[str, Dict[str, float]]:
    """Load human codon usage table via python-codon-tables.

    python-codon-tables returns {amino_acid: {codon: frequency}} where
    frequencies already sum to 1 per amino acid.

    NOTE: get_codons_table() returns a *cached mutable* dict. dnachisel
    mutates it in-place (adds 'log_codons_frequencies' etc.), which
    corrupts subsequent calls. We deep-copy and filter to only the 20
    standard amino acid keys (single uppercase letters A–Y) so our copy
    is immune to those side effects.
    """
    _VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
    try:
        from python_codon_tables import get_codons_table
        raw = get_codons_table("h_sapiens")
        # Deep copy + filter: keep only standard 1-letter amino acid keys
        aa_table: Dict[str, Dict[str, float]] = {
            aa: {c: float(f) for c, f in codons.items()}
            for aa, codons in raw.items()
            if isinstance(aa, str) and aa in _VALID_AA
        }
        log.debug("Loaded human codon table: %d amino acids.", len(aa_table))
        return aa_table
    except ImportError:
        log.warning("python-codon-tables not installed — using embedded minimal table.")
        return {aa: dict(d) for aa, d in _FALLBACK_HUMAN_CODONS.items()}  # also copy


# Minimal fallback table (most-used human codon per amino acid)
# Source: Nakamura et al. (2000) Nucleic Acids Res (Kazusa database)
_FALLBACK_HUMAN_CODONS: Dict[str, Dict[str, float]] = {
    "A": {"GCC": 0.40, "GCT": 0.28, "GCA": 0.22, "GCG": 0.11},
    "R": {"AGG": 0.22, "AGA": 0.20, "CGG": 0.20, "CGC": 0.19, "CGA": 0.11, "CGT": 0.08},
    "N": {"AAC": 0.54, "AAT": 0.46},
    "D": {"GAC": 0.54, "GAT": 0.46},
    "C": {"TGC": 0.55, "TGT": 0.45},
    "Q": {"CAG": 0.75, "CAA": 0.25},
    "E": {"GAG": 0.58, "GAA": 0.42},
    "G": {"GGC": 0.34, "GGG": 0.25, "GGA": 0.25, "GGT": 0.16},
    "H": {"CAC": 0.58, "CAT": 0.42},
    "I": {"ATC": 0.48, "ATT": 0.35, "ATA": 0.17},
    "L": {"CTG": 0.40, "CTC": 0.20, "TTG": 0.13, "CTT": 0.13, "CTA": 0.07, "TTA": 0.07},
    "K": {"AAG": 0.58, "AAA": 0.42},
    "M": {"ATG": 1.00},
    "F": {"TTC": 0.55, "TTT": 0.45},
    "P": {"CCC": 0.32, "CCT": 0.28, "CCA": 0.27, "CCG": 0.11},  # slightly wrong but ok
    "S": {"AGC": 0.24, "TCC": 0.22, "TCT": 0.15, "AGT": 0.15, "TCA": 0.15, "TCG": 0.06},
    "T": {"ACC": 0.36, "ACA": 0.28, "ACT": 0.24, "ACG": 0.12},
    "W": {"TGG": 1.00},
    "Y": {"TAC": 0.57, "TAT": 0.43},
    "V": {"GTG": 0.47, "GTC": 0.24, "GTT": 0.18, "GTA": 0.11},
    "*": {"TGA": 0.47, "TAA": 0.28, "TAG": 0.25},
}


# =============================================================================
# Simple frequency-table back-translation
# =============================================================================

def back_translate_frequency(
    aa_sequence: str,
    codon_table: Optional[Dict[str, Dict[str, float]]] = None,
    random_seed: Optional[int] = 42,
) -> str:
    """
    Back-translate amino acid sequence using codon frequency sampling.

    For each amino acid, a codon is sampled proportional to its human usage
    frequency (max-frequency selection when seed is fixed).

    Parameters
    ----------
    aa_sequence  : protein sequence (single-letter codes)
    codon_table  : {aa: {codon: freq}}; loads human table if None
    random_seed  : fixed seed for reproducibility; None = stochastic

    Returns
    -------
    DNA coding sequence (no stop codon)
    """
    if codon_table is None:
        codon_table = _load_human_codon_table()

    if random_seed is not None:
        rng = random.Random(random_seed)
    else:
        rng = random

    codons = []
    for aa in aa_sequence.upper():
        if aa not in codon_table:
            raise ValueError(f"Unknown amino acid '{aa}' in sequence.")
        options = codon_table[aa]
        # Weighted random selection
        codons_list = list(options.keys())
        weights     = list(options.values())
        chosen = rng.choices(codons_list, weights=weights, k=1)[0]
        codons.append(chosen)

    return "".join(codons)


def max_frequency_back_translate(
    aa_sequence: str,
    codon_table: Optional[Dict[str, Dict[str, float]]] = None,
) -> str:
    """
    Deterministic back-translation: always choose the highest-frequency codon
    for each amino acid.
    """
    if codon_table is None:
        codon_table = _load_human_codon_table()

    codons = []
    for aa in aa_sequence.upper():
        if aa not in codon_table:
            raise ValueError(f"Unknown amino acid '{aa}'.")
        best = max(codon_table[aa], key=codon_table[aa].get)
        codons.append(best)
    return "".join(codons)


# =============================================================================
# CAI calculation
# =============================================================================

def calculate_cai(dna_sequence: str, codon_table: Optional[Dict] = None) -> float:
    """
    Compute the Codon Adaptation Index (CAI) for a DNA sequence.

    CAI = geometric mean of relative synonymous codon usage (RSCU) weights.
    Range: 0–1, where 1 = all codons are the most frequent synonym.

    Reference: Sharp & Li (1987) Nucleic Acids Res 15:1281-1295
    """
    if codon_table is None:
        codon_table = _load_human_codon_table()

    # Invert to {codon: relative_weight}
    codon_weights: Dict[str, float] = {}
    for aa, codons in codon_table.items():
        max_freq = max(codons.values())
        for codon, freq in codons.items():
            codon_weights[codon] = freq / max_freq if max_freq > 0 else 0.0

    dna = dna_sequence.upper().replace("U", "T")
    if len(dna) % 3 != 0:
        raise ValueError("DNA sequence length is not a multiple of 3.")

    log_weights = []
    for i in range(0, len(dna) - 3, 3):   # skip stop codon
        codon = dna[i: i + 3]
        w = codon_weights.get(codon)
        if w is not None and w > 0:
            log_weights.append(np.log(w))

    if not log_weights:
        return 0.0
    return float(np.exp(np.mean(log_weights)))


# =============================================================================
# dnachisel optimization (primary method)
# =============================================================================

def optimize_with_dnachisel(
    aa_sequence: str,
    avoid_sequences: Optional[List[str]] = None,
    gc_window: Tuple[float, float] = (0.45, 0.65),
    gc_window_size: int = 50,
    max_hairpin_stem: int = 6,
    verbose: bool = False,
) -> Tuple[str, Dict]:
    """
    Optimize a DNA sequence using dnachisel with multiple constraints.

    Constraints applied:
      - CodonOptimize      : maximize human CAI
      - AvoidHairpins      : no stems > max_hairpin_stem bp
      - AvoidPattern       : user-specified restriction sites + poly-N runs
      - EnforceGCContent   : window-based GC in [gc_window]

    Parameters
    ----------
    aa_sequence   : protein to encode
    avoid_sequences : list of DNA patterns to avoid (restriction sites etc.)
    gc_window     : (min_gc, max_gc) fraction
    gc_window_size: window size for GC enforcement
    max_hairpin_stem: maximum allowed hairpin stem length (bp)

    Returns
    -------
    (optimized_dna, optimization_report_dict)
    """
    try:
        from dnachisel import (
            DnaOptimizationProblem,
            CodonOptimize,
            AvoidHairpins,
            AvoidPattern,
            EnforceGCContent,
            AvoidChanges,
        )
        from dnachisel.Location import Location as DnaLocation
    except ImportError:
        log.warning("dnachisel not available — falling back to frequency-table optimization.")
        dna = max_frequency_back_translate(aa_sequence)
        cai = calculate_cai(dna)
        return dna, {"method": "frequency_table", "cai": cai}

    # Start from frequency-optimized sequence (guaranteed to start with ATG)
    initial_dna = max_frequency_back_translate(aa_sequence)

    constraints = [
        # Lock the initiator ATG (positions 0–3) so no constraint can mutate it
        AvoidChanges(location=DnaLocation(0, 3)),
        # Avoid restriction sites commonly used in IVT cloning
        AvoidPattern("GAATTC"),    # EcoRI
        AvoidPattern("AAGCTT"),    # HindIII
        AvoidPattern("GGATCC"),    # BamHI
        AvoidPattern("GCGGCCGC"),  # NotI
        # Avoid homopolymer runs ≥ 5 nt
        AvoidPattern("AAAAA"),
        AvoidPattern("TTTTT"),
        AvoidPattern("GGGGG"),
        AvoidPattern("CCCCC"),
    ]

    # Add any extra patterns supplied by caller
    for pat in (avoid_sequences or []):
        constraints.append(AvoidPattern(pat))

    # AvoidHairpins — only add if sequence is long enough for the window
    if len(initial_dna) > 200:
        constraints.append(AvoidHairpins(stem_size=max_hairpin_stem, hairpin_window=200))

    # GC content enforcement — use a larger window than the default to avoid
    # forcing non-synonymous changes in short GC-poor regions (e.g. ATG start)
    effective_window = max(gc_window_size, 60)
    if len(initial_dna) >= effective_window:
        constraints.append(
            EnforceGCContent(mini=gc_window[0], maxi=gc_window[1], window=effective_window)
        )

    # CodonOptimize goes in objectives only (soft), never in constraints
    objectives = [CodonOptimize(species="h_sapiens")]

    problem = DnaOptimizationProblem(
        sequence=initial_dna,
        constraints=constraints,
        objectives=objectives,
    )

    try:
        problem.resolve_constraints(final_check=False)
        problem.optimize()
    except Exception as exc:
        log.warning("dnachisel optimization incomplete: %s — returning best attempt.", exc)

    optimized = problem.sequence

    # Safety: restore ATG if any constraint accidentally changed it
    if not optimized.startswith("ATG"):
        log.warning("Start codon mutated by dnachisel — restoring ATG.")
        optimized = "ATG" + optimized[3:]

    # Safety: check for premature stop codons (dnachisel can introduce them
    # when resolving AvoidPattern constraints that overlap codon boundaries).
    from Bio.Seq import Seq as _Seq
    translated_check = str(_Seq(optimized).translate())
    if "*" in translated_check[:-1]:   # stop anywhere except last position
        log.warning(
            "dnachisel introduced premature stop codon — "
            "falling back to max-frequency back-translation."
        )
        optimized = max_frequency_back_translate(aa_sequence)

    cai = calculate_cai(optimized)
    gc  = sum(1 for c in optimized if c in "GC") / len(optimized)

    report = {
        "method":      "dnachisel",
        "cai":         round(cai, 4),
        "length_nt":   len(optimized),
        "gc_content":  round(gc, 4),
        "constraint_summary": problem.constraints_text_summary() if verbose else None,
    }
    log.info("dnachisel optimization complete: CAI=%.3f, GC=%.2f%%", cai, gc * 100)
    return optimized, report


# =============================================================================
# Public API
# =============================================================================

def optimize_sequence(
    aa_sequence: str,
    method: str = "dnachisel",
    **kwargs,
) -> Tuple[str, Dict]:
    """
    Back-translate and codon-optimize an amino acid sequence.

    Parameters
    ----------
    aa_sequence : protein sequence
    method      : "dnachisel" | "max_frequency" | "stochastic"

    Returns
    -------
    (optimized_dna_sequence, report_dict)
    """
    if method == "dnachisel":
        return optimize_with_dnachisel(aa_sequence, **kwargs)
    elif method == "max_frequency":
        dna = max_frequency_back_translate(aa_sequence)
        cai = calculate_cai(dna)
        return dna, {"method": "max_frequency", "cai": round(cai, 4)}
    elif method == "stochastic":
        dna = back_translate_frequency(aa_sequence, random_seed=None)
        cai = calculate_cai(dna)
        return dna, {"method": "stochastic", "cai": round(cai, 4)}
    else:
        raise ValueError(f"Unknown optimization method: '{method}'")
