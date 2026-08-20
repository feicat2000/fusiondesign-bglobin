"""
epitope_analysis.py — MHC-I binding prediction and epitope ranking.

Workflow per neoantigen:
  1. Slide windows of length 8, 9, 10 over the full neoantigen sequence.
  2. Score each peptide with mhcflurry (affinity + presentation score).
  3. Optionally re-score with NetMHCpan binary if available.
  4. Return ranked DataFrame with recommended peptides per HLA allele.

References:
  O'Donnell et al. (2020) Cell Systems — mhcflurry 2.0
  Jurtz et al. (2017) J Immunol — NetMHCpan-4.0
"""

from __future__ import annotations

import itertools
import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Force CPU-only TensorFlow for mhcflurry: avoids missing-libdevice GPU errors
# on systems where the CUDA toolkit is incomplete. Remove these lines if you
# have a fully configured CUDA environment and want GPU acceleration.
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")    # suppress TF info/warnings

import pandas as pd

try:
    from mhcflurry import Class1PresentationPredictor
    _MHCFLURRY_AVAILABLE = True
except ImportError:
    _MHCFLURRY_AVAILABLE = False
    logging.warning("mhcflurry not available — MHC binding scoring disabled.")

from config import (
    DEFAULT_HLA_ALLELES,
    MHC1_PEPTIDE_LENGTHS,
    STRONG_BINDER_nM,
    WEAK_BINDER_nM,
    PRESENTATION_SCORE,
)

log = logging.getLogger(__name__)


# =============================================================================
# Peptide window generation
# =============================================================================

def slide_windows(sequence: str, lengths: List[int]) -> List[Tuple[str, int, int]]:
    """
    Generate all sub-peptides of given lengths from `sequence`.

    Returns
    -------
    list of (peptide, start, end)  — 0-based, end exclusive
    """
    peptides = []
    for length in lengths:
        for start in range(len(sequence) - length + 1):
            peptide = sequence[start: start + length]
            peptides.append((peptide, start, start + length))
    return peptides


# =============================================================================
# mhcflurry prediction
# =============================================================================

def predict_mhcflurry(
    peptides: List[str],
    alleles: List[str],
    n_flanks: Optional[List[str]] = None,
    c_flanks: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Score a list of peptides against HLA alleles using mhcflurry.

    Parameters
    ----------
    peptides : list of peptide strings (8–10-mers)
    alleles  : list of HLA allele strings, e.g. "HLA-A*02:01"
    n_flanks : N-terminal flanking sequences for each peptide (optional)
    c_flanks : C-terminal flanking sequences for each peptide (optional)

    Returns
    -------
    DataFrame with columns:
        peptide, allele, affinity [nM], affinity_percentile,
        presentation_score, binding_class
    """
    if not _MHCFLURRY_AVAILABLE:
        raise RuntimeError("mhcflurry is not installed. Run: pip install mhcflurry && mhcflurry-downloads fetch")

    predictor = Class1PresentationPredictor.load()

    # mhcflurry 2.x Class1PresentationPredictor.predict() API:
    #   alleles must be a dict: {sample_name: [allele, ...]}
    #   Each "sample" is evaluated independently against its allele list.
    #   Passing each allele as its own single-allele sample gives us
    #   per-allele predictions for every peptide.
    alleles_dict = {allele: [allele] for allele in alleles}

    kwargs: dict = dict(peptides=peptides, alleles=alleles_dict)
    if n_flanks:
        kwargs["n_flanks"] = n_flanks
        kwargs["c_flanks"] = c_flanks

    result = predictor.predict(**kwargs)
    # result has columns: peptide, sample_name, allele, affinity,
    #   affinity_percentile, presentation_score, processing_score, best_allele

    # Rename mhcflurry columns to pipeline-standard names
    df = result.rename(columns={
        "sample_name":           "allele",
        "affinity":              "affinity_nM",
        "affinity_percentile":   "affinity_percentile",
        "presentation_score":    "presentation_score",
        "processing_score":      "processing_score",
    })

    # Classify binders
    df["binding_class"] = "non-binder"
    df.loc[df["affinity_nM"] < WEAK_BINDER_nM,   "binding_class"] = "weak_binder"
    df.loc[df["affinity_nM"] < STRONG_BINDER_nM,  "binding_class"] = "strong_binder"

    return df.sort_values(["allele", "affinity_nM"])


# =============================================================================
# NetMHCpan binary (optional, higher accuracy)
# =============================================================================

def predict_netmhcpan(
    peptides: List[str],
    alleles: List[str],
    netmhcpan_path: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """
    Run NetMHCpan binary (if available) for cross-validation.

    Returns None if binary not found (graceful degradation).
    """
    binary = netmhcpan_path or os.environ.get("NETMHCPAN_PATH", "netMHCpan")
    try:
        subprocess.run([binary, "-version"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        log.info("NetMHCpan binary not found — skipping.")
        return None

    rows = []
    with tempfile.TemporaryDirectory() as tmpdir:
        pep_file = Path(tmpdir) / "peptides.txt"
        pep_file.write_text("\n".join(peptides))

        for allele in alleles:
            cmd = [
                binary,
                "-p",            # peptide input
                "-f", str(pep_file),
                "-a", allele,
                "-xls",          # tabular output
                "-xlsfile", str(Path(tmpdir) / "out.xls"),
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                log.warning("NetMHCpan failed for %s: %s", allele, result.stderr[:200])
                continue

            out_file = Path(tmpdir) / "out.xls"
            if out_file.exists():
                df_allele = pd.read_csv(out_file, sep="\t", skiprows=1)
                df_allele["allele"] = allele
                rows.append(df_allele)

    if not rows:
        return None

    df = pd.concat(rows, ignore_index=True)
    # Standardize columns
    df = df.rename(columns={
        "Peptide": "peptide",
        "nM":      "affinity_nM",
        "%Rank_EL": "affinity_percentile",
    })
    return df[["peptide", "allele", "affinity_nM", "affinity_percentile"]].copy()


# =============================================================================
# Per-neoantigen epitope scan
# =============================================================================

def scan_neoantigen(
    name: str,
    sequence: str,
    alleles: List[str] = DEFAULT_HLA_ALLELES,
    peptide_lengths: List[int] = MHC1_PEPTIDE_LENGTHS,
    netmhcpan_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Full MHC-I epitope scan for a single neoantigen sequence.

    Returns
    -------
    DataFrame sorted by presentation_score (desc), with columns:
        neoantigen, peptide, start, end, length, allele,
        affinity_nM, affinity_percentile, presentation_score, binding_class
    """
    windows = slide_windows(sequence, peptide_lengths)
    peptides = [w[0] for w in windows]
    starts   = [w[1] for w in windows]
    ends     = [w[2] for w in windows]

    log.info("Scanning %d peptide windows for '%s' against %d alleles …",
             len(peptides), name, len(alleles))

    df = predict_mhcflurry(peptides, alleles)

    # Add positional metadata
    pos_map = {p: (s, e) for p, s, e in zip(peptides, starts, ends)}
    df["start"]  = df["peptide"].map(lambda p: pos_map.get(p, (None, None))[0])
    df["end"]    = df["peptide"].map(lambda p: pos_map.get(p, (None, None))[1])
    df["length"] = df["peptide"].str.len()
    df["neoantigen"] = name

    # Optionally merge NetMHCpan scores
    nmhc = predict_netmhcpan(peptides, alleles, netmhcpan_path)
    if nmhc is not None:
        df = df.merge(
            nmhc[["peptide", "allele", "affinity_nM"]].rename(
                columns={"affinity_nM": "netmhcpan_nM"}),
            on=["peptide", "allele"],
            how="left",
        )
        log.info("NetMHCpan scores merged.")

    cols_order = [
        "neoantigen", "peptide", "start", "end", "length",
        "allele", "affinity_nM", "affinity_percentile",
        "presentation_score", "binding_class",
    ]
    cols_order = [c for c in cols_order if c in df.columns]
    return df[cols_order].sort_values("presentation_score", ascending=False).reset_index(drop=True)


def scan_all_neoantigens(
    neoantigen_dict: Dict[str, str],
    alleles: List[str] = DEFAULT_HLA_ALLELES,
    peptide_lengths: List[int] = MHC1_PEPTIDE_LENGTHS,
    netmhcpan_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Run `scan_neoantigen` for every entry in `neoantigen_dict`.

    Parameters
    ----------
    neoantigen_dict : {name: amino_acid_sequence}

    Returns
    -------
    Concatenated DataFrame for all neoantigens.
    """
    frames = []
    for name, seq in neoantigen_dict.items():
        df = scan_neoantigen(name, seq, alleles, peptide_lengths, netmhcpan_path)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def select_top_epitopes(
    scan_df: pd.DataFrame,
    top_n: int = 3,
    min_presentation_score: float = PRESENTATION_SCORE,
) -> pd.DataFrame:
    """
    Select top-N epitopes per neoantigen × allele combination.

    Criteria:
      - presentation_score >= threshold
      - binding_class in {strong_binder, weak_binder}
      - highest presentation_score first
    """
    filtered = scan_df[
        (scan_df["presentation_score"] >= min_presentation_score) &
        (scan_df["binding_class"].isin(["strong_binder", "weak_binder"]))
    ].copy()

    top = (
        filtered
        .sort_values("presentation_score", ascending=False)
        .groupby(["neoantigen", "allele"])
        .head(top_n)
        .reset_index(drop=True)
    )
    log.info("Selected %d top epitopes (threshold=%.2f).", len(top), min_presentation_score)
    return top
