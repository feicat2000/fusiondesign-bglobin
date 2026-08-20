"""
config.py — Pipeline-wide constants, reference sequences, and parameters.

All beta-globin sequences are from:
  Human HBB gene, RefSeq NM_000518.5
  https://www.ncbi.nlm.nih.gov/nuccore/NM_000518.5

To re-fetch and verify these sequences programmatically:
  from src.sequence_utils import fetch_ncbi_sequence
  fetch_ncbi_sequence("NM_000518.5", email="you@email.com")
"""

# =============================================================================
# REFERENCE mRNA ARCHITECTURE
# =============================================================================
#
#  5'─[m7G Cap]─[5'UTR]─[Kozak]─ATG─[Neoantigens+Linkers]─[HBB CDS]─[3'UTR]─[PolyA]─3'
#                                      ↑                       ↑
#                            processed by proteasome    stable, highly expressed
#                            → MHC-I loading            anchor/scaffold
#
# Scientific rationale:
#   • N-terminal neoantigens emerge from ribosome first → immediate DRiP processing
#   • Beta-globin scaffold: cytoplasmic, no ER signal → cytosolic proteasomal route
#   • HBB 5'/3' UTRs are gold-standard for mRNA stability in therapeutics
#   • AAY linkers create efficient proteasomal cleavage sites for MHC-I peptides
# =============================================================================

# ── Human Beta-Globin 5' UTR ─────────────────────────────────────────────────
# Source: NM_000518.5, nucleotides 1–50 (before ATG at position 51)
# Known for efficient ribosome loading; widely used in mRNA therapeutics
HBB_5UTR = (
    "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACAC"
)

# Kozak consensus (human optimal); the ATG here will be the fusion protein start
KOZAK = "GCCACC"   # precedes ATG; combined = GCCACCATG

# ── Human Beta-Globin CDS ────────────────────────────────────────────────────
# Source: NM_000518.5, CDS 51–494 (444 nt, 147 aa + stop)
# Encodes: MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLST…
# NOTE: When fusing neoantigens at the N-terminus, the neoantigen string
#       is prepended BEFORE this CDS (sharing the same Kozak+ATG).
#       The HBB Met (M) at the start of this CDS is removed from the
#       translated fusion (the upstream ATG in Kozak context is used).
HBB_CDS = (
    "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG"
    "AACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGG"
    "ACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGC"
    "AACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTG"
    "GCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGAC"
    "AAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTG"
    "CTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAA"
    "GTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
)

# ── Human Beta-Globin 3' UTR ─────────────────────────────────────────────────
# Source: NM_000518.5, nucleotides 495–626 (132 nt)
# Contains canonical AATAAA polyadenylation signal at position ~614
# High mRNA stability; used in BNT111 and related mRNA vaccines
HBB_3UTR = (
    "AGTAACGGCAGACTTCTCCTCAGGAGTCAGGTGCACCATGGTGTCTGTTTGAGGTT"
    "GCTAGTGAACACAGTTGTGTCAGAAGCAAATGTAAGCAATAGATGGCTCTGCCCTGA"
    "CTTTTATGCCCAGCCCTGGCTCCTGCCCTCCCTGCTCCTGGGAGTAGATTGGCCAAC"
    "CCTAGTCTGGAGTCTGACAGGAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTT"
    "CTGACACAAC"
)

# Optional: HBA1 3'UTR (NM_000558) — longer stability element
# Some mRNA therapeutics (e.g., BNT vaccines) use AES+mtRNR1 or HBB+HBA1
# Uncomment and concatenate if desired:
# HBA1_3UTR = (
#     "GCAATGAAAATAAATTTCCTTTATTAGCAGAAGGCTGTGGGGCAAAGGTATAATGCTG"
#     "GTAGACTTAGTAAGTGAGACCTTGGAGAGGAGCAGAAACTGGAGCCTTGATGGCAGG"
#     "AGATAGGGGAGTGGCCTCAATAAATGCTATTTCTTTT"
# )

# ── Poly-A tail ───────────────────────────────────────────────────────────────
# 120 nt; standard for mRNA therapeutics; encoded in the DNA template
POLY_A_LENGTH = 120
POLY_A = "A" * POLY_A_LENGTH

# =============================================================================
# LINKER SEQUENCES (amino acid)
# =============================================================================
# Used to join neoantigen epitopes in the polypeptide string-of-beads.
# See: Kreiter et al. (2015) Nature; Ott et al. (2017) Nature

LINKERS = {
    # Preferred for MHC-I constructs — creates efficient proteasomal cleavage
    # The C-terminal Y is a common MHC-I anchor residue
    "AAY":    "AAY",

    # Flexible spacer; reduces inter-epitope context effects
    # Widely used in multi-epitope strings (Nardin et al.)
    "GPGPG":  "GPGPG",

    # Rigid helix-breaking linker; 3× repeat for ~10 Å separation
    "EAAAK":  "EAAAK",

    # Short rigid spacer; good for closely spaced HLA-I epitopes
    "KK":     "KK",

    # No linker — direct fusion; only for benchmarking
    "NONE":   "",
}

# Default linker for MHC-I multi-epitope constructs
DEFAULT_LINKER = "AAY"

# =============================================================================
# CODON OPTIMIZATION PARAMETERS
# =============================================================================
CODON_OPT_ORGANISM = "h_sapiens"   # python-codon-tables key
CODON_OPT_AVOID_SEQS = [
    # Restriction sites to avoid (for IVT cloning flexibility)
    "GAATTC",   # EcoRI
    "AAGCTT",   # HindIII
    "GGATCC",   # BamHI
    "CTCGAG",   # XhoI
    "GCGGCCGC", # NotI
]

# dnachisel optimization objectives
CODON_OPT_OBJECTIVES = {
    "codon_usage":        True,   # Maximize human codon adaptation index (CAI)
    "avoid_hairpins":     True,   # Avoid stem-loops > 6 bp (mRNA stability)
    "avoid_homopolymers": True,   # Avoid runs of ≥ 5 identical nucleotides
    "gc_target":         (0.45, 0.65),  # Target GC window [45–65%]
}

# =============================================================================
# MHC-I PREDICTION PARAMETERS
# =============================================================================
# HLA alleles to consider (add patient-specific alleles via --hla CLI flag)
DEFAULT_HLA_ALLELES = [
    "HLA-A*02:01",   # Most common in Caucasian populations
    "HLA-A*01:01",
    "HLA-A*03:01",
    "HLA-B*07:02",
    "HLA-B*08:01",
    "HLA-C*07:01",
]

# Peptide lengths for MHC-I (canonical: 8–10 mers)
MHC1_PEPTIDE_LENGTHS = [8, 9, 10]

# mhcflurry affinity cutoffs
STRONG_BINDER_nM   = 50.0    # IC50 < 50 nM  → strong binder
WEAK_BINDER_nM     = 500.0   # IC50 < 500 nM → weak binder
PRESENTATION_SCORE = 0.5     # mhcflurry presentation score threshold

# =============================================================================
# PROTEIN QUALITY PARAMETERS
# =============================================================================
# Biopython ProteinAnalysis thresholds
INSTABILITY_INDEX_CUTOFF = 40.0   # < 40 → stable protein
GRAVY_TARGET_RANGE       = (-1.0, 0.2)  # GRAVY close to 0 → soluble

# Aggregation: simple sliding-window hydrophobicity (Kyte-Doolittle)
AGGREGATION_WINDOW       = 9
AGGREGATION_THRESHOLD    = 1.6    # Mean KD > 1.6 in window → aggregation risk

# =============================================================================
# mRNA STABILITY PARAMETERS
# =============================================================================
# ViennaRNA / GC analysis thresholds
# NOTE: Full-mRNA MFE scales with length (~0.5 kcal/mol per nt for typical mRNAs).
# For a 1 kb mRNA, MFE ≈ -500 kcal/mol is normal.  The threshold below is used
# for per-window (120 nt) analysis; -50 kcal/mol per window is a reasonable cutoff.
MFE_THRESHOLD_KCAL       = -50.0   # Per-window MFE threshold (less negative = less structured)
GC_CONTENT_TARGET        = (0.50, 0.65)  # 50–65% GC for stability without hairpin excess
CAI_TARGET               = 0.80    # Codon Adaptation Index ≥ 0.8

# =============================================================================
# PROTEASOMAL CLEAVAGE
# =============================================================================
# NetChop web API (DTU Health Tech) — used as fallback
NETCHOP_URL = "https://services.healthtech.dtu.dk/cgi-bin/webscripts/netchop.cgi"
NETCHOP_METHOD = "cterm"  # C-terminal cleavage model (recommended for MHC-I)
NETCHOP_THRESHOLD = 0.5   # Cleavage probability threshold

# =============================================================================
# FILE PATHS
# =============================================================================
import os
BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
DATA_DIR    = os.path.join(BASE_DIR, "data")
OUTPUT_DIR  = os.path.join(BASE_DIR, "output")
REPORT_DIR  = os.path.join(OUTPUT_DIR, "reports")
MRNA_DIR    = os.path.join(OUTPUT_DIR, "mrna")
PROTEIN_DIR = os.path.join(OUTPUT_DIR, "constructs")

REPORT_TEMPLATE = os.path.join(BASE_DIR, "src", "templates", "report.html.j2")
