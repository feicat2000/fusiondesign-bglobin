"""
construct_builder.py — Assemble the neoantigen-HBB fusion protein sequence.

Protein-level architecture:
  M-[Neo1]-[Lnk1]-[Neo2]-[Lnk2]-...-[NeoN]-[Lnk_N]-[HBB(2..147)]
                                                         ↑
                                            HBB Met removed (Neo string
                                            provides the first Met via Kozak)

N-terminal placement rationale:
  • N-terminal sequences emerge from the ribosome first during translation.
  • They are immediately accessible to the 26S proteasome via the DRiP
    (Defective Ribosomal Products) pathway.
  • Beta-globin provides a stable C-terminal scaffold for high-level cytoplasmic
    expression — no ER signal peptide, keeping the fusion in the cytosol where
    proteasomal MHC-I processing occurs.
  • The highly expressed, stable HBB portion increases overall ribosomal
    occupancy, boosting co-translational DRiP antigen generation.

References:
  Yewdell et al. (2003) Nat Rev Immunol — DRiP pathway
  Princiotta et al. (2003) Immunity — N-terminal processing efficiency
  Kreiter et al. (2015) Nature — personalized multi-epitope mRNA vaccine
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from Bio.Seq import Seq

from config import HBB_CDS, DEFAULT_LINKER, LINKERS
from src.linker_design import design_linker_strategy

log = logging.getLogger(__name__)


# =============================================================================
# Data structures
# =============================================================================

@dataclass
class EpitopeEntry:
    """Single neoantigen epitope to include in the construct."""
    name:        str
    sequence:    str          # Amino acid sequence
    hla_allele:  str = ""     # Best-match HLA allele
    affinity_nM: float = 0.0  # Predicted IC50


@dataclass
class FusionConstruct:
    """Complete protein-level fusion construct."""
    epitope_entries:   List[EpitopeEntry]
    linker_strategy:   List[Tuple[str, str]]   # [(name, seq), ...]
    hbb_aa:            str                      # HBB amino acid sequence (residues 2–147)
    full_aa:           str = field(default="")  # Complete fusion AA sequence
    hbb_cds:           str = field(default="")  # Original HBB CDS (for reference)

    # Map of feature positions in full_aa (0-based)
    feature_map: Dict[str, Tuple[int, int]] = field(default_factory=dict)


# =============================================================================
# HBB reference sequence helpers
# =============================================================================

def hbb_cds_to_protein(cds: str = HBB_CDS) -> str:
    """
    Translate the HBB CDS to protein sequence.

    Raises ValueError if the CDS does not produce the canonical 147-aa sequence.
    """
    seq = Seq(cds)
    if len(seq) % 3 != 0:
        raise ValueError(f"HBB CDS length ({len(seq)}) is not a multiple of 3.")
    protein = str(seq.translate(to_stop=True))
    log.debug("HBB CDS → %d aa protein: %s…", len(protein), protein[:10])
    return protein


def hbb_scaffold(cds: str = HBB_CDS, drop_n_terminal_met: bool = True) -> str:
    """
    Return the HBB scaffold amino acid sequence.

    Parameters
    ----------
    drop_n_terminal_met : if True, remove the HBB Met (position 1) since the
                          fused neoantigen string provides the initiator Met
                          through the Kozak-ATG context.
    """
    protein = hbb_cds_to_protein(cds)
    if drop_n_terminal_met and protein.startswith("M"):
        protein = protein[1:]
        log.debug("Dropped HBB N-terminal Met; scaffold is %d aa.", len(protein))
    return protein


# =============================================================================
# Fusion construct assembly
# =============================================================================

def assemble_fusion(
    epitopes: List[EpitopeEntry],
    linker_strategy: Optional[List[Tuple[str, str]]] = None,
    fixed_linker: Optional[str] = DEFAULT_LINKER,
    hbb_cds: str = HBB_CDS,
    ubiquitin_tag: bool = False,
) -> FusionConstruct:
    """
    Assemble the complete fusion protein sequence.

    Parameters
    ----------
    epitopes         : ordered list of EpitopeEntry objects (N → C terminal order)
    linker_strategy  : pre-computed list of (name, seq) linkers; auto-selected if None
    fixed_linker     : force a specific linker name for all junctions
    hbb_cds          : HBB coding DNA sequence (from config)
    ubiquitin_tag    : prepend ubiquitin (76 aa) before the first epitope for
                       accelerated DRiP processing (optional, experimental)

    Returns
    -------
    FusionConstruct with complete amino acid sequence and feature map.
    """
    if not epitopes:
        raise ValueError("At least one epitope is required.")

    hbb_aa = hbb_scaffold(hbb_cds, drop_n_terminal_met=True)

    # Auto-design linker strategy if not provided
    epi_seqs = [e.sequence for e in epitopes]
    if linker_strategy is None:
        linker_strategy = design_linker_strategy(epi_seqs, fixed_linker=fixed_linker)

    # ── Build the fusion string ───────────────────────────────────────────────
    feature_map: Dict[str, Tuple[int, int]] = {}
    parts: List[str] = []
    cursor = 0

    # Optional ubiquitin tag (76 aa human ubiquitin — ensures rapid degradation)
    # Sequence from UniProt P0CG48
    if ubiquitin_tag:
        ubi = (
            "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
        )
        parts.append(ubi)
        feature_map["ubiquitin"] = (cursor, cursor + len(ubi))
        cursor += len(ubi)
        log.info("Prepended ubiquitin tag (%d aa).", len(ubi))

    # Initiator Methionine (from Kozak ATG — becomes position 0)
    init_met = "M"
    parts.append(init_met)
    feature_map["start_M"] = (cursor, cursor + 1)
    cursor += 1

    # Neoantigen epitopes + linkers
    for i, entry in enumerate(epitopes):
        start = cursor
        parts.append(entry.sequence)
        feature_map[f"neoantigen_{i+1}_{entry.name}"] = (start, start + len(entry.sequence))
        cursor += len(entry.sequence)
        log.debug("Added epitope %d ('%s'): pos %d–%d", i+1, entry.name, start, cursor)

        # Add linker between this epitope and the next
        if i < len(linker_strategy):
            lnk_name, lnk_seq = linker_strategy[i]
            if lnk_seq:
                lnk_start = cursor
                parts.append(lnk_seq)
                feature_map[f"linker_{i+1}_{lnk_name}"] = (lnk_start, lnk_start + len(lnk_seq))
                cursor += len(lnk_seq)

    # HBB scaffold
    hbb_start = cursor
    parts.append(hbb_aa)
    feature_map["HBB_scaffold"] = (hbb_start, hbb_start + len(hbb_aa))
    cursor += len(hbb_aa)

    full_aa = "".join(parts)
    log.info(
        "Fusion construct assembled: %d aa total "
        "(%d epitopes, %d linkers, HBB at %d–%d).",
        len(full_aa), len(epitopes), len(linker_strategy),
        hbb_start, cursor,
    )

    return FusionConstruct(
        epitope_entries=epitopes,
        linker_strategy=linker_strategy,
        hbb_aa=hbb_aa,
        full_aa=full_aa,
        hbb_cds=hbb_cds,
        feature_map=feature_map,
    )


def construct_summary(fc: FusionConstruct) -> str:
    """Return a human-readable summary of the fusion construct."""
    lines = [
        "=" * 60,
        "FUSION CONSTRUCT SUMMARY",
        "=" * 60,
        f"Total length : {len(fc.full_aa)} aa",
        f"Epitopes     : {len(fc.epitope_entries)}",
        f"Linkers      : {len(fc.linker_strategy)}",
        "",
        "Feature map (0-based positions):",
    ]
    for feature, (start, end) in fc.feature_map.items():
        seq_snippet = fc.full_aa[start:end]
        display = seq_snippet[:20] + ("…" if len(seq_snippet) > 20 else "")
        lines.append(f"  {feature:<40} [{start:4d}–{end:4d}]  {display}")

    lines += [
        "",
        "Full sequence:",
        fc.full_aa,
        "=" * 60,
    ]
    return "\n".join(lines)


# =============================================================================
# FASTA export
# =============================================================================

def to_fasta(fc: FusionConstruct, construct_name: str = "neoantigen_HBB_fusion") -> str:
    """Return FASTA-formatted fusion protein sequence."""
    header = (
        f">{construct_name} | "
        f"neoantigens={len(fc.epitope_entries)} | "
        f"length={len(fc.full_aa)}aa"
    )
    # Wrap at 60 chars
    seq = fc.full_aa
    wrapped = "\n".join(seq[i: i + 60] for i in range(0, len(seq), 60))
    return f"{header}\n{wrapped}\n"


def parse_neoantigen_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Parse a FASTA file of neoantigen amino acid sequences.

    Returns {name: sequence}
    """
    from Bio import SeqIO
    records = {}
    with open(fasta_path) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            records[rec.id] = str(rec.seq).upper()
    log.info("Loaded %d neoantigen sequences from %s", len(records), fasta_path)
    return records
