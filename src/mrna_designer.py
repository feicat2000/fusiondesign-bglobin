"""
mrna_designer.py — Assemble the complete mRNA sequence.

Full mRNA architecture (5'→3'):
  [m7G Cap] ─ [5'UTR] ─ [Kozak] ─ ATG ─ [Optimized Fusion CDS] ─ [Stop] ─ [3'UTR] ─ [PolyA]
                                    ↑
                          Provided by Kozak prefix (GCCACCATG)
                          The ATG from Kozak IS the initiator Met of the fusion protein.

Notes on cap:
  - m7GpppN cap is added enzymatically during IVT — not encoded in the DNA template.
  - For CleanCap® AG (TriLink): use AGATCGGAAGAGC... as the first two nucleotides.
  - For ARCA cap analog: template starts with G.
  - This file outputs the IVT template DNA (T7 promoter not included).

Notes on UTR choice:
  - HBB 5'/3' UTRs: gold-standard, used in BNT111, Moderna mRNA-1273 variants.
  - Alpha-globin 3'UTR (HBA1): additional stability option.
  - AES+mtRNR1 3'UTR: used in some BioNTech constructs (Orlandini von Niessen 2019).

References:
  Karikó et al. (2008) Mol Ther — pseudouridine-modified mRNA stability
  Orlandini von Niessen et al. (2019) Nat Commun — UTR optimization
  Sahin et al. (2014) Nat Rev Drug Discov — mRNA therapeutics
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, Tuple

from config import (
    HBB_5UTR,
    HBB_3UTR,
    KOZAK,
    POLY_A,
    POLY_A_LENGTH,
)

log = logging.getLogger(__name__)

# T7 RNA polymerase promoter (used for IVT; not present in the final mRNA)
T7_PROMOTER = "TAATACGACTCACTATA"

# Nucleotide modification flags (for documentation only — modifications applied during IVT)
MODIFICATION_NOTES = {
    "N1-methylpseudouridine (m1Ψ)": "Replace all U with m1Ψ during IVT for immune evasion and stability.",
    "5-methylcytosine (m5C)":       "Optional; reduces innate immune activation.",
    "CleanCap AG":                  "Cap-1 structure; use TriLink CleanCap AG in IVT reaction.",
}


# =============================================================================
# Data structure
# =============================================================================

@dataclass
class MRNAConstruct:
    """Complete mRNA construct at DNA template level."""
    name:           str
    five_utr:       str
    kozak:          str
    cds:            str          # Optimized coding sequence (no start/stop — see below)
    stop_codon:     str
    three_utr:      str
    poly_a:         str

    # Assembled sequences
    mrna_dna:       str = ""     # Full DNA template (5'→3', sense strand)
    mrna_rna:       str = ""     # RNA sequence (T→U)

    # Positions (0-based in mrna_dna)
    pos_utr5_start:  int = 0
    pos_cds_start:   int = 0
    pos_cds_end:     int = 0
    pos_utr3_start:  int = 0
    pos_polya_start: int = 0

    def gc_content(self) -> float:
        s = self.mrna_dna
        return sum(1 for c in s if c in "GC") / len(s) if s else 0.0

    def length(self) -> int:
        return len(self.mrna_dna)


# =============================================================================
# mRNA assembly
# =============================================================================

def assemble_mrna(
    construct_name:   str,
    optimized_cds:    str,      # Optimized CDS from codon_optimizer (no ATG prefix, no stop)
    fusion_aa_seq:    str,      # Used only for validation
    five_utr:         str  = HBB_5UTR,
    three_utr:        str  = HBB_3UTR,
    kozak:            str  = KOZAK,
    poly_a_length:    int  = POLY_A_LENGTH,
    stop_codon:       str  = "TGA",       # TGA preferred (most common human)
    include_t7:       bool = False,        # Prepend T7 promoter for IVT template
) -> MRNAConstruct:
    """
    Assemble the complete mRNA sequence.

    The `optimized_cds` should encode the full fusion protein INCLUDING the
    initiator Met but NOT the stop codon.  The Kozak ATG is provided by
    the `kozak` parameter and is NOT duplicated.

    Sequence order:
      [T7_PROMOTER] - 5'UTR - KOZAK(ATG) - CDS(minus ATG) - STOP - 3'UTR - PolyA

    Parameters
    ----------
    optimized_cds : DNA sequence for the entire fusion ORF (with leading ATG)
    fusion_aa_seq : expected protein (for stop-codon validation)
    include_t7    : add T7 promoter upstream (for IVT cloning template)
    """
    # Validate: CDS must start with ATG
    cds_upper = optimized_cds.upper()
    if not cds_upper.startswith("ATG"):
        raise ValueError("optimized_cds must start with ATG (the initiator Met codon).")

    # Kozak provides GCCACC prefix; the ATG is part of the CDS
    # So the final order is: 5'UTR + KOZAK(GCCACC) + CDS(ATG...) + STOP + 3'UTR + PolyA
    poly_a = "A" * poly_a_length

    parts = []
    positions = {}
    cursor = 0

    if include_t7:
        parts.append(T7_PROMOTER)
        positions["t7_promoter"] = (cursor, cursor + len(T7_PROMOTER))
        cursor += len(T7_PROMOTER)

    # 5' UTR
    positions["utr5_start"] = cursor
    parts.append(five_utr)
    cursor += len(five_utr)

    # Kozak (GCCACC) immediately upstream of ATG
    parts.append(kozak)
    cursor += len(kozak)

    # CDS (starts with ATG, includes all codons for fusion protein)
    positions["cds_start"] = cursor
    parts.append(cds_upper)
    cursor += len(cds_upper)

    # Stop codon
    parts.append(stop_codon)
    positions["cds_end"] = cursor + len(stop_codon)
    cursor += len(stop_codon)

    # 3' UTR
    positions["utr3_start"] = cursor
    parts.append(three_utr)
    cursor += len(three_utr)

    # Poly-A tail
    positions["polya_start"] = cursor
    parts.append(poly_a)
    cursor += len(poly_a)

    mrna_dna = "".join(parts)
    mrna_rna = mrna_dna.upper().replace("T", "U")

    mc = MRNAConstruct(
        name          = construct_name,
        five_utr      = five_utr,
        kozak         = kozak,
        cds           = cds_upper,
        stop_codon    = stop_codon,
        three_utr     = three_utr,
        poly_a        = poly_a,
        mrna_dna      = mrna_dna,
        mrna_rna      = mrna_rna,
        pos_utr5_start  = positions.get("utr5_start", 0),
        pos_cds_start   = positions["cds_start"],
        pos_cds_end     = positions["cds_end"],
        pos_utr3_start  = positions["utr3_start"],
        pos_polya_start = positions["polya_start"],
    )

    log.info(
        "mRNA assembled: total %d nt | 5'UTR %d nt | CDS %d nt | "
        "3'UTR %d nt | PolyA %d nt | GC %.1f%%",
        mc.length(),
        len(five_utr),
        len(cds_upper) + len(stop_codon),
        len(three_utr),
        len(poly_a),
        mc.gc_content() * 100,
    )
    return mc


# =============================================================================
# Validation helpers
# =============================================================================

def validate_orf(mc: MRNAConstruct, expected_aa: str) -> Tuple[bool, str]:
    """
    Confirm the CDS in the mRNA construct translates to the expected protein.

    Returns (is_valid, message).
    """
    from Bio.Seq import Seq

    # Extract CDS region (from pos_cds_start to pos_cds_end)
    cds_with_stop = mc.mrna_dna[mc.pos_cds_start: mc.pos_cds_end]
    translated = str(Seq(cds_with_stop).translate(to_stop=True))

    # The expected_aa starts with M (Kozak ATG provides it, CDS starts with ATG)
    if translated == expected_aa:
        return True, f"ORF validates: {len(translated)} aa."
    else:
        return False, (
            f"ORF mismatch! Expected {len(expected_aa)} aa, got {len(translated)} aa.\n"
            f"Expected[:20]: {expected_aa[:20]}\n"
            f"Got[:20]:      {translated[:20]}"
        )


def mrna_summary(mc: MRNAConstruct) -> str:
    """Print a human-readable mRNA construct summary."""
    lines = [
        "=" * 60,
        f"mRNA CONSTRUCT: {mc.name}",
        "=" * 60,
        f"Total length      : {mc.length()} nt",
        f"GC content        : {mc.gc_content()*100:.1f}%",
        f"5' UTR            : {len(mc.five_utr)} nt",
        f"Kozak             : {mc.kozak}",
        f"CDS               : {len(mc.cds)} nt ({len(mc.cds)//3} codons)",
        f"Stop codon        : {mc.stop_codon}",
        f"3' UTR            : {len(mc.three_utr)} nt",
        f"Poly-A tail       : {len(mc.poly_a)} nt",
        "",
        "Position map (0-based, DNA template):",
        f"  5'UTR  : 0 – {mc.pos_utr5_start + len(mc.five_utr)}",
        f"  CDS    : {mc.pos_cds_start} – {mc.pos_cds_end}",
        f"  3'UTR  : {mc.pos_utr3_start} – {mc.pos_utr3_start + len(mc.three_utr)}",
        f"  PolyA  : {mc.pos_polya_start} – {mc.pos_polya_start + len(mc.poly_a)}",
        "",
        "Modification notes (apply during IVT):",
    ]
    for mod, note in MODIFICATION_NOTES.items():
        lines.append(f"  [{mod}] {note}")
    lines.append("=" * 60)
    return "\n".join(lines)


# =============================================================================
# File export
# =============================================================================

def export_fasta(mc: MRNAConstruct, path: str, mode: str = "rna") -> None:
    """
    Export the mRNA (RNA or DNA) to a FASTA file.

    Parameters
    ----------
    mode : "rna" (U-based) or "dna" (T-based template)
    """
    seq = mc.mrna_rna if mode == "rna" else mc.mrna_dna
    header = f">{mc.name} | {len(seq)} nt | GC={mc.gc_content()*100:.1f}% | mode={mode}"
    wrapped = "\n".join(seq[i: i+70] for i in range(0, len(seq), 70))
    with open(path, "w") as fh:
        fh.write(f"{header}\n{wrapped}\n")
    log.info("Exported mRNA FASTA (%s) to %s", mode, path)


def export_genbank(mc: MRNAConstruct, path: str) -> None:
    """Export annotated GenBank file with feature annotations."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    record = SeqRecord(
        Seq(mc.mrna_rna),
        id=mc.name,
        name=mc.name[:16],
        description=f"mRNA neoantigen-HBB fusion, {mc.length()} nt",
    )
    record.annotations["molecule_type"] = "RNA"

    features = [
        ("5UTR",  0,                mc.pos_cds_start, "five_prime_UTR"),
        ("CDS",   mc.pos_cds_start, mc.pos_cds_end,   "CDS"),
        ("3UTR",  mc.pos_utr3_start, mc.pos_polya_start, "three_prime_UTR"),
        ("polyA", mc.pos_polya_start, mc.length(),    "polyA_signal"),
    ]
    for label, start, end, ftype in features:
        record.features.append(
            SeqFeature(FeatureLocation(start, end, strand=1),
                       type=ftype,
                       qualifiers={"label": [label]})
        )

    with open(path, "w") as fh:
        SeqIO.write(record, fh, "genbank")
    log.info("Exported GenBank file to %s", path)
