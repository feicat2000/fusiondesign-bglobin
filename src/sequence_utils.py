"""
sequence_utils.py — NCBI fetch helpers and sequence verification utilities.

Use these to verify/update the reference sequences in config.py against
the current NCBI database.

Example:
  python -c "
  from src.sequence_utils import verify_hbb_sequences
  verify_hbb_sequences('your@email.com')
  "
"""

from __future__ import annotations

import logging
import time
from typing import Dict, Optional, Tuple

log = logging.getLogger(__name__)


def fetch_ncbi_sequence(
    accession: str,
    email: str,
    rettype: str = "fasta",
    retmode: str = "text",
    start: Optional[int] = None,
    stop: Optional[int] = None,
) -> str:
    """
    Fetch a sequence from NCBI Entrez.

    Parameters
    ----------
    accession : NCBI accession number (e.g. "NM_000518.5")
    email     : required by NCBI API policy
    start/stop: 0-based nucleotide range (for CDS extraction)

    Returns
    -------
    Raw FASTA string
    """
    from Bio import Entrez, SeqIO
    import io

    Entrez.email = email
    kwargs = {"db": "nucleotide", "id": accession, "rettype": rettype, "retmode": retmode}
    if start is not None and stop is not None:
        kwargs["seq_start"] = start + 1   # NCBI uses 1-based
        kwargs["seq_stop"]  = stop

    log.info("Fetching %s from NCBI Entrez …", accession)
    time.sleep(0.34)   # Respect NCBI rate limit (3 req/sec for unregistered)

    with Entrez.efetch(**kwargs) as handle:
        data = handle.read()
    return data


def fetch_hbb_components(email: str) -> Dict[str, str]:
    """
    Fetch human beta-globin (HBB, NM_000518.5) UTRs and CDS from NCBI.

    Verifies and returns:
      {
        "full_mrna": <full mRNA sequence>,
        "utr5":      <5' UTR>,
        "cds":       <coding sequence>,
        "utr3":      <3' UTR>,
      }
    """
    from Bio import Entrez, SeqIO
    import io

    Entrez.email = email
    log.info("Fetching NM_000518.5 (HBB) from NCBI …")
    time.sleep(0.5)

    with Entrez.efetch(
        db="nucleotide", id="NM_000518.5",
        rettype="gb", retmode="text"
    ) as handle:
        record = SeqIO.read(handle, "genbank")

    full_seq = str(record.seq).upper()

    # Extract features
    utr5_seq = cds_seq = utr3_seq = ""
    for feat in record.features:
        if feat.type == "5'UTR":
            utr5_seq = str(feat.extract(record.seq)).upper()
        elif feat.type == "CDS":
            cds_seq  = str(feat.extract(record.seq)).upper()
        elif feat.type == "3'UTR":
            utr3_seq = str(feat.extract(record.seq)).upper()

    return {
        "full_mrna": full_seq,
        "utr5":      utr5_seq,
        "cds":       cds_seq,
        "utr3":      utr3_seq,
    }


def verify_hbb_sequences(email: str) -> None:
    """
    Compare embedded config.py HBB sequences against NCBI NM_000518.5.

    Prints a diff summary; does NOT modify config.py automatically.
    """
    from config import HBB_5UTR, HBB_CDS, HBB_3UTR

    components = fetch_hbb_components(email)

    checks = [
        ("5'UTR",  HBB_5UTR.upper(), components["utr5"]),
        ("CDS",    HBB_CDS.upper(),  components["cds"]),
        ("3'UTR",  HBB_3UTR.upper(), components["utr3"]),
    ]

    all_ok = True
    for name, embedded, fetched in checks:
        if not fetched:
            print(f"  {name}: Could not extract from GenBank record.")
            continue
        if embedded == fetched:
            print(f"  {name}: ✓ Match ({len(embedded)} nt)")
        else:
            all_ok = False
            print(f"  {name}: ✗ MISMATCH")
            print(f"    Embedded length: {len(embedded)}")
            print(f"    NCBI length:     {len(fetched)}")
            # Show first diff
            for i, (a, b) in enumerate(zip(embedded, fetched)):
                if a != b:
                    print(f"    First diff at position {i}: embedded={a}, NCBI={b}")
                    break

    if all_ok:
        print("\nAll HBB sequences match NM_000518.5. config.py is up-to-date.")
    else:
        print("\nUpdate config.py with the fetched sequences above.")


def translate_cds(dna: str) -> str:
    """Translate a DNA CDS to protein (requires Biopython)."""
    from Bio.Seq import Seq
    return str(Seq(dna.upper()).translate(to_stop=True))


def reverse_complement(dna: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    from Bio.Seq import Seq
    return str(Seq(dna.upper()).reverse_complement())
