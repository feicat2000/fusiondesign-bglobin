# bglobin — mRNA Neoantigen Vaccine Design Pipeline

A pipeline for designing personalized mRNA neoantigen cancer vaccines using a
beta-globin (HBB) fusion scaffold.

## Idea

Neoantigen peptides are fused N-terminally to beta-globin, a highly-expressed
cytoplasmic scaffold protein. This exploits the Defective Ribosomal Products
(DRiP) pathway — N-terminal sequences are exposed to the proteasome
co-translationally, before the fusion protein finishes folding, which
improves MHC-I antigen presentation. HBB's 5'/3' UTRs (used in BioNTech's
BNT111/BNT112) are included for mRNA stability and translation efficiency.

```
5'─[m7G Cap]─[5'UTR]─[Kozak]─[Neo1-Lnk-Neo2-Lnk-...-NeoN-Lnk]─[HBB aa 2-147]─[Stop]─[3'UTR]─[PolyA]─3'
```

See [DESIGN.md](DESIGN.md) for the full design rationale, linker options,
codon-optimization strategy, and references.

## What it does

- `src/epitope_analysis.py` / `src/protein_analyzer.py` — MHC-I binding
  prediction (mhcflurry) and proteasomal cleavage prediction (pepsickle) over
  candidate neoantigens
- `src/linker_design.py` — selects proteasome-friendly linkers (AAY, GPGPG,
  EAAAK, KK) between epitopes
- `src/construct_builder.py` — assembles the full fusion construct
- `src/codon_optimizer.py` — constrained codon optimization via `dnachisel`
  (CAI, hairpin avoidance, GC windowing, restriction-site avoidance)
- `src/mrna_designer.py` — builds the final mRNA (UTRs, Kozak, poly-A) and
  checks secondary structure with ViennaRNA
- `src/stability_analyzer.py`, `src/report_generator.py` — QC and an HTML
  report per construct

## Usage

```bash
./setup.sh              # creates venv, installs requirements
python pipeline.py       # run the pipeline (see config.py for inputs)
```

Inputs go in `data/`, outputs land in `output/` (constructs, mRNA sequences,
analysis, HTML reports).

## Dependencies

See `requirements.txt` — key packages: biopython, mhcflurry, dnachisel,
python-codon-tables, ViennaRNA, pepsickle.
