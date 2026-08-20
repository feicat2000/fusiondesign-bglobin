# mRNA Neoantigen Vaccine — Design Rationale

## Construct Architecture

```
5'─[m7G Cap]─[5'UTR (HBB 50 nt)]─[Kozak: GCCACCATG]─[M-Neo1-Lnk-Neo2-Lnk-...-NeoN-Lnk]─[HBB aa 2-147]─[Stop]─[3'UTR (HBB 132 nt)]─[PolyA 120 nt]─3'
                                                         ↑                                    ↑
                                               N-terminal neoantigens                 Stable scaffold
                                               (proteasome-accessible)                (drives expression)
```

## Why N-Terminal Fusion to Beta-Globin?

### 1. DRiP Pathway Efficiency
Defective Ribosomal Products (DRiPs) are the primary source of MHC-I peptides during active translation. N-terminal sequences emerge first from the ribosome exit tunnel, making them immediately available to the 26S proteasome **before** the full protein folds. Fusion of neoantigens at the N-terminus maximizes their exposure to co-translational proteasomal surveillance.

> Reference: Yewdell et al. (2003) Nat Rev Immunol 3:952–961

### 2. Cytoplasmic Localization → Correct Antigen Processing Route
Beta-globin is a **cytoplasmic protein** (no ER signal peptide). This ensures:
- Processing by the **cytosolic 26S proteasome** → 8-10 mer peptides
- **TAP1/2** transport into the ER lumen
- Loading onto **nascent MHC-I/β2m** complexes in the ER
- Surface presentation to CD8+ cytotoxic T cells

### 3. High Expression Levels
HBB is one of the most highly expressed genes in erythroid cells and is efficiently translated in other cell types. As a scaffold, it ensures **high ribosomal occupancy** → more DRiPs → more antigen presentation per transfected cell.

### 4. HBB UTRs for mRNA Stability
The HBB 5'/3' UTRs are gold-standard for mRNA therapeutics:
- **5'UTR**: Efficient ribosome loading (m7G cap recognition)
- **3'UTR**: Contains canonical AATAAA polyadenylation signal; high mRNA half-life
- Used in BNT111, BNT112, and related BioNTech mRNA vaccines

> Reference: Orlandini von Niessen et al. (2019) Nat Commun 10:3872

## Linker Design for Proteasomal Processing

| Linker | Sequence | Function |
|--------|----------|----------|
| **AAY** | Ala-Ala-Tyr | Creates C-terminal Tyr — hydrophobic anchor for proteasomal cleavage; leaves 9-mer epitopes |
| **GPGPG** | Gly-Pro-Gly-Pro-Gly | Flexible spacer; prevents inter-epitope conformational interference; ~15 Å extension |
| **EAAAK** | Glu-Ala-Ala-Ala-Lys | Rigid α-helix spacer; for structural separation when needed |
| **KK** | Lys-Lys | Short cleavage site; cathepsin B accessible |

**Default (AAY)**: Creates an optimal proteasomal C-terminal cleavage site. The Tyr at position P1 is a common anchor for HLA-A*02:01 and HLA-B*07:02, and the preceding Ala residues are good P2/P3 donors for proteasomal cleavage.

## Codon Optimization Strategy

Applied in layers using **dnachisel**:
1. **CAI maximization** — human codon usage (Kazusa database)
2. **Hairpin avoidance** — no stems > 6 bp (prevents ribosomal stalling)
3. **GC windowing** — 45–65% GC per 50-nt window (stability without excessive structure)
4. **Restriction site avoidance** — EcoRI, BamHI, NotI, HindIII (cloning flexibility)
5. **Homopolymer avoidance** — no runs of ≥ 5 identical nucleotides

## mRNA Modification Recommendations

For clinical/therapeutic use, apply during IVT:
- **N1-methylpseudouridine (m1Ψ)**: Replace all U → m1Ψ. Reduces TLR7/8 innate immune activation; increases translation efficiency ~10×
- **CleanCap® AG**: Provides Cap-1 (2'-O-methylation on +1 nucleotide) — required for optimal translation and immune evasion
- **HPLC purification**: Remove dsRNA byproducts from IVT (major immunostimulant)

> Reference: Karikó et al. (2008) Mol Ther; Karikó et al. (2012) Mol Ther

## Pipeline Dependencies

| Package | Version | Role |
|---------|---------|------|
| biopython | ≥1.81 | Sequence manipulation, NCBI Entrez, ProtParam |
| mhcflurry | ≥2.1 | Pan-allele MHC-I binding (neural network) |
| dnachisel | ≥3.2 | Constrained codon optimization |
| python-codon-tables | ≥0.1.12 | Human codon usage tables (Kazusa) |
| ViennaRNA | ≥2.5 | mRNA secondary structure (MFE, ensemble) |
| pepsickle | any | Proteasomal cleavage prediction |
| pandas/numpy | ≥2.0/1.24 | Data handling |
| matplotlib/seaborn | ≥3.7/0.12 | Visualization |
| jinja2 | ≥3.1 | Report templating |
| rich | ≥13.0 | CLI output |

## References

1. Kreiter S et al. (2015) Mutant MHC class II epitopes drive therapeutic immune responses to cancer. *Nature* 520:692–696.
2. Ott PA et al. (2017) An immunogenic personal neoantigen vaccine for patients with melanoma. *Nature* 547:217–221.
3. Sahin U et al. (2017) Personalized RNA mutanome vaccines mobilize poly-specific therapeutic immunity against cancer. *Nature* 547:222–226.
4. O'Donnell TJ et al. (2020) MHCflurry 2.0: Improved pan-allele prediction of MHC class I-presented peptides. *Cell Systems* 11:42–48.
5. Lorenz R et al. (2011) ViennaRNA Package 2.0. *Algorithms Mol Biol* 6:26.
6. Zulkower V & Rosser S (2020) DNA Chisel, a versatile sequence optimizer. *Bioinformatics* 36:4508–4509.
7. Yewdell JW & Nicchitta CV (2006) The DRiP hypothesis decennial: support, controversy, refinement and extension. *Trends Immunol* 27:368–373.
8. Karikó K et al. (2008) Incorporation of pseudouridine into mRNA yields superior nonimmunogenic vector with increased translational capacity. *Mol Ther* 16:1833–1840.
