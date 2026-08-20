#!/usr/bin/env python3
"""
pipeline.py — Main orchestration CLI for the mRNA Neoantigen Vaccine Pipeline.

Full pipeline flow:
  1. Parse neoantigen FASTA + HLA alleles
  2. MHC-I epitope scan (mhcflurry / NetMHCpan)
  3. Select top epitopes per neoantigen
  4. Design linker strategy (pepsickle / NetChop / heuristic)
  5. Assemble fusion protein (N-terminal neoantigens → HBB scaffold)
  6. Codon optimize (dnachisel with human usage + structural constraints)
  7. Assemble full mRNA (5'UTR + Kozak + CDS + 3'UTR + PolyA)
  8. mRNA stability analysis (ViennaRNA)
  9. Protein quality analysis (solubility, aggregation, cleavage)
 10. Generate HTML/PDF report + figures

Usage:
  python pipeline.py \\
    --fasta data/example_neoantigens.fasta \\
    --hla HLA-A*02:01,HLA-B*07:02 \\
    --linker AAY \\
    --name MyVaccine \\
    --output output/

Advanced:
  --epitope-select top1     # Take only the best epitope per neoantigen
  --no-hbb-truncate         # Keep full HBB (don't drop N-terminal Met)
  --codon-method dnachisel  # Codon optimization backend
  --ubiquitin-tag           # Prepend ubiquitin for accelerated DRiP
  --netmhcpan /path/to/bin  # Use NetMHCpan binary for cross-validation
  --t7-promoter             # Include T7 promoter in IVT template output
  --no-stability            # Skip ViennaRNA (fast mode)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

# ── Configure logging ──────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(markup=True, rich_tracebacks=True)],
)
log = logging.getLogger("pipeline")
console = Console()

# ── Pipeline imports ───────────────────────────────────────────────────────────
from config import (
    DEFAULT_HLA_ALLELES,
    MHC1_PEPTIDE_LENGTHS,
    DEFAULT_LINKER,
    OUTPUT_DIR,
    REPORT_DIR,
    MRNA_DIR,
    PROTEIN_DIR,
    HBB_CDS,
    HBB_5UTR,
    HBB_3UTR,
    POLY_A_LENGTH,
)
from src.construct_builder import (
    EpitopeEntry,
    FusionConstruct,
    assemble_fusion,
    construct_summary,
    to_fasta as fusion_to_fasta,
    parse_neoantigen_fasta,
)
from src.epitope_analysis import scan_all_neoantigens, select_top_epitopes
from src.codon_optimizer import optimize_sequence, calculate_cai
from src.mrna_designer import assemble_mrna, validate_orf, mrna_summary, export_fasta, export_genbank
from src.stability_analyzer import analyze_mrna_stability
from src.protein_analyzer import analyze_protein
from src.report_generator import generate_html_report


# =============================================================================
# CLI argument parser
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pipeline.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── Input ──────────────────────────────────────────────────────────────────
    inp = p.add_argument_group("Input")
    inp.add_argument(
        "--fasta", required=True, metavar="FILE",
        help="FASTA file with neoantigen amino acid sequences (one record per neoantigen).",
    )
    inp.add_argument(
        "--hla", default=None, metavar="ALLELES",
        help="Comma-separated HLA alleles, e.g. 'HLA-A*02:01,HLA-B*07:02'. "
             "Defaults to a common population panel.",
    )
    inp.add_argument(
        "--peptide-lengths", default="8,9,10", metavar="LENGTHS",
        help="Comma-separated peptide lengths for MHC-I scanning (default: 8,9,10).",
    )

    # ── Epitope selection ──────────────────────────────────────────────────────
    epi = p.add_argument_group("Epitope Selection")
    epi.add_argument(
        "--top-n", type=int, default=1, metavar="N",
        help="Number of top epitopes to include per neoantigen (default: 1).",
    )
    epi.add_argument(
        "--min-presentation", type=float, default=0.5, metavar="SCORE",
        help="Minimum mhcflurry presentation score (default: 0.5).",
    )
    epi.add_argument(
        "--use-full-sequence", action="store_true",
        help="Use the full input neoantigen sequence instead of MHC-predicted sub-peptides.",
    )

    # ── Construct design ───────────────────────────────────────────────────────
    des = p.add_argument_group("Construct Design")
    des.add_argument(
        "--linker", default=DEFAULT_LINKER,
        choices=["AAY", "GPGPG", "EAAAK", "KK", "NONE", "auto"],
        help="Linker between epitopes. 'auto' selects best per junction (default: AAY).",
    )
    des.add_argument(
        "--ubiquitin-tag", action="store_true",
        help="Prepend N-terminal ubiquitin for enhanced DRiP processing.",
    )
    des.add_argument(
        "--no-hbb-truncate", action="store_true",
        help="Keep full HBB sequence (include HBB Met); default removes it.",
    )

    # ── Codon optimization ─────────────────────────────────────────────────────
    copt = p.add_argument_group("Codon Optimization")
    copt.add_argument(
        "--codon-method", default="dnachisel",
        choices=["dnachisel", "max_frequency", "stochastic"],
        help="Codon optimization backend (default: dnachisel).",
    )

    # ── mRNA design ────────────────────────────────────────────────────────────
    mrna = p.add_argument_group("mRNA Design")
    mrna.add_argument("--t7-promoter", action="store_true",
                      help="Include T7 promoter in IVT DNA template output.")
    mrna.add_argument("--stop-codon", default="TGA", choices=["TGA", "TAA", "TAG"],
                      help="Stop codon to use (default: TGA).")
    mrna.add_argument("--poly-a-length", type=int, default=POLY_A_LENGTH,
                      help=f"Poly-A tail length (default: {POLY_A_LENGTH}).")

    # ── Analysis ───────────────────────────────────────────────────────────────
    ana = p.add_argument_group("Analysis")
    ana.add_argument("--no-stability", action="store_true",
                     help="Skip ViennaRNA stability analysis (faster).")
    ana.add_argument("--netmhcpan", default=None, metavar="PATH",
                     help="Path to netMHCpan binary for cross-validation.")

    # ── Output ─────────────────────────────────────────────────────────────────
    out = p.add_argument_group("Output")
    out.add_argument("--name", default="NeoAg_HBB", metavar="NAME",
                     help="Construct name used in filenames and report (default: NeoAg_HBB).")
    out.add_argument("--output", default=OUTPUT_DIR, metavar="DIR",
                     help=f"Output directory (default: {OUTPUT_DIR}).")
    out.add_argument("--verbose", action="store_true",
                     help="Enable verbose/debug logging.")

    return p


# =============================================================================
# Utility
# =============================================================================

def parse_alleles(hla_str: Optional[str]) -> List[str]:
    if not hla_str:
        return DEFAULT_HLA_ALLELES
    return [a.strip() for a in hla_str.split(",") if a.strip()]


def parse_lengths(lengths_str: str) -> List[int]:
    return [int(x.strip()) for x in lengths_str.split(",")]


def print_banner():
    console.print("""
[bold blue]╔══════════════════════════════════════════════════════════════╗[/]
[bold blue]║  mRNA Neoantigen Vaccine Design Pipeline  v1.0.0             ║[/]
[bold blue]║  N-terminal Beta-Globin Fusion for MHC-I Presentation        ║[/]
[bold blue]╚══════════════════════════════════════════════════════════════╝[/]
""")


# =============================================================================
# Main pipeline
# =============================================================================

def run_pipeline(args: argparse.Namespace) -> Dict:
    """
    Execute the full pipeline and return a results dict.
    """
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    print_banner()
    results: Dict = {"construct_name": args.name, "timestamp": datetime.now().isoformat()}

    # ── Ensure output dirs ─────────────────────────────────────────────────────
    out_root   = Path(args.output)
    mrna_dir   = out_root / "mrna"
    prot_dir   = out_root / "constructs"
    report_dir = out_root / "reports" / args.name
    for d in [mrna_dir, prot_dir, report_dir]:
        d.mkdir(parents=True, exist_ok=True)

    alleles        = parse_alleles(args.hla)
    peptide_lengths = parse_lengths(args.peptide_lengths)
    linker_arg     = None if args.linker == "auto" else args.linker

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:

        # ── STEP 1: Load neoantigens ───────────────────────────────────────────
        t1 = progress.add_task("[cyan]Step 1/8  Loading neoantigen sequences …", total=None)
        neoantigen_seqs = parse_neoantigen_fasta(args.fasta)
        console.print(f"  [green]✓[/] Loaded {len(neoantigen_seqs)} neoantigen sequences.")
        progress.remove_task(t1)

        # ── STEP 2: MHC-I epitope scan ─────────────────────────────────────────
        t2 = progress.add_task("[cyan]Step 2/8  Scanning MHC-I epitopes …", total=None)

        if args.use_full_sequence:
            console.print("  [yellow]Using full neoantigen sequences (--use-full-sequence).[/]")
            epitope_df      = None
            epitope_entries = [
                EpitopeEntry(name=name, sequence=seq)
                for name, seq in neoantigen_seqs.items()
            ]
        else:
            epitope_df = scan_all_neoantigens(
                neoantigen_seqs,
                alleles=alleles,
                peptide_lengths=peptide_lengths,
                netmhcpan_path=args.netmhcpan,
            )
            epitope_df.to_csv(str(out_root / "analysis" / "epitope_scan.csv"), index=False)

            top_df = select_top_epitopes(
                epitope_df,
                top_n=args.top_n,
                min_presentation_score=args.min_presentation,
            )

            # Build EpitopeEntry list — one per neoantigen, best peptide
            epitope_entries = []
            for name in neoantigen_seqs:
                subset = top_df[top_df["neoantigen"] == name]
                if subset.empty:
                    console.print(f"  [yellow]⚠ No qualifying epitope for '{name}' — using full sequence.[/]")
                    epitope_entries.append(EpitopeEntry(name=name, sequence=neoantigen_seqs[name]))
                else:
                    best = subset.iloc[0]
                    epitope_entries.append(EpitopeEntry(
                        name=name,
                        sequence=best["peptide"],
                        hla_allele=best.get("allele", ""),
                        affinity_nM=best.get("affinity_nM", 0.0),
                    ))
                    console.print(
                        f"  [green]✓[/] {name}: {best['peptide']} "
                        f"({best.get('allele','?')}, "
                        f"{best.get('affinity_nM',0):.0f} nM)"
                    )
        progress.remove_task(t2)

        # ── STEP 3: Assemble fusion protein ────────────────────────────────────
        t3 = progress.add_task("[cyan]Step 3/8  Assembling fusion construct …", total=None)
        fusion = assemble_fusion(
            epitopes=epitope_entries,
            fixed_linker=linker_arg,
            hbb_cds=HBB_CDS,
            ubiquitin_tag=args.ubiquitin_tag,
        )
        console.print(f"  [green]✓[/] Fusion protein: {len(fusion.full_aa)} aa")

        # Save fusion protein
        fasta_path = str(prot_dir / f"{args.name}_fusion.fasta")
        with open(fasta_path, "w") as fh:
            fh.write(fusion_to_fasta(fusion, args.name))
        console.print(f"  [green]✓[/] Fusion FASTA: {fasta_path}")
        console.print(construct_summary(fusion))
        progress.remove_task(t3)

        # ── STEP 4: Codon optimization ─────────────────────────────────────────
        t4 = progress.add_task("[cyan]Step 4/8  Codon optimizing …", total=None)
        optimized_cds, codon_report = optimize_sequence(
            fusion.full_aa,
            method=args.codon_method,
        )
        results["codon_report"] = {k: v for k, v in codon_report.items() if k != "constraint_breach"}
        console.print(
            f"  [green]✓[/] Codon optimization: CAI={codon_report.get('cai', 0):.3f}, "
            f"GC={codon_report.get('gc_content', 0)*100:.1f}%"
        )
        progress.remove_task(t4)

        # ── STEP 5: Assemble mRNA ──────────────────────────────────────────────
        t5 = progress.add_task("[cyan]Step 5/8  Assembling mRNA …", total=None)
        mc = assemble_mrna(
            construct_name=args.name,
            optimized_cds=optimized_cds,
            fusion_aa_seq=fusion.full_aa,
            poly_a_length=args.poly_a_length,
            stop_codon=args.stop_codon,
            include_t7=args.t7_promoter,
        )

        # Validate ORF
        valid, msg = validate_orf(mc, fusion.full_aa)
        if valid:
            console.print(f"  [green]✓[/] ORF validated: {msg}")
        else:
            console.print(f"  [red]✗ ORF mismatch: {msg}[/]")
            log.warning("ORF validation failed — check codon optimizer output.")

        console.print(mrna_summary(mc))

        # Export mRNA
        rna_path = str(mrna_dir / f"{args.name}_mrna.fasta")
        dna_path = str(mrna_dir / f"{args.name}_ivt_template.fasta")
        gb_path  = str(mrna_dir / f"{args.name}_mrna.gb")
        export_fasta(mc, rna_path, mode="rna")
        export_fasta(mc, dna_path, mode="dna")
        export_genbank(mc, gb_path)
        console.print(f"  [green]✓[/] mRNA FASTA (RNA): {rna_path}")
        console.print(f"  [green]✓[/] IVT template DNA: {dna_path}")
        console.print(f"  [green]✓[/] GenBank: {gb_path}")
        progress.remove_task(t5)

        # ── STEP 6: mRNA stability ─────────────────────────────────────────────
        t6 = progress.add_task("[cyan]Step 6/8  mRNA stability analysis …", total=None)
        if args.no_stability:
            console.print("  [yellow]Stability analysis skipped (--no-stability).[/]")
            stability_report = {}
        else:
            stability_report = analyze_mrna_stability(
                mrna_rna  = mc.mrna_rna,
                utr5_rna  = mc.five_utr.replace("T", "U"),
                utr3_seq  = mc.three_utr,
                cds_start = mc.pos_cds_start,
                cds_end   = mc.pos_cds_end,
                run_full_fold = len(mc.mrna_rna) < 5000,  # full fold only for short constructs
            )
            mfe = stability_report.get("global_mfe", {})
            gc  = stability_report.get("gc_profile", {})
            console.print(
                f"  [green]✓[/] MFE: {mfe.get('mfe_kcal', 'N/A')} kcal/mol | "
                f"GC: {gc.get('global_gc_pct', '?')}% | "
                f"{mfe.get('flag', '')}"
            )
            # Serialize (drop large window lists for JSON)
            stab_json = {
                k: {kk: vv for kk, vv in v.items() if kk != "data"}
                if isinstance(v, dict) else v
                for k, v in stability_report.items()
            }
            with open(str(out_root / "analysis" / f"{args.name}_stability.json"), "w") as fh:
                json.dump(stab_json, fh, indent=2, default=str)
        progress.remove_task(t6)

        # ── STEP 7: Protein analysis ───────────────────────────────────────────
        t7 = progress.add_task("[cyan]Step 7/8  Protein quality analysis …", total=None)
        protein_report = analyze_protein(
            aa_sequence=fusion.full_aa,
            feature_map=fusion.feature_map,
            hla_alleles=alleles,
        )
        pc = protein_report.get("physicochemical", {})
        console.print(
            f"  [green]✓[/] Stability: {'stable' if pc.get('is_stable') else 'UNSTABLE'} "
            f"(II={pc.get('instability_index', '?'):.1f}) | "
            f"GRAVY={pc.get('gravy', '?'):.3f} | "
            f"pI={pc.get('isoelectric_point', '?'):.2f}"
        )
        sol = protein_report.get("solubility", {})
        console.print(f"  {sol.get('flag', '')}")
        progress.remove_task(t7)

        # ── STEP 8: Generate report ────────────────────────────────────────────
        t8 = progress.add_task("[cyan]Step 8/8  Generating report …", total=None)
        html_path = generate_html_report(
            construct_name   = args.name,
            mrna_construct   = mc,
            fusion_construct = fusion,
            epitope_df       = epitope_df if epitope_df is not None else
                               __import__("pandas").DataFrame(),
            protein_report   = protein_report,
            stability_report = stability_report,
            codon_report     = codon_report,
            output_dir       = str(report_dir),
        )
        progress.remove_task(t8)

    # ── Summary table ──────────────────────────────────────────────────────────
    table = Table(title="Pipeline Complete", border_style="blue")
    table.add_column("Metric",   style="bold")
    table.add_column("Value",    style="green")
    table.add_column("File",     style="dim")

    table.add_row("Fusion protein",  f"{len(fusion.full_aa)} aa",         fasta_path)
    table.add_row("mRNA length",     f"{mc.length()} nt",                 rna_path)
    table.add_row("GC content",      f"{mc.gc_content()*100:.1f}%",       "")
    table.add_row("CAI",             f"{codon_report.get('cai', 0):.3f}", "")
    table.add_row("Epitopes used",   str(len(epitope_entries)),            "")
    table.add_row("Report",          "HTML" + ("+PDF" if (report_dir / "report.pdf").exists() else ""),
                  html_path)

    console.print(table)

    results.update({
        "fusion_length_aa": len(fusion.full_aa),
        "mrna_length_nt":   mc.length(),
        "gc_content":       mc.gc_content(),
        "cai":              codon_report.get("cai"),
        "html_report":      html_path,
    })
    return results


# =============================================================================
# Entry point
# =============================================================================

def main():
    parser = build_parser()
    args   = parser.parse_args()

    # Ensure analysis subdir exists
    (Path(args.output) / "analysis").mkdir(parents=True, exist_ok=True)

    try:
        results = run_pipeline(args)
        sys.exit(0)
    except KeyboardInterrupt:
        console.print("\n[red]Interrupted by user.[/]")
        sys.exit(1)
    except Exception as exc:
        log.exception("Pipeline failed: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
