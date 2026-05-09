#!/usr/bin/env python3
"""
ClusterCatcher Configuration Generator
======================================

Generate pipeline configuration YAML file with all necessary parameters.

This script creates a complete config.yaml for the ClusterCatcher pipeline,
supporting:
- Cell Ranger alignment
- QC and annotation (Scanpy + popV)
- Dysregulation detection (CytoTRACE2 + inferCNV)
- Viral detection (Kraken2)
- Mutation calling (SComatic)
- Signature analysis (COSMIC NNLS)

Usage:
    python create_config.py --output-dir ./results --sample-pickle samples.pkl \
        --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A
"""

import argparse
import os
import sys
import yaml
from pathlib import Path


def _repo_root():
    """Return the ClusterCatcher repo root (two levels up from this file)."""
    return Path(__file__).parent.parent


def _default_scomatic_scripts():
    """Return the bundled SComatic scripts path from the git submodule."""
    return str(_repo_root() / 'external' / 'SComatic' / 'scripts')


def validate_path(path, name, must_exist=True, create_dir=False):
    """Validate and return absolute path."""
    if path is None:
        return None

    abs_path = os.path.abspath(path)

    if create_dir and not os.path.exists(abs_path):
        os.makedirs(abs_path, exist_ok=True)
        print(f"  Created directory: {abs_path}")
    elif must_exist and not os.path.exists(abs_path):
        raise FileNotFoundError(f"{name} not found: {abs_path}")

    return abs_path


def auto_derive_reference_paths(cellranger_ref):
    """
    Auto-derive FASTA and GTF paths from Cell Ranger reference directory.

    Standard Cell Ranger reference layout:
    {cellranger_ref}/
    ├── fasta/
    │   └── genome.fa
    ├── genes/
    │   └── genes.gtf
    └── star/
        └── ...
    """
    fasta_path = None
    gtf_path = None

    # Try standard locations
    potential_fasta = os.path.join(cellranger_ref, "fasta", "genome.fa")
    if os.path.exists(potential_fasta):
        fasta_path = potential_fasta

    potential_gtf = os.path.join(cellranger_ref, "genes", "genes.gtf")
    if os.path.exists(potential_gtf):
        gtf_path = potential_gtf

    return fasta_path, gtf_path


def create_config(args):
    """Create configuration dictionary and write to YAML file."""

    print("="*70)
    print("CLUSTERCATCHER - CREATE CONFIGURATION")
    print("="*70)

    # Validate and process paths
    print("\nValidating input paths...")

    output_dir = validate_path(args.output_dir, "Output directory",
                               must_exist=False, create_dir=True)

    # Sample information
    sample_info_path = None
    sample_ids = []

    if args.sample_pickle:
        sample_info_path = validate_path(args.sample_pickle, "Sample pickle file")
        import pickle
        with open(sample_info_path, 'rb') as f:
            sample_dict = pickle.load(f)
        sample_ids = list(sample_dict.keys())
        print(f"  Loaded {len(sample_ids)} samples from pickle")
    elif args.sample_ids:
        sample_ids = args.sample_ids
        print(f"  Using {len(sample_ids)} sample IDs from command line")
    else:
        raise ValueError("Must provide either --sample-pickle or --sample-ids")

    # Cell Ranger reference (REQUIRED)
    cellranger_ref = validate_path(args.cellranger_reference, "Cell Ranger reference")
    print(f"  Cell Ranger reference: {cellranger_ref}")

    # Auto-derive fasta and gtf from Cell Ranger reference
    print("\nProcessing reference files...")
    auto_fasta, auto_gtf = auto_derive_reference_paths(cellranger_ref)

    # Use provided paths or auto-derived
    if args.reference_fasta:
        reference_fasta = validate_path(args.reference_fasta, "Reference FASTA")
        print(f"  Reference FASTA (provided): {reference_fasta}")
    elif auto_fasta:
        reference_fasta = auto_fasta
        print(f"  Reference FASTA (auto-derived): {reference_fasta}")
    else:
        reference_fasta = None
        print("  Reference FASTA: Not found (SComatic will not work without it)")

    if args.gtf_file:
        gtf_file = validate_path(args.gtf_file, "GTF annotation file")
        print(f"  GTF annotation (provided): {gtf_file}")
    elif auto_gtf:
        gtf_file = auto_gtf
        print(f"  GTF annotation (auto-derived): {gtf_file}")
    else:
        gtf_file = None
        print("  GTF annotation: Not found (inferCNV may have issues)")

    # Optional modules
    print("\nConfiguring optional modules...")

    # Viral detection
    viral_enabled = args.enable_viral
    kraken_db = validate_path(args.kraken_db, "Kraken2 database") if viral_enabled and args.kraken_db else None
    viral_db = validate_path(args.viral_db, "Viral database") if viral_enabled and args.viral_db else None
    print(f"  Viral detection: {'ENABLED' if viral_enabled else 'DISABLED'}")

    # SComatic mutation calling
    scomatic_enabled = args.enable_scomatic
    scomatic_config = {}
    if scomatic_enabled:
        if reference_fasta is None:
            print("  WARNING: SComatic requires reference FASTA - it was not found or provided")
        scomatic_scripts = args.scomatic_scripts_dir or _default_scomatic_scripts()
        scomatic_config = {
            'scripts_dir': validate_path(scomatic_scripts, "SComatic scripts directory"),
            'editing_sites': validate_path(args.scomatic_editing_sites, "SComatic RNA editing sites"),
            'pon_file': validate_path(args.scomatic_pon_file, "SComatic Panel of Normals"),
            'bed_file': validate_path(args.scomatic_bed_file, "SComatic BED file"),
            'min_cov': args.scomatic_min_cov,
            'min_cells': args.scomatic_min_cells,
            'cell_types': None,
            'min_base_quality': args.scomatic_min_base_quality,
            'min_map_quality': args.scomatic_min_map_quality,
        }
        print(f"  SComatic mutation calling: ENABLED")
    else:
        print(f"  SComatic mutation calling: DISABLED")

    # Signature analysis
    signatures_enabled = args.enable_signatures
    signatures_config = {}
    if signatures_enabled:
        core_signatures = args.core_signatures if args.core_signatures else ['SBS2', 'SBS13', 'SBS5']
        candidate_order = args.candidate_order if args.candidate_order else None

        signatures_config = {
            'cosmic_file': validate_path(args.cosmic_file, "COSMIC signature database"),
            'use_scree_plot': args.use_scree_plot,
            'core_signatures': core_signatures,
            'candidate_order': candidate_order,
            'mutation_threshold': args.mutation_threshold,
            'max_signatures': args.max_signatures,
            'hnscc_only': args.hnscc_only,
            'hnscc_signatures': [
                'SBS1', 'SBS2', 'SBS4', 'SBS5', 'SBS7a', 'SBS7b', 'SBS13',
                'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS29', 'SBS39', 'SBS40', 'SBS44'
            ],
        }
        print(f"  Signature analysis: ENABLED")
        print(f"    Core signatures: {core_signatures}")
        print(f"    Use scree plot: {args.use_scree_plot}")
    else:
        print(f"  Signature analysis: DISABLED")

    # popV annotation settings
    print("\nConfiguring cell type annotation (popV)...")
    print(f"  Hugging Face repo: {args.popv_huggingface_repo}")
    print(f"  Prediction mode: {args.popv_prediction_mode}")

    # QC mode information
    print("\nConfiguring QC filtering...")
    print(f"  QC mode: {args.qc_mode}")
    if args.qc_mode == 'adaptive':
        print(f"    MAD thresholds: counts={args.mad_n_counts}, genes={args.mad_n_genes}, "
              f"top_genes={args.mad_top_genes}, mito={args.mad_mito}")
    else:
        print(f"    Fixed thresholds: min_genes={args.min_genes}, max_genes={args.max_genes}, "
              f"min_counts={args.min_counts}, max_counts={args.max_counts}")
    print(f"  Common: min_cells={args.min_cells}, max_mito_pct={args.max_mito_pct}%")

    # Determine fastq_base_dir
    fastq_base_dir = args.fastq_base_dir
    if not fastq_base_dir and sample_info_path:
        # Try to auto-derive from sample_info path
        # SRAscraper structure: {output_dir}/metadata/dictionary_file.pkl
        #                       {output_dir}/fastq/{GSE_ID}/{SRR_ID}/
        sample_info_dir = os.path.dirname(os.path.abspath(sample_info_path))
        potential_fastq_dir = os.path.join(os.path.dirname(sample_info_dir), "fastq")
        if os.path.exists(potential_fastq_dir):
            fastq_base_dir = potential_fastq_dir
            print(f"\n  Auto-derived fastq_base_dir: {fastq_base_dir}")
        else:
            print(f"\n  WARNING: Could not auto-derive fastq_base_dir")
            print(f"    Checked: {potential_fastq_dir}")
            print(f"    You may need to set --fastq-base-dir manually")

    # Build configuration
    print("\nBuilding configuration...")

    config = {
        # Output directories
        'output_dir': output_dir,
        'log_dir': None,  # Defaults to {output_dir}/logs

        # Resource settings
        'threads': args.threads,
        'memory_gb': args.memory_gb,

        # Sample specification
        'sample_info': sample_info_path,
        'sample_ids': sample_ids,
        'samples': {},

        # FASTQ base directory for SRAscraper format (top level)
        # Structure: {fastq_base_dir}/{GSE_ID}/{SRR_ID}/*.fastq.gz
        'fastq_base_dir': fastq_base_dir,

        # Reference files (simplified)
        'reference': {
            'cellranger': cellranger_ref,
            'fasta': reference_fasta,
            'gtf': gtf_file,
            'genome': args.genome,
        },

        # Cell Ranger settings
        'cellranger': {
            'chemistry': args.chemistry,
            'expect_cells': args.expect_cells,
            'include_introns': args.include_introns,
            'localcores': args.threads,
            'localmem': args.memory_gb,
            'create_bam': True,
        },

        # QC settings (all parameters exposed)
        # Supports two modes: "adaptive" (MAD-based) and "fixed" (threshold-based)
        'qc': {
            # QC mode selection
            'qc_mode': args.qc_mode,

            # Common parameters (used in both modes)
            'min_genes': args.min_genes,
            'max_mito_pct': args.max_mito_pct,
            'min_cells': args.min_cells,
            'doublet_removal': args.doublet_removal,
            'doublet_rate': args.doublet_rate,

            # Adaptive mode parameters (MAD-based outlier detection)
            'mad_n_counts': args.mad_n_counts,
            'mad_n_genes': args.mad_n_genes,
            'mad_top_genes': args.mad_top_genes,
            'mad_mito': args.mad_mito,

            # Fixed mode parameters (threshold-based filtering)
            'max_genes': args.max_genes,
            'min_counts': args.min_counts,
            'max_counts': args.max_counts,
        },

        # Preprocessing settings (post-annotation)
        'preprocessing': {
            'target_sum': args.target_sum,
            'n_pcs': args.n_pcs,
            'n_neighbors': args.n_neighbors,
            'leiden_resolution': args.leiden_resolution,
            'run_bbknn': args.run_bbknn,
            'bbknn_batch_key': args.bbknn_batch_key,
        },

        # Annotation settings (popV only)
        'annotation': {
            'popv_huggingface_repo': args.popv_huggingface_repo,
            'popv_prediction_mode': args.popv_prediction_mode,
            'popv_gene_symbol_key': args.popv_gene_symbol_key,
            'popv_cache_dir': args.popv_cache_dir,
            'batch_key': args.batch_key,
        },

        # Module enable/disable flags
        'modules': {
            'cellranger': True,
            'qc': True,
            'annotation': True,
            'viral': viral_enabled,
            'dysregulation': args.enable_dysregulation,
            'scomatic': scomatic_enabled,
            'signatures': signatures_enabled,
        },

        # Dysregulation settings (order before viral to match workflow)
        'dysregulation': {
            'cytotrace2': {
                'enabled': args.cytotrace2_enabled,
                'species': args.species,
                'max_cells_per_chunk': args.max_cells_chunk,
            },
            'infercnv': {
                'enabled': args.infercnv_enabled,
                'window_size': args.infercnv_window,
            },
            'infercnv_reference_groups': args.infercnv_reference_groups,
            'agreement': {
                'alpha': args.agreement_alpha,
                'min_correlation': args.min_correlation,
            },
        },

        # Viral detection settings
        'viral': {
            'kraken_db': kraken_db,
            'viral_db': viral_db,
            'confidence': args.viral_confidence,
            'include_organisms': None,
            'exclude_organisms': None,
            'organisms_of_interest': [],
        },

        # SComatic settings
        'scomatic': scomatic_config if scomatic_enabled else {
            'scripts_dir': None,
            'editing_sites': None,
            'pon_file': None,
            'bed_file': None,
            'min_cov': 5,
            'min_cells': 5,
            'cell_types': None,
            'min_base_quality': 30,
            'min_map_quality': 30,
        },

        # Signature analysis settings
        'signatures': signatures_config if signatures_enabled else {
            'cosmic_file': None,
            'core_signatures': ['SBS2', 'SBS13', 'SBS5'],
            'use_scree_plot': True,
            'candidate_order': None,
            'mutation_threshold': 0,
            'max_signatures': 15,
            'hnscc_only': False,
            'hnscc_signatures': [
                'SBS1', 'SBS2', 'SBS4', 'SBS5', 'SBS7a', 'SBS7b', 'SBS13',
                'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS29', 'SBS39', 'SBS40', 'SBS44'
            ],
        },

        # Advanced settings
        'random_seed': 42,
        'keep_intermediate': False,

        'reports': {
            'generate_html': True,
            'generate_pdf': False,
            'figure_format': 'png',
            'figure_dpi': 150,
        },
    }

    # Write configuration file
    config_path = os.path.join(output_dir, 'config.yaml')
    print(f"\nWriting configuration to: {config_path}")

    # Custom YAML representer to handle None values nicely
    def represent_none(dumper, data):
        return dumper.represent_scalar('tag:yaml.org,2002:null', 'null')

    yaml.add_representer(type(None), represent_none)

    with open(config_path, 'w') as f:
        f.write("# =============================================================================\n")
        f.write("# ClusterCatcher Pipeline Configuration\n")
        f.write("# =============================================================================\n")
        f.write("# Generated by: create_config.py\n")
        f.write("# Documentation: https://github.com/JakeLehle/ClusterCatcher\n")
        f.write("# =============================================================================\n\n")
        yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

    print("\n" + "="*70)
    print("CONFIGURATION CREATED SUCCESSFULLY")
    print("="*70)
    print(f"\nConfiguration file: {config_path}")
    print(f"Samples: {len(sample_ids)}")
    print(f"\nReference files:")
    print(f"  Cell Ranger: {cellranger_ref}")
    print(f"  FASTA: {reference_fasta or 'Not found'}")
    print(f"  GTF: {gtf_file or 'Not found'}")
    print(f"\nQC Mode: {args.qc_mode}")
    if args.qc_mode == 'adaptive':
        print(f"  Using MAD-based outlier detection (recommended)")
    else:
        print(f"  Using fixed threshold filtering")
    print(f"\nAnnotation: popV ({args.popv_huggingface_repo})")
    print(f"\nTo run the pipeline:")
    print(f"  ClusterCatcher run-config {config_path}")

    return config_path


def main():
    parser = argparse.ArgumentParser(
        description='Create ClusterCatcher pipeline configuration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with adaptive QC (recommended)
  python create_config.py \\
    --output-dir /path/to/output \\
    --sample-pickle samples.pkl \\
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A

  # Use fixed threshold QC mode
  python create_config.py \\
    --output-dir /path/to/output \\
    --sample-pickle samples.pkl \\
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \\
    --qc-mode fixed \\
    --max-genes 5000 --min-counts 500 --max-counts 50000
        """)

    # Required arguments
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('--output-dir', required=True, help='Output directory')
    required.add_argument('--cellranger-reference', required=True,
                          help='Cell Ranger reference directory (contains fasta/, genes/, star/)')

    # Reference files
    refs = parser.add_argument_group('Reference Files (Optional - Auto-derived from Cell Ranger reference)')
    refs.add_argument('--reference-fasta',
                      help='Reference FASTA file (auto-derived as {cellranger}/fasta/genome.fa if not provided)')
    refs.add_argument('--gtf-file',
                      help='GTF annotation file (auto-derived as {cellranger}/genes/genes.gtf if not provided)')
    refs.add_argument('--genome', default='GRCh38', choices=['GRCh38', 'GRCh37', 'mm10', 'mm39'],
                      help='Genome build (default: GRCh38)')

    # Sample specification
    samples = parser.add_argument_group('Sample Specification (one required)')
    samples.add_argument('--sample-pickle', help='Pickle file with sample dictionary')
    samples.add_argument('--sample-ids', nargs='+', help='List of sample IDs')

    # Resource settings
    resources = parser.add_argument_group('Resource Settings')
    resources.add_argument('--threads', type=int, default=8, help='Number of threads (default: 8)')
    resources.add_argument('--memory-gb', type=int, default=64, help='Memory in GB (default: 64)')

    # Cell Ranger settings
    cellranger = parser.add_argument_group('Cell Ranger Settings')
    cellranger.add_argument('--chemistry', default='auto', help='Chemistry (default: auto)')
    cellranger.add_argument('--expect-cells', type=int, default=10000, help='Expected cells (default: 10000)')
    cellranger.add_argument('--include-introns', action='store_true', help='Include introns')
    cellranger.add_argument('--fastq-base-dir',
                            help='Base directory containing FASTQ files.')

    # QC settings
    qc = parser.add_argument_group('QC Settings')
    qc.add_argument('--qc-mode', default='adaptive', choices=['adaptive', 'fixed'],
                    help='QC filtering mode (default: adaptive)')
    qc.add_argument('--min-genes', type=int, default=200,
                    help='Min genes per cell (default: 200)')
    qc.add_argument('--max-mito-pct', type=float, default=20,
                    help='Max mitochondrial %% (default: 20)')
    qc.add_argument('--min-cells', type=int, default=20,
                    help='Min cells expressing a gene (default: 20)')
    qc.add_argument('--doublet-removal', action='store_true', default=True,
                    help='Enable doublet removal (default: True)')
    qc.add_argument('--no-doublet-removal', action='store_false', dest='doublet_removal',
                    help='Disable doublet removal')
    qc.add_argument('--doublet-rate', type=float, default=0.08,
                    help='Expected doublet rate (default: 0.08)')
    qc.add_argument('--mad-n-counts', type=float, default=5,
                    help='MADs for log1p_total_counts outliers (default: 5)')
    qc.add_argument('--mad-n-genes', type=float, default=5,
                    help='MADs for log1p_n_genes_by_counts outliers (default: 5)')
    qc.add_argument('--mad-top-genes', type=float, default=5,
                    help='MADs for pct_counts_in_top_20_genes outliers (default: 5)')
    qc.add_argument('--mad-mito', type=float, default=3,
                    help='MADs for pct_counts_mt outliers (default: 3)')
    qc.add_argument('--max-genes', type=int, default=None,
                    help='Max genes per cell - fixed mode only')
    qc.add_argument('--min-counts', type=int, default=None,
                    help='Min counts per cell - fixed mode only')
    qc.add_argument('--max-counts', type=int, default=None,
                    help='Max counts per cell - fixed mode only')

    # Preprocessing settings
    preproc = parser.add_argument_group('Preprocessing Settings (Post-Annotation)')
    preproc.add_argument('--target-sum', type=int, default=1000000,
                         help='Normalization target sum (default: 1000000 for CPM)')
    preproc.add_argument('--n-pcs', type=int, default=None,
                         help='Number of PCs for neighbors (default: None = CPU count)')
    preproc.add_argument('--n-neighbors', type=int, default=15,
                         help='Number of neighbors for graph (default: 15)')
    preproc.add_argument('--leiden-resolution', type=float, default=1.0,
                         help='Leiden clustering resolution (default: 1.0)')
    preproc.add_argument('--run-bbknn', action='store_true',
                         help='Enable BBKNN batch correction')
    preproc.add_argument('--bbknn-batch-key', default='sample_id',
                         help='Batch key for BBKNN (default: sample_id)')

    # Annotation settings
    annot = parser.add_argument_group('Cell Annotation Settings (popV)')
    annot.add_argument('--popv-huggingface-repo', default='popV/tabula_sapiens_All_Cells',
                       help='Hugging Face repo for popV model')
    annot.add_argument('--popv-prediction-mode', default='inference', choices=['inference', 'fast'],
                       help='popV prediction mode (default: inference)')
    annot.add_argument('--popv-gene-symbol-key', default='feature_name',
                       help='Gene symbol key in adata.var (default: feature_name)')
    annot.add_argument('--popv-cache-dir', default='tmp/popv_models',
                       help='Cache directory for popV models (default: tmp/popv_models)')
    annot.add_argument('--batch-key', default='sample_id',
                       help='Batch key for annotation (default: sample_id)')

    # Dysregulation settings
    dysreg = parser.add_argument_group('Dysregulation Detection')
    dysreg.add_argument('--enable-dysregulation', action='store_true', default=True,
                        help='Enable dysregulation detection (default: True)')
    dysreg.add_argument('--no-dysregulation', action='store_false', dest='enable_dysregulation',
                        help='Disable dysregulation detection')
    dysreg.add_argument('--cytotrace2-enabled', action='store_true', default=True,
                        help='Enable CytoTRACE2 (default: True)')
    dysreg.add_argument('--no-cytotrace2', action='store_false', dest='cytotrace2_enabled',
                        help='Disable CytoTRACE2')
    dysreg.add_argument('--infercnv-enabled', action='store_true', default=True,
                        help='Enable inferCNV (default: True)')
    dysreg.add_argument('--no-infercnv', action='store_false', dest='infercnv_enabled',
                        help='Disable inferCNV')
    dysreg.add_argument('--infercnv-reference-groups', nargs='+',
                        help='Reference cell types for inferCNV (e.g., "T cells" "B cells")')
    dysreg.add_argument('--species', default='human', choices=['human', 'mouse'],
                        help='Species for CytoTRACE2 (default: human)')
    dysreg.add_argument('--max-cells-chunk', type=int, default=200000,
                        help='Max cells per CytoTRACE2 chunk (default: 200000)')
    dysreg.add_argument('--infercnv-window', type=int, default=250,
                        help='InferCNV window size (default: 250)')
    dysreg.add_argument('--agreement-alpha', type=float, default=0.5,
                        help='Agreement weight between rank and value (default: 0.5)')
    dysreg.add_argument('--min-correlation', type=float, default=0.5,
                        help='Minimum Spearman correlation for quartile selection (default: 0.5)')

    # Viral detection
    viral = parser.add_argument_group('Viral Detection')
    viral.add_argument('--enable-viral', action='store_true', help='Enable viral detection')
    viral.add_argument('--kraken-db', help='Kraken2 database path')
    viral.add_argument('--viral-db', help='Viral database path')
    viral.add_argument('--viral-confidence', type=float, default=0.1,
                       help='Viral detection confidence (default: 0.1)')

    # SComatic settings
    scomatic = parser.add_argument_group('SComatic Mutation Calling')
    scomatic.add_argument('--enable-scomatic', action='store_true', help='Enable SComatic mutation calling')
    scomatic.add_argument('--scomatic-scripts-dir',
                          help='SComatic scripts directory (default: external/SComatic/scripts)')
    scomatic.add_argument('--scomatic-editing-sites', help='RNA editing sites file')
    scomatic.add_argument('--scomatic-pon-file', help='Panel of Normals file')
    scomatic.add_argument('--scomatic-bed-file', help='Mappable regions BED file')
    scomatic.add_argument('--scomatic-min-cov', type=int, default=5, help='Min coverage (default: 5)')
    scomatic.add_argument('--scomatic-min-cells', type=int, default=5, help='Min cells with variant (default: 5)')
    scomatic.add_argument('--scomatic-min-base-quality', type=int, default=30,
                          help='Min base quality (default: 30)')
    scomatic.add_argument('--scomatic-min-map-quality', type=int, default=30,
                          help='Min mapping quality (default: 30)')

    # Signature analysis settings
    sigs = parser.add_argument_group('Signature Analysis')
    sigs.add_argument('--enable-signatures', action='store_true', help='Enable signature analysis')
    sigs.add_argument('--cosmic-file', help='COSMIC signature database file')
    sigs.add_argument('--core-signatures', nargs='+', default=['SBS2', 'SBS13', 'SBS5'],
                     help='Core signatures to always include (default: SBS2 SBS13 SBS5)')
    sigs.add_argument('--candidate-order', nargs='+', help='Candidate signature ranking order')
    sigs.add_argument('--use-scree-plot', action='store_true', help='Use scree plot for signature selection')
    sigs.add_argument('--mutation-threshold', type=int, default=0, help='Min mutations per cell (default: 0)')
    sigs.add_argument('--max-signatures', type=int, default=15, help='Max signatures to test (default: 15)')
    sigs.add_argument('--hnscc-only', action='store_true', help='Use only HNSCC-relevant signatures')

    args = parser.parse_args()

    try:
        create_config(args)
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)


import click
from types import SimpleNamespace


@click.command('create-config')
# ---- Required ----
@click.option('--output-dir', required=True, help='Output directory for results and config.yaml')
@click.option('--cellranger-reference', required=True,
              help='Cell Ranger reference directory (contains fasta/, genes/, star/)')
# ---- Reference files ----
@click.option('--reference-fasta', default=None,
              help='Reference FASTA (auto-derived as {cellranger}/fasta/genome.fa if omitted)')
@click.option('--gtf-file', default=None,
              help='GTF annotation file (auto-derived as {cellranger}/genes/genes.gtf if omitted)')
@click.option('--genome', default='GRCh38', type=click.Choice(['GRCh38', 'GRCh37', 'mm10', 'mm39']),
              help='Genome build (default: GRCh38)')
# ---- Sample specification ----
@click.option('--sample-pickle', default=None, help='Pickle file with sample dictionary')
@click.option('--sample-ids', multiple=True, help='Sample IDs (repeat flag for each ID)')
# ---- Resources ----
@click.option('--threads', default=8, type=int, help='Number of threads (default: 8)')
@click.option('--memory-gb', default=64, type=int, help='Memory in GB (default: 64)')
# ---- Cell Ranger ----
@click.option('--chemistry', default='auto', help='10X chemistry (default: auto)')
@click.option('--expect-cells', default=10000, type=int, help='Expected cells (default: 10000)')
@click.option('--include-introns', is_flag=True, default=False, help='Include intronic reads')
@click.option('--fastq-base-dir', default=None,
              help='Base directory for FASTQ files; auto-derived from --sample-pickle if omitted')
# ---- QC ----
@click.option('--qc-mode', default='adaptive', type=click.Choice(['adaptive', 'fixed']),
              help='QC filtering mode: adaptive (MAD-based, default) or fixed (threshold-based)')
@click.option('--min-genes', default=200, type=int, help='Min genes per cell (default: 200)')
@click.option('--max-mito-pct', default=20.0, type=float, help='Max mitochondrial %% (default: 20)')
@click.option('--min-cells', default=20, type=int, help='Min cells expressing a gene (default: 20)')
@click.option('--doublet-removal/--no-doublet-removal', default=True, help='Enable doublet removal (default: on)')
@click.option('--doublet-rate', default=0.08, type=float, help='Expected doublet rate (default: 0.08)')
@click.option('--mad-n-counts', default=5.0, type=float, help='MADs for total_counts outliers (default: 5)')
@click.option('--mad-n-genes', default=5.0, type=float, help='MADs for n_genes outliers (default: 5)')
@click.option('--mad-top-genes', default=5.0, type=float, help='MADs for top-gene pct outliers (default: 5)')
@click.option('--mad-mito', default=3.0, type=float, help='MADs for mito pct outliers (default: 3)')
@click.option('--max-genes', default=None, type=int, help='Max genes per cell - fixed mode only')
@click.option('--min-counts', default=None, type=int, help='Min counts per cell - fixed mode only')
@click.option('--max-counts', default=None, type=int, help='Max counts per cell - fixed mode only')
# ---- Preprocessing ----
@click.option('--target-sum', default=1000000, type=int, help='Normalization target sum (default: 1e6 CPM)')
@click.option('--n-pcs', default=None, type=int, help='PCs for neighbor graph (default: CPU count)')
@click.option('--n-neighbors', default=15, type=int, help='Neighbors for graph (default: 15)')
@click.option('--leiden-resolution', default=1.0, type=float, help='Leiden resolution (default: 1.0)')
@click.option('--run-bbknn', is_flag=True, default=False, help='Enable BBKNN batch correction')
@click.option('--bbknn-batch-key', default='sample_id', help='Batch key for BBKNN (default: sample_id)')
# ---- Annotation ----
@click.option('--popv-huggingface-repo', default='popV/tabula_sapiens_All_Cells',
              help='Hugging Face repo for popV model')
@click.option('--popv-prediction-mode', default='inference', type=click.Choice(['inference', 'fast']),
              help='popV prediction mode (default: inference)')
@click.option('--popv-gene-symbol-key', default='feature_name', help='Gene symbol key in adata.var')
@click.option('--popv-cache-dir', default='tmp/popv_models', help='Cache directory for popV models')
@click.option('--batch-key', default='sample_id', help='Batch key for annotation (default: sample_id)')
# ---- Dysregulation ----
@click.option('--enable-dysregulation/--no-dysregulation', default=True,
              help='Enable dysregulation detection (default: on)')
@click.option('--cytotrace2-enabled/--no-cytotrace2', default=True, help='Enable CytoTRACE2 (default: on)')
@click.option('--infercnv-enabled/--no-infercnv', default=True, help='Enable inferCNV (default: on)')
@click.option('--infercnv-reference-groups', multiple=True,
              help='Reference cell types for inferCNV (repeat flag for each)')
@click.option('--species', default='human', type=click.Choice(['human', 'mouse']),
              help='Species for CytoTRACE2 (default: human)')
@click.option('--max-cells-chunk', default=200000, type=int,
              help='Max cells per CytoTRACE2 chunk (default: 200000)')
@click.option('--infercnv-window', default=250, type=int, help='InferCNV window size (default: 250)')
@click.option('--agreement-alpha', default=0.5, type=float,
              help='Agreement weight between rank and value (default: 0.5)')
@click.option('--min-correlation', default=0.5, type=float,
              help='Min Spearman correlation for quartile selection (default: 0.5)')
# ---- Viral detection ----
@click.option('--enable-viral', is_flag=True, default=False, help='Enable viral detection')
@click.option('--kraken-db', default=None, help='Kraken2 database path')
@click.option('--viral-db', default=None, help='Viral database path')
@click.option('--viral-confidence', default=0.1, type=float,
              help='Viral detection confidence threshold (default: 0.1)')
# ---- SComatic ----
@click.option('--enable-scomatic', is_flag=True, default=False, help='Enable SComatic mutation calling')
@click.option('--scomatic-scripts-dir', default=None,
              help='SComatic scripts directory (default: external/SComatic/scripts from bundled submodule)')
@click.option('--scomatic-editing-sites', default=None, help='RNA editing sites file')
@click.option('--scomatic-pon-file', default=None, help='Panel of Normals file')
@click.option('--scomatic-bed-file', default=None, help='Mappable regions BED file')
@click.option('--scomatic-min-cov', default=5, type=int, help='Min coverage (default: 5)')
@click.option('--scomatic-min-cells', default=5, type=int, help='Min cells with variant (default: 5)')
@click.option('--scomatic-min-base-quality', default=30, type=int, help='Min base quality (default: 30)')
@click.option('--scomatic-min-map-quality', default=30, type=int, help='Min mapping quality (default: 30)')
# ---- Signatures ----
@click.option('--enable-signatures', is_flag=True, default=False, help='Enable mutational signature analysis')
@click.option('--cosmic-file', default=None, help='COSMIC signature database file')
@click.option('--core-signatures', multiple=True, default=['SBS2', 'SBS13', 'SBS5'],
              help='Core signatures to always include (repeat flag for each, default: SBS2 SBS13 SBS5)')
@click.option('--candidate-order', multiple=True, help='Candidate signature ranking order (repeat flag for each)')
@click.option('--use-scree-plot', is_flag=True, default=False, help='Use scree plot for signature selection')
@click.option('--mutation-threshold', default=0, type=int, help='Min mutations per cell (default: 0)')
@click.option('--max-signatures', default=15, type=int, help='Max signatures to test (default: 15)')
@click.option('--hnscc-only', is_flag=True, default=False, help='Use only HNSCC-relevant signatures')
def create_config_cmd(
    output_dir, cellranger_reference, reference_fasta, gtf_file, genome,
    sample_pickle, sample_ids, threads, memory_gb,
    chemistry, expect_cells, include_introns, fastq_base_dir,
    qc_mode, min_genes, max_mito_pct, min_cells, doublet_removal, doublet_rate,
    mad_n_counts, mad_n_genes, mad_top_genes, mad_mito,
    max_genes, min_counts, max_counts,
    target_sum, n_pcs, n_neighbors, leiden_resolution, run_bbknn, bbknn_batch_key,
    popv_huggingface_repo, popv_prediction_mode, popv_gene_symbol_key, popv_cache_dir, batch_key,
    enable_dysregulation, cytotrace2_enabled, infercnv_enabled, infercnv_reference_groups,
    species, max_cells_chunk, infercnv_window, agreement_alpha, min_correlation,
    enable_viral, kraken_db, viral_db, viral_confidence,
    enable_scomatic, scomatic_scripts_dir, scomatic_editing_sites, scomatic_pon_file,
    scomatic_bed_file, scomatic_min_cov, scomatic_min_cells,
    scomatic_min_base_quality, scomatic_min_map_quality,
    enable_signatures, cosmic_file, core_signatures, candidate_order,
    use_scree_plot, mutation_threshold, max_signatures, hnscc_only,
):
    """
    Generate pipeline configuration file.

    \b
    Required arguments:
      --output-dir            Directory where config.yaml and results will be written
      --cellranger-reference  Cell Ranger reference directory

    \b
    One of the following must also be provided:
      --sample-pickle         Pickle file produced by sample-information (or SRAscraper)
      --sample-ids            One or more sample IDs (repeat flag: --sample-ids S1 --sample-ids S2)

    \b
    Minimal example:
      ClusterCatcher create-config \\
          --output-dir ./results \\
          --sample-pickle samples.pkl \\
          --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A

    \b
    With optional modules:
      ClusterCatcher create-config \\
          --output-dir ./results \\
          --sample-pickle samples.pkl \\
          --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \\
          --enable-viral --kraken-db /path/to/kraken2_db \\
          --enable-scomatic \\
          --enable-signatures --cosmic-file /path/to/COSMIC_v3.4_SBS_GRCh38.txt
    """
    args = SimpleNamespace(
        output_dir=output_dir,
        cellranger_reference=cellranger_reference,
        reference_fasta=reference_fasta,
        gtf_file=gtf_file,
        genome=genome,
        sample_pickle=sample_pickle,
        sample_ids=list(sample_ids) if sample_ids else None,
        threads=threads,
        memory_gb=memory_gb,
        chemistry=chemistry,
        expect_cells=expect_cells,
        include_introns=include_introns,
        fastq_base_dir=fastq_base_dir,
        qc_mode=qc_mode,
        min_genes=min_genes,
        max_mito_pct=max_mito_pct,
        min_cells=min_cells,
        doublet_removal=doublet_removal,
        doublet_rate=doublet_rate,
        mad_n_counts=mad_n_counts,
        mad_n_genes=mad_n_genes,
        mad_top_genes=mad_top_genes,
        mad_mito=mad_mito,
        max_genes=max_genes,
        min_counts=min_counts,
        max_counts=max_counts,
        target_sum=target_sum,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        leiden_resolution=leiden_resolution,
        run_bbknn=run_bbknn,
        bbknn_batch_key=bbknn_batch_key,
        popv_huggingface_repo=popv_huggingface_repo,
        popv_prediction_mode=popv_prediction_mode,
        popv_gene_symbol_key=popv_gene_symbol_key,
        popv_cache_dir=popv_cache_dir,
        batch_key=batch_key,
        enable_dysregulation=enable_dysregulation,
        cytotrace2_enabled=cytotrace2_enabled,
        infercnv_enabled=infercnv_enabled,
        infercnv_reference_groups=list(infercnv_reference_groups) if infercnv_reference_groups else None,
        species=species,
        max_cells_chunk=max_cells_chunk,
        infercnv_window=infercnv_window,
        agreement_alpha=agreement_alpha,
        min_correlation=min_correlation,
        enable_viral=enable_viral,
        kraken_db=kraken_db,
        viral_db=viral_db,
        viral_confidence=viral_confidence,
        enable_scomatic=enable_scomatic,
        scomatic_scripts_dir=scomatic_scripts_dir,
        scomatic_editing_sites=scomatic_editing_sites,
        scomatic_pon_file=scomatic_pon_file,
        scomatic_bed_file=scomatic_bed_file,
        scomatic_min_cov=scomatic_min_cov,
        scomatic_min_cells=scomatic_min_cells,
        scomatic_min_base_quality=scomatic_min_base_quality,
        scomatic_min_map_quality=scomatic_min_map_quality,
        enable_signatures=enable_signatures,
        cosmic_file=cosmic_file,
        core_signatures=list(core_signatures) if core_signatures else ['SBS2', 'SBS13', 'SBS5'],
        candidate_order=list(candidate_order) if candidate_order else None,
        use_scree_plot=use_scree_plot,
        mutation_threshold=mutation_threshold,
        max_signatures=max_signatures,
        hnscc_only=hnscc_only,
    )

    try:
        create_config(args)
    except Exception as e:
        click.echo(f"\nERROR: {e}", err=True)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
