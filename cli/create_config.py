#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
create_config.py
================

CLI command to create the ClusterCatcher pipeline configuration file.

Usage:
    python create_config.py [OPTIONS]
    clustercatcher create-config [OPTIONS]

Author: Jake Lehle
Date: 2025
"""

import os
import sys
import yaml
import argparse
from pathlib import Path


def validate_path(path, name, must_exist=True, create_dir=False):
    """Validate a file or directory path."""
    if path is None:
        return None
    
    path = Path(path).resolve()
    
    if create_dir and not path.exists():
        path.mkdir(parents=True, exist_ok=True)
        print(f"  Created directory: {path}")
    
    if must_exist and not path.exists():
        raise FileNotFoundError(f"{name} not found: {path}")
    
    return str(path)


def auto_derive_reference_paths(cellranger_ref):
    """
    Auto-derive fasta and gtf paths from Cell Ranger reference directory.
    
    Standard 10x Genomics reference layout:
        {cellranger}/
        ├── fasta/
        │   └── genome.fa
        ├── genes/
        │   └── genes.gtf
        └── star/
            └── ...
    
    Parameters
    ----------
    cellranger_ref : str
        Path to Cell Ranger reference directory
        
    Returns
    -------
    tuple
        (fasta_path, gtf_path) - paths if found, None if not
    """
    cellranger_path = Path(cellranger_ref)
    
    # Try to find FASTA
    fasta_path = None
    potential_fasta_files = [
        cellranger_path / "fasta" / "genome.fa",
        cellranger_path / "fasta" / "genome.fasta",
        cellranger_path / "fasta" / "reference.fa",
    ]
    for fasta_candidate in potential_fasta_files:
        if fasta_candidate.exists():
            fasta_path = str(fasta_candidate)
            break
    
    # Try to find GTF
    gtf_path = None
    potential_gtf_files = [
        cellranger_path / "genes" / "genes.gtf",
        cellranger_path / "genes" / "annotations.gtf",
        cellranger_path / "genes" / "reference.gtf",
    ]
    for gtf_candidate in potential_gtf_files:
        if gtf_candidate.exists():
            gtf_path = str(gtf_candidate)
            break
    
    return fasta_path, gtf_path


def create_config(args):
    """Create the pipeline configuration file."""
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
        scomatic_config = {
            'scripts_dir': validate_path(args.scomatic_scripts_dir, "SComatic scripts directory"),
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
        'qc': {
            'min_genes': args.min_genes,
            'max_genes': args.max_genes,
            'min_counts': args.min_counts,
            'max_counts': args.max_counts,
            'max_mito_pct': args.max_mito_pct,
            'min_cells': args.min_cells,
            'doublet_removal': args.doublet_removal,
            'doublet_rate': args.doublet_rate,
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
    print(f"\nAnnotation: popV ({args.popv_huggingface_repo})")
    print(f"\nTo run the pipeline:")
    print(f"  cd snakemake_wrapper")
    print(f"  snakemake --configfile {config_path} --cores {args.threads} --use-conda")
    
    return config_path


def main():
    parser = argparse.ArgumentParser(
        description='Create ClusterCatcher pipeline configuration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - only cellranger-reference is required for references
  # FASTA and GTF are auto-derived from standard 10x reference layout
  python create_config.py \\
    --output-dir /path/to/output \\
    --sample-pickle samples.pkl \\
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A

  # Override auto-derived FASTA if using non-standard layout
  python create_config.py \\
    --output-dir /path/to/output \\
    --sample-pickle samples.pkl \\
    --cellranger-reference /path/to/cellranger_ref \\
    --reference-fasta /custom/path/to/genome.fa

  # Enable all modules with custom popV model
  python create_config.py \\
    --output-dir /path/to/output \\
    --sample-ids SAMPLE1 SAMPLE2 SAMPLE3 \\
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \\
    --popv-huggingface-repo "popV/tabula_sapiens_immune" \\
    --enable-viral --kraken-db /path/to/kraken_db \\
    --enable-scomatic \\
      --scomatic-scripts-dir /path/to/SComatic/scripts \\
      --scomatic-editing-sites /path/to/editing_sites.txt \\
      --scomatic-pon-file /path/to/pon.tsv \\
      --scomatic-bed-file /path/to/mappable.bed \\
    --enable-signatures \\
      --cosmic-file /path/to/COSMIC_v3.4_SBS_GRCh38.txt \\
      --core-signatures SBS2 SBS13 SBS5 \\
      --use-scree-plot
        """)
    
    # ==========================================================================
    # Required arguments
    # ==========================================================================
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('--output-dir', required=True, help='Output directory')
    required.add_argument('--cellranger-reference', required=True, 
                          help='Cell Ranger reference directory (contains fasta/, genes/, star/)')
    
    # ==========================================================================
    # Reference files (optional - auto-derived from cellranger-reference)
    # ==========================================================================
    refs = parser.add_argument_group('Reference Files (Optional - Auto-derived from Cell Ranger reference)')
    refs.add_argument('--reference-fasta', 
                      help='Reference FASTA file (auto-derived as {cellranger}/fasta/genome.fa if not provided)')
    refs.add_argument('--gtf-file', 
                      help='GTF annotation file (auto-derived as {cellranger}/genes/genes.gtf if not provided)')
    refs.add_argument('--genome', default='GRCh38', choices=['GRCh38', 'GRCh37', 'mm10', 'mm39'],
                      help='Genome build (default: GRCh38)')
    
    # ==========================================================================
    # Sample specification
    # ==========================================================================
    samples = parser.add_argument_group('Sample Specification (one required)')
    samples.add_argument('--sample-pickle', help='Pickle file with sample dictionary')
    samples.add_argument('--sample-ids', nargs='+', help='List of sample IDs')
    
    # ==========================================================================
    # Resource settings
    # ==========================================================================
    resources = parser.add_argument_group('Resource Settings')
    resources.add_argument('--threads', type=int, default=8, help='Number of threads (default: 8)')
    resources.add_argument('--memory-gb', type=int, default=64, help='Memory in GB (default: 64)')
    
    # ==========================================================================
    # Cell Ranger settings
    # ==========================================================================
    cellranger = parser.add_argument_group('Cell Ranger Settings')
    cellranger.add_argument('--chemistry', default='auto', help='Chemistry (default: auto)')
    cellranger.add_argument('--expect-cells', type=int, default=10000, help='Expected cells (default: 10000)')
    cellranger.add_argument('--include-introns', action='store_true', help='Include introns')
    
    # ==========================================================================
    # QC settings
    # ==========================================================================
    qc = parser.add_argument_group('QC Settings')
    qc.add_argument('--min-genes', type=int, default=200, help='Min genes per cell (default: 200)')
    qc.add_argument('--max-genes', type=int, default=5000, help='Max genes per cell (default: 5000)')
    qc.add_argument('--min-counts', type=int, default=500, help='Min counts per cell (default: 500)')
    qc.add_argument('--max-counts', type=int, default=50000, help='Max counts per cell (default: 50000)')
    qc.add_argument('--max-mito-pct', type=float, default=20, help='Max mitochondrial %% (default: 20)')
    qc.add_argument('--min-cells', type=int, default=3, help='Min cells expressing a gene (default: 3)')
    qc.add_argument('--doublet-removal', action='store_true', default=True, help='Enable doublet removal (default: True)')
    qc.add_argument('--no-doublet-removal', action='store_false', dest='doublet_removal', help='Disable doublet removal')
    qc.add_argument('--doublet-rate', type=float, default=0.08, help='Expected doublet rate (default: 0.08)')
    
    # ==========================================================================
    # Preprocessing settings (post-annotation)
    # ==========================================================================
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
    
    # ==========================================================================
    # Annotation settings (popV only)
    # ==========================================================================
    annotation = parser.add_argument_group('Annotation Settings (popV)')
    annotation.add_argument('--popv-huggingface-repo', default='popV/tabula_sapiens_All_Cells',
                           help='Hugging Face repository for popV model (default: popV/tabula_sapiens_All_Cells)')
    annotation.add_argument('--popv-prediction-mode', default='inference', choices=['inference', 'fast'],
                           help='popV prediction mode: inference (full) or fast (default: inference)')
    annotation.add_argument('--popv-gene-symbol-key', default='feature_name',
                           help='Gene symbol column in adata.var (default: feature_name)')
    annotation.add_argument('--popv-cache-dir', default='tmp/popv_models',
                           help='Cache directory for popV models (default: tmp/popv_models)')
    annotation.add_argument('--batch-key', default='sample_id',
                           help='Batch key for annotation (default: sample_id)')
    
    # ==========================================================================
    # Dysregulation settings
    # ==========================================================================
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
    
    # ==========================================================================
    # Viral detection
    # ==========================================================================
    viral = parser.add_argument_group('Viral Detection')
    viral.add_argument('--enable-viral', action='store_true', help='Enable viral detection')
    viral.add_argument('--kraken-db', help='Kraken2 database path')
    viral.add_argument('--viral-db', help='Viral database path')
    viral.add_argument('--viral-confidence', type=float, default=0.1, 
                       help='Viral detection confidence (default: 0.1)')
    
    # ==========================================================================
    # SComatic settings
    # ==========================================================================
    scomatic = parser.add_argument_group('SComatic Mutation Calling')
    scomatic.add_argument('--enable-scomatic', action='store_true', help='Enable SComatic mutation calling')
    scomatic.add_argument('--scomatic-scripts-dir', help='SComatic scripts directory')
    scomatic.add_argument('--scomatic-editing-sites', help='RNA editing sites file')
    scomatic.add_argument('--scomatic-pon-file', help='Panel of Normals file')
    scomatic.add_argument('--scomatic-bed-file', help='Mappable regions BED file')
    scomatic.add_argument('--scomatic-min-cov', type=int, default=5, help='Min coverage (default: 5)')
    scomatic.add_argument('--scomatic-min-cells', type=int, default=5, help='Min cells with variant (default: 5)')
    scomatic.add_argument('--scomatic-min-base-quality', type=int, default=30, 
                          help='Min base quality (default: 30)')
    scomatic.add_argument('--scomatic-min-map-quality', type=int, default=30, 
                          help='Min mapping quality (default: 30)')
    
    # ==========================================================================
    # Signature analysis settings
    # ==========================================================================
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


if __name__ == '__main__':
    main()
