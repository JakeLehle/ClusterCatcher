#!/usr/bin/env python3
"""
create-config command
=====================

Generate the master configuration YAML file for the ClusterCatcher pipeline.

This command takes sample information (either from sample-information command
or directly from SRAscraper output) and user-specified parameters to create
a comprehensive configuration file for the Snakemake pipeline.
"""

import click
import os
import sys
import pickle
import yaml
from pathlib import Path
from datetime import datetime


def get_default_config():
    """
    Return default configuration values.
    
    Returns
    -------
    dict
        Default configuration dictionary
    """
    return {
        # Pipeline info
        'pipeline_version': '0.1.0',
        'created': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        
        # Output settings
        'output_dir': './results',
        'temp_dir': './tmp',
        'log_dir': './logs',
        
        # Computing resources
        'threads': 8,
        'memory_gb': 32,
        
        # Reference data
        'reference': {
            'genome': 'GRCh38',
            'transcriptome': None,  # Path to Cell Ranger transcriptome reference (REQUIRED)
            'fasta': None,          # Path to reference FASTA (for SComatic)
            'gtf': None,            # Path to GTF annotation file (for inferCNV)
        },
        
        # Cell Ranger settings
        # Note: Cell Ranger must be installed separately and available in PATH
        # See: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
        'cellranger': {
            'chemistry': 'auto',      # auto, SC3Pv2, SC3Pv3, SC3Pv3HT, SC3Pv4, threeprime, fiveprime, etc.
            'expect_cells': None,     # None for auto-detection
            'force_cells': None,      # Force specific cell count
            'include_introns': True,  # Include intronic reads (recommended for most analyses)
            'create_bam': True,       # Create BAM file (needed for mutation calling)
            'localcores': 8,          # Cores per Cell Ranger job
            'localmem': 32,           # Memory (GB) per Cell Ranger job
            'no_bam': False,          # Skip BAM creation entirely
        },
        
        # QC and filtering (Scanpy)
        'qc': {
            'min_genes': 200,
            'min_cells': 3,
            'max_genes': 5000,
            'min_counts': 500,
            'max_counts': 50000,
            'max_pct_mito': 20,
            'doublet_removal': True,
            'doublet_rate': 0.06,
        },
        
        # Preprocessing parameters
        'preprocessing': {
            'target_sum': 1e4,       # Normalization target
            'n_top_genes': 2000,     # Number of HVGs
            'n_pcs': 50,             # PCA components
            'n_neighbors': 15,       # For neighbor graph
            'leiden_resolution': 1.0,
        },
        
        # Cell annotation
        'annotation': {
            'method': 'popv',  # popv or leiden
            'popv_model': 'popv_immune_All_human_umap_from_cellxgene',  # HuggingFace model
            'batch_key': 'sample_id',
            'min_confidence': 0.5,
        },
        
        # Dysregulation detection (Cancer cell identification)
        # Uses CytoTRACE2 (stemness) and inferCNV (copy number variation)
        'dysregulation': {
            'enabled': True,
            'cytotrace2': {
                'species': 'human',           # human or mouse
                'max_cells_per_chunk': 200000,  # Chunk size for large datasets
                'seed': 42,
            },
            'infercnv': {
                'window_size': 250,           # Genomic window size
                'step': 10,
                'dynamic_threshold': 1.5,
            },
            'agreement': {
                'alpha': 0.5,                 # Weight: 0=value only, 1=rank only, 0.5=both
                'min_correlation': 0.5,       # Min Spearman rho for quartile selection
            },
        },
        
        # Viral detection (Kraken2)
        # Requires: Kraken2 database setup (see README for instructions)
        'viral_detection': {
            'enabled': True,
            'kraken2_db': None,       # Path to Kraken2 database (REQUIRED if enabled)
            'human_viral_db': None,   # Path to human-specific viral Kraken2 inspect.txt (for filtering)
            'confidence': 0.0,        # Kraken2 confidence threshold (0.0-1.0)
            'include_organisms': None,  # List of organism patterns to include (None = all)
            'exclude_organisms': None,  # List of organism patterns to exclude
            'organisms_of_interest': [  # Organisms to highlight in reports
                'Human papillomavirus',
                'Epstein-Barr virus',
                'Hepatitis B virus',
                'Hepatitis C virus',
                'Human T-lymphotropic virus',
                'Human immunodeficiency virus',
                'Merkel cell polyomavirus',
                'Kaposi sarcoma-associated herpesvirus',
            ],
        },
        
        # Somatic mutation calling
        'scomatic': {
            'enabled': True,
            'min_cells': 5,
            'min_cov': 5,
            'min_cc': 5,
            'min_ar': 0.1,
            # reference_fasta comes from reference.fasta
        },
        
        # Mutational signature analysis
        'signatures': {
            'enabled': True,
            'method': 'semi_supervised_nmf',
            'cosmic_version': '3.4',
            'relevant_signatures': [
                'SBS1', 'SBS2', 'SBS5', 'SBS13', 'SBS18',
                'SBS29', 'SBS40', 'SBS44', 'SBS88', 'SBS89',
            ],
            'min_mutations': 10,
            'similarity_metrics': ['cosine', 'euclidean', 'dot_product'],
        },
        
        # Target files for the pipeline
        'targets': [
            'cellranger_counts',
            'qc_report',
            'cell_annotations',
            'dysregulation_scores',
            'viral_detection_report',
            'somatic_mutations',
            'signature_profiles',
            'master_summary',
        ],
    }


def validate_config(config):
    """
    Validate configuration values.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary to validate
        
    Returns
    -------
    tuple
        (is_valid, error_messages, warnings)
    """
    errors = []
    warnings = []
    
    # Check required paths
    if config['reference']['transcriptome'] is None:
        errors.append("Cell Ranger transcriptome reference path is required (--transcriptome)")
        errors.append("  Download from: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest")
        
    if config['scomatic']['enabled'] and config['reference'].get('fasta') is None:
        errors.append("SComatic requires reference FASTA (--reference-fasta)")
        
    if config['viral_detection']['enabled'] and config['viral_detection']['kraken2_db'] is None:
        errors.append("Viral detection requires Kraken2 database path (--kraken2-db)")
        
    # Check sample information
    if 'samples' not in config or not config['samples']:
        errors.append("No samples provided")
        
    # Validate sample FASTQ paths
    samples = config.get('samples', {})
    for sample_id, sample_info in samples.items():
        fastq_r1 = sample_info.get('fastq_r1', [])
        if not fastq_r1:
            errors.append(f"Sample {sample_id}: No FASTQ R1 files specified")
            
    # Warnings for suboptimal settings
    if config['qc']['max_pct_mito'] > 30:
        warnings.append(f"High mitochondrial percentage threshold ({config['qc']['max_pct_mito']}%) may include low-quality cells")
        
    if config['threads'] > 32:
        warnings.append(f"High thread count ({config['threads']}) - ensure your system can handle this")
        
    if config['cellranger']['chemistry'] != 'auto':
        warnings.append(f"Using fixed chemistry '{config['cellranger']['chemistry']}' - will fallback to auto if this fails")
        
    if not config['cellranger']['create_bam'] and config['scomatic']['enabled']:
        errors.append("BAM creation is disabled but SComatic (mutation calling) is enabled - BAM files are required")
        
    return len(errors) == 0, errors, warnings


def merge_configs(default, user):
    """
    Recursively merge user config into default config.
    
    Parameters
    ----------
    default : dict
        Default configuration
    user : dict
        User-specified configuration
        
    Returns
    -------
    dict
        Merged configuration
    """
    result = default.copy()
    for key, value in user.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = merge_configs(result[key], value)
        else:
            result[key] = value
    return result


@click.command('create-config')
@click.option(
    '--samples', '-s', 'samples_pkl',
    required=True,
    type=click.Path(exists=True),
    help='Path to sample dictionary pickle (from sample-information or SRAscraper)'
)
@click.option(
    '--output', '-o', 'output_yaml',
    required=True,
    type=click.Path(),
    help='Path for output config YAML file'
)
@click.option(
    '--results-dir', '-r',
    default='./results',
    type=click.Path(),
    help='Directory for pipeline output (default: ./results)'
)
@click.option(
    '--threads', '-t',
    default=8,
    type=int,
    help='Number of threads per job (default: 8)'
)
@click.option(
    '--memory', '-m',
    default=32,
    type=int,
    help='Memory in GB per job (default: 32)'
)
@click.option(
    '--transcriptome',
    type=click.Path(),
    help='Path to Cell Ranger transcriptome reference'
)
@click.option(
    '--genome',
    default='GRCh38',
    type=click.Choice(['GRCh38', 'GRCh37', 'mm10', 'mm39']),
    help='Reference genome build (default: GRCh38)'
)
@click.option(
    '--reference-fasta',
    type=click.Path(),
    help='Path to reference FASTA for SComatic'
)
@click.option(
    '--gtf',
    type=click.Path(),
    help='Path to GTF annotation file (for inferCNV chromosomal locations)'
)
@click.option(
    '--kraken2-db',
    type=click.Path(),
    help='Path to Kraken2 database for viral detection'
)
@click.option(
    '--human-viral-db',
    type=click.Path(),
    help='Path to human-specific viral Kraken2 database inspect.txt (for filtering)'
)
@click.option(
    '--chemistry',
    default='auto',
    help='10X chemistry version (default: auto)'
)
@click.option(
    '--expect-cells',
    type=int,
    help='Expected number of cells per sample'
)
@click.option(
    '--force-cells',
    type=int,
    help='Force this exact number of cells (use with caution)'
)
@click.option(
    '--no-introns',
    is_flag=True,
    default=False,
    help='Exclude intronic reads from Cell Ranger counting'
)
@click.option(
    '--no-bam',
    is_flag=True,
    default=False,
    help='Skip BAM file creation (disables mutation calling)'
)
@click.option(
    '--no-viral-detection',
    is_flag=True,
    default=False,
    help='Disable viral detection step'
)
@click.option(
    '--kraken2-confidence',
    default=0.0,
    type=float,
    help='Kraken2 confidence threshold (0.0-1.0, default: 0.0)'
)
@click.option(
    '--include-organisms',
    type=str,
    help='Comma-separated organism patterns to include in viral detection'
)
@click.option(
    '--exclude-organisms',
    type=str,
    help='Comma-separated organism patterns to exclude from viral detection'
)
@click.option(
    '--no-scomatic',
    is_flag=True,
    default=False,
    help='Disable somatic mutation calling'
)
@click.option(
    '--no-cancer-detection',
    is_flag=True,
    default=False,
    help='Disable cancer cell detection (CytoTRACE2 + inferCNV)'
)
@click.option(
    '--cytotrace-species',
    default='human',
    type=click.Choice(['human', 'mouse']),
    help='Species for CytoTRACE2 (default: human)'
)
@click.option(
    '--cytotrace-chunk-size',
    default=200000,
    type=int,
    help='Max cells per CytoTRACE2 chunk (default: 200000)'
)
@click.option(
    '--agreement-alpha',
    default=0.5,
    type=float,
    help='Weight for rank vs value agreement (0-1, default: 0.5)'
)
@click.option(
    '--no-signatures',
    is_flag=True,
    default=False,
    help='Disable mutational signature analysis'
)
@click.option(
    '--cosmic-signatures',
    default='SBS1,SBS2,SBS5,SBS13,SBS18,SBS29,SBS40,SBS44,SBS88,SBS89',
    help='Comma-separated list of COSMIC signatures to use'
)
@click.option(
    '--annotation-method',
    default='popv',
    type=click.Choice(['popv', 'leiden']),
    help='Cell type annotation method (default: popv)'
)
@click.option(
    '--popv-model',
    default='popv_immune_All_human_umap_from_cellxgene',
    help='popV model from HuggingFace (default: popv_immune_All_human_umap_from_cellxgene)'
)
@click.option(
    '--leiden-resolution',
    default=1.0,
    type=float,
    help='Resolution for Leiden clustering (default: 1.0)'
)
@click.option(
    '--n-hvgs',
    default=2000,
    type=int,
    help='Number of highly variable genes for preprocessing (default: 2000)'
)
@click.option(
    '--min-genes',
    default=200,
    type=int,
    help='Minimum genes per cell for QC (default: 200)'
)
@click.option(
    '--max-mito',
    default=20,
    type=float,
    help='Maximum mitochondrial percentage for QC (default: 20)'
)
@click.option(
    '--force',
    is_flag=True,
    default=False,
    help='Overwrite existing config file'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    default=False,
    help='Print verbose output'
)
def create_config(samples_pkl, output_yaml, results_dir, threads, memory,
                  transcriptome, genome, reference_fasta, gtf, kraken2_db, human_viral_db,
                  chemistry, expect_cells, force_cells, no_introns, no_bam,
                  no_viral_detection, kraken2_confidence, include_organisms, exclude_organisms,
                  no_scomatic, no_cancer_detection, cytotrace_species, cytotrace_chunk_size,
                  agreement_alpha,
                  no_signatures, cosmic_signatures, annotation_method,
                  popv_model, leiden_resolution, n_hvgs,
                  min_genes, max_mito, force, verbose):
    """
    Generate master configuration YAML for the pipeline.
    
    This command creates a comprehensive configuration file that controls
    all aspects of the ClusterCatcher pipeline. The configuration can be
    manually edited after creation if needed.
    
    \b
    The pipeline steps controlled by this config:
      1. Cell Ranger alignment and counting
      2. Scanpy QC and filtering
      3. Cell type annotation (popV/CellTypist)
      4. Dysregulation detection (CytoTRACE2 + inferCNV)
      5. Viral detection in unmapped reads (Kraken2)
      6. Somatic mutation calling (SComatic)
      7. Mutational signature deconvolution (semi-supervised NMF)
    
    \b
    Example:
      ClusterCatcher create-config \\
        --samples samples.pkl \\
        --output config.yaml \\
        --transcriptome /path/to/refdata-gex-GRCh38-2020-A \\
        --reference-fasta /path/to/GRCh38.fa \\
        --kraken2-db /path/to/k2_viral_db \\
        --threads 16 \\
        --memory 64
    """
    
    click.echo(f"\n{'='*60}")
    click.echo("ClusterCatcher: create-config")
    click.echo(f"{'='*60}\n")
    
    # Check if output exists
    if os.path.exists(output_yaml) and not force:
        click.echo(f"ERROR: Output file exists: {output_yaml}", err=True)
        click.echo("Use --force to overwrite", err=True)
        sys.exit(1)
        
    # Load sample information
    click.echo(f"Loading sample information from: {samples_pkl}")
    try:
        with open(samples_pkl, 'rb') as f:
            samples = pickle.load(f)
    except Exception as e:
        click.echo(f"ERROR: Failed to load sample pickle: {e}", err=True)
        sys.exit(1)
        
    click.echo(f"Found {len(samples)} samples")
    
    if verbose:
        click.echo("\nSamples:")
        for sample_id in samples:
            click.echo(f"  - {sample_id}")
            
    # Build configuration
    click.echo("\nBuilding configuration...")
    config = get_default_config()
    
    # Update with user options
    config['output_dir'] = results_dir
    config['threads'] = threads
    config['memory_gb'] = memory
    
    config['reference']['genome'] = genome
    config['reference']['transcriptome'] = transcriptome
    config['reference']['fasta'] = reference_fasta
    config['reference']['gtf'] = gtf
    
    config['cellranger']['chemistry'] = chemistry
    config['cellranger']['localcores'] = threads
    config['cellranger']['localmem'] = memory
    config['cellranger']['include_introns'] = not no_introns
    config['cellranger']['create_bam'] = not no_bam
    config['cellranger']['no_bam'] = no_bam
    if expect_cells:
        config['cellranger']['expect_cells'] = expect_cells
    if force_cells:
        config['cellranger']['force_cells'] = force_cells
        
    config['qc']['min_genes'] = min_genes
    config['qc']['max_pct_mito'] = max_mito
    
    config['preprocessing']['n_top_genes'] = n_hvgs
    config['preprocessing']['leiden_resolution'] = leiden_resolution
    
    config['annotation']['method'] = annotation_method
    config['annotation']['popv_model'] = popv_model
    
    config['viral_detection']['enabled'] = not no_viral_detection
    config['viral_detection']['kraken2_db'] = kraken2_db
    config['viral_detection']['human_viral_db'] = human_viral_db
    config['viral_detection']['confidence'] = kraken2_confidence
    if include_organisms:
        config['viral_detection']['include_organisms'] = include_organisms.split(',')
    if exclude_organisms:
        config['viral_detection']['exclude_organisms'] = exclude_organisms.split(',')
    
    config['scomatic']['enabled'] = not no_scomatic
    config['scomatic']['reference_fasta'] = reference_fasta
    
    config['dysregulation']['enabled'] = not no_cancer_detection
    config['dysregulation']['cytotrace2']['species'] = cytotrace_species
    config['dysregulation']['cytotrace2']['max_cells_per_chunk'] = cytotrace_chunk_size
    config['dysregulation']['agreement']['alpha'] = agreement_alpha
    
    config['signatures']['enabled'] = not no_signatures
    config['signatures']['relevant_signatures'] = cosmic_signatures.split(',')
    
    # Add sample information
    config['samples'] = samples
    
    # Determine targets based on enabled steps
    targets = ['cellranger_counts', 'qc_report', 'cell_annotations', 'dysregulation_scores']
    if not no_viral_detection:
        targets.append('viral_detection_report')
    if not no_scomatic:
        targets.append('somatic_mutations')
    if not no_signatures:
        targets.append('signature_profiles')
    targets.append('master_summary')
    config['targets'] = targets
    
    # Validate configuration
    click.echo("\nValidating configuration...")
    is_valid, errors, warnings = validate_config(config)
    
    if warnings:
        click.echo("\nWarnings:")
        for w in warnings:
            click.echo(f"  ⚠ {w}")
            
    if not is_valid:
        click.echo("\nERRORS:", err=True)
        for e in errors:
            click.echo(f"  ✗ {e}", err=True)
        sys.exit(1)
        
    click.echo("  ✓ Configuration validated")
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_yaml)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    # Write configuration
    click.echo(f"\nWriting configuration to: {output_yaml}")
    
    # Custom YAML representer for cleaner output
    def str_representer(dumper, data):
        if '\n' in data:
            return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
        return dumper.represent_scalar('tag:yaml.org,2002:str', data)
    
    yaml.add_representer(str, str_representer)
    
    with open(output_yaml, 'w') as f:
        f.write("# ClusterCatcher Pipeline Configuration\n")
        f.write(f"# Generated: {config['created']}\n")
        f.write(f"# Pipeline version: {config['pipeline_version']}\n")
        f.write("#\n")
        f.write("# Edit this file to customize pipeline parameters.\n")
        f.write("# See documentation for detailed parameter descriptions.\n")
        f.write("#\n\n")
        yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
        
    click.echo(f"\n{'='*60}")
    click.echo("SUCCESS: Configuration created")
    click.echo(f"{'='*60}")
    click.echo(f"\nOutput: {output_yaml}")
    click.echo(f"Samples: {len(samples)}")
    click.echo(f"Targets: {len(targets)}")
    click.echo("\nEnabled steps:")
    click.echo(f"  ✓ Cell Ranger counting")
    click.echo(f"  ✓ Scanpy QC")
    click.echo(f"  ✓ Cell annotation ({annotation_method})")
    click.echo(f"  ✓ Dysregulation detection (CytoTRACE2 + inferCNV)")
    click.echo(f"  {'✓' if not no_viral_detection else '✗'} Viral detection")
    click.echo(f"  {'✓' if not no_scomatic else '✗'} Somatic mutation calling")
    click.echo(f"  {'✓' if not no_signatures else '✗'} Mutational signatures")
    click.echo("\nNext step:")
    click.echo(f"  ClusterCatcher run-config {output_yaml}")
    click.echo("")


if __name__ == '__main__':
    create_config()
