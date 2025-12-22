#!/usr/bin/env python3
"""
sample-information command
==========================

Process a CSV file containing sample information and create a Python dictionary
pickle file that can be used by create-config.

This command is intended for users who generated their own single-cell FASTQ files
and did not use SRAscraper to download from GEO.

Expected CSV columns:
- sample_id: Unique identifier for each sample
- fastq_r1: Path to R1 FASTQ file (can be comma-separated for multiple lanes)
- fastq_r2: Path to R2 FASTQ file (can be comma-separated for multiple lanes)
- sample_name: Human-readable sample name (optional)
- condition: Experimental condition/group (optional)
- donor: Donor/patient ID for multi-sample studies (optional)

Additional optional columns will be preserved in the metadata.
"""

import click
import os
import sys
import pickle
import pandas as pd
from pathlib import Path


def validate_fastq_paths(df, check_existence=True):
    """
    Validate that FASTQ paths exist and are properly formatted.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with fastq_r1 and fastq_r2 columns
    check_existence : bool
        Whether to check if files actually exist
        
    Returns
    -------
    tuple
        (is_valid, error_messages)
    """
    errors = []
    
    for idx, row in df.iterrows():
        sample_id = row['sample_id']
        
        # Handle comma-separated paths for multiple lanes
        r1_paths = [p.strip() for p in str(row['fastq_r1']).split(',') if p.strip()]
        r2_paths = [p.strip() for p in str(row['fastq_r2']).split(',') if p.strip()]
        
        # Check for empty paths
        if not r1_paths:
            errors.append(f"Sample {sample_id}: No R1 paths provided")
            continue
        if not r2_paths:
            errors.append(f"Sample {sample_id}: No R2 paths provided")
            continue
        
        # Check matching number of R1 and R2 files
        if len(r1_paths) != len(r2_paths):
            errors.append(
                f"Sample {sample_id}: Mismatched number of R1 ({len(r1_paths)}) "
                f"and R2 ({len(r2_paths)}) files"
            )
            continue
        
        if check_existence:
            for r1, r2 in zip(r1_paths, r2_paths):
                if not os.path.exists(r1):
                    errors.append(f"Sample {sample_id}: R1 file not found: {r1}")
                if not os.path.exists(r2):
                    errors.append(f"Sample {sample_id}: R2 file not found: {r2}")
                    
    return len(errors) == 0, errors


def create_sample_dict(df):
    """
    Convert DataFrame to nested dictionary structure.
    
    The structure matches what the Snakefile and create_config.py expect:
    {
        'sample_id': {
            'fastq_r1': [list of R1 paths],
            'fastq_r2': [list of R2 paths],
            'lanes': number of lanes,
            ... additional metadata
        }
    }
    
    Parameters
    ----------
    df : pd.DataFrame
        Validated sample information DataFrame
        
    Returns
    -------
    dict
        Nested dictionary with sample information
    """
    sample_dict = {}
    
    for idx, row in df.iterrows():
        sample_id = str(row['sample_id'])
        
        # Parse FASTQ paths
        r1_paths = [p.strip() for p in str(row['fastq_r1']).split(',') if p.strip()]
        r2_paths = [p.strip() for p in str(row['fastq_r2']).split(',') if p.strip()]
        
        # Build sample entry
        sample_entry = {
            'fastq_r1': r1_paths,
            'fastq_r2': r2_paths,
            'R1': r1_paths,  # Alternative key for compatibility
            'R2': r2_paths,  # Alternative key for compatibility
            'lanes': len(r1_paths),
        }
        
        # Add optional fields
        optional_fields = ['sample_name', 'condition', 'donor', 'tissue', 
                          'species', 'chemistry', 'notes', 'fastq_dir', 'path']
        for field in optional_fields:
            if field in row and pd.notna(row[field]):
                sample_entry[field] = row[field]
        
        # Add fastq_dir if not present (derive from first R1 path)
        if 'fastq_dir' not in sample_entry and r1_paths:
            sample_entry['fastq_dir'] = os.path.dirname(r1_paths[0])
            sample_entry['path'] = os.path.dirname(r1_paths[0])
        
        # Add any additional columns as metadata
        standard_cols = ['sample_id', 'fastq_r1', 'fastq_r2'] + optional_fields
        for col in df.columns:
            if col not in standard_cols and pd.notna(row[col]):
                sample_entry[col] = row[col]
        
        sample_dict[sample_id] = sample_entry
    
    return sample_dict


def create_srascraper_compatible_dict(df):
    """
    Create a dictionary structure compatible with SRAscraper output.
    
    This creates a nested structure: {series_id: DataFrame}
    which is what some parts of the pipeline expect from SRAscraper.
    
    Parameters
    ----------
    df : pd.DataFrame
        Sample information DataFrame
        
    Returns
    -------
    dict
        SRAscraper-compatible dictionary
    """
    # If there's a series_id column, group by it
    if 'series_id' in df.columns:
        grouped = {}
        for series_id, group_df in df.groupby('series_id'):
            # Rename columns to match SRAscraper format
            sra_df = group_df.copy()
            if 'sample_id' in sra_df.columns:
                sra_df = sra_df.rename(columns={'sample_id': 'run_accession'})
            grouped[series_id] = sra_df
        return grouped
    else:
        # Create a single group with a generic series ID
        sra_df = df.copy()
        if 'sample_id' in sra_df.columns:
            sra_df = sra_df.rename(columns={'sample_id': 'run_accession'})
        return {'LOCAL': sra_df}


@click.command('sample-information')
@click.option(
    '--input', '-i', 'input_csv',
    required=True,
    type=click.Path(exists=True),
    help='Path to CSV file containing sample information'
)
@click.option(
    '--output', '-o', 'output_pkl',
    required=True,
    type=click.Path(),
    help='Path for output pickle file (sample dictionary)'
)
@click.option(
    '--skip-validation', '-s',
    is_flag=True,
    default=False,
    help='Skip validation of FASTQ file existence (useful for cluster environments)'
)
@click.option(
    '--srascraper-format',
    is_flag=True,
    default=False,
    help='Output in SRAscraper-compatible format (nested dict with DataFrames)'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    default=False,
    help='Print verbose output'
)
def sample_information(input_csv, output_pkl, skip_validation, srascraper_format, verbose):
    """
    Process sample CSV and create sample dictionary pickle.
    
    This command reads a CSV file containing information about your single-cell
    FASTQ files and creates a standardized pickle file that can be used with
    the create-config command.
    
    \b
    Required CSV columns:
      - sample_id: Unique identifier for each sample
      - fastq_r1: Path(s) to R1 FASTQ file(s) - comma-separated for multiple lanes
      - fastq_r2: Path(s) to R2 FASTQ file(s) - comma-separated for multiple lanes
    
    \b
    Optional CSV columns:
      - sample_name: Human-readable sample name
      - condition: Experimental condition/group
      - donor: Donor/patient ID
      - tissue: Tissue type
      - species: Species (default: human)
      - chemistry: 10X chemistry version (e.g., SC3Pv3)
      - series_id: GEO series ID (for grouping)
      - notes: Additional notes
    
    \b
    Example CSV format:
      sample_id,fastq_r1,fastq_r2,condition,donor
      S1,/path/S1_R1.fq.gz,/path/S1_R2.fq.gz,tumor,P001
      S2,/path/S2_R1.fq.gz,/path/S2_R2.fq.gz,normal,P001
    
    \b
    For multi-lane samples (comma-separated paths):
      sample_id,fastq_r1,fastq_r2
      S1,/path/S1_L001_R1.fq.gz;/path/S1_L002_R1.fq.gz,/path/S1_L001_R2.fq.gz;/path/S1_L002_R2.fq.gz
    
    \b
    Usage examples:
    
      # Basic usage
      ClusterCatcher sample-information -i samples.csv -o samples.pkl
    
      # Skip file existence check (for cluster)
      ClusterCatcher sample-information -i samples.csv -o samples.pkl --skip-validation
    
      # Create SRAscraper-compatible format
      ClusterCatcher sample-information -i samples.csv -o samples.pkl --srascraper-format
    """
    
    click.echo(f"\n{'='*60}")
    click.echo("ClusterCatcher: sample-information")
    click.echo(f"{'='*60}\n")
    
    # Read CSV
    click.echo(f"Reading sample information from: {input_csv}")
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        click.echo(f"ERROR: Failed to read CSV file: {e}", err=True)
        sys.exit(1)
    
    # Check required columns
    required_cols = ['sample_id', 'fastq_r1', 'fastq_r2']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        click.echo(f"ERROR: Missing required columns: {missing_cols}", err=True)
        click.echo(f"Found columns: {list(df.columns)}", err=True)
        click.echo("\nExpected CSV format:", err=True)
        click.echo("  sample_id,fastq_r1,fastq_r2[,optional columns...]", err=True)
        sys.exit(1)
    
    # Check for duplicate sample IDs
    if df['sample_id'].duplicated().any():
        dups = df[df['sample_id'].duplicated()]['sample_id'].tolist()
        click.echo(f"ERROR: Duplicate sample IDs found: {dups}", err=True)
        sys.exit(1)
    
    click.echo(f"Found {len(df)} samples")
    
    if verbose:
        click.echo("\nSample IDs:")
        for sid in df['sample_id']:
            click.echo(f"  - {sid}")
        click.echo(f"\nColumns found: {list(df.columns)}")
    
    # Validate FASTQ paths
    if not skip_validation:
        click.echo("\nValidating FASTQ paths...")
        is_valid, errors = validate_fastq_paths(df, check_existence=True)
        if not is_valid:
            click.echo("ERROR: FASTQ validation failed:", err=True)
            for err in errors[:10]:  # Show first 10 errors
                click.echo(f"  - {err}", err=True)
            if len(errors) > 10:
                click.echo(f"  ... and {len(errors) - 10} more errors", err=True)
            click.echo("\nUse --skip-validation to skip file existence checks", err=True)
            sys.exit(1)
        click.echo("  ✓ All FASTQ paths validated successfully")
    else:
        click.echo("\nSkipping FASTQ path validation (--skip-validation)")
        # Still check format
        is_valid, errors = validate_fastq_paths(df, check_existence=False)
        if not is_valid:
            click.echo("ERROR: FASTQ path format validation failed:", err=True)
            for err in errors:
                click.echo(f"  - {err}", err=True)
            sys.exit(1)
    
    # Create sample dictionary
    click.echo("\nCreating sample dictionary...")
    if srascraper_format:
        sample_dict = create_srascraper_compatible_dict(df)
        click.echo("  Using SRAscraper-compatible format")
    else:
        sample_dict = create_sample_dict(df)
        click.echo("  Using standard format")
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_pkl)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Save pickle
    click.echo(f"Saving to: {output_pkl}")
    with open(output_pkl, 'wb') as f:
        pickle.dump(sample_dict, f)
    
    click.echo(f"\n{'='*60}")
    click.echo("SUCCESS: Sample information processed")
    click.echo(f"{'='*60}")
    click.echo(f"\nOutput: {output_pkl}")
    click.echo(f"Samples: {len(df)}")
    
    if verbose:
        click.echo("\nSample dictionary structure:")
        for sid in list(sample_dict.keys())[:3]:
            if isinstance(sample_dict[sid], dict):
                click.echo(f"  {sid}:")
                for k, v in list(sample_dict[sid].items())[:4]:
                    click.echo(f"    {k}: {v}")
    
    click.echo("\n" + "="*60)
    click.echo("Next step: Create pipeline configuration")
    click.echo("="*60)
    click.echo(f"\n  python snakemake_wrapper/create_config.py \\")
    click.echo(f"    --output-dir ./results \\")
    click.echo(f"    --sample-pickle {output_pkl} \\")
    click.echo(f"    --reference-fasta /path/to/GRCh38.fa \\")
    click.echo(f"    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A")
    click.echo("")


if __name__ == '__main__':
    sample_information()
