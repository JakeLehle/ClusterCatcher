#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cellranger_count.py
===================

Run Cell Ranger count on single-cell FASTQ files with automatic chemistry detection.

This script is adapted from the standalone Cell Ranger processing script to work
within the ClusterCatcher Snakemake pipeline. It handles:
- Automatic chemistry detection with intelligent fallbacks
- Multiplexed sample detection (I1/I2 files)
- Failed run cleanup and retry logic
- Tracking of successful runs

Usage:
    Called via Snakemake rule with snakemake.input/output/params
    
    Or standalone:
    python cellranger_count.py --config config.yaml --sample SAMPLE_ID

Requirements:
    - Cell Ranger must be installed and available in PATH
    - Reference transcriptome must be downloaded
"""

import os
import sys
import yaml
import pickle
import shutil
import subprocess
import traceback
import logging
import argparse
from pathlib import Path
from os import cpu_count

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Chemistry Options
# =============================================================================

# List of all possible chemistry options to try (in order of likelihood)
CHEMISTRY_OPTIONS = [
    "auto",           # Try auto first
    "threeprime",
    "fiveprime", 
    "SC3Pv4-polyA",
    "SC3Pv4",
    "SC3Pv3-polyA",
    "SC3Pv3",
    "SC3Pv2",
    "SC3Pv3HT",
    "SC3Pv3HT-polyA",
    "SC5P-PE-v3",
    "SC5P-PE",
    "SC5P-R2-v3",
    "SC5P-R2",
    "SC3Pv1",         # These last two can't be auto-detected
    "ARC-v1"
]


# =============================================================================
# Helper Functions
# =============================================================================

def detect_sample_pattern(directory):
    """
    Dynamically detect sample pattern from FASTQ files.
    
    Parameters
    ----------
    directory : str
        Path to directory containing FASTQ files
        
    Returns
    -------
    str
        Detected sample pattern (e.g., "_S1_L001_")
    """
    try:
        fastq_files = [f for f in os.listdir(directory) if f.endswith('.fastq.gz')]
        if not fastq_files:
            return "_S1_L001_"  # Default fallback
        
        # Extract sample pattern from first file
        first_file = fastq_files[0]
        parts = first_file.split('_')
        if len(parts) >= 4:
            # Reconstruct pattern like _S1_L001_
            return '_' + '_'.join(parts[1:3]) + '_'
        return "_S1_L001_"  # Default fallback
    except Exception as e:
        logger.warning(f"Could not detect sample pattern in {directory}: {e}")
        return "_S1_L001_"  # Default fallback


def check_if_multiplexed(directory):
    """
    Check if sample needs demultiplexing by looking for index files.
    
    Parameters
    ----------
    directory : str
        Path to directory containing FASTQ files
        
    Returns
    -------
    bool
        True if I1 and I2 index files are present
    """
    try:
        files = os.listdir(directory)
        i1_files = [f for f in files if '_I1_' in f]
        i2_files = [f for f in files if '_I2_' in f]
        return len(i1_files) > 0 and len(i2_files) > 0
    except Exception as e:
        logger.warning(f"Could not check multiplexing status for {directory}: {e}")
        return False


def suggest_chemistry_based_on_files(directory):
    """
    Suggest chemistry options based on file patterns.
    
    Parameters
    ----------
    directory : str
        Path to directory containing FASTQ files
        
    Returns
    -------
    list
        Ordered list of suggested chemistries to try
    """
    if check_if_multiplexed(directory):
        logger.info("Sample appears to be multiplexed (has I1/I2 files).")
        # For multiplexed data, try these chemistries first
        return ["SC3Pv3", "SC3Pv3HT", "threeprime", "auto"]
    else:
        # For non-multiplexed data
        return ["auto", "threeprime", "fiveprime", "SC3Pv3"]


def format_error_output(result):
    """
    Format detailed error output from subprocess result.
    
    Parameters
    ----------
    result : subprocess.CompletedProcess
        Result from subprocess.run
        
    Returns
    -------
    str
        Formatted error message
    """
    error_lines = []
    if result.stdout:
        error_lines.append("=== STDOUT ===")
        error_lines.extend(result.stdout.split('\n')[-20:])  # Last 20 lines of stdout
    if result.stderr:
        error_lines.append("\n=== STDERR ===")
        error_lines.extend(result.stderr.split('\n')[-20:])  # Last 20 lines of stderr
    return "\n".join(error_lines)


def cleanup_directory(directory_path):
    """
    Safely remove a directory if it exists.
    
    Parameters
    ----------
    directory_path : str
        Path to directory to remove
        
    Returns
    -------
    bool
        True if cleanup succeeded or directory doesn't exist
    """
    if os.path.exists(directory_path):
        try:
            logger.info(f"Cleaning up previous run: {directory_path}")
            shutil.rmtree(directory_path)
            return True
        except Exception as e:
            logger.warning(f"Could not remove directory {directory_path}: {e}")
            return False
    return True  # Directory doesn't exist, so considered clean


def verify_cellranger_installed():
    """
    Verify that Cell Ranger is installed and accessible.
    
    Returns
    -------
    bool
        True if cellranger is found in PATH
    """
    try:
        result = subprocess.run(
            ["cellranger", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            logger.info(f"Cell Ranger version: {version}")
            return True
        return False
    except FileNotFoundError:
        return False


def get_fastq_dir_for_sample(sample_info):
    """
    Get the FASTQ directory for a sample from sample info dict.
    
    Parameters
    ----------
    sample_info : dict
        Sample information dictionary
        
    Returns
    -------
    str
        Path to FASTQ directory
    """
    # Handle both list and string formats for fastq paths
    fastq_r1 = sample_info.get('fastq_r1', [])
    if isinstance(fastq_r1, list) and len(fastq_r1) > 0:
        return os.path.dirname(fastq_r1[0])
    elif isinstance(fastq_r1, str):
        return os.path.dirname(fastq_r1)
    return None


# =============================================================================
# Main Cell Ranger Function
# =============================================================================

def run_cellranger_count(
    sample_id,
    fastq_dir,
    transcriptome_ref,
    output_dir,
    chemistry='auto',
    localcores=None,
    localmem=None,
    create_bam=True,
    expect_cells=None,
    force_cells=None,
    include_introns=True,
    no_bam=False
):
    """
    Run Cell Ranger count for a single sample with automatic chemistry detection.
    
    Parameters
    ----------
    sample_id : str
        Sample identifier (must match FASTQ file prefix)
    fastq_dir : str
        Path to directory containing FASTQ files
    transcriptome_ref : str
        Path to Cell Ranger transcriptome reference
    output_dir : str
        Directory where output will be written
    chemistry : str
        Chemistry type or 'auto' for automatic detection
    localcores : int, optional
        Number of cores to use (default: all available)
    localmem : int, optional
        Memory in GB (default: system limit)
    create_bam : bool
        Whether to create BAM file
    expect_cells : int, optional
        Expected number of cells
    force_cells : int, optional
        Force this number of cells
    include_introns : bool
        Include intronic reads
    no_bam : bool
        Skip BAM generation
        
    Returns
    -------
    dict
        Result dictionary with status, chemistry used, and output paths
    """
    
    # Verify Cell Ranger is installed
    if not verify_cellranger_installed():
        logger.error("Cell Ranger is not installed or not in PATH")
        logger.error("Please install Cell Ranger from 10X Genomics and add to PATH")
        return {'success': False, 'error': 'Cell Ranger not found'}
    
    # Set defaults
    if localcores is None:
        localcores = cpu_count()
    if localmem is None:
        localmem = 64  # Default to 64GB
        
    # Verify inputs
    if not os.path.exists(fastq_dir):
        logger.error(f"FASTQ directory not found: {fastq_dir}")
        return {'success': False, 'error': f'FASTQ directory not found: {fastq_dir}'}
        
    if not os.path.exists(transcriptome_ref):
        logger.error(f"Transcriptome reference not found: {transcriptome_ref}")
        return {'success': False, 'error': f'Reference not found: {transcriptome_ref}'}
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    original_dir = os.getcwd()
    os.chdir(output_dir)
    
    # Detect sample pattern
    sample_pattern = detect_sample_pattern(fastq_dir)
    logger.info(f"Detected sample pattern for {sample_id}: {sample_pattern}")
    
    # Check if multiplexed
    if check_if_multiplexed(fastq_dir):
        logger.info(f"Sample {sample_id} has index files (I1/I2). May require special handling.")
    
    # Determine chemistries to try
    if chemistry == 'auto':
        suggested_chemistries = suggest_chemistry_based_on_files(fastq_dir)
        trial_chemistries = suggested_chemistries + [c for c in CHEMISTRY_OPTIONS if c not in suggested_chemistries]
    else:
        # User specified chemistry - try that first, then fallback to auto-detection
        trial_chemistries = [chemistry] + CHEMISTRY_OPTIONS
        trial_chemistries = list(dict.fromkeys(trial_chemistries))  # Remove duplicates, preserve order
    
    success = False
    last_errors = {}
    successful_chemistry = None
    result_paths = {}
    
    for chem in trial_chemistries:
        run_id = f"{sample_id}{sample_pattern}"
        run_dir = os.path.join(output_dir, run_id)
        
        # Clean up any previous attempts
        if not cleanup_directory(run_dir):
            logger.warning(f"Could not clean up {run_dir}. Skipping chemistry {chem}.")
            last_errors[chem] = "Could not clean up previous run directory"
            continue
        
        try:
            logger.info(f"Trying chemistry {chem} for sample {sample_id}")
            
            # Build command
            cmd = [
                "cellranger", "count",
                f"--id={run_id}",
                f"--fastqs={fastq_dir}",
                f"--sample={sample_id}",
                f"--transcriptome={transcriptome_ref}",
                f"--localcores={localcores}",
                f"--localmem={localmem}",
                f"--chemistry={chem}"
            ]
            
            # Optional arguments
            if create_bam and not no_bam:
                cmd.append("--create-bam=true")
            elif no_bam:
                cmd.append("--create-bam=false")
                
            if expect_cells:
                cmd.append(f"--expect-cells={expect_cells}")
            if force_cells:
                cmd.append(f"--force-cells={force_cells}")
            if include_introns:
                cmd.append("--include-introns=true")
            
            logger.info(f"Running: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            if result.returncode == 0:
                logger.info(f"SUCCESS with chemistry {chem} for sample {sample_id}")
                success = True
                successful_chemistry = chem
                
                # Record output paths
                outs_dir = os.path.join(run_dir, 'outs')
                result_paths = {
                    'run_dir': run_dir,
                    'outs_dir': outs_dir,
                    'matrix_h5': os.path.join(outs_dir, 'filtered_feature_bc_matrix.h5'),
                    'matrix_dir': os.path.join(outs_dir, 'filtered_feature_bc_matrix'),
                    'bam': os.path.join(outs_dir, 'possorted_genome_bam.bam'),
                    'bai': os.path.join(outs_dir, 'possorted_genome_bam.bam.bai'),
                    'summary': os.path.join(outs_dir, 'web_summary.html'),
                    'metrics': os.path.join(outs_dir, 'metrics_summary.csv'),
                }
                break
                
            else:
                error_output = format_error_output(result)
                last_errors[chem] = error_output
                logger.warning(f"Chemistry {chem} failed for sample {sample_id}")
                logger.debug(f"Error details:\n{error_output}")
                cleanup_directory(run_dir)
                
        except Exception as e:
            error_msg = f"Exception with chemistry {chem}:\n{traceback.format_exc()}"
            last_errors[chem] = error_msg
            logger.error(error_msg)
            cleanup_directory(run_dir)
            continue
    
    os.chdir(original_dir)
    
    if success:
        return {
            'success': True,
            'sample_id': sample_id,
            'chemistry': successful_chemistry,
            'paths': result_paths
        }
    else:
        logger.error(f"All chemistry options failed for sample {sample_id}")
        logger.error("Last errors encountered:")
        for chem, error in list(last_errors.items())[-3:]:
            first_line = error.splitlines()[0] if error else 'Unknown error'
            logger.error(f"  Chemistry {chem}: {first_line}")
        return {
            'success': False,
            'sample_id': sample_id,
            'errors': last_errors
        }


# =============================================================================
# Snakemake Integration
# =============================================================================

def run_from_snakemake():
    """Run Cell Ranger from Snakemake rule."""
    
    # Get Snakemake variables
    sample_id = snakemake.wildcards.sample
    output_matrix = snakemake.output.matrix
    output_bam = snakemake.output.bam
    output_summary = snakemake.output.summary
    log_file = snakemake.log[0]
    
    params = snakemake.params
    transcriptome_ref = params.transcriptome
    chemistry = params.chemistry
    localcores = params.localcores
    localmem = params.localmem
    expect_cells = params.get('expect_cells')
    include_introns = params.get('include_introns', True)
    
    # Get FASTQ directory from sample info
    samples = snakemake.config.get('samples', {})
    sample_info = samples.get(sample_id, {})
    fastq_dir = get_fastq_dir_for_sample(sample_info)
    
    if not fastq_dir:
        logger.error(f"Could not determine FASTQ directory for sample {sample_id}")
        sys.exit(1)
    
    # Determine output directory (parent of where outputs should go)
    output_dir = os.path.dirname(os.path.dirname(output_matrix))
    
    # Set up file logging
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    logger.info(f"Starting Cell Ranger count for sample: {sample_id}")
    logger.info(f"FASTQ directory: {fastq_dir}")
    logger.info(f"Transcriptome reference: {transcriptome_ref}")
    logger.info(f"Output directory: {output_dir}")
    
    # Run Cell Ranger
    result = run_cellranger_count(
        sample_id=sample_id,
        fastq_dir=fastq_dir,
        transcriptome_ref=transcriptome_ref,
        output_dir=output_dir,
        chemistry=chemistry,
        localcores=localcores,
        localmem=localmem,
        expect_cells=expect_cells,
        include_introns=include_introns
    )
    
    if not result['success']:
        logger.error(f"Cell Ranger failed for sample {sample_id}")
        sys.exit(1)
    
    logger.info(f"Cell Ranger completed successfully for sample {sample_id}")
    logger.info(f"Chemistry used: {result['chemistry']}")


# =============================================================================
# Standalone CLI
# =============================================================================

def main():
    """Main function for standalone CLI usage."""
    
    parser = argparse.ArgumentParser(
        description='Run Cell Ranger count with automatic chemistry detection'
    )
    parser.add_argument('--config', required=True, help='Path to config YAML file')
    parser.add_argument('--sample', help='Specific sample ID to process (optional)')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--transcriptome', required=True, help='Path to transcriptome reference')
    parser.add_argument('--chemistry', default='auto', help='Chemistry type (default: auto)')
    parser.add_argument('--cores', type=int, default=None, help='Number of cores')
    parser.add_argument('--memory', type=int, default=64, help='Memory in GB')
    
    args = parser.parse_args()
    
    # Load config
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    samples = config.get('samples', {})
    
    if args.sample:
        # Process single sample
        sample_ids = [args.sample]
    else:
        # Process all samples
        sample_ids = list(samples.keys())
    
    results = {}
    for sample_id in sample_ids:
        sample_info = samples.get(sample_id, {})
        fastq_dir = get_fastq_dir_for_sample(sample_info)
        
        if not fastq_dir:
            logger.error(f"Could not determine FASTQ directory for sample {sample_id}")
            continue
        
        result = run_cellranger_count(
            sample_id=sample_id,
            fastq_dir=fastq_dir,
            transcriptome_ref=args.transcriptome,
            output_dir=args.output_dir,
            chemistry=args.chemistry,
            localcores=args.cores,
            localmem=args.memory
        )
        results[sample_id] = result
    
    # Summary
    successful = sum(1 for r in results.values() if r.get('success'))
    logger.info(f"\n=== Processing Summary ===")
    logger.info(f"Successful: {successful}/{len(results)}")
    
    for sample_id, result in results.items():
        status = "SUCCESS" if result.get('success') else "FAILED"
        chem = result.get('chemistry', 'N/A')
        logger.info(f"  {sample_id}: {status} (chemistry: {chem})")
    
    # Return exit code based on results
    if successful == len(results):
        sys.exit(0)
    elif successful > 0:
        sys.exit(1)  # Partial success
    else:
        sys.exit(2)  # Complete failure


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == '__main__':
    # Check if running from Snakemake
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        # Running standalone
        main()
