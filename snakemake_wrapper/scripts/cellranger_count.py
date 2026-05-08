#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cellranger_count.py
===================

Run Cell Ranger count on single-cell FASTQ files with automatic chemistry detection.

This script is designed to work with SRAscraper output structure where:
- Pickle file contains: {GSE_ID: DataFrame with run_accession column}
- FASTQs are organized as: {fastq_base_dir}/{GSE_ID}/{SRR_ID}/*.fastq.gz
- Each SRR is processed separately (one Cell Ranger run per SRR)

The script processes all SRRs within a GSE and tracks successful runs,
matching the behavior of the standalone SRAscraper cellranger script.

Snakemake Integration:
    The declared Snakemake output is a status JSON file (not H5/BAM/etc).
    Cell Ranger outputs (H5, BAM, web_summary) are side effects on disk.
    A downstream checkpoint rule validates which samples produced usable
    outputs and builds the dynamic DAG for QC, annotation, and beyond.
    This ensures one failed sample never blocks the entire pipeline.

Usage:
    Called via Snakemake rule with snakemake.input/output/params
    
    Or standalone:
    python cellranger_count.py --gse-id GSE123456 --fastq-base-dir /path/to/fastq ...

Requirements:
    - Cell Ranger must be installed and available in PATH
    - Reference transcriptome must be downloaded
"""

import os
import sys
import json
import yaml
import pickle
import shutil
import subprocess
import traceback
import logging
import argparse
import pandas as pd
from pathlib import Path
from os import cpu_count
from datetime import datetime

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
# Matches the order from original script
CHEMISTRY_OPTIONS = [
    "auto",           # Try auto first
    "SC3Pv4",
    "SC3Pv3",
    "SC3Pv2",
    "SC3Pv3HT",
    "SC3Pv3LT",
    "threeprime",
    "fiveprime",
    "SC5P-PE-v3",
    "SC5P-PE",
    "SC5P-R2-v3",
    "SC5P-R2",
    "SC3Pv1",         # These last two can't be auto-detected
    "ARC-v1"
]

# Default sample pattern for Cell Ranger --id
SAMPLE_PATTERN = "_S1_L001_"


# =============================================================================
# Status JSON Helper
# =============================================================================

def write_status_json(output_path, sample_id, success, chemistry=None,
                      error=None, paths=None, series_id=None):
    """
    Write a status JSON file for the checkpoint to consume.
    
    This is the ONLY declared Snakemake output for the cellranger_count rule.
    It is always written regardless of success or failure, ensuring the
    Snakemake DAG is never broken by a single sample failure.
    
    Parameters
    ----------
    output_path : str
        Path to write the status JSON
    sample_id : str
        SRR or GSE sample identifier
    success : bool
        Whether Cell Ranger completed successfully
    chemistry : str, optional
        Chemistry that succeeded (None if failed)
    error : str, optional
        Error message (None if succeeded)
    paths : dict, optional
        Dict of output file paths (None if failed)
    series_id : str, optional
        GSE series ID for this SRR
    """
    status = {
        'sample_id': sample_id,
        'success': success,
        'chemistry': chemistry,
        'error': error,
        'series_id': series_id,
        'timestamp': datetime.now().isoformat(),
        'paths': paths or {},
    }
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(status, f, indent=2)
    
    if success:
        logger.info(f"Status JSON written: {output_path} (SUCCESS)")
    else:
        logger.warning(f"Status JSON written: {output_path} (FAILED: {error})")


# =============================================================================
# Helper Functions
# =============================================================================

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
        error_lines.extend(result.stdout.split('\n')[-50:])  # Last 50 lines of stdout
    if result.stderr:
        error_lines.append("\n=== STDERR ===")
        error_lines.extend(result.stderr.split('\n')[-50:])  # Last 50 lines of stderr
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


def get_srr_ids_from_sample_info(sample_info, gse_id):
    """
    Extract SRR IDs from sample info (DataFrame or dict).
    
    Parameters
    ----------
    sample_info : DataFrame or dict
        Sample information from pickle file
    gse_id : str
        GSE series ID
        
    Returns
    -------
    list
        List of SRR run accession IDs
    """
    if isinstance(sample_info, pd.DataFrame):
        if len(sample_info) == 0:
            logger.error(f"Empty DataFrame for {gse_id} - no runs available")
            return []
        
        if 'run_accession' in sample_info.columns:
            return sample_info['run_accession'].tolist()
        elif 'run_alias' in sample_info.columns:
            # run_alias might have format like "SRR14340883_GSM5234567"
            return [r.split('_')[0] for r in sample_info['run_alias'].tolist()]
        else:
            logger.error(f"Could not find run accession column in DataFrame for {gse_id}")
            logger.error(f"Available columns: {list(sample_info.columns)}")
            return []
    
    elif isinstance(sample_info, dict):
        # Handle dict format if present
        if 'run_accession' in sample_info:
            runs = sample_info['run_accession']
            return runs if isinstance(runs, list) else [runs]
        if 'srr_ids' in sample_info:
            return sample_info['srr_ids']
    
    return []


# =============================================================================
# Cell Ranger Execution Functions
# =============================================================================

def run_cellranger_for_srr(
    srr_id,
    fastq_dir,
    transcriptome_ref,
    output_base_dir,
    chemistry='auto',
    localcores=None,
    localmem=None,
    create_bam=True,
    expect_cells=None,
    include_introns=True,
):
    """
    Run Cell Ranger count for a single SRR with automatic chemistry detection.
    
    This matches custom individual script behavior I used before making this pipeline:
    - --id={srr_id}_S1_L001_
    - --fastqs={fastq_dir}
    - --sample={srr_id}
    
    Parameters
    ----------
    srr_id : str
        SRR run accession ID
    fastq_dir : str
        Path to directory containing FASTQ files for this SRR
    transcriptome_ref : str
        Path to Cell Ranger transcriptome reference
    output_base_dir : str
        Base directory where Cell Ranger output will be written
    chemistry : str
        Chemistry type or 'auto' for automatic detection
    localcores : int, optional
        Number of cores to use (default: all available)
    localmem : int, optional
        Memory in GB (default: 64)
    create_bam : bool
        Whether to create BAM file
    expect_cells : int, optional
        Expected number of cells
    include_introns : bool
        Include intronic reads
        
    Returns
    -------
    dict
        Result dictionary with status, chemistry used, and output paths
    """
    # Set defaults
    if localcores is None:
        localcores = cpu_count()
    if localmem is None:
        localmem = 64
    
    # Verify FASTQ directory exists
    if not os.path.exists(fastq_dir):
        logger.error(f"FASTQ directory not found: {fastq_dir}")
        return {'success': False, 'srr_id': srr_id, 'error': f'FASTQ directory not found: {fastq_dir}'}
    
    # Check for FASTQ files
    fastq_files = [f for f in os.listdir(fastq_dir) if f.endswith('.fastq.gz')]
    if not fastq_files:
        logger.error(f"No FASTQ files found in: {fastq_dir}")
        return {'success': False, 'srr_id': srr_id, 'error': f'No FASTQ files in {fastq_dir}'}
    
    logger.info(f"Found {len(fastq_files)} FASTQ files in {fastq_dir}")
    
    # Cell Ranger output ID
    cellranger_id = f"{srr_id}{SAMPLE_PATTERN}"
    
    # Cell Ranger creates output in: {cwd}/{cellranger_id}/outs/
    # We run from the output base directory
    os.makedirs(output_base_dir, exist_ok=True)
    original_dir = os.getcwd()
    os.chdir(output_base_dir)
    
    success = False
    last_errors = {}
    successful_chemistry = None
    result_paths = {}
    
    for chem in CHEMISTRY_OPTIONS:
        # Cell Ranger output directory
        run_dir = os.path.join(output_base_dir, cellranger_id)
        
        # Clean up any previous attempts
        if not cleanup_directory(run_dir):
            logger.warning(f"Could not clean up {run_dir}. Skipping chemistry {chem}.")
            last_errors[chem] = "Could not clean up previous run directory"
            continue
        
        try:
            logger.info(f"Trying chemistry {chem} for SRR {srr_id}")
            
            # Build command
            cmd = [
                "cellranger", "count",
                f"--id={cellranger_id}",
                f"--fastqs={fastq_dir}",
                f"--sample={srr_id}",
                f"--transcriptome={transcriptome_ref}",
                f"--localcores={localcores}",
                f"--chemistry={chem}"
            ]
            
            # Optional arguments
            if create_bam:
                cmd.append("--create-bam=true")
            else:
                cmd.append("--create-bam=false")
                
            if expect_cells:
                cmd.append(f"--expect-cells={expect_cells}")
            if include_introns:
                cmd.append("--include-introns=true")
            if localmem:
                cmd.append(f"--localmem={localmem}")
            
            logger.info(f"Running: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            if result.returncode == 0:
                logger.info(f"SUCCESS with chemistry {chem} for SRR {srr_id}")
                success = True
                successful_chemistry = chem
                
                # Record output paths
                outs_dir = os.path.join(run_dir, 'outs')
                result_paths = {
                    'run_dir': run_dir,
                    'outs_dir': outs_dir,
                    'cellranger_id': cellranger_id,
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
                logger.warning(f"Chemistry {chem} failed for SRR {srr_id}")
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
            'srr_id': srr_id,
            'chemistry': successful_chemistry,
            'paths': result_paths
        }
    else:
        logger.error(f"All chemistry options failed for SRR {srr_id}")
        return {
            'success': False,
            'srr_id': srr_id,
            'errors': last_errors
        }


def run_cellranger_for_gse(
    gse_id,
    sample_info,
    fastq_base_dir,
    transcriptome_ref,
    output_dir,
    chemistry='auto',
    localcores=None,
    localmem=None,
    create_bam=True,
    expect_cells=None,
    include_introns=True,
):
    """
    Run Cell Ranger for all SRRs within a GSE series.
    
    Matches previous workflow:
    - Iterates over each run_accession in the GSE DataFrame
    - Runs Cell Ranger separately for each SRR
    - Tracks successful runs
    
    Parameters
    ----------
    gse_id : str
        GSE series ID
    sample_info : DataFrame
        Sample information DataFrame from pickle file
    fastq_base_dir : str
        Base directory containing FASTQs: {fastq_base_dir}/{GSE_ID}/{SRR_ID}/
    transcriptome_ref : str
        Path to Cell Ranger transcriptome reference
    output_dir : str
        Output directory for Cell Ranger results
    chemistry : str
        Chemistry type or 'auto'
    localcores : int, optional
        Number of cores
    localmem : int, optional
        Memory in GB
    create_bam : bool
        Whether to create BAM
    expect_cells : int, optional
        Expected cells
    include_introns : bool
        Include intronic reads
        
    Returns
    -------
    dict
        Results with successful_runs DataFrame and failed_runs list
    """
    # Verify Cell Ranger is installed
    if not verify_cellranger_installed():
        logger.error("Cell Ranger is not installed or not in PATH")
        return {'success': False, 'error': 'Cell Ranger not found'}
    
    # Get SRR IDs from sample info
    srr_ids = get_srr_ids_from_sample_info(sample_info, gse_id)
    
    if not srr_ids:
        logger.error(f"No SRR IDs found for {gse_id}")
        return {'success': False, 'error': 'No SRR IDs found'}
    
    logger.info(f"Processing {len(srr_ids)} SRR runs for {gse_id}")
    
    # Output directory for this GSE's Cell Ranger runs
    gse_output_dir = os.path.join(output_dir, 'cellranger', gse_id)
    os.makedirs(gse_output_dir, exist_ok=True)
    
    # Track results
    successful_runs = []
    failed_runs = []
    all_results = {}
    
    for srr_id in srr_ids:
        logger.info("="*60)
        logger.info(f"Processing SRR: {srr_id}")
        logger.info("="*60)
        
        # FASTQ directory for this SRR: {fastq_base_dir}/{GSE_ID}/{SRR_ID}/
        srr_fastq_dir = os.path.join(fastq_base_dir, gse_id, srr_id)
        
        if not os.path.exists(srr_fastq_dir):
            logger.warning(f"FASTQ directory not found for {srr_id}: {srr_fastq_dir}")
            failed_runs.append({'srr_id': srr_id, 'error': 'FASTQ directory not found'})
            continue
        
        # Run Cell Ranger for this SRR
        result = run_cellranger_for_srr(
            srr_id=srr_id,
            fastq_dir=srr_fastq_dir,
            transcriptome_ref=transcriptome_ref,
            output_base_dir=gse_output_dir,
            chemistry=chemistry,
            localcores=localcores,
            localmem=localmem,
            create_bam=create_bam,
            expect_cells=expect_cells,
            include_introns=include_introns,
        )
        
        all_results[srr_id] = result
        
        if result['success']:
            successful_runs.append(srr_id)
            logger.info(f"SUCCESS: {srr_id} (chemistry: {result['chemistry']})")
        else:
            failed_runs.append({'srr_id': srr_id, 'error': result.get('error', 'Unknown error')})
            logger.warning(f"FAILED: {srr_id}")
    
    # Create successful runs DataFrame
    successful_df = pd.DataFrame()
    if isinstance(sample_info, pd.DataFrame) and successful_runs:
        successful_df = sample_info[sample_info['run_accession'].isin(successful_runs)].copy()
    
    # Summary
    logger.info("="*60)
    logger.info(f"SUMMARY for {gse_id}")
    logger.info("="*60)
    logger.info(f"Total SRRs: {len(srr_ids)}")
    logger.info(f"Successful: {len(successful_runs)}")
    logger.info(f"Failed: {len(failed_runs)}")
    
    if failed_runs:
        logger.warning("Failed runs:")
        for fail in failed_runs:
            logger.warning(f"  {fail['srr_id']}: {fail['error']}")
    
    return {
        'success': len(successful_runs) > 0,
        'gse_id': gse_id,
        'total_runs': len(srr_ids),
        'successful_runs': successful_runs,
        'failed_runs': failed_runs,
        'successful_df': successful_df,
        'all_results': all_results,
        'output_dir': gse_output_dir,
    }


# =============================================================================
# Snakemake Integration
# =============================================================================

def run_from_snakemake():
    """
    Run Cell Ranger from Snakemake rule.
    
    The ONLY declared Snakemake output is a status JSON file.
    Cell Ranger outputs (H5, BAM, etc.) are side effects on disk — the
    downstream checkpoint rule validates which ones are usable.
    
    This function NEVER calls sys.exit(1) for sample-level failures.
    Instead it writes a failure status JSON so the pipeline continues.
    """
    
    # Get Snakemake variables
    sample_id = snakemake.params.sample_id
    status_output = snakemake.output.status
    log_file = snakemake.log[0] if snakemake.log else None
    
    # Get parameters from Snakemake
    params = snakemake.params
    transcriptome_ref = params.transcriptome
    chemistry = params.chemistry
    localcores = params.localcores
    localmem = params.localmem
    expect_cells = params.expect_cells
    include_introns = params.include_introns
    create_bam = params.create_bam
    output_dir = params.output_dir
    fastq_base_dir = getattr(params, 'fastq_base_dir', None)
    
    # Get sample information from params
    samples_dict = params.samples_dict
    sample_info = samples_dict.get(sample_id)
    
    # Set up file logging
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    logger.info("="*60)
    logger.info(f"CELL RANGER COUNT - GSE: {sample_id}")
    logger.info("="*60)
    
    # =========================================================================
    # Validate inputs — write failure status JSON instead of sys.exit(1)
    # =========================================================================
    if sample_info is None:
        logger.error(f"No sample info found for {sample_id}")
        write_status_json(status_output, sample_id, success=False,
                          error=f"No sample info found for {sample_id}")
        return
    
    if isinstance(sample_info, pd.DataFrame):
        logger.info(f"Sample info: DataFrame with {len(sample_info)} runs")
        if len(sample_info) == 0:
            logger.error(f"Empty DataFrame for {sample_id} - no runs to process")
            write_status_json(status_output, sample_id, success=False,
                              error=f"Empty DataFrame for {sample_id}")
            return
    else:
        logger.info(f"Sample info type: {type(sample_info)}")
    
    if not fastq_base_dir:
        logger.error("fastq_base_dir not provided in config")
        write_status_json(status_output, sample_id, success=False,
                          error="fastq_base_dir not provided in config")
        return
    
    if not os.path.exists(fastq_base_dir):
        logger.error(f"FASTQ base directory not found: {fastq_base_dir}")
        write_status_json(status_output, sample_id, success=False,
                          error=f"FASTQ base directory not found: {fastq_base_dir}")
        return
    
    logger.info(f"FASTQ base directory: {fastq_base_dir}")
    logger.info(f"Transcriptome reference: {transcriptome_ref}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Chemistry: {chemistry}")
    logger.info(f"Cores: {localcores}, Memory: {localmem}GB")
    
    # =========================================================================
    # Determine processing mode: SRR-level or GSE-level
    # =========================================================================
    is_srr_level = isinstance(sample_info, dict) and 'run_accession' in sample_info
    
    if is_srr_level:
        # =================================================================
        # SRR-level processing (new behavior with Snakefile fix)
        # =================================================================
        srr_id = sample_info.get('run_accession', sample_id)
        series_id = sample_info.get('series_id', sample_id)
        
        logger.info(f"SRR-level mode: {srr_id} (series: {series_id})")
        
        # FASTQ directory: {fastq_base_dir}/{series_id}/{srr_id}/
        srr_fastq_dir = os.path.join(fastq_base_dir, series_id, srr_id)
        
        if not os.path.exists(srr_fastq_dir):
            logger.error(f"FASTQ directory not found: {srr_fastq_dir}")
            write_status_json(status_output, sample_id, success=False,
                              error=f"FASTQ directory not found: {srr_fastq_dir}",
                              series_id=series_id)
            return
        
        # Output goes under: cellranger/{series_id}/{srr}_S1_L001_/outs/
        gse_output_dir = os.path.join(output_dir, 'cellranger', series_id)
        os.makedirs(gse_output_dir, exist_ok=True)
        
        result_srr = run_cellranger_for_srr(
            srr_id=srr_id,
            fastq_dir=srr_fastq_dir,
            transcriptome_ref=transcriptome_ref,
            output_base_dir=gse_output_dir,
            chemistry=chemistry,
            localcores=localcores,
            localmem=localmem,
            create_bam=create_bam,
            expect_cells=expect_cells,
            include_introns=include_introns,
        )
        
        if not result_srr['success']:
            logger.error(f"Cell Ranger failed for {srr_id}")
            write_status_json(status_output, sample_id, success=False,
                              error="All chemistry options failed",
                              series_id=series_id)
            logger.info("="*60)
            logger.info(f"Cell Ranger FAILED for {sample_id}")
            logger.info(f"Sample will be excluded by checkpoint validation")
            logger.info("="*60)
            return
        
        # Create per-SRR symlink:
        #   cellranger/{srr_id}/outs → cellranger/{series_id}/{srr}_S1_L001_/outs
        actual_outs = os.path.join(gse_output_dir, f"{srr_id}{SAMPLE_PATTERN}", 'outs')
        expected_outs_dir = os.path.join(output_dir, 'cellranger', srr_id, 'outs')
        
        if os.path.exists(actual_outs):
            os.makedirs(os.path.dirname(expected_outs_dir), exist_ok=True)
            
            if os.path.exists(expected_outs_dir):
                if os.path.islink(expected_outs_dir):
                    os.unlink(expected_outs_dir)
                else:
                    shutil.rmtree(expected_outs_dir)
            
            os.symlink(actual_outs, expected_outs_dir)
            logger.info(f"Created per-SRR symlink: {expected_outs_dir} -> {actual_outs}")
        
        # Write success status JSON
        write_status_json(
            status_output, sample_id, success=True,
            chemistry=result_srr['chemistry'],
            series_id=series_id,
            paths=result_srr.get('paths', {}),
        )
        
        logger.info("="*60)
        logger.info(f"Cell Ranger completed for {sample_id}")
        logger.info(f"Successful runs: 1/1")
        logger.info("="*60)
        
    else:
        # =================================================================
        # GSE-level processing (original behavior for backward compatibility)
        # =================================================================
        result = run_cellranger_for_gse(
            gse_id=sample_id,
            sample_info=sample_info,
            fastq_base_dir=fastq_base_dir,
            transcriptome_ref=transcriptome_ref,
            output_dir=output_dir,
            chemistry=chemistry,
            localcores=localcores,
            localmem=localmem,
            create_bam=create_bam,
            expect_cells=expect_cells,
            include_introns=include_introns,
        )
        
        if not result['success']:
            logger.error(f"Cell Ranger failed for all runs in {sample_id}")
            write_status_json(status_output, sample_id, success=False,
                              error="All runs in GSE failed")
            return
        
        # Create per-SRR symlinks for ALL successful runs (not just the first)
        # This ensures downstream modules can access each sample individually
        for srr_id in result['successful_runs']:
            actual_outs = os.path.join(
                result['output_dir'], f"{srr_id}{SAMPLE_PATTERN}", 'outs'
            )
            expected_srr_dir = os.path.join(output_dir, 'cellranger', srr_id, 'outs')
            
            if os.path.exists(actual_outs):
                os.makedirs(os.path.dirname(expected_srr_dir), exist_ok=True)
                
                if os.path.exists(expected_srr_dir):
                    if os.path.islink(expected_srr_dir):
                        os.unlink(expected_srr_dir)
                    else:
                        shutil.rmtree(expected_srr_dir)
                
                os.symlink(actual_outs, expected_srr_dir)
                logger.info(f"Created per-SRR symlink: {expected_srr_dir} -> {actual_outs}")
        
        # Also create the GSE-level symlink for backward compatibility
        # (points to first successful run, same as before)
        expected_outs_dir = os.path.join(output_dir, 'cellranger', sample_id, 'outs')
        if result['successful_runs']:
            first_srr = result['successful_runs'][0]
            first_run_dir = os.path.join(
                result['output_dir'], f"{first_srr}{SAMPLE_PATTERN}", 'outs'
            )
            
            if os.path.exists(first_run_dir):
                os.makedirs(os.path.dirname(expected_outs_dir), exist_ok=True)
                
                if os.path.exists(expected_outs_dir):
                    if os.path.islink(expected_outs_dir):
                        os.unlink(expected_outs_dir)
                    else:
                        shutil.rmtree(expected_outs_dir)
                
                os.symlink(first_run_dir, expected_outs_dir)
                logger.info(f"Created GSE-level symlink: {expected_outs_dir} -> {first_run_dir}")
        
        # Write success status JSON
        write_status_json(
            status_output, sample_id, success=True,
            chemistry='mixed',
            paths={'successful_runs': result['successful_runs'],
                   'failed_runs': [f['srr_id'] for f in result['failed_runs']],
                   'output_dir': result['output_dir']},
        )
        
        logger.info("="*60)
        logger.info(f"Cell Ranger completed for {sample_id}")
        logger.info(f"Successful runs: {len(result['successful_runs'])}/{result['total_runs']}")
        logger.info("="*60)


# =============================================================================
# Standalone CLI
# =============================================================================

def main():
    """Main function for standalone CLI usage."""
    
    parser = argparse.ArgumentParser(
        description='Run Cell Ranger count for SRAscraper output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run for a single GSE with pickle file
    python cellranger_count.py \\
        --gse-id GSE173468 \\
        --pickle-file /path/to/dictionary_file.pkl \\
        --fastq-base-dir /path/to/fastq \\
        --transcriptome /path/to/GRCh38 \\
        --output-dir /path/to/output

    # Run for a single SRR directly
    python cellranger_count.py \\
        --srr-id SRR14340883 \\
        --fastq-dir /path/to/fastq/GSE173468/SRR14340883 \\
        --transcriptome /path/to/GRCh38 \\
        --output-dir /path/to/output
        """
    )
    
    # GSE-level options
    gse_group = parser.add_argument_group('GSE-level processing')
    gse_group.add_argument('--gse-id', help='GSE series ID')
    gse_group.add_argument('--pickle-file', help='Path to SRAscraper pickle file')
    gse_group.add_argument('--fastq-base-dir', help='Base directory containing FASTQs ({dir}/{GSE}/{SRR}/)')
    
    # SRR-level options (direct processing)
    srr_group = parser.add_argument_group('SRR-level processing (single run)')
    srr_group.add_argument('--srr-id', help='SRR run accession ID')
    srr_group.add_argument('--fastq-dir', help='Path to FASTQ directory for single SRR')
    
    # Common options
    common = parser.add_argument_group('Common options')
    common.add_argument('--transcriptome', required=True, help='Path to transcriptome reference')
    common.add_argument('--output-dir', required=True, help='Output directory')
    common.add_argument('--chemistry', default='auto', help='Chemistry type (default: auto)')
    common.add_argument('--cores', type=int, default=None, help='Number of cores')
    common.add_argument('--memory', type=int, default=64, help='Memory in GB')
    common.add_argument('--expect-cells', type=int, default=None, help='Expected number of cells')
    common.add_argument('--include-introns', action='store_true', default=True, help='Include intronic reads')
    common.add_argument('--no-bam', action='store_true', help='Skip BAM generation')
    
    args = parser.parse_args()
    
    # Determine mode: GSE-level or SRR-level
    if args.gse_id and args.pickle_file:
        # GSE-level processing
        logger.info(f"Running in GSE mode for {args.gse_id}")
        
        # Load pickle file
        with open(args.pickle_file, 'rb') as f:
            gse_dict = pickle.load(f)
        
        if args.gse_id not in gse_dict:
            logger.error(f"GSE {args.gse_id} not found in pickle file")
            sys.exit(1)
        
        sample_info = gse_dict[args.gse_id]
        
        result = run_cellranger_for_gse(
            gse_id=args.gse_id,
            sample_info=sample_info,
            fastq_base_dir=args.fastq_base_dir,
            transcriptome_ref=args.transcriptome,
            output_dir=args.output_dir,
            chemistry=args.chemistry,
            localcores=args.cores,
            localmem=args.memory,
            create_bam=not args.no_bam,
            expect_cells=args.expect_cells,
            include_introns=args.include_introns,
        )
        
        if result['success']:
            logger.info(f"SUCCESS: {len(result['successful_runs'])}/{result['total_runs']} runs completed")
            sys.exit(0)
        else:
            logger.error("FAILED: No successful runs")
            sys.exit(1)
            
    elif args.srr_id and args.fastq_dir:
        # SRR-level processing
        logger.info(f"Running in SRR mode for {args.srr_id}")
        
        result = run_cellranger_for_srr(
            srr_id=args.srr_id,
            fastq_dir=args.fastq_dir,
            transcriptome_ref=args.transcriptome,
            output_base_dir=args.output_dir,
            chemistry=args.chemistry,
            localcores=args.cores,
            localmem=args.memory,
            create_bam=not args.no_bam,
            expect_cells=args.expect_cells,
            include_introns=args.include_introns,
        )
        
        if result['success']:
            logger.info(f"SUCCESS - Chemistry: {result['chemistry']}")
            sys.exit(0)
        else:
            logger.error("FAILED")
            sys.exit(1)
    else:
        parser.error("Must provide either (--gse-id and --pickle-file and --fastq-base-dir) or (--srr-id and --fastq-dir)")


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
