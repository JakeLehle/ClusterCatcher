#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_cellranger.py
======================

Checkpoint validation script for Cell Ranger outputs.

This script runs AFTER all cellranger_count jobs complete and BEFORE any
downstream analysis (QC, annotation, viral detection, etc.). It serves as
the bridge between Cell Ranger's per-sample processing and the aggregate
analysis steps, mirroring the dictionary_file_filtered.pkl pattern from
the standalone Cell Ranger processing script.

What it does:
1. Reads every cellranger_status.json to determine which runs succeeded
2. For samples that reported success, validates outputs on disk:
   - barcodes.tsv.gz > 100 bytes  (catches empty barcode whitelists)
   - filtered_feature_bc_matrix.h5 > 1000 bytes  (catches empty H5 files)
   - possorted_genome_bam.bam > 10,000 bytes  (catches empty BAMs)
3. Writes validated_samples.txt — the dynamic sample list for all downstream rules
4. Backs up the original pickle and writes a filtered version
5. Writes validation_summary.json with full details and skipped_samples.tsv

The Snakemake checkpoint mechanism re-evaluates the DAG after this script
completes, so downstream rules only see validated samples.

Usage:
    Called via Snakemake checkpoint rule — not intended for standalone use.
"""

import os
import sys
import json
import pickle
import shutil
import logging
from pathlib import Path
from datetime import datetime

import pandas as pd

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Validation Thresholds
# =============================================================================
# These thresholds are based on observed failure modes:
# - barcodes.tsv.gz = 33 bytes when Cell Ranger finds no valid barcodes
#   (chemistry mismatch that still exits 0)
# - H5 file can exist but contain only headers with no actual data
# - BAM can exist but be essentially empty (no aligned reads)

THRESHOLDS = {
    'barcodes_tsv_gz': 100,      # bytes — empty gzip is ~33 bytes
    'filtered_matrix_h5': 1000,  # bytes — empty H5 header is ~few hundred
    'bam': 10000,                # bytes — empty BAM header is ~few thousand
}


# =============================================================================
# Validation Functions
# =============================================================================

def validate_sample_outputs(cellranger_dir, sample_id):
    """
    Validate Cell Ranger outputs for a single sample on disk.
    
    Checks file existence and size thresholds to catch cases where
    Cell Ranger reports success but produces empty/unusable outputs
    (e.g., chemistry mismatch with no barcode matches).
    
    Parameters
    ----------
    cellranger_dir : str
        Base cellranger output directory (contains {sample}/outs/)
    sample_id : str
        Sample identifier (SRR ID)
        
    Returns
    -------
    dict
        Validation result with 'valid' bool, file sizes, and any error
    """
    outs_dir = os.path.join(cellranger_dir, sample_id, 'outs')
    matrix_dir = os.path.join(outs_dir, 'filtered_feature_bc_matrix')
    
    result = {
        'valid': False,
        'checks': {},
        'error': None,
    }
    
    # ---- Check 1: barcodes.tsv.gz > 100 bytes ----
    barcodes_path = os.path.join(matrix_dir, 'barcodes.tsv.gz')
    if os.path.exists(barcodes_path):
        size = os.path.getsize(barcodes_path)
        result['checks']['barcodes_tsv_gz'] = size
        if size <= THRESHOLDS['barcodes_tsv_gz']:
            result['error'] = (f"barcodes.tsv.gz is {size} bytes "
                             f"(threshold: >{THRESHOLDS['barcodes_tsv_gz']})")
            return result
    else:
        result['checks']['barcodes_tsv_gz'] = 0
        result['error'] = "barcodes.tsv.gz not found"
        return result
    
    # ---- Check 2: filtered_feature_bc_matrix.h5 > 1000 bytes ----
    h5_path = os.path.join(outs_dir, 'filtered_feature_bc_matrix.h5')
    if os.path.exists(h5_path):
        size = os.path.getsize(h5_path)
        result['checks']['filtered_matrix_h5'] = size
        if size <= THRESHOLDS['filtered_matrix_h5']:
            result['error'] = (f"filtered_feature_bc_matrix.h5 is {size} bytes "
                             f"(threshold: >{THRESHOLDS['filtered_matrix_h5']})")
            return result
    else:
        result['checks']['filtered_matrix_h5'] = 0
        result['error'] = "filtered_feature_bc_matrix.h5 not found"
        return result
    
    # ---- Check 3: possorted_genome_bam.bam > 10,000 bytes ----
    bam_path = os.path.join(outs_dir, 'possorted_genome_bam.bam')
    if os.path.exists(bam_path):
        size = os.path.getsize(bam_path)
        result['checks']['bam'] = size
        if size <= THRESHOLDS['bam']:
            result['error'] = (f"possorted_genome_bam.bam is {size} bytes "
                             f"(threshold: >{THRESHOLDS['bam']})")
            return result
    else:
        result['checks']['bam'] = 0
        result['error'] = "possorted_genome_bam.bam not found"
        return result
    
    # ---- All checks passed ----
    result['valid'] = True
    return result


def update_sample_pickle(sample_info_path, validated_ids, failed_details):
    """
    Back up the original pickle and write a filtered version with only
    validated samples. Mirrors the standalone script's
    dictionary_file_filtered.pkl pattern.
    
    Parameters
    ----------
    sample_info_path : str
        Path to the sample info pickle file
    validated_ids : set
        Set of sample IDs that passed validation
    failed_details : dict
        Dict of sample_id -> failure info for logging
    """
    if not sample_info_path or not os.path.exists(sample_info_path):
        logger.warning("No sample_info pickle path provided — skipping pickle update")
        return

    # Load original to inspect its keying before touching anything on disk.
    with open(sample_info_path, 'rb') as f:
        raw_dict = pickle.load(f)

    # In multi (patient-level) mode the validated IDs are multi_ids (e.g. P18),
    # not SRR run accessions, so filtering the SRAscraper pickle by run_accession
    # would match nothing for every group and clobber the pickle with an empty
    # dict. Detect that keying mismatch and skip the update cleanly — the library
    # sheet governs multi-mode sample identity and validated_samples.txt is the
    # downstream gate, so the pickle is not consulted after the checkpoint.
    pickle_ids = set()
    for key, value in raw_dict.items():
        if hasattr(value, 'empty') and not value.empty and 'run_accession' in getattr(value, 'columns', []):
            pickle_ids.update(value['run_accession'].astype(str).tolist())
        elif isinstance(value, dict):
            pickle_ids.add(str(value.get('run_accession', key)))
        pickle_ids.add(str(key))

    if not ({str(v) for v in validated_ids} & pickle_ids):
        logger.info("Validated IDs are not run-accession keys in the sample pickle "
                    "(multi / patient-level mode); skipping pickle filtering to avoid "
                    "overwriting it. validated_samples.txt is the downstream gate.")
        return

    # Back up original (only if backup doesn't already exist)
    backup_path = sample_info_path.replace('.pkl', '_original.pkl')
    if not os.path.exists(backup_path):
        shutil.copy2(sample_info_path, backup_path)
        logger.info(f"Backed up original pickle to: {backup_path}")
    else:
        logger.info(f"Backup already exists: {backup_path}")
    
    # Load original
    with open(backup_path, 'rb') as f:
        raw_dict = pickle.load(f)
    
    # Filter to only validated samples
    filtered_dict = {}
    removed_count = 0
    
    for key, value in raw_dict.items():
        if hasattr(value, 'empty') and not value.empty:
            # DataFrame (SRAscraper format)
            if 'run_accession' in value.columns:
                filtered = value[value['run_accession'].isin(validated_ids)]
                if not filtered.empty:
                    filtered_dict[key] = filtered
                    removed = len(value) - len(filtered)
                    if removed > 0:
                        removed_count += removed
                        logger.info(f"  {key}: kept {len(filtered)}/{len(value)} samples")
                else:
                    removed_count += len(value)
                    logger.warning(f"  {key}: ALL samples removed (0/{len(value)} passed)")
            else:
                filtered_dict[key] = value
        elif isinstance(value, dict):
            # Dict format — check if this sample's run_accession is validated
            run_acc = value.get('run_accession', key)
            if run_acc in validated_ids or key in validated_ids:
                filtered_dict[key] = value
            else:
                removed_count += 1
        else:
            filtered_dict[key] = value

    # Safety guard: never overwrite the pickle with an empty dict. If filtering
    # collapsed everything (keying mismatch or all-failed), keep the original.
    if not filtered_dict:
        logger.warning("Filtered pickle would be empty; keeping the original pickle "
                       "unchanged (backup preserved).")
        return

    # Write filtered pickle (overwrites current, backup is safe)
    with open(sample_info_path, 'wb') as f:
        pickle.dump(filtered_dict, f)
    
    logger.info(f"Updated pickle: removed {removed_count} failed samples")
    logger.info(f"  Original backup: {backup_path}")
    logger.info(f"  Filtered pickle: {sample_info_path}")


# =============================================================================
# Main Validation Logic
# =============================================================================

def run_validation(status_files, cellranger_dir, sample_info_path,
                   output_validated, output_summary, log_file=None):
    """
    Run the full validation pipeline.
    
    Parameters
    ----------
    status_files : list
        Paths to all cellranger_status.json files
    cellranger_dir : str
        Base cellranger output directory
    sample_info_path : str
        Path to sample info pickle (for updating)
    output_validated : str
        Path to write validated_samples.txt
    output_summary : str
        Path to write validation_summary.json
    log_file : str, optional
        Path to log file
    """
    # Set up file logging
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    logger.info("="*60)
    logger.info("CELL RANGER OUTPUT VALIDATION (Checkpoint)")
    logger.info("="*60)
    logger.info(f"Validating {len(status_files)} samples")
    logger.info(f"Thresholds: barcodes>{THRESHOLDS['barcodes_tsv_gz']}B, "
                f"H5>{THRESHOLDS['filtered_matrix_h5']}B, "
                f"BAM>{THRESHOLDS['bam']}B")
    
    # =========================================================================
    # Phase 1: Read status JSONs
    # =========================================================================
    validated = []
    failed_cellranger = []
    failed_validation = []
    details = {}
    
    for status_file in status_files:
        try:
            with open(status_file, 'r') as f:
                status = json.load(f)
        except (json.JSONDecodeError, FileNotFoundError) as e:
            sample_id = Path(status_file).parent.name
            logger.error(f"  {sample_id}: Could not read status JSON — {e}")
            failed_cellranger.append(sample_id)
            details[sample_id] = {
                'status': 'failed_cellranger',
                'reason': f'Could not read status JSON: {e}',
            }
            continue
        
        sample_id = status.get('sample_id', Path(status_file).parent.name)
        
        # Check if Cell Ranger itself reported failure
        if not status.get('success', False):
            logger.warning(f"  {sample_id}: Cell Ranger failed — {status.get('error', 'unknown')}")
            failed_cellranger.append(sample_id)
            details[sample_id] = {
                'status': 'failed_cellranger',
                'reason': status.get('error', 'Cell Ranger reported failure'),
                'chemistry': status.get('chemistry'),
            }
            continue
        
        # =================================================================
        # Phase 2: Validate outputs on disk
        # =================================================================
        validation = validate_sample_outputs(cellranger_dir, sample_id)
        
        if validation['valid']:
            logger.info(f"  {sample_id}: PASSED "
                       f"(barcodes={validation['checks'].get('barcodes_tsv_gz', 0):,}B, "
                       f"H5={validation['checks'].get('filtered_matrix_h5', 0):,}B, "
                       f"BAM={validation['checks'].get('bam', 0):,}B)")
            validated.append(sample_id)
            details[sample_id] = {
                'status': 'passed',
                'chemistry': status.get('chemistry'),
                'file_sizes': validation['checks'],
            }
        else:
            logger.warning(f"  {sample_id}: FAILED VALIDATION — {validation['error']}")
            failed_validation.append(sample_id)
            details[sample_id] = {
                'status': 'failed_validation',
                'reason': validation['error'],
                'chemistry': status.get('chemistry'),
                'file_sizes': validation['checks'],
            }
    
    # =========================================================================
    # Phase 3: Write outputs
    # =========================================================================
    
    # Validated samples list (one per line — consumed by checkpoint input functions)
    os.makedirs(os.path.dirname(output_validated), exist_ok=True)
    with open(output_validated, 'w') as f:
        for sid in sorted(validated):
            f.write(f"{sid}\n")
    
    logger.info(f"Wrote {len(validated)} validated samples to: {output_validated}")
    
    # Validation summary JSON
    summary = {
        'timestamp': datetime.now().isoformat(),
        'total_samples': len(status_files),
        'passed': len(validated),
        'failed_cellranger': len(failed_cellranger),
        'failed_validation': len(failed_validation),
        'thresholds': THRESHOLDS,
        'validated_samples': sorted(validated),
        'failed_cellranger_samples': sorted(failed_cellranger),
        'failed_validation_samples': sorted(failed_validation),
        'details': details,
    }
    
    with open(output_summary, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Wrote validation summary to: {output_summary}")
    
    # Skipped samples TSV (convenient for review)
    all_failed = (
        [(sid, 'failed_cellranger', details[sid].get('reason', ''))
         for sid in failed_cellranger] +
        [(sid, 'failed_validation', details[sid].get('reason', ''))
         for sid in failed_validation]
    )
    
    if all_failed:
        skip_log = os.path.join(os.path.dirname(output_summary), 'skipped_samples.tsv')
        pd.DataFrame(all_failed, columns=['sample_id', 'failure_type', 'reason']).to_csv(
            skip_log, sep='\t', index=False)
        logger.info(f"Wrote skipped samples log to: {skip_log}")
    
    # =========================================================================
    # Phase 4: Update sample pickle
    # =========================================================================
    if sample_info_path:
        update_sample_pickle(sample_info_path, set(validated), details)
    
    # =========================================================================
    # Final Summary
    # =========================================================================
    logger.info("")
    logger.info("="*60)
    logger.info("VALIDATION SUMMARY")
    logger.info("="*60)
    logger.info(f"  Total samples:       {len(status_files)}")
    logger.info(f"  Passed validation:   {len(validated)}")
    logger.info(f"  Failed Cell Ranger:  {len(failed_cellranger)}")
    logger.info(f"  Failed validation:   {len(failed_validation)}")
    logger.info("")
    
    if failed_cellranger:
        logger.info("  Failed Cell Ranger samples:")
        for sid in failed_cellranger:
            logger.info(f"    {sid}: {details[sid].get('reason', 'unknown')}")
    
    if failed_validation:
        logger.info("  Failed validation samples (Cell Ranger 'succeeded' but outputs unusable):")
        for sid in failed_validation:
            logger.info(f"    {sid}: {details[sid].get('reason', 'unknown')}")
    
    if len(validated) == 0:
        logger.error("NO SAMPLES PASSED VALIDATION — pipeline cannot continue")
        raise ValueError("No samples passed Cell Ranger validation. "
                        "Check logs and validation_summary.json for details.")
    
    logger.info("")
    logger.info(f"Pipeline will continue with {len(validated)} validated samples")
    logger.info("="*60)


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == '__main__':
    # Check if running from Snakemake
    try:
        snakemake
        
        run_validation(
            status_files=snakemake.input.statuses,
            cellranger_dir=snakemake.params.cellranger_dir,
            sample_info_path=snakemake.params.sample_info_path,
            output_validated=snakemake.output.validated,
            output_summary=snakemake.output.summary,
            log_file=snakemake.log[0] if snakemake.log else None,
        )
        
    except NameError:
        print("This script is intended to be called from Snakemake.")
        print("Use: checkpoint validate_cellranger in your Snakefile")
        sys.exit(1)
