#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kraken2_viral_detection.py
==========================

Detect viral/microbial sequences in unmapped single-cell reads using Kraken2.

This script:
1. Extracts unmapped reads from Cell Ranger BAM files
2. Runs Kraken2 classification on unmapped reads
3. Builds single-cell count matrices linking barcodes to detected organisms
4. Generates summary reports

The output is a sparse matrix similar to Cell Ranger's filtered_feature_bc_matrix
but with organisms instead of genes, allowing integration with scanpy.

Requirements:
    - Kraken2 installed and database built (see README for setup)
    - samtools
    - pysam

Usage:
    Called via Snakemake rule with snakemake.input/output/params
    
    Or standalone:
    python kraken2_viral_detection.py --bam input.bam --db /path/to/kraken2_db --output outdir
"""

import os
import sys
import gzip
import shutil
import logging
import argparse
import subprocess
import collections.abc
from pathlib import Path

import pysam
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import csv

# Try to import regex, fall back to re
try:
    import regex as re
except ImportError:
    import re

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Core Functions
# =============================================================================

def extract_unmapped_reads(bam_path, output_bam, output_fastq, threads=8):
    """
    Extract unmapped reads from a BAM file.
    
    Parameters
    ----------
    bam_path : str
        Path to input BAM file (possorted_genome_bam.bam)
    output_bam : str
        Path for output BAM with unmapped reads
    output_fastq : str
        Path for output FASTQ file
    threads : int
        Number of threads for samtools
        
    Returns
    -------
    int
        Number of unmapped reads extracted
    """
    logger.info(f"Extracting unmapped reads from {bam_path}")
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_bam), exist_ok=True)
    
    # Extract unmapped reads to BAM
    cmd1 = [
        "samtools", "view",
        "-@", str(threads),
        "-b", "-f", "4",  # -f 4 means unmapped
        bam_path,
        "-o", output_bam
    ]
    
    logger.info(f"  Running: {' '.join(cmd1)}")
    result = subprocess.run(cmd1, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"samtools view failed: {result.stderr}")
        raise RuntimeError(f"Failed to extract unmapped reads: {result.stderr}")
    
    # Convert to FASTQ
    cmd2 = [
        "samtools", "fastq",
        "-@", str(threads),
        "-n",  # Don't append /1 or /2 to read names
        output_bam
    ]
    
    logger.info(f"  Running: {' '.join(cmd2[:4])}...")
    with open(output_fastq, 'w') as f:
        result = subprocess.run(cmd2, stdout=f, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        logger.error(f"samtools fastq failed: {result.stderr}")
        raise RuntimeError(f"Failed to convert to FASTQ: {result.stderr}")
    
    # Count reads
    n_reads = sum(1 for _ in pysam.AlignmentFile(output_bam, "rb"))
    logger.info(f"  Extracted {n_reads:,} unmapped reads")
    
    return n_reads


def run_kraken2(fastq_path, db_path, output_file, report_file, threads=8, confidence=0.0):
    """
    Run Kraken2 classification on FASTQ file.
    
    Parameters
    ----------
    fastq_path : str
        Path to input FASTQ file
    db_path : str
        Path to Kraken2 database
    output_file : str
        Path for Kraken2 output file
    report_file : str
        Path for Kraken2 report file
    threads : int
        Number of threads
    confidence : float
        Confidence threshold for Kraken2
        
    Returns
    -------
    dict
        Summary statistics from the run
    """
    logger.info(f"Running Kraken2 with database: {db_path}")
    
    cmd = [
        "kraken2",
        "--threads", str(threads),
        "--db", db_path,
        "--confidence", str(confidence),
        "--report", report_file,
        "--output", output_file,
        fastq_path
    ]
    
    logger.info(f"  Running: kraken2 --db {db_path} ...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Kraken2 failed: {result.stderr}")
        raise RuntimeError(f"Kraken2 classification failed: {result.stderr}")
    
    # Parse summary from stderr (Kraken2 outputs stats there)
    stats = {}
    for line in result.stderr.split('\n'):
        if 'sequences classified' in line.lower():
            try:
                parts = line.strip().split()
                stats['classified'] = int(parts[0])
            except:
                pass
        elif 'sequences unclassified' in line.lower():
            try:
                parts = line.strip().split()
                stats['unclassified'] = int(parts[0])
            except:
                pass
    
    logger.info(f"  Kraken2 completed")
    if stats:
        logger.info(f"  Classified: {stats.get('classified', 'N/A')}, Unclassified: {stats.get('unclassified', 'N/A')}")
    
    return stats


def extract_ids(bamfile, krakenfile):
    """
    Build nested dictionary with Kraken2 taxonomy code for each transcript and cell.
    
    Parameters
    ----------
    bamfile : str
        BAM file with unmapped reads (has CB and UB tags)
    krakenfile : str
        Kraken2 output file
        
    Returns
    -------
    dict
        Nested dictionary {cellbarcode: {transcriptbarcode: [taxonomy_ids]}}
    """
    logger.info("Extracting cell barcode and taxonomy associations...")
    
    line = 0
    skipped = 0
    nested_dict = {}
    
    # Iterate simultaneously through bam and kraken file
    bam_iter = pysam.AlignmentFile(bamfile, "rb")
    
    with open(krakenfile, "r") as kf:
        for sread, kread in zip(bam_iter, kf):
            line += 1
            
            # Check that read names match
            kread_parts = kread.split('\t')
            if sread.query_name != kread_parts[1]:
                skipped += 1
                logger.debug(f"Read name mismatch: BAM={sread.query_name}, Kraken={kread_parts[1]}")
                continue
            
            # Get cell barcode and UMI from BAM file
            try:
                sread_CB = sread.get_tag('CB')
                sread_UB = sread.get_tag('UB')
            except KeyError:
                # Some reads don't have cell barcode or UMI
                skipped += 1
                continue
            
            # Get taxonomy ID from kraken file
            kread_taxid = kread_parts[2]
            if not kread_taxid.isdigit():
                try:
                    # Sometimes format is "name (taxid #)"
                    match = re.search(r'\(taxid\s*(\d+)\)', kread_taxid)
                    if match:
                        kread_taxid = match.group(1)
                    else:
                        # Try to extract just numbers
                        match = re.search(r'(\d+)', kread_taxid)
                        if match:
                            kread_taxid = match.group(1)
                        else:
                            skipped += 1
                            continue
                except:
                    logger.debug(f"Could not parse taxonomy ID: {kread_taxid}")
                    skipped += 1
                    continue
            
            # Build nested dictionary
            if sread_CB not in nested_dict:
                nested_dict[sread_CB] = {}
            
            if sread_UB not in nested_dict[sread_CB]:
                nested_dict[sread_CB][sread_UB] = []
            
            nested_dict[sread_CB][sread_UB].append(kread_taxid)
    
    bam_iter.close()
    
    logger.info(f"  Total reads: {line:,}, Skipped: {skipped:,}")
    logger.info(f"  Unique cells with hits: {len(nested_dict):,}")
    
    return nested_dict


def most_frequent(lst):
    """Find the most frequent element in a list."""
    if not lst:
        return '0'
    return max(set(lst), key=lst.count)


def map_nested_dicts(ob, func):
    """Apply a function to innermost items of nested dictionaries."""
    for k, v in ob.items():
        if isinstance(v, collections.abc.Mapping):
            map_nested_dicts(v, func)
        else:
            ob[k] = func(v)


def twist_dict(nested):
    """
    Convert to count dictionary format.
    
    Parameters
    ----------
    nested : dict
        {cellbarcode: {transcriptbarcode: taxonomy_id}}
        
    Returns
    -------
    dict
        {cellbarcode: {taxonomy_id: count}}
    """
    newdict = {}
    for ckey, tdict in nested.items():
        for tkey, kvalue in tdict.items():
            if ckey not in newdict:
                newdict[ckey] = {}
            
            if kvalue not in newdict[ckey]:
                newdict[ckey][kvalue] = 0
            newdict[ckey][kvalue] += 1
    
    return newdict


def dict2lists(nested):
    """
    Convert nested dictionary to sparse matrix format.
    
    Returns
    -------
    tuple
        (rows, columns, values, cell_list, taxid_list)
    """
    rows = []      # taxonomy coordinate
    columns = []   # cell coordinate
    values = []    # count
    
    cell_list = []
    taxid_list = []
    
    j = 0
    for ckey, taxdict in nested.items():
        for taxkey, count in taxdict.items():
            try:
                k = taxid_list.index(taxkey)
            except ValueError:
                taxid_list.append(taxkey)
                k = len(taxid_list) - 1
            
            rows.append(k)
            columns.append(j)
            values.append(count)
        
        cell_list.append(ckey)
        j += 1
    
    return rows, columns, values, cell_list, taxid_list


def krakenID2dict(dbfile, taxid_list):
    """
    Get taxonomy names from Kraken2 database inspect file.
    
    Parameters
    ----------
    dbfile : str
        Path to Kraken2 inspect.txt file
    taxid_list : list
        List of taxonomy IDs to look up
        
    Returns
    -------
    dict
        {taxid: taxonomy_name}
    """
    taxdict = {'0': 'unclassified'}
    taxid_set = set(taxid_list)
    
    if not os.path.exists(dbfile):
        logger.warning(f"Database inspect file not found: {dbfile}")
        # Return just IDs as names
        return {tid: f"taxid_{tid}" for tid in taxid_list}
    
    with open(dbfile) as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                taxid_db = parts[4]
                taxname = parts[5].strip()
                
                if taxid_db in taxid_set:
                    taxdict[taxid_db] = taxname
    
    return taxdict


def filter_organisms(count_dict, taxdict, include_list=None, exclude_list=None):
    """
    Filter organisms based on inclusion/exclusion lists.
    
    Parameters
    ----------
    count_dict : dict
        {cellbarcode: {taxonomy_id: count}}
    taxdict : dict
        {taxid: taxonomy_name}
    include_list : list, optional
        List of organism name patterns to include (if provided, only these are kept)
    exclude_list : list, optional
        List of organism name patterns to exclude
        
    Returns
    -------
    dict
        Filtered count dictionary
    """
    if include_list is None and exclude_list is None:
        return count_dict
    
    # Build set of taxonomy IDs to keep
    keep_taxids = set()
    
    for taxid, taxname in taxdict.items():
        keep = True
        taxname_lower = taxname.lower()
        
        if include_list:
            # Only keep if matches any include pattern
            keep = any(pattern.lower() in taxname_lower for pattern in include_list)
        
        if exclude_list and keep:
            # Exclude if matches any exclude pattern
            keep = not any(pattern.lower() in taxname_lower for pattern in exclude_list)
        
        if keep:
            keep_taxids.add(taxid)
    
    # Filter count dictionary
    filtered_dict = {}
    for cell, tax_counts in count_dict.items():
        filtered_counts = {tid: cnt for tid, cnt in tax_counts.items() if tid in keep_taxids}
        if filtered_counts:
            filtered_dict[cell] = filtered_counts
    
    logger.info(f"  Filtered to {len(keep_taxids)} organisms, {len(filtered_dict)} cells with hits")
    
    return filtered_dict


def build_sparse_matrix(bamfile, krakenfile, dbpath, outdir, 
                        include_organisms=None, exclude_organisms=None):
    """
    Build sparse matrix with transcript counts per organism for each cell.
    
    Parameters
    ----------
    bamfile : str
        BAM file with unmapped reads
    krakenfile : str
        Kraken2 output file
    dbpath : str
        Path to Kraken2 database directory
    outdir : str
        Output directory for matrix files
    include_organisms : list, optional
        List of organism patterns to include
    exclude_organisms : list, optional
        List of organism patterns to exclude
        
    Returns
    -------
    dict
        Summary statistics
    """
    logger.info("Building single-cell sparse matrix...")
    
    # Create output directory - this is the matrix_dir output
    matrix_dir = outdir
    os.makedirs(matrix_dir, exist_ok=True)
    
    # Output file paths
    matrixfile = os.path.join(matrix_dir, 'matrix.mtx')
    cellfile = os.path.join(matrix_dir, 'barcodes.tsv')
    taxfile = os.path.join(matrix_dir, 'features.tsv')
    
    # Database inspect file
    inspectfile = os.path.join(dbpath, 'inspect.txt')
    
    # Extract taxonomy IDs for each transcript
    mg_dict = extract_ids(bamfile, krakenfile)
    
    if not mg_dict:
        logger.warning("No valid cell-taxonomy associations found!")
        # Create empty files
        with open(matrixfile, 'w') as f:
            f.write("%%MatrixMarket matrix coordinate integer general\n%\n0 0 0\n")
        with open(cellfile, 'w') as f:
            pass
        with open(taxfile, 'w') as f:
            pass
        # Gzip empty files
        for filename in ['matrix.mtx', 'barcodes.tsv', 'features.tsv']:
            filepath = os.path.join(matrix_dir, filename)
            with open(filepath, 'rb') as f_in:
                with gzip.open(filepath + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(filepath)
        return {'cells': 0, 'organisms': 0, 'total_counts': 0}
    
    # Find most frequent taxonomy for each transcript (UMI deduplication)
    map_nested_dicts(mg_dict, most_frequent)
    
    # Convert to count dictionary
    count_dict = twist_dict(mg_dict)
    
    # Get all unique taxonomy IDs
    all_taxids = set()
    for tax_counts in count_dict.values():
        all_taxids.update(tax_counts.keys())
    
    # Get taxonomy names
    taxdict = krakenID2dict(inspectfile, list(all_taxids))
    
    # Filter organisms if specified
    if include_organisms or exclude_organisms:
        count_dict = filter_organisms(count_dict, taxdict, include_organisms, exclude_organisms)
    
    if not count_dict:
        logger.warning("No data remaining after filtering!")
        return {'cells': 0, 'organisms': 0, 'total_counts': 0}
    
    # Convert to sparse matrix format
    rows, cols, vals, cell_list, taxid_list = dict2lists(count_dict)
    
    # Get taxonomy names in order
    taxname_list = [taxdict.get(k, f"taxid_{k}") for k in taxid_list]
    
    # Create and save sparse matrix
    sparsematrix = csr_matrix((vals, (rows, cols)), 
                               shape=(len(taxid_list), len(cell_list)))
    
    mmwrite(matrixfile, sparsematrix)
    
    # Save cell barcodes
    with open(cellfile, 'w') as f:
        for cell in cell_list:
            f.write(f"{cell}\n")
    
    # Save features (taxonomy IDs and names)
    with open(taxfile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for taxid, taxname in zip(taxid_list, taxname_list):
            writer.writerow([taxid, taxname, 'Taxonomy'])
    
    # Copy database hierarchy if available
    if os.path.exists(inspectfile):
        shutil.copy(inspectfile, os.path.join(matrix_dir, 'hierarchy.txt'))
    
    # Gzip the files to match Cell Ranger format
    for filename in ['matrix.mtx', 'barcodes.tsv', 'features.tsv']:
        filepath = os.path.join(matrix_dir, filename)
        if os.path.exists(filepath):
            with open(filepath, 'rb') as f_in:
                with gzip.open(filepath + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(filepath)
    
    stats = {
        'cells': len(cell_list),
        'organisms': len(taxid_list),
        'total_counts': sum(vals),
        'matrix_dir': matrix_dir
    }
    
    logger.info(f"  Created matrix: {stats['cells']} cells x {stats['organisms']} organisms")
    logger.info(f"  Total UMI counts: {stats['total_counts']:,}")
    
    return stats


def generate_summary_report(kraken_report, taxdict, output_file, top_n=50):
    """
    Generate a summary report of detected organisms.
    
    Parameters
    ----------
    kraken_report : str
        Path to Kraken2 report file
    taxdict : dict
        {taxid: taxonomy_name}
    output_file : str
        Path for output summary TSV
    top_n : int
        Number of top organisms to include in summary
    """
    logger.info("Generating summary report...")
    
    # Parse Kraken2 report
    report_data = []
    
    with open(kraken_report, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                pct = float(parts[0])
                num_clade = int(parts[1])
                num_taxon = int(parts[2])
                rank = parts[3]
                taxid = parts[4]
                name = parts[5].strip()
                
                if num_taxon > 0:  # Only include taxa with direct hits
                    report_data.append({
                        'taxonomy_id': taxid,
                        'name': name,
                        'rank': rank,
                        'reads_direct': num_taxon,
                        'reads_clade': num_clade,
                        'percentage': pct
                    })
    
    # Sort by direct reads and take top N
    report_data.sort(key=lambda x: x['reads_direct'], reverse=True)
    top_organisms = report_data[:top_n]
    
    # Save summary
    df = pd.DataFrame(top_organisms)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"  Summary saved to {output_file}")
    logger.info(f"  Top organisms detected:")
    for org in top_organisms[:10]:
        logger.info(f"    {org['name']}: {org['reads_direct']:,} reads")


def process_sample(
    bam_path,
    db_path,
    output_dir,
    matrix_dir,
    summary_file,
    sample_id,
    threads=8,
    confidence=0.0,
    include_organisms=None,
    exclude_organisms=None
):
    """
    Process a single sample through the complete Kraken2 pipeline.
    
    Parameters
    ----------
    bam_path : str
        Path to Cell Ranger BAM file
    db_path : str
        Path to Kraken2 database
    output_dir : str
        Base output directory for this sample
    matrix_dir : str
        Output directory for matrix files (Snakemake output)
    summary_file : str
        Path for summary TSV file (Snakemake output)
    sample_id : str
        Sample identifier
    threads : int
        Number of threads
    confidence : float
        Kraken2 confidence threshold
    include_organisms : list, optional
        Organisms to include
    exclude_organisms : list, optional
        Organisms to exclude
        
    Returns
    -------
    dict
        Processing results and statistics
    """
    logger.info(f"\n{'='*60}")
    logger.info(f"Processing sample: {sample_id}")
    logger.info(f"{'='*60}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Define intermediate file paths
    unmapped_bam = os.path.join(output_dir, f"{sample_id}_unmapped.bam")
    unmapped_fq = os.path.join(output_dir, f"{sample_id}_unmapped.fq")
    kraken_output = os.path.join(output_dir, f"{sample_id}_kraken.out")
    kraken_report = os.path.join(output_dir, f"{sample_id}_kraken_report.txt")
    
    results = {
        'sample_id': sample_id,
        'success': False
    }
    
    try:
        # Step 1: Extract unmapped reads
        n_unmapped = extract_unmapped_reads(bam_path, unmapped_bam, unmapped_fq, threads)
        results['unmapped_reads'] = n_unmapped
        
        if n_unmapped == 0:
            logger.warning(f"No unmapped reads found for {sample_id}")
            # Create empty outputs
            os.makedirs(matrix_dir, exist_ok=True)
            with open(os.path.join(matrix_dir, 'matrix.mtx'), 'w') as f:
                f.write("%%MatrixMarket matrix coordinate integer general\n%\n0 0 0\n")
            for filename in ['barcodes.tsv', 'features.tsv']:
                with open(os.path.join(matrix_dir, filename), 'w') as f:
                    pass
            for filename in ['matrix.mtx', 'barcodes.tsv', 'features.tsv']:
                filepath = os.path.join(matrix_dir, filename)
                with open(filepath, 'rb') as f_in:
                    with gzip.open(filepath + '.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(filepath)
            pd.DataFrame().to_csv(summary_file, sep='\t', index=False)
            results['success'] = True
            results['cells'] = 0
            results['organisms'] = 0
            return results
        
        # Step 2: Run Kraken2
        kraken_stats = run_kraken2(unmapped_fq, db_path, kraken_output, kraken_report, 
                                    threads, confidence)
        results.update(kraken_stats)
        
        # Step 3: Build sparse matrix (output to matrix_dir)
        matrix_stats = build_sparse_matrix(
            unmapped_bam, kraken_output, db_path, matrix_dir,
            include_organisms, exclude_organisms
        )
        results.update(matrix_stats)
        
        # Step 4: Generate summary report
        inspectfile = os.path.join(db_path, 'inspect.txt')
        taxdict = {}
        if os.path.exists(inspectfile):
            taxdict = krakenID2dict(inspectfile, [])
        generate_summary_report(kraken_report, taxdict, summary_file)
        
        results['success'] = True
        
        # Cleanup intermediate files (optional)
        # os.remove(unmapped_fq)
        
    except Exception as e:
        logger.error(f"Error processing {sample_id}: {e}")
        import traceback
        traceback.print_exc()
        results['error'] = str(e)
    
    return results


# =============================================================================
# Snakemake Integration
# =============================================================================

def run_from_snakemake():
    """Run from Snakemake rule."""
    
    # Get inputs
    bam_path = snakemake.input.bam
    
    # Get outputs - these are the required Snakemake outputs
    summary_file = snakemake.output.summary
    matrix_dir = snakemake.output.matrix_dir
    
    # Get params
    sample_id = snakemake.params.sample_id
    db_path = snakemake.params.db_path
    output_dir = snakemake.params.output_dir
    threads = snakemake.threads
    confidence = snakemake.params.get('confidence', 0.0)
    include_organisms = snakemake.params.get('include_organisms')
    exclude_organisms = snakemake.params.get('exclude_organisms')
    
    # Build sample-specific output directory
    sample_output_dir = os.path.join(output_dir, sample_id)
    
    # Set up logging
    if snakemake.log:
        os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
        file_handler = logging.FileHandler(snakemake.log[0])
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Run processing
    results = process_sample(
        bam_path=bam_path,
        db_path=db_path,
        output_dir=sample_output_dir,
        matrix_dir=matrix_dir,
        summary_file=summary_file,
        sample_id=sample_id,
        threads=threads,
        confidence=confidence,
        include_organisms=include_organisms,
        exclude_organisms=exclude_organisms
    )
    
    if not results['success']:
        sys.exit(1)


# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Detect viral/microbial sequences in unmapped single-cell reads'
    )
    parser.add_argument('--bam', required=True, help='Input BAM file (Cell Ranger output)')
    parser.add_argument('--db', required=True, help='Path to Kraken2 database')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--sample-id', default='sample', help='Sample identifier')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads')
    parser.add_argument('--confidence', type=float, default=0.0, help='Kraken2 confidence threshold')
    parser.add_argument('--include', nargs='*', help='Organism patterns to include')
    parser.add_argument('--exclude', nargs='*', help='Organism patterns to exclude')
    
    args = parser.parse_args()
    
    # For CLI, create standard output structure
    matrix_dir = os.path.join(args.output, 'kraken2_filtered_feature_bc_matrix')
    summary_file = os.path.join(args.output, f'{args.sample_id}_organism_summary.tsv')
    
    results = process_sample(
        bam_path=args.bam,
        db_path=args.db,
        output_dir=args.output,
        matrix_dir=matrix_dir,
        summary_file=summary_file,
        sample_id=args.sample_id,
        threads=args.threads,
        confidence=args.confidence,
        include_organisms=args.include,
        exclude_organisms=args.exclude
    )
    
    if results['success']:
        logger.info("\nProcessing completed successfully!")
        logger.info(f"  Cells with detections: {results.get('cells', 0)}")
        logger.info(f"  Organisms detected: {results.get('organisms', 0)}")
    else:
        logger.error(f"\nProcessing failed: {results.get('error', 'Unknown error')}")
        sys.exit(1)


if __name__ == '__main__':
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        main()
