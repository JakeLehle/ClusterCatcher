#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scomatic_mutation_calling.py
============================

Comprehensive SComatic pipeline for single-cell somatic mutation detection.

This script orchestrates the SComatic workflow:
1. Filter BAM to annotated cells only
2. Split BAM by cell type (SplitBamCellTypes.py)
3. Count bases per cell (BaseCellCounter.py)
4. Merge counts across cell types (MergeBaseCellCounts.py)
5. Call variants - Step 1 (BaseCellCalling.step1.py)
6. Call variants - Step 2 with filtering (BaseCellCalling.step2.py)
7. Filter with BED file for mappable regions
8. Compute callable sites per cell
9. Generate single-cell genotypes
10. Filter and annotate genotypes with trinucleotide context

Usage:
    Called by Snakemake or run standalone with arguments.

Author: Jake Lehle
Date: 2025
"""

import os
import sys
import yaml
import pandas as pd
import scanpy as sc
from tqdm import tqdm
import pysam
import subprocess
import shutil
import glob
import pickle
import threading
import time
import psutil
import numpy as np
import argparse
import logging
from multiprocessing.pool import ThreadPool
from math import ceil
from time import perf_counter
from os import cpu_count
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def create_filtered_bam(input_bam, output_bam, valid_barcodes):
    """
    Create filtered BAM file containing only specified barcodes.
    Handles both raw barcodes (from adata) and CB tags with suffix (e.g., '-1')
    """
    logger.info(f"Creating filtered BAM: {output_bam}")
    logger.info(f"Valid barcodes count: {len(valid_barcodes)}")
    
    kept_reads = 0
    total_reads = 0
    
    try:
        with pysam.AlignmentFile(input_bam, "rb") as inbam:
            with pysam.AlignmentFile(output_bam, "wb", template=inbam) as outbam:
                for read in tqdm(inbam, desc="Filtering BAM", miniters=1000):
                    total_reads += 1
                    
                    if not read.has_tag('CB'):
                        continue
                    
                    cb_tag = read.get_tag('CB')
                    raw_barcode = cb_tag.split('-')[0]
                    
                    if raw_barcode in valid_barcodes:
                        outbam.write(read)
                        kept_reads += 1
        
        logger.info(f"Total reads processed: {total_reads:,}")
        logger.info(f"Reads kept: {kept_reads:,}")
        logger.info(f"Reads filtered out: {total_reads - kept_reads:,}")
        
        if kept_reads == 0:
            logger.warning("No reads were kept in the filtered BAM!")
            return False
        else:
            logger.info("Indexing filtered BAM...")
            pysam.index(output_bam)
            return True
            
    except Exception as e:
        logger.error(f"CRITICAL ERROR in create_filtered_bam: {str(e)}")
        return False


def run_scomatic(args):
    """Run SComatic SplitBam for a single sample"""
    filtered_bam, annotation_file, run_acc, scomatic_dir, scomatic_scripts_dir = args
    cmd = [
        "python",
        os.path.join(scomatic_scripts_dir, "SplitBam", "SplitBamCellTypes.py"),
        "--bam", filtered_bam,
        "--meta", annotation_file,
        "--id", run_acc,
        "--max_nM", "5",
        "--max_NH", "1",
        "--min_MQ", "255",
        "--outdir", scomatic_dir
    ]
    try:
        logger.info(f"Running SComatic for {run_acc}...")
        result = subprocess.run(
            cmd, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True
        )
        logger.info(f"Completed SComatic for {run_acc}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"SComatic failed for {run_acc}: {e.returncode}")
        logger.error(f"Error message: {e.stderr[:500]}...")
        return False
    except Exception as e:
        logger.error(f"UNEXPECTED ERROR in run_scomatic: {str(e)}")
        return False


def run_basecell_counter(args):
    """Run BaseCellCounter with comprehensive validation"""
    bam_file, ref_genome, output_dir, scomatic_scripts_dir = args
    cell_type = os.path.basename(bam_file).split('.')[-2]
    temp_dir = os.path.join(output_dir, f'temp_{cell_type}')
    
    try:
        if not os.path.exists(bam_file) or os.path.getsize(bam_file) == 0:
            logger.error(f"BAM file missing or empty: {bam_file}")
            return False
            
        if not os.path.exists(ref_genome) or os.path.getsize(ref_genome) == 0:
            logger.error(f"Reference genome missing or empty: {ref_genome}")
            return False
        
        os.makedirs(temp_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "BaseCellCounter", "BaseCellCounter.py"),
            "--bam", bam_file,
            "--ref", ref_genome,
            "--chrom", "all",
            "--out_folder", output_dir,
            "--min_bq", "30",
            "--tmp_dir", temp_dir,
            "--nprocs", "1"
        ]
        
        logger.info(f"Running BaseCellCounter for {cell_type}...")
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        output_files = glob.glob(os.path.join(output_dir, "*.tsv"))
        if not output_files:
            logger.error(f"No output files created for {cell_type}")
            return False
            
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"BaseCellCounter failed for {cell_type}: {e.returncode}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error processing {cell_type}: {str(e)}")
        return False
    finally:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)


def run_basecell_counter_parallel(basecell_counter_args, n_workers):
    """Run BaseCellCounter in parallel"""
    if not basecell_counter_args:
        logger.warning("No BAM files to process in BaseCellCounter")
        return []
    
    logger.info("="*80)
    logger.info("Running BaseCellCounter")
    logger.info(f"Processing {len(basecell_counter_args)} BAM files with {n_workers} workers")
    logger.info("="*80)
    
    start = perf_counter()
    
    # Validate all BAM files exist
    missing_bams = [bam for bam, *_ in basecell_counter_args if not os.path.exists(bam)]
    if missing_bams:
        logger.error(f"Missing BAM files: {len(missing_bams)}")
        for bam in missing_bams[:5]:
            logger.error(f"- {bam}")
        raise FileNotFoundError(f"{len(missing_bams)} BAM files not found")
    
    with ThreadPool(n_workers) as pool:
        chunksize = ceil(len(basecell_counter_args) / n_workers)
        results = list(pool.map(run_basecell_counter, basecell_counter_args, chunksize=chunksize))
    
    end = perf_counter()
    
    success_count = sum(results)
    failure_count = len(results) - success_count
    duration = end - start
    
    logger.info("="*80)
    logger.info("BASE CELL COUNTER SUMMARY:")
    logger.info(f"Total BAMs processed: {len(results)}")
    logger.info(f"Successfully processed: {success_count}")
    logger.info(f"Failures: {failure_count}")
    logger.info(f"Time elapsed: {duration:.2f} seconds")
    logger.info("="*80)
    
    return results


def run_merge_counts(tsv_folder, output_file, scomatic_scripts_dir, sample_id):
    """Run MergeBaseCellCounts.py"""
    try:
        tsv_files = glob.glob(os.path.join(tsv_folder, "*.tsv"))
        if not tsv_files:
            logger.warning(f"No TSV files found in {tsv_folder}")
            return False
            
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "MergeCounts", "MergeBaseCellCounts.py"),
            "--tsv_folder", tsv_folder,
            "--outfile", output_file
        ]
        logger.info(f"Merging {len(tsv_files)} count files from {tsv_folder}...")
        
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            logger.error(f"Merged file not created or empty: {output_file}")
            return False
            
        logger.info(f"Merged counts saved to {output_file}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to merge counts: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error in merge_counts: {str(e)}")
        return False


def run_variant_calling_step1(args):
    """Run BaseCellCalling step 1"""
    merged_counts, output_prefix, ref_genome, scomatic_scripts_dir = args
    step1_output = f"{output_prefix}.calling.step1.tsv"
    
    try:
        if not os.path.exists(merged_counts) or os.path.getsize(merged_counts) == 0:
            raise ValueError(f"Input file missing or empty: {merged_counts}")
            
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "BaseCellCalling", "BaseCellCalling.step1.py"),
            "--infile", merged_counts,
            "--outfile", output_prefix,
            "--ref", ref_genome
        ]
        logger.info(f"Running variant calling step1 for {output_prefix}...")
        
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if not os.path.exists(step1_output) or os.path.getsize(step1_output) == 0:
            raise RuntimeError(f"Output file missing or empty: {step1_output}")
            
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Variant calling step1 failed: {e.returncode}")
        return False
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return False


def run_variant_calling_step2(args):
    """Run BaseCellCalling step 2"""
    step1_output, output_prefix, editing_sites, pon_file, scomatic_scripts_dir = args
    step2_output = f"{output_prefix}.calling.step2.tsv"
    
    try:
        if not os.path.exists(step1_output) or os.path.getsize(step1_output) == 0:
            raise ValueError(f"Step1 output missing or empty: {step1_output}")
            
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "BaseCellCalling", "BaseCellCalling.step2.py"),
            "--infile", step1_output,
            "--outfile", output_prefix,
            "--editing", editing_sites,
            "--pon", pon_file
        ]
        logger.info(f"Running variant calling step2 for {output_prefix}...")
        
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if not os.path.exists(step2_output) or os.path.getsize(step2_output) == 0:
            raise RuntimeError(f"Output file missing or empty: {step2_output}")
            
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Variant calling step2 failed: {e.returncode}")
        return False
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return False


def filter_with_bed(input_tsv, bed_file, output_tsv):
    """Filter variants using bedtools intersect"""
    try:
        if not os.path.exists(input_tsv) or os.path.getsize(input_tsv) == 0:
            logger.error(f"Input file missing or empty: {input_tsv}")
            return False
            
        if not os.path.exists(bed_file) or os.path.getsize(bed_file) == 0:
            logger.error(f"BED file missing or empty: {bed_file}")
            return False
            
        cmd = f"bedtools intersect -header -a {input_tsv} -b {bed_file} | " \
              f"awk '$1 ~ /^#/ || $6 == \"PASS\"' > {output_tsv}"
              
        logger.info(f"Filtering {input_tsv} with BED file...")
        subprocess.run(cmd, shell=True, check=True)
        
        if not os.path.exists(output_tsv) or os.path.getsize(output_tsv) == 0:
            logger.warning(f"Output file empty after filtering: {output_tsv}")
            with open(output_tsv, 'w') as f:
                f.write("## Empty after filtering\n")
                
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"BED filtering failed: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error in filter_with_bed: {str(e)}")
        return False


def run_callable_sites(args):
    """Compute callable sites per cell type"""
    step1_file, output_prefix, max_cov, min_cell_types, scomatic_scripts_dir = args
    try:
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "GetCallableSites", "GetAllCallableSites.py"),
            "--infile", step1_file,
            "--outfile", output_prefix,
            "--max_cov", str(max_cov),
            "--min_cell_types", str(min_cell_types)
        ]
        logger.info(f"Computing callable sites for {output_prefix}...")
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Callable sites computation failed: {e}")
        return False


def run_sites_per_cell(args):
    """Compute callable sites per cell"""
    bam_file, step1_file, ref_genome, output_dir, temp_dir, scomatic_scripts_dir = args
    cell_type = os.path.basename(bam_file).split('.')[-2]
    output_file = os.path.join(output_dir, f"{os.path.basename(bam_file)}.SitesPerCell.tsv")
    
    try:
        if not os.path.exists(bam_file) or os.path.getsize(bam_file) == 0:
            logger.warning(f"BAM file missing or empty: {bam_file}")
            return False
            
        if not os.path.exists(step1_file) or os.path.getsize(step1_file) == 0:
            logger.warning(f"Step1 file missing or empty: {step1_file}")
            return False

        os.makedirs(temp_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

        if os.path.exists(output_file):
            os.remove(output_file)
        
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "SitesPerCell", "SitesPerCell.py"),
            "--bam", bam_file,
            "--infile", step1_file,
            "--ref", ref_genome,
            "--out_folder", output_dir,
            "--tmp_dir", temp_dir,
            "--nprocs", "1"
        ]
        logger.info(f"Processing sites per cell for {cell_type}...")
        
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if not os.path.exists(output_file):
            logger.error(f"Output file not created: {output_file}")
            return False
            
        if os.path.getsize(output_file) < 20:
            logger.warning(f"Output file empty, creating placeholder: {output_file}")
            with open(output_file, 'w') as f:
                f.write("CB,SitesPerCell\n")
                
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed for {cell_type}: {e.returncode}")
        return False
    except Exception as e:
        logger.error(f"Processing failed for {cell_type}: {str(e)}")
        return False
    finally:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)


def run_single_cell_genotype(args):
    """Compute genotypes for each cell"""
    bam_file, variant_file, meta_file, ref_genome, output_file, temp_dir, scomatic_scripts_dir, custom_script = args
    cell_type = os.path.basename(bam_file).split('.')[-2]
    
    try:
        if os.path.exists(output_file):
            logger.info(f"Output file already exists, skipping: {output_file}")
            return True
            
        if not os.path.exists(bam_file) or os.path.getsize(bam_file) == 0:
            logger.error(f"BAM file missing or empty: {bam_file}")
            return False
            
        if not os.path.exists(variant_file) or os.path.getsize(variant_file) == 0:
            logger.error(f"Variant file missing or empty: {variant_file}")
            return False
            
        if not os.path.exists(meta_file) or os.path.getsize(meta_file) == 0:
            logger.error(f"Meta file missing or empty: {meta_file}")
            return False
            
        os.makedirs(temp_dir, exist_ok=True)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Use custom script if provided, otherwise use default
        if custom_script and os.path.exists(custom_script):
            genotype_script = custom_script
        else:
            genotype_script = os.path.join(scomatic_scripts_dir, "SingleCellGenotype", "SingleCellGenotype.py")
        
        cmd = [
            "python",
            genotype_script,
            "--bam", bam_file,
            "--infile", variant_file,
            "--nprocs", "1",
            "--meta", meta_file,
            "--outfile", output_file,
            "--tmp_dir", temp_dir,
            "--ref", ref_genome
        ]
        logger.info(f"Computing genotypes for {cell_type}...")
        
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            logger.error(f"Output file not created or empty: {output_file}")
            return False
            
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Single cell genotype failed for {bam_file}: {e}")
        return False
    finally:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)


def filter_and_annotate_sc_genotypes(input_dir, variant_file, output_root, sample_id):
    """Filter single cell genotypes and add trinucleotide context"""
    final_output_dir = os.path.join(output_root, 'FilteredSingleCellAlleles')
    os.makedirs(final_output_dir, exist_ok=True)
    
    # Load variant information
    logger.info(f"Loading variant context information from {variant_file}...")
    variants_dict = {}
    context_issues = 0
    
    try:
        with open(variant_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or line.startswith('Chr'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    context_issues += 1
                    continue
                    
                chrom = parts[0]
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]
                up_context = parts[7]
                down_context = parts[8]
                
                # Extract flanking bases
                up_base = up_context[1] if len(up_context) >= 2 else 'N'
                down_base = down_context[1] if len(down_context) >= 2 else 'N'
                
                key = f"{chrom}_{pos}"
                variants_dict[key] = {
                    'ref': ref,
                    'alt': alt,
                    'up_base': up_base,
                    'down_base': down_base
                }
                
    except Exception as e:
        logger.error(f"Error loading variant file {variant_file}: {str(e)}")
        return [], None
    
    if context_issues > 0:
        logger.warning(f"{context_issues} variants had incomplete context information")

    # Process each cell type file
    logger.info(f"Processing cell type files in {input_dir}...")
    input_files = glob.glob(os.path.join(input_dir, '*.single_cell_genotype.tsv'))
    
    if not input_files:
        logger.warning("No input files found matching pattern *.single_cell_genotype.tsv")
        return [], None
    
    processed_files = []
    all_cells_data = []

    # Tracking statistics
    total_rows = 0
    filtered_by_base_mismatch = 0
    filtered_by_alt_reads = 0
    filtered_by_total_depth = 0
    rows_kept = 0
    
    for input_file in input_files:
        try:
            if os.path.getsize(input_file) == 0:
                logger.warning(f"Skipping empty file: {input_file}")
                continue
                
            df = pd.read_csv(input_file, sep='\t')
            if df.empty:
                logger.warning(f"Empty dataframe in {input_file}")
                continue
            
            initial_count = len(df)
            total_rows += initial_count
            
            def get_alt_bases(alt_str):
                """Extract unique bases from ALT_expected"""
                if pd.isna(alt_str) or alt_str == '.':
                    return set()
                alt_str = str(alt_str).replace('|', ',')
                bases = [base.strip() for base in alt_str.split(',') if base.strip()]
                return set(bases)
            
            df['ALT_bases_set'] = df['ALT_expected'].apply(get_alt_bases)
            
            def base_matches_alt(row):
                """Check if Base_observed matches any of the ALT_expected bases"""
                if not row['ALT_bases_set']:
                    return False
                return row['Base_observed'] in row['ALT_bases_set']
            
            mask_base = df.apply(base_matches_alt, axis=1)
            filtered_by_base_mismatch += (~mask_base).sum()
            
            # Filter 2: Must have >= 3 reads supporting ALT
            mask_alt_reads = df['Num_reads'] >= 3
            filtered_by_alt_reads += (~mask_alt_reads & mask_base).sum()
            
            # Filter 3: Must have >= 5 total depth
            mask_total_depth = df['Total_depth'] >= 5
            filtered_by_total_depth += (~mask_total_depth & mask_base & mask_alt_reads).sum()
            
            # Combined filter
            df_filtered = df[mask_base & mask_alt_reads & mask_total_depth].copy()
            rows_kept += len(df_filtered)
            
            df_filtered = df_filtered.drop(columns=['ALT_bases_set'])
            
            if df_filtered.empty:
                logger.info(f"No matching variants after filtering in {input_file}")
                continue
                
            # Add barcode suffix
            df_filtered['CB'] = df_filtered['CB'] + f"-1-{sample_id}"
            
            # Add trinucleotide context
            def get_trinucleotide(row):
                chrom_pos = f"{row['#CHROM']}_{row['Start']}"
                variant_info = variants_dict.get(chrom_pos, {})
                ref_tri = variant_info.get('up_base', 'N') + row['REF'] + variant_info.get('down_base', 'N')
                alt_tri = variant_info.get('up_base', 'N') + row['Base_observed'] + variant_info.get('down_base', 'N')
                return pd.Series([ref_tri, alt_tri])
            
            df_filtered[['REF_TRI', 'ALT_TRI']] = df_filtered.apply(get_trinucleotide, axis=1)
            
            # Save cell type file
            cell_type = os.path.basename(input_file).split('.')[0]
            output_file = os.path.join(final_output_dir, f"{cell_type}.single_cell_genotype.filtered.tsv")
            df_filtered.to_csv(output_file, sep='\t', index=False)
            processed_files.append(output_file)
            
            all_cells_data.append(df_filtered)
            
        except Exception as e:
            logger.error(f"ERROR processing {input_file}: {str(e)}")
            continue
    
    # Print filtering statistics
    logger.info("="*60)
    logger.info("FILTERING STATISTICS")
    logger.info("="*60)
    logger.info(f"Total rows processed: {total_rows}")
    logger.info(f"Filtered by base mismatch: {filtered_by_base_mismatch}")
    logger.info(f"Filtered by insufficient ALT reads (<3): {filtered_by_alt_reads}")
    logger.info(f"Filtered by insufficient total depth (<5X): {filtered_by_total_depth}")
    logger.info(f"Rows kept after all filters: {rows_kept}")
    if total_rows > 0:
        logger.info(f"Percentage retained: {100 * rows_kept / total_rows:.2f}%")
    logger.info("="*60)
    
    # Create combined output
    if all_cells_data:
        combined_cells = pd.concat(all_cells_data, ignore_index=True)
        combined_output = os.path.join(final_output_dir, "all_cell.single_cell_genotype.filtered.tsv")
        combined_cells.to_csv(combined_output, sep='\t', index=False)
        processed_files.append(combined_output)
        return processed_files, combined_cells
    return processed_files, None


def run_trinucleotide_context(args):
    """Compute trinucleotide context background"""
    step1_files_list, output_file, scomatic_scripts_dir = args
    try:
        temp_dir = os.path.join(os.path.dirname(output_file), "temp_trinuc")
        os.makedirs(temp_dir, exist_ok=True)
        
        input_list_file = os.path.join(temp_dir, "trinuc_input_list.txt")
        
        with open(input_list_file, 'w') as f:
            for file in step1_files_list:
                if os.path.exists(file):
                    f.write(f"{file}\n")
        
        if os.path.getsize(input_list_file) == 0:
            logger.warning("No valid input files for trinucleotide context")
            return False
        
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "TrinucleotideBackground", "TrinucleotideContextBackground.py"),
            "--in_tsv", input_list_file,
            "--out_file", output_file
        ]
        logger.info("Computing trinucleotide context background...")
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Trinucleotide context failed: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error in trinucleotide context: {str(e)}")
        return False
    finally:
        if 'temp_dir' in locals() and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)


def generate_complete_callable_sites(output_dir, valid_samples, adata_pp, cell_annotations, sample_metadata):
    """Generate complete callable sites file including all cells"""
    combined_callable_dir = os.path.join(output_dir, 'CombinedCallableSites')
    os.makedirs(combined_callable_dir, exist_ok=True)
    
    all_cells = pd.DataFrame({'CB': cell_annotations['cell_barcodes']})
    all_cells['SitesPerCell'] = 0
    
    stats = {
        'total_cells': len(all_cells),
        'processed_cells': 0,
        'missing_cells': 0,
        'samples_processed': 0,
        'samples_missing': 0,
        'files_processed': 0,
        'files_skipped': 0
    }
    
    for run_acc in tqdm(valid_samples, desc="Processing callable sites"):
        mask = adata_pp.obs['run_accession'] == run_acc
        subset = adata_pp[mask, :]
        series_id = subset.obs['series_id'].iloc[0] if 'series_id' in subset.obs.columns else run_acc
        
        callable_sites_dir = os.path.join(
            output_dir, 'scomatic', series_id, run_acc, 'UniqueCellCallableSites'
        )
        
        if not os.path.exists(callable_sites_dir):
            stats['samples_missing'] += 1
            continue
            
        stats['samples_processed'] += 1
        
        callable_files = glob.glob(os.path.join(callable_sites_dir, "*.SitesPerCell.tsv"))
        
        if not callable_files:
            callable_files = glob.glob(os.path.join(callable_sites_dir, "*.tsv"))
        
        for file_path in callable_files:
            try:
                if os.path.getsize(file_path) < 20:
                    stats['files_skipped'] += 1
                    continue
                
                with open(file_path, 'r') as f:
                    first_line = f.readline()
                
                if ',' in first_line:
                    df = pd.read_csv(file_path, sep=',')
                elif '\t' in first_line:
                    df = pd.read_csv(file_path, sep='\t')
                else:
                    stats['files_skipped'] += 1
                    continue
                
                if df.empty:
                    stats['files_skipped'] += 1
                    continue
                
                if 'CB' not in df.columns or 'SitesPerCell' not in df.columns:
                    stats['files_skipped'] += 1
                    continue
                    
                df['CB'] = df['CB'] + f"-1-{run_acc}"
                
                cells_updated = 0
                for _, row in df.iterrows():
                    if row['CB'] in all_cells['CB'].values:
                        idx = all_cells.index[all_cells['CB'] == row['CB']][0]
                        all_cells.at[idx, 'SitesPerCell'] = row['SitesPerCell']
                        cells_updated += 1
                
                stats['processed_cells'] += cells_updated
                stats['files_processed'] += 1
                    
            except Exception as e:
                logger.error(f"Error processing {file_path}: {str(e)}")
                stats['files_skipped'] += 1
                continue
            
    stats['missing_cells'] = stats['total_cells'] - stats['processed_cells']
    
    output_path = os.path.join(combined_callable_dir, 'complete_callable_sites.tsv')
    all_cells.to_csv(output_path, sep='\t', index=False)
    
    stats_path = os.path.join(combined_callable_dir, 'callable_sites_stats.txt')
    with open(stats_path, 'w') as f:
        f.write("Callable Sites Processing Statistics\n")
        f.write("="*50 + "\n")
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")
    
    logger.info(f"Saved complete callable sites to: {output_path}")
    logger.info(f"Statistics: {stats}")
    
    return output_path


def prepare_sample_args(run_acc, adata_pp, output_dir, scomatic_scripts_dir, ref_genome, custom_genotype_script=None):
    """Prepare arguments for sample processing with phase completion validation"""
    if not all([adata_pp is not None, output_dir, scomatic_scripts_dir, ref_genome]):
        raise ValueError("Missing required arguments")
    
    mask = adata_pp.obs['run_accession'] == run_acc
    if not any(mask):
        logger.info(f"SKIPPING {run_acc}: Not found in AnnData")
        return None
        
    subset = adata_pp[mask, :]
    series_id = subset.obs['series_id'].iloc[0] if 'series_id' in subset.obs.columns else run_acc
    
    scomatic_dir = os.path.join(output_dir, 'scomatic', series_id, run_acc)
    variant_calling_dir = os.path.join(scomatic_dir, 'VariantCalling')
    
    # Check required files
    meta_file = os.path.join(scomatic_dir, 'cell_barcode_annotation.tsv')
    if not os.path.exists(meta_file):
        logger.info(f"SKIPPING {run_acc}: Missing annotation file")
        return None
    
    bam_files = glob.glob(os.path.join(scomatic_dir, f"{run_acc}.*.bam"))
    if not bam_files:
        logger.info(f"SKIPPING {run_acc}: No split BAM files found")
        return None
    
    basecell_counts_dir = os.path.join(scomatic_dir, 'BaseCellCounts')
    if not os.path.exists(basecell_counts_dir):
        logger.info(f"SKIPPING {run_acc}: BaseCellCounts directory missing")
        return None
    
    merged_counts_file = os.path.join(scomatic_dir, 'MergedCounts', f"{run_acc}.BaseCellCounts.AllCellTypes.tsv")
    if not os.path.exists(merged_counts_file) or os.path.getsize(merged_counts_file) == 0:
        logger.info(f"SKIPPING {run_acc}: Merged counts file missing or empty")
        return None
    
    step2_file = os.path.join(variant_calling_dir, f"{run_acc}.calling.step2.tsv")
    if not os.path.exists(step2_file) or os.path.getsize(step2_file) == 0:
        logger.info(f"SKIPPING {run_acc}: Variant calling step2 file missing")
        return None
    
    filtered_dir = os.path.join(scomatic_dir, 'FilteredVariants')
    filtered_files = glob.glob(os.path.join(filtered_dir, "*.filtered.tsv"))
    if not filtered_files:
        logger.info(f"SKIPPING {run_acc}: No filtered variant files found")
        return None
    
    # Prepare arguments
    callable_sites_dir = os.path.join(scomatic_dir, 'CellTypeCallableSites')
    step1_file = os.path.join(variant_calling_dir, f"{run_acc}.calling.step1.tsv")
    
    callable_args = None
    if os.path.exists(step1_file) and os.path.getsize(step1_file) > 0:
        os.makedirs(callable_sites_dir, exist_ok=True)
        callable_args = (step1_file, os.path.join(callable_sites_dir, run_acc), "150", "2", scomatic_scripts_dir)
    
    sites_per_cell_args = []
    sites_per_cell_dir = os.path.join(scomatic_dir, 'UniqueCellCallableSites')
    
    if os.path.exists(step1_file) and os.path.getsize(step1_file) > 0:
        for bam in bam_files:
            if os.path.exists(bam) and os.path.getsize(bam) > 0:
                temp_dir = os.path.join(sites_per_cell_dir, f"temp_{os.path.basename(bam).split('.')[-2]}")
                sites_per_cell_args.append((bam, step1_file, ref_genome, sites_per_cell_dir, temp_dir, scomatic_scripts_dir))
    
    genotype_args = []
    single_cell_dir = os.path.join(scomatic_dir, 'SingleCellAlleles')
    
    if os.path.exists(step2_file) and os.path.getsize(step2_file) > 0:
        for bam in bam_files:
            cell_type = os.path.basename(bam).split('.')[-2]
            output_file = os.path.join(single_cell_dir, f"{cell_type}.single_cell_genotype.tsv")
            temp_dir = os.path.join(single_cell_dir, f"temp_{cell_type}")
            genotype_args.append((bam, step2_file, meta_file, ref_genome, output_file, temp_dir, scomatic_scripts_dir, custom_genotype_script))
    
    filtered_genotype_files = []
    if os.path.exists(step2_file) and os.path.getsize(step2_file) > 0:
        filtered_genotype_files = [(step2_file, single_cell_dir)]
    
    return {
        'run_acc': run_acc,
        'callable_args': callable_args,
        'sites_per_cell_args': sites_per_cell_args,
        'genotype_args': genotype_args,
        'filtered_genotype_files': filtered_genotype_files,
        'step1_file': step1_file,
        'directories': {
            'root': scomatic_dir,
            'callable_sites': callable_sites_dir,
            'sites_per_cell': sites_per_cell_dir,
            'single_cell': single_cell_dir,
            'variant_calling': variant_calling_dir
        }
    }


def process_sample(sample_args):
    """Process all steps for a single sample"""
    run_acc = sample_args['run_acc']
    results = {'run_acc': run_acc, 'success': True}
    
    try:
        for dir_path in sample_args['directories'].values():
            if dir_path:
                os.makedirs(dir_path, exist_ok=True)
        
        # 1. Callable sites per cell type
        if sample_args['callable_args']:
            try:
                success = run_callable_sites(sample_args['callable_args'])
                results['callable_sites'] = success
                if not success:
                    results['success'] = False
            except Exception as e:
                logger.error(f"Error in callable sites for {run_acc}: {str(e)}")
                results['callable_sites'] = False
                results['success'] = False
        
        # 2. Sites per cell
        if sample_args['sites_per_cell_args']:
            try:
                with ThreadPool(min(cpu_count(), len(sample_args['sites_per_cell_args']))) as pool:
                    sites_results = list(pool.imap(run_sites_per_cell, sample_args['sites_per_cell_args']))
                results['sites_per_cell'] = sum(sites_results)
            except Exception as e:
                logger.error(f"Error in sites per cell for {run_acc}: {str(e)}")
                results['sites_per_cell'] = 0
                results['success'] = False
        
        # 3. Single cell genotypes
        if sample_args['genotype_args']:
            try:
                with ThreadPool(min(cpu_count(), len(sample_args['genotype_args']))) as pool:
                    genotype_results = list(pool.imap(run_single_cell_genotype, sample_args['genotype_args']))
                
                results['single_cell_genotypes'] = sum(genotype_results)
                
                if sum(genotype_results) < len(genotype_results):
                    results['success'] = False
                    sample_args['filtered_genotype_files'] = []
                    results['filtered_genotypes'] = False
                    
            except Exception as e:
                logger.error(f"Error in genotype processing for {run_acc}: {str(e)}")
                results['single_cell_genotypes'] = 0
                results['success'] = False
                sample_args['filtered_genotype_files'] = []
                results['filtered_genotypes'] = False
        
        # 4. Filtered genotypes
        if sample_args['filtered_genotype_files'] and results.get('single_cell_genotypes', 0) > 0:
            try:
                variant_file, input_dir = sample_args['filtered_genotype_files'][0]
                
                if not os.path.exists(variant_file) or os.path.getsize(variant_file) == 0:
                    raise FileNotFoundError(f"Variant file not found: {variant_file}")
                
                if not os.path.exists(input_dir):
                    raise FileNotFoundError(f"Input directory not found: {input_dir}")
                
                processed_files, combined_data = filter_and_annotate_sc_genotypes(
                    input_dir, variant_file, sample_args['directories']['root'], run_acc
                )
                results['filtered_genotypes'] = bool(processed_files)
                results['combined_data'] = combined_data
                    
            except Exception as e:
                logger.error(f"Error in genotype filtering for {run_acc}: {str(e)}")
                results['filtered_genotypes'] = False
                results['success'] = False
        else:
            results['filtered_genotypes'] = False
                
        return results
        
    except Exception as e:
        logger.error(f"CRITICAL ERROR processing sample {run_acc}: {str(e)}")
        return {'run_acc': run_acc, 'success': False, 'error': str(e)}


def monitor_resources(interval=300, stop_event=None):
    """Monitor system resources during processing"""
    while stop_event is None or not stop_event.is_set():
        cpu = psutil.cpu_percent()
        mem = psutil.virtual_memory().percent
        disk = psutil.disk_usage('/').percent
        logger.info(f"Resource Monitor | CPU: {cpu}% | Mem: {mem}% | Disk: {disk}%")
        time.sleep(interval)
        if stop_event and stop_event.wait(timeout=0):
            break


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_scomatic_pipeline(
    adata_path,
    output_dir,
    sample_ids,
    scomatic_scripts_dir,
    ref_genome,
    editing_sites,
    pon_file,
    bed_file,
    cellranger_dir,
    custom_genotype_script=None,
    min_cov=5,
    min_cells=5,
    n_workers=None
):
    """
    Run the complete SComatic pipeline.
    
    Parameters
    ----------
    adata_path : str
        Path to the input AnnData file
    output_dir : str
        Output directory for results
    sample_ids : list
        List of sample IDs to process
    scomatic_scripts_dir : str
        Path to SComatic scripts directory
    ref_genome : str
        Path to reference genome FASTA
    editing_sites : str
        Path to RNA editing sites file
    pon_file : str
        Path to Panel of Normals file
    bed_file : str
        Path to BED file for mappable regions
    cellranger_dir : str
        Path to Cell Ranger output directory
    custom_genotype_script : str, optional
        Path to custom SingleCellGenotype.py script
    min_cov : int, optional
        Minimum coverage threshold (default: 5)
    min_cells : int, optional
        Minimum cells threshold (default: 5)
    n_workers : int, optional
        Number of workers for parallel processing
    
    Returns
    -------
    dict
        Dictionary with pipeline results and statistics
    """
    if n_workers is None:
        n_workers = cpu_count()
    
    # Create output directories
    mutations_dir = os.path.join(output_dir, 'mutations')
    os.makedirs(mutations_dir, exist_ok=True)
    
    # Start resource monitoring
    stop_event = threading.Event()
    monitor_thread = threading.Thread(target=monitor_resources, args=(300, stop_event))
    monitor_thread.daemon = True
    monitor_thread.start()
    
    try:
        # Load AnnData
        logger.info("Loading AnnData object...")
        adata_pp = sc.read_h5ad(adata_path)
        
        # Process barcodes
        adata_pp.obs['original_barcode'] = adata_pp.obs_names.str.split('-').str[0]
        adata_pp.obs['cell_barcodes'] = adata_pp.obs.index
        
        # Create cell annotations file
        adata_annotation_df = adata_pp.obs[['cell_barcodes', 'final_annotation']]
        cell_annotations_path = os.path.join(mutations_dir, 'cell_annotations.tsv')
        adata_annotation_df.to_csv(cell_annotations_path, sep="\t", header=True, index=False)
        
        # Get run accessions
        if 'run_accession' not in adata_pp.obs.columns:
            # Try to extract from barcode
            adata_pp.obs['run_accession'] = adata_pp.obs_names.str.split('-').str[-1]
        
        run_accessions = adata_pp.obs['run_accession'].unique()
        valid_samples = [s for s in run_accessions if s in sample_ids]
        
        logger.info(f"Found {len(valid_samples)} valid samples to process")
        
        # ======================
        # PHASE 1: Prepare and Split BAMs
        # ======================
        logger.info("="*80)
        logger.info("PHASE 1: Preparing inputs and running SplitBam")
        logger.info("="*80)
        
        scomatic_args = []
        for run_acc in tqdm(valid_samples, desc='Preparing samples'):
            mask = adata_pp.obs['run_accession'] == run_acc
            subset = adata_pp[mask, :]
            series_id = subset.obs['series_id'].iloc[0] if 'series_id' in subset.obs.columns else run_acc
            
            # Find BAM file
            sample_cellranger_dir = os.path.join(cellranger_dir, run_acc, "outs")
            original_bam = os.path.join(sample_cellranger_dir, "possorted_genome_bam.bam")
            
            if not os.path.exists(original_bam):
                logger.warning(f"Skipping {run_acc} - BAM not found")
                continue
            
            scomatic_sample_dir = os.path.join(mutations_dir, 'scomatic', series_id, run_acc)
            os.makedirs(scomatic_sample_dir, exist_ok=True)
            
            filtered_bam = os.path.join(scomatic_sample_dir, f"filtered_{run_acc}.bam")
            valid_barcodes = set(subset.obs['original_barcode'])
            
            if not os.path.exists(filtered_bam):
                if not create_filtered_bam(original_bam, filtered_bam, valid_barcodes):
                    continue
            
            annotation_file = os.path.join(scomatic_sample_dir, 'cell_barcode_annotation.tsv')
            if not os.path.exists(annotation_file):
                annotation_df = pd.DataFrame({
                    'Index': subset.obs['original_barcode'],
                    'Cell_type': subset.obs['final_annotation'].str.replace(r'[~.`!@#$%^&*(){|}/\\:;"\'<>?,=+\s]', '_', regex=True)
                })
                annotation_df.to_csv(annotation_file, sep='\t', index=False)
            
            scomatic_args.append((filtered_bam, annotation_file, run_acc, scomatic_sample_dir, scomatic_scripts_dir))
        
        # Run SplitBam
        if scomatic_args:
            with ThreadPool(min(n_workers, len(scomatic_args))) as pool:
                split_results = list(tqdm(
                    pool.imap(run_scomatic, scomatic_args),
                    total=len(scomatic_args),
                    desc="Running SplitBam"
                ))
            
            valid_samples = [args[2] for args, success in zip(scomatic_args, split_results) if success]
            logger.info(f"SplitBam completed: {len(valid_samples)}/{len(scomatic_args)} successful")
        
        # ======================
        # PHASE 2: BaseCellCounter
        # ======================
        logger.info("="*80)
        logger.info("PHASE 2: BaseCellCounter Processing")
        logger.info("="*80)
        
        basecell_counter_args = []
        for run_acc in valid_samples:
            mask = adata_pp.obs['run_accession'] == run_acc
            subset = adata_pp[mask, :]
            series_id = subset.obs['series_id'].iloc[0] if 'series_id' in subset.obs.columns else run_acc
            
            scomatic_sample_dir = os.path.join(mutations_dir, 'scomatic', series_id, run_acc)
            
            bam_files = glob.glob(os.path.join(scomatic_sample_dir, f"{run_acc}.*.bam"))
            if not bam_files:
                bam_files = [f for f in glob.glob(os.path.join(scomatic_sample_dir, "*.bam")) 
                            if not f.endswith(f"filtered_{run_acc}.bam")]
            
            if not bam_files:
                continue
                
            basecell_counts_dir = os.path.join(scomatic_sample_dir, 'BaseCellCounts')
            os.makedirs(basecell_counts_dir, exist_ok=True)
            
            for bam in bam_files:
                basecell_counter_args.append((bam, ref_genome, basecell_counts_dir, scomatic_scripts_dir))
        
        if basecell_counter_args:
            basecell_results = run_basecell_counter_parallel(basecell_counter_args, n_workers)
        
        # ======================
        # PHASE 3: Merge Counts
        # ======================
        logger.info("="*80)
        logger.info("PHASE 3: Merging Counts")
        logger.info("="*80)
        
        for run_acc in valid_samples:
            mask = adata_pp.obs['run_accession'] == run_acc
            subset = adata_pp[mask, :]
            series_id = subset.obs['series_id'].iloc[0] if 'series_id' in subset.obs.columns else run_acc
            
            scomatic_sample_dir = os.path.join(mutations_dir, 'scomatic', series_id, run_acc)
            basecell_counts_dir = os.path.join(scomatic_sample_dir, 'BaseCellCounts')
            
            if not os.path.exists(basecell_counts_dir):
                continue
                
            merged_counts_dir = os.path.join(scomatic_sample_dir, 'MergedCounts')
            os.makedirs(merged_counts_dir, exist_ok=True)
            
            output_file = os.path.join(merged_counts_dir, f"{run_acc}.BaseCellCounts.AllCellTypes.tsv")
            run_merge_counts(basecell_counts_dir, output_file, scomatic_scripts_dir, run_acc)
        
        # ======================
        # PHASE 4: Variant Calling
        # ======================
        logger.info("="*80)
        logger.info("PHASE 4: Variant Calling")
        logger.info("="*80)
        
        variant_step1_args = []
        variant_step2_args = []
        variant_valid_samples = []
        
        for run_acc in valid_samples:
            mask = adata_pp.obs['run_accession'] == run_acc
            subset = adata_pp[mask, :]
            series_id = subset.obs['series_id'].iloc[0] if 'series_id' in subset.obs.columns else run_acc
            
            scomatic_sample_dir = os.path.join(mutations_dir, 'scomatic', series_id, run_acc)
            merged_counts_file = os.path.join(scomatic_sample_dir, 'MergedCounts', f"{run_acc}.BaseCellCounts.AllCellTypes.tsv")
            
            if not os.path.exists(merged_counts_file) or os.path.getsize(merged_counts_file) == 0:
                continue
                
            variant_calling_dir = os.path.join(scomatic_sample_dir, 'VariantCalling')
            os.makedirs(variant_calling_dir, exist_ok=True)
            
            output_prefix_step1 = os.path.join(variant_calling_dir, run_acc)
            variant_step1_args.append((merged_counts_file, output_prefix_step1, ref_genome, scomatic_scripts_dir))
            variant_step2_args.append((f"{output_prefix_step1}.calling.step1.tsv", output_prefix_step1, editing_sites, pon_file, scomatic_scripts_dir))
            variant_valid_samples.append(run_acc)
        
        if variant_step1_args:
            with ThreadPool(min(n_workers, len(variant_step1_args))) as pool:
                step1_results = list(tqdm(
                    pool.imap(run_variant_calling_step1, variant_step1_args),
                    total=len(variant_step1_args),
                    desc="Variant Calling Step 1"
                ))
                
                step2_args = [args for args, success in zip(variant_step2_args, step1_results) if success]
                step2_results = list(tqdm(
                    pool.imap(run_variant_calling_step2, step2_args),
                    total=len(step2_args),
                    desc="Variant Calling Step 2"
                ))
        
        # ======================
        # PHASE 5: BED Filtering
        # ======================
        logger.info("="*80)
        logger.info("PHASE 5: BED Filtering")
        logger.info("="*80)
        
        final_valid_samples = []
        for run_acc in variant_valid_samples:
            mask = adata_pp.obs['run_accession'] == run_acc
            subset = adata_pp[mask, :]
            series_id = subset.obs['series_id'].iloc[0] if 'series_id' in subset.obs.columns else run_acc
            
            scomatic_sample_dir = os.path.join(mutations_dir, 'scomatic', series_id, run_acc)
            variant_calling_dir = os.path.join(scomatic_sample_dir, 'VariantCalling')
            filtered_dir = os.path.join(scomatic_sample_dir, 'FilteredVariants')
            os.makedirs(filtered_dir, exist_ok=True)
            
            step2_files = glob.glob(os.path.join(variant_calling_dir, "*.calling.step2.tsv"))
            for variant_file in step2_files:
                filename = os.path.basename(variant_file)
                filtered_file = os.path.join(filtered_dir, filename.replace('.step2.tsv', '.filtered.tsv'))
                if filter_with_bed(variant_file, bed_file, filtered_file):
                    final_valid_samples.append(run_acc)
        
        final_valid_samples = list(set(final_valid_samples))
        
        # ======================
        # PHASE 6: Final Processing
        # ======================
        logger.info("="*80)
        logger.info("PHASE 6: Final Processing")
        logger.info("="*80)
        
        all_sample_args = [
            prepare_sample_args(
                run_acc=run_acc,
                adata_pp=adata_pp,
                output_dir=mutations_dir,
                scomatic_scripts_dir=scomatic_scripts_dir,
                ref_genome=ref_genome,
                custom_genotype_script=custom_genotype_script
            ) for run_acc in final_valid_samples
        ]
        all_sample_args = [args for args in all_sample_args if args is not None]
        
        if all_sample_args:
            with ThreadPool(min(n_workers, len(all_sample_args))) as pool:
                final_results = list(tqdm(
                    pool.imap(process_sample, all_sample_args),
                    total=len(all_sample_args),
                    desc="Phase 6 Processing"
                ))
        else:
            final_results = []
        
        # Combine mutation data
        logger.info("Combining mutation data from all samples...")
        all_samples_data = []
        for result in final_results:
            if result.get('combined_data') is not None and not result['combined_data'].empty:
                all_samples_data.append(result['combined_data'])
        
        if all_samples_data:
            final_combined = pd.concat(all_samples_data, ignore_index=True)
            final_output_path = os.path.join(mutations_dir, 'all_samples.single_cell_genotype.filtered.tsv')
            final_combined.to_csv(final_output_path, sep='\t', index=False)
            logger.info(f"Combined mutations saved to: {final_output_path}")
        else:
            final_output_path = os.path.join(mutations_dir, 'all_samples.single_cell_genotype.filtered.tsv')
            # Create empty file with header for Snakemake to find
            pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'Cell', 'Sample']).to_csv(
                final_output_path, sep='\t', index=False
            )
            logger.warning("No mutation data collected - created empty output file")
        
        # Compute trinucleotide context
        trinuc_step1_files = []
        for sample_args in all_sample_args:
            step1_file = sample_args['step1_file']
            if os.path.exists(step1_file) and os.path.getsize(step1_file) > 0:
                trinuc_step1_files.append(step1_file)
        
        output_trinuc = os.path.join(mutations_dir, 'trinucleotide_background.tsv')
        if trinuc_step1_files:
            trinuc_args = (trinuc_step1_files, output_trinuc, scomatic_scripts_dir)
            run_trinucleotide_context(trinuc_args)
        else:
            # Create empty trinucleotide file for Snakemake
            pd.DataFrame(columns=['Context', 'Count']).to_csv(output_trinuc, sep='\t', index=False)
            logger.warning("No trinucleotide data - created empty output file")
        
        # Generate callable sites
        cell_annotations = pd.read_csv(cell_annotations_path, sep='\t')
        callable_path = generate_complete_callable_sites(
            mutations_dir, final_valid_samples, adata_pp, cell_annotations, None
        )
        
        # Summary
        logger.info("="*80)
        logger.info("SCOMATIC PIPELINE COMPLETE")
        logger.info("="*80)
        logger.info(f"Samples processed: {len(final_valid_samples)}")
        logger.info(f"Mutations file: {final_output_path}")
        logger.info(f"Callable sites: {callable_path}")
        logger.info(f"Cell annotations: {cell_annotations_path}")
        logger.info(f"Trinucleotide background: {output_trinuc}")
        
        return {
            'mutations_file': final_output_path,
            'callable_sites': callable_path,
            'cell_annotations': cell_annotations_path,
            'trinucleotide_background': output_trinuc,
            'samples_processed': len(final_valid_samples),
            'final_results': final_results
        }
        
    except Exception as e:
        logger.error(f"CRITICAL ERROR in pipeline: {str(e)}")
        import traceback
        traceback.print_exc()
        raise
    finally:
        stop_event.set()
        monitor_thread.join(timeout=5)


# =============================================================================
# SNAKEMAKE/CLI INTERFACE
# =============================================================================

def run_from_snakemake():
    """Run pipeline from Snakemake context"""
    run_scomatic_pipeline(
        adata_path=snakemake.input.adata,
        output_dir=snakemake.params.output_dir,
        sample_ids=snakemake.params.sample_ids,
        scomatic_scripts_dir=snakemake.params.scomatic_scripts_dir,
        ref_genome=snakemake.params.ref_genome,
        editing_sites=snakemake.params.editing_sites,
        pon_file=snakemake.params.pon_file,
        bed_file=snakemake.params.bed_file,
        cellranger_dir=snakemake.params.cellranger_dir,
        custom_genotype_script=getattr(snakemake.params, 'custom_genotype_script', None),
        min_cov=getattr(snakemake.params, 'min_cov', 5),
        min_cells=getattr(snakemake.params, 'min_cells', 5),
        n_workers=snakemake.threads
    )


def main():
    """CLI entry point"""
    parser = argparse.ArgumentParser(
        description="SComatic mutation calling pipeline for ClusterCatcher"
    )
    parser.add_argument('--adata', required=True, help='Path to input AnnData file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--sample-ids', required=True, nargs='+', help='Sample IDs to process')
    parser.add_argument('--scomatic-scripts-dir', required=True, help='Path to SComatic scripts')
    parser.add_argument('--ref-genome', required=True, help='Path to reference genome FASTA')
    parser.add_argument('--editing-sites', required=True, help='Path to RNA editing sites file')
    parser.add_argument('--pon-file', required=True, help='Path to Panel of Normals file')
    parser.add_argument('--bed-file', required=True, help='Path to BED file for mappable regions')
    parser.add_argument('--cellranger-dir', required=True, help='Path to Cell Ranger output directory')
    parser.add_argument('--custom-genotype-script', help='Path to custom SingleCellGenotype.py')
    parser.add_argument('--min-cov', type=int, default=5, help='Minimum coverage threshold')
    parser.add_argument('--min-cells', type=int, default=5, help='Minimum cells threshold')
    parser.add_argument('--workers', type=int, default=None, help='Number of workers')
    
    args = parser.parse_args()
    
    run_scomatic_pipeline(
        adata_path=args.adata,
        output_dir=args.output_dir,
        sample_ids=args.sample_ids,
        scomatic_scripts_dir=args.scomatic_scripts_dir,
        ref_genome=args.ref_genome,
        editing_sites=args.editing_sites,
        pon_file=args.pon_file,
        bed_file=args.bed_file,
        cellranger_dir=args.cellranger_dir,
        custom_genotype_script=args.custom_genotype_script,
        min_cov=args.min_cov,
        min_cells=args.min_cells,
        n_workers=args.workers
    )


if __name__ == '__main__':
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        main()

