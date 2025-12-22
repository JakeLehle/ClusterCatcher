#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_viral_detection.py
============================

Summarize viral/microbial detection results across all samples.

This script:
1. Aggregates organism counts across all samples
2. Identifies samples with specific organisms of interest
3. Generates summary tables and visualizations
4. Creates a combined AnnData object for integration with cell annotations

Usage:
    Called via Snakemake rule with snakemake.input/output/params
"""

import os
import sys
import gzip
import logging
import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import vstack, csr_matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_kraken_matrix(matrix_dir):
    """
    Load Kraken2 sparse matrix output.
    
    Parameters
    ----------
    matrix_dir : str
        Path to kraken2_filtered_feature_bc_matrix directory
        
    Returns
    -------
    tuple
        (sparse_matrix, barcodes, features)
    """
    # Load matrix
    matrix_file = os.path.join(matrix_dir, 'matrix.mtx.gz')
    if os.path.exists(matrix_file):
        with gzip.open(matrix_file, 'rb') as f:
            matrix = mmread(f).tocsr()
    else:
        matrix_file = os.path.join(matrix_dir, 'matrix.mtx')
        if os.path.exists(matrix_file):
            matrix = mmread(matrix_file).tocsr()
        else:
            return None, [], []
    
    # Load barcodes
    barcodes = []
    barcodes_file = os.path.join(matrix_dir, 'barcodes.tsv.gz')
    if os.path.exists(barcodes_file):
        with gzip.open(barcodes_file, 'rt') as f:
            barcodes = [line.strip() for line in f]
    else:
        barcodes_file = os.path.join(matrix_dir, 'barcodes.tsv')
        if os.path.exists(barcodes_file):
            with open(barcodes_file, 'r') as f:
                barcodes = [line.strip() for line in f]
    
    # Load features
    features = []
    features_file = os.path.join(matrix_dir, 'features.tsv.gz')
    if os.path.exists(features_file):
        with gzip.open(features_file, 'rt') as f:
            features = [line.strip().split('\t') for line in f]
    else:
        # Try alternative names
        for alt_name in ['features.tsv', 'genes.tsv.gz', 'genes.tsv']:
            alt_file = os.path.join(matrix_dir, alt_name)
            if os.path.exists(alt_file):
                opener = gzip.open if alt_name.endswith('.gz') else open
                mode = 'rt' if alt_name.endswith('.gz') else 'r'
                with opener(alt_file, mode) as f:
                    features = [line.strip().split('\t') for line in f]
                break
    
    return matrix, barcodes, features


def aggregate_sample_summaries(sample_dirs, sample_ids):
    """
    Aggregate organism summaries across all samples.
    
    Parameters
    ----------
    sample_dirs : list
        List of sample output directories
    sample_ids : list
        List of sample identifiers
        
    Returns
    -------
    pd.DataFrame
        Aggregated summary with organism counts per sample
    """
    logger.info("Aggregating sample summaries...")
    
    all_data = []
    
    for sample_dir, sample_id in zip(sample_dirs, sample_ids):
        summary_file = os.path.join(sample_dir, f"{sample_id}_organism_summary.tsv")
        
        if os.path.exists(summary_file):
            try:
                df = pd.read_csv(summary_file, sep='\t')
                if len(df) > 0:
                    df['sample_id'] = sample_id
                    all_data.append(df)
                else:
                    logger.warning(f"Empty summary file for {sample_id}")
            except Exception as e:
                logger.warning(f"Could not read summary for {sample_id}: {e}")
        else:
            logger.warning(f"Summary file not found for {sample_id}")
    
    if not all_data:
        logger.warning("No summary data found!")
        return pd.DataFrame()
    
    combined = pd.concat(all_data, ignore_index=True)
    
    return combined


def identify_organisms_of_interest(combined_df, organisms_of_interest):
    """
    Filter for specific organisms of interest.
    
    Parameters
    ----------
    combined_df : pd.DataFrame
        Combined summary dataframe
    organisms_of_interest : list
        List of organism name patterns to search for
        
    Returns
    -------
    pd.DataFrame
        Filtered dataframe with organisms of interest
    """
    if not organisms_of_interest or combined_df.empty:
        return combined_df
    
    logger.info(f"Filtering for organisms of interest: {organisms_of_interest}")
    
    mask = combined_df['name'].str.lower().apply(
        lambda x: any(org.lower() in str(x) for org in organisms_of_interest)
    )
    
    filtered = combined_df[mask].copy()
    logger.info(f"  Found {len(filtered)} entries matching organisms of interest")
    
    return filtered


def create_organism_sample_matrix(combined_df):
    """
    Create organism x sample count matrix.
    
    Parameters
    ----------
    combined_df : pd.DataFrame
        Combined summary dataframe
        
    Returns
    -------
    pd.DataFrame
        Pivot table with organisms as rows and samples as columns
    """
    if combined_df.empty:
        return pd.DataFrame()
    
    # Pivot to get organism x sample matrix
    pivot = combined_df.pivot_table(
        index='name',
        columns='sample_id',
        values='reads_direct',
        aggfunc='sum',
        fill_value=0
    )
    
    # Sort by total counts
    pivot['total'] = pivot.sum(axis=1)
    pivot = pivot.sort_values('total', ascending=False)
    pivot = pivot.drop('total', axis=1)
    
    return pivot


def generate_summary_plots(organism_matrix, output_dir, top_n=20):
    """
    Generate summary visualization plots.
    
    Parameters
    ----------
    organism_matrix : pd.DataFrame
        Organism x sample matrix
    output_dir : str
        Output directory for plots
    top_n : int
        Number of top organisms to show
    """
    if organism_matrix.empty:
        logger.warning("No data for plotting")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot 1: Heatmap of top organisms
    top_organisms = organism_matrix.head(top_n)
    
    if not top_organisms.empty and top_organisms.shape[1] > 0:
        try:
            fig, ax = plt.subplots(figsize=(12, max(8, len(top_organisms) * 0.4)))
            sns.heatmap(
                np.log10(top_organisms + 1),
                cmap='YlOrRd',
                annot=False,
                ax=ax
            )
            ax.set_title(f'Top {top_n} Detected Organisms (log10 reads)')
            ax.set_xlabel('Sample')
            ax.set_ylabel('Organism')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'organism_heatmap.pdf'), bbox_inches='tight')
            plt.savefig(os.path.join(output_dir, 'organism_heatmap.png'), dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            logger.warning(f"Could not generate heatmap: {e}")
    
    # Plot 2: Bar chart of total reads per organism
    total_counts = organism_matrix.sum(axis=1).head(top_n)
    
    if not total_counts.empty:
        try:
            fig, ax = plt.subplots(figsize=(10, max(6, len(total_counts) * 0.3)))
            total_counts.plot(kind='barh', ax=ax)
            ax.set_xlabel('Total Reads')
            ax.set_title(f'Top {top_n} Organisms by Total Reads')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'organism_barplot.pdf'), bbox_inches='tight')
            plt.savefig(os.path.join(output_dir, 'organism_barplot.png'), dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            logger.warning(f"Could not generate barplot: {e}")
    
    # Plot 3: Samples with detections
    samples_with_detections = (organism_matrix > 0).sum(axis=0)
    
    if not samples_with_detections.empty:
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            samples_with_detections.plot(kind='bar', ax=ax)
            ax.set_xlabel('Sample')
            ax.set_ylabel('Number of Organisms Detected')
            ax.set_title('Organism Diversity per Sample')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'sample_diversity.pdf'), bbox_inches='tight')
            plt.savefig(os.path.join(output_dir, 'sample_diversity.png'), dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            logger.warning(f"Could not generate diversity plot: {e}")
    
    logger.info(f"Plots saved to {output_dir}")


def create_combined_anndata(sample_dirs, sample_ids, output_file):
    """
    Create combined AnnData object from all sample matrices.
    
    Parameters
    ----------
    sample_dirs : list
        List of sample output directories
    sample_ids : list
        List of sample identifiers
    output_file : str
        Path for output H5AD file
        
    Returns
    -------
    AnnData
        Combined AnnData object
    """
    logger.info("Creating combined AnnData object...")
    
    adatas = []
    
    for sample_dir, sample_id in zip(sample_dirs, sample_ids):
        matrix_dir = os.path.join(sample_dir, 'kraken2_filtered_feature_bc_matrix')
        
        if not os.path.exists(matrix_dir):
            logger.warning(f"Matrix directory not found for {sample_id}")
            continue
        
        try:
            matrix, barcodes, features = load_kraken_matrix(matrix_dir)
            
            if matrix is None or len(barcodes) == 0:
                logger.warning(f"Empty matrix for {sample_id}")
                continue
            
            # Create AnnData
            adata = sc.AnnData(matrix.T)  # Transpose: cells x organisms
            adata.obs_names = [f"{sample_id}_{bc}" for bc in barcodes]
            adata.var_names = [f[1] if len(f) > 1 else f[0] for f in features]
            adata.var['taxonomy_id'] = [f[0] for f in features]
            adata.obs['sample_id'] = sample_id
            
            adatas.append(adata)
            logger.info(f"  {sample_id}: {adata.n_obs} cells, {adata.n_vars} organisms")
            
        except Exception as e:
            logger.warning(f"Could not load matrix for {sample_id}: {e}")
    
    if not adatas:
        logger.warning("No AnnData objects created!")
        # Create empty AnnData
        combined = sc.AnnData()
        combined.write(output_file)
        return combined
    
    # Concatenate
    if len(adatas) == 1:
        combined = adatas[0]
    else:
        combined = sc.concat(adatas, join='outer')
    
    combined.obs['sample_id'] = combined.obs['sample_id'].astype('category')
    
    # Save
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    combined.write(output_file)
    logger.info(f"Combined AnnData saved to {output_file}")
    logger.info(f"  Total: {combined.n_obs} cells, {combined.n_vars} organisms")
    
    return combined


def generate_final_summary(combined_df, organism_matrix, output_file):
    """
    Generate final summary TSV file.
    
    Parameters
    ----------
    combined_df : pd.DataFrame
        Combined summary dataframe
    organism_matrix : pd.DataFrame
        Organism x sample matrix
    output_file : str
        Path for output TSV
    """
    logger.info(f"Generating final summary: {output_file}")
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    if organism_matrix.empty:
        # Create empty summary
        pd.DataFrame(columns=['organism', 'total_reads', 'n_samples']).to_csv(
            output_file, sep='\t', index=False
        )
        return
    
    # Summary statistics
    summary = pd.DataFrame({
        'organism': organism_matrix.index,
        'total_reads': organism_matrix.sum(axis=1).values,
        'n_samples': (organism_matrix > 0).sum(axis=1).values,
        'max_reads_single_sample': organism_matrix.max(axis=1).values,
        'mean_reads_when_present': organism_matrix.replace(0, np.nan).mean(axis=1).values
    })
    
    summary = summary.sort_values('total_reads', ascending=False)
    summary.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"  {len(summary)} organisms summarized")


def run_summarize_pipeline(
    sample_dirs,
    sample_ids,
    summary_file,
    adata_file,
    organisms_of_interest=None
):
    """
    Run the complete summarization pipeline.
    
    Parameters
    ----------
    sample_dirs : list
        List of sample viral detection directories
    sample_ids : list
        List of sample identifiers
    summary_file : str
        Path for output summary TSV
    adata_file : str
        Path for output AnnData H5AD
    organisms_of_interest : list, optional
        Organisms to highlight
    """
    logger.info("="*60)
    logger.info("Viral Detection Summary Pipeline")
    logger.info("="*60)
    
    output_dir = os.path.dirname(summary_file)
    
    # Aggregate summaries
    combined_df = aggregate_sample_summaries(sample_dirs, sample_ids)
    
    # Filter for organisms of interest
    if organisms_of_interest:
        filtered_df = identify_organisms_of_interest(combined_df, organisms_of_interest)
    else:
        filtered_df = combined_df
    
    # Create organism x sample matrix
    organism_matrix = create_organism_sample_matrix(filtered_df)
    
    # Generate plots
    generate_summary_plots(organism_matrix, os.path.join(output_dir, 'plots'))
    
    # Generate final summary
    generate_final_summary(combined_df, organism_matrix, summary_file)
    
    # Create combined AnnData
    create_combined_anndata(sample_dirs, sample_ids, adata_file)
    
    logger.info("\n" + "="*60)
    logger.info("Summary pipeline completed!")
    logger.info("="*60)


# =============================================================================
# Snakemake Integration
# =============================================================================

def run_from_snakemake():
    """Run from Snakemake rule."""
    
    # Get inputs
    sample_reports = snakemake.input.reports
    
    # Get outputs
    summary_file = snakemake.output.summary
    adata_file = snakemake.output.adata
    
    # Get params
    sample_ids = snakemake.params.sample_ids
    organisms_of_interest = snakemake.params.get('organisms_of_interest', [])
    viral_base_dir = snakemake.params.viral_base_dir
    
    # Build sample directories from base directory
    sample_dirs = [os.path.join(viral_base_dir, sid) for sid in sample_ids]
    
    # Set up logging
    if snakemake.log:
        os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
        file_handler = logging.FileHandler(snakemake.log[0])
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Run pipeline
    run_summarize_pipeline(
        sample_dirs=sample_dirs,
        sample_ids=sample_ids,
        summary_file=summary_file,
        adata_file=adata_file,
        organisms_of_interest=organisms_of_interest
    )


# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Summarize viral detection results across samples'
    )
    parser.add_argument('--input-dirs', nargs='+', required=True,
                        help='Sample output directories')
    parser.add_argument('--sample-ids', nargs='+', required=True,
                        help='Sample identifiers')
    parser.add_argument('--output', required=True, help='Output summary file')
    parser.add_argument('--adata-output', help='Output AnnData file (default: same dir as summary)')
    parser.add_argument('--organisms', nargs='*', 
                        help='Organisms of interest to highlight')
    
    args = parser.parse_args()
    
    # Determine AnnData output path
    if args.adata_output:
        adata_file = args.adata_output
    else:
        adata_file = os.path.join(os.path.dirname(args.output), 'viral_counts.h5ad')
    
    run_summarize_pipeline(
        sample_dirs=args.input_dirs,
        sample_ids=args.sample_ids,
        summary_file=args.output,
        adata_file=adata_file,
        organisms_of_interest=args.organisms
    )


if __name__ == '__main__':
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        main()
