#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
viral_integration.py
====================

Post-processing script for viral detection that integrates Kraken2 results
with annotated single-cell data.

This script:
1. Loads Kraken2 viral detection results for each sample
2. Filters for human-specific viruses using a reference database
3. Joins viral counts with gene expression data
4. Calculates differential expression of viruses across cell types
5. Creates comprehensive visualizations
6. Saves integrated AnnData with viral information

Requirements:
    - Kraken2 results from viral detection step
    - Annotated AnnData (adata_pp) from QC/annotation step
    - Human virus hierarchy file (from Kraken2 human viral database)

Usage:
    Called via Snakemake rule with snakemake.input/output/params
    
    Or standalone:
    python viral_integration.py --adata-pp adata_pp.h5ad --viral-dir viral/ \
                                --human-viral-db /path/to/human_viral/inspect.txt \
                                --output output_dir
"""

import os
import sys
import re
import gzip
import csv
import logging
import argparse
import warnings
from pathlib import Path
from tqdm import tqdm

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings
warnings.filterwarnings('ignore')
matplotlib.use('Agg')

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Utility Functions
# =============================================================================

def clean_column_names(df):
    """
    Clean column names to be HDF5-compatible.
    
    Replaces problematic characters that can cause issues when saving H5AD files.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with potentially problematic column names
        
    Returns
    -------
    pd.DataFrame
        DataFrame with cleaned column names
    """
    logger.info("Cleaning column names for HDF5 compatibility...")
    rename_map = {}
    
    for col in df.columns:
        new_col = str(col)
        # Replace problematic characters
        new_col = new_col.replace('/', '_')
        new_col = new_col.replace('\\', '_')
        new_col = new_col.replace(' ', '_')
        new_col = new_col.replace('(', '')
        new_col = new_col.replace(')', '')
        new_col = new_col.replace('[', '')
        new_col = new_col.replace(']', '')
        new_col = new_col.replace('{', '')
        new_col = new_col.replace('}', '')
        new_col = new_col.replace("'", '')
        new_col = new_col.replace('"', '')
        
        if new_col != col:
            rename_map[col] = new_col
            logger.debug(f"  Renamed: '{col}' -> '{new_col}'")
    
    if rename_map:
        df = df.rename(columns=rename_map)
        logger.info(f"  Cleaned {len(rename_map)} column names")
    else:
        logger.info("  No problematic column names found")
    
    return df


def clean_obs_for_saving(adata):
    """
    Clean AnnData obs columns for HDF5 saving.
    
    Converts all obs columns to strings to avoid serialization issues.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object to clean
        
    Returns
    -------
    AnnData
        Cleaned AnnData
    """
    for col in adata.obs.columns:
        adata.obs[col] = [str(element) for element in adata.obs[col]]
    
    if 'gene_ids' in adata.var.columns:
        adata.var['gene_ids'] = [str(item) for item in adata.var['gene_ids']]
    
    return adata


def parse_kraken_hierarchy(hierarchy_file):
    """
    Parse Kraken2 hierarchy/inspect file.
    
    Parameters
    ----------
    hierarchy_file : str
        Path to Kraken2 inspect.txt or hierarchy.txt file
        
    Returns
    -------
    dict
        Dictionary mapping tax_id to virus information
    list
        List of species-level entries
    """
    logger.info(f"Parsing Kraken2 hierarchy file: {hierarchy_file}")
    
    v_hierarchy = {}
    species_list = []
    
    with open(hierarchy_file, 'rt') as f:
        lines = f.readlines()
    
    for line in lines:
        if not line.strip():
            continue
        
        parts = line.strip().split('\t')
        
        if len(parts) >= 6:
            # Standard inspect.txt format with tabs
            percentage = parts[0].strip()
            read_count = parts[1].strip()
            tax_count = parts[2].strip()
            rank_code = parts[3].strip()
            tax_id = parts[4].strip()
            name = parts[5].strip()
        else:
            # Try space-separated format
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            percentage = parts[0]
            read_count = parts[1]
            tax_count = parts[2]
            rank_code = parts[3]
            tax_id = parts[4]
            name = ' '.join(parts[5:])
        
        v_hierarchy[tax_id] = {
            'rank_code': rank_code,
            'name': name,
            'read_count': read_count,
            'tax_count': tax_count,
            'percentage': percentage
        }
        
        # Add to species list if it's a species (starts with S)
        if rank_code.startswith('S'):
            species_list.append({
                'tax_id': tax_id,
                'name': name,
                'read_count': read_count
            })
    
    logger.info(f"  Found {len(v_hierarchy)} total entries, {len(species_list)} species")
    
    return v_hierarchy, species_list


def parse_sample_hierarchy(hierarchy_file):
    """
    Parse sample-level hierarchy file (simpler format).
    
    Parameters
    ----------
    hierarchy_file : str
        Path to hierarchy.txt from sample Kraken2 output
        
    Returns
    -------
    dict
        Dictionary mapping tax_id to virus name
    """
    v_hierarchy = {}
    
    with open(hierarchy_file, 'rt') as f:
        lines = f.readlines()
    
    for line in lines:
        if not line.strip() or line.startswith('#'):
            continue
        
        parts = line.strip().split('\t')
        if len(parts) >= 6:
            tax_id = parts[4].strip()
            # Extract just the name, removing leading spaces
            name_match = re.search(r'.*?([a-zA-Z].*)', parts[5])
            if name_match:
                name = name_match.group(1).strip()
            else:
                name = parts[5].strip()
            v_hierarchy[tax_id] = name
    
    return v_hierarchy


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_viral_data_for_samples(sample_dirs, sample_ids, sample_hierarchy):
    """
    Load Kraken2 viral detection matrices for all samples.
    
    Parameters
    ----------
    sample_dirs : list
        List of sample viral detection output directories
    sample_ids : list
        List of sample identifiers
    sample_hierarchy : dict
        Virus name mapping from reference hierarchy
        
    Returns
    -------
    AnnData
        Combined viral counts AnnData
    """
    logger.info("Loading viral detection data for all samples...")
    
    adata_v = None
    
    for sample_dir, sample_id in tqdm(zip(sample_dirs, sample_ids), 
                                       total=len(sample_ids), 
                                       desc="Loading samples"):
        matrix_dir = os.path.join(sample_dir, 'kraken2_filtered_feature_bc_matrix')
        
        if not os.path.exists(matrix_dir):
            logger.warning(f"  Skipping {sample_id}: matrix directory not found")
            continue
        
        # Check required files
        required_files = ['barcodes.tsv.gz', 'matrix.mtx.gz']
        genes_file = 'features.tsv.gz' if os.path.exists(os.path.join(matrix_dir, 'features.tsv.gz')) else 'genes.tsv.gz'
        required_files.append(genes_file)
        
        missing = [f for f in required_files if not os.path.exists(os.path.join(matrix_dir, f))]
        if missing:
            logger.warning(f"  Skipping {sample_id}: missing files {missing}")
            continue
        
        try:
            # Decompress files temporarily
            temp_files = []
            for gz_file in required_files:
                gz_path = os.path.join(matrix_dir, gz_file)
                out_path = os.path.join(matrix_dir, gz_file[:-3])
                
                with gzip.open(gz_path, 'rt') as f_in, open(out_path, 'wt') as f_out:
                    f_out.writelines(f_in)
                temp_files.append(out_path)
            
            # Read the 10X matrix
            adata_sample = sc.read_10x_mtx(matrix_dir)
            
            # Add sample metadata
            adata_sample.obs['sample_id'] = sample_id
            
            # Create unique cell barcodes
            adata_sample.obs_names = [f"{bc}-{sample_id}" for bc in adata_sample.obs_names]
            
            # Add gene symbol column
            adata_sample.var['gene_symbol'] = adata_sample.var_names
            
            # Concatenate
            if adata_v is None:
                adata_v = adata_sample
            else:
                adata_v = ad.concat([adata_v, adata_sample], join='outer', merge='same')
            
            # Cleanup temp files
            for temp_file in temp_files:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                    
        except Exception as e:
            logger.error(f"  Error processing {sample_id}: {e}")
            continue
    
    if adata_v is None:
        logger.warning("No viral data was loaded!")
        return ad.AnnData()
    
    logger.info(f"  Loaded viral data: {adata_v.n_obs} cells, {adata_v.n_vars} organisms")
    
    return adata_v


def ensure_all_viruses_present(sample_dirs, sample_ids, reference_hierarchy):
    """
    Ensure all samples have entries for all viruses in the reference.
    
    This is necessary for proper concatenation of sparse matrices.
    
    Parameters
    ----------
    sample_dirs : list
        List of sample viral detection output directories
    sample_ids : list
        List of sample identifiers
    reference_hierarchy : dict
        Reference virus name mapping
    """
    logger.info("Ensuring consistent virus features across samples...")
    
    for sample_dir, sample_id in zip(sample_dirs, sample_ids):
        matrix_dir = os.path.join(sample_dir, 'kraken2_filtered_feature_bc_matrix')
        
        if not os.path.exists(matrix_dir):
            continue
        
        genes_file = os.path.join(matrix_dir, 'genes.tsv.gz')
        if not os.path.exists(genes_file):
            genes_file = os.path.join(matrix_dir, 'features.tsv.gz')
            if not os.path.exists(genes_file):
                continue
        
        try:
            # Read current genes
            existing_viruses = []
            with gzip.open(genes_file, 'rt') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if len(row) >= 2:
                        existing_viruses.append(row[1])
            
            # Add missing viruses
            viruses_to_add = []
            for tax_id, virus_name in reference_hierarchy.items():
                if virus_name not in existing_viruses:
                    viruses_to_add.append((tax_id, virus_name))
            
            if viruses_to_add:
                with gzip.open(genes_file, 'at') as f:
                    writer = csv.writer(f, delimiter='\t', lineterminator='\n')
                    for tax_id, virus_name in viruses_to_add:
                        writer.writerow([tax_id, virus_name])
                
                # Update matrix header with new feature count
                matrix_file = os.path.join(matrix_dir, 'matrix.mtx.gz')
                new_feature_count = len(existing_viruses) + len(viruses_to_add)
                
                with gzip.open(matrix_file, 'rb') as f:
                    lines = f.readlines()
                
                # Update the header line (line 3, 0-indexed as 2)
                if len(lines) >= 3:
                    header_parts = lines[2].decode().strip().split()
                    if len(header_parts) >= 3:
                        new_header = f"{new_feature_count} {header_parts[1]} {len(lines) - 3}\n"
                        lines[2] = new_header.encode()
                
                with gzip.open(matrix_file, 'wb') as f:
                    for line in lines:
                        f.write(line)
                
                logger.debug(f"  Added {len(viruses_to_add)} viruses to {sample_id}")
                
        except Exception as e:
            logger.warning(f"  Error updating {sample_id}: {e}")


# =============================================================================
# Integration Functions
# =============================================================================

def filter_human_viruses(adata_v, human_species_names):
    """
    Filter viral AnnData to only human-associated viruses.
    
    Parameters
    ----------
    adata_v : AnnData
        Viral counts AnnData
    human_species_names : tuple
        Tuple of human virus species names
        
    Returns
    -------
    AnnData
        Filtered AnnData with only human viruses
    """
    logger.info("Filtering for human-associated viruses...")
    
    human_virus_names = tuple(name.strip() for name in human_species_names if name.strip())
    human_virus_mask = adata_v.var_names.isin(human_virus_names)
    
    n_before = adata_v.n_vars
    adata_v = adata_v[:, human_virus_mask].copy()
    n_after = adata_v.n_vars
    
    logger.info(f"  Filtered from {n_before} to {n_after} organisms")
    
    return adata_v


def integrate_viral_with_expression(adata_pp, adata_v):
    """
    Integrate viral counts with gene expression data.
    
    Parameters
    ----------
    adata_pp : AnnData
        Preprocessed gene expression data (log-normalized)
    adata_v : AnnData
        Viral counts data (raw counts)
        
    Returns
    -------
    AnnData
        Integrated AnnData with genes and viruses
    """
    logger.info("Integrating viral counts with gene expression data...")
    
    # Store original gene names for later filtering
    gene_names = list(adata_pp.var_names)
    
    # Revert log1p transformation: exp(x) - 1
    logger.info("  Reverting log normalization...")
    adata_pp_denorm = adata_pp.copy()
    adata_pp_denorm.X = np.expm1(adata_pp_denorm.X)
    
    # Filter viral data to cells present in adata_pp
    common_cells = list(set(adata_pp.obs_names) & set(adata_v.obs_names))
    logger.info(f"  Common cells: {len(common_cells)}")
    
    if len(common_cells) == 0:
        logger.warning("No common cells found between expression and viral data!")
        return adata_pp
    
    adata_v_filtered = adata_v[adata_v.obs_names.isin(common_cells)].copy()
    
    # Filter to cells with viral reads
    sc.pp.filter_cells(adata_v_filtered, min_counts=1)
    logger.info(f"  Cells with viral reads: {adata_v_filtered.n_obs}")
    
    # Filter expression data to cells with viral reads
    adata_pp_filtered = adata_pp_denorm[adata_pp_denorm.obs_names.isin(adata_v_filtered.obs_names)].copy()
    
    # Clear varm to avoid concatenation issues
    adata_pp_filtered.varm = {}
    adata_v_filtered.varm = {}
    
    # Join the data
    adata_joined = ad.concat(
        [adata_pp_filtered, adata_v_filtered],
        join='outer',
        axis=1,
        merge='same'
    )
    
    # Reorder to match adata_pp
    common_cells_ordered = [c for c in adata_pp.obs_names if c in adata_joined.obs_names]
    adata_joined = adata_joined[common_cells_ordered, :].copy()
    
    # Transfer metadata
    adata_joined.obs = adata_pp.obs.loc[common_cells_ordered].copy()
    
    # Re-normalize the joined data
    logger.info("  Re-normalizing joined dataset...")
    sc.pp.normalize_total(adata_joined, target_sum=1e4)
    sc.pp.log1p(adata_joined)
    
    logger.info(f"  Integrated data: {adata_joined.n_obs} cells, {adata_joined.n_vars} features")
    
    # Store gene names for filtering markers later
    adata_joined.uns['original_gene_names'] = gene_names
    
    return adata_joined


# =============================================================================
# Analysis Functions
# =============================================================================

def calculate_viral_markers(adata_joined, human_virus_names, groupby='final_annotation'):
    """
    Calculate differential expression of viruses across cell types.
    
    Parameters
    ----------
    adata_joined : AnnData
        Integrated AnnData with genes and viruses
    human_virus_names : tuple
        Tuple of human virus names
    groupby : str
        Column in obs to group by
        
    Returns
    -------
    pd.DataFrame
        DataFrame with viral markers per cell type
    """
    logger.info(f"Calculating viral markers by {groupby}...")
    
    # Rank genes
    sc.tl.rank_genes_groups(adata_joined, groupby=groupby, 
                            key_added='rank_genes', method='wilcoxon')
    
    results = adata_joined.uns['rank_genes']
    
    # Get original gene names to filter out
    remove_list = adata_joined.uns.get('original_gene_names', [])
    
    # Build markers dataframe
    out = []
    for group in results['names'].dtype.names:
        for i in range(len(results['names'][group])):
            out.append({
                'virus': results['names'][group][i],
                'scores': results['scores'][group][i],
                'pval_adj': results['pvals_adj'][group][i],
                'lfc': results['logfoldchanges'][group][i],
                'cluster': group
            })
    
    markers = pd.DataFrame(out)
    
    # Filter to significant markers
    markers = markers[markers['scores'] > 0.0001]
    
    # Remove gene markers (keep only viruses)
    if remove_list:
        mask = markers['virus'].isin(remove_list)
        markers = markers[~mask]
    
    # Format p-values
    markers['pval_adj'] = markers['pval_adj'].apply(lambda x: f'{x:.2e}')
    
    # Filter to human viruses
    virus_markers = markers[markers['virus'].isin(human_virus_names)]
    
    logger.info(f"  Found {len(virus_markers)} viral marker entries")
    
    return markers, virus_markers


def aggregate_virus_scores(virus_markers):
    """
    Aggregate virus scores across all cell types.
    
    Parameters
    ----------
    virus_markers : pd.DataFrame
        Viral markers dataframe
        
    Returns
    -------
    pd.Series
        Aggregated scores per virus, sorted descending
    """
    if len(virus_markers) == 0:
        return pd.Series(dtype=float)
    
    subset = virus_markers[['virus', 'scores']].copy()
    subset['scores'] = pd.to_numeric(subset['scores'])
    aggregated = subset.groupby('virus')['scores'].sum().sort_values(ascending=False)
    
    return aggregated


# =============================================================================
# Visualization Functions
# =============================================================================

def generate_viral_plots(adata_joined, adata_pp, virus_scores, figures_dir, top_n=10):
    """
    Generate comprehensive viral detection plots.
    
    Parameters
    ----------
    adata_joined : AnnData
        Integrated AnnData
    adata_pp : AnnData
        Original expression AnnData
    virus_scores : pd.Series
        Aggregated virus scores
    figures_dir : str
        Output directory for figures
    top_n : int
        Number of top viruses to show
    """
    os.makedirs(figures_dir, exist_ok=True)
    
    if len(virus_scores) == 0:
        logger.warning("No viruses detected - skipping visualization")
        return
    
    logger.info(f"Generating viral detection plots in {figures_dir}...")
    
    # Ensure gene_symbol index
    if 'gene_symbol' in adata_joined.var.columns:
        adata_joined.var.index = adata_joined.var['gene_symbol']
    
    top_viruses = list(virus_scores.head(top_n).index)
    
    # 1. Matrix plot - Top viruses by cell type
    logger.info("  Creating matrix plot...")
    try:
        plt.rcParams['figure.figsize'] = (10, 8)
        sc.settings.set_figure_params(scanpy=True, fontsize=14)
        
        fig = sc.pl.matrixplot(
            adata_joined,
            top_viruses,
            'final_annotation',
            dendrogram=True,
            var_group_rotation=30,
            cmap='plasma',
            log=True,
            return_fig=True,
            show=False
        )
        fig.savefig(os.path.join(figures_dir, 'virus_matrix_plot.pdf'), 
                    bbox_inches='tight', dpi=300)
        plt.close()
    except Exception as e:
        logger.warning(f"  Matrix plot failed: {e}")
    
    # 2. UMAP colored by top virus
    logger.info("  Creating UMAP plots...")
    top_virus = top_viruses[0] if top_viruses else None
    
    if top_virus and top_virus in adata_joined.var_names:
        try:
            sc.set_figure_params(scanpy=True, fontsize=16)
            
            fig = sc.pl.umap(
                adata_joined,
                color=top_virus,
                use_raw=False,
                size=5,
                color_map='plasma',
                frameon=False,
                title=f'Top Virus: {top_virus}',
                return_fig=True,
                show=False
            )
            fig.savefig(os.path.join(figures_dir, f'umap_top_virus.pdf'),
                        bbox_inches='tight', dpi=300)
            plt.close()
        except Exception as e:
            logger.warning(f"  UMAP plot failed: {e}")
    
    # 3. Bar plot of virus scores
    logger.info("  Creating score bar plot...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        virus_scores.head(top_n).plot(kind='barh', ax=ax)
        ax.set_xlabel('Aggregate Score', fontsize=14)
        ax.set_ylabel('Virus', fontsize=14)
        ax.set_title(f'Top {top_n} Detected Viruses', fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'virus_score_barplot.pdf'),
                    bbox_inches='tight', dpi=300)
        plt.close()
    except Exception as e:
        logger.warning(f"  Bar plot failed: {e}")
    
    # 4. Violin plot for top viruses
    logger.info("  Creating violin plots...")
    for virus in top_viruses[:3]:
        if virus in adata_joined.var_names:
            try:
                fig = sc.pl.violin(
                    adata_joined,
                    virus,
                    groupby='final_annotation',
                    rotation=90,
                    return_fig=True,
                    show=False
                )
                safe_name = virus.replace(' ', '_').replace('/', '_')
                fig.savefig(os.path.join(figures_dir, f'violin_{safe_name}.pdf'),
                            bbox_inches='tight', dpi=300)
                plt.close()
            except Exception as e:
                logger.warning(f"  Violin plot for {virus} failed: {e}")
    
    plt.rcdefaults()
    logger.info("  Plots saved successfully")


def add_top_virus_to_adata(adata_pp, adata_joined, virus_scores):
    """
    Add top virus counts to the main expression AnnData.
    
    Parameters
    ----------
    adata_pp : AnnData
        Main expression AnnData
    adata_joined : AnnData
        Integrated AnnData with viral counts
    virus_scores : pd.Series
        Aggregated virus scores
        
    Returns
    -------
    AnnData
        Updated adata_pp with top virus counts
    """
    if len(virus_scores) == 0:
        logger.info("No viruses to add to adata_pp")
        return adata_pp
    
    top_virus = virus_scores.index[0]
    logger.info(f"Adding top virus ({top_virus}) counts to adata_pp...")
    
    if top_virus in adata_joined.var_names:
        # Get virus counts
        virus_counts = adata_joined[:, top_virus].X.toarray().flatten()
        virus_series = pd.Series(virus_counts, index=adata_joined.obs_names, name=top_virus)
        
        # Map to adata_pp
        safe_name = top_virus.replace(' ', '_').replace('/', '_')
        adata_pp.obs[safe_name] = adata_pp.obs_names.map(virus_series).fillna(0)
        
        n_positive = (adata_pp.obs[safe_name] > 0).sum()
        logger.info(f"  Added column '{safe_name}'")
        logger.info(f"  Cells with virus: {n_positive} ({100*n_positive/len(adata_pp):.1f}%)")
    
    return adata_pp


# =============================================================================
# Main Pipeline Function
# =============================================================================

def run_viral_integration(
    adata_pp_path,
    viral_sample_dirs,
    sample_ids,
    human_viral_hierarchy_path,
    output_dir,
):
    """
    Run the complete viral integration pipeline.
    
    Parameters
    ----------
    adata_pp_path : str
        Path to preprocessed expression H5AD
    viral_sample_dirs : list
        List of sample viral detection directories
    sample_ids : list
        List of sample identifiers
    human_viral_hierarchy_path : str
        Path to human viral database inspect.txt file
    output_dir : str
        Output directory
        
    Returns
    -------
    tuple
        (adata_pp_with_virus, adata_joined, summary)
    """
    os.makedirs(output_dir, exist_ok=True)
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    logger.info("="*60)
    logger.info("Viral Integration Pipeline")
    logger.info("="*60)
    
    # Load expression data
    logger.info(f"\nLoading expression data: {adata_pp_path}")
    adata_pp = sc.read_h5ad(adata_pp_path)
    logger.info(f"  Loaded: {adata_pp.n_obs} cells, {adata_pp.n_vars} genes")
    
    # Parse human viral hierarchy
    v_hierarchy_human, species_list = parse_kraken_hierarchy(human_viral_hierarchy_path)
    human_virus_names = tuple(entry['name'] for entry in species_list)
    logger.info(f"  Human virus species: {len(human_virus_names)}")
    
    # Get reference hierarchy from first sample
    first_sample_dir = viral_sample_dirs[0]
    hierarchy_file = os.path.join(first_sample_dir, 'kraken2_filtered_feature_bc_matrix', 'hierarchy.txt')
    
    if os.path.exists(hierarchy_file):
        sample_hierarchy = parse_sample_hierarchy(hierarchy_file)
        # Ensure all samples have all viruses
        ensure_all_viruses_present(viral_sample_dirs, sample_ids, sample_hierarchy)
    else:
        sample_hierarchy = {}
    
    # Load viral data
    adata_v = load_viral_data_for_samples(viral_sample_dirs, sample_ids, sample_hierarchy)
    
    if adata_v.n_obs == 0:
        logger.warning("No viral data loaded - saving empty results")
        summary = {
            'status': 'no_viral_data',
            'total_cells': adata_pp.n_obs,
            'cells_with_virus': 0,
            'viruses_detected': 0
        }
        pd.DataFrame([summary]).to_csv(
            os.path.join(output_dir, 'viral_integration_summary.tsv'),
            sep='\t', index=False
        )
        adata_pp = clean_obs_for_saving(adata_pp)
        adata_pp.write(os.path.join(output_dir, 'adata_with_virus.h5ad'))
        return adata_pp, None, summary
    
    # Save raw viral data
    adata_v = clean_obs_for_saving(adata_v)
    adata_v.write(os.path.join(output_dir, 'adata_viral_raw.h5ad'))
    
    # Filter to human viruses
    adata_v = filter_human_viruses(adata_v, human_virus_names)
    
    if adata_v.n_vars == 0:
        logger.warning("No human viruses detected after filtering")
        summary = {
            'status': 'no_human_viruses',
            'total_cells': adata_pp.n_obs,
            'cells_with_virus': 0,
            'viruses_detected': 0
        }
        pd.DataFrame([summary]).to_csv(
            os.path.join(output_dir, 'viral_integration_summary.tsv'),
            sep='\t', index=False
        )
        adata_pp = clean_obs_for_saving(adata_pp)
        adata_pp.write(os.path.join(output_dir, 'adata_with_virus.h5ad'))
        return adata_pp, None, summary
    
    # Integrate viral with expression
    adata_joined = integrate_viral_with_expression(adata_pp, adata_v)
    
    # Calculate neighborhood and clusters for integrated data
    logger.info("Computing UMAP for integrated data...")
    sc.pp.neighbors(adata_joined, n_pcs=min(50, adata_joined.n_obs - 1))
    sc.tl.umap(adata_joined)
    sc.tl.leiden(adata_joined, key_added='clusters', resolution=1, random_state=42)
    
    # Apply batch correction if sample_id exists
    if 'sample_id' in adata_joined.obs.columns:
        try:
            sc.external.pp.bbknn(adata_joined, batch_key='sample_id')
            sc.tl.umap(adata_joined)
        except Exception as e:
            logger.warning(f"BBKNN batch correction failed: {e}")
    
    # Calculate viral markers
    markers, virus_markers = calculate_viral_markers(
        adata_joined, human_virus_names, groupby='final_annotation'
    )
    adata_joined.uns['markers'] = markers
    
    # Aggregate scores
    virus_scores = aggregate_virus_scores(virus_markers)
    
    if len(virus_scores) > 0:
        logger.info(f"\nTop 10 detected viruses:")
        for i, (virus, score) in enumerate(virus_scores.head(10).items(), 1):
            logger.info(f"  {i}. {virus}: {score:.4f}")
    
    # Generate plots
    generate_viral_plots(adata_joined, adata_pp, virus_scores, figures_dir)
    
    # Add top virus to adata_pp
    adata_pp = add_top_virus_to_adata(adata_pp, adata_joined, virus_scores)
    
    # Save results
    logger.info("\nSaving results...")
    adata_joined = clean_obs_for_saving(adata_joined)
    adata_joined.write(os.path.join(output_dir, 'adata_viral_integrated.h5ad'))
    
    adata_pp = clean_obs_for_saving(adata_pp)
    adata_pp.write(os.path.join(output_dir, 'adata_with_virus.h5ad'))
    
    # Summary
    summary = {
        'status': 'success',
        'total_cells': adata_pp.n_obs,
        'cells_with_viral_data': adata_joined.n_obs,
        'human_viruses_detected': len(virus_scores),
        'top_virus': virus_scores.index[0] if len(virus_scores) > 0 else None,
        'top_virus_score': virus_scores.iloc[0] if len(virus_scores) > 0 else 0,
    }
    
    pd.DataFrame([summary]).to_csv(
        os.path.join(output_dir, 'viral_integration_summary.tsv'),
        sep='\t', index=False
    )
    
    # Save virus detection summary
    if len(virus_scores) > 0:
        virus_scores.to_csv(
            os.path.join(output_dir, 'virus_scores.tsv'),
            sep='\t', header=['score']
        )
    
    logger.info("\n" + "="*60)
    logger.info("Pipeline completed successfully!")
    logger.info("="*60)
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Cells with viral data: {adata_joined.n_obs}")
    logger.info(f"Human viruses detected: {len(virus_scores)}")
    
    return adata_pp, adata_joined, summary


# =============================================================================
# Snakemake Integration
# =============================================================================

def run_from_snakemake():
    """Run from Snakemake rule."""
    
    # Get inputs
    adata_pp_path = snakemake.input.adata_pp
    viral_summary = snakemake.input.viral_summary  # This ensures viral detection completed
    
    # Get outputs
    output_dir = os.path.dirname(snakemake.output.adata)
    
    # Get params
    sample_ids = snakemake.params.sample_ids
    human_viral_db = snakemake.params.human_viral_db
    viral_base_dir = snakemake.params.viral_base_dir
    
    # Build sample directories
    viral_sample_dirs = [os.path.join(viral_base_dir, sid) for sid in sample_ids]
    
    # Set up logging
    if snakemake.log:
        file_handler = logging.FileHandler(snakemake.log[0])
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Run pipeline
    run_viral_integration(
        adata_pp_path=adata_pp_path,
        viral_sample_dirs=viral_sample_dirs,
        sample_ids=sample_ids,
        human_viral_hierarchy_path=human_viral_db,
        output_dir=output_dir,
    )


# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Integrate viral detection with expression data'
    )
    parser.add_argument('--adata-pp', required=True, help='Preprocessed expression H5AD')
    parser.add_argument('--viral-dir', required=True, help='Base directory with viral detection results')
    parser.add_argument('--sample-ids', nargs='+', required=True, help='Sample identifiers')
    parser.add_argument('--human-viral-db', required=True, 
                        help='Path to human viral Kraken2 database inspect.txt')
    parser.add_argument('--output', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    viral_sample_dirs = [os.path.join(args.viral_dir, sid) for sid in args.sample_ids]
    
    run_viral_integration(
        adata_pp_path=args.adata_pp,
        viral_sample_dirs=viral_sample_dirs,
        sample_ids=args.sample_ids,
        human_viral_hierarchy_path=args.human_viral_db,
        output_dir=args.output,
    )


if __name__ == '__main__':
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        main()
