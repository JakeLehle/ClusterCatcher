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

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings
warnings.filterwarnings('ignore')

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
    
    if rename_map:
        df = df.rename(columns=rename_map)
        logger.debug(f"Cleaned {len(rename_map)} column names")
    
    return df


def clean_obs_for_saving(adata):
    """
    Clean AnnData obs columns for HDF5 saving.
    
    Converts problematic columns to strings to avoid serialization issues.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object to clean
        
    Returns
    -------
    AnnData
        Cleaned AnnData
    """
    adata.obs = clean_column_names(adata.obs)
    
    for col in adata.obs.columns:
        try:
            # Try to convert to string if it's a problematic type
            if adata.obs[col].dtype == object:
                adata.obs[col] = adata.obs[col].astype(str)
        except:
            pass
    
    if 'gene_ids' in adata.var.columns:
        adata.var['gene_ids'] = adata.var['gene_ids'].astype(str)
    
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
    
    if not os.path.exists(hierarchy_file):
        logger.warning(f"Hierarchy file not found: {hierarchy_file}")
        return v_hierarchy, species_list
    
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


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_viral_adata(viral_adata_path):
    """
    Load the combined viral counts AnnData from summarize step.
    
    Parameters
    ----------
    viral_adata_path : str
        Path to viral_counts.h5ad from summarize_viral step
        
    Returns
    -------
    AnnData
        Viral counts AnnData
    """
    logger.info(f"Loading viral counts AnnData: {viral_adata_path}")
    
    if not os.path.exists(viral_adata_path):
        logger.warning(f"Viral AnnData not found: {viral_adata_path}")
        return ad.AnnData()
    
    adata_v = sc.read_h5ad(viral_adata_path)
    logger.info(f"  Loaded: {adata_v.n_obs} cells, {adata_v.n_vars} organisms")
    
    return adata_v


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
    human_species_names : tuple or list
        Human virus species names
        
    Returns
    -------
    AnnData
        Filtered AnnData with only human viruses
    """
    logger.info("Filtering for human-associated viruses...")
    
    if adata_v.n_vars == 0:
        return adata_v
    
    human_virus_names = tuple(name.strip() for name in human_species_names if name and name.strip())
    
    if not human_virus_names:
        logger.warning("No human virus names provided for filtering")
        return adata_v
    
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
    
    if adata_v.n_obs == 0 or adata_v.n_vars == 0:
        logger.warning("No viral data to integrate")
        return adata_pp.copy()
    
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
        return adata_pp.copy()
    
    adata_v_filtered = adata_v[adata_v.obs_names.isin(common_cells)].copy()
    
    # Filter to cells with viral reads
    cell_sums = np.array(adata_v_filtered.X.sum(axis=1)).flatten()
    cells_with_reads = cell_sums > 0
    adata_v_filtered = adata_v_filtered[cells_with_reads].copy()
    logger.info(f"  Cells with viral reads: {adata_v_filtered.n_obs}")
    
    if adata_v_filtered.n_obs == 0:
        logger.warning("No cells have viral reads")
        return adata_pp.copy()
    
    # Filter expression data to cells with viral reads
    adata_pp_filtered = adata_pp_denorm[adata_pp_denorm.obs_names.isin(adata_v_filtered.obs_names)].copy()
    
    # Clear varm to avoid concatenation issues
    adata_pp_filtered.varm = {}
    adata_v_filtered.varm = {}
    
    # Ensure obs columns are compatible
    for col in adata_pp_filtered.obs.columns:
        if col not in adata_v_filtered.obs.columns:
            adata_v_filtered.obs[col] = adata_pp_filtered.obs.loc[adata_v_filtered.obs_names, col].values
    
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
    human_virus_names : tuple or list
        Human virus names
    groupby : str
        Column in obs to group by
        
    Returns
    -------
    tuple
        (all_markers_df, virus_markers_df)
    """
    logger.info(f"Calculating viral markers by {groupby}...")
    
    if groupby not in adata_joined.obs.columns:
        logger.warning(f"Column '{groupby}' not found in obs")
        return pd.DataFrame(), pd.DataFrame()
    
    # Check for sufficient groups
    n_groups = adata_joined.obs[groupby].nunique()
    if n_groups < 2:
        logger.warning(f"Need at least 2 groups for marker analysis, found {n_groups}")
        return pd.DataFrame(), pd.DataFrame()
    
    try:
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
        
    except Exception as e:
        logger.warning(f"Marker calculation failed: {e}")
        return pd.DataFrame(), pd.DataFrame()


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
    subset['scores'] = pd.to_numeric(subset['scores'], errors='coerce')
    aggregated = subset.groupby('virus')['scores'].sum().sort_values(ascending=False)
    
    return aggregated


# =============================================================================
# Visualization Functions
# =============================================================================

def generate_viral_plots(adata_joined, virus_scores, figures_dir, top_n=10):
    """
    Generate comprehensive viral detection plots.
    
    Parameters
    ----------
    adata_joined : AnnData
        Integrated AnnData
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
    top_viruses = [v for v in top_viruses if v in adata_joined.var_names]
    
    if not top_viruses:
        logger.warning("No top viruses found in adata_joined")
        return
    
    # 1. Matrix plot - Top viruses by cell type
    if 'final_annotation' in adata_joined.obs.columns:
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
            fig.savefig(os.path.join(figures_dir, 'virus_matrix_plot.png'), 
                        dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            logger.warning(f"  Matrix plot failed: {e}")
    
    # 2. UMAP colored by top virus
    logger.info("  Creating UMAP plots...")
    top_virus = top_viruses[0] if top_viruses else None
    
    if top_virus and 'X_umap' in adata_joined.obsm:
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
            fig.savefig(os.path.join(figures_dir, 'umap_top_virus.pdf'),
                        bbox_inches='tight', dpi=300)
            fig.savefig(os.path.join(figures_dir, 'umap_top_virus.png'),
                        dpi=150, bbox_inches='tight')
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
        plt.savefig(os.path.join(figures_dir, 'virus_score_barplot.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()
    except Exception as e:
        logger.warning(f"  Bar plot failed: {e}")
    
    # 4. Violin plot for top viruses
    if 'final_annotation' in adata_joined.obs.columns:
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
                    fig.savefig(os.path.join(figures_dir, f'violin_{safe_name}.png'),
                                dpi=150, bbox_inches='tight')
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
        virus_data = adata_joined[:, top_virus].X
        if hasattr(virus_data, 'toarray'):
            virus_counts = virus_data.toarray().flatten()
        else:
            virus_counts = np.array(virus_data).flatten()
        
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
    viral_adata_path,
    human_viral_hierarchy_path,
    output_adata_path,
    output_integrated_path,
    output_summary_path,
    figures_dir,
):
    """
    Run the complete viral integration pipeline.
    
    Parameters
    ----------
    adata_pp_path : str
        Path to preprocessed expression H5AD
    viral_adata_path : str
        Path to viral_counts.h5ad from summarize step
    human_viral_hierarchy_path : str
        Path to human viral database inspect.txt file
    output_adata_path : str
        Path for output adata_with_virus.h5ad
    output_integrated_path : str
        Path for output adata_viral_integrated.h5ad
    output_summary_path : str
        Path for output summary TSV
    figures_dir : str
        Output directory for figures
        
    Returns
    -------
    tuple
        (adata_pp_with_virus, adata_joined, summary)
    """
    os.makedirs(os.path.dirname(output_adata_path), exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    
    logger.info("="*60)
    logger.info("Viral Integration Pipeline")
    logger.info("="*60)
    
    # Load expression data
    logger.info(f"\nLoading expression data: {adata_pp_path}")
    adata_pp = sc.read_h5ad(adata_pp_path)
    logger.info(f"  Loaded: {adata_pp.n_obs} cells, {adata_pp.n_vars} genes")
    
    # Parse human viral hierarchy
    human_virus_names = []
    if human_viral_hierarchy_path and os.path.exists(human_viral_hierarchy_path):
        v_hierarchy_human, species_list = parse_kraken_hierarchy(human_viral_hierarchy_path)
        human_virus_names = [entry['name'] for entry in species_list]
        logger.info(f"  Human virus species: {len(human_virus_names)}")
    else:
        logger.warning("Human viral hierarchy not provided or not found")
    
    # Load viral counts AnnData
    adata_v = load_viral_adata(viral_adata_path)
    
    if adata_v.n_obs == 0:
        logger.warning("No viral data loaded - saving empty results")
        summary = {
            'status': 'no_viral_data',
            'total_cells': adata_pp.n_obs,
            'cells_with_virus': 0,
            'viruses_detected': 0
        }
        pd.DataFrame([summary]).to_csv(output_summary_path, sep='\t', index=False)
        adata_pp = clean_obs_for_saving(adata_pp)
        adata_pp.write(output_adata_path)
        # Create empty integrated file
        ad.AnnData().write(output_integrated_path)
        return adata_pp, None, summary
    
    # Filter to human viruses if we have the list
    if human_virus_names:
        adata_v = filter_human_viruses(adata_v, human_virus_names)
    
    if adata_v.n_vars == 0:
        logger.warning("No human viruses detected after filtering")
        summary = {
            'status': 'no_human_viruses',
            'total_cells': adata_pp.n_obs,
            'cells_with_virus': 0,
            'viruses_detected': 0
        }
        pd.DataFrame([summary]).to_csv(output_summary_path, sep='\t', index=False)
        adata_pp = clean_obs_for_saving(adata_pp)
        adata_pp.write(output_adata_path)
        ad.AnnData().write(output_integrated_path)
        return adata_pp, None, summary
    
    # Integrate viral with expression
    adata_joined = integrate_viral_with_expression(adata_pp, adata_v)
    
    # Compute UMAP if not present
    if 'X_umap' not in adata_joined.obsm:
        logger.info("Computing UMAP for integrated data...")
        try:
            n_pcs = min(50, adata_joined.n_obs - 1, adata_joined.n_vars - 1)
            sc.pp.pca(adata_joined, n_comps=n_pcs)
            sc.pp.neighbors(adata_joined, n_pcs=n_pcs)
            sc.tl.umap(adata_joined)
        except Exception as e:
            logger.warning(f"UMAP computation failed: {e}")
    
    # Calculate viral markers
    markers, virus_markers = calculate_viral_markers(
        adata_joined, human_virus_names, groupby='final_annotation'
    )
    if len(markers) > 0:
        adata_joined.uns['markers'] = markers
    
    # Aggregate scores
    virus_scores = aggregate_virus_scores(virus_markers)
    
    if len(virus_scores) > 0:
        logger.info(f"\nTop 10 detected viruses:")
        for i, (virus, score) in enumerate(virus_scores.head(10).items(), 1):
            logger.info(f"  {i}. {virus}: {score:.4f}")
    
    # Generate plots
    generate_viral_plots(adata_joined, virus_scores, figures_dir)
    
    # Add top virus to adata_pp
    adata_pp = add_top_virus_to_adata(adata_pp, adata_joined, virus_scores)
    
    # Save results
    logger.info("\nSaving results...")
    adata_joined = clean_obs_for_saving(adata_joined)
    adata_joined.write(output_integrated_path)
    
    adata_pp = clean_obs_for_saving(adata_pp)
    adata_pp.write(output_adata_path)
    
    # Summary
    summary = {
        'status': 'success',
        'total_cells': adata_pp.n_obs,
        'cells_with_viral_data': adata_joined.n_obs,
        'human_viruses_detected': len(virus_scores),
        'top_virus': virus_scores.index[0] if len(virus_scores) > 0 else None,
        'top_virus_score': float(virus_scores.iloc[0]) if len(virus_scores) > 0 else 0,
    }
    
    pd.DataFrame([summary]).to_csv(output_summary_path, sep='\t', index=False)
    
    # Save virus detection summary
    if len(virus_scores) > 0:
        virus_scores.to_csv(
            os.path.join(os.path.dirname(output_summary_path), 'virus_scores.tsv'),
            sep='\t', header=['score']
        )
    
    logger.info("\n" + "="*60)
    logger.info("Pipeline completed successfully!")
    logger.info("="*60)
    logger.info(f"Output directory: {os.path.dirname(output_adata_path)}")
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
    viral_adata_path = snakemake.input.viral_adata
    
    # Get outputs
    output_adata = snakemake.output.adata
    output_integrated = snakemake.output.integrated
    output_summary = snakemake.output.summary
    output_figures = snakemake.output.figures
    
    # Get params
    human_viral_db = snakemake.params.human_viral_db
    
    # Set up logging
    if snakemake.log:
        os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
        file_handler = logging.FileHandler(snakemake.log[0])
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Run pipeline
    run_viral_integration(
        adata_pp_path=adata_pp_path,
        viral_adata_path=viral_adata_path,
        human_viral_hierarchy_path=human_viral_db,
        output_adata_path=output_adata,
        output_integrated_path=output_integrated,
        output_summary_path=output_summary,
        figures_dir=output_figures,
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
    parser.add_argument('--viral-adata', required=True, help='Viral counts H5AD from summarize step')
    parser.add_argument('--human-viral-db', help='Path to human viral Kraken2 database inspect.txt')
    parser.add_argument('--output', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    output_dir = args.output
    
    run_viral_integration(
        adata_pp_path=args.adata_pp,
        viral_adata_path=args.viral_adata,
        human_viral_hierarchy_path=args.human_viral_db,
        output_adata_path=os.path.join(output_dir, 'adata_with_virus.h5ad'),
        output_integrated_path=os.path.join(output_dir, 'adata_viral_integrated.h5ad'),
        output_summary_path=os.path.join(output_dir, 'viral_integration_summary.tsv'),
        figures_dir=os.path.join(output_dir, 'figures'),
    )


if __name__ == '__main__':
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        main()

