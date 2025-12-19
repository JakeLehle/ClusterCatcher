#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scanpy_qc_annotation.py
=======================

Comprehensive QC, preprocessing, and cell type annotation pipeline using
Scanpy and popV.

This script performs:
1. Loading Cell Ranger count matrices
2. Quality control filtering (genes, UMIs, mitochondrial content)
3. Doublet detection and removal (Scrublet)
4. Normalization and log transformation
5. Highly variable gene selection
6. Dimensionality reduction (PCA, UMAP)
7. Clustering (Leiden)
8. Automated cell type annotation (popV)
9. Visualization and QC reports

Usage:
    Called via Snakemake rule with snakemake.input/output/params
    
    Or standalone:
    python scanpy_qc_annotation.py --config config.yaml

Requirements:
    - scanpy
    - scrublet
    - popv
    - pandas, numpy, matplotlib, seaborn
"""

import os
import sys
import yaml
import pickle
import logging
import argparse
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration and Defaults
# =============================================================================

DEFAULT_QC_PARAMS = {
    'min_genes': 200,
    'max_genes': 5000,
    'min_cells': 3,
    'max_pct_mito': 20,
    'min_counts': 500,
    'max_counts': 50000,
    'doublet_removal': True,
    'doublet_rate': 0.06,  # Expected doublet rate
}

DEFAULT_PREPROCESSING_PARAMS = {
    'target_sum': 1e4,      # Normalization target
    'n_top_genes': 2000,    # Number of HVGs
    'n_pcs': 50,            # PCA components
    'n_neighbors': 15,      # For neighbor graph
    'leiden_resolution': 1.0,
}

DEFAULT_ANNOTATION_PARAMS = {
    'method': 'popv',
    'popv_model': 'popv_immune_All_human_umap_from_cellxgene',  # Default model
    'min_confidence': 0.5,
    'batch_key': 'sample_id',
}


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_cellranger_outputs(cellranger_dirs, sample_ids):
    """
    Load Cell Ranger count matrices from multiple samples.
    
    Parameters
    ----------
    cellranger_dirs : list
        List of paths to Cell Ranger output directories
    sample_ids : list
        List of sample identifiers
        
    Returns
    -------
    AnnData
        Concatenated AnnData object with all samples
    """
    logger.info(f"Loading {len(sample_ids)} samples...")
    
    adatas = []
    for sample_id, cr_dir in zip(sample_ids, cellranger_dirs):
        # Find the matrix file
        matrix_h5 = os.path.join(cr_dir, 'outs', 'filtered_feature_bc_matrix.h5')
        matrix_dir = os.path.join(cr_dir, 'outs', 'filtered_feature_bc_matrix')
        
        if os.path.exists(matrix_h5):
            logger.info(f"  Loading {sample_id} from H5 file...")
            adata = sc.read_10x_h5(matrix_h5)
        elif os.path.exists(matrix_dir):
            logger.info(f"  Loading {sample_id} from matrix directory...")
            adata = sc.read_10x_mtx(matrix_dir)
        else:
            logger.warning(f"  Could not find matrix for {sample_id}, skipping...")
            continue
        
        # Make var names unique
        adata.var_names_make_unique()
        
        # Add sample metadata
        adata.obs['sample_id'] = sample_id
        adata.obs['original_barcode'] = adata.obs.index
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
        
        # Store gene symbols if available
        if 'gene_ids' in adata.var.columns:
            adata.var['gene_ids'] = adata.var['gene_ids']
        if 'feature_types' in adata.var.columns:
            adata.var['feature_types'] = adata.var['feature_types']
        
        # Add gene_symbol column for compatibility
        adata.var['gene_symbol'] = adata.var_names
        
        logger.info(f"    {sample_id}: {adata.n_obs} cells, {adata.n_vars} genes")
        adatas.append(adata)
    
    if len(adatas) == 0:
        raise ValueError("No samples were loaded successfully!")
    
    # Concatenate all samples
    logger.info("Concatenating samples...")
    adata = sc.concat(adatas, join='outer', label='sample_id', keys=[a.obs['sample_id'].iloc[0] for a in adatas])
    
    # Fix the sample_id column after concat
    adata.obs['sample_id'] = adata.obs['sample_id'].astype('category')
    
    logger.info(f"Total: {adata.n_obs} cells, {adata.n_vars} genes from {len(adatas)} samples")
    
    return adata


def load_from_srascraper_dict(gse_dict, cellranger_base_dir):
    """
    Load Cell Ranger outputs using SRAscraper dictionary structure.
    
    Parameters
    ----------
    gse_dict : dict
        SRAscraper-style dictionary with sample metadata
    cellranger_base_dir : str
        Base directory containing Cell Ranger outputs
        
    Returns
    -------
    AnnData
        Concatenated AnnData object
    """
    logger.info("Loading from SRAscraper dictionary structure...")
    
    adatas = []
    for gse_key in gse_dict.keys():
        for idx, row in gse_dict[gse_key].iterrows():
            accession = row['run_accession']
            
            # Find Cell Ranger output directory
            # Try multiple possible directory structures
            possible_dirs = [
                os.path.join(cellranger_base_dir, gse_key, accession),
                os.path.join(cellranger_base_dir, accession),
            ]
            
            # Also check for Cell Ranger's typical output naming
            for base in possible_dirs:
                if os.path.exists(base):
                    # Look for subdirectories that might be the actual CR output
                    for item in os.listdir(base):
                        item_path = os.path.join(base, item)
                        if os.path.isdir(item_path) and os.path.exists(os.path.join(item_path, 'outs')):
                            possible_dirs.append(item_path)
            
            cr_dir = None
            for d in possible_dirs:
                outs_dir = os.path.join(d, 'outs')
                if os.path.exists(outs_dir):
                    cr_dir = d
                    break
            
            if cr_dir is None:
                logger.warning(f"Could not find Cell Ranger output for {accession}")
                continue
            
            # Load the matrix
            matrix_h5 = os.path.join(cr_dir, 'outs', 'filtered_feature_bc_matrix.h5')
            if os.path.exists(matrix_h5):
                logger.info(f"  Loading {accession}...")
                adata = sc.read_10x_h5(matrix_h5)
                adata.var_names_make_unique()
                
                # Add metadata
                adata.obs['sample_id'] = accession
                adata.obs['series_id'] = gse_key
                adata.obs['original_barcode'] = adata.obs.index
                adata.obs_names = [f"{accession}_{bc}" for bc in adata.obs_names]
                adata.var['gene_symbol'] = adata.var_names
                
                # Add any additional metadata from the dictionary
                for col in row.index:
                    if col not in ['run_accession']:
                        adata.obs[col] = row[col]
                
                logger.info(f"    {accession}: {adata.n_obs} cells")
                adatas.append(adata)
    
    if len(adatas) == 0:
        raise ValueError("No samples were loaded!")
    
    logger.info("Concatenating samples...")
    adata = sc.concat(adatas, join='outer')
    adata.obs['sample_id'] = adata.obs['sample_id'].astype('category')
    
    logger.info(f"Total: {adata.n_obs} cells from {len(adatas)} samples")
    
    return adata


# =============================================================================
# QC Functions
# =============================================================================

def calculate_qc_metrics(adata):
    """
    Calculate QC metrics for all cells.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object
        
    Returns
    -------
    AnnData
        AnnData with QC metrics added to obs
    """
    logger.info("Calculating QC metrics...")
    
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    
    # Ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL', 'Rps', 'Rpl'))
    
    # Hemoglobin genes
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]', case=False, regex=True)
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo', 'hb'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    logger.info(f"  Mean genes/cell: {adata.obs['n_genes_by_counts'].mean():.1f}")
    logger.info(f"  Mean counts/cell: {adata.obs['total_counts'].mean():.1f}")
    logger.info(f"  Mean % mito: {adata.obs['pct_counts_mt'].mean():.1f}%")
    
    return adata


def filter_cells_and_genes(adata, params):
    """
    Filter cells and genes based on QC thresholds.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData with QC metrics
    params : dict
        QC parameters
        
    Returns
    -------
    AnnData
        Filtered AnnData
    """
    logger.info("Filtering cells and genes...")
    
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=params['min_genes'])
    
    if params.get('min_counts'):
        adata = adata[adata.obs['total_counts'] >= params['min_counts'], :]
    
    if params.get('max_counts'):
        adata = adata[adata.obs['total_counts'] <= params['max_counts'], :]
    
    adata = adata[adata.obs['n_genes_by_counts'] <= params['max_genes'], :]
    adata = adata[adata.obs['pct_counts_mt'] <= params['max_pct_mito'], :]
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=params['min_cells'])
    
    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars
    
    logger.info(f"  Cells: {n_cells_before} -> {n_cells_after} ({n_cells_before - n_cells_after} removed)")
    logger.info(f"  Genes: {n_genes_before} -> {n_genes_after} ({n_genes_before - n_genes_after} removed)")
    
    return adata


def remove_doublets(adata, expected_doublet_rate=0.06):
    """
    Detect and remove doublets using Scrublet.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData
    expected_doublet_rate : float
        Expected doublet rate
        
    Returns
    -------
    AnnData
        AnnData with doublets removed
    """
    logger.info("Detecting doublets with Scrublet...")
    
    try:
        import scrublet as scr
    except ImportError:
        logger.warning("Scrublet not installed, skipping doublet detection")
        return adata
    
    n_before = adata.n_obs
    
    # Run Scrublet per sample for better accuracy
    doublet_scores = []
    predicted_doublets = []
    
    for sample in adata.obs['sample_id'].unique():
        sample_mask = adata.obs['sample_id'] == sample
        sample_adata = adata[sample_mask].copy()
        
        if sample_adata.n_obs < 100:
            logger.warning(f"  Sample {sample} has <100 cells, skipping doublet detection")
            doublet_scores.extend([0.0] * sample_adata.n_obs)
            predicted_doublets.extend([False] * sample_adata.n_obs)
            continue
        
        try:
            scrub = scr.Scrublet(sample_adata.X, expected_doublet_rate=expected_doublet_rate)
            scores, predictions = scrub.scrub_doublets(min_counts=2, min_cells=3, verbose=False)
            doublet_scores.extend(scores)
            predicted_doublets.extend(predictions)
            n_doublets = sum(predictions)
            logger.info(f"  {sample}: {n_doublets} doublets detected ({100*n_doublets/len(predictions):.1f}%)")
        except Exception as e:
            logger.warning(f"  Scrublet failed for {sample}: {e}")
            doublet_scores.extend([0.0] * sample_adata.n_obs)
            predicted_doublets.extend([False] * sample_adata.n_obs)
    
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets
    
    # Remove doublets
    adata = adata[~adata.obs['predicted_doublet'], :].copy()
    
    n_after = adata.n_obs
    logger.info(f"  Total doublets removed: {n_before - n_after}")
    
    return adata


# =============================================================================
# Preprocessing Functions
# =============================================================================

def preprocess_data(adata, params):
    """
    Normalize, find HVGs, and reduce dimensionality.
    
    Parameters
    ----------
    adata : AnnData
        QC-filtered AnnData
    params : dict
        Preprocessing parameters
        
    Returns
    -------
    AnnData
        Preprocessed AnnData
    """
    logger.info("Preprocessing data...")
    
    # Store raw counts
    adata.layers['counts'] = adata.X.copy()
    
    # Normalize
    logger.info("  Normalizing...")
    sc.pp.normalize_total(adata, target_sum=params['target_sum'])
    sc.pp.log1p(adata)
    
    # Store normalized data
    adata.raw = adata.copy()
    
    # Find highly variable genes
    logger.info(f"  Finding {params['n_top_genes']} highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=params['n_top_genes'],
        batch_key='sample_id' if 'sample_id' in adata.obs.columns else None,
        flavor='seurat_v3',
        layer='counts'
    )
    
    # Subset to HVGs for downstream analysis
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    
    # Scale
    logger.info("  Scaling...")
    sc.pp.scale(adata_hvg, max_value=10)
    
    # PCA
    logger.info(f"  Running PCA ({params['n_pcs']} components)...")
    sc.tl.pca(adata_hvg, n_comps=params['n_pcs'])
    
    # Copy PCA results back to full adata
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.uns['pca'] = adata_hvg.uns['pca']
    adata.varm['PCs'] = np.zeros((adata.n_vars, params['n_pcs']))
    adata.varm['PCs'][adata.var['highly_variable'], :] = adata_hvg.varm['PCs']
    
    # Neighbors
    logger.info(f"  Computing neighbors (k={params['n_neighbors']})...")
    sc.pp.neighbors(adata, n_neighbors=params['n_neighbors'], n_pcs=params['n_pcs'])
    
    # UMAP
    logger.info("  Computing UMAP...")
    sc.tl.umap(adata)
    
    # Leiden clustering
    logger.info(f"  Clustering (resolution={params['leiden_resolution']})...")
    sc.tl.leiden(adata, resolution=params['leiden_resolution'])
    
    return adata


# =============================================================================
# Annotation Functions
# =============================================================================

def annotate_with_popv(adata, model_name, batch_key='sample_id'):
    """
    Annotate cell types using popV.
    
    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData
    model_name : str
        Name of the popV model to use
    batch_key : str
        Key in obs for batch correction
        
    Returns
    -------
    AnnData
        AnnData with cell type annotations
    """
    logger.info(f"Running popV annotation with model: {model_name}")
    
    try:
        import popv
    except ImportError:
        logger.error("popV not installed! Install with: pip install popv")
        raise
    
    # Prepare data for popV
    # popV expects raw counts in a specific format
    adata_popv = adata.copy()
    
    # Make sure we have the counts layer
    if 'counts' in adata_popv.layers:
        adata_popv.X = adata_popv.layers['counts'].copy()
    
    # Run popV annotation
    logger.info("  Running popV prediction (this may take a while)...")
    
    try:
        # Query the model
        popv.annotation.annotate_data(
            adata_popv,
            model_name,
            batch_key=batch_key if batch_key in adata_popv.obs.columns else None,
        )
        
        # Transfer annotations back
        adata.obs['popv_prediction'] = adata_popv.obs['popv_prediction']
        
        if 'popv_prediction_score' in adata_popv.obs.columns:
            adata.obs['popv_prediction_score'] = adata_popv.obs['popv_prediction_score']
        
        # Use popV prediction as final annotation
        adata.obs['final_annotation'] = adata.obs['popv_prediction']
        
        # Summary statistics
        annotation_counts = adata.obs['final_annotation'].value_counts()
        logger.info("  Cell type distribution:")
        for ct, count in annotation_counts.head(10).items():
            logger.info(f"    {ct}: {count} ({100*count/adata.n_obs:.1f}%)")
        
    except Exception as e:
        logger.error(f"popV annotation failed: {e}")
        logger.info("Falling back to leiden clusters as annotations")
        adata.obs['final_annotation'] = adata.obs['leiden'].astype(str)
    
    return adata


# =============================================================================
# Visualization Functions
# =============================================================================

def generate_qc_plots(adata, output_dir):
    """
    Generate QC visualization plots.
    
    Parameters
    ----------
    adata : AnnData
        AnnData with QC metrics
    output_dir : str
        Directory to save plots
    """
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Generating QC plots in {output_dir}...")
    
    sc.settings.figdir = output_dir
    sc.settings.set_figure_params(dpi=150, frameon=False)
    
    # Violin plots of QC metrics
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    sc.pl.violin(adata, 'n_genes_by_counts', groupby='sample_id', ax=axes[0], show=False)
    axes[0].set_title('Genes per Cell')
    axes[0].tick_params(axis='x', rotation=90)
    
    sc.pl.violin(adata, 'total_counts', groupby='sample_id', ax=axes[1], show=False)
    axes[1].set_title('UMI Counts per Cell')
    axes[1].tick_params(axis='x', rotation=90)
    
    sc.pl.violin(adata, 'pct_counts_mt', groupby='sample_id', ax=axes[2], show=False)
    axes[2].set_title('% Mitochondrial')
    axes[2].tick_params(axis='x', rotation=90)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'qc_violin_plots.pdf'), bbox_inches='tight')
    plt.close()
    
    # Scatter plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt',
                  ax=axes[0], show=False)
    axes[0].set_title('Counts vs Genes (colored by % mito)')
    
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color='n_genes_by_counts',
                  ax=axes[1], show=False)
    axes[1].set_title('Counts vs % Mito (colored by genes)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'qc_scatter_plots.pdf'), bbox_inches='tight')
    plt.close()


def generate_annotation_plots(adata, output_dir):
    """
    Generate annotation visualization plots.
    
    Parameters
    ----------
    adata : AnnData
        Annotated AnnData
    output_dir : str
        Directory to save plots
    """
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Generating annotation plots in {output_dir}...")
    
    sc.settings.figdir = output_dir
    
    # UMAP with cell types
    fig = sc.pl.umap(
        adata,
        color='final_annotation',
        legend_loc='on data',
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        title='Cell Type Annotations',
        return_fig=True,
        show=False
    )
    plt.savefig(os.path.join(output_dir, 'umap_cell_types.pdf'), bbox_inches='tight')
    plt.close()
    
    # UMAP with samples
    sc.pl.umap(
        adata,
        color='sample_id',
        frameon=False,
        title='Samples',
        show=False
    )
    plt.savefig(os.path.join(output_dir, 'umap_samples.pdf'), bbox_inches='tight')
    plt.close()
    
    # UMAP with clusters
    sc.pl.umap(
        adata,
        color='leiden',
        legend_loc='on data',
        legend_fontsize=10,
        legend_fontoutline=2,
        frameon=False,
        title='Leiden Clusters',
        show=False
    )
    plt.savefig(os.path.join(output_dir, 'umap_leiden.pdf'), bbox_inches='tight')
    plt.close()
    
    # Cell type proportions per sample
    ct_props = pd.crosstab(
        adata.obs['sample_id'],
        adata.obs['final_annotation'],
        normalize='index'
    ) * 100
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ct_props.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
    ax.set_ylabel('Percentage')
    ax.set_xlabel('Sample')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cell_type_proportions.pdf'), bbox_inches='tight')
    plt.close()


def export_qc_metrics(adata, output_path):
    """
    Export QC metrics summary to TSV.
    
    Parameters
    ----------
    adata : AnnData
        AnnData with QC metrics
    output_path : str
        Output file path
    """
    logger.info(f"Exporting QC metrics to {output_path}")
    
    # Per-sample summary
    qc_summary = adata.obs.groupby('sample_id').agg({
        'n_genes_by_counts': ['count', 'mean', 'median', 'std'],
        'total_counts': ['mean', 'median', 'std'],
        'pct_counts_mt': ['mean', 'median', 'std'],
    })
    
    qc_summary.columns = ['_'.join(col).strip() for col in qc_summary.columns.values]
    qc_summary = qc_summary.rename(columns={'n_genes_by_counts_count': 'n_cells'})
    
    qc_summary.to_csv(output_path, sep='\t')


def export_annotation_summary(adata, output_path):
    """
    Export annotation summary to TSV.
    
    Parameters
    ----------
    adata : AnnData
        Annotated AnnData
    output_path : str
        Output file path
    """
    logger.info(f"Exporting annotation summary to {output_path}")
    
    # Cell type counts per sample
    summary = pd.crosstab(adata.obs['sample_id'], adata.obs['final_annotation'])
    summary['total_cells'] = summary.sum(axis=1)
    
    summary.to_csv(output_path, sep='\t')


# =============================================================================
# Main Pipeline Function
# =============================================================================

def run_qc_annotation_pipeline(
    cellranger_dirs,
    sample_ids,
    output_dir,
    qc_params=None,
    preprocessing_params=None,
    annotation_params=None,
    figures_dir=None,
):
    """
    Run the complete QC and annotation pipeline.
    
    Parameters
    ----------
    cellranger_dirs : list
        List of Cell Ranger output directories
    sample_ids : list
        List of sample identifiers
    output_dir : str
        Output directory for results
    qc_params : dict, optional
        QC parameters (uses defaults if not provided)
    preprocessing_params : dict, optional
        Preprocessing parameters
    annotation_params : dict, optional
        Annotation parameters
    figures_dir : str, optional
        Directory for figures (default: output_dir/figures)
        
    Returns
    -------
    AnnData
        Final annotated AnnData object
    """
    # Set defaults
    qc_params = {**DEFAULT_QC_PARAMS, **(qc_params or {})}
    preprocessing_params = {**DEFAULT_PREPROCESSING_PARAMS, **(preprocessing_params or {})}
    annotation_params = {**DEFAULT_ANNOTATION_PARAMS, **(annotation_params or {})}
    
    if figures_dir is None:
        figures_dir = os.path.join(output_dir, 'figures')
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    
    logger.info("="*60)
    logger.info("Starting QC and Annotation Pipeline")
    logger.info("="*60)
    
    # Step 1: Load data
    logger.info("\n[Step 1/6] Loading Cell Ranger outputs...")
    adata = load_cellranger_outputs(cellranger_dirs, sample_ids)
    
    # Save raw data
    adata.write(os.path.join(output_dir, 'adata_raw.h5ad'))
    
    # Step 2: Calculate QC metrics
    logger.info("\n[Step 2/6] Calculating QC metrics...")
    adata = calculate_qc_metrics(adata)
    
    # Generate pre-filter QC plots
    generate_qc_plots(adata, os.path.join(figures_dir, 'pre_filter'))
    
    # Step 3: Filter cells and genes
    logger.info("\n[Step 3/6] Filtering cells and genes...")
    adata = filter_cells_and_genes(adata, qc_params)
    
    # Step 4: Remove doublets
    if qc_params.get('doublet_removal', True):
        logger.info("\n[Step 4/6] Removing doublets...")
        adata = remove_doublets(adata, qc_params.get('doublet_rate', 0.06))
    else:
        logger.info("\n[Step 4/6] Skipping doublet removal...")
    
    # Generate post-filter QC plots
    generate_qc_plots(adata, os.path.join(figures_dir, 'post_filter'))
    
    # Export QC metrics
    export_qc_metrics(adata, os.path.join(output_dir, 'qc_metrics.tsv'))
    
    # Save QC-filtered data
    adata.write(os.path.join(output_dir, 'adata_qc.h5ad'))
    
    # Step 5: Preprocessing
    logger.info("\n[Step 5/6] Preprocessing data...")
    adata = preprocess_data(adata, preprocessing_params)
    
    # Save preprocessed data
    adata.write(os.path.join(output_dir, 'adata_pp.h5ad'))
    
    # Step 6: Annotation
    logger.info("\n[Step 6/6] Annotating cell types...")
    if annotation_params['method'] == 'popv':
        adata = annotate_with_popv(
            adata,
            annotation_params['popv_model'],
            annotation_params.get('batch_key', 'sample_id')
        )
    else:
        logger.info("Using leiden clusters as annotations")
        adata.obs['final_annotation'] = adata.obs['leiden'].astype(str)
    
    # Generate annotation plots
    generate_annotation_plots(adata, figures_dir)
    
    # Export annotation summary
    export_annotation_summary(adata, os.path.join(output_dir, 'annotation_summary.tsv'))
    
    # Save final annotated data
    adata.write(os.path.join(output_dir, 'adata_annotated.h5ad'))
    
    logger.info("\n" + "="*60)
    logger.info("Pipeline completed successfully!")
    logger.info("="*60)
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Final cells: {adata.n_obs}")
    logger.info(f"Final genes: {adata.n_vars}")
    logger.info(f"Cell types identified: {adata.obs['final_annotation'].nunique()}")
    
    return adata


# =============================================================================
# Snakemake Integration
# =============================================================================

def run_from_snakemake():
    """Run pipeline from Snakemake rule."""
    
    # Get inputs
    cellranger_dirs = snakemake.input.cellranger_dirs
    
    # Get outputs
    output_adata = snakemake.output.adata
    output_qc_metrics = snakemake.output.qc_metrics
    output_annotation_summary = snakemake.output.annotation_summary
    output_figures = snakemake.output.figures
    
    # Get params from config
    config = snakemake.config
    sample_ids = list(config.get('samples', {}).keys())
    
    qc_params = config.get('qc', {})
    preprocessing_params = config.get('preprocessing', {})
    annotation_params = config.get('annotation', {})
    
    output_dir = os.path.dirname(output_adata)
    
    # Set up logging to file
    log_file = snakemake.log[0] if snakemake.log else None
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Run pipeline
    adata = run_qc_annotation_pipeline(
        cellranger_dirs=cellranger_dirs,
        sample_ids=sample_ids,
        output_dir=output_dir,
        qc_params=qc_params,
        preprocessing_params=preprocessing_params,
        annotation_params=annotation_params,
        figures_dir=output_figures,
    )


# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    """Main function for standalone CLI usage."""
    
    parser = argparse.ArgumentParser(
        description='Single-cell QC and annotation pipeline'
    )
    parser.add_argument('--config', required=True, help='Path to config YAML file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--cellranger-dir', help='Base directory containing Cell Ranger outputs')
    parser.add_argument('--srascraper-dict', help='Path to SRAscraper dictionary pickle')
    
    args = parser.parse_args()
    
    # Load config
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    # Determine data source
    if args.srascraper_dict:
        # Load from SRAscraper dictionary
        with open(args.srascraper_dict, 'rb') as f:
            gse_dict = pickle.load(f)
        
        cellranger_base = args.cellranger_dir or config.get('output_dir', '.')
        adata = load_from_srascraper_dict(gse_dict, cellranger_base)
        
        # Then run the rest of the pipeline on the loaded data
        # ... (simplified for now)
        
    else:
        # Load from config sample information
        samples = config.get('samples', {})
        sample_ids = list(samples.keys())
        
        cellranger_base = args.cellranger_dir or config.get('output_dir', '.')
        cellranger_dirs = [
            os.path.join(cellranger_base, 'cellranger', sid)
            for sid in sample_ids
        ]
        
        adata = run_qc_annotation_pipeline(
            cellranger_dirs=cellranger_dirs,
            sample_ids=sample_ids,
            output_dir=args.output_dir,
            qc_params=config.get('qc', {}),
            preprocessing_params=config.get('preprocessing', {}),
            annotation_params=config.get('annotation', {}),
        )


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == '__main__':
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        main()
