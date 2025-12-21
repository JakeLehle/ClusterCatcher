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
    python scanpy_qc_annotation.py --cellranger-dir /path/to/cellranger --output-dir /path/to/output

Requirements:
    - scanpy
    - scrublet
    - popv (optional, for cell type annotation)
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
    'max_mito_pct': 20,      # Matches create_config.py default
    'min_counts': 500,
    'max_counts': 50000,
    'doublet_removal': True,
    'doublet_rate': 0.08,    # Matches create_config.py default
}

DEFAULT_PREPROCESSING_PARAMS = {
    'target_sum': 1e4,       # Normalization target
    'n_top_genes': 2000,     # Number of HVGs
    'n_pcs': 50,             # PCA components
    'n_neighbors': 15,       # For neighbor graph
    'leiden_resolution': 1.0,
}

DEFAULT_ANNOTATION_PARAMS = {
    'method': 'popv',        # Matches create_config.py default
    'popv_model': 'popv_immune_All_human_umap_from_cellxgene',
    'reference': None,       # Matches create_config.py field name
    'min_confidence': 0.5,
    'batch_key': 'sample_id',
}


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_cellranger_outputs(cellranger_base_dir, sample_ids):
    """
    Load Cell Ranger count matrices from multiple samples.
    
    Parameters
    ----------
    cellranger_base_dir : str
        Base directory containing Cell Ranger outputs (e.g., results/cellranger)
    sample_ids : list
        List of sample identifiers
        
    Returns
    -------
    AnnData
        Concatenated AnnData object with all samples
    """
    logger.info(f"Loading {len(sample_ids)} samples from {cellranger_base_dir}...")
    
    adatas = []
    for sample_id in sample_ids:
        # Construct path to Cell Ranger output
        cr_dir = os.path.join(cellranger_base_dir, sample_id)
        
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
            logger.warning(f"  Could not find matrix for {sample_id} at {cr_dir}, skipping...")
            continue
        
        # Make var names unique
        adata.var_names_make_unique()
        
        # Add sample metadata
        adata.obs['sample_id'] = sample_id
        adata.obs['original_barcode'] = adata.obs.index.tolist()
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
        
        # Store gene info if available
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
    if len(adatas) == 1:
        adata = adatas[0]
    else:
        adata = sc.concat(
            adatas, 
            join='outer', 
            label='sample_id', 
            keys=[a.obs['sample_id'].iloc[0] for a in adatas],
            index_unique=None
        )
    
    # Fix the sample_id column after concat
    adata.obs['sample_id'] = adata.obs['sample_id'].astype('category')
    
    logger.info(f"Total: {adata.n_obs} cells, {adata.n_vars} genes from {len(adatas)} samples")
    
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
        QC parameters (from config.yaml qc section)
        
    Returns
    -------
    AnnData
        Filtered AnnData
    """
    logger.info("Filtering cells and genes...")
    
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Get thresholds with defaults
    min_genes = params.get('min_genes', DEFAULT_QC_PARAMS['min_genes'])
    max_genes = params.get('max_genes', DEFAULT_QC_PARAMS['max_genes'])
    min_counts = params.get('min_counts', DEFAULT_QC_PARAMS['min_counts'])
    max_counts = params.get('max_counts', DEFAULT_QC_PARAMS['max_counts'])
    max_mito_pct = params.get('max_mito_pct', params.get('max_pct_mito', DEFAULT_QC_PARAMS['max_mito_pct']))
    min_cells = params.get('min_cells', DEFAULT_QC_PARAMS['min_cells'])
    
    logger.info(f"  Thresholds: min_genes={min_genes}, max_genes={max_genes}, "
                f"min_counts={min_counts}, max_mito={max_mito_pct}%")
    
    # Filter cells by minimum genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # Filter by counts
    if min_counts:
        adata = adata[adata.obs['total_counts'] >= min_counts, :].copy()
    if max_counts:
        adata = adata[adata.obs['total_counts'] <= max_counts, :].copy()
    
    # Filter by max genes
    adata = adata[adata.obs['n_genes_by_counts'] <= max_genes, :].copy()
    
    # Filter by mitochondrial content
    adata = adata[adata.obs['pct_counts_mt'] <= max_mito_pct, :].copy()
    
    # Filter genes by minimum cells
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars
    
    logger.info(f"  Cells: {n_cells_before} -> {n_cells_after} ({n_cells_before - n_cells_after} removed)")
    logger.info(f"  Genes: {n_genes_before} -> {n_genes_after} ({n_genes_before - n_genes_after} removed)")
    
    return adata


def remove_doublets(adata, expected_doublet_rate=0.08):
    """
    Detect and remove doublets using Scrublet.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData
    expected_doublet_rate : float
        Expected doublet rate (matches create_config.py default of 0.08)
        
    Returns
    -------
    AnnData
        AnnData with doublets removed
    """
    logger.info(f"Detecting doublets with Scrublet (expected rate: {expected_doublet_rate})...")
    
    try:
        import scrublet as scr
    except ImportError:
        logger.warning("Scrublet not installed, skipping doublet detection")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        return adata
    
    n_before = adata.n_obs
    
    # Run Scrublet per sample for better accuracy
    doublet_scores = []
    predicted_doublets = []
    
    samples = adata.obs['sample_id'].unique()
    
    for sample in samples:
        sample_mask = adata.obs['sample_id'] == sample
        sample_idx = np.where(sample_mask)[0]
        sample_adata = adata[sample_mask].copy()
        
        if sample_adata.n_obs < 100:
            logger.warning(f"  Sample {sample} has <100 cells, skipping doublet detection")
            doublet_scores.extend([0.0] * sample_adata.n_obs)
            predicted_doublets.extend([False] * sample_adata.n_obs)
            continue
        
        try:
            scrub = scr.Scrublet(sample_adata.X, expected_doublet_rate=expected_doublet_rate)
            scores, predictions = scrub.scrub_doublets(min_counts=2, min_cells=3, verbose=False)
            doublet_scores.extend(scores.tolist())
            predicted_doublets.extend(predictions.tolist())
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
    
    # Merge with defaults
    params = {**DEFAULT_PREPROCESSING_PARAMS, **params}
    
    # Store raw counts
    adata.layers['counts'] = adata.X.copy()
    
    # Normalize
    logger.info(f"  Normalizing (target_sum={params['target_sum']})...")
    sc.pp.normalize_total(adata, target_sum=params['target_sum'])
    sc.pp.log1p(adata)
    
    # Store normalized data
    adata.raw = adata.copy()
    
    # Find highly variable genes
    logger.info(f"  Finding {params['n_top_genes']} highly variable genes...")
    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=params['n_top_genes'],
            batch_key='sample_id' if 'sample_id' in adata.obs.columns and adata.obs['sample_id'].nunique() > 1 else None,
            flavor='seurat_v3',
            layer='counts'
        )
    except Exception as e:
        logger.warning(f"  seurat_v3 HVG failed: {e}, trying cell_ranger flavor...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=params['n_top_genes'],
            flavor='cell_ranger'
        )
    
    # Subset to HVGs for dimensionality reduction
    n_hvg = adata.var['highly_variable'].sum()
    logger.info(f"  Found {n_hvg} highly variable genes")
    
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    
    # Scale
    logger.info("  Scaling...")
    sc.pp.scale(adata_hvg, max_value=10)
    
    # PCA
    n_pcs = min(params['n_pcs'], adata_hvg.n_vars - 1, adata_hvg.n_obs - 1)
    logger.info(f"  Running PCA ({n_pcs} components)...")
    sc.tl.pca(adata_hvg, n_comps=n_pcs)
    
    # Copy PCA results back to full adata
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.uns['pca'] = adata_hvg.uns['pca']
    adata.varm['PCs'] = np.zeros((adata.n_vars, n_pcs))
    adata.varm['PCs'][adata.var['highly_variable'], :] = adata_hvg.varm['PCs']
    
    # Neighbors
    logger.info(f"  Computing neighbors (k={params['n_neighbors']})...")
    sc.pp.neighbors(adata, n_neighbors=params['n_neighbors'], n_pcs=n_pcs)
    
    # UMAP
    logger.info("  Computing UMAP...")
    sc.tl.umap(adata)
    
    # Leiden clustering
    logger.info(f"  Clustering (resolution={params['leiden_resolution']})...")
    sc.tl.leiden(adata, resolution=params['leiden_resolution'])
    
    n_clusters = adata.obs['leiden'].nunique()
    logger.info(f"  Found {n_clusters} clusters")
    
    return adata


# =============================================================================
# Annotation Functions
# =============================================================================

def annotate_with_popv(adata, model_name=None, reference=None, batch_key='sample_id'):
    """
    Annotate cell types using popV.
    
    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData
    model_name : str, optional
        Name of the popV model to use
    reference : str, optional
        Reference dataset or model path
    batch_key : str
        Key in obs for batch correction
        
    Returns
    -------
    AnnData
        AnnData with cell type annotations
    """
    # Use reference as model name if provided, otherwise use default
    if reference:
        model_name = reference
    elif model_name is None:
        model_name = DEFAULT_ANNOTATION_PARAMS['popv_model']
    
    logger.info(f"Running popV annotation with model: {model_name}")
    
    try:
        import popv
    except ImportError:
        logger.warning("popV not installed! Falling back to leiden clusters.")
        logger.warning("Install with: pip install popv")
        adata.obs['final_annotation'] = 'Cluster_' + adata.obs['leiden'].astype(str)
        return adata
    
    # Prepare data for popV
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
        adata.obs['final_annotation'] = 'Cluster_' + adata.obs['leiden'].astype(str)
    
    return adata


def annotate_cells(adata, params):
    """
    Annotate cells using the specified method.
    
    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData
    params : dict
        Annotation parameters from config
        
    Returns
    -------
    AnnData
        Annotated AnnData
    """
    # Merge with defaults
    params = {**DEFAULT_ANNOTATION_PARAMS, **params}
    
    method = params.get('method', 'popv')
    
    if method == 'popv':
        adata = annotate_with_popv(
            adata,
            model_name=params.get('popv_model'),
            reference=params.get('reference'),
            batch_key=params.get('batch_key', 'sample_id')
        )
    elif method == 'celltypist':
        logger.info("CellTypist annotation requested...")
        try:
            import celltypist
            from celltypist import models
            
            model_name = params.get('reference', 'Immune_All_Low.pkl')
            logger.info(f"  Using model: {model_name}")
            
            model = models.Model.load(model=model_name)
            predictions = celltypist.annotate(adata, model=model, majority_voting=True)
            adata.obs['celltypist_prediction'] = predictions.predicted_labels['majority_voting']
            adata.obs['final_annotation'] = adata.obs['celltypist_prediction']
            
        except ImportError:
            logger.warning("CellTypist not installed! Falling back to leiden clusters.")
            adata.obs['final_annotation'] = 'Cluster_' + adata.obs['leiden'].astype(str)
    else:
        logger.info(f"Unknown annotation method '{method}', using leiden clusters")
        adata.obs['final_annotation'] = 'Cluster_' + adata.obs['leiden'].astype(str)
    
    return adata


# =============================================================================
# Visualization Functions
# =============================================================================

def generate_qc_plots(adata, output_dir, prefix=''):
    """
    Generate QC visualization plots.
    
    Parameters
    ----------
    adata : AnnData
        AnnData with QC metrics
    output_dir : str
        Directory to save plots
    prefix : str
        Prefix for output files (e.g., 'pre_filter_' or 'post_filter_')
    """
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Generating QC plots in {output_dir}...")
    
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
    plt.savefig(os.path.join(output_dir, f'{prefix}qc_violin_plots.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f'{prefix}qc_violin_plots.png'), dpi=150, bbox_inches='tight')
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
    plt.savefig(os.path.join(output_dir, f'{prefix}qc_scatter_plots.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f'{prefix}qc_scatter_plots.png'), dpi=150, bbox_inches='tight')
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
    
    # UMAP with cell types
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(
        adata,
        color='final_annotation',
        legend_loc='on data',
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        title='Cell Type Annotations',
        ax=ax,
        show=False
    )
    plt.savefig(os.path.join(output_dir, 'umap_cell_types.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'umap_cell_types.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # UMAP with samples
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color='sample_id',
        frameon=False,
        title='Samples',
        ax=ax,
        show=False
    )
    plt.savefig(os.path.join(output_dir, 'umap_samples.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'umap_samples.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # UMAP with clusters
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color='leiden',
        legend_loc='on data',
        legend_fontsize=10,
        legend_fontoutline=2,
        frameon=False,
        title='Leiden Clusters',
        ax=ax,
        show=False
    )
    plt.savefig(os.path.join(output_dir, 'umap_leiden.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'umap_leiden.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Cell type proportions per sample
    try:
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
        plt.savefig(os.path.join(output_dir, 'cell_type_proportions.png'), dpi=150, bbox_inches='tight')
        plt.close()
    except Exception as e:
        logger.warning(f"Could not generate cell type proportions plot: {e}")


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
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
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
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    logger.info(f"Exporting annotation summary to {output_path}")
    
    # Cell type counts per sample
    summary = pd.crosstab(adata.obs['sample_id'], adata.obs['final_annotation'])
    summary['total_cells'] = summary.sum(axis=1)
    
    summary.to_csv(output_path, sep='\t')


# =============================================================================
# Main Pipeline Function
# =============================================================================

def run_qc_annotation_pipeline(
    cellranger_dir,
    sample_ids,
    output_adata_path,
    output_qc_path,
    output_annotation_path,
    figures_dir,
    qc_params=None,
    preprocessing_params=None,
    annotation_params=None,
):
    """
    Run the complete QC and annotation pipeline.
    
    Parameters
    ----------
    cellranger_dir : str
        Base directory containing Cell Ranger outputs
    sample_ids : list
        List of sample identifiers
    output_adata_path : str
        Path for output AnnData file
    output_qc_path : str
        Path for QC metrics TSV
    output_annotation_path : str
        Path for annotation summary TSV
    figures_dir : str
        Directory for figures
    qc_params : dict, optional
        QC parameters (uses defaults if not provided)
    preprocessing_params : dict, optional
        Preprocessing parameters
    annotation_params : dict, optional
        Annotation parameters
        
    Returns
    -------
    AnnData
        Final annotated AnnData object
    """
    # Set defaults
    qc_params = {**DEFAULT_QC_PARAMS, **(qc_params or {})}
    preprocessing_params = {**DEFAULT_PREPROCESSING_PARAMS, **(preprocessing_params or {})}
    annotation_params = {**DEFAULT_ANNOTATION_PARAMS, **(annotation_params or {})}
    
    # Create output directories
    os.makedirs(os.path.dirname(output_adata_path), exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    
    logger.info("="*60)
    logger.info("Starting QC and Annotation Pipeline")
    logger.info("="*60)
    logger.info(f"Cell Ranger directory: {cellranger_dir}")
    logger.info(f"Samples: {len(sample_ids)}")
    logger.info(f"Output: {output_adata_path}")
    
    # Step 1: Load data
    logger.info("\n[Step 1/6] Loading Cell Ranger outputs...")
    adata = load_cellranger_outputs(cellranger_dir, sample_ids)
    
    # Step 2: Calculate QC metrics
    logger.info("\n[Step 2/6] Calculating QC metrics...")
    adata = calculate_qc_metrics(adata)
    
    # Generate pre-filter QC plots
    generate_qc_plots(adata, figures_dir, prefix='pre_filter_')
    
    # Step 3: Filter cells and genes
    logger.info("\n[Step 3/6] Filtering cells and genes...")
    adata = filter_cells_and_genes(adata, qc_params)
    
    # Step 4: Remove doublets
    if qc_params.get('doublet_removal', True):
        logger.info("\n[Step 4/6] Removing doublets...")
        doublet_rate = qc_params.get('doublet_rate', DEFAULT_QC_PARAMS['doublet_rate'])
        adata = remove_doublets(adata, expected_doublet_rate=doublet_rate)
    else:
        logger.info("\n[Step 4/6] Skipping doublet removal...")
    
    # Generate post-filter QC plots
    generate_qc_plots(adata, figures_dir, prefix='post_filter_')
    
    # Export QC metrics
    export_qc_metrics(adata, output_qc_path)
    
    # Step 5: Preprocessing
    logger.info("\n[Step 5/6] Preprocessing data...")
    adata = preprocess_data(adata, preprocessing_params)
    
    # Step 6: Annotation
    logger.info("\n[Step 6/6] Annotating cell types...")
    adata = annotate_cells(adata, annotation_params)
    
    # Generate annotation plots
    generate_annotation_plots(adata, figures_dir)
    
    # Export annotation summary
    export_annotation_summary(adata, output_annotation_path)
    
    # Save final annotated data
    logger.info(f"\nSaving annotated data to {output_adata_path}...")
    adata.write(output_adata_path)
    
    logger.info("\n" + "="*60)
    logger.info("Pipeline completed successfully!")
    logger.info("="*60)
    logger.info(f"Final cells: {adata.n_obs}")
    logger.info(f"Final genes: {adata.n_vars}")
    logger.info(f"Cell types identified: {adata.obs['final_annotation'].nunique()}")
    
    return adata


# =============================================================================
# Snakemake Integration
# =============================================================================

def run_from_snakemake():
    """Run pipeline from Snakemake rule."""
    
    # Get outputs (these define the expected output paths)
    output_adata = snakemake.output.adata
    output_qc_metrics = snakemake.output.qc_metrics
    output_annotation_summary = snakemake.output.annotation_summary
    output_figures = snakemake.output.figures
    
    # Get params
    params = snakemake.params
    sample_ids = params.sample_ids
    cellranger_dir = params.cellranger_dir
    qc_params = params.qc_params
    preprocessing_params = params.preprocessing_params
    annotation_params = params.annotation_params
    
    # Set up logging to file
    log_file = snakemake.log[0] if snakemake.log else None
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Run pipeline
    adata = run_qc_annotation_pipeline(
        cellranger_dir=cellranger_dir,
        sample_ids=sample_ids,
        output_adata_path=output_adata,
        output_qc_path=output_qc_metrics,
        output_annotation_path=output_annotation_summary,
        figures_dir=output_figures,
        qc_params=qc_params,
        preprocessing_params=preprocessing_params,
        annotation_params=annotation_params,
    )


# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    """Main function for standalone CLI usage."""
    
    parser = argparse.ArgumentParser(
        description='Single-cell QC and annotation pipeline'
    )
    parser.add_argument('--cellranger-dir', required=True, 
                        help='Base directory containing Cell Ranger outputs')
    parser.add_argument('--sample-ids', nargs='+', required=True,
                        help='List of sample IDs')
    parser.add_argument('--output-dir', required=True, 
                        help='Output directory')
    parser.add_argument('--min-genes', type=int, default=200,
                        help='Minimum genes per cell')
    parser.add_argument('--min-counts', type=int, default=500,
                        help='Minimum counts per cell')
    parser.add_argument('--max-mito-pct', type=float, default=20,
                        help='Maximum mitochondrial percentage')
    parser.add_argument('--doublet-rate', type=float, default=0.08,
                        help='Expected doublet rate')
    parser.add_argument('--annotation-method', default='popv',
                        choices=['popv', 'celltypist', 'leiden'],
                        help='Cell type annotation method')
    
    args = parser.parse_args()
    
    # Build params from CLI
    qc_params = {
        'min_genes': args.min_genes,
        'min_counts': args.min_counts,
        'max_mito_pct': args.max_mito_pct,
        'doublet_rate': args.doublet_rate,
    }
    
    annotation_params = {
        'method': args.annotation_method,
    }
    
    # Define output paths
    output_adata = os.path.join(args.output_dir, 'annotation', 'adata_annotated.h5ad')
    output_qc = os.path.join(args.output_dir, 'qc', 'qc_metrics.tsv')
    output_annotation = os.path.join(args.output_dir, 'annotation', 'annotation_summary.tsv')
    figures_dir = os.path.join(args.output_dir, 'figures', 'qc')
    
    run_qc_annotation_pipeline(
        cellranger_dir=args.cellranger_dir,
        sample_ids=args.sample_ids,
        output_adata_path=output_adata,
        output_qc_path=output_qc,
        output_annotation_path=output_annotation,
        figures_dir=figures_dir,
        qc_params=qc_params,
        annotation_params=annotation_params,
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
