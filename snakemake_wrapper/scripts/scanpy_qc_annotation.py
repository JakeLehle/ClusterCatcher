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
   - Supports both fixed threshold and MAD-based adaptive filtering
3. Doublet detection and removal (Scrublet)
4. Automated cell type annotation using popV (HubModel from Hugging Face)
5. Post-annotation normalization and dimensionality reduction
6. Cluster-based annotation refinement
7. Visualization and QC reports

All parameters are read from the config file - no hardcoded defaults in this script.

The popV annotation workflow follows the structure:
1. Run popV on raw counts (no normalization)
2. Merge annotations back to original object
3. Normalize and run dimensionality reduction
4. Run BBKNN batch correction (optional)
5. Assign cluster-based final annotations using weighted scoring
6. Generate publication-quality plots

QC Filtering Modes:
- "fixed": Traditional threshold-based filtering (min/max genes, counts, mito%)
- "adaptive": MAD-based outlier detection (recommended for heterogeneous samples)

Usage:
    Called via Snakemake rule with snakemake.input/output/params
    
    Or standalone:
    python scanpy_qc_annotation.py --config config.yaml --output-dir /path/to/output

Requirements:
    - scanpy
    - scrublet
    - popv
    - bbknn (optional, for batch correction)
    - adjustText (for non-overlapping labels)
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
from textwrap import wrap

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

# Set matplotlib backend for headless operation (must be before pyplot import)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Patch
from matplotlib.colors import to_hex
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
# MAD-Based Outlier Detection Function
# =============================================================================

def is_outlier(adata, metric, nmads):
    """
    Detect outliers using Median Absolute Deviation (MAD).
    
    Identifies cells where the metric value is more than `nmads` MADs
    away from the median. This approach adapts to each dataset's distribution,
    making it more robust than fixed thresholds for heterogeneous samples.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    metric : str
        Column name in adata.obs to use for outlier detection.
        Common metrics: 'log1p_total_counts', 'log1p_n_genes_by_counts', 
        'pct_counts_mt', 'pct_counts_in_top_20_genes'
    nmads : float
        Number of MADs to use as threshold. Cells beyond this many MADs
        from the median are flagged as outliers.
        Typical values: 5 for counts/genes, 3 for mitochondrial
        
    Returns
    -------
    pd.Series
        Boolean series indicating outlier status (True = outlier)
        
    Notes
    -----
    The MAD (Median Absolute Deviation) is a robust measure of variability:
        MAD = median(|X_i - median(X)|)
    
    This function flags a cell as an outlier if:
        |X_i - median(X)| > nmads * MAD
    
    If MAD is 0 (can happen with very homogeneous data), falls back to
    standard deviation.
    
    References
    ----------
    This approach follows best practices from:
    - Luecken & Theis (2019) Current best practices in single-cell RNA-seq analysis
    - Heumos et al. (2023) Best practices for single-cell analysis across modalities
    """
    M = adata.obs[metric]
    median_val = M.median()
    mad = stats.median_abs_deviation(M, nan_policy='omit')
    
    # Handle edge case where MAD is 0 (very homogeneous data)
    if mad == 0:
        logger.warning(f"  MAD for {metric} is 0, using std instead")
        mad = M.std()
        if mad == 0:
            logger.warning(f"  Std for {metric} is also 0, no outliers detected")
            return pd.Series([False] * len(M), index=M.index)
    
    outlier_mask = (np.abs(M - median_val) > nmads * mad)
    
    n_outliers = outlier_mask.sum()
    logger.info(f"  {metric}: {n_outliers} outliers ({100*n_outliers/len(M):.1f}%) "
                f"[median={median_val:.2f}, MAD={mad:.2f}, threshold=±{nmads*mad:.2f}]")
    
    return outlier_mask


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_cellranger_outputs(cellranger_base_dir, sample_ids):
    """Load Cell Ranger count matrices from multiple samples."""
    logger.info(f"Loading {len(sample_ids)} samples from {cellranger_base_dir}...")
    
    adatas = []
    for sample_id in sample_ids:
        cr_dir = os.path.join(cellranger_base_dir, sample_id)
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
        
        adata.var_names_make_unique()
        adata.obs['sample_id'] = sample_id
        adata.obs['original_barcode'] = adata.obs.index.tolist()
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
        
        if 'gene_ids' in adata.var.columns:
            adata.var['gene_ids'] = adata.var['gene_ids']
        if 'feature_types' in adata.var.columns:
            adata.var['feature_types'] = adata.var['feature_types']
        
        adata.var['gene_symbol'] = adata.var_names
        adata.var['feature_name'] = adata.var_names
        
        logger.info(f"    {sample_id}: {adata.n_obs} cells, {adata.n_vars} genes")
        adatas.append(adata)
    
    if len(adatas) == 0:
        raise ValueError("No samples were loaded successfully!")
    
    logger.info("Concatenating samples...")
    if len(adatas) == 1:
        adata = adatas[0]
    else:
        adata = sc.concat(adatas, join='outer', label='sample_id', 
                         keys=[a.obs['sample_id'].iloc[0] for a in adatas], index_unique=None)
    
    adata.obs['sample_id'] = adata.obs['sample_id'].astype('category')
    logger.info(f"Total: {adata.n_obs} cells, {adata.n_vars} genes from {len(adatas)} samples")
    
    return adata


# =============================================================================
# QC Functions
# =============================================================================

def calculate_qc_metrics(adata, qc_params=None):
    """
    Calculate QC metrics for all cells.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with raw counts
    qc_params : dict, optional
        QC parameters from config. If provided and qc_mode is 'adaptive',
        calculates additional metrics needed for MAD-based filtering.
        
    Returns
    -------
    AnnData
        Updated AnnData with QC metrics in .obs and gene annotations in .var
        
    Notes
    -----
    For adaptive mode, calculates:
    - log1p_total_counts: log1p-transformed UMI counts (for MAD filtering)
    - log1p_n_genes_by_counts: log1p-transformed gene counts (for MAD filtering)
    - pct_counts_in_top_20_genes: Concentration metric (for MAD filtering)
    
    For fixed mode (or if qc_params not provided):
    - Standard metrics without log1p transformation
    """
    logger.info("Calculating QC metrics...")
    
    # Identify gene categories
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL', 'Rps', 'Rpl'))
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]', case=False, regex=True)
    
    # Check if adaptive mode needs percent_top and log1p
    qc_mode = qc_params.get('qc_mode', 'fixed') if qc_params else 'fixed'
    
    if qc_mode == 'adaptive':
        # For adaptive mode, calculate with log1p and percent_top for MAD-based filtering
        logger.info("  Using adaptive mode metrics (log1p=True, percent_top=[20])")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mt', 'ribo', 'hb'], 
            percent_top=[20], log1p=True, inplace=True
        )
    else:
        # Original behavior for fixed mode
        logger.info("  Using fixed mode metrics (log1p=False)")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mt', 'ribo', 'hb'], 
            percent_top=None, log1p=False, inplace=True
        )
    
    logger.info(f"  Mean genes/cell: {adata.obs['n_genes_by_counts'].mean():.1f}")
    logger.info(f"  Mean counts/cell: {adata.obs['total_counts'].mean():.1f}")
    logger.info(f"  Mean % mito: {adata.obs['pct_counts_mt'].mean():.1f}%")
    
    return adata


def filter_cells_and_genes(adata, qc_params):
    """
    Filter cells and genes based on QC thresholds or MAD-based outlier detection.
    
    Supports two filtering modes controlled by qc_params['qc_mode']:
    
    1. "fixed" (default): Traditional threshold-based filtering
       - Uses min_genes, max_genes, min_counts, max_counts, max_mito_pct
       - Good for well-characterized datasets with known quality ranges
    
    2. "adaptive": MAD-based outlier detection
       - Uses median absolute deviation to identify outliers
       - Adapts to each dataset's distribution
       - Recommended for heterogeneous samples or when optimal thresholds are unknown
       - Based on best practices from Luecken & Theis (2019)
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with QC metrics calculated
    qc_params : dict
        QC parameters from config file. Expected keys depend on mode:
        
        Common (both modes):
        - qc_mode: "fixed" or "adaptive" (default: "fixed")
        - min_genes: Minimum genes per cell (default: 200)
        - max_mito_pct: Maximum mitochondrial % (default: 20)
        - min_cells: Minimum cells expressing a gene (default: 20)
        
        Fixed mode specific:
        - max_genes: Maximum genes per cell (default: None)
        - min_counts: Minimum UMI counts (default: None)
        - max_counts: Maximum UMI counts (default: None)
        
        Adaptive mode specific:
        - mad_n_counts: MADs for total counts outliers (default: 5)
        - mad_n_genes: MADs for gene count outliers (default: 5)
        - mad_top_genes: MADs for top-20 gene concentration outliers (default: 5)
        - mad_mito: MADs for mitochondrial outliers (default: 3)
        
    Returns
    -------
    AnnData
        Filtered AnnData object
    """
    logger.info("Filtering cells and genes...")
    
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Get common parameters with defaults
    qc_mode = qc_params.get('qc_mode', 'fixed')
    min_genes = qc_params.get('min_genes', 200)
    max_mito_pct = qc_params.get('max_mito_pct', 20)
    min_cells = qc_params.get('min_cells', 20)
    
    logger.info(f"  QC mode: {qc_mode}")
    
    if qc_mode == 'adaptive':
        # =====================================================================
        # ADAPTIVE MODE (MAD-based) - matches Jake's interactive workflow
        # =====================================================================
        mad_n_counts = qc_params.get('mad_n_counts', 5)
        mad_n_genes = qc_params.get('mad_n_genes', 5)
        mad_top_genes = qc_params.get('mad_top_genes', 5)
        mad_mito = qc_params.get('mad_mito', 3)
        
        logger.info(f"  MAD thresholds: counts={mad_n_counts}, genes={mad_n_genes}, "
                    f"top_genes={mad_top_genes}, mito={mad_mito}")
        logger.info(f"  Hard caps: min_genes={min_genes}, max_mito_pct={max_mito_pct}%")
        
        # Verify required columns exist
        required_cols = ['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_in_top_20_genes']
        missing_cols = [col for col in required_cols if col not in adata.obs.columns]
        if missing_cols:
            raise ValueError(
                f"Adaptive QC mode requires columns {missing_cols} which are missing. "
                "Ensure calculate_qc_metrics() was called with qc_params containing qc_mode='adaptive'."
            )
        
        # Detect outliers using MAD on log1p-transformed metrics
        logger.info("  Detecting outliers...")
        adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", mad_n_counts)
            | is_outlier(adata, "log1p_n_genes_by_counts", mad_n_genes)
            | is_outlier(adata, "pct_counts_in_top_20_genes", mad_top_genes)
        )
        
        # Mitochondrial outliers: MAD-based OR hard cap
        # This follows the approach: cells are outliers if they're statistical outliers
        # OR if they exceed a hard biological threshold
        adata.obs["mt_outlier"] = (
            is_outlier(adata, "pct_counts_mt", mad_mito) 
            | (adata.obs["pct_counts_mt"] > max_mito_pct)
        )
        
        n_outliers = adata.obs["outlier"].sum()
        n_mt_outliers = adata.obs["mt_outlier"].sum()
        n_both = (adata.obs["outlier"] & adata.obs["mt_outlier"]).sum()
        logger.info(f"  General outliers: {n_outliers} cells")
        logger.info(f"  Mito outliers: {n_mt_outliers} cells")
        logger.info(f"  Both: {n_both} cells")
        
        # Filter based on outlier flags
        adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
        
        # Apply minimum gene threshold (always applied as a sanity check)
        sc.pp.filter_cells(adata, min_genes=min_genes)
        
    else:
        # =====================================================================
        # FIXED MODE (threshold-based) - original behavior
        # =====================================================================
        max_genes = qc_params.get('max_genes')
        min_counts = qc_params.get('min_counts')
        max_counts = qc_params.get('max_counts')
        
        logger.info(f"  Thresholds: min_genes={min_genes}, max_genes={max_genes}, "
                    f"min_counts={min_counts}, max_counts={max_counts}, "
                    f"max_mito={max_mito_pct}%")
        
        # Apply filters sequentially
        sc.pp.filter_cells(adata, min_genes=min_genes)
        
        if min_counts is not None:
            adata = adata[adata.obs['total_counts'] >= min_counts, :].copy()
        if max_counts is not None:
            adata = adata[adata.obs['total_counts'] <= max_counts, :].copy()
        if max_genes is not None:
            adata = adata[adata.obs['n_genes_by_counts'] <= max_genes, :].copy()
        if max_mito_pct is not None:
            adata = adata[adata.obs['pct_counts_mt'] <= max_mito_pct, :].copy()
    
    # Filter genes by min_cells (applied in both modes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    logger.info(f"  Cells: {n_cells_before} -> {adata.n_obs} ({n_cells_before - adata.n_obs} removed)")
    logger.info(f"  Genes: {n_genes_before} -> {adata.n_vars} ({n_genes_before - adata.n_vars} removed)")
    
    return adata


def remove_doublets(adata, qc_params):
    """Detect and remove doublets using Scrublet."""
    doublet_removal = qc_params.get('doublet_removal', True)
    expected_doublet_rate = qc_params.get('doublet_rate', 0.08)
    
    if not doublet_removal:
        logger.info("Doublet removal disabled, skipping...")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        return adata
    
    logger.info(f"Detecting doublets with Scrublet (expected rate: {expected_doublet_rate})...")
    
    try:
        import scrublet as scr
    except ImportError:
        logger.warning("Scrublet not installed, skipping doublet detection")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        return adata
    
    n_before = adata.n_obs
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
            doublet_scores.extend(scores.tolist())
            predicted_doublets.extend(predictions.tolist())
            n_doublets = sum(predictions)
            logger.info(f"  {sample}: {n_doublets} doublets ({100*n_doublets/len(predictions):.1f}%)")
        except Exception as e:
            logger.warning(f"  Scrublet failed for {sample}: {e}")
            doublet_scores.extend([0.0] * sample_adata.n_obs)
            predicted_doublets.extend([False] * sample_adata.n_obs)
    
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets
    adata = adata[~adata.obs['predicted_doublet'], :].copy()
    
    logger.info(f"  Total doublets removed: {n_before - adata.n_obs}")
    return adata


# =============================================================================
# popV Annotation Functions
# =============================================================================

def run_popv_annotation(adata, annotation_params):
    """
    Run popV annotation using HubModel from Hugging Face.
    popV needs raw counts (no normalization) in adata.X.
    """
    huggingface_repo = annotation_params['popv_huggingface_repo']
    query_batch_key = annotation_params['batch_key']
    prediction_mode = annotation_params['popv_prediction_mode']
    gene_symbol_key = annotation_params['popv_gene_symbol_key']
    cache_dir = annotation_params.get('popv_cache_dir', 'tmp/popv_models')
    
    logger.info("="*60)
    logger.info("Running popV Annotation")
    logger.info("="*60)
    logger.info(f"  Hugging Face repo: {huggingface_repo}")
    logger.info(f"  Prediction mode: {prediction_mode}")
    logger.info(f"  Batch key: {query_batch_key}")
    logger.info(f"  Gene symbol key: {gene_symbol_key}")
    
    try:
        import popv
    except ImportError:
        raise ImportError("popV is required. Install with: pip install popv")
    
    if gene_symbol_key not in adata.var.columns:
        adata.var[gene_symbol_key] = adata.var_names
    
    if query_batch_key not in adata.obs.columns:
        if 'sample_id' in adata.obs.columns:
            query_batch_key = 'sample_id'
            logger.info(f"  Using 'sample_id' as batch key")
        else:
            raise ValueError("No valid batch key found")
    
    original_obs_columns = list(adata.obs.columns)
    
    logger.info("  Pulling model from Hugging Face...")
    os.makedirs(cache_dir, exist_ok=True)
    hmo = popv.hub.HubModel.pull_from_huggingface_hub(huggingface_repo, cache_dir=cache_dir)
    
    logger.info("  Running popV annotation (this may take a while)...")
    adata_annotated = hmo.annotate_data(
        adata,
        query_batch_key=query_batch_key,
        prediction_mode=prediction_mode,
        gene_symbols=gene_symbol_key
    )
    
    # Filter to query cells only
    if '_dataset' in adata_annotated.obs.columns:
        logger.info("  Filtering to query cells only...")
        adata_annotated = adata_annotated[adata_annotated.obs["_dataset"] == "query"].copy()
    
    # Merge annotations back
    logger.info("  Merging popV annotations...")
    popv_obs = adata_annotated.obs.drop(columns=original_obs_columns, errors='ignore')
    merged_df = pd.merge(adata.obs, popv_obs, left_index=True, right_index=True, how='inner')
    
    adata = adata[adata.obs_names.isin(merged_df.index)].copy()
    merged_df = merged_df.loc[adata.obs_names]
    adata.obs = merged_df
    
    assert all(adata.obs_names == merged_df.index), "Indices do not match after merge!"
    
    logger.info(f"  Merged {len(popv_obs.columns)} popV columns, {adata.n_obs} cells")
    
    if 'popv_prediction' in adata.obs.columns:
        counts = adata.obs['popv_prediction'].value_counts()
        logger.info("  Cell type distribution (top 10):")
        for ct, count in counts.head(10).items():
            logger.info(f"    {ct}: {count} ({100*count/adata.n_obs:.1f}%)")
    
    return adata


def post_annotation_processing(adata, preprocessing_params):
    """Process data after popV annotation: normalize, UMAP, cluster."""
    logger.info("Post-annotation processing...")
    
    target_sum = preprocessing_params['target_sum']
    n_pcs = preprocessing_params['n_pcs']
    n_neighbors = preprocessing_params['n_neighbors']
    leiden_resolution = preprocessing_params['leiden_resolution']
    run_bbknn = preprocessing_params.get('run_bbknn', False)
    bbknn_batch_key = preprocessing_params.get('bbknn_batch_key', 'sample_id')
    
    adata.layers['counts'] = adata.X.copy()
    
    logger.info(f"  Normalizing (target_sum={target_sum})...")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    
    if n_pcs is None:
        from os import cpu_count
        n_pcs = cpu_count()
        logger.info(f"  Using CPU count for n_pcs: {n_pcs}")
    
    n_pcs_actual = min(n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    
    logger.info(f"  Computing neighbors (n_pcs={n_pcs_actual})...")
    sc.pp.neighbors(adata, n_pcs=n_pcs_actual)
    
    logger.info("  Computing UMAP...")
    sc.tl.umap(adata)
    
    logger.info(f"  Running Leiden clustering (resolution={leiden_resolution})...")
    sc.tl.leiden(adata, key_added='clusters', resolution=leiden_resolution, random_state=42)
    logger.info(f"  Found {adata.obs['clusters'].nunique()} clusters")
    
    if run_bbknn:
        logger.info(f"  Running BBKNN batch correction...")
        try:
            sc.external.pp.bbknn(adata, batch_key=bbknn_batch_key, n_pcs=n_pcs_actual)
            sc.tl.umap(adata)
            logger.info("  BBKNN complete, UMAP recomputed")
        except Exception as e:
            logger.warning(f"  BBKNN failed: {e}")
    
    return adata


def assign_cluster_based_annotation(adata):
    """Assign final annotations using cluster-level weighted scoring."""
    logger.info("Assigning cluster-based annotations...")
    
    for col in ['clusters', 'popv_prediction']:
        if col not in adata.obs.columns:
            raise ValueError(f"Required column '{col}' not found")
    
    has_scores = 'popv_prediction_score' in adata.obs.columns
    
    if has_scores:
        df = adata.obs[['clusters', 'popv_prediction', 'popv_prediction_score']].copy()
        df['popv_prediction_score'] = pd.to_numeric(df['popv_prediction_score'], errors='coerce')
        
        def min_max_normalize(x):
            x_min, x_range = x.min(), x.max() - x.min()
            return np.zeros(len(x)) if x_range == 0 else (x - x_min) / x_range
        
        def linear_weights(x):
            total = x.sum()
            return np.ones(len(x)) / len(x) if total == 0 else x / total
        
        df['normalized_score'] = df.groupby('clusters')['popv_prediction_score'].transform(min_max_normalize)
        df['weight'] = df.groupby('clusters')['normalized_score'].transform(linear_weights)
        df['weighted_score'] = df['popv_prediction_score'] * df['weight']
        
        cluster_type_scores = df.groupby(['clusters', 'popv_prediction'])['weighted_score'].sum().reset_index()
        dominant_types = cluster_type_scores.loc[
            cluster_type_scores.groupby('clusters')['weighted_score'].idxmax()
        ].set_index('clusters')['popv_prediction']
    else:
        logger.warning("  No popv_prediction_score found, using majority voting")
        cluster_type_counts = adata.obs.groupby(['clusters', 'popv_prediction']).size().reset_index(name='count')
        dominant_types = cluster_type_counts.loc[
            cluster_type_counts.groupby('clusters')['count'].idxmax()
        ].set_index('clusters')['popv_prediction']
    
    adata.obs['final_annotation'] = adata.obs['clusters'].map(dominant_types)
    
    final_counts = adata.obs['final_annotation'].value_counts()
    logger.info("  Final annotation distribution:")
    for ct, count in final_counts.items():
        logger.info(f"    {ct}: {count} ({100*count/adata.n_obs:.1f}%)")
    
    return adata


# =============================================================================
# Visualization Functions
# =============================================================================

def generate_qc_plots(adata, output_dir, prefix=''):
    """Generate QC visualization plots."""
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Generating QC plots in {output_dir}...")
    
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
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', ax=axes[0], show=False)
    axes[0].set_title('Counts vs Genes (colored by % mito)')
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color='n_genes_by_counts', ax=axes[1], show=False)
    axes[1].set_title('Counts vs % Mito (colored by genes)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{prefix}qc_scatter_plots.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f'{prefix}qc_scatter_plots.png'), dpi=150, bbox_inches='tight')
    plt.close()


def gen_mpl_labels(adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None, color_by_group=False):
    """Generate non-overlapping labels for UMAP plots."""
    try:
        from adjustText import adjust_text
    except ImportError:
        logger.warning("adjustText not installed, skipping label adjustment")
        return
    
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}
    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)

    text_colors = {group: 'black' for group in adata.obs[groupby].cat.categories}
    if color_by_group and groupby + "_colors" in adata.uns:
        for i, group in enumerate(adata.obs[groupby].cat.categories):
            if group not in exclude:
                text_colors[group] = adata.uns[groupby + "_colors"][i]

    if ax is None:
        texts = [plt.text(x=x, y=y, s=k, color=text_colors.get(k, 'black'), **text_kwargs) 
                 for k, (x, y) in medians.items()]
    else:
        texts = [ax.text(x=x, y=y, s=k, color=text_colors.get(k, 'black'), **text_kwargs) 
                 for k, (x, y) in medians.items()]
    adjust_text(texts, **adjust_kwargs)


def plot_umap_popv_prediction(adata, output_dir):
    """Generate UMAP plot colored by popV prediction."""
    os.makedirs(output_dir, exist_ok=True)
    
    if 'popv_prediction' not in adata.obs.columns:
        logger.warning("popv_prediction not found, skipping")
        return
    
    logger.info("  Generating popV prediction UMAP...")
    
    cmap = plt.get_cmap('turbo')
    value_cat = pd.Categorical(adata.obs['popv_prediction'])
    values = np.linspace(0, 1, len(value_cat.categories))
    palette = [cmap(value) for value in values]
    
    with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": 150, "figure.frameon": False}):
        fig = sc.pl.umap(adata, color=['popv_prediction'], palette=palette, 
                        frameon=False, size=5, return_fig=True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'UMAP_popv_prediction.pdf'), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, 'UMAP_popv_prediction.png'), dpi=150, bbox_inches='tight')
        plt.close()


def plot_umap_popv_score(adata, output_dir):
    """Generate UMAP plot colored by popV prediction score."""
    os.makedirs(output_dir, exist_ok=True)
    
    if 'popv_prediction_score' not in adata.obs.columns:
        logger.warning("popv_prediction_score not found, skipping")
        return
    
    logger.info("  Generating popV prediction score UMAP...")
    
    with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": 150, "figure.frameon": False}):
        fig = sc.pl.umap(adata, color=['popv_prediction_score'], cmap='magma',
                        frameon=False, size=5, return_fig=True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'UMAP_popv_prediction_score.pdf'), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, 'UMAP_popv_prediction_score.png'), dpi=150, bbox_inches='tight')
        plt.close()


def plot_umap_final_annotation(adata, output_dir):
    """Generate UMAP plot with non-overlapping labels for final annotation."""
    os.makedirs(output_dir, exist_ok=True)
    
    if 'final_annotation' not in adata.obs.columns:
        logger.warning("final_annotation not found, skipping")
        return
    
    logger.info("  Generating final annotation UMAP with labels...")
    
    adata.obs['final_annotation'] = adata.obs['final_annotation'].astype('category')
    adata.obs['final_annotation'] = adata.obs['final_annotation'].cat.remove_unused_categories()
    
    cmap = plt.get_cmap('turbo')
    value_cat = pd.Categorical(adata.obs['final_annotation'])
    values = np.linspace(0, 1, len(value_cat.categories))
    palette = [cmap(value) for value in values]
    
    combined_effects = [
        pe.withStroke(linewidth=6, foreground="white"),
        pe.withStroke(linewidth=1, foreground="black"),
        pe.Normal()
    ]
    
    with plt.rc_context({"figure.figsize": (12, 12), "figure.dpi": 150, "figure.frameon": False}):
        ax = sc.pl.umap(adata, color='final_annotation', show=False, legend_loc=None, 
                       frameon=False, size=5, palette=palette)
        gen_mpl_labels(
            adata, 'final_annotation', exclude=("None", "Unknown"), ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
            text_kwargs=dict(fontsize=14, path_effects=combined_effects),
            color_by_group=True
        )
        fig = ax.get_figure()
        fig.tight_layout()
        plt.savefig(os.path.join(output_dir, 'UMAP_final_annotation.pdf'), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, 'UMAP_final_annotation.png'), dpi=150, bbox_inches='tight')
        plt.close()


def plot_umap_samples_clusters(adata, output_dir):
    """Generate UMAP plots colored by samples and clusters."""
    os.makedirs(output_dir, exist_ok=True)
    logger.info("  Generating sample and cluster UMAPs...")
    
    with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": 150, "figure.frameon": False}):
        fig = sc.pl.umap(adata, color=['sample_id'], frameon=False, size=5, return_fig=True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'UMAP_samples.pdf'), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, 'UMAP_samples.png'), dpi=150, bbox_inches='tight')
        plt.close()
    
    with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": 150, "figure.frameon": False}):
        fig = sc.pl.umap(adata, color=['clusters'], frameon=False, size=5, return_fig=True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'UMAP_clusters.pdf'), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, 'UMAP_clusters.png'), dpi=150, bbox_inches='tight')
        plt.close()


def get_dominant_clusters(adata, min_cells=50):
    """Get clusters with at least min_cells cells, sorted by size."""
    cluster_counts = adata.obs['clusters'].value_counts()
    dominant = cluster_counts[cluster_counts >= min_cells].index.tolist()
    return sorted(dominant, key=lambda x: int(x) if x.isdigit() else x)


def plot_stacked_bar_cluster_composition(adata, output_dir):
    """Generate stacked bar plot showing cell type composition per cluster."""
    os.makedirs(output_dir, exist_ok=True)
    
    if 'clusters' not in adata.obs.columns or 'popv_prediction' not in adata.obs.columns:
        logger.warning("Missing required columns for stacked bar plot")
        return
    
    logger.info("  Generating stacked bar plot of cluster composition...")
    
    CLUSTERS = get_dominant_clusters(adata)
    MIN_PERCENT = 5
    BAR_WIDTH = 0.6
    PLOT_SPACING = 1.0
    
    fig, ax = plt.subplots(figsize=(max(24, len(CLUSTERS) * 2), 8))
    
    all_celltypes = adata.obs['popv_prediction'].astype('category').cat.categories
    cmap = plt.get_cmap('turbo')
    values = np.linspace(0, 1, len(all_celltypes))
    color_map = dict(zip(all_celltypes, [to_hex(cmap(value)) for value in values]))
    
    for cluster_idx, cluster_num in enumerate(CLUSTERS):
        sub_clust_cells = adata.obs[adata.obs['clusters'] == cluster_num]
        celltype_counts = sub_clust_cells['popv_prediction'].value_counts()
        percentages = (celltype_counts / celltype_counts.sum() * 100).round(2)
        
        nonzero_mask = percentages > 0
        filtered_labels = percentages.index[nonzero_mask].tolist()
        filtered_percentages = percentages[nonzero_mask].tolist()
        filtered_colors = [color_map.get(label, '#888888') for label in filtered_labels]
        
        x_pos = cluster_idx * PLOT_SPACING
        
        bottom = 0
        for percent, color in zip(filtered_percentages, filtered_colors):
            alpha = 0.4 if percent < MIN_PERCENT else 1.0
            ax.bar(x_pos, percent, width=BAR_WIDTH, bottom=bottom, 
                   color=color, edgecolor='white', alpha=alpha)
            bottom += percent
    
    ax.set_xticks([i * PLOT_SPACING for i in range(len(CLUSTERS))])
    ax.set_xticklabels([f"Cluster {c}" for c in CLUSTERS], rotation=45, ha='right')
    ax.set_ylabel('Percentage')
    ax.set_title('Cell Type Composition per Cluster')
    ax.set_ylim(0, 100)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'Stacked_Bar_Cluster_Composition.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'Stacked_Bar_Cluster_Composition.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_cell_type_proportions(adata, output_dir):
    """Generate cell type proportions per sample."""
    os.makedirs(output_dir, exist_ok=True)
    
    if 'final_annotation' not in adata.obs.columns:
        logger.warning("final_annotation not found, skipping proportions plot")
        return
    
    logger.info("  Generating cell type proportions plot...")
    
    proportions = pd.crosstab(adata.obs['sample_id'], adata.obs['final_annotation'], normalize='index') * 100
    
    fig, ax = plt.subplots(figsize=(12, 8))
    proportions.plot(kind='bar', stacked=True, ax=ax, colormap='turbo')
    ax.set_ylabel('Percentage')
    ax.set_xlabel('Sample')
    ax.set_title('Cell Type Proportions per Sample')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cell_type_proportions.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'cell_type_proportions.png'), dpi=150, bbox_inches='tight')
    plt.close()


def generate_annotation_plots(adata, output_dir):
    """Generate all annotation-related plots."""
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Generating annotation plots in {output_dir}...")
    
    if 'X_umap' not in adata.obsm:
        logger.warning("No UMAP found, skipping UMAP-based plots")
        plot_cell_type_proportions(adata, output_dir)
        return
    
    plot_umap_popv_prediction(adata, output_dir)
    plot_umap_popv_score(adata, output_dir)
    plot_umap_final_annotation(adata, output_dir)
    plot_umap_samples_clusters(adata, output_dir)
    plot_stacked_bar_cluster_composition(adata, output_dir)
    plot_cell_type_proportions(adata, output_dir)


def export_qc_metrics(adata, output_path):
    """Export QC metrics summary to TSV."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    logger.info(f"Exporting QC metrics to {output_path}")
    
    qc_summary = adata.obs.groupby('sample_id').agg({
        'n_genes_by_counts': ['count', 'mean', 'median', 'std'],
        'total_counts': ['mean', 'median', 'std'],
        'pct_counts_mt': ['mean', 'median', 'std'],
    })
    qc_summary.columns = ['_'.join(col).strip() for col in qc_summary.columns.values]
    qc_summary = qc_summary.rename(columns={'n_genes_by_counts_count': 'n_cells'})
    qc_summary.to_csv(output_path, sep='\t')


def export_annotation_summary(adata, output_path):
    """Export annotation summary to TSV."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    logger.info(f"Exporting annotation summary to {output_path}")
    
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
    qc_params,
    preprocessing_params,
    annotation_params,
):
    """
    Run the complete QC and annotation pipeline.
    
    Workflow:
    1. Load Cell Ranger outputs
    2. Calculate QC metrics
    3. Filter cells and genes
    4. Remove doublets
    5. Run popV annotation (on raw counts)
    6. Post-annotation processing (normalize, UMAP, cluster)
    7. Assign cluster-based final annotations
    8. Generate plots
    9. Export summaries
    """
    os.makedirs(os.path.dirname(output_adata_path), exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    
    logger.info("="*60)
    logger.info("ClusterCatcher QC and Annotation Pipeline")
    logger.info("="*60)
    logger.info(f"Cell Ranger directory: {cellranger_dir}")
    logger.info(f"Samples: {len(sample_ids)}")
    logger.info(f"Output: {output_adata_path}")
    
    logger.info("\nQC Parameters:")
    for k, v in qc_params.items():
        logger.info(f"  {k}: {v}")
    logger.info("\nPreprocessing Parameters:")
    for k, v in preprocessing_params.items():
        logger.info(f"  {k}: {v}")
    logger.info("\nAnnotation Parameters:")
    for k, v in annotation_params.items():
        logger.info(f"  {k}: {v}")
    
    # Step 1: Load data
    logger.info("\n" + "="*60)
    logger.info("[Step 1/7] Loading Cell Ranger outputs...")
    logger.info("="*60)
    adata = load_cellranger_outputs(cellranger_dir, sample_ids)
    
    # Step 2: Calculate QC metrics - pass qc_params for mode detection
    logger.info("\n" + "="*60)
    logger.info("[Step 2/7] Calculating QC metrics...")
    logger.info("="*60)
    adata = calculate_qc_metrics(adata, qc_params)
    generate_qc_plots(adata, figures_dir, prefix='pre_filter_')
    
    # Step 3: Filter cells and genes
    logger.info("\n" + "="*60)
    logger.info("[Step 3/7] Filtering cells and genes...")
    logger.info("="*60)
    adata = filter_cells_and_genes(adata, qc_params)
    
    # Step 4: Remove doublets
    logger.info("\n" + "="*60)
    logger.info("[Step 4/7] Doublet detection...")
    logger.info("="*60)
    adata = remove_doublets(adata, qc_params)
    generate_qc_plots(adata, figures_dir, prefix='post_filter_')
    export_qc_metrics(adata, output_qc_path)
    
    # Step 5: popV Annotation (on raw counts - NO normalization yet)
    logger.info("\n" + "="*60)
    logger.info("[Step 5/7] Running popV annotation...")
    logger.info("="*60)
    adata = run_popv_annotation(adata, annotation_params)
    
    # Step 6: Post-annotation processing
    logger.info("\n" + "="*60)
    logger.info("[Step 6/7] Post-annotation processing...")
    logger.info("="*60)
    adata = post_annotation_processing(adata, preprocessing_params)
    
    # Step 7: Assign cluster-based final annotations
    logger.info("\n" + "="*60)
    logger.info("[Step 7/7] Assigning cluster-based annotations...")
    logger.info("="*60)
    adata = assign_cluster_based_annotation(adata)
    
    # Generate annotation plots
    generate_annotation_plots(adata, figures_dir)
    export_annotation_summary(adata, output_annotation_path)
    
    # Save final annotated data
    logger.info(f"\nSaving annotated data to {output_adata_path}...")
    adata.var.index = adata.var.index.astype(str)
    if adata.raw is not None:
        adata.raw.var.index = adata.raw.var.index.astype(str)
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].fillna('').astype(str)
    adata.write(output_adata_path)
    
    logger.info("\n" + "="*60)
    logger.info("Pipeline completed successfully!")
    logger.info("="*60)
    logger.info(f"Final cells: {adata.n_obs}")
    logger.info(f"Final genes: {adata.n_vars}")
    logger.info(f"Clusters: {adata.obs['clusters'].nunique()}")
    logger.info(f"Cell types (final): {adata.obs['final_annotation'].nunique()}")
    
    return adata


# =============================================================================
# Snakemake Integration
# =============================================================================

def run_from_snakemake():
    """Run pipeline from Snakemake rule."""
    output_adata = snakemake.output.adata
    output_qc_metrics = snakemake.output.qc_metrics
    output_annotation_summary = snakemake.output.annotation_summary
    output_figures = snakemake.output.figures
    
    params = snakemake.params
    sample_ids = params.sample_ids
    cellranger_dir = params.cellranger_dir
    qc_params = params.qc_params
    preprocessing_params = params.preprocessing_params
    annotation_params = params.annotation_params
    
    log_file = snakemake.log[0] if snakemake.log else None
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    run_qc_annotation_pipeline(
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
        description='ClusterCatcher QC and annotation pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python scanpy_qc_annotation.py \\
        --config config.yaml \\
        --cellranger-dir /path/to/cellranger \\
        --sample-ids SAMPLE1 SAMPLE2 \\
        --output-dir /path/to/output
        """
    )
    parser.add_argument('--config', required=True, help='Path to config.yaml file')
    parser.add_argument('--cellranger-dir', required=True, help='Base directory containing Cell Ranger outputs')
    parser.add_argument('--sample-ids', nargs='+', required=True, help='List of sample IDs')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    qc_params = config.get('qc', {})
    preprocessing_params = config.get('preprocessing', {})
    annotation_params = config.get('annotation', {})
    
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
        preprocessing_params=preprocessing_params,
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
