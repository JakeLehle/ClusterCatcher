#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cancer_cell_detection.py
========================

Detect cancer cells using dual-model approach with CytoTRACE2 (stemness) and
inferCNV (copy number variation). Cells are classified as cancer if they show
high scores on both models with strong agreement.

Pipeline:
1. Run CytoTRACE2 for stemness/potency scoring
2. Use CytoTRACE2 threshold to define reference normal cells
3. Run inferCNV using normal cells as reference
4. Calculate weighted agreement between models
5. Classify cells as Cancer/Normal based on agreement threshold

Key outputs:
- CytoTRACE2 scores per cell
- CNV scores per cell  
- Final cancer cell classification
- Comprehensive visualizations

Requirements:
    - cytotrace2-py
    - infercnvpy
    - gffutils (for GTF parsing)
    - scanpy, matplotlib, seaborn

Usage:
    Called via Snakemake rule with snakemake.input/output/params
    
    Or standalone:
    python cancer_cell_detection.py --adata input.h5ad --gtf genes.gtf --output output_dir
"""

import os
import sys
import time
import shutil
import logging
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Patch, Rectangle
from matplotlib.text import Text
import seaborn as sns
from scipy.stats import spearmanr, pearsonr, norm
from scipy.signal import convolve, find_peaks
from adjustText import adjust_text

# Suppress warnings
warnings.filterwarnings('ignore')
matplotlib.use('Agg')  # Non-interactive backend for cluster

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================

DEFAULT_CONFIG = {
    'cytotrace2': {
        'species': 'human',
        'max_cells_per_chunk': 200000,
        'seed': 42,
    },
    'infercnv': {
        'window_size': 250,
        'step': 10,
        'dynamic_threshold': 1.5,
    },
    'agreement': {
        'alpha': 0.5,  # Weight between rank and value agreement
        'min_correlation': 0.5,  # Minimum Spearman rho for quartile selection
    },
}


# =============================================================================
# Visualization Helper Functions
# =============================================================================

def gen_mpl_labels(adata, groupby, exclude=(), ax=None, adjust_kwargs=None, 
                   text_kwargs=None, color_by_group=False):
    """Generate non-overlapping labels for UMAP plots."""
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}
    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)

    text_colors = {group: None for group in adata.obs[groupby].cat.categories}

    if color_by_group and groupby + "_colors" in adata.uns:
        for i, group in enumerate(adata.obs[groupby].cat.categories):
            if group in exclude:
                continue
            text_colors[group] = adata.uns[groupby + "_colors"][i]

    if ax is None:
        texts = [plt.text(x=x, y=y, s=k, color=text_colors[k], **text_kwargs) 
                 for k, (x, y) in medians.items()]
    else:
        texts = [ax.text(x=x, y=y, s=k, color=text_colors[k], **text_kwargs) 
                 for k, (x, y) in medians.items()]

    adjust_text(texts, **adjust_kwargs)


def nonoverlapping_UMAP(adata_obj, group_name, output_path):
    """Create UMAP with non-overlapping text labels."""
    cmap = plt.get_cmap('turbo')
    adata_obj.obs[group_name] = adata_obj.obs[group_name].cat.remove_unused_categories()
    value_cat = pd.Categorical(adata_obj.obs[group_name])
    values = np.linspace(0, 1, len(value_cat.categories))
    palette = [cmap(value) for value in values]
    
    combined_effects = [
        pe.withStroke(linewidth=6, foreground="white"),
        pe.withStroke(linewidth=1, foreground="black"),
        pe.Normal()
    ]
    
    with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": 150, "figure.frameon": False}):
        ax = sc.pl.umap(adata_obj, color=group_name, show=False, legend_loc=None, 
                        frameon=False, title='', size=5, palette=palette)
        gen_mpl_labels(
            adata_obj,
            group_name,
            exclude=("None",),
            ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
            text_kwargs=dict(fontsize=32, path_effects=combined_effects),
            color_by_group=True
        )
        fig = ax.get_figure()
        fig.tight_layout()
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()


# =============================================================================
# CytoTRACE2 Functions
# =============================================================================

def safe_makedirs(dir_path):
    """Create directory, removing existing one if it exists."""
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)


def process_cytotrace_chunk(adata_chunk, chunk_name, working_dir, species='human', seed=42):
    """Process a single chunk of cells through CytoTRACE2."""
    from cytotrace2_py.cytotrace2_py import cytotrace2
    
    cytotrace2_results_dir = os.path.join(working_dir, "cytotrace2_results")
    
    # Prepare cell annotations
    adata_chunk.obs['cell_barcodes'] = adata_chunk.obs.index
    adata_annotation_df = adata_chunk.obs[['cell_barcodes', 'final_annotation']]
    
    annotation_path = os.path.join(working_dir, f"cell_annotations_{chunk_name}.txt")
    adata_annotation_df.to_csv(annotation_path, sep="\t", header=True, index=False)
    
    # Check for duplicates
    assert adata_chunk.var_names.is_unique, "Duplicate gene names detected!"
    assert adata_chunk.obs_names.is_unique, "Duplicate cell barcodes detected!"
    
    # Prepare gene expression matrix
    if 'gene_symbol' in adata_chunk.var.columns:
        gene_names = adata_chunk.var['gene_symbol'].values
    else:
        gene_names = adata_chunk.var_names.values
    
    adata_X_df_T = pd.DataFrame.sparse.from_spmatrix(
        scipy.sparse.csr_matrix(adata_chunk.X),
        index=adata_chunk.obs.index,
        columns=gene_names
    ).T
    
    expression_path = os.path.join(working_dir, f"gene_expression_matrix_{chunk_name}.txt")
    adata_X_df_T.to_csv(expression_path, sep="\t", header=True, index=True, chunksize=10000)
    
    # Clean up any existing results
    if os.path.exists(cytotrace2_results_dir):
        shutil.rmtree(cytotrace2_results_dir)
    
    # Run CytoTRACE2
    results = cytotrace2(
        expression_path,
        annotation_path=annotation_path,
        species=species,
        max_cores=os.cpu_count(),
        seed=seed
    )
    
    # Read results
    project_results_path = os.path.join(cytotrace2_results_dir, "cytotrace2_results.txt")
    if not os.path.exists(project_results_path):
        raise FileNotFoundError(f"CytoTRACE2 output missing at {project_results_path}")
    
    cytotrace_txt = pd.read_csv(project_results_path, sep='\t')
    
    # Clean up temporary files
    for path in [annotation_path, expression_path]:
        if os.path.exists(path):
            os.remove(path)
    
    return cytotrace_txt


def run_cytotrace2(adata, working_dir, config):
    """
    Run CytoTRACE2 on all cells, processing in chunks if necessary.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data with 'final_annotation' in obs and 'gene_symbol' in var
    working_dir : str
        Working directory for temporary files
    config : dict
        Configuration parameters
        
    Returns
    -------
    pd.DataFrame
        CytoTRACE2 results for all cells
    """
    logger.info("Running CytoTRACE2 analysis...")
    
    # Ensure gene names are unique
    adata.var_names_make_unique()
    if 'gene_symbol' in adata.var.columns:
        adata.var.index = adata.var['gene_symbol']
        adata.var_names_make_unique()
    
    MAX_CELLS_PER_CHUNK = config.get('max_cells_per_chunk', 200000)
    RANDOM_SEED = config.get('seed', 42)
    species = config.get('species', 'human')
    
    # Set up directories
    chunk_results_dir = os.path.join(working_dir, "cytotrace_chunk_results")
    final_results_dir = os.path.join(working_dir, "cytotrace_final_results")
    cytotrace2_results_dir = os.path.join(working_dir, "cytotrace2_results")
    
    safe_makedirs(chunk_results_dir)
    safe_makedirs(final_results_dir)
    
    all_cytotrace_results = []
    
    # Get unique sample/project IDs if available
    if 'sample_id' in adata.obs.columns:
        project_ids = adata.obs['sample_id'].unique()
        project_key = 'sample_id'
    elif 'series_id' in adata.obs.columns:
        project_ids = adata.obs['series_id'].unique()
        project_key = 'series_id'
    else:
        # Process all cells together
        project_ids = ['all']
        project_key = None
    
    np.random.seed(RANDOM_SEED)
    
    for project_id in project_ids:
        logger.info(f"\n{'='*50}")
        logger.info(f"Processing: {project_id}")
        project_start = time.time()
        
        if project_key:
            project_mask = adata.obs[project_key] == project_id
            adata_project = adata[project_mask].copy()
        else:
            adata_project = adata.copy()
        
        n_cells = adata_project.n_obs
        logger.info(f"Total cells: {n_cells}")
        
        cell_barcodes = adata_project.obs.index.values.copy()
        shuffled_positions = np.arange(n_cells)
        np.random.shuffle(shuffled_positions)
        
        n_chunks = max(1, (n_cells + MAX_CELLS_PER_CHUNK - 1) // MAX_CELLS_PER_CHUNK)
        logger.info(f"Splitting into {n_chunks} chunks")
        
        for chunk_idx in range(n_chunks):
            chunk_start = time.time()
            chunk_name = f"{project_id}_chunk{chunk_idx+1}"
            logger.info(f"\nProcessing {chunk_name}...")
            
            start_idx = chunk_idx * MAX_CELLS_PER_CHUNK
            end_idx = min((chunk_idx + 1) * MAX_CELLS_PER_CHUNK, n_cells)
            chunk_positions = shuffled_positions[start_idx:end_idx]
            chunk_cell_barcodes = cell_barcodes[chunk_positions]
            
            adata_chunk = adata_project[adata_project.obs.index.isin(chunk_cell_barcodes)].copy()
            
            chunk_results = process_cytotrace_chunk(
                adata_chunk, chunk_name, working_dir, species, RANDOM_SEED
            )
            chunk_results = chunk_results.rename(columns={'Unnamed: 0': 'cell_barcodes'})
            
            if project_key:
                chunk_results[project_key] = project_id
            chunk_results['chunk_id'] = chunk_name
            all_cytotrace_results.append(chunk_results)
            
            # Move results to chunk-specific directory
            chunk_output_dir = os.path.join(chunk_results_dir, chunk_name)
            safe_makedirs(chunk_output_dir)
            
            if os.path.exists(cytotrace2_results_dir):
                unique_dir_name = f"cytotrace2_results_{chunk_name}"
                shutil.move(cytotrace2_results_dir, os.path.join(chunk_output_dir, unique_dir_name))
            
            logger.info(f"Completed {chunk_name} in {time.time()-chunk_start:.2f}s")
    
    # Combine all results
    final_cytotrace_results = pd.concat(all_cytotrace_results, axis=0)
    
    # Merge with original annotations
    Cell_type_anno_df = pd.DataFrame(adata.obs[['final_annotation']])
    Cell_type_anno_df['cell_barcodes'] = Cell_type_anno_df.index
    
    final_cytotrace_results = pd.merge(
        final_cytotrace_results,
        Cell_type_anno_df,
        on=['cell_barcodes'],
        how='outer'
    )
    
    # Save final results
    final_output_path = os.path.join(final_results_dir, "combined_cytotrace_results.txt")
    final_cytotrace_results.to_csv(final_output_path, sep="\t", header=True, index=False)
    
    logger.info(f"\nCytoTRACE2 completed: {len(final_cytotrace_results)} cells processed")
    
    return final_cytotrace_results


def detect_cytotrace_threshold(scores, figures_dir=None, sigma=3):
    """
    Detect threshold between normal and cancer populations using peak detection.
    
    Uses Gaussian convolution to smooth the histogram and finds the valley
    between the first two peaks.
    
    Parameters
    ----------
    scores : array-like
        CytoTRACE2 scores
    figures_dir : str, optional
        Directory to save threshold plot
    sigma : float
        Gaussian kernel sigma for smoothing
        
    Returns
    -------
    float
        Threshold value
    """
    logger.info("Detecting CytoTRACE2 threshold...")
    
    counts, edges = np.histogram(scores, bins=100, density=False)
    counts_data = np.asarray(counts)
    kernel_size = int(2 * np.ceil(3 * sigma) + 0.1)
    
    # Create Gaussian kernel
    x = np.linspace(-2.5*sigma, 2.5*sigma, kernel_size)
    kernel = norm.pdf(x, loc=0, scale=sigma)
    kernel = kernel / np.sum(kernel)
    convolved = convolve(counts_data, kernel, mode='same')
    peaks, _ = find_peaks(convolved)
    
    arr = np.linspace(0, 1, 100)
    
    # Find valley between first two peaks
    if len(peaks) >= 2:
        valley_idx = np.argmin(convolved[peaks[0]:peaks[1]]) + peaks[0]
        threshold = arr[valley_idx]
    else:
        # Fallback to median
        threshold = np.median(scores)
        logger.warning(f"Could not detect peaks, using median threshold: {threshold:.3f}")
    
    logger.info(f"Detected threshold: {threshold:.3f}")
    
    # Generate plot if directory provided
    if figures_dir:
        plt.figure(figsize=(14, 4))
        
        # Original data
        plt.subplot(1, 2, 1)
        plt.stem(arr, counts_data, linefmt='C0-', markerfmt='C0o', basefmt='C0-')
        plt.ylabel('Cell count', fontsize=18)
        plt.xlabel('CytoTRACE2 Score', fontsize=18)
        plt.title('Original Data', fontsize=20)
        
        # Smoothed data
        ax = plt.subplot(1, 2, 2)
        plt.stem(arr, convolved, linefmt='C1-', markerfmt='C1o', basefmt='C1-')
        ax.axvspan(0, threshold, facecolor='lightblue', alpha=0.3, zorder=0)
        ax.axvspan(threshold, 1, facecolor='mistyrose', alpha=0.3, zorder=0)
        plt.axvline(x=threshold, color='black', linestyle='--', linewidth=2, 
                    label=f'Threshold: {threshold:.3f}')
        
        if len(peaks) >= 2:
            plt.scatter(arr[peaks], convolved[peaks], color='black', s=100, zorder=3, label='Peaks')
        
        plt.text(threshold/2, 0.5*max(convolved), 'Normal', 
                 ha='center', va='center', fontsize=18, weight='bold')
        plt.text((threshold+1)/2, 0.5*max(convolved), 'Cancer', 
                 ha='center', va='center', fontsize=18, weight='bold')
        
        plt.ylabel('Cell count (smoothed)', fontsize=18)
        plt.xlabel('CytoTRACE2 Score', fontsize=18)
        plt.title('Smoothed Data with Population Separation', fontsize=18)
        plt.legend(loc='upper right', fontsize=14)
        
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'cytotrace_threshold_detection.pdf'), 
                    bbox_inches='tight')
        plt.close()
    
    return threshold


# =============================================================================
# InferCNV Functions
# =============================================================================

def add_chromosomal_info(adata, gtf_path, working_dir):
    """
    Add chromosomal information to adata.var using GTF file.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data
    gtf_path : str
        Path to GTF annotation file
    working_dir : str
        Working directory for database cache
        
    Returns
    -------
    AnnData
        AnnData with chromosomal info in var
    """
    import gffutils
    
    logger.info("Adding chromosomal information from GTF...")
    
    db_path = os.path.join(working_dir, "genes.db")
    
    # Create or load database
    if not os.path.exists(db_path):
        logger.info("Creating gene database (this may take a few minutes)...")
        db = gffutils.create_db(
            gtf_path,
            dbfn=db_path,
            force=True,
            keep_order=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
        logger.info(f"Database created at: {db_path}")
    else:
        logger.info(f"Using existing database at: {db_path}")
    
    db = gffutils.FeatureDB(db_path)
    
    # Get gene IDs from adata
    if 'gene_ids' in adata.var.columns:
        gene_ids = adata.var['gene_ids'].values
    else:
        gene_ids = adata.var_names.values
    
    # Extract gene info
    gene_data = []
    found = 0
    for gene_id in gene_ids:
        try:
            gene = db[gene_id]
            gene_data.append({
                'chromosome': gene.chrom,
                'start': gene.start,
                'end': gene.end
            })
            found += 1
        except:
            gene_data.append({
                'chromosome': None,
                'start': None,
                'end': None
            })
    
    logger.info(f"Found chromosomal info for {found}/{len(gene_ids)} genes")
    
    # Add to adata.var
    gene_df = pd.DataFrame(gene_data, index=adata.var.index)
    
    for col in ['chromosome', 'start', 'end']:
        adata.var[col] = gene_df[col].values
    
    return adata


def run_infercnv(adata, reference_key, reference_cat, config, figures_dir=None):
    """
    Run inferCNV analysis.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data with chromosomal info in var
    reference_key : str
        Column in obs containing reference labels
    reference_cat : list
        Categories to use as reference (normal cells)
    config : dict
        InferCNV configuration
    figures_dir : str, optional
        Directory for plots
        
    Returns
    -------
    AnnData
        AnnData with CNV scores
    """
    import infercnvpy as cnv
    
    logger.info("Running inferCNV analysis...")
    
    # Ensure gene_ids is the index for inferCNV
    if 'gene_ids' in adata.var.columns:
        adata.var = adata.var.set_index('gene_ids')
    
    # Run inferCNV
    cnv.tl.infercnv(
        adata,
        reference_key=reference_key,
        reference_cat=reference_cat,
        window_size=config.get('window_size', 250),
    )
    
    # Generate heatmaps
    if figures_dir:
        cnv.pl.chromosome_heatmap(adata, groupby=reference_key, show=False)
        plt.savefig(os.path.join(figures_dir, 'chromosome_heatmap_reference.pdf'), 
                    bbox_inches='tight')
        plt.close()
        
        if 'final_annotation' in adata.obs.columns:
            cnv.pl.chromosome_heatmap(adata, groupby="final_annotation", show=False)
            plt.savefig(os.path.join(figures_dir, 'chromosome_heatmap_celltypes.pdf'), 
                        bbox_inches='tight')
            plt.close()
    
    # CNV-based clustering
    logger.info("Performing CNV-based clustering...")
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata)
    cnv.tl.umap(adata)
    cnv.tl.cnv_score(adata)
    
    if figures_dir:
        cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True, show=False)
        plt.savefig(os.path.join(figures_dir, 'chromosome_heatmap_cnv_clusters.pdf'), 
                    bbox_inches='tight')
        plt.close()
    
    logger.info(f"InferCNV completed. CNV scores range: {adata.obs['cnv_score'].min():.3f} - {adata.obs['cnv_score'].max():.3f}")
    
    return adata


# =============================================================================
# Agreement Scoring Functions
# =============================================================================

def calculate_agreement_score(cytotrace_scores, cnv_scores, alpha=0.5):
    """
    Calculate weighted agreement between CytoTRACE2 and CNV scores.
    
    Parameters
    ----------
    cytotrace_scores : array-like
        CytoTRACE2 scores (0-1)
    cnv_scores : array-like
        CNV scores
    alpha : float
        Weight between rank agreement (alpha) and value agreement (1-alpha)
        
    Returns
    -------
    np.ndarray
        Agreement scores (0-1, higher = better agreement)
    """
    array1 = np.array(cytotrace_scores)
    array2 = np.array(cnv_scores)
    
    # Rank-based agreement
    rank1 = np.argsort(np.argsort(array1))
    rank2 = np.argsort(np.argsort(array2))
    normalized_rank1 = rank1 / len(array1)
    normalized_rank2 = rank2 / len(array2)
    rank_agreement = 1 - np.abs(normalized_rank1 - normalized_rank2)
    
    # Value-based agreement (normalize to 0-1)
    norm_array1 = (array1 - np.min(array1)) / (np.max(array1) - np.min(array1) + 1e-10)
    norm_array2 = (array2 - np.min(array2)) / (np.max(array2) - np.min(array2) + 1e-10)
    value_agreement = 1 - np.abs(norm_array1 - norm_array2)
    
    # Combined agreement
    agreement_score = (alpha * rank_agreement) + ((1 - alpha) * value_agreement)
    
    return agreement_score


def find_agreement_threshold(agreement_scores, cytotrace_scores, cnv_scores, 
                             min_correlation=0.5, figures_dir=None):
    """
    Find optimal agreement threshold using quartile-based correlation analysis.
    
    Parameters
    ----------
    agreement_scores : array-like
        Agreement scores between models
    cytotrace_scores : array-like
        CytoTRACE2 scores
    cnv_scores : array-like
        CNV scores
    min_correlation : float
        Minimum Spearman correlation to accept quartile
    figures_dir : str, optional
        Directory for plots
        
    Returns
    -------
    tuple
        (selected_threshold, quartile_info)
    """
    logger.info("Finding optimal agreement threshold...")
    
    quartiles = np.percentile(agreement_scores, [25, 50, 75])
    q1, q2, q3 = quartiles
    
    quartile_ranges = [(0, q1), (q1, q2), (q2, q3), (q3, 1)]
    quartile_names = ['Q1', 'Q2', 'Q3', 'Q4']
    
    selected_threshold = None
    selected_quartile = None
    quartile_info = []
    
    for q_range, q_name in zip(quartile_ranges, quartile_names):
        mask = (agreement_scores >= q_range[0]) & (agreement_scores < q_range[1])
        cyto_q = cytotrace_scores[mask]
        cnv_q = cnv_scores[mask]
        
        if len(cyto_q) > 1:
            rho, _ = spearmanr(cyto_q, cnv_q)
            quartile_info.append({
                'quartile': q_name,
                'range': q_range,
                'n_cells': len(cyto_q),
                'spearman_rho': rho
            })
            
            logger.info(f"{q_name} Spearman ρ: {rho:.3f} (n={len(cyto_q)})")
            
            if rho > min_correlation and selected_threshold is None:
                selected_threshold = q_range[0]
                selected_quartile = q_name
                logger.info(f"SELECTED: {q_name} with threshold > {q_range[0]:.3f}")
    
    # Fallback
    if selected_threshold is None:
        selected_threshold = q1
        logger.warning(f"No quartile met correlation threshold. Using Q1: {q1:.3f}")
    
    # Generate plots
    if figures_dir:
        # Quartile correlation bar plot
        plt.figure(figsize=(4, 6))
        for i, info in enumerate(quartile_info):
            plt.bar(i, info['spearman_rho'], label=f"{info['quartile']} (n={info['n_cells']})")
        plt.xticks(range(4), ['Q1\n(25%)', 'Q2\n(50%)', 'Q3\n(75%)', 'Q4\n(100%)'], fontsize=12)
        plt.ylabel('Spearman Correlation (ρ)', fontsize=14)
        plt.axhline(y=min_correlation, color='r', linestyle='--', label=f'Threshold: {min_correlation}')
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'quartile_correlation_analysis.pdf'), 
                    bbox_inches='tight')
        plt.close()
        
        # Agreement histogram
        colors = ['#ff0000', '#ff9900', '#ffff00', '#00cc00']
        plt.figure(figsize=(8, 6))
        ax = plt.gca()
        
        n, bins, patches = plt.hist(agreement_scores, bins=50, edgecolor='white', 
                                     alpha=0.8, range=(0.3, 1))
        
        for i in range(len(patches)):
            bin_center = (bins[i] + bins[i+1]) / 2
            if bin_center <= q1:
                patches[i].set_facecolor(colors[0])
            elif bin_center <= q2:
                patches[i].set_facecolor(colors[1])
            elif bin_center <= q3:
                patches[i].set_facecolor(colors[2])
            else:
                patches[i].set_facecolor(colors[3])
        
        for q, label in zip(quartiles, ['Q1', 'Q2', 'Q3']):
            plt.axvline(q, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
            plt.text(q, plt.ylim()[1]*0.9, f'{q:.2f}', fontsize=12, ha='center')
        
        plt.xlabel('Agreement Score', fontsize=14)
        plt.ylabel('Number of cells', fontsize=14)
        
        legend_elements = [
            Patch(facecolor=colors[0], label='Q1 (0-25%)'),
            Patch(facecolor=colors[1], label='Q2 (25-50%)'),
            Patch(facecolor=colors[2], label='Q3 (50-75%)'),
            Patch(facecolor=colors[3], label='Q4 (75-100%)')
        ]
        plt.legend(handles=legend_elements, loc='upper left')
        
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'agreement_score_distribution.pdf'), 
                    bbox_inches='tight')
        plt.close()
    
    return selected_threshold, quartile_info


# =============================================================================
# Final Classification Functions
# =============================================================================

def classify_cancer_cells(adata, cytotrace_threshold, agreement_threshold, figures_dir=None):
    """
    Final classification of cancer vs normal cells.
    
    A cell is classified as cancer if:
    1. CytoTRACE2 score > cytotrace_threshold AND
    2. Agreement score > agreement_threshold
    
    Parameters
    ----------
    adata : AnnData
        AnnData with CytoTRACE2_Score, cnv_score, and agreement_score
    cytotrace_threshold : float
        CytoTRACE2 threshold for cancer
    agreement_threshold : float
        Minimum agreement between models
    figures_dir : str, optional
        Directory for plots
        
    Returns
    -------
    AnnData
        AnnData with Final_cancer_cell_status
    """
    logger.info("Classifying cancer cells...")
    
    # Get cells with high agreement
    high_agreement = adata.obs['agreement_score'] > agreement_threshold
    
    # Get cells with high CytoTRACE2 score
    high_cytotrace = adata.obs['CytoTRACE2_Score'] > cytotrace_threshold
    
    # Cancer cells must meet both criteria
    cancer_mask = high_agreement & high_cytotrace
    
    adata.obs['Final_cancer_cell_status'] = np.where(
        cancer_mask, 'Cancer cell', 'Normal cell'
    )
    adata.obs['Final_cancer_cell_status'] = adata.obs['Final_cancer_cell_status'].astype('category')
    
    n_cancer = cancer_mask.sum()
    n_total = len(adata.obs)
    logger.info(f"Classified {n_cancer} cancer cells ({100*n_cancer/n_total:.1f}%)")
    
    # Generate final visualization
    if figures_dir:
        generate_final_plots(adata, figures_dir)
    
    return adata


def generate_final_plots(adata, figures_dir):
    """Generate comprehensive final visualization plots."""
    logger.info("Generating final visualization plots...")
    
    from matplotlib import rcParams
    rcParams['font.size'] = 14
    rcParams['axes.titlesize'] = 16
    
    # Multi-panel UMAP
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 13),
                                                   gridspec_kw={"wspace": 0.2, "hspace": 0.2})
    
    # CytoTRACE2 Score
    sc.pl.umap(adata, color="CytoTRACE2_Score", ax=ax1, title="CytoTRACE2 Score",
               show=False, frameon=False, color_map='plasma', size=5)
    
    # CNV Score
    sc.pl.umap(adata, color="cnv_score", ax=ax2, title="CNV Score",
               show=False, frameon=False, color_map='plasma', size=5)
    
    # Cancer Status
    cancer_palette = {'Cancer cell': '#FF7F0E', 'Normal cell': '#1F77B4'}
    sc.pl.umap(adata, color="Final_cancer_cell_status", ax=ax3, title="Cancer Cell Status",
               show=False, frameon=False, palette=cancer_palette, size=5)
    
    # Cell Types
    sc.pl.umap(adata, color="final_annotation", ax=ax4, title="Cell Types",
               show=False, frameon=False, size=5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'cancer_detection_summary.pdf'), 
                bbox_inches='tight', dpi=300)
    plt.close()
    
    # Correlation plot
    fig, ax = plt.subplots(figsize=(6, 6))
    
    sc = ax.scatter(
        adata.obs['CytoTRACE2_Score'],
        adata.obs['cnv_score'],
        c=adata.obs['agreement_score'],
        cmap="RdYlGn",
        alpha=0.5,
        s=5
    )
    
    rho, p = spearmanr(adata.obs['CytoTRACE2_Score'], adata.obs['cnv_score'])
    ax.text(0.05, 0.95, f"Spearman ρ = {rho:.2f}\n(p = {p:.2e})",
            transform=ax.transAxes, fontsize=12, va='top')
    
    ax.set_xlabel("CytoTRACE2 Score", fontsize=14)
    ax.set_ylabel("CNV Score", fontsize=14)
    
    cbar = plt.colorbar(sc)
    cbar.set_label("Agreement Score", fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'score_correlation.pdf'), 
                bbox_inches='tight')
    plt.close()
    
    # Cancer cell type distribution
    if 'final_annotation' in adata.obs.columns:
        cancer_cells = adata.obs[adata.obs['Final_cancer_cell_status'] == 'Cancer cell']
        if len(cancer_cells) > 0:
            celltype_counts = cancer_cells['final_annotation'].value_counts()
            
            fig, ax = plt.subplots(figsize=(10, 6))
            celltype_counts.plot(kind='bar', ax=ax)
            ax.set_xlabel("Cell Type", fontsize=14)
            ax.set_ylabel("Number of Cancer Cells", fontsize=14)
            ax.set_title("Cancer Cells by Cell Type", fontsize=16)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(figures_dir, 'cancer_celltype_distribution.pdf'), 
                        bbox_inches='tight')
            plt.close()


# =============================================================================
# Main Pipeline Function
# =============================================================================

def run_cancer_detection_pipeline(
    adata_path,
    gtf_path,
    output_dir,
    config=None,
):
    """
    Run the complete cancer cell detection pipeline.
    
    Parameters
    ----------
    adata_path : str
        Path to annotated H5AD file (must have 'final_annotation' in obs)
    gtf_path : str
        Path to GTF annotation file
    output_dir : str
        Output directory
    config : dict, optional
        Configuration parameters
        
    Returns
    -------
    AnnData
        AnnData with cancer cell classifications
    """
    if config is None:
        config = DEFAULT_CONFIG
    
    os.makedirs(output_dir, exist_ok=True)
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    logger.info("="*60)
    logger.info("Cancer Cell Detection Pipeline")
    logger.info("="*60)
    
    # Load data
    logger.info(f"\nLoading data from {adata_path}")
    adata = sc.read_h5ad(adata_path)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Ensure required columns exist
    if 'final_annotation' not in adata.obs.columns:
        raise ValueError("AnnData must have 'final_annotation' column in obs")
    
    # Step 1: Run CytoTRACE2
    logger.info("\n" + "="*60)
    logger.info("Step 1: CytoTRACE2 Analysis")
    logger.info("="*60)
    
    cytotrace_results = run_cytotrace2(adata, output_dir, config.get('cytotrace2', {}))
    
    # Merge results with adata
    cytotrace_txt = cytotrace_results[['cell_barcodes', 'CytoTRACE2_Score', 'CytoTRACE2_Potency']].copy()
    cytotrace_txt = cytotrace_txt.drop_duplicates(subset=['cell_barcodes'])
    cytotrace_txt = cytotrace_txt.set_index('cell_barcodes')
    
    adata.obs['cell_barcodes'] = adata.obs.index
    adata.obs = adata.obs.join(cytotrace_txt, on='cell_barcodes', how='left')
    
    # Detect threshold
    cytotrace_threshold = detect_cytotrace_threshold(
        adata.obs['CytoTRACE2_Score'].dropna(),
        figures_dir
    )
    
    # Create initial cancer/normal labels for inferCNV reference
    adata.obs['infer_cnv_reference'] = np.where(
        adata.obs['CytoTRACE2_Score'] > cytotrace_threshold, 'Cancer', 'Normal'
    )
    
    # Save intermediate result
    adata.write(os.path.join(output_dir, 'adata_cytotrace.h5ad'))
    
    # Step 2: Add chromosomal info and run inferCNV
    logger.info("\n" + "="*60)
    logger.info("Step 2: InferCNV Analysis")
    logger.info("="*60)
    
    adata = add_chromosomal_info(adata, gtf_path, output_dir)
    
    adata = run_infercnv(
        adata,
        reference_key='infer_cnv_reference',
        reference_cat=['Normal'],
        config=config.get('infercnv', {}),
        figures_dir=figures_dir
    )
    
    # Save intermediate result
    adata.write(os.path.join(output_dir, 'adata_infercnv.h5ad'))
    
    # Step 3: Calculate agreement and final classification
    logger.info("\n" + "="*60)
    logger.info("Step 3: Agreement Analysis and Final Classification")
    logger.info("="*60)
    
    # Calculate agreement
    agreement_config = config.get('agreement', {})
    adata.obs['agreement_score'] = calculate_agreement_score(
        adata.obs['CytoTRACE2_Score'].values,
        adata.obs['cnv_score'].values,
        alpha=agreement_config.get('alpha', 0.5)
    )
    
    # Find agreement threshold
    agreement_threshold, quartile_info = find_agreement_threshold(
        adata.obs['agreement_score'].values,
        adata.obs['CytoTRACE2_Score'].values,
        adata.obs['cnv_score'].values,
        min_correlation=agreement_config.get('min_correlation', 0.5),
        figures_dir=figures_dir
    )
    
    # Final classification
    adata = classify_cancer_cells(
        adata, cytotrace_threshold, agreement_threshold, figures_dir
    )
    
    # Save final result
    adata.write(os.path.join(output_dir, 'adata_cancer_detected.h5ad'))
    
    # Save summary
    summary = {
        'total_cells': adata.n_obs,
        'cancer_cells': (adata.obs['Final_cancer_cell_status'] == 'Cancer cell').sum(),
        'normal_cells': (adata.obs['Final_cancer_cell_status'] == 'Normal cell').sum(),
        'cytotrace_threshold': float(cytotrace_threshold),
        'agreement_threshold': float(agreement_threshold),
        'spearman_correlation': float(spearmanr(
            adata.obs['CytoTRACE2_Score'], adata.obs['cnv_score']
        )[0]),
    }
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(os.path.join(output_dir, 'cancer_detection_summary.tsv'), 
                      sep='\t', index=False)
    
    logger.info("\n" + "="*60)
    logger.info("Pipeline completed successfully!")
    logger.info("="*60)
    logger.info(f"Total cells: {summary['total_cells']}")
    logger.info(f"Cancer cells: {summary['cancer_cells']} ({100*summary['cancer_cells']/summary['total_cells']:.1f}%)")
    logger.info(f"Output directory: {output_dir}")
    
    return adata


# =============================================================================
# Snakemake Integration
# =============================================================================

def run_from_snakemake():
    """Run from Snakemake rule."""
    
    # Get inputs
    adata_path = snakemake.input.adata
    gtf_path = snakemake.params.gtf_path
    
    # Get outputs
    output_dir = os.path.dirname(snakemake.output.adata)
    
    # Get config
    config = {
        'cytotrace2': snakemake.config.get('dysregulation', {}).get('cytotrace2', {}),
        'infercnv': snakemake.config.get('dysregulation', {}).get('infercnv', {}),
        'agreement': snakemake.config.get('dysregulation', {}).get('agreement', {}),
    }
    
    # Set up logging
    if snakemake.log:
        file_handler = logging.FileHandler(snakemake.log[0])
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Run pipeline
    run_cancer_detection_pipeline(
        adata_path=adata_path,
        gtf_path=gtf_path,
        output_dir=output_dir,
        config=config,
    )


# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Detect cancer cells using CytoTRACE2 and inferCNV'
    )
    parser.add_argument('--adata', required=True, help='Input H5AD file')
    parser.add_argument('--gtf', required=True, help='GTF annotation file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--species', default='human', help='Species (default: human)')
    parser.add_argument('--max-cells-chunk', type=int, default=200000,
                        help='Max cells per CytoTRACE2 chunk')
    
    args = parser.parse_args()
    
    config = {
        'cytotrace2': {
            'species': args.species,
            'max_cells_per_chunk': args.max_cells_chunk,
        },
        'infercnv': {},
        'agreement': {},
    }
    
    run_cancer_detection_pipeline(
        adata_path=args.adata,
        gtf_path=args.gtf,
        output_dir=args.output,
        config=config,
    )


if __name__ == '__main__':
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        main()
