#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
signature_analysis.py
=====================

Semi-supervised COSMIC signature analysis for single-cell mutation data.

This script:
1. Converts mutations to 96-trinucleotide context matrix
2. Extracts relevant COSMIC signatures (all or HNSCC-specific)
3. Uses scree plot elbow detection for signature selection (optional)
4. Fits signatures using Non-Negative Least Squares (NNLS)
5. Evaluates reconstruction quality using Frobenius norm
6. Adds signature weights to AnnData
7. Generates comprehensive visualizations

Author: Jake Lehle
Date: 2025
"""

import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import logging
from pathlib import Path
from datetime import datetime
from scipy.optimize import nnls
from scipy.stats import pearsonr, spearmanr
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

# COSMIC signature colors
COSMIC_COLORS = {
    'C>A': '#1EBFF0',
    'C>G': '#050708',
    'C>T': '#E62725',
    'T>A': '#CBCACB',
    'T>C': '#A1CE63',
    'T>G': '#EDB6C2'
}

# Standard 96 trinucleotide contexts (used for empty matrix creation)
ALL_MUTATION_TYPES = [
    f"{five}[{ref}>{alt}]{three}"
    for ref in ['C', 'T']
    for alt in (['A', 'G', 'T'] if ref == 'C' else ['A', 'C', 'G'])
    for five in ['A', 'C', 'G', 'T']
    for three in ['A', 'C', 'G', 'T']
]


def reverse_complement(sequence):
    """Generate reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in sequence[::-1]])


def get_mutation_type(ref, alt, context):
    """Determine mutation type in SigProfiler format"""
    if len(ref) != 1 or len(alt) != 1:
        return None
    if ref in ['C', 'T']:
        return f"{context[0]}[{ref}>{alt}]{context[2]}"
    return None


def create_empty_mutation_matrix(output_file):
    """Create an empty mutation matrix with standard 96 contexts."""
    result_df = pd.DataFrame(index=ALL_MUTATION_TYPES, columns=[], dtype=int)
    result_df.to_csv(output_file, sep='\t')
    logger.info(f"Created empty mutation matrix at {output_file}")
    return result_df


def process_mutations(input_file, output_file):
    """Process mutations into 96-trinucleotide context matrix."""
    logger.info(f"Processing mutations from {input_file}...")
    
    # Check if file exists
    if not os.path.exists(input_file):
        logger.warning(f"Mutations file not found: {input_file}")
        return create_empty_mutation_matrix(output_file)
    
    # Check if file is empty
    file_size = os.path.getsize(input_file)
    if file_size == 0:
        logger.warning(f"Mutations file is empty (0 bytes): {input_file}")
        return create_empty_mutation_matrix(output_file)
    
    # Try to find header line starting with #
    header = None
    line_count = 0
    with open(input_file, 'r') as f:
        for line in f:
            line_count += 1
            if line.startswith('#'):
                header = line.lstrip('#').strip().split('\t')
                break
    
    # Handle missing header or empty file
    if header is None:
        logger.warning(f"No header line (starting with #) found in {input_file}")
        logger.warning(f"File has {line_count} lines but no valid header")
        
        # Check if file has any content at all
        if line_count == 0:
            logger.warning("File appears to be empty")
            return create_empty_mutation_matrix(output_file)
        
        # Try to read without header comment
        try:
            df = pd.read_csv(input_file, sep='\t')
            if df.empty or len(df.columns) == 0:
                logger.warning("File has no data columns")
                return create_empty_mutation_matrix(output_file)
            # Use the DataFrame's columns as header
            header = list(df.columns)
            logger.info(f"Using column names from file: {header[:5]}...")
        except Exception as e:
            logger.warning(f"Could not parse mutations file: {e}")
            return create_empty_mutation_matrix(output_file)
    else:
        # Read with the header we found
        try:
            df = pd.read_csv(input_file, sep='\t', comment='#', names=header)
        except Exception as e:
            logger.warning(f"Error reading mutations file with header: {e}")
            return create_empty_mutation_matrix(output_file)
    
    # Check if DataFrame is empty
    if df.empty:
        logger.warning("Mutations DataFrame is empty after reading")
        return create_empty_mutation_matrix(output_file)
    
    # Check for required columns
    required_cols = ['REF', 'ALT_expected', 'REF_TRI', 'ALT_TRI', 'CB']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.warning(f"Missing required columns: {missing_cols}")
        logger.warning(f"Available columns: {list(df.columns)}")
        return create_empty_mutation_matrix(output_file)
    
    # Filter to SNVs only
    df = df[(df['REF'].str.len() == 1) & (df['ALT_expected'].str.len() == 1)]
    
    if df.empty:
        logger.warning("No SNVs found after filtering")
        return create_empty_mutation_matrix(output_file)
    
    # Filter out N-containing contexts
    df = df[~df['REF_TRI'].str.contains('N', na=True) & ~df['ALT_TRI'].str.contains('N', na=True)]
    
    if df.empty:
        logger.warning("No mutations remaining after filtering N-containing contexts")
        return create_empty_mutation_matrix(output_file)
    
    # Convert to pyrimidine context
    for idx, row in df.iterrows():
        ref = row['REF']
        alt = row['ALT_expected']
        if ref in ['G', 'A']:
            df.at[idx, 'REF'] = reverse_complement(ref)
            df.at[idx, 'ALT_expected'] = reverse_complement(alt)
            df.at[idx, 'REF_TRI'] = reverse_complement(row['REF_TRI'])
            df.at[idx, 'ALT_TRI'] = reverse_complement(row['ALT_TRI'])
    
    # Count mutations per cell per context
    mutation_counts = defaultdict(lambda: defaultdict(int))
    for _, row in df.iterrows():
        cb = row['CB']
        ref = row['REF']
        alt = row['ALT_expected']
        context = row['REF_TRI']
        if ref not in ['C', 'T']:
            continue
        if f"{ref}>{alt}" not in ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']:
            continue
        mutation_type = get_mutation_type(ref, alt, context)
        if mutation_type:
            mutation_counts[cb][mutation_type] += 1
    
    # Check if we have any mutations
    if not mutation_counts:
        logger.warning("No valid mutations found after processing")
        return create_empty_mutation_matrix(output_file)
    
    # Build result DataFrame
    cell_barcodes = list(mutation_counts.keys())
    data_dict = {mut_type: [mutation_counts[cb].get(mut_type, 0) for cb in cell_barcodes] for mut_type in ALL_MUTATION_TYPES}
    result_df = pd.DataFrame(data_dict, index=cell_barcodes).T
    result_df = result_df.reindex(ALL_MUTATION_TYPES, fill_value=0).astype(int)
    result_df.to_csv(output_file, sep='\t')
    
    logger.info(f"Matrix shape: {result_df.shape} ({result_df.shape[1]} cells with mutations)")
    return result_df


def extract_cosmic_signatures(cosmic_file, output_dir=None, hnscc_only=False):
    """Extract COSMIC signatures."""
    logger.info(f"Loading COSMIC signatures from: {cosmic_file}")
    
    if cosmic_file is None or not os.path.exists(cosmic_file):
        raise FileNotFoundError(f"COSMIC signatures file not found: {cosmic_file}")
    
    cosmic_all = pd.read_csv(cosmic_file, sep='\t', index_col=0)
    
    if hnscc_only:
        hnscc_sigs = ["SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS13", "SBS16", 
                      "SBS17a", "SBS17b", "SBS18", "SBS33", "SBS40a", "SBS40b", "SBS40c"]
        available = [s for s in hnscc_sigs if s in cosmic_all.columns]
        sigs = cosmic_all[available].copy()
        logger.info(f"Using HNSCC-specific signatures: {available}")
    else:
        sigs = cosmic_all.copy()
        logger.info(f"Using all {len(sigs.columns)} COSMIC signatures")
    
    # Normalize columns to sum to 1
    col_sums = sigs.sum(axis=0)
    if not np.allclose(col_sums, 1.0, atol=1e-5):
        sigs = sigs / col_sums
    
    if output_dir:
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        sigs.to_csv(Path(output_dir) / "cosmic_signatures_used.txt", sep='\t', float_format='%.6f')
    
    return sigs


def fit_signatures_nnls(mutation_matrix, signature_matrix, verbose=True):
    """Fit signatures using NNLS."""
    if verbose:
        logger.info("Fitting signatures using NNLS...")
    
    # Align indices
    if not all(mutation_matrix.index == signature_matrix.index):
        common_contexts = mutation_matrix.index.intersection(signature_matrix.index)
        mutation_matrix = mutation_matrix.loc[common_contexts]
        signature_matrix = signature_matrix.loc[common_contexts]
    
    n_cells = mutation_matrix.shape[1]
    n_sigs = signature_matrix.shape[1]
    
    if n_cells == 0:
        logger.warning("No cells in mutation matrix")
        return {
            'weights': pd.DataFrame(index=signature_matrix.columns, columns=[]),
            'residuals': np.array([]),
            'reconstruction': pd.DataFrame(index=mutation_matrix.index, columns=[])
        }
    
    X = mutation_matrix.values
    H = signature_matrix.values
    W = np.zeros((n_sigs, n_cells))
    residuals = np.zeros(n_cells)
    
    for i in range(n_cells):
        try:
            weights, residual = nnls(H, X[:, i])
            W[:, i] = weights
            residuals[i] = residual
        except Exception as e:
            logger.warning(f"NNLS failed for cell {i}: {e}")
            W[:, i] = 0
            residuals[i] = np.inf
    
    weights_df = pd.DataFrame(W, index=signature_matrix.columns, columns=mutation_matrix.columns)
    reconstruction = H @ W
    reconstruction_df = pd.DataFrame(reconstruction, index=mutation_matrix.index, columns=mutation_matrix.columns)
    
    if verbose:
        logger.info(f"Fitted {n_sigs} signatures to {n_cells} cells")
    
    return {'weights': weights_df, 'residuals': residuals, 'reconstruction': reconstruction_df}


def evaluate_reconstruction(original_matrix, reconstructed_matrix, verbose=True):
    """Evaluate reconstruction quality."""
    if original_matrix.shape[1] == 0:
        logger.warning("Empty matrix - skipping evaluation")
        return {
            'frobenius_norm_original': 0,
            'frobenius_norm_reconstructed': 0,
            'frobenius_error': 0,
            'relative_frobenius_error': 0,
            'optimality_ratio': 1.0,
            'pearson_per_cell': np.array([]),
            'cosine_per_cell': np.array([]),
            'cell_frobenius_errors': np.array([]),
            'mean_pearson': np.nan,
            'median_pearson': np.nan,
            'std_pearson': np.nan,
            'mean_cosine': np.nan,
            'median_cosine': np.nan,
            'mean_cell_error': 0,
            'quality': 'N/A',
            'quality_correlation': 'N/A'
        }
    
    X = original_matrix.values
    X_recon = reconstructed_matrix.values
    
    frobenius_norm_original = np.linalg.norm(X, 'fro')
    frobenius_error = np.linalg.norm(X - X_recon, 'fro')
    relative_error = frobenius_error / frobenius_norm_original if frobenius_norm_original > 0 else 0
    
    U, sv, Vt = np.linalg.svd(X, full_matrices=False)
    _, sv_recon, _ = np.linalg.svd(X_recon, full_matrices=False)
    decomp_rank = np.sum(sv_recon > 1e-10)
    
    theoretical_min = np.sqrt(np.sum(sv[decomp_rank:]**2)) if decomp_rank < len(sv) else 0
    optimality = theoretical_min / frobenius_error if frobenius_error > 0 else 1.0
    
    pearson_corrs = []
    cosine_sims = []
    cell_errors = []
    
    for i in range(X.shape[1]):
        cell_errors.append(np.linalg.norm(X[:, i] - X_recon[:, i]))
        if X[:, i].sum() > 0 and X_recon[:, i].sum() > 0:
            r, _ = pearsonr(X[:, i], X_recon[:, i])
            pearson_corrs.append(r)
            norm_x = np.linalg.norm(X[:, i])
            norm_r = np.linalg.norm(X_recon[:, i])
            if norm_x > 0 and norm_r > 0:
                cosine_sims.append(np.dot(X[:, i], X_recon[:, i]) / (norm_x * norm_r))
            else:
                cosine_sims.append(np.nan)
        else:
            pearson_corrs.append(np.nan)
            cosine_sims.append(np.nan)
    
    pearson_corrs = np.array(pearson_corrs)
    cosine_sims = np.array(cosine_sims)
    
    quality = "EXCELLENT" if relative_error < 0.1 else "GOOD" if relative_error < 0.2 else "MODERATE" if relative_error < 0.3 else "POOR"
    quality_corr = "EXCELLENT" if np.nanmean(pearson_corrs) > 0.8 else "GOOD" if np.nanmean(pearson_corrs) > 0.6 else "MODERATE" if np.nanmean(pearson_corrs) > 0.4 else "POOR"
    
    if verbose:
        logger.info(f"Frobenius error: {frobenius_error:.2f}, Relative: {100*relative_error:.2f}%")
        logger.info(f"Mean Pearson: {np.nanmean(pearson_corrs):.4f}, Quality: {quality}")
    
    return {
        'frobenius_norm_original': frobenius_norm_original,
        'frobenius_norm_reconstructed': np.linalg.norm(X_recon, 'fro'),
        'frobenius_error': frobenius_error,
        'relative_frobenius_error': relative_error,
        'optimality_ratio': optimality,
        'pearson_per_cell': pearson_corrs,
        'cosine_per_cell': cosine_sims,
        'cell_frobenius_errors': np.array(cell_errors),
        'mean_pearson': np.nanmean(pearson_corrs),
        'median_pearson': np.nanmedian(pearson_corrs),
        'std_pearson': np.nanstd(pearson_corrs),
        'mean_cosine': np.nanmean(cosine_sims),
        'median_cosine': np.nanmedian(cosine_sims),
        'mean_cell_error': np.mean(cell_errors),
        'quality': quality,
        'quality_correlation': quality_corr
    }


def select_signatures_scree(mutation_matrix, signature_pool, core_sigs, candidates, output_dir, max_sigs=15):
    """Select signatures via scree plot elbow detection."""
    logger.info("Selecting signatures via scree plot...")
    
    if mutation_matrix.shape[1] == 0:
        logger.warning("Empty mutation matrix - returning core signatures only")
        available_core = [s for s in core_sigs if s in signature_pool.columns]
        return {
            'selected_signatures': available_core,
            'signature_matrix': signature_pool[available_core] if available_core else signature_pool.iloc[:, :1],
            'n_signatures': len(available_core),
            'scree_data': []
        }
    
    available = set(signature_pool.columns)
    core = [s for s in core_sigs if s in available]
    cands = [s for s in candidates if s in available and s not in core]
    
    X = mutation_matrix.values
    X_norm = np.linalg.norm(X, 'fro')
    
    if X_norm == 0:
        logger.warning("Mutation matrix has zero norm - returning core signatures")
        return {
            'selected_signatures': core,
            'signature_matrix': signature_pool[core] if core else signature_pool.iloc[:, :1],
            'n_signatures': len(core),
            'scree_data': []
        }
    
    scores = []
    for sig in cands:
        fit = fit_signatures_nnls(mutation_matrix, signature_pool[[sig]], verbose=False)
        res = X - fit['reconstruction'].values
        exp_var = 1 - (np.linalg.norm(res, 'fro') / X_norm)**2
        scores.append({'signature': sig, 'explained_variance': exp_var})
    
    scores = sorted(scores, key=lambda x: x['explained_variance'], reverse=True)
    ordered_cands = [s['signature'] for s in scores]
    all_sigs = core + ordered_cands
    
    scree_data = []
    for n in range(len(core), min(max_sigs + 1, len(all_sigs) + 1)):
        test_sigs = all_sigs[:n]
        fit = fit_signatures_nnls(mutation_matrix, signature_pool[test_sigs], verbose=False)
        ev = evaluate_reconstruction(mutation_matrix, fit['reconstruction'], verbose=False)
        scree_data.append({
            'n': n,
            'signatures': test_sigs.copy(),
            'error': ev['frobenius_error'],
            'rel_error': ev['relative_frobenius_error'],
            'exp_var': 1 - ev['relative_frobenius_error']**2
        })
    
    if len(scree_data) < 2:
        logger.warning("Not enough data points for scree plot - using core signatures")
        return {
            'selected_signatures': core,
            'signature_matrix': signature_pool[core] if core else signature_pool.iloc[:, :1],
            'n_signatures': len(core),
            'scree_data': scree_data
        }
    
    errors = np.array([s['error'] for s in scree_data])
    d2y = np.gradient(np.gradient(errors))
    elbow_idx = np.argmax(d2y)
    final_n = scree_data[elbow_idx]['n']
    selected = scree_data[elbow_idx]['signatures']
    
    logger.info(f"Selected {final_n} signatures via scree plot")
    
    return {
        'selected_signatures': selected,
        'signature_matrix': signature_pool[selected],
        'n_signatures': final_n,
        'scree_data': scree_data
    }


def plot_results(weights_df, evaluation, output_dir):
    """Generate visualization plots."""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    if weights_df.shape[1] == 0:
        logger.warning("No cells to plot - skipping visualization")
        return
    
    # Weights summary
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    mean_weights = weights_df.mean(axis=1).sort_values(ascending=False)
    axes[0].barh(range(len(mean_weights)), mean_weights.values)
    axes[0].set_yticks(range(len(mean_weights)))
    axes[0].set_yticklabels(mean_weights.index, fontsize=8)
    axes[0].set_xlabel('Mean Weight')
    axes[0].set_title('Signature Activity')
    axes[0].invert_yaxis()
    
    pearson = evaluation.get('pearson_per_cell', np.array([]))
    if len(pearson) > 0:
        valid_p = pearson[~np.isnan(pearson)]
        if len(valid_p) > 0:
            axes[1].hist(valid_p, bins=50, color='steelblue', alpha=0.7)
            axes[1].axvline(np.nanmean(pearson), color='red', linestyle='--', label=f'Mean: {np.nanmean(pearson):.3f}')
            axes[1].set_xlabel('Pearson Correlation')
            axes[1].set_ylabel('Count')
            axes[1].set_title('Reconstruction Quality')
            axes[1].legend()
        else:
            axes[1].text(0.5, 0.5, 'No valid correlations', ha='center', va='center', transform=axes[1].transAxes)
            axes[1].set_title('Reconstruction Quality')
    else:
        axes[1].text(0.5, 0.5, 'No data', ha='center', va='center', transform=axes[1].transAxes)
        axes[1].set_title('Reconstruction Quality')
    
    plt.tight_layout()
    plt.savefig(output_path / 'signature_analysis_summary.png', dpi=200)
    plt.close()
    
    logger.info(f"Saved plots to {output_path}")


def plot_signature_umaps(adata, sig_cols, output_dir):
    """Generate UMAPs for signatures."""
    output_path = Path(output_dir) / 'signature_UMAPs'
    output_path.mkdir(exist_ok=True, parents=True)
    
    if 'X_umap' not in adata.obsm:
        logger.warning("No UMAP found in AnnData - skipping signature UMAPs")
        return
    
    for sig in sig_cols:
        if sig not in adata.obs.columns:
            continue
            
        if adata.obs[sig].dtype == 'object':
            adata.obs[sig] = pd.to_numeric(adata.obs[sig], errors='coerce').fillna(0)
        
        vals = adata.obs[sig].values.astype(float)
        vmax = max(np.percentile(vals[vals > 0], 95) if (vals > 0).any() else 1, 0.1)
        
        fig, ax = plt.subplots(figsize=(8, 8))
        sc.pl.umap(adata, color=sig, size=5, frameon=False, cmap='plasma', ax=ax, show=False, vmax=vmax, title=sig)
        plt.tight_layout()
        fig.savefig(output_path / f'UMAP_{sig}.png', dpi=150)
        plt.close()
    
    logger.info(f"Saved {len(sig_cols)} UMAPs")


def run_signature_analysis(
    mutations_file, adata_path, cosmic_file, output_dir,
    callable_sites_file=None, use_scree=False, core_sigs=None,
    candidate_order=None, mut_threshold=0, max_sigs=15, hnscc_only=False
):
    """Run complete signature analysis."""
    logger.info("="*60)
    logger.info("SIGNATURE ANALYSIS PIPELINE")
    logger.info("="*60)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    figures_dir = output_path / 'figures'
    figures_dir.mkdir(exist_ok=True)
    
    # Process mutations
    matrix_file = output_path / "mutation_matrix_96contexts.txt"
    mut_matrix = process_mutations(mutations_file, matrix_file)
    
    # Load AnnData early so we can save it even if there are no mutations
    logger.info(f"Loading AnnData from {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Handle empty mutations case
    if mut_matrix.shape[1] == 0:
        logger.warning("No mutations found - creating output files with zero mutation counts")
        
        # Create empty weights file
        empty_weights = pd.DataFrame(index=[], columns=['Cell'])
        empty_weights.to_csv(output_path / "signature_weights_per_cell.txt", sep='\t', index=True)
        
        # Add zero mutation count to adata
        adata.obs['total_mutations'] = 0
        
        # Handle column types for HDF5
        for col in adata.obs.columns:
            if adata.obs[col].dtype not in ['float64', 'float32', 'int64', 'int32', 'bool']:
                adata.obs[col] = adata.obs[col].astype(str)
        
        final_path = output_path / "adata_final.h5ad"
        adata.write(final_path)
        logger.info(f"Saved: {final_path}")
        
        logger.info("="*60)
        logger.info("COMPLETE (no mutations to analyze)")
        logger.info("="*60)
        
        return {'weights': empty_weights, 'evaluation': None, 'adata_path': str(final_path)}
    
    # Apply mutation threshold filter
    if mut_threshold > 0:
        cells_keep = mut_matrix.sum(axis=0) >= mut_threshold
        n_before = mut_matrix.shape[1]
        mut_matrix = mut_matrix.loc[:, cells_keep]
        logger.info(f"After filtering (>={mut_threshold} mutations): {mut_matrix.shape[1]} cells (removed {n_before - mut_matrix.shape[1]})")
    
    # Handle callable sites
    if callable_sites_file and os.path.exists(callable_sites_file):
        try:
            callable_df = pd.read_csv(callable_sites_file, sep='\t')
            if 'CB' in callable_df.columns:
                callable_bcs = set(callable_df['CB'])
                missing = set(mut_matrix.columns) - callable_bcs
                if missing:
                    logger.info(f"Setting {len(missing)} cells not in callable sites to 0")
                    mut_matrix.loc[:, list(missing)] = 0
        except Exception as e:
            logger.warning(f"Could not process callable sites file: {e}")
    
    # Load COSMIC signatures
    if cosmic_file is None:
        raise ValueError("COSMIC signature file is required but not provided")
    
    cosmic_sigs = extract_cosmic_signatures(cosmic_file, output_dir, hnscc_only)
    
    # Signature selection
    if use_scree and core_sigs:
        # Handle SBS40 variants
        if 'SBS40a' in cosmic_sigs.columns and 'SBS40' in (core_sigs or []):
            core_sigs = [s if s != 'SBS40' else 'SBS40a' for s in core_sigs]
        if candidate_order is None:
            candidate_order = [s for s in cosmic_sigs.columns if s not in core_sigs]
        
        selection = select_signatures_scree(mut_matrix, cosmic_sigs, core_sigs, candidate_order, output_dir, max_sigs)
        final_sigs = selection['signature_matrix']
        logger.info(f"Selected signatures: {list(final_sigs.columns)}")
    else:
        final_sigs = cosmic_sigs
        selection = None
        logger.info(f"Using all {len(final_sigs.columns)} signatures")
    
    # Fit and evaluate
    fitting = fit_signatures_nnls(mut_matrix, final_sigs)
    evaluation = evaluate_reconstruction(mut_matrix, fitting['reconstruction'])
    
    # Plot results
    plot_results(fitting['weights'], evaluation, figures_dir)
    
    # Save weights
    fitting['weights'].to_csv(output_path / "signature_weights_per_cell.txt", sep='\t', float_format='%.6f')
    logger.info(f"Saved weights for {fitting['weights'].shape[1]} cells")
    
    # Save evaluation metrics
    eval_df = pd.DataFrame([{
        'frobenius_error': evaluation['frobenius_error'],
        'relative_error': evaluation['relative_frobenius_error'],
        'mean_pearson': evaluation['mean_pearson'],
        'median_pearson': evaluation['median_pearson'],
        'quality': evaluation['quality'],
        'n_cells': mut_matrix.shape[1],
        'n_signatures': final_sigs.shape[1]
    }])
    eval_df.to_csv(output_path / "reconstruction_evaluation.txt", sep='\t', index=False)
    
    # Add to adata
    muts_per_cell = mut_matrix.sum(axis=0)
    muts_reindex = muts_per_cell.reindex(adata.obs.index).fillna(0)
    adata.obs['total_mutations'] = muts_reindex.values.astype(int)
    
    for sig in fitting['weights'].index:
        sig_vals = fitting['weights'].loc[sig].reindex(adata.obs.index).fillna(0)
        adata.obs[sig] = sig_vals.values.astype(float)
    
    # Handle column types for HDF5 compatibility
    for col in adata.obs.columns:
        if adata.obs[col].dtype not in ['float64', 'float32', 'int64', 'int32', 'bool']:
            adata.obs[col] = adata.obs[col].astype(str)
    
    # Generate signature UMAPs
    sig_cols = [c for c in adata.obs.columns if c.startswith('SBS')]
    if sig_cols:
        plot_signature_umaps(adata, sig_cols, figures_dir)
    
    # Save final AnnData
    final_path = output_path / "adata_final.h5ad"
    adata.write(final_path)
    logger.info(f"Saved: {final_path}")
    
    logger.info("="*60)
    logger.info("COMPLETE")
    logger.info("="*60)
    logger.info(f"Cells with mutations: {mut_matrix.shape[1]}")
    logger.info(f"Signatures used: {list(final_sigs.columns)}")
    logger.info(f"Reconstruction quality: {evaluation['quality']}")
    
    return {'weights': fitting['weights'], 'evaluation': evaluation, 'adata_path': str(final_path)}


def run_from_snakemake():
    """Run from Snakemake context."""
    # Set up file logging if log path provided
    if hasattr(snakemake, 'log') and snakemake.log:
        log_path = snakemake.log[0]
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    run_signature_analysis(
        mutations_file=snakemake.input.mutations,
        adata_path=snakemake.input.adata,
        cosmic_file=snakemake.params.cosmic_file,
        output_dir=snakemake.params.output_dir,
        callable_sites_file=snakemake.input.callable_sites,
        use_scree=getattr(snakemake.params, 'use_scree_plot', False),
        core_sigs=getattr(snakemake.params, 'core_signatures', None),
        candidate_order=getattr(snakemake.params, 'candidate_order', None),
        mut_threshold=getattr(snakemake.params, 'mutation_threshold', 0),
        max_sigs=getattr(snakemake.params, 'max_signatures', 15),
        hnscc_only=getattr(snakemake.params, 'hnscc_only', False)
    )


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Signature analysis for single-cell mutations")
    parser.add_argument('--mutations', required=True, help='Input mutations TSV file')
    parser.add_argument('--adata', required=True, help='Input AnnData H5AD file')
    parser.add_argument('--cosmic-file', required=True, help='COSMIC signatures file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--callable-sites', help='Callable sites TSV file')
    parser.add_argument('--use-scree', action='store_true', help='Use scree plot for signature selection')
    parser.add_argument('--core-signatures', nargs='+', help='Core signatures to always include')
    parser.add_argument('--candidate-order', nargs='+', help='Order of candidate signatures')
    parser.add_argument('--mutation-threshold', type=int, default=0, help='Minimum mutations per cell')
    parser.add_argument('--max-signatures', type=int, default=15, help='Maximum signatures to test')
    parser.add_argument('--hnscc-only', action='store_true', help='Use HNSCC-specific signatures only')
    
    args = parser.parse_args()
    
    run_signature_analysis(
        mutations_file=args.mutations,
        adata_path=args.adata,
        cosmic_file=args.cosmic_file,
        output_dir=args.output_dir,
        callable_sites_file=args.callable_sites,
        use_scree=args.use_scree,
        core_sigs=args.core_signatures,
        candidate_order=args.candidate_order,
        mut_threshold=args.mutation_threshold,
        max_sigs=args.max_signatures,
        hnscc_only=args.hnscc_only
    )


if __name__ == '__main__':
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        main()
