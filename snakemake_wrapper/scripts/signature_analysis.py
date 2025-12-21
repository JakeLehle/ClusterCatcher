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


def process_mutations(input_file, output_file):
    """Process mutations into 96-trinucleotide context matrix."""
    logger.info(f"Processing mutations from {input_file}...")
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header = line.lstrip('#').strip().split('\t')
                break
    
    df = pd.read_csv(input_file, sep='\t', comment='#', names=header)
    df = df[(df['REF'].str.len() == 1) & (df['ALT_expected'].str.len() == 1)]
    df = df[~df['REF_TRI'].str.contains('N') & ~df['ALT_TRI'].str.contains('N')]
    
    for idx, row in df.iterrows():
        ref = row['REF']
        alt = row['ALT_expected']
        if ref in ['G', 'A']:
            df.at[idx, 'REF'] = reverse_complement(ref)
            df.at[idx, 'ALT_expected'] = reverse_complement(alt)
            df.at[idx, 'REF_TRI'] = reverse_complement(row['REF_TRI'])
            df.at[idx, 'ALT_TRI'] = reverse_complement(row['ALT_TRI'])
    
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
    
    all_mutation_types = [
        f"{five}[{ref}>{alt}]{three}"
        for ref in ['C', 'T']
        for alt in (['A', 'G', 'T'] if ref == 'C' else ['A', 'C', 'G'])
        for five in ['A', 'C', 'G', 'T']
        for three in ['A', 'C', 'G', 'T']
    ]
    
    cell_barcodes = list(mutation_counts.keys())
    data_dict = {mut_type: [mutation_counts[cb].get(mut_type, 0) for cb in cell_barcodes] for mut_type in all_mutation_types}
    result_df = pd.DataFrame(data_dict, index=cell_barcodes).T
    result_df = result_df.reindex(all_mutation_types, fill_value=0).astype(int)
    result_df.to_csv(output_file, sep='\t')
    
    logger.info(f"Matrix shape: {result_df.shape}")
    return result_df


def extract_cosmic_signatures(cosmic_file, output_dir=None, hnscc_only=False):
    """Extract COSMIC signatures."""
    logger.info(f"Loading COSMIC signatures from: {cosmic_file}")
    cosmic_all = pd.read_csv(cosmic_file, sep='\t', index_col=0)
    
    if hnscc_only:
        hnscc_sigs = ["SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS13", "SBS16", 
                      "SBS17a", "SBS17b", "SBS18", "SBS33", "SBS40a", "SBS40b", "SBS40c"]
        available = [s for s in hnscc_sigs if s in cosmic_all.columns]
        sigs = cosmic_all[available].copy()
    else:
        sigs = cosmic_all.copy()
    
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
    
    if not all(mutation_matrix.index == signature_matrix.index):
        signature_matrix = signature_matrix.loc[mutation_matrix.index]
    
    n_cells = mutation_matrix.shape[1]
    n_sigs = signature_matrix.shape[1]
    
    X = mutation_matrix.values
    H = signature_matrix.values
    W = np.zeros((n_sigs, n_cells))
    residuals = np.zeros(n_cells)
    
    for i in range(n_cells):
        try:
            weights, residual = nnls(H, X[:, i])
            W[:, i] = weights
            residuals[i] = residual
        except:
            W[:, i] = 0
            residuals[i] = np.inf
    
    weights_df = pd.DataFrame(W, index=signature_matrix.columns, columns=mutation_matrix.columns)
    reconstruction = H @ W
    reconstruction_df = pd.DataFrame(reconstruction, index=mutation_matrix.index, columns=mutation_matrix.columns)
    
    return {'weights': weights_df, 'residuals': residuals, 'reconstruction': reconstruction_df}


def evaluate_reconstruction(original_matrix, reconstructed_matrix, verbose=True):
    """Evaluate reconstruction quality."""
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
    
    available = set(signature_pool.columns)
    core = [s for s in core_sigs if s in available]
    cands = [s for s in candidates if s in available and s not in core]
    
    X = mutation_matrix.values
    X_norm = np.linalg.norm(X, 'fro')
    
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
    
    # Weights summary
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    mean_weights = weights_df.mean(axis=1).sort_values(ascending=False)
    axes[0].barh(range(len(mean_weights)), mean_weights.values)
    axes[0].set_yticks(range(len(mean_weights)))
    axes[0].set_yticklabels(mean_weights.index, fontsize=8)
    axes[0].set_xlabel('Mean Weight')
    axes[0].set_title('Signature Activity')
    axes[0].invert_yaxis()
    
    pearson = evaluation['pearson_per_cell']
    valid_p = pearson[~np.isnan(pearson)]
    axes[1].hist(valid_p, bins=50, color='steelblue', alpha=0.7)
    axes[1].axvline(np.nanmean(pearson), color='red', linestyle='--', label=f'Mean: {np.nanmean(pearson):.3f}')
    axes[1].set_xlabel('Pearson Correlation')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Reconstruction Quality')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(output_path / 'signature_analysis_summary.png', dpi=200)
    plt.close()
    
    logger.info(f"Saved plots to {output_path}")


def plot_signature_umaps(adata, sig_cols, output_dir):
    """Generate UMAPs for signatures."""
    output_path = Path(output_dir) / 'signature_UMAPs'
    output_path.mkdir(exist_ok=True, parents=True)
    
    for sig in sig_cols:
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
    
    if mut_threshold > 0:
        cells_keep = mut_matrix.sum(axis=0) >= mut_threshold
        mut_matrix = mut_matrix.loc[:, cells_keep]
        logger.info(f"After filtering: {mut_matrix.shape[1]} cells")
    
    # Load AnnData
    adata = sc.read_h5ad(adata_path)
    
    # Handle callable sites
    if callable_sites_file and os.path.exists(callable_sites_file):
        callable_df = pd.read_csv(callable_sites_file, sep='\t')
        callable_bcs = set(callable_df['CB'])
        missing = set(mut_matrix.columns) - callable_bcs
        if missing:
            mut_matrix.loc[:, list(missing)] = 0
    
    # Load signatures
    cosmic_sigs = extract_cosmic_signatures(cosmic_file, output_dir, hnscc_only)
    
    # Signature selection
    if use_scree and core_sigs:
        if 'SBS40a' in cosmic_sigs.columns and 'SBS40' in core_sigs:
            core_sigs = [s if s != 'SBS40' else 'SBS40a' for s in core_sigs]
        if candidate_order is None:
            candidate_order = [s for s in cosmic_sigs.columns if s not in core_sigs]
        
        selection = select_signatures_scree(mut_matrix, cosmic_sigs, core_sigs, candidate_order, output_dir, max_sigs)
        final_sigs = selection['signature_matrix']
    else:
        final_sigs = cosmic_sigs
        selection = None
    
    # Fit and evaluate
    fitting = fit_signatures_nnls(mut_matrix, final_sigs)
    evaluation = evaluate_reconstruction(mut_matrix, fitting['reconstruction'])
    
    # Plot
    plot_results(fitting['weights'], evaluation, figures_dir)
    
    # Save weights
    fitting['weights'].to_csv(output_path / "signature_weights_per_cell.txt", sep='\t', float_format='%.6f')
    
    # Add to adata
    muts_per_cell = mut_matrix.sum(axis=0)
    muts_reindex = muts_per_cell.reindex(adata.obs.index).fillna(0)
    adata.obs['total_mutations'] = muts_reindex.values
    
    for sig in fitting['weights'].index:
        sig_vals = fitting['weights'].loc[sig].reindex(adata.obs.index).fillna(0)
        adata.obs[sig] = sig_vals.values
    
    # Handle column types for HDF5
    for col in adata.obs.columns:
        if adata.obs[col].dtype not in ['float64', 'float32', 'int64', 'int32', 'bool']:
            adata.obs[col] = adata.obs[col].astype(str)
    
    # Signature UMAPs
    sig_cols = [c for c in adata.obs.columns if c.startswith('SBS')]
    if sig_cols:
        plot_signature_umaps(adata, sig_cols, figures_dir)
    
    # Save
    final_path = output_path / "adata_final.h5ad"
    adata.write(final_path)
    logger.info(f"Saved: {final_path}")
    
    logger.info("="*60)
    logger.info("COMPLETE")
    logger.info("="*60)
    
    return {'weights': fitting['weights'], 'evaluation': evaluation, 'adata_path': str(final_path)}


def run_from_snakemake():
    """Run from Snakemake context."""
    run_signature_analysis(
        mutations_file=snakemake.input.mutations,
        adata_path=snakemake.input.adata,
        cosmic_file=snakemake.params.cosmic_file,
        output_dir=snakemake.params.output_dir,
        callable_sites_file=snakemake.input.get('callable_sites'),
        use_scree=snakemake.params.get('use_scree_plot', False),
        core_sigs=snakemake.params.get('core_signatures'),
        candidate_order=snakemake.params.get('candidate_order'),
        mut_threshold=snakemake.params.get('mutation_threshold', 0),
        max_sigs=snakemake.params.get('max_signatures', 15),
        hnscc_only=snakemake.params.get('hnscc_only', False)
    )


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Signature analysis")
    parser.add_argument('--mutations', required=True)
    parser.add_argument('--adata', required=True)
    parser.add_argument('--cosmic-file', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--callable-sites')
    parser.add_argument('--use-scree', action='store_true')
    parser.add_argument('--core-signatures', nargs='+')
    parser.add_argument('--candidate-order', nargs='+')
    parser.add_argument('--mutation-threshold', type=int, default=0)
    parser.add_argument('--max-signatures', type=int, default=15)
    parser.add_argument('--hnscc-only', action='store_true')
    
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
