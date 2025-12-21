#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_summary.py
===================

Generate master summary YAML file combining all pipeline results.

This script aggregates metrics from all pipeline steps into a single
comprehensive summary file.

Author: Jake Lehle
Date: 2025
"""

import os
import sys
import yaml
import pandas as pd
from datetime import datetime
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def safe_read_tsv(filepath):
    """Safely read a TSV file, returning empty dict if file doesn't exist."""
    if filepath and os.path.exists(filepath):
        try:
            df = pd.read_csv(filepath, sep='\t')
            return df.to_dict(orient='records')[0] if len(df) > 0 else {}
        except Exception as e:
            logger.warning(f"Could not read {filepath}: {e}")
    return {}


def count_file_lines(filepath, skip_header=True):
    """Count lines in a file."""
    if not filepath or not os.path.exists(filepath):
        return 0
    try:
        with open(filepath, 'r') as f:
            lines = sum(1 for _ in f)
        return lines - 1 if skip_header else lines
    except:
        return 0


def generate_summary(
    qc_file,
    annotation_file,
    dysregulation_file,
    viral_file,
    mutations_file,
    signatures_file,
    config,
    sample_ids,
    output_dir,
    output_file
):
    """Generate comprehensive pipeline summary."""
    
    logger.info("="*60)
    logger.info("GENERATING MASTER SUMMARY")
    logger.info("="*60)
    
    summary = {
        'pipeline': {
            'name': 'ClusterCatcher',
            'version': '1.0.0',
            'completion_time': datetime.now().isoformat(),
        },
        'configuration': {
            'output_dir': output_dir,
            'n_samples': len(sample_ids),
            'sample_ids': sample_ids,
            'threads': config.get('threads', 8),
            'memory_gb': config.get('memory_gb', 64),
        },
        'modules': {
            'cellranger': True,
            'qc_annotation': True,
            'dysregulation': config.get('modules', {}).get('dysregulation', True),
            'viral_detection': config.get('modules', {}).get('viral', False),
            'scomatic': config.get('modules', {}).get('scomatic', False),
            'signatures': config.get('modules', {}).get('signatures', False),
        },
        'results': {}
    }
    
    # QC metrics
    qc_data = safe_read_tsv(qc_file)
    if qc_data:
        summary['results']['qc'] = {
            'status': 'completed',
            'metrics': qc_data
        }
    
    # Annotation metrics
    annotation_data = safe_read_tsv(annotation_file)
    if annotation_data:
        summary['results']['annotation'] = {
            'status': 'completed',
            'metrics': annotation_data
        }
    
    # Dysregulation metrics
    if dysregulation_file:
        dysreg_data = safe_read_tsv(dysregulation_file)
        if dysreg_data:
            summary['results']['dysregulation'] = {
                'status': 'completed',
                'metrics': dysreg_data
            }
    
    # Viral detection metrics
    if viral_file:
        viral_data = safe_read_tsv(viral_file)
        if viral_data:
            summary['results']['viral_detection'] = {
                'status': 'completed',
                'metrics': viral_data
            }
    
    # Mutation calling metrics
    if mutations_file and os.path.exists(mutations_file):
        n_mutations = count_file_lines(mutations_file)
        summary['results']['mutations'] = {
            'status': 'completed',
            'n_mutations': n_mutations,
            'output_file': mutations_file
        }
    
    # Signature analysis metrics
    if signatures_file and os.path.exists(signatures_file):
        try:
            sig_df = pd.read_csv(signatures_file, sep='\t', index_col=0)
            n_signatures = sig_df.shape[0]
            n_cells = sig_df.shape[1]
            mean_weights = sig_df.mean(axis=1).to_dict()
            
            summary['results']['signatures'] = {
                'status': 'completed',
                'n_signatures': n_signatures,
                'n_cells_analyzed': n_cells,
                'mean_signature_weights': {k: float(v) for k, v in mean_weights.items()}
            }
        except Exception as e:
            logger.warning(f"Could not parse signature file: {e}")
            summary['results']['signatures'] = {
                'status': 'completed',
                'output_file': signatures_file
            }
    
    # Output files
    summary['output_files'] = {
        'master_summary': output_file,
        'qc_metrics': qc_file,
        'annotation_summary': annotation_file,
    }
    
    if dysregulation_file:
        summary['output_files']['dysregulation_summary'] = dysregulation_file
        summary['output_files']['adata_cancer_detected'] = os.path.join(
            output_dir, 'dysregulation', 'adata_cancer_detected.h5ad'
        )
    
    if viral_file:
        summary['output_files']['viral_summary'] = viral_file
    
    if mutations_file:
        summary['output_files']['mutations'] = mutations_file
        summary['output_files']['callable_sites'] = os.path.join(
            output_dir, 'mutations', 'CombinedCallableSites', 'complete_callable_sites.tsv'
        )
    
    if signatures_file:
        summary['output_files']['signature_weights'] = signatures_file
        summary['output_files']['adata_final'] = os.path.join(
            output_dir, 'signatures', 'adata_final.h5ad'
        )
    
    # Write summary
    with open(output_file, 'w') as f:
        yaml.dump(summary, f, default_flow_style=False, sort_keys=False)
    
    logger.info(f"Summary written to: {output_file}")
    
    return summary


def run_from_snakemake():
    """Run from Snakemake context."""
    # Get inputs - handle both string and list inputs
    qc_file = snakemake.input.qc
    annotation_file = snakemake.input.annotation
    
    # Optional inputs may be empty lists
    dysregulation_file = snakemake.input.get('dysregulation')
    if isinstance(dysregulation_file, list):
        dysregulation_file = dysregulation_file[0] if dysregulation_file else None
    
    viral_file = snakemake.input.get('viral')
    if isinstance(viral_file, list):
        viral_file = viral_file[0] if viral_file else None
    
    mutations_file = snakemake.input.get('mutations')
    if isinstance(mutations_file, list):
        mutations_file = mutations_file[0] if mutations_file else None
    
    signatures_file = snakemake.input.get('signatures')
    if isinstance(signatures_file, list):
        signatures_file = signatures_file[0] if signatures_file else None
    
    generate_summary(
        qc_file=qc_file,
        annotation_file=annotation_file,
        dysregulation_file=dysregulation_file,
        viral_file=viral_file,
        mutations_file=mutations_file,
        signatures_file=signatures_file,
        config=snakemake.params.config,
        sample_ids=snakemake.params.sample_ids,
        output_dir=snakemake.params.output_dir,
        output_file=snakemake.output.summary
    )


if __name__ == '__main__':
    try:
        snakemake
        run_from_snakemake()
    except NameError:
        # CLI usage
        import argparse
        parser = argparse.ArgumentParser(description='Generate pipeline summary')
        parser.add_argument('--qc', required=True)
        parser.add_argument('--annotation', required=True)
        parser.add_argument('--dysregulation')
        parser.add_argument('--viral')
        parser.add_argument('--mutations')
        parser.add_argument('--signatures')
        parser.add_argument('--output-dir', required=True)
        parser.add_argument('--output', required=True)
        parser.add_argument('--sample-ids', nargs='+', default=[])
        
        args = parser.parse_args()
        
        generate_summary(
            qc_file=args.qc,
            annotation_file=args.annotation,
            dysregulation_file=args.dysregulation,
            viral_file=args.viral,
            mutations_file=args.mutations,
            signatures_file=args.signatures,
            config={},
            sample_ids=args.sample_ids,
            output_dir=args.output_dir,
            output_file=args.output
        )
