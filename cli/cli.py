#!/usr/bin/env python3
"""
ClusterCatcher CLI
==================

Single-cell sequencing analysis pipeline for:
- Cell Ranger alignment and counting
- Cell QC and annotation (Scanpy + popV)
- Dysregulated cell detection (CytoTRACE2 + inferCNV)
- Viral detection in unmapped reads (Kraken2)
- Somatic mutation calling (SComatic)
- Mutational signature deconvolution (semi-supervised NNLS)

Three main commands:
1. sample-information: Process sample CSV and create sample dictionary pickle
2. create-config: Generate master config YAML for pipeline
3. run-config: Execute the Snakemake pipeline
"""

import click
import os
import sys

# Import subcommands
from cli.sample_information import sample_information
from cli.run_config import run_config


@click.group()
@click.version_option(version='1.0.3', prog_name='ClusterCatcher')
def main():
    """
    ClusterCatcher: Single-cell sequencing analysis pipeline
    
    A comprehensive pipeline for analyzing single-cell RNA sequencing data,
    detecting mutational signatures at single-cell resolution, and identifying
    dysregulated cells through multiple complementary approaches.
    
    \b
    Typical workflow:
    1. ClusterCatcher sample-information --input samples.csv --output samples.pkl
    2. ClusterCatcher create-config --samples samples.pkl --output-dir ./results [options]
    3. ClusterCatcher run-config ./results/config.yaml
    
    \b
    For SRAscraper users:
    - Skip step 1, use the metadata pickle from SRAscraper directly
    - Provide the SRAscraper pickle to create-config with --sample-pickle
    
    \b
    Quick start with existing data:
    ClusterCatcher create-config \\
        --output-dir ./results \\
        --sample-pickle samples.pkl \\
        --reference-fasta /path/to/GRCh38.fa \\
        --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A
    
    \b
    Full pipeline with all modules:
    ClusterCatcher create-config \\
        --output-dir ./results \\
        --sample-pickle samples.pkl \\
        --reference-fasta /path/to/GRCh38.fa \\
        --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \\
        --gtf-file /path/to/genes.gtf \\
        --enable-viral --kraken-db /path/to/kraken2_db \\
        --enable-scomatic --scomatic-scripts-dir /path/to/SComatic/scripts \\
        --enable-signatures --cosmic-file /path/to/COSMIC_v3.4_SBS_GRCh38.txt
    """
    pass


@main.command('create-config')
@click.pass_context
def create_config_cmd(ctx):
    """
    Generate pipeline configuration file.
    
    This is a wrapper that calls the create_config.py script directly.
    For full options, run: python snakemake_wrapper/create_config.py --help
    """
    click.echo("For create-config, please run directly:")
    click.echo("  python snakemake_wrapper/create_config.py --help")
    click.echo("\nExample:")
    click.echo("  python snakemake_wrapper/create_config.py \\")
    click.echo("    --output-dir ./results \\")
    click.echo("    --sample-pickle samples.pkl \\")
    click.echo("    --reference-fasta /path/to/GRCh38.fa \\")
    click.echo("    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A")


# Register subcommands
main.add_command(sample_information)
main.add_command(run_config)


if __name__ == '__main__':
    main()
