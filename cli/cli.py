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
- Mutational signature deconvolution (semi-supervised NMF)

Three main commands:
1. sample-information: Process sample CSV and create sample dictionary pickle
2. create-config: Generate master config YAML for pipeline
3. run-config: Execute the Snakemake pipeline

"""

import click
import os
import sys

from cli.sample_information import sample_information
from cli.create_config import create_config
from cli.run_config import run_config


@click.group()
@click.version_option(version='0.1.0', prog_name='ClusterCatcher')
def main():
    """
    ClusterCatcher: Single-cell sequencing analysis pipeline
    
    A comprehensive pipeline for analyzing single-cell RNA sequencing data,
    detecting mutational signatures at single-cell resolution, and identifying
    dysregulated cells through multiple complementary approaches.
    
    \b
    Typical workflow:
    1. ClusterCatcher sample-information --input samples.csv --output samples.pkl
    2. ClusterCatcher create-config --samples samples.pkl --output config.yaml [options]
    3. ClusterCatcher run-config config.yaml
    
    \b
    For SRAscraper users:
    - Skip step 1, use the metadata pickle from SRAscraper directly
    - Provide the SRAscraper pickle to create-config with --samples
    """
    pass


# Register subcommands
main.add_command(sample_information)
main.add_command(create_config)
main.add_command(run_config)


if __name__ == '__main__':
    main()
