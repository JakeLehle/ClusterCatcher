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
Four main commands:
1. inspect: Build a cellranger multi library sheet from SRAscraper metadata
2. sample-information: Process sample CSV and create sample dictionary pickle
3. create-config: Generate master config YAML for pipeline
4. run-config: Execute the Snakemake pipeline
"""
import click
import os
import sys
# Import subcommands
from cli.sample_information import sample_information
from cli.run_config import run_config
from cli.create_config import create_config_cmd
from cli.inspect_libraries import inspect_libraries
@click.group()
@click.version_option(version='1.4.2', prog_name='ClusterCatcher')
def main():
    """
    ClusterCatcher: Single-cell sequencing analysis pipeline
    
    A comprehensive pipeline for analyzing single-cell RNA sequencing data,
    detecting mutational signatures at single-cell resolution, and identifying
    dysregulated cells through multiple complementary approaches.
    
    \b
    Typical workflow:
    1. ClusterCatcher sample-information --input samples.csv --output samples.pkl
    2. ClusterCatcher create-config --sample-pickle samples.pkl --output-dir ./results [options]
    3. ClusterCatcher run-config ./results/config.yaml
    
    \b
    For SRAscraper users:
    - Skip step 1, use the metadata pickle from SRAscraper directly
    - Provide the SRAscraper pickle to create-config with --sample-pickle
    
    \b
    For paired 5' GEX + VDJ (TCR/BCR) data headed for `cellranger multi`:
    - First run `ClusterCatcher inspect` on the SRAscraper metadata to build and
      review a library sheet (patient grouping + feature types) before create-config.
    - ClusterCatcher inspect --metadata-dir /path/to/GSE/metadata \\
          --fastq-base /path/to/GSE/fastq --out ./library_sheet_draft.csv
    
    \b
    Quick start with existing data:
    ClusterCatcher create-config \\
        --output-dir ./results \\
        --sample-pickle samples.pkl \\
        --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A
    
    \b
    Full pipeline with all modules:
    ClusterCatcher create-config \\
        --output-dir ./results \\
        --sample-pickle samples.pkl \\
        --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \\
        --gtf-file /path/to/genes.gtf \\
        --enable-viral --kraken-db /path/to/kraken2_db \\
        --enable-scomatic \\
        --enable-signatures --cosmic-file /path/to/COSMIC_v3.4_SBS_GRCh38.txt
    """
    pass
# Register subcommands
main.add_command(inspect_libraries)
main.add_command(sample_information)
main.add_command(create_config_cmd)
main.add_command(run_config)
if __name__ == '__main__':
    main()
