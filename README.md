# ClusterCatcher

Single-cell sequencing analysis pipeline for mutation signature detection and cell annotation.

## Overview

ClusterCatcher is a comprehensive Snakemake-based pipeline for analyzing single-cell RNA sequencing data. It integrates multiple analysis steps:

1. **Cell Ranger alignment and counting** - Align FASTQ files to reference transcriptome
2. **QC and filtering** - Scanpy-based quality control with automatic doublet removal
3. **Cell type annotation** - Consensus annotation using popV and/or CellTypist
4. **Dysregulation detection** - CytoTRACE2 (stemness) + inferCNV (chromosomal instability)
5. **Viral detection** - Kraken2-based detection of viral sequences in unmapped reads
6. **Somatic mutation calling** - SComatic for single-cell mutation detection
7. **Mutational signature deconvolution** - Semi-supervised NMF with COSMIC signatures

## Requirements

### External Software (must be installed separately)

- **Cell Ranger** (v7.0+): Required for alignment and counting
  - Download from: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
  - Must be added to PATH

### Conda/Mamba

- Conda or Mamba package manager
- All other dependencies are managed via conda environments

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/JakeLehle/ClusterCatcher.git
cd ClusterCatcher
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate ClusterCatcher
```

### 3. Install ClusterCatcher

```bash
pip install .
```

### 4. Install Cell Ranger (if not already installed)

Follow instructions at: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

Ensure `cellranger` is in your PATH:
```bash
cellranger --version
```

### 5. Download reference data

Download the appropriate Cell Ranger reference:
```bash
# For human (GRCh38)
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz
```

## Quick Start

### For new users (with their own FASTQ files)

```bash
# Step 1: Create sample information pickle from CSV
ClusterCatcher sample-information \
    --input samples.csv \
    --output samples.pkl

# Step 2: Generate pipeline configuration
ClusterCatcher create-config \
    --samples samples.pkl \
    --output config.yaml \
    --transcriptome /path/to/refdata-gex-GRCh38-2020-A \
    --reference-fasta /path/to/GRCh38.fa \
    --threads 16 \
    --memory 64

# Step 3: Run the pipeline
ClusterCatcher run-config config.yaml --cores 32
```

### For SRAscraper users

If you used SRAscraper to download data, you can use its pickle file directly:

```bash
# Use SRAscraper's dictionary_file.pkl directly
ClusterCatcher create-config \
    --samples /path/to/SRAscraper_output/metadata/dictionary_file.pkl \
    --output config.yaml \
    --transcriptome /path/to/refdata-gex-GRCh38-2020-A \
    --reference-fasta /path/to/GRCh38.fa

ClusterCatcher run-config config.yaml --cores 32
```

## Sample CSV Format

For the `sample-information` command, prepare a CSV file with the following format:

| sample_id | fastq_r1 | fastq_r2 | condition | donor |
|-----------|----------|----------|-----------|-------|
| S1 | /path/S1_R1.fq.gz | /path/S1_R2.fq.gz | tumor | P001 |
| S2 | /path/S2_R1.fq.gz | /path/S2_R2.fq.gz | normal | P001 |

For multi-lane samples, use comma-separated paths:
```
sample_id,fastq_r1,fastq_r2
S1,/path/S1_L1_R1.fq.gz,/path/S1_L2_R1.fq.gz,/path/S1_L1_R2.fq.gz,/path/S1_L2_R2.fq.gz
```

## CLI Commands

### `ClusterCatcher sample-information`

Convert a sample CSV file to a pickle dictionary.

```bash
ClusterCatcher sample-information --input samples.csv --output samples.pkl
```

Options:
- `--input, -i`: Path to CSV file (required)
- `--output, -o`: Path for output pickle (required)
- `--skip-validation, -s`: Skip FASTQ file existence checks
- `--verbose, -v`: Print verbose output

### `ClusterCatcher create-config`

Generate the master configuration YAML file.

```bash
ClusterCatcher create-config \
    --samples samples.pkl \
    --output config.yaml \
    --transcriptome /path/to/reference \
    [options]
```

Key options:
- `--samples, -s`: Path to sample pickle (required)
- `--output, -o`: Path for config YAML (required)
- `--transcriptome`: Path to Cell Ranger reference (required)
- `--reference-fasta`: Path to genome FASTA (for mutation calling)
- `--chemistry`: 10X chemistry (default: auto)
- `--threads, -t`: Threads per job (default: 8)
- `--memory, -m`: Memory in GB (default: 32)
- `--no-viral-detection`: Disable viral detection
- `--no-scomatic`: Disable mutation calling
- `--no-signatures`: Disable signature analysis

### `ClusterCatcher run-config`

Execute the Snakemake pipeline.

```bash
ClusterCatcher run-config config.yaml --cores 32
```

Options:
- `--cores, -c`: Total cores for Snakemake (default: 1)
- `--jobs, -j`: Concurrent jobs (default: 1)
- `--dryrun, -n`: Show what would be done
- `--profile`: Path to Snakemake profile for cluster execution
- `--cluster`: Cluster submission command

## Cluster Execution

### SLURM example

```bash
ClusterCatcher run-config config.yaml \
    --cores 100 \
    --jobs 10 \
    --cluster "sbatch --partition=normal --nodes=1 --ntasks-per-node {threads} --time 04:00:00"
```

### Using a Snakemake profile

```bash
ClusterCatcher run-config config.yaml --profile /path/to/slurm_profile
```

## Output Structure

```
results/
├── cellranger/           # Cell Ranger output per sample
│   └── {sample}/
│       └── outs/
│           ├── filtered_feature_bc_matrix.h5
│           ├── possorted_genome_bam.bam
│           └── web_summary.html
├── qc/                   # Quality control results
│   ├── qc_filtered.h5ad
│   ├── qc_metrics.tsv
│   ├── plots/
│   └── multiqc_report.html
├── annotation/           # Cell type annotations
│   ├── cell_annotations.h5ad
│   └── annotation_summary.tsv
├── dysregulation/        # Dysregulation scores
│   ├── cytotrace2_scores.tsv
│   ├── infercnv_scores.tsv
│   └── dysregulation_summary.tsv
├── viral/                # Viral detection results
│   └── viral_detection_summary.tsv
├── mutations/            # Somatic mutations
│   └── somatic_mutations.vcf.gz
├── signatures/           # Mutational signatures
│   └── signature_weights.tsv
└── master_summary.yaml   # Pipeline summary
```

## Chemistry Options

Cell Ranger will attempt automatic chemistry detection. If needed, you can specify:

- `auto` (default): Automatic detection
- `threeprime`: 3' gene expression
- `fiveprime`: 5' gene expression
- `SC3Pv2`, `SC3Pv3`, `SC3Pv3HT`, `SC3Pv4`: Specific 3' versions
- `SC5P-PE`, `SC5P-R2`: 5' paired-end or R2-only

The pipeline will automatically try multiple chemistry options if the specified one fails.

## Configuration Reference

See the generated `config.yaml` for all available parameters. Key sections:

- `cellranger`: Cell Ranger settings (chemistry, expect_cells, etc.)
- `qc`: Filtering thresholds (min_genes, max_pct_mito, etc.)
- `annotation`: Cell type annotation settings
- `dysregulation`: CytoTRACE2 and inferCNV settings
- `viral_detection`: Kraken2 settings
- `scomatic`: Mutation calling settings
- `signatures`: Signature deconvolution settings

## Troubleshooting

### Cell Ranger not found
Ensure Cell Ranger is installed and in your PATH:
```bash
export PATH=/path/to/cellranger:$PATH
```

### Chemistry detection fails
Try specifying a chemistry explicitly:
```bash
ClusterCatcher create-config ... --chemistry SC3Pv3
```

### Memory issues
Increase memory allocation:
```bash
ClusterCatcher create-config ... --memory 128
```

## Citation

Please cite ClusterCatcher if you use it in your research:

> Lehle, J. (2025). ClusterCatcher: Single-cell sequencing analysis pipeline for mutation signature detection. GitHub. https://github.com/JakeLehle/ClusterCatcher

## License

GPL-3.0

## Contact

For issues or feature requests, please open an issue on GitHub.
