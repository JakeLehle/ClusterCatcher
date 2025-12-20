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

### 6. Set up Kraken2 database (optional, for viral detection)

If you want to use the viral/microbial detection feature, you need to set up a Kraken2 database. This is a one-time setup.

#### Option A: Download pre-built database

```bash
# Create database directory
mkdir -p ~/kraken2_db

# Download a pre-built viral database (smallest option)
# See: https://benlangmead.github.io/aws-indexes/k2 for options
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20231009.tar.gz
tar -xzf k2_viral_20231009.tar.gz -C ~/kraken2_db/viral
```

#### Option B: Build custom database

For more control, build your own database:

```bash
# Install Kraken2 (if not already installed)
conda create -n kraken2_build -c bioconda kraken2 blast
conda activate kraken2_build

# Create database directory
mkdir -p ~/kraken2_db/custom_viral
cd ~/kraken2_db/custom_viral

# Download taxonomy
kraken2-build --download-taxonomy --db .

# Download viral library
kraken2-build --download-library viral --db .

# Or add custom sequences (e.g., human disease-associated viruses only)
# kraken2-build --add-to-library your_sequences.fasta --db .

# Build the database (adjust threads as needed)
kraken2-build --build --db . --threads 16

# Generate inspect file (required for organism names)
kraken2-inspect --db . > inspect.txt

# Clean up intermediate files (optional, saves space)
kraken2-build --clean --db .
```

#### Option C: Human disease-associated viruses only

For cancer research, you may want only human viruses known to cause disease:

```bash
# Create filtered database
mkdir -p ~/kraken2_db/human_viral

# Download NCBI viral sequences
datasets download virus genome taxon "human virus" --filename human_viruses.zip
unzip human_viruses.zip

# Add to Kraken2 database
kraken2-build --download-taxonomy --db ~/kraken2_db/human_viral
kraken2-build --add-to-library ncbi_dataset/data/genomic.fna --db ~/kraken2_db/human_viral
kraken2-build --build --db ~/kraken2_db/human_viral --threads 16
kraken2-inspect --db ~/kraken2_db/human_viral > ~/kraken2_db/human_viral/inspect.txt
```

#### Verify database setup

```bash
# Check database files exist
ls -la ~/kraken2_db/viral/
# Should contain: hash.k2d, opts.k2d, taxo.k2d, inspect.txt

# Test with a small query
echo ">test" > test.fa
echo "ATCGATCGATCGATCG" >> test.fa
kraken2 --db ~/kraken2_db/viral test.fa
rm test.fa
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

## Cancer Cell Detection

The pipeline uses a dual-model approach combining CytoTRACE2 (stemness scoring) and inferCNV (copy number variation) to identify cancer cells with high confidence.

### How it works

1. **CytoTRACE2**: Scores each cell's developmental potential/stemness (0-1 scale)
2. **Threshold Detection**: Automatically detects the bimodal threshold separating normal and cancer populations
3. **InferCNV**: Uses normal cells as reference to detect copy number variations
4. **Agreement Scoring**: Calculates weighted agreement between both models
5. **Final Classification**: Cells are classified as cancer only if both models agree

### CytoTRACE2 Setup

CytoTRACE2 requires manual installation from GitHub:

```bash
# Clone the repository
git clone https://github.com/digitalcytometry/cytotrace2.git
cd cytotrace2/cytotrace2_python

# Install in your conda environment
conda activate ClusterCatcher
pip install .
```

### GTF File Requirement

InferCNV requires a GTF annotation file to map genes to chromosomal positions. You should provide the same GTF file used for Cell Ranger reference:

```bash
# For Cell Ranger GRCh38 reference
gunzip refdata-gex-GRCh38-2020-A/genes/genes.gtf.gz

# Then specify in config
ClusterCatcher create-config \
    --samples samples.pkl \
    --gtf /path/to/genes.gtf \
    ...
```

### Configuration Options

```yaml
dysregulation:
  enabled: true
  cytotrace2:
    species: human          # human or mouse
    max_cells_per_chunk: 200000  # Reduce for memory issues
    seed: 42
  infercnv:
    window_size: 250        # Genomic window size
  agreement:
    alpha: 0.5              # 0.5 = equal weight to rank and value agreement
    min_correlation: 0.5    # Minimum Spearman rho for quartile selection
```

### Outputs

The cancer detection step produces:
- `adata_cancer_detected.h5ad`: AnnData with all scores and classifications
- `cancer_detection_summary.tsv`: Summary statistics
- `figures/`: Comprehensive visualization plots including:
  - Score distribution histograms
  - Threshold detection plots
  - Chromosome heatmaps
  - Agreement analysis
  - Final UMAP visualizations

## Viral Detection

The viral detection module uses Kraken2 to identify viral/microbial sequences in unmapped reads from Cell Ranger. This creates a sparse matrix similar to Cell Ranger's gene expression matrix, but with organisms instead of genes.

### How it works

1. **Extract unmapped reads**: Reads that didn't map to the human genome are extracted from Cell Ranger's BAM file
2. **Kraken2 classification**: Unmapped reads are classified against a Kraken2 database
3. **Single-cell assignment**: Using cell barcodes and UMIs from the original BAM, each organism detection is linked back to specific cells
4. **Matrix generation**: A sparse matrix (cells × organisms) is created for integration with downstream analysis
5. **Viral integration**: Optionally filter to human-specific viruses and integrate with gene expression

### Two Kraken2 Databases

The pipeline can use two Kraken2 databases:

1. **Primary Database** (`--kraken2-db`): Used for classification. Can be broad (all viruses) or specific.
2. **Human Viral Database** (`--human-viral-db`): Optional. Filters results to human-associated viruses only.

### Setting up the Human Viral Database

```bash
# Create human-specific viral database
mkdir -p ~/kraken2_db/human_viral
kraken2-build --download-taxonomy --db ~/kraken2_db/human_viral

# Download human viral sequences
datasets download virus genome taxon "human virus" --filename human_viruses.zip
unzip human_viruses.zip

# Build database
kraken2-build --add-to-library ncbi_dataset/data/genomic.fna --db ~/kraken2_db/human_viral
kraken2-build --build --db ~/kraken2_db/human_viral --threads 16

# Generate inspect file (REQUIRED for viral integration)
kraken2-inspect --db ~/kraken2_db/human_viral > ~/kraken2_db/human_viral/inspect.txt
```

### Configuration options

```yaml
viral_detection:
  enabled: true
  kraken2_db: /path/to/kraken2_db          # Primary database
  human_viral_db: /path/to/inspect.txt     # For human virus filtering
  confidence: 0.0                           # Kraken2 confidence (0.0-1.0)
  include_organisms:                        # Only include these (optional)
    - "Human papillomavirus"
    - "Epstein-Barr"
  exclude_organisms:                        # Exclude these (optional)
    - "Bacteriophage"
  organisms_of_interest:                    # Highlight in reports
    - "Human papillomavirus"
    - "Epstein-Barr virus"
```

### Filtering organisms

You can filter the results in two ways:

1. **At detection time**: Use `--include-organisms` or `--exclude-organisms` to filter during Kraken2 processing
2. **At integration time**: Provide `--human-viral-db` to filter to human-specific viruses

Example:
```bash
ClusterCatcher create-config \
    --samples samples.pkl \
    --kraken2-db ~/kraken2_db/viral \
    --human-viral-db ~/kraken2_db/human_viral/inspect.txt \
    --output config.yaml
```

### Viral Integration Outputs

- `adata_with_virus.h5ad`: Expression data with top virus counts added
- `adata_viral_integrated.h5ad`: Full integrated data (genes + viruses)
- `viral_integration_summary.tsv`: Summary statistics
- `virus_scores.tsv`: Aggregated scores by cell type
- `figures/`: Matrix plots, UMAPs, violin plots

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
