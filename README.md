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
7. **Mutational signature deconvolution** - Semi-supervised NNLS with COSMIC signatures

## Requirements

### External Software (must be installed separately)

- **Cell Ranger** (v7.0+): Required for alignment and counting
  - Download from: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
  - Must be added to PATH

- **SComatic** (for mutation calling): Required if using the somatic mutation module
  - Clone from: https://github.com/cortes-ciriano-lab/SComatic
  - See [SComatic Setup](#scomatic-setup) section below

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
python snakemake_wrapper/create_config.py \
    --output-dir ./results \
    --sample-pickle samples.pkl \
    --reference-fasta /path/to/GRCh38.fa \
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \
    --threads 16 \
    --memory-gb 64

# Step 3: Run the pipeline
cd snakemake_wrapper
snakemake --configfile ../results/config.yaml --cores 32 --use-conda
```

### For SRAscraper users

If you used SRAscraper to download data, you can use its pickle file directly:

```bash
python snakemake_wrapper/create_config.py \
    --output-dir ./results \
    --sample-pickle /path/to/SRAscraper_output/metadata/dictionary_file.pkl \
    --reference-fasta /path/to/GRCh38.fa \
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A

snakemake --configfile results/config.yaml --cores 32 --use-conda
```

### Enable all modules (mutations + signatures)

```bash
python snakemake_wrapper/create_config.py \
    --output-dir ./results \
    --sample-pickle samples.pkl \
    --reference-fasta /path/to/GRCh38.fa \
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \
    --enable-viral \
        --kraken-db /path/to/kraken2_db \
    --enable-scomatic \
        --scomatic-scripts-dir /path/to/SComatic/scripts \
        --scomatic-editing-sites /path/to/RNA_editing_sites.txt \
        --scomatic-pon-file /path/to/PoN.tsv \
        --scomatic-bed-file /path/to/mappable_regions.bed \
    --enable-signatures \
        --cosmic-file /path/to/COSMIC_v3.4_SBS_GRCh38.txt \
        --core-signatures SBS2 SBS13 SBS5 \
        --use-scree-plot
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
S1,/path/S1_L1_R1.fq.gz;/path/S1_L2_R1.fq.gz,/path/S1_L1_R2.fq.gz;/path/S1_L2_R2.fq.gz
```

## Pipeline Modules

### Module Overview

| Module | Required | Description |
|--------|----------|-------------|
| Cell Ranger | Yes | Alignment and counting |
| QC & Annotation | Yes | Quality control and cell typing |
| Dysregulation | Yes (default) | Cancer cell detection |
| Viral Detection | Optional | Kraken2-based pathogen detection |
| SComatic | Optional | Somatic mutation calling |
| Signatures | Optional | Mutational signature analysis |

### Enabling/Disabling Modules

```bash
# Enable viral detection
--enable-viral --kraken-db /path/to/db

# Enable mutation calling
--enable-scomatic --scomatic-scripts-dir /path/to/SComatic/scripts ...

# Enable signature analysis (requires SComatic)
--enable-signatures --cosmic-file /path/to/COSMIC_signatures.txt
```

## SComatic Setup

SComatic is required for somatic mutation calling. The pipeline uses a custom wrapper around SComatic scripts.

### 1. Clone SComatic

```bash
git clone https://github.com/cortes-ciriano-lab/SComatic.git
cd SComatic
```

### 2. Required Files

SComatic requires several reference files:

| File | Description | Source |
|------|-------------|--------|
| RNA editing sites | Known RNA editing positions | SComatic repository or RADAR database |
| Panel of Normals (PoN) | Common germline variants | Generate from matched normals or use provided |
| BED file | Mappable genomic regions | Generate or download for your genome build |

### 3. Generate reference files (if needed)

```bash
# Example: Generate mappable regions BED
# (Use appropriate tool like mappability or download pre-computed)
wget https://example.com/GRCh38_mappable_regions.bed.gz
gunzip GRCh38_mappable_regions.bed.gz
```

### 4. Pipeline configuration

```bash
python snakemake_wrapper/create_config.py \
    --enable-scomatic \
    --scomatic-scripts-dir /path/to/SComatic/scripts \
    --scomatic-editing-sites /path/to/RNA_editing_sites.txt \
    --scomatic-pon-file /path/to/PoN.tsv \
    --scomatic-bed-file /path/to/mappable_regions.bed \
    --scomatic-min-cov 5 \
    --scomatic-min-cells 5
```

### SComatic Pipeline Steps

The ClusterCatcher SComatic module performs a comprehensive 10-phase workflow:

1. **BAM Filtering** - Filter to annotated cells with valid barcodes
2. **Cell Type Splitting** - Split BAMs by cell type annotation
3. **Base Counting** - Count bases per cell (parallel processing)
4. **Count Merging** - Merge counts across cell types
5. **Variant Calling Step 1** - Initial variant detection
6. **Variant Calling Step 2** - Filter with PoN and editing sites
7. **BED Filtering** - Filter to mappable regions
8. **Callable Sites (Cell Type)** - Compute callable sites per cell type
9. **Callable Sites (Per Cell)** - Compute callable sites per individual cell
10. **Single-Cell Genotypes** - Generate filtered mutations with trinucleotide context

### SComatic Outputs

```
results/mutations/
├── all_samples.single_cell_genotype.filtered.tsv  # Final filtered mutations
├── CombinedCallableSites/
│   └── complete_callable_sites.tsv                # Per-cell callable sites
├── cell_annotations.tsv                           # Cell barcode to type mapping
├── trinucleotide_background.tsv                   # Background trinucleotide frequencies
└── {sample}/                                      # Per-sample intermediate files
    ├── SplitBams/
    ├── BaseCellCounts/
    └── VariantCalling/
```

## Signature Analysis

The signature analysis module uses semi-supervised Non-Negative Least Squares (NNLS) to decompose single-cell mutations into COSMIC mutational signatures.

### COSMIC Signatures

Download COSMIC signatures (v3.4 recommended):

```bash
# From COSMIC website (requires registration)
# https://cancer.sanger.ac.uk/signatures/downloads/

# Or use the GRCh38 version
wget https://cancer.sanger.ac.uk/signatures/documents/2123/COSMIC_v3.4_SBS_GRCh38.txt
```

### Configuration Options

```bash
python snakemake_wrapper/create_config.py \
    --enable-signatures \
    --cosmic-file /path/to/COSMIC_v3.4_SBS_GRCh38.txt \
    --core-signatures SBS2 SBS13 SBS5 \      # Always include these
    --candidate-order SBS1 SBS18 SBS40 \     # Try these in order
    --use-scree-plot \                        # Use elbow detection
    --mutation-threshold 0 \                  # Min mutations per cell
    --max-signatures 15 \                     # Max signatures to test
    --hnscc-only                              # Use HNSCC-specific signature set
```

### Signature Selection Methods

1. **All signatures** (default): Fit all COSMIC signatures
2. **HNSCC-specific** (`--hnscc-only`): Use curated set of 110 HNSCC-relevant signatures
3. **Scree plot** (`--use-scree-plot`): Automatic selection using elbow detection

### Signature Outputs

```
results/signatures/
├── signature_weights_per_cell.txt      # Per-cell signature weights matrix
├── adata_final.h5ad                    # Final AnnData with all annotations
├── mutation_matrix_96contexts.txt      # 96-trinucleotide context matrix
├── cosmic_signatures_used.txt          # Signatures used for fitting
├── reconstruction_evaluation.txt       # Quality metrics
├── per_cell_quality_metrics.txt        # Per-cell reconstruction quality
├── signature_weight_summary.txt        # Summary statistics
└── figures/
    ├── signature_analysis_summary.png
    ├── reconstruction_quality.png
    └── signature_UMAPs/                # Individual signature UMAPs
        ├── UMAP_SBS2.png
        ├── UMAP_SBS13.png
        └── ...
```

## Kraken2 Viral Detection Setup

### Option A: Download pre-built database

```bash
mkdir -p ~/kraken2_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20231009.tar.gz
tar -xzf k2_viral_20231009.tar.gz -C ~/kraken2_db/viral
```

### Option B: Build custom human viral database

```bash
mkdir -p ~/kraken2_db/human_viral
kraken2-build --download-taxonomy --db ~/kraken2_db/human_viral

# Download human viral sequences
datasets download virus genome taxon "human virus" --filename human_viruses.zip
unzip human_viruses.zip

kraken2-build --add-to-library ncbi_dataset/data/genomic.fna --db ~/kraken2_db/human_viral
kraken2-build --build --db ~/kraken2_db/human_viral --threads 16
kraken2-inspect --db ~/kraken2_db/human_viral > ~/kraken2_db/human_viral/inspect.txt
```

### Configuration

```bash
python snakemake_wrapper/create_config.py \
    --enable-viral \
    --kraken-db ~/kraken2_db/viral \
    --viral-db ~/kraken2_db/human_viral/inspect.txt \
    --viral-confidence 0.1
```

## Cancer Cell Detection

The pipeline uses a dual-model approach combining CytoTRACE2 (stemness scoring) and inferCNV (copy number variation).

### How it works

1. **CytoTRACE2**: Scores each cell's developmental potential/stemness (0-1 scale)
2. **Threshold Detection**: Automatically detects bimodal threshold
3. **InferCNV**: Uses normal cells as reference to detect CNVs
4. **Agreement Scoring**: Weighted agreement between models
5. **Final Classification**: Cells classified as cancer only if both models agree

### CytoTRACE2 Installation

```bash
git clone https://github.com/digitalcytometry/cytotrace2.git
cd cytotrace2/cytotrace2_python
pip install .
```

### Configuration

```bash
python snakemake_wrapper/create_config.py \
    --enable-dysregulation \
    --cytotrace2-enabled \
    --infercnv-enabled \
    --infercnv-reference-groups "T cells" "B cells" "Fibroblasts"
```

## Output Structure

```
results/
├── cellranger/                    # Cell Ranger output per sample
│   └── {sample}/outs/
│       ├── filtered_feature_bc_matrix.h5
│       ├── possorted_genome_bam.bam
│       └── web_summary.html
├── qc/                            # Quality control
│   ├── qc_metrics.tsv
│   └── multiqc_report.html
├── annotation/                    # Cell type annotations
│   ├── adata_annotated.h5ad
│   └── annotation_summary.tsv
├── dysregulation/                 # Cancer detection
│   ├── adata_cancer_detected.h5ad
│   ├── cancer_detection_summary.tsv
│   ├── dysregulation_summary.tsv
│   └── figures/
├── viral/                         # Viral detection (if enabled)
│   ├── viral_detection_summary.tsv
│   └── viral_counts.h5ad
├── viral_integration/             # Viral integration (if enabled)
│   ├── adata_with_virus.h5ad
│   └── viral_integration_summary.tsv
├── mutations/                     # SComatic output (if enabled)
│   ├── all_samples.single_cell_genotype.filtered.tsv
│   ├── CombinedCallableSites/
│   │   └── complete_callable_sites.tsv
│   └── cell_annotations.tsv
├── signatures/                    # Signature analysis (if enabled)
│   ├── signature_weights_per_cell.txt
│   ├── adata_final.h5ad
│   └── figures/
├── figures/                       # All visualization outputs
├── logs/                          # Pipeline logs
├── adata_final.h5ad              # Final annotated AnnData (copied from signatures/)
└── master_summary.yaml            # Pipeline summary
```

## Cluster Execution

### SLURM example

```bash
snakemake --configfile results/config.yaml \
    --cores 100 \
    --jobs 10 \
    --use-conda \
    --cluster "sbatch --partition=normal --nodes=1 --ntasks-per-node={threads} --time=04:00:00"
```

### Using a Snakemake profile

```bash
snakemake --configfile results/config.yaml --profile slurm_profile --use-conda
```

## Configuration Reference

### create_config.py Options

| Category | Option | Default | Description |
|----------|--------|---------|-------------|
| **Required** | `--output-dir` | - | Output directory |
| | `--reference-fasta` | - | Reference genome FASTA |
| | `--cellranger-reference` | - | Cell Ranger reference |
| **Samples** | `--sample-pickle` | - | Sample info pickle file |
| | `--sample-ids` | - | List of sample IDs |
| **Resources** | `--threads` | 8 | Threads per job |
| | `--memory-gb` | 64 | Memory in GB |
| **Cell Ranger** | `--chemistry` | auto | 10X chemistry |
| | `--expect-cells` | 10000 | Expected cells |
| **QC** | `--min-genes` | 200 | Min genes per cell |
| | `--min-counts` | 500 | Min counts per cell |
| | `--max-mito-pct` | 20 | Max mitochondrial % |
| **Viral** | `--enable-viral` | False | Enable viral detection |
| | `--kraken-db` | - | Kraken2 database path |
| **SComatic** | `--enable-scomatic` | False | Enable mutation calling |
| | `--scomatic-scripts-dir` | - | SComatic scripts path |
| | `--scomatic-editing-sites` | - | RNA editing sites file |
| | `--scomatic-pon-file` | - | Panel of Normals |
| | `--scomatic-bed-file` | - | Mappable regions BED |
| **Signatures** | `--enable-signatures` | False | Enable signature analysis |
| | `--cosmic-file` | - | COSMIC signatures file |
| | `--core-signatures` | SBS2,SBS13,SBS5 | Core signatures |
| | `--use-scree-plot` | False | Use elbow detection |
| | `--hnscc-only` | False | HNSCC signature set |

## Troubleshooting

### Cell Ranger not found
```bash
export PATH=/path/to/cellranger:$PATH
```

### Chemistry detection fails
```bash
python create_config.py ... --chemistry SC3Pv3
```

### Memory issues
```bash
python create_config.py ... --memory-gb 128
```

### SComatic script errors
Ensure SComatic is properly installed:
```bash
ls /path/to/SComatic/scripts/
# Should contain: SplitBamCellTypes.py, BaseCellCounter.py, etc.
```

### Signature analysis - no mutations detected
- Check mutation filtering thresholds
- Verify BAM files contain cell barcodes (CB tag)
- Check callable sites output

## Citation

Please cite ClusterCatcher if you use it in your research:

> Lehle, J. (2025). ClusterCatcher: Single-cell sequencing analysis pipeline for mutation signature detection. GitHub. https://github.com/JakeLehle/ClusterCatcher

If you use specific modules, please also cite:
- **Cell Ranger**: 10x Genomics
- **CytoTRACE2**: Gulati et al., Science 2020
- **SComatic**: Muyas et al., Nature Biotechnology 2024
- **COSMIC Signatures**: Alexandrov et al., Nature 2020

## License

GPL-3.0

## Contact

For issues or feature requests, please open an issue on GitHub.
