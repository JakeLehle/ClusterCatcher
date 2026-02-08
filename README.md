# ClusterCatcher

**Single-cell RNA sequencing analysis pipeline for mutational signature detection, cell annotation, and cancer cell identification.**

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%203.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io)

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Quick Start](#quick-start)
6. [Detailed Usage](#detailed-usage)
7. [Pipeline Modules](#pipeline-modules)
8. [Configuration Reference](#configuration-reference)
9. [Output Structure](#output-structure)
10. [Troubleshooting](#troubleshooting)
11. [Citation](#citation)

---

## Overview

ClusterCatcher is a comprehensive Snakemake-based pipeline designed for end-to-end analysis of single-cell RNA sequencing data with a focus on:

- **Mutational signature detection** at single-cell resolution
- **Cancer/dysregulated cell identification** using dual-model consensus
- **Viral pathogen detection** in unmapped reads
- **Automated cell type annotation** using popV with cluster-based refinement

The pipeline integrates seven major analysis modules that can be flexibly enabled or disabled based on your research needs.

### Key Features

- **Modular design**: Enable only the modules you need
- **Reproducible**: Snakemake workflow with conda environment management
- **Scalable**: Supports local execution and HPC cluster submission
- **Comprehensive**: From raw FASTQs to annotated single-cell mutations
- **Well-documented**: Extensive logging and summary reports
- **Publication-ready plots**: Non-overlapping UMAP labels, stacked bar plots
- **Adaptive QC filtering**: MAD-based outlier detection that adapts to each dataset

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           ClusterCatcher Pipeline                           │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌──────────────────┐                                                       │
│  │   FASTQ Files    │                                                       │
│  └────────┬─────────┘                                                       │
│           │                                                                 │
│           ▼                                                                 │
│  ┌──────────────────┐     ┌──────────────────┐                              │
│  │  1. Cell Ranger  │────▶│  Unmapped Reads  │                              │
│  │  (Alignment)     │     └────────┬─────────┘                              │
│  └────────┬─────────┘              │                                        │
│           │                        ▼                                        │
│           │              ┌──────────────────┐                               │
│           │              │  5. Kraken2      │                               │
│           │              │  Viral Detection │                               │
│           │              └────────┬─────────┘                               │
│           ▼                       │                                         │
│  ┌──────────────────┐             │                                         │
│  │  2. QC/Filter    │             │                                         │
│  │  (Adaptive MAD)  │             │                                         │
│  └────────┬─────────┘             │                                         │
│           │                       │                                         │
│           ▼                       │                                         │
│  ┌──────────────────┐             │                                         │
│  │  3. popV         │             │                                         │
│  │  Annotation      │             │                                         │
│  └────────┬─────────┘             │                                         │
│           │                       │                                         │
│           ▼                       │                                         │
│  ┌──────────────────┐             │                                         │
│  │  4. Dysregulation│◀────────────┘                                         │
│  │  (CytoTRACE2 +   │                                                       │
│  │   inferCNV)      │                                                       │
│  └────────┬─────────┘                                                       │
│           │                                                                 │
│           ▼                                                                 │
│  ┌──────────────────┐                                                       │
│  │  6. SComatic     │                                                       │
│  │  Mutations       │                                                       │
│  └────────┬─────────┘                                                       │
│           │                                                                 │
│           ▼                                                                 │
│  ┌──────────────────┐                                                       │
│  │  7. Mutational   │                                                       │
│  │  Signatures      │                                                       │
│  │  (NNLS/COSMIC)   │                                                       │
│  └────────┬─────────┘                                                       │
│           │                                                                 │
│           ▼                                                                 │
│  ┌──────────────────┐                                                       │
│  │  Final AnnData   │                                                       │
│  │  + Summaries     │                                                       │
│  └──────────────────┘                                                       │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Pipeline Modules Summary

| # | Module | Script | Required | Description |
|---|--------|--------|----------|-------------|
| 1 | Cell Ranger | `cellranger_count.py` | Yes | FASTQ alignment, counting, BAM generation |
| 2 | QC & Filtering | `scanpy_qc_annotation.py` | Yes | Quality control, doublet removal, adaptive filtering |
| 3 | Cell Annotation | `scanpy_qc_annotation.py` | Yes | popV annotation with cluster-based refinement |
| 4 | Dysregulation | `cancer_cell_detection.py` | Default | Cancer cell detection (CytoTRACE2 + inferCNV) |
| 5 | Viral Detection | `kraken2_viral_detection.py` | Optional | Pathogen detection in unmapped reads |
| 6 | Mutation Calling | `scomatic_mutation_calling.py` | Optional | Somatic mutation calling (SComatic) |
| 7 | Signatures | `signature_analysis.py` | Optional | Mutational signature deconvolution |

---

## Requirements

### System Requirements

- **Linux** (tested on Ubuntu 20.04+, CentOS 7+)
- **Python** >= 3.9
- **Memory**: Minimum 64GB RAM recommended (128GB for large datasets)
- **Storage**: ~100GB per sample (including intermediates)

### External Software

These must be installed separately and available in PATH:

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **Cell Ranger** | >=7.0 | Alignment & counting | [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) |
| **Snakemake** | >=7.0 | Pipeline orchestration | `conda install -c bioconda snakemake` |
| **samtools** | >=1.15 | BAM processing | Installed via conda |
| **Kraken2** | >=2.1 | Viral detection (optional) | `conda install -c bioconda kraken2` |

### Optional External Dependencies

| Software | Purpose | Required for |
|----------|---------|--------------|
| **SComatic** | Mutation calling | `--enable-scomatic` |
| **CytoTRACE2** | Stemness scoring | `--enable-dysregulation` |

---

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/JakeLehle/ClusterCatcher.git
cd ClusterCatcher
```

### 2. Create Conda Environment

```bash
# Using mamba (recommended, faster)
mamba env create -f environment.yml

# Or using conda
conda env create -f environment.yml

# Activate the environment
conda activate ClusterCatcher
```

### 3. Install ClusterCatcher Package

```bash
pip install -e .
```

### 4. Install Cell Ranger

Download from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest):

```bash
# Extract Cell Ranger
tar -xzf cellranger-7.2.0.tar.gz

# Add to PATH
export PATH=/path/to/cellranger-7.2.0:$PATH

# Verify installation
cellranger --version
```

### 5. Download Reference Data

```bash
# Download Cell Ranger reference (GRCh38)
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz
```

**Reference directory structure** (standard 10x Genomics layout):

```
refdata-gex-GRCh38-2020-A/
├── fasta/
│   └── genome.fa          # Reference FASTA (auto-derived)
├── genes/
│   └── genes.gtf          # GTF annotation (auto-derived)
├── star/
│   └── ...                # STAR index
└── reference.json
```

> **Note**: ClusterCatcher only requires you to specify the Cell Ranger reference directory (`--cellranger-reference`). The FASTA and GTF paths are **automatically derived** from this directory using the standard 10x layout. You can override them with `--reference-fasta` or `--gtf-file` if using a non-standard layout.

### 6. (Optional) Install CytoTRACE2

Required for dysregulation detection:

```bash
git clone https://github.com/digitalcytometry/cytotrace2.git
cd cytotrace2/cytotrace2_python
pip install .
```

### 7. (Optional) Clone SComatic

Required for mutation calling:

```bash
git clone https://github.com/cortes-ciriano-lab/SComatic.git
```

---

## Quick Start

### Step 1: Prepare Sample Information

Create a CSV file with your sample information:

```csv
sample_id,fastq_r1,fastq_r2,condition,donor
SAMPLE1,/path/to/SAMPLE1_S1_L001_R1_001.fastq.gz,/path/to/SAMPLE1_S1_L001_R2_001.fastq.gz,tumor,P001
SAMPLE2,/path/to/SAMPLE2_S1_L001_R1_001.fastq.gz,/path/to/SAMPLE2_S1_L001_R2_001.fastq.gz,normal,P001
```

Convert to pickle format:

```bash
ClusterCatcher sample-information \
    --input samples.csv \
    --output samples.pkl
```

### Step 2: Generate Configuration

**Simplified usage** - only `--cellranger-reference` is required for references:

```bash
python snakemake_wrapper/create_config.py \
    --output-dir ./results \
    --sample-pickle samples.pkl \
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \
    --threads 16 \
    --memory-gb 64
```

The FASTA and GTF files are automatically derived from the Cell Ranger reference directory.

### Step 3: Run the Pipeline

```bash
cd snakemake_wrapper
snakemake --configfile ../results/config.yaml \
    --cores 16 \
    --use-conda
```

---

## Detailed Usage

### SRAscraper Integration

If using [SRAscraper](https://github.com/JakeLehle/SRAscraper) to download data:

```bash
# Use SRAscraper pickle directly - no conversion needed
python snakemake_wrapper/create_config.py \
    --output-dir ./results \
    --sample-pickle /path/to/srascraper/metadata/samples.pkl \
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A
```

### Sample CSV Format

#### Required Columns

| Column | Description |
|--------|-------------|
| `sample_id` | Unique identifier for the sample |
| `fastq_r1` | Path to R1 FASTQ file(s) |
| `fastq_r2` | Path to R2 FASTQ file(s) |

#### Optional Columns

| Column | Description |
|--------|-------------|
| `condition` | Experimental condition (tumor, normal, etc.) |
| `donor` | Donor/patient ID |
| `tissue` | Tissue type |
| `chemistry` | 10X chemistry version (SC3Pv2, SC3Pv3, etc.) |

#### Multi-Lane Samples

For samples with multiple sequencing lanes, use comma-separated paths:

```csv
sample_id,fastq_r1,fastq_r2
SAMPLE1,/path/L001_R1.fq.gz,/path/L002_R1.fq.gz,/path/L001_R2.fq.gz,/path/L002_R2.fq.gz
```

### Enable All Modules

```bash
python snakemake_wrapper/create_config.py \
    --output-dir ./results \
    --sample-pickle samples.pkl \
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \
    \
    # Viral detection
    --enable-viral \
    --kraken-db /path/to/kraken2_db \
    --viral-db /path/to/human_viral_db/inspect.txt \
    \
    # Dysregulation (enabled by default)
    --enable-dysregulation \
    --infercnv-reference-groups "T cells" "B cells" "Fibroblasts" \
    \
    # SComatic mutation calling
    --enable-scomatic \
    --scomatic-scripts-dir /path/to/SComatic/scripts \
    --scomatic-editing-sites /path/to/RNA_editing_sites.txt \
    --scomatic-pon-file /path/to/PoN.tsv \
    --scomatic-bed-file /path/to/mappable_regions.bed \
    \
    # Signature analysis
    --enable-signatures \
    --cosmic-file /path/to/COSMIC_v3.4_SBS_GRCh38.txt \
    --core-signatures SBS2 SBS13 SBS5 \
    --use-scree-plot \
    --hnscc-only
```

### Cluster Execution

#### SLURM Example

```bash
snakemake --configfile results/config.yaml \
    --cores 100 \
    --jobs 10 \
    --use-conda \
    --latency-wait 60 \
    --cluster "sbatch --partition=normal --nodes=1 --ntasks-per-node={threads} --mem={resources.mem_mb}M --time=08:00:00"
```

#### Using Snakemake Profiles

```bash
snakemake --configfile results/config.yaml \
    --profile /path/to/slurm_profile \
    --use-conda
```

---

## Pipeline Modules

### Module 1: Cell Ranger Alignment

**Script**: `scripts/cellranger_count.py`

Aligns FASTQ files to the reference transcriptome and generates gene expression matrices.

#### Inputs
- FASTQ files (R1 and R2)
- Cell Ranger reference transcriptome

#### Outputs
- `cellranger/{sample}/outs/filtered_feature_bc_matrix.h5` - Gene expression matrix
- `cellranger/{sample}/outs/possorted_genome_bam.bam` - Aligned reads
- `cellranger/{sample}/outs/web_summary.html` - QC report

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `chemistry` | `auto` | 10X chemistry version |
| `expect_cells` | `10000` | Expected number of cells |
| `include_introns` | `true` | Include intronic reads |
| `localcores` | `8` | CPU cores for Cell Ranger |
| `localmem` | `64` | Memory (GB) for Cell Ranger |

---

### Module 2: QC and Filtering

**Script**: `scripts/scanpy_qc_annotation.py`

Performs quality control, doublet detection, and cell filtering using Scanpy. Supports two filtering modes: **adaptive** (MAD-based, recommended) and **fixed** (threshold-based).

#### QC Metrics Calculated
- Number of genes per cell (`n_genes_by_counts`)
- Total UMI counts per cell (`total_counts`)
- Mitochondrial gene percentage (`pct_counts_mt`)
- Ribosomal gene percentage (`pct_counts_ribo`)
- Hemoglobin gene percentage (`pct_counts_hb`)
- Percent counts in top 20 genes (`pct_counts_in_top_20_genes`)
- Doublet scores (Scrublet)

#### QC Filtering Modes

ClusterCatcher supports two QC filtering approaches:

##### Adaptive Mode (Default, Recommended)

Uses **Median Absolute Deviation (MAD)** to identify statistical outliers. This approach:
- Adapts automatically to each dataset's distribution
- Preserves biological heterogeneity while removing technical artifacts
- Works well for heterogeneous samples (e.g., tumor tissues, mixed populations)
- Based on best practices from Luecken & Theis (2019) and Heumos et al. (2023)

**How it works:**
1. Calculates log1p-transformed metrics: `log1p_total_counts`, `log1p_n_genes_by_counts`
2. Computes `pct_counts_in_top_20_genes` to detect cells dominated by few genes
3. For each metric, calculates median and MAD
4. Flags cells where `|value - median| > nmads × MAD` as outliers
5. Combines outlier flags with OR logic for general QC and mitochondrial filtering

```
Outlier = (log1p_total_counts outlier) OR 
          (log1p_n_genes_by_counts outlier) OR 
          (pct_counts_in_top_20_genes outlier)

MT_Outlier = (pct_counts_mt outlier) OR (pct_counts_mt > max_mito_pct)

Final Filter = NOT(Outlier) AND NOT(MT_Outlier) AND (n_genes >= min_genes)
```

##### Fixed Mode (Original Behavior)

Uses hard thresholds for filtering. Good for:
- Well-characterized datasets with known quality ranges
- Reproducibility with specific threshold values
- Backward compatibility with existing workflows

**Filtering Steps (Fixed Mode):**
1. Remove cells with < `min_genes` genes detected
2. Remove cells with > `max_genes` genes detected (if set)
3. Remove cells with < `min_counts` total counts (if set)
4. Remove cells with > `max_counts` total counts (if set)
5. Remove cells with > `max_mito_pct` mitochondrial reads
6. Remove genes expressed in < `min_cells` cells
7. Remove predicted doublets

#### Key Parameters

##### Common Parameters (Both Modes)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `qc_mode` | `adaptive` | Filtering mode: "adaptive" or "fixed" |
| `min_genes` | `200` | Minimum genes per cell (always applied) |
| `max_mito_pct` | `20` | Maximum mitochondrial % (hard cap in adaptive, threshold in fixed) |
| `min_cells` | `20` | Minimum cells expressing a gene |
| `doublet_removal` | `true` | Enable Scrublet doublet detection |
| `doublet_rate` | `0.08` | Expected doublet rate for Scrublet |

##### Adaptive Mode Parameters (MAD-based)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mad_n_counts` | `5` | MADs for `log1p_total_counts` outliers |
| `mad_n_genes` | `5` | MADs for `log1p_n_genes_by_counts` outliers |
| `mad_top_genes` | `5` | MADs for `pct_counts_in_top_20_genes` outliers |
| `mad_mito` | `3` | MADs for `pct_counts_mt` outliers |

> **Note:** Higher MAD values = more permissive (fewer cells filtered). 5 MADs is permissive for counts/genes to preserve biological heterogeneity; 3 MADs is more stringent for mitochondrial content to remove dying cells.

##### Fixed Mode Parameters (Threshold-based)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_genes` | `null` | Maximum genes per cell |
| `min_counts` | `null` | Minimum UMI counts per cell |
| `max_counts` | `null` | Maximum UMI counts per cell |

> **Note:** Defaults are `null` to prevent over-filtering. When using fixed mode, set these explicitly based on your dataset's characteristics.

#### Outputs
- `qc/qc_metrics.tsv` - Per-sample QC statistics
- `qc/multiqc_report.html` - MultiQC summary
- `figures/qc/pre_filter_*.png` - Pre-filter QC plots
- `figures/qc/post_filter_*.png` - Post-filter QC plots

---

### Module 3: Cell Type Annotation (popV)

**Script**: `scripts/scanpy_qc_annotation.py` (combined with QC)

Annotates cells with cell type labels using popV with cluster-based refinement.

#### Workflow

The annotation follows a 7-step pipeline:

```
┌──────────────────┐
│ 1. Load Cell     │
│    Ranger Data   │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ 2. QC Metrics    │
│    Calculation   │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ 3. Cell/Gene     │
│    Filtering     │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ 4. Doublet       │
│    Detection     │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ 5. popV Annot.   │◄── RAW COUNTS (no normalization)
│  (HubModel)      │    Pulls from Hugging Face
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ 6. Post-Annot.   │◄── Normalize (CPM), UMAP, Leiden
│    Processing    │    Optional BBKNN batch correction
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ 7. Cluster-Based │◄── Weighted scoring per cluster
│    Refinement    │    Assigns final_annotation
└──────────────────┘
```

#### popV Annotation Details

ClusterCatcher uses **popV** (Population-level Voting) with HubModel from Hugging Face for automated cell type annotation.

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `popv_huggingface_repo` | `popV/tabula_sapiens_All_Cells` | Model repository |
| `popv_prediction_mode` | `inference` | Prediction mode |
| `popv_gene_symbol_key` | `feature_name` | Gene symbol column |
| `batch_key` | `sample_id` | Batch column for correction |

#### Outputs
- `annotation/adata_annotated.h5ad` - Final annotated AnnData
- `annotation/annotation_summary.tsv` - Cell type counts per sample
- `figures/qc/UMAP_*.pdf` - Annotation visualizations

---

### Module 4: Dysregulation Detection

**Script**: `scripts/cancer_cell_detection.py`

Identifies dysregulated/cancer cells using dual-model consensus (CytoTRACE2 + inferCNV).

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cytotrace2.enabled` | `true` | Enable CytoTRACE2 |
| `infercnv.enabled` | `true` | Enable inferCNV |
| `infercnv.window_size` | `250` | Window size for CNV |
| `infercnv_reference_groups` | `null` | Reference cell types |

#### Outputs
- `dysregulation/adata_cancer_detected.h5ad` - AnnData with cancer scores
- `dysregulation/dysregulation_summary.tsv` - Summary statistics

---

### Module 5: Viral Detection

**Script**: `scripts/kraken2_viral_detection.py`

Detects viral pathogens in unmapped reads using Kraken2.

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `kraken_db` | required | Kraken2 database path |
| `viral_db` | optional | Viral database for integration |
| `confidence` | `0.1` | Classification confidence |

#### Outputs
- `viral/viral_detection_summary.tsv` - Viral detection results
- `viral/viral_counts.h5ad` - Viral counts AnnData

---

### Module 6: SComatic Mutation Calling

**Script**: `scripts/scomatic_mutation_calling.py`

Calls somatic mutations at single-cell resolution using SComatic.

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `scripts_dir` | required | SComatic scripts directory |
| `editing_sites` | required | RNA editing sites file |
| `pon_file` | required | Panel of Normals |
| `bed_file` | required | Mappable regions BED |
| `min_cov` | `5` | Minimum coverage |
| `min_cells` | `5` | Minimum cells with variant |

#### Outputs
- `mutations/all_samples.single_cell_genotype.filtered.tsv` - Filtered mutations
- `mutations/CombinedCallableSites/` - Callable sites

---

### Module 7: Mutational Signature Analysis

**Script**: `scripts/signature_analysis.py`

Deconvolves single-cell mutations into COSMIC mutational signatures.

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cosmic_file` | required | COSMIC signatures file |
| `core_signatures` | `SBS2,SBS13,SBS5` | Always include |
| `use_scree_plot` | `true` | Use elbow detection |
| `hnscc_only` | `false` | HNSCC signature set |

#### Outputs
- `signatures/signature_weights_per_cell.txt` - Per-cell weights
- `signatures/adata_final.h5ad` - Final AnnData

---

## Configuration Reference

### Reference Files (Simplified)

The reference configuration has been **simplified** to reduce redundancy:

```yaml
reference:
  # Cell Ranger reference directory - THE ONLY REQUIRED FIELD
  # Must contain fasta/, genes/, star/ subdirectories
  cellranger: "/path/to/refdata-gex-GRCh38-2020-A"
  
  # Auto-derived from cellranger: {cellranger}/fasta/genome.fa
  # Override only if using non-standard layout
  fasta: null
  
  # Auto-derived from cellranger: {cellranger}/genes/genes.gtf
  # Override only if using non-standard layout
  gtf: null
  
  # Genome build (GRCh38, GRCh37, mm10, mm39)
  genome: "GRCh38"
```

### Complete Configuration Example

```yaml
# =============================================================================
# ClusterCatcher Configuration
# =============================================================================

output_dir: "/path/to/results"
threads: 8
memory_gb: 64

# Sample specification
sample_info: "/path/to/samples.pkl"

# Reference (only cellranger required)
reference:
  cellranger: "/path/to/refdata-gex-GRCh38-2020-A"
  genome: "GRCh38"

# QC parameters - Adaptive mode (recommended)
qc:
  qc_mode: "adaptive"           # "adaptive" (MAD-based) or "fixed" (thresholds)
  
  # Common parameters
  min_genes: 200                # Always applied
  max_mito_pct: 20              # Hard cap in adaptive, threshold in fixed
  min_cells: 20                 # Min cells expressing a gene
  doublet_removal: true
  doublet_rate: 0.08
  
  # Adaptive mode: MAD thresholds
  mad_n_counts: 5               # MADs for log1p_total_counts
  mad_n_genes: 5                # MADs for log1p_n_genes_by_counts
  mad_top_genes: 5              # MADs for pct_counts_in_top_20_genes
  mad_mito: 3                   # MADs for pct_counts_mt
  
  # Fixed mode: Hard thresholds (only used when qc_mode: "fixed")
  max_genes: null               # Set explicitly if using fixed mode
  min_counts: null              # Set explicitly if using fixed mode
  max_counts: null              # Set explicitly if using fixed mode

# Preprocessing (post-annotation)
preprocessing:
  target_sum: 1000000  # CPM
  n_pcs: null          # Uses CPU count
  n_neighbors: 15
  leiden_resolution: 1.0
  run_bbknn: false
  bbknn_batch_key: "sample_id"

# Annotation (popV)
annotation:
  popv_huggingface_repo: "popV/tabula_sapiens_All_Cells"
  popv_prediction_mode: "inference"
  popv_gene_symbol_key: "feature_name"
  popv_cache_dir: "tmp/popv_models"
  batch_key: "sample_id"

# Module flags
modules:
  cellranger: true
  qc: true
  annotation: true
  viral: false
  dysregulation: true
  scomatic: false
  signatures: false
```

### Complete `create_config.py` Options

#### Required Arguments

| Option | Description |
|--------|-------------|
| `--output-dir` | Output directory for results |
| `--cellranger-reference` | Cell Ranger reference directory |

#### Reference Files (Optional - Auto-derived)

| Option | Description |
|--------|-------------|
| `--reference-fasta` | Override auto-derived FASTA path |
| `--gtf-file` | Override auto-derived GTF path |
| `--genome` | Genome build (default: GRCh38) |

#### Sample Specification (one required)

| Option | Description |
|--------|-------------|
| `--sample-pickle` | Pickle file with sample dictionary |
| `--sample-ids` | Space-separated list of sample IDs |

#### Resource Settings

| Option | Default | Description |
|--------|---------|-------------|
| `--threads` | `8` | Threads per job |
| `--memory-gb` | `64` | Memory in GB |

#### Cell Ranger Settings

| Option | Default | Description |
|--------|---------|-------------|
| `--chemistry` | `auto` | 10X chemistry version |
| `--expect-cells` | `10000` | Expected cell count |
| `--include-introns` | `false` | Include intronic reads |

#### QC Settings

| Option | Default | Description |
|--------|---------|-------------|
| `--qc-mode` | `adaptive` | QC filtering mode: "adaptive" (MAD-based) or "fixed" (thresholds) |
| `--min-genes` | `200` | Min genes per cell (always applied) |
| `--max-mito-pct` | `20` | Max mitochondrial % |
| `--min-cells` | `20` | Min cells expressing gene |
| `--doublet-rate` | `0.08` | Expected doublet rate |

##### Adaptive Mode Options (MAD-based outlier detection)

| Option | Default | Description |
|--------|---------|-------------|
| `--mad-n-counts` | `5` | MADs for log1p_total_counts outliers |
| `--mad-n-genes` | `5` | MADs for log1p_n_genes_by_counts outliers |
| `--mad-top-genes` | `5` | MADs for pct_counts_in_top_20_genes outliers |
| `--mad-mito` | `3` | MADs for pct_counts_mt outliers |

##### Fixed Mode Options (threshold-based filtering)

| Option | Default | Description |
|--------|---------|-------------|
| `--max-genes` | `null` | Max genes per cell |
| `--min-counts` | `null` | Min UMI counts |
| `--max-counts` | `null` | Max UMI counts |

#### Module Enable/Disable

| Option | Default | Description |
|--------|---------|-------------|
| `--enable-viral` | `false` | Enable viral detection |
| `--enable-dysregulation` | `true` | Enable dysregulation |
| `--enable-scomatic` | `false` | Enable mutation calling |
| `--enable-signatures` | `false` | Enable signature analysis |

---

## Output Structure

```
results/
├── cellranger/                          # Cell Ranger outputs
│   └── {sample}/outs/
│       ├── filtered_feature_bc_matrix.h5
│       ├── possorted_genome_bam.bam
│       ├── possorted_genome_bam.bam.bai
│       └── web_summary.html
│
├── qc/                                  # Quality control
│   ├── qc_metrics.tsv
│   └── multiqc_report.html
│
├── figures/                             # Visualizations
│   └── qc/
│       ├── pre_filter_qc_violin_plots.pdf
│       ├── pre_filter_qc_scatter_plots.pdf
│       ├── post_filter_qc_violin_plots.pdf
│       ├── post_filter_qc_scatter_plots.pdf
│       ├── UMAP_popv_prediction.pdf
│       ├── UMAP_popv_prediction_score.pdf
│       ├── UMAP_final_annotation.pdf
│       ├── UMAP_samples.pdf
│       ├── UMAP_clusters.pdf
│       ├── Stacked_Bar_Cluster_Composition.pdf
│       └── cell_type_proportions.pdf
│
├── annotation/                          # Cell type annotation
│   ├── adata_annotated.h5ad
│   └── annotation_summary.tsv
│
├── dysregulation/                       # Cancer cell detection
│   ├── adata_cancer_detected.h5ad
│   ├── cancer_detection_summary.tsv
│   ├── dysregulation_summary.tsv
│   └── figures/
│
├── viral/                               # Viral detection (if enabled)
│   ├── {sample}/
│   │   ├── kraken2_filtered_feature_bc_matrix/
│   │   └── {sample}_organism_summary.tsv
│   ├── viral_detection_summary.tsv
│   └── viral_counts.h5ad
│
├── mutations/                           # SComatic output (if enabled)
│   ├── all_samples.single_cell_genotype.filtered.tsv
│   └── {sample}/
│
├── signatures/                          # Signature analysis (if enabled)
│   ├── signature_weights_per_cell.txt
│   ├── adata_final.h5ad
│   └── figures/
│
├── logs/                                # Pipeline logs
│
├── config.yaml                          # Pipeline configuration
└── master_summary.yaml                  # Pipeline summary
```

---

## Troubleshooting

### Common Issues

#### Cell Ranger not found

```bash
# Add Cell Ranger to PATH
export PATH=/path/to/cellranger-7.2.0:$PATH

# Verify
which cellranger
cellranger --version
```

#### Reference FASTA not found (auto-derivation failed)

If using a non-standard Cell Ranger reference layout:

```bash
python create_config.py ... \
    --reference-fasta /custom/path/to/genome.fa \
    --gtf-file /custom/path/to/genes.gtf
```

#### Chemistry detection fails

Specify chemistry explicitly:

```bash
python create_config.py ... --chemistry SC3Pv3
```

#### Memory issues

Increase memory allocation:

```bash
python create_config.py ... --memory-gb 128
```

#### popV model download fails

Check internet connection and try manually downloading:

```bash
# Clear cache and retry
rm -rf tmp/popv_models
# Run pipeline again
```

Or specify a different model:

```yaml
annotation:
  popv_huggingface_repo: "popV/tabula_sapiens_immune"
```

### QC Filtering Issues

#### Too many cells filtered

If adaptive mode is removing too many cells:

1. Check your QC violin plots to understand your data distribution
2. Increase MAD thresholds (more permissive):
   ```yaml
   qc:
     mad_n_counts: 6    # Default: 5
     mad_n_genes: 6     # Default: 5
     mad_mito: 4        # Default: 3
   ```

3. Or switch to fixed mode with custom thresholds based on your plots:
   ```bash
   python create_config.py ... \
       --qc-mode fixed \
       --max-genes 8000 \
       --min-counts 300 \
       --max-counts 80000
   ```

#### Too few cells filtered

If you suspect low-quality cells are passing through:

1. Decrease MAD thresholds (more stringent):
   ```yaml
   qc:
     mad_n_counts: 3
     mad_n_genes: 3
     mad_mito: 2
   ```

2. Or increase `min_genes`:
   ```yaml
   qc:
     min_genes: 500
   ```

#### Comparing Adaptive vs Fixed Mode

To compare filtering approaches on your data:

```python
# In a Jupyter notebook after running the pipeline
import scanpy as sc

# Check outlier flags (adaptive mode adds these)
adata = sc.read_h5ad("results/annotation/adata_annotated.h5ad")

# If adaptive mode was used, these columns exist:
if 'outlier' in adata.obs.columns:
    print(f"General outliers: {adata.obs['outlier'].sum()}")
    print(f"MT outliers: {adata.obs['mt_outlier'].sum()}")
```

### SComatic Issues

#### SComatic script errors

Verify SComatic installation:

```bash
ls /path/to/SComatic/scripts/
# Should contain: SplitBamCellTypes.py, BaseCellCounter.py, etc.
```

#### No mutations detected

1. Check BAM files contain cell barcodes (CB tag)
2. Verify reference files are correct
3. Check callable sites output
4. Lower `min_cov` and `min_cells` thresholds

### Log Files

All log files are stored in `{output_dir}/logs/`:

```bash
# View Cell Ranger log
cat results/logs/cellranger/SAMPLE1.log

# View QC/annotation log
cat results/logs/qc_annotation/qc_annotation.log

# View Snakemake log
snakemake ... 2>&1 | tee snakemake.log
```

### Debugging

Run with verbose output:

```bash
snakemake --configfile config.yaml \
    --cores 8 \
    --verbose \
    --printshellcmds
```

Dry run to check workflow:

```bash
snakemake --configfile config.yaml --dryrun
```

---

## Citation

If you use ClusterCatcher in your research, please cite:

> Lehle, J. (2025). ClusterCatcher: Single-cell sequencing analysis pipeline for mutation signature detection. GitHub. https://github.com/JakeLehle/ClusterCatcher

Please also cite the tools used in each module:

| Module | Citation |
|--------|----------|
| Cell Ranger | 10x Genomics |
| Scanpy | Wolf et al., Genome Biology 2018 |
| popV | https://github.com/YosefLab/popV |
| CytoTRACE2 | Gulati et al., Science 2020 |
| inferCNV | https://github.com/broadinstitute/infercnv |
| SComatic | Muyas et al., Nature Biotechnology 2024 |
| Kraken2 | Wood et al., Genome Biology 2019 |
| COSMIC | Alexandrov et al., Nature 2020 |

### QC Best Practices References

- Luecken, M.D. & Theis, F.J. (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology*, 15(6), e8746.
- Heumos, L. et al. (2023). Best practices for single-cell analysis across modalities. *Nature Reviews Genetics*, 24, 550–572.

---

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

---

## Contact

- **Author**: Jake Lehle
- **Institution**: Texas Biomedical Research Institute
- **Issues**: [GitHub Issues](https://github.com/JakeLehle/ClusterCatcher/issues)

For bug reports or feature requests, please open an issue on GitHub.
