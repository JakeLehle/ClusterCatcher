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
│  │  1. Cell Ranger  │───▶│   BAM + Matrix   │                              │
│  │   (Alignment)    │     │                  │                              │
│  └──────────────────┘     └────────┬─────────┘                              │
│                                    │                                        │
│           ┌────────────────────────┼────────────────────────┐               │
│           │                        │                        │               │
│           ▼                        ▼                        ▼               │
│  ┌──────────────────┐     ┌──────────────────┐     ┌──────────────────┐     │
│  │ 2. QC & Filter   │     │ 5. Viral Detect  │     │ 6. SComatic      │     │
│  │   (Scanpy)       │     │   (Kraken2)      │     │ (Mutations)      │     │
│  └────────┬─────────┘     └────────┬─────────┘     └────────┬─────────┘     │
│           │                        │                        │               │
│           ▼                        ▼                        │               │
│  ┌──────────────────┐     ┌──────────────────┐              │               │
│  │ 3. Cell Annot.   │     │ Viral Integration│              │               │
│  │   (popV)         │     │                  │              │               │
│  └────────┬─────────┘     └──────────────────┘              │               │
│           │                                                 │               │
│           ▼                                                 │               │
│  ┌──────────────────┐                                       │               │
│  │ 4. Dysregulation │                                       │               │
│  │ (CytoTRACE2+CNV) │                                       │               │
│  └────────┬─────────┘                                       │               │
│           │                                                 │               │
│           └──────────────────────────┬──────────────────────┘               │
│                                      │                                      │
│                                      ▼                                      │
│                             ┌──────────────────┐                            │
│                             │ 7. Signatures    │                            │
│                             │   (NNLS/COSMIC)  │                            │
│                             └────────┬─────────┘                            │
│                                      │                                      │
│                                      ▼                                      │
│                             ┌──────────────────┐                            │
│                             │  Final AnnData   │                            │
│                             │  + Summaries     │                            │
│                             └──────────────────┘                            │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Pipeline Modules Summary

| # | Module | Script | Required | Description |
|---|--------|--------|----------|-------------|
| 1 | Cell Ranger | `cellranger_count.py` | Yes | FASTQ alignment, counting, BAM generation |
| 2 | QC & Filtering | `scanpy_qc_annotation.py` | Yes | Quality control, doublet removal, filtering |
| 3 | Cell Annotation | `scanpy_qc_annotation.py` | Yes | popV annotation with cluster-based refinement |
| 4 | Dysregulation | `cancer_cell_detection.py` | Default | Cancer cell detection (CytoTRACE2 + inferCNV) |
| 5 | Viral Detection | `kraken2_viral_detection.py`, `summarize_viral_detection.py`, `viral_integration.py` | Optional | Pathogen detection in unmapped reads |
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

**Override if needed** (for non-standard reference layouts):

```bash
python snakemake_wrapper/create_config.py \
    --output-dir ./results \
    --sample-pickle samples.pkl \
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A \
    --reference-fasta /custom/path/to/genome.fa \
    --gtf-file /custom/path/to/genes.gtf \
    --threads 16 \
    --memory-gb 64
```

### Step 3: Run Pipeline

```bash
# Dry run first (recommended)
cd snakemake_wrapper
snakemake --configfile ../results/config.yaml --dryrun

# Full run
snakemake --configfile ../results/config.yaml \
    --cores 32 \
    --use-conda
```

### For SRAscraper Users

If you used [SRAscraper](https://github.com/JakeLehle/SRAscraper) to download data:

```bash
python snakemake_wrapper/create_config.py \
    --output-dir ./results \
    --sample-pickle /path/to/SRAscraper_output/metadata/sample_dict.pkl \
    --cellranger-reference /path/to/refdata-gex-GRCh38-2020-A
```

---

## Detailed Usage

### Sample CSV Format

#### Required Columns

| Column | Description | Example |
|--------|-------------|---------|
| `sample_id` | Unique sample identifier | `SAMPLE1` |
| `fastq_r1` | Path to R1 FASTQ file(s) | `/path/to/R1.fastq.gz` |
| `fastq_r2` | Path to R2 FASTQ file(s) | `/path/to/R2.fastq.gz` |

#### Optional Columns

| Column | Description |
|--------|-------------|
| `sample_name` | Human-readable name |
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

Performs quality control, doublet detection, and cell filtering using Scanpy.

#### QC Metrics Calculated
- Number of genes per cell
- Total UMI counts per cell
- Mitochondrial gene percentage
- Ribosomal gene percentage
- Doublet scores (Scrublet)

#### Filtering Steps
1. Remove cells with < `min_genes` genes detected
2. Remove cells with > `max_genes` genes detected
3. Remove cells with < `min_counts` total counts
4. Remove cells with > `max_counts` total counts
5. Remove cells with > `max_mito_pct` mitochondrial reads
6. Remove predicted doublets

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_genes` | `200` | Minimum genes per cell |
| `max_genes` | `5000` | Maximum genes per cell |
| `min_counts` | `500` | Minimum UMI counts |
| `max_counts` | `50000` | Maximum UMI counts |
| `max_mito_pct` | `20` | Maximum mitochondrial % |
| `min_cells` | `3` | Minimum cells expressing a gene |
| `doublet_removal` | `true` | Enable doublet detection |
| `doublet_rate` | `0.08` | Expected doublet rate |

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

**Key features:**
- Runs on **raw counts** (before normalization)
- Downloads pre-trained models from Hugging Face
- Uses batch correction during annotation
- Transfers predictions back to original object

**Available models** (set via `popv_huggingface_repo`):

| Model | Description |
|-------|-------------|
| `popV/tabula_sapiens_All_Cells` | All cell types (default, recommended) |
| `popV/tabula_sapiens_immune` | Immune cells only |
| `popV/human_lung` | Lung-specific |

See [Hugging Face popV](https://huggingface.co/popV) for more models.

#### Cluster-Based Annotation Refinement

After popV annotation, ClusterCatcher refines predictions using cluster-level weighted scoring:

1. **Min-Max normalization** of `popv_prediction_score` within each cluster
2. **Linear weighting** to sum to 1 within each cluster
3. **Weighted aggregation** of scores per cell type per cluster
4. **Assignment** of dominant cell type to all cells in cluster

This approach:
- Reduces noise from low-confidence predictions
- Leverages cluster structure for more robust assignments
- Produces cleaner cell type boundaries

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `popv_huggingface_repo` | `popV/tabula_sapiens_All_Cells` | Hugging Face model repository |
| `popv_prediction_mode` | `inference` | `inference` (full) or `fast` |
| `popv_gene_symbol_key` | `feature_name` | Gene symbol column in adata.var |
| `popv_cache_dir` | `tmp/popv_models` | Cache for downloaded models |
| `batch_key` | `sample_id` | Batch correction key |

#### Post-Annotation Processing Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `target_sum` | `1000000` | Normalization target (1e6 = CPM) |
| `n_pcs` | `null` | PCs for neighbors (null = CPU count) |
| `n_neighbors` | `15` | Neighbors for graph construction |
| `leiden_resolution` | `1.0` | Clustering resolution |
| `run_bbknn` | `false` | Enable BBKNN batch correction |
| `bbknn_batch_key` | `sample_id` | Batch key for BBKNN |

#### Outputs

- `annotation/adata_annotated.h5ad` - Annotated AnnData with:
  - `popv_prediction`: Raw popV cell type calls
  - `popv_prediction_score`: Confidence scores
  - `clusters`: Leiden cluster assignments
  - `final_annotation`: Cluster-refined cell types
- `annotation/annotation_summary.tsv` - Cell type counts per sample
- `figures/qc/UMAP_popv_prediction.pdf` - UMAP colored by popV calls
- `figures/qc/UMAP_popv_prediction_score.pdf` - UMAP colored by confidence
- `figures/qc/UMAP_final_annotation.pdf` - UMAP with non-overlapping labels
- `figures/qc/UMAP_samples.pdf` - UMAP colored by sample
- `figures/qc/UMAP_clusters.pdf` - UMAP colored by cluster
- `figures/qc/Stacked_Bar_Cluster_Composition.pdf` - Cell type composition per cluster
- `figures/qc/cell_type_proportions.pdf` - Cell type proportions per sample

---

### Module 4: Dysregulation Detection (Cancer Cell Identification)

**Script**: `scripts/cancer_cell_detection.py`

Identifies cancer/dysregulated cells using a dual-model consensus approach.

#### How It Works

```
┌─────────────────┐     ┌─────────────────┐
│   CytoTRACE2    │     │    InferCNV     │
│  (Stemness)     │     │  (CNV Scores)   │
└────────┬────────┘     └────────┬────────┘
         │                       │
         ▼                       ▼
    ┌─────────┐            ┌─────────┐
    │ Potency │            │   CNV   │
    │  Score  │            │  Score  │
    └────┬────┘            └────┬────┘
         │                      │
         └──────────┬───────────┘
                    ▼
           ┌──────────────┐
           │  Agreement   │
           │    Score     │
           └──────┬───────┘
                  ▼
           ┌──────────────┐
           │ Cancer/Normal│
           │Classification│
           └──────────────┘
```

#### CytoTRACE2 Component
- Scores cell developmental potential (0-1)
- Higher scores = more stem-like/cancer-like
- Automatic bimodal threshold detection

#### InferCNV Component
- Uses normal cell types as reference
- Detects copy number variations
- Scores chromosomal instability

#### Agreement Scoring
- Combines both models using weighted correlation
- `alpha` parameter controls rank vs. value weighting
- Cells classified as cancer only if both models agree

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cytotrace2.enabled` | `true` | Enable CytoTRACE2 |
| `cytotrace2.species` | `human` | Species for model |
| `infercnv.enabled` | `true` | Enable inferCNV |
| `infercnv.window_size` | `250` | CNV sliding window |
| `infercnv_reference_groups` | `null` | Normal cell types for reference |
| `agreement.alpha` | `0.5` | Rank (0) vs. value (1) weight |
| `agreement.min_correlation` | `0.5` | Minimum Spearman rho |

#### Outputs
- `dysregulation/adata_cancer_detected.h5ad` - AnnData with cancer labels
- `dysregulation/cancer_detection_summary.tsv` - Classification summary
- `dysregulation/figures/` - CytoTRACE2, CNV, and agreement plots

---

### Module 5: Viral Detection

**Scripts**: 
- `scripts/kraken2_viral_detection.py` - Per-sample viral detection
- `scripts/summarize_viral_detection.py` - Cross-sample summary
- `scripts/viral_integration.py` - Integration with expression data

Detects viral/microbial sequences in unmapped reads using Kraken2.

#### Workflow

```
┌──────────────────┐
│ Cell Ranger BAM  │
│ (unmapped reads) │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Extract Unmapped │
│    (samtools)    │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│    Kraken2       │
│ Classification   │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Build SC Matrix  │
│(organism x cell) │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Filter Human     │
│   Viruses Only   │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Integrate with   │
│ Gene Expression  │
└──────────────────┘
```

#### Outputs
- `viral/{sample}/kraken2_filtered_feature_bc_matrix/` - Per-sample viral count matrices
- `viral/{sample}/{sample}_organism_summary.tsv` - Per-sample organism summary
- `viral/viral_detection_summary.tsv` - Cross-sample summary
- `viral/viral_counts.h5ad` - Combined viral AnnData
- `viral_integration/adata_with_virus.h5ad` - Expression + viral data
- `viral_integration/figures/` - Viral detection visualizations

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `kraken_db` | required | Kraken2 database path |
| `viral_db` | optional | Human viral inspect.txt for filtering |
| `confidence` | `0.1` | Kraken2 confidence threshold |
| `include_organisms` | `null` | Organism patterns to include |
| `exclude_organisms` | `null` | Organism patterns to exclude |

#### Kraken2 Database Setup

```bash
# Option A: Download pre-built viral database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20231009.tar.gz
tar -xzf k2_viral_20231009.tar.gz -C ~/kraken2_db/

# Option B: Build custom human viral database
kraken2-build --download-taxonomy --db ~/kraken2_db/human_viral
datasets download virus genome taxon "human virus" --filename human_viruses.zip
unzip human_viruses.zip
kraken2-build --add-to-library ncbi_dataset/data/genomic.fna --db ~/kraken2_db/human_viral
kraken2-build --build --db ~/kraken2_db/human_viral --threads 16
kraken2-inspect --db ~/kraken2_db/human_viral > ~/kraken2_db/human_viral/inspect.txt
```

---

### Module 6: SComatic Mutation Calling

**Script**: `scripts/scomatic_mutation_calling.py`

Comprehensive 10-phase somatic mutation calling pipeline using SComatic.

#### Workflow Phases

| Phase | Description | Output |
|-------|-------------|--------|
| 1 | Filter BAM to annotated cells | Filtered BAM |
| 2 | Split BAM by cell type | Per-cell-type BAMs |
| 3 | Count bases per cell | Base count tables |
| 4 | Merge counts across samples | Combined counts |
| 5 | Variant calling (Step 1) | Raw variants |
| 6 | Filter with PoN/editing sites | Filtered variants |
| 7 | BED file filtering | Mappable variants |
| 8 | Callable sites (cell type) | Cell type callable |
| 9 | Callable sites (per cell) | Per-cell callable |
| 10 | Single-cell genotypes | Final mutations |

#### Required Reference Files

| File | Description | How to Obtain |
|------|-------------|---------------|
| `editing_sites` | Known RNA editing positions | SComatic repo or RADAR database |
| `pon_file` | Panel of Normals | Generate from matched normals |
| `bed_file` | Mappable genomic regions | Download or generate |

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `scripts_dir` | required | Path to SComatic/scripts |
| `editing_sites` | required | RNA editing sites file |
| `pon_file` | required | Panel of Normals |
| `bed_file` | required | Mappable regions BED |
| `min_cov` | `5` | Minimum coverage |
| `min_cells` | `5` | Minimum cells with variant |
| `min_base_quality` | `30` | Base quality threshold |
| `min_map_quality` | `30` | Mapping quality threshold |

#### Outputs
- `mutations/all_samples.single_cell_genotype.filtered.tsv` - Final mutations
- `mutations/CombinedCallableSites/complete_callable_sites.tsv` - Callable sites
- `mutations/cell_annotations.tsv` - Cell barcode to type mapping
- `mutations/trinucleotide_background.tsv` - Background frequencies
- `mutations/{sample}/` - Per-sample intermediate files

---

### Module 7: Mutational Signature Analysis

**Script**: `scripts/signature_analysis.py`

Deconvolves single-cell mutations into COSMIC mutational signatures using semi-supervised NNLS (Non-Negative Least Squares).

#### Workflow

```
┌──────────────────┐
│ Single-Cell      │
│ Mutations (TSV)  │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Build 96-Context │
│ Mutation Matrix  │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Aggregate Cell   │
│ Profiles         │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ NNLS Fitting to  │
│ COSMIC Signatures│
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Add Signature    │
│ Weights to       │
│ AnnData          │
└──────────────────┘
```

#### Signature Selection Methods

1. **All Signatures**: Fit all provided COSMIC signatures
2. **HNSCC-Specific** (`--hnscc-only`): Use curated set of 15 HNSCC-relevant signatures
3. **Scree Plot** (`--use-scree-plot`): Automatic selection using elbow detection

#### HNSCC Signature Set

When `--hnscc-only` is enabled, these signatures are used:

| Signature | Proposed Etiology |
|-----------|-------------------|
| SBS1 | Spontaneous deamination (age) |
| SBS2 | APOBEC activity |
| SBS4 | Tobacco smoking |
| SBS5 | Clock-like (unknown) |
| SBS7a/b | UV exposure |
| SBS13 | APOBEC activity |
| SBS16 | Unknown (liver-associated) |
| SBS17a/b | Unknown |
| SBS18 | Reactive oxygen species |
| SBS29 | Tobacco chewing |
| SBS39 | Unknown |
| SBS40 | Unknown |
| SBS44 | Defective DNA mismatch repair |

#### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cosmic_file` | required | COSMIC signatures file |
| `core_signatures` | `SBS2,SBS13,SBS5` | Always include these |
| `use_scree_plot` | `false` | Use elbow detection |
| `mutation_threshold` | `0` | Min mutations per cell |
| `max_signatures` | `15` | Max signatures to test |
| `hnscc_only` | `false` | Use HNSCC signature set |

#### Outputs
- `signatures/signature_weights_per_cell.txt` - Per-cell signature weights
- `signatures/adata_final.h5ad` - Final AnnData with all annotations
- `signatures/mutation_matrix_96contexts.txt` - 96-context mutation matrix
- `signatures/cosmic_signatures_used.txt` - Signatures used
- `signatures/reconstruction_evaluation.txt` - Quality metrics
- `signatures/figures/` - Visualization plots

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

**How it works**:
1. You only need to specify `cellranger` (the Cell Ranger reference directory)
2. The pipeline automatically derives:
   - `fasta` from `{cellranger}/fasta/genome.fa`
   - `gtf` from `{cellranger}/genes/genes.gtf`
3. Override these only if your reference has a non-standard layout

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

# QC parameters
qc:
  min_genes: 200
  max_genes: 5000
  min_counts: 500
  max_counts: 50000
  max_mito_pct: 20
  min_cells: 3
  doublet_removal: true
  doublet_rate: 0.08

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
| `--min-genes` | `200` | Min genes per cell |
| `--max-genes` | `5000` | Max genes per cell |
| `--min-counts` | `500` | Min UMI counts |
| `--max-counts` | `50000` | Max UMI counts |
| `--max-mito-pct` | `20` | Max mitochondrial % |
| `--min-cells` | `3` | Min cells expressing gene |
| `--doublet-rate` | `0.08` | Expected doublet rate |

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

---

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

---

## Contact

- **Author**: Jake Lehle
- **Institution**: Texas Biomedical Research Institute
- **Issues**: [GitHub Issues](https://github.com/JakeLehle/ClusterCatcher/issues)

For bug reports or feature requests, please open an issue on GitHub.
