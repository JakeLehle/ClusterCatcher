# =============================================================================
# CLUSTERCATCHER SNAKEMAKE RULES - SCOMATIC & SIGNATURE ANALYSIS
# =============================================================================
# These rules should be included in the main Snakefile
# Include with: include: "rules/scomatic_signatures.smk"
# =============================================================================

# ----------------------------------------------------------------------------
# Somatic Mutation Rules (SComatic)
# ----------------------------------------------------------------------------

rule scomatic_mutation_calling:
    """
    Run SComatic pipeline for somatic mutation detection.
    
    This rule:
    1. Filters BAMs to annotated cells
    2. Splits BAMs by cell type
    3. Counts bases per cell
    4. Calls variants with filtering
    5. Generates single-cell genotypes
    6. Creates filtered mutation matrix with trinucleotide context
    
    Requires:
    - SComatic installed and scripts accessible
    - Reference genome FASTA
    - RNA editing sites file
    - Panel of Normals (PoN) file
    - BED file for mappable regions
    """
    input:
        adata=f"{OUTPUT_DIR}/dysregulation/adata_cancer_detected.h5ad",
        bams=expand(
            f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/possorted_genome_bam.bam",
            sample=SAMPLE_IDS
        ),
    output:
        mutations=f"{OUTPUT_DIR}/mutations/all_samples.single_cell_genotype.filtered.tsv",
        callable_sites=f"{OUTPUT_DIR}/mutations/CombinedCallableSites/complete_callable_sites.tsv",
        annotations=f"{OUTPUT_DIR}/mutations/cell_annotations.tsv",
        trinuc=f"{OUTPUT_DIR}/mutations/trinucleotide_background.tsv",
    params:
        output_dir=OUTPUT_DIR,
        sample_ids=SAMPLE_IDS,
        scomatic_scripts_dir=config.get("scomatic", {}).get("scripts_dir"),
        ref_genome=REFERENCE_FASTA,
        editing_sites=config.get("scomatic", {}).get("editing_sites"),
        pon_file=config.get("scomatic", {}).get("pon_file"),
        bed_file=config.get("scomatic", {}).get("bed_file"),
        cellranger_dir=f"{OUTPUT_DIR}/cellranger",
        custom_genotype_script=workflow.source_path("scripts/SingleCellGenotype.py"),
    log:
        f"{LOG_DIR}/mutations/scomatic.log"
    threads: THREADS
    conda:
        "envs/scomatic.yaml"
    script:
        "scripts/scomatic_mutation_calling.py"


# Simplified callable sites path alias
rule callable_sites_link:
    """Create simplified link to callable sites file."""
    input:
        f"{OUTPUT_DIR}/mutations/CombinedCallableSites/complete_callable_sites.tsv"
    output:
        f"{OUTPUT_DIR}/mutations/callable_sites.tsv"
    shell:
        "ln -sf CombinedCallableSites/complete_callable_sites.tsv {output}"


# ----------------------------------------------------------------------------
# Mutational Signature Rules
# ----------------------------------------------------------------------------

rule signature_analysis:
    """
    Run semi-supervised COSMIC signature analysis.
    
    This rule:
    1. Converts mutations to 96-trinucleotide context matrix
    2. Extracts relevant COSMIC signatures
    3. Uses scree plot elbow detection for signature selection (optional)
    4. Fits signatures using NNLS
    5. Evaluates reconstruction quality
    6. Adds signature weights to AnnData
    7. Generates comprehensive visualizations (UMAPs, etc.)
    
    Requires:
    - COSMIC signature database (v3.4 recommended)
    """
    input:
        mutations=f"{OUTPUT_DIR}/mutations/all_samples.single_cell_genotype.filtered.tsv",
        adata=f"{OUTPUT_DIR}/dysregulation/adata_cancer_detected.h5ad",
        callable_sites=f"{OUTPUT_DIR}/mutations/CombinedCallableSites/complete_callable_sites.tsv",
    output:
        weights=f"{OUTPUT_DIR}/signatures/signature_weights_per_cell.txt",
        adata_final=f"{OUTPUT_DIR}/signatures/adata_final.h5ad",
    params:
        output_dir=f"{OUTPUT_DIR}/signatures",
        cosmic_file=config.get("signatures", {}).get("cosmic_file"),
        use_scree_plot=config.get("signatures", {}).get("use_scree_plot", False),
        core_signatures=config.get("signatures", {}).get("core_signatures", ["SBS2", "SBS13", "SBS5"]),
        candidate_order=config.get("signatures", {}).get("candidate_order"),
        mutation_threshold=config.get("signatures", {}).get("mutation_threshold", 0),
        max_signatures=config.get("signatures", {}).get("max_signatures", 15),
        hnscc_only=config.get("signatures", {}).get("hnscc_only", False),
    log:
        f"{LOG_DIR}/signatures/signature_analysis.log"
    threads: THREADS
    conda:
        "envs/signatures.yaml"
    script:
        "scripts/signature_analysis.py"


# ----------------------------------------------------------------------------
# Helper rule to copy final AnnData to main output
# ----------------------------------------------------------------------------

rule copy_final_adata:
    """Copy final annotated AnnData to main output directory."""
    input:
        f"{OUTPUT_DIR}/signatures/adata_final.h5ad"
    output:
        f"{OUTPUT_DIR}/adata_final.h5ad"
    shell:
        "cp {input} {output}"

