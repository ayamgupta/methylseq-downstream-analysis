#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
METHYLATION ANALYSIS PIPELINE
========================================================================================
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --samplesheet samplesheet.csv --input_dir bismark_results/ --gtf genes.gtf

    Mandatory arguments:
      --samplesheet             Path to sample sheet CSV file
      --input_dir               Path to input directory containing methylation files
      --gtf                     Path to GTF annotation file

    Optional arguments:
      --file_type               Type of input files: 'bismark' or 'methyldeckyl' (default: 'bismark')
      --outdir                  Output directory (default: './results')
      --min_coverage            Minimum coverage threshold (default: 20)
      --high_percentile         High percentile for coverage filtering (default: 99.9)
      --meth_diff_threshold     Methylation difference threshold (default: 25)
      --qval_threshold          Q-value threshold (default: 0.05)
      --top_gene_count          Number of top genes to select (default: 10)
      --min_genes_for_enrichment Minimum genes required for enrichment (default: 3)
      --enrichr_databases       Comma-separated list of EnrichR databases

    Other options:
      --help                    Show this help message and exit
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
PARAMETER VALIDATION
========================================================================================
*/

// Check mandatory parameters
if (!params.samplesheet) {
    log.error "Please provide a samplesheet with --samplesheet"
    exit 1
}

if (!params.input_dir) {
    log.error "Please provide input directory with --input_dir"
    exit 1
}

if (!params.gtf) {
    log.error "Please provide GTF file with --gtf"
    exit 1
}

/*
========================================================================================
PARAMETER DEFAULTS
========================================================================================
*/

params.file_type = 'bismark'
params.outdir = './results'
params.min_coverage = 20
params.high_percentile = 99.9
params.meth_diff_threshold = 25
params.qval_threshold = 0.05
params.top_gene_count = 10
params.min_genes_for_enrichment = 3
params.enrichr_databases = 'GO_Molecular_Function_2023,GO_Biological_Process_2023,GO_Cellular_Component_2023,Reactome_2022,WikiPathway_2023_Human,DGIdb_Drug_Targets_2024'

/*
========================================================================================
WORKFLOWS
========================================================================================
*/

workflow {
    
    // Create channels
    samplesheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
    input_dir_ch = Channel.fromPath(params.input_dir, type: 'dir', checkIfExists: true)
    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)
    
    // Run methylation analysis
    METHYLATION_ANALYSIS(
        samplesheet_ch,
        input_dir_ch,
        gtf_ch
    )
}

/*
========================================================================================
PROCESSES
========================================================================================
*/

process METHYLATION_ANALYSIS {
    tag "methylation_analysis"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path samplesheet
    path input_dir
    path gtf
    
    output:
    path "*.csv", emit: csv_files
    path "*.png", optional: true, emit: plots
    path "enrichment_*.csv", optional: true, emit: enrichment
    
    script:
    def r_script = "${projectDir}/bin/methylation_analysis.R"
    
    """
    Rscript ${r_script} --metadata_csv ${samplesheet} --base_dir ${input_dir} --gtf_file ${gtf} --file_type ${params.file_type} --min_coverage ${params.min_coverage} --high_percentile ${params.high_percentile} --meth_diff_threshold ${params.meth_diff_threshold} --qval_threshold ${params.qval_threshold} --top_gene_count ${params.top_gene_count} --min_genes_for_enrichment ${params.min_genes_for_enrichment} --enrichr_databases "${params.enrichr_databases}"
    """
    
    stub:
    """
    touch methylkit_pca_scores.csv
    touch volcano_input.csv
    touch heatmap_input.csv
    touch enrichment_hyper_combined.csv
    touch enrichment_hypo_combined.csv
    """
}

/*
========================================================================================
WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    if (workflow.success) {
        log.info "Pipeline completed successfully!"
        log.info "Results are in: ${params.outdir}"
    } else {
        log.error "Pipeline failed!"
    }
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}