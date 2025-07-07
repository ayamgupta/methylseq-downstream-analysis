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
    nextflow run main.nf --samplesheet samplesheet.csv --gtf genes.gtf

    Mandatory arguments:
      --samplesheet             Path to sample sheet CSV file (must contain 'path' column)
      --gtf                     Path to GTF annotation file

    Optional arguments:
      --input_dir               Optional base directory for relative file paths in samplesheet
      --file_type               Type of input files: 'bismark' or 'methyldackel' (default: 'bismark')
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

    Note: The samplesheet must contain 'path' column with individual file paths.
    If paths are relative, use --input_dir to specify the base directory.
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

if (!params.gtf) {
    log.error "Please provide GTF file with --gtf"
    exit 1
}

/*
========================================================================================
WORKFLOWS
========================================================================================
*/

workflow {
    
    // Create channels - ensure files exist
    samplesheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)
    r_script_ch = Channel.fromPath("${projectDir}/bin/methylation_analysis.R", checkIfExists: true)
    
    // Run methylation analysis
    METHYLATION_ANALYSIS(
        samplesheet_ch,
        gtf_ch,
        r_script_ch
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
    path gtf
    path r_script
    
    output:
    path "*.csv", emit: csv_files
    path "*.png", optional: true, emit: plots
    path "enrichment_*.csv", optional: true, emit: enrichment
    
    script:
    def base_dir_arg = params.input_dir ? "--base_dir ${params.input_dir}" : ""
    
    """
    Rscript ${r_script} \\
        --metadata_csv ${samplesheet} \\
        ${base_dir_arg} \\
        --gtf_file ${gtf} \\
        --file_type ${params.file_type} \\
        --min_coverage ${params.min_coverage} \\
        --high_percentile ${params.high_percentile} \\
        --meth_diff_threshold ${params.meth_diff_threshold} \\
        --qval_threshold ${params.qval_threshold} \\
        --top_gene_count ${params.top_gene_count} \\
        --min_genes_for_enrichment ${params.min_genes_for_enrichment} \\
        --output_dir . \\
        --enrichr_databases "${params.enrichr_databases}"
    """
    
    stub:
    """
    touch pca.csv
    touch volcano.csv
    touch grapher.csv
    touch cohort_to_condition_mapping.csv
    touch condition_comparison.csv
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