#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Single-Cell QC Metrics Processing Pipeline
 * 
 * This pipeline ingests metrics CSVs and data files,
 * processes them with a Python script, and generates output CSVs.
 */

// Define parameters
params.num_input_cells = null          // Number of input cells (integer)
params.data_dir = null                 // Path to directory containing all data files
params.file_prefix = null              // Prefix for all input files
params.run_guide_assignment = true     // Whether to run guide assignment
params.outdir = "results"              // Output directory
params.help = false

// Help message
def helpMessage() {
    log.info"""
    Usage:
      nextflow run main.nf --num_input_cells <int> --data_dir <path> --file_prefix <prefix> [options]

    Required arguments:
      --num_input_cells          Number of input cells (integer)
      --data_dir                 Path to directory containing all data files
      --file_prefix              Common prefix for all input files

    Expected files in data_dir:
      <prefix>.scRNA_metrics.csv
      <prefix>.filtered.matrix.mtx.gz
      <prefix>.filtered.barcodes.tsv.gz
      <prefix>.filtered.features.tsv.gz

    Optional arguments:
      --run_guide_assignment     Whether to run CRISPR guide assignment (default: ${params.run_guide_assignment})
      --outdir                   Output directory (default: ${params.outdir})
      --help                     Show this help message
    """.stripIndent()
}

// Import modules
include { GENERATE_REPORT_DATA } from './modules/generate_report_data'
include { EXTRACT_CRISPR_FEATURES } from './modules/extract_crispr_features'
include { GUIDE_ASSIGNMENT } from './modules/guide_assignment'

/*
 * Main workflow
 */
workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate inputs
    if (!params.num_input_cells) {
        log.error "ERROR: --num_input_cells is required"
        helpMessage()
        exit 1
    }

    if (!params.data_dir) {
        log.error "ERROR: --data_dir is required"
        helpMessage()
        exit 1
    }

    if (!params.file_prefix) {
        log.error "ERROR: --file_prefix is required"
        helpMessage()
        exit 1
    }

    // Construct file paths using data_dir and file_prefix
    def metrics_path = "${params.data_dir}/${params.file_prefix}.scRNA_metrics.csv"
    def matrix_path = "${params.data_dir}/${params.file_prefix}.filtered.matrix.mtx.gz"
    def barcodes_path = "${params.data_dir}/${params.file_prefix}.filtered.barcodes.tsv.gz"
    def features_path = "${params.data_dir}/${params.file_prefix}.filtered.features.tsv.gz"

    // Log pipeline parameters
    log.info """
        ========================================
        Single-Cell QC Metrics Processing
        ========================================
        Num input cells      : ${params.num_input_cells}
        Data directory       : ${params.data_dir}
        File prefix          : ${params.file_prefix}
        scRNA metrics        : ${metrics_path}
        Filtered matrix      : ${matrix_path}
        Filtered barcodes    : ${barcodes_path}
        Filtered features    : ${features_path}
        Run guide assignment : ${params.run_guide_assignment}
        Output dir           : ${params.outdir}
        ========================================
        """.stripIndent()

    // Create channels from input files
    num_cells_ch = Channel.value(params.num_input_cells)
    metrics_ch = Channel.fromPath(metrics_path, checkIfExists: true)
    matrix_ch = Channel.fromPath(matrix_path, checkIfExists: true)
    barcodes_ch = Channel.fromPath(barcodes_path, checkIfExists: true)
    features_ch = Channel.fromPath(features_path, checkIfExists: true)
    
    // Conditionally run guide assignment if requested
    if (params.run_guide_assignment) {
        log.info "Running CRISPR feature extraction and guide assignment..."
        
        // Step 1: Extract CRISPR features and create h5ad file
        crispr_input_ch = matrix_ch
            .combine(barcodes_ch)
            .combine(features_ch)
        
        EXTRACT_CRISPR_FEATURES(crispr_input_ch)
        
        EXTRACT_CRISPR_FEATURES.out.crispr_adata.view { "Generated CRISPR h5ad: $it" }
        
        // Step 2: Run CRISPAT guide assignment using the h5ad file
        guide_input_ch = EXTRACT_CRISPR_FEATURES.out.crispr_adata
            .combine(metrics_ch)
        
        GUIDE_ASSIGNMENT(guide_input_ch)
        
        // Use guide assignment output for metrics processing
        assignments_ch = GUIDE_ASSIGNMENT.out.assignments
        
        GUIDE_ASSIGNMENT.out.assignments.view { "Generated guide assignments: $it" }
    } else {
        log.info "Skipping guide assignment"
        // Create empty channel for assignments
        assignments_ch = Channel.fromPath('NO_FILE')
    }
    
    // Combine all inputs for metrics processing
    input_ch = num_cells_ch
        .combine(metrics_ch)
        .combine(matrix_ch)
        .combine(barcodes_ch)
        .combine(features_ch)
        .combine(assignments_ch)
    
    // Process metrics
    GENERATE_REPORT_DATA(input_ch)
    
    // Emit results
    GENERATE_REPORT_DATA.out.results.view { "Generated output: $it" }
}
