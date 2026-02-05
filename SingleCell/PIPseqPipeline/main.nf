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
params.dragen_file_prefix = null              // Prefix for all input files
params.sample_id = null                // Sample identifier (defaults to file_prefix if not provided)
params.output_basename = null          // Output basename (defaults to file_prefix if not provided)
params.run_guide_assignment = true     // Whether to run guide assignment
params.outdir = "results"              // Output directory
params.help = false

// Help message
def helpMessage() {
    log.info"""
    Usage:
      nextflow run main.nf --num_input_cells <int> --data_dir <path> --dragen_file_prefix <prefix> [options]

    Required arguments:
      --num_input_cells          Number of input cells (integer)
      --data_dir                 Path to directory containing all data files
      --dragen_file_prefix       Common prefix for all input files
      --sample_id                Sample identifier
      --output_basename          Output basename for generated files

    Expected files in data_dir:
      <dragen_file_prefix>.scRNA_metrics.csv
      <dragen_file_prefix>.filtered.matrix.mtx.gz
      <dragen_file_prefix>.filtered.barcodes.tsv.gz
      <dragen_file_prefix>.filtered.features.tsv.gz

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

    if (!params.dragen_file_prefix) {
        log.error "ERROR: --dragen_file_prefix is required"
        helpMessage()
        exit 1
    }

    if (!params.sample_id) {
        log.error "ERROR: --sample_id is required"
        helpMessage()
        exit 1
    }

    if (!params.output_basename) {
        log.error "ERROR: --output_basename is required"
        helpMessage()
        exit 1
    }

    // Force all parameters to strings to preserve leading zeros
    // (Nextflow may parse numeric-looking strings as numbers)
    def dragen_file_prefix = "${params.dragen_file_prefix}"
    def sample_id = "${params.sample_id}"
    def output_basename = "${params.output_basename}"
    
    // Construct file paths using data_dir and file_prefix
    def metrics_path = "${params.data_dir}/${dragen_file_prefix}.scRNA_metrics.csv"
    def barcode_summary_path = "${params.data_dir}/${dragen_file_prefix}.scRNA.barcodeSummary.tsv"
    def matrix_path = "${params.data_dir}/${dragen_file_prefix}.filtered.matrix.mtx.gz"
    def barcodes_path = "${params.data_dir}/${dragen_file_prefix}.filtered.barcodes.tsv.gz"
    def features_path = "${params.data_dir}/${dragen_file_prefix}.filtered.features.tsv.gz"
    // Log pipeline parameters
    log.info """
        ========================================
        Single-Cell QC Metrics Processing
        ========================================
        Num input cells      : ${params.num_input_cells}
        Data directory       : ${params.data_dir}
        File prefix          : ${dragen_file_prefix}
        Sample ID            : ${sample_id}
        Output basename      : ${output_basename}
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
    sample_id_ch = Channel.value(sample_id)
    output_basename_ch = Channel.value(output_basename)
    metrics_ch = Channel.fromPath(metrics_path, checkIfExists: true)
    barcode_summary_ch = Channel.fromPath(barcode_summary_path, checkIfExists: true)
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
            .combine(output_basename_ch)
        
        EXTRACT_CRISPR_FEATURES(crispr_input_ch)
        
        EXTRACT_CRISPR_FEATURES.out.crispr_adata.view { "Generated CRISPR h5ad: $it" }
        
        // Step 2: Run CRISPAT guide assignment using the h5ad file
        guide_input_ch = EXTRACT_CRISPR_FEATURES.out.crispr_adata
            .combine(output_basename_ch)
        
        GUIDE_ASSIGNMENT(guide_input_ch)
        
        // Use guide assignment output for metrics processing
        assignments_ch = GUIDE_ASSIGNMENT.out.assignments
        
        GUIDE_ASSIGNMENT.out.assignments.view { "Generated guide assignments: $it" }
    } else {
        log.info "Skipping guide assignment"
        // Create empty channel for assignments
        assignments_ch = Channel.fromPath('NO_FILE')
    }
    
    // Combine all inputs for report generation
    input_ch = num_cells_ch
        .combine(metrics_ch)
        .combine(barcode_summary_ch)
        .combine(sample_id_ch)
        .combine(output_basename_ch)
        .combine(assignments_ch)
    
    // Generate report data
    GENERATE_REPORT_DATA(input_ch)
    
    // Emit results
    GENERATE_REPORT_DATA.out.results.view { "Generated output: $it" }
}
