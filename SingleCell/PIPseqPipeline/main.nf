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
params.samplesheet = null              // CSV file with samples
params.run_guide_assignment = true     // Whether to run guide assignment
params.outdir = "results"              // Output directory
params.help = false

// Help message
def helpMessage() {
    log.info"""
    Usage:
      nextflow run main.nf --num_input_cells <int> --samplesheet <samplesheet.csv> --supersample_id <id> --supersample_basename <name> [options]

    Required arguments:
      --num_input_cells          Number of input cells (integer)
      --samplesheet              CSV file with columns: data_dir,dragen_file_prefix,subsample_id,subsample_basename
      --supersample_id           Supersample identifier
      --supersample_basename     Supersample basename for output organization
      --min_valid_guides         Minimum number of valid guides for guide assignment QC (integer)
      --max_valid_guides         Maximum number of valid guides for guide assignment QC (integer)

    Samplesheet format:
      CSV file with columns: data_dir, dragen_file_prefix, subsample_id, subsample_basename
      See samplesheet.csv.example for an example
      - Each row represents a subsample belonging to the supersample specified in --supersample_id/--supersample_basename

    Expected files in each data_dir:
      <dragen_file_prefix>.scRNA_metrics.csv
      <dragen_file_prefix>.scRNA.barcodeSummary.tsv
      <dragen_file_prefix>.scRNA.filtered.matrix.mtx.gz
      <dragen_file_prefix>.scRNA.filtered.barcodes.tsv.gz
      <dragen_file_prefix>.scRNA.filtered.features.tsv.gz

    Optional arguments:
      --run_guide_assignment     Whether to run CRISPR guide assignment (default: ${params.run_guide_assignment})
      --outdir                   Output directory (default: ${params.outdir})
      --help                     Show this help message
    
    Note: If any parameter values in the samplesheet contain leading zeros (e.g., 001234),
          ensure they are properly formatted as strings in your CSV.
    
    Behavior:
      - Concatenates all subsamples (handles single subsample case automatically)
      - Runs GUIDE_ASSIGNMENT on concatenated CRISPR features (if enabled)
      - Per-subsample QC reports are always generated in outdir/<supersample_basename>/<subsample_basename>/qc/
      - Concatenated AnnData outputs to outdir/<supersample_basename>/adata/
      - Guide assignments are output to outdir/<supersample_basename>/crispat_ga/
    """.stripIndent()
}

// Import modules
include { GENERATE_REPORT_DATA } from './modules/generate_report_data'
include { GENERATE_SUPERSAMPLE_QC } from './modules/generate_supersample_qc'
include { GUIDE_ASSIGNMENT } from './modules/guide_assignment'
include { CONCATENATE } from './modules/concatenate'

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

    if (!params.samplesheet) {
        log.error "ERROR: --samplesheet is required"
        helpMessage()
        exit 1
    }

    if (!params.supersample_id) {
        log.error "ERROR: --supersample_id is required"
        helpMessage()
        exit 1
    }

    if (!params.supersample_basename) {
        log.error "ERROR: --supersample_basename is required"
        helpMessage()
        exit 1
    }

    if (!params.min_valid_guides || !params.max_valid_guides) {
        log.error "ERROR: --min_valid_guides and --max_valid_guides are required"
        helpMessage()
        exit 1
    }

    log.info "Reading subsamples from samplesheet..."
    
    // Read samplesheet and create channels
    subsample_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // Force strings and strip quotes for each field
            def data_dir = "${row.data_dir ?: ''}".replaceAll(/^['"]|['"]$/, '')
            def dragen_file_prefix = "${row.dragen_file_prefix ?: ''}".replaceAll(/^['"]|['"]$/, '')
            def subsample_id = "${row.subsample_id ?: ''}".replaceAll(/^['"]|['"]$/, '')
            def subsample_basename = "${row.subsample_basename ?: ''}".replaceAll(/^['"]|['"]$/, '')
            
            // Validate that all required columns are present
            if (!data_dir || !dragen_file_prefix || !subsample_id || !subsample_basename) {
                error "ERROR: Samplesheet is missing required column(s). Required columns: data_dir, dragen_file_prefix, subsample_id, subsample_basename\nFound columns: ${row.keySet().join(', ')}"
            }
            
            tuple(data_dir, dragen_file_prefix, subsample_id, subsample_basename)
        }
    
    // Create channels for all files across all subsamples
    all_subsamples = subsample_ch.map { data_dir, dragen_file_prefix, subsample_id, subsample_basename ->
        def metrics_path = "${data_dir}/${dragen_file_prefix}.scRNA_metrics.csv"
        def barcode_summary_path = "${data_dir}/${dragen_file_prefix}.scRNA.barcodeSummary.tsv"
        def matrix_path = "${data_dir}/${dragen_file_prefix}.scRNA.filtered.matrix.mtx.gz"
        def barcodes_path = "${data_dir}/${dragen_file_prefix}.scRNA.filtered.barcodes.tsv.gz"
        def features_path = "${data_dir}/${dragen_file_prefix}.scRNA.filtered.features.tsv.gz"
        
        tuple(
            file(metrics_path),
            file(barcode_summary_path),
            file(matrix_path),
            file(barcodes_path),
            file(features_path),
            subsample_id,
            subsample_basename
        )
    }
    
    // Generate per-subsample QC reports
    qc_input_ch = all_subsamples.map { metrics, barcode_summary, _matrix, _barcodes, _features, subsample_id, subsample_basename ->
        tuple(
            params.num_input_cells,
            metrics,
            barcode_summary,
            subsample_id,
            subsample_basename
        )
    }
    
    GENERATE_REPORT_DATA(qc_input_ch)
    
    if (params.run_guide_assignment) {
        log.info "Running CRISPR feature extraction and guide assignment..."
        
        // Collect all subsample data for concatenation
        concatenate_input_ch = all_subsamples
            .toList()
            .map { subsamples ->
                tuple(
                    subsamples.collect { it[2] },  // matrices
                    subsamples.collect { it[3] },  // barcodes
                    subsamples.collect { it[4] },  // features
                    subsamples.collect { it[5] }   // subsample_ids
                )
            }
        
        CONCATENATE(concatenate_input_ch)
        
        GUIDE_ASSIGNMENT(CONCATENATE.out.concatenated_crispr_adata)
        
        GUIDE_ASSIGNMENT.out.guide_assignments.view { "Generated guide assignments: $it" }
        
        // Set guide assignments channel
        guide_assignments_ch = GUIDE_ASSIGNMENT.out.guide_assignments
    } else {
        // Use placeholder for guide assignments
        guide_assignments_ch = Channel.of(file('NO_FILE'))
    }
    
    // Generate supersample QC (always runs)
    supersample_qc_input = GENERATE_REPORT_DATA.out.qc_files
        .collect()
        .combine(guide_assignments_ch)
        .map { qc_files, guide_assignments ->
            tuple(
                params.num_input_cells,
                qc_files,
                params.supersample_basename,
                params.supersample_id,
                guide_assignments,
                params.min_valid_guides,
                params.max_valid_guides
            )
        }
    
    GENERATE_SUPERSAMPLE_QC(supersample_qc_input)
}
