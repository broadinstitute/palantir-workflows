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

    Samplesheet format:
      CSV file with columns: data_dir, dragen_file_prefix, subsample_id, subsample_basename, supersample_id, supersample_basename
      See samplesheet.csv.example for an example
      - Each row represents a subsample
      - All rows should have the same supersample_id and supersample_basename (one supersample per run)

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
      - Single subsample: Runs EXTRACT_CRISPR_FEATURES then GUIDE_ASSIGNMENT (if enabled)
      - Multiple subsamples: Runs CONCATENATE on all subsamples, then GUIDE_ASSIGNMENT (if enabled)
      - Per-subsample QC reports are always generated in outdir/<supersample_basename>/<subsample_basename>/qc/
      - Guide assignments are output to outdir/<supersample_basename>/crispat_ga/
    """.stripIndent()
}

// Import modules
include { GENERATE_REPORT_DATA } from './modules/generate_report_data'
include { EXTRACT_CRISPR_FEATURES } from './modules/extract_crispr_features'
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
        
        // Collect all subsample data with count
        collected_data = all_subsamples
            .map { it[2..6] }  // matrix, barcodes, features, subsample_id, subsample_basename
            .toList()
            .map { subsamples ->
                [
                    subsamples.collect { it[0] },  // matrices
                    subsamples.collect { it[1] },  // barcodes
                    subsamples.collect { it[2] },  // features
                    subsamples.collect { it[3] },  // subsample_ids
                    subsamples.collect { it[4] },  // subsample_basenames
                    subsamples.size()              // count
                ]
            }
        
        // Single subsample: extract CRISPR features
        single_subsample_ch = collected_data
            .filter { matrices, barcodes, features, _subsample_ids, _subsample_basenames, count -> count == 1 }
            .map { matrices, barcodes, features, _subsample_ids, subsample_basenames, _count ->
                tuple(matrices[0], barcodes[0], features[0], subsample_basenames[0])
            }
        
        EXTRACT_CRISPR_FEATURES(single_subsample_ch)
        
        // Multiple subsamples: concatenate
        multi_subsample_ch = collected_data
            .filter { matrices, barcodes, features, _subsample_ids, _subsample_basenames, count -> count > 1 }
            .map { matrices, barcodes, features, subsample_ids, _subsample_basenames, _count ->
                tuple(matrices, barcodes, features, subsample_ids)
            }
        
        CONCATENATE(multi_subsample_ch)
        
        // Combine outputs for guide assignment
        guide_input_ch = EXTRACT_CRISPR_FEATURES.out.crispr_adata
            .mix(CONCATENATE.out.concatenated_adata)
        
        GUIDE_ASSIGNMENT(guide_input_ch)
        
        GUIDE_ASSIGNMENT.out.guide_assignments.view { "Generated guide assignments: $it" }
    }
}
