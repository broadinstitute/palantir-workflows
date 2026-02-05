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
      nextflow run main.nf --num_input_cells <int> --samplesheet <samplesheet.csv> [options]

    Required arguments:
      --num_input_cells          Number of input cells (integer)
      --samplesheet              CSV file with columns: data_dir,dragen_file_prefix,sample_id,sample_basename

    Samplesheet format:
      CSV file with columns: data_dir, dragen_file_prefix, sample_id, sample_basename
      See samplesheet.csv.example for an example

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
      - Single sample: Runs EXTRACT_CRISPR_FEATURES then GUIDE_ASSIGNMENT (if enabled)
      - Multiple samples: Runs CONCATENATE on all samples, then GUIDE_ASSIGNMENT (if enabled)
      - Per-sample QC reports are always generated in outdir/<sample_basename>/qc/
      - Guide assignments are output to outdir/ga_crispat/
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

    log.info "Reading samples from samplesheet..."
    
    // Read samplesheet and create channels
    sample_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // Force strings and strip quotes for each field
            def data_dir = "${row.data_dir}".replaceAll(/^['"]|['"]$/, '')
            def dragen_file_prefix = "${row.dragen_file_prefix}".replaceAll(/^['"]|['"]$/, '')
            def sample_id = "${row.sample_id}".replaceAll(/^['"]|['"]$/, '')
            def sample_basename = "${row.sample_basename}".replaceAll(/^['"]|['"]$/, '')
            
            tuple(data_dir, dragen_file_prefix, sample_id, sample_basename)
        }
    
    // Create channels for all files across all samples
    all_samples = sample_ch.map { data_dir, dragen_file_prefix, sample_id, sample_basename ->
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
            sample_id,
            sample_basename
        )
    }
    
    // Generate per-sample QC reports
    qc_input_ch = all_samples.map { metrics, barcode_summary, _matrix, _barcodes, _features, sample_id, sample_basename ->
        tuple(
            params.num_input_cells,
            metrics,
            barcode_summary,
            sample_id,
            sample_basename
        )
    }
    
    GENERATE_REPORT_DATA(qc_input_ch)
    
    if (params.run_guide_assignment) {
        log.info "Running CRISPR feature extraction and guide assignment..."
        
        // Collect all sample data and branch based on count
        all_samples
            .map { it[2..6] }  // matrix, barcodes, features, sample_id, sample_basename
            .toList()
            .map { samples ->
                [
                    samples.collect { it[0] },  // matrices
                    samples.collect { it[1] },  // barcodes
                    samples.collect { it[2] },  // features
                    samples.collect { it[3] },  // sample_ids
                    samples.collect { it[4] },  // sample_basenames
                    samples.size()              // count
                ]
            }
            .branch {
                single: it[5] == 1
                multiple: it[5] > 1
            }
            .set { sample_data }
        
        // Single sample: extract CRISPR features
        single_sample_ch = sample_data.single
            .map { matrices, barcodes, features, _sample_ids, sample_basenames, _count ->
                tuple(matrices[0], barcodes[0], features[0], sample_basenames[0])
            }
        
        EXTRACT_CRISPR_FEATURES(single_sample_ch)
        
        // Multiple samples: concatenate
        multi_sample_ch = sample_data.multiple
            .map { matrices, barcodes, features, sample_ids, _sample_basenames, _count ->
                tuple(matrices, barcodes, features, sample_ids)
            }
        
        CONCATENATE(multi_sample_ch)
        
        // Combine outputs for guide assignment
        guide_input_ch = EXTRACT_CRISPR_FEATURES.out.crispr_adata
            .mix(CONCATENATE.out.concatenated_adata)
        
        GUIDE_ASSIGNMENT(guide_input_ch)
        
        GUIDE_ASSIGNMENT.out.guide_assignments.view { "Generated guide assignments: $it" }
    }
}
