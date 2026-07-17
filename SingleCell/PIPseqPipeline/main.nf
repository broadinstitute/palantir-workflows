#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Single-Cell QC Metrics Processing Pipeline
 *
 * Production entrypoint: describes every subsample via a --fastq_list CSV (potentially
 * many subsamples per run). For one-off single-subsample runs where hand-writing that CSV
 * is unnecessary friction, see main_simple.nf, which takes flat expression/feature/hashing
 * FASTQ params instead. Both entrypoints share the same engine -- see workflows/pipseq_core.nf.
 */

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { PIPSEQ_CORE; writeOutputManifest } from './workflows/pipseq_core'

// Define parameters
params.num_input_cells = null          // Number of input cells (integer)
params.fastq_list = null               // CSV file with FASTQ information (RGID, RGSM, RGTY, Read1File, Read2File)
params.ref_tar = null                  // DRAGEN reference tar file
params.annotation_file = null          // Gene annotation file for DRAGEN
params.scrna_feature_barcode_reference = null  // Feature barcode reference for DRAGEN
params.scrna_barcode_sequence_list = null      // Optional barcode sequence list for DRAGEN
params.scrna_cell_hashing_reference = null     // Optional cell hashing reference for DRAGEN
params.run_guide_assignment = true     // Whether to run guide assignment
params.outdir = "out"              // Output directory
params.help = false
params.dragen_container = null         // DRAGEN container image
params.qc_container = null             // QC container image
params.use_direct_capture_mode = true       // Whether to use direct capture mode in DRAGEN
params.guide_assignment_num_processes = null  // Number of processes for guide assignment (default: all available cores)
params.additional_dragen_args = null        // Optional additional arguments to pass to DRAGEN command line

// Recognized RGTY values in --fastq_list (see helpMessage() below)
def VALID_RGTY = ['expression', 'feature', 'hashing']

// Help message
def helpMessage() {
    log.info"""
    Usage:
      nextflow run main.nf --num_input_cells <int> --fastq_list <fastq_list.csv> --supersample_id <id> --supersample_basename <name> [options]

    For one-off single-subsample runs, see main_simple.nf instead -- it takes flat
    expression/feature/hashing FASTQ params and doesn't require writing a fastq_list CSV.

    Required arguments:
      --num_input_cells          Number of input cells (integer)
      --fastq_list               CSV file with FASTQ information (columns: RGID, RGSM, RGTY, Read1File, Read2File)
      --supersample_id           Supersample identifier
      --supersample_basename     Supersample basename for output organization
      --min_valid_guides         Minimum number of valid guides for guide assignment QC (integer)
      --max_valid_guides         Maximum number of valid guides for guide assignment QC (integer)
      --ref_tar                  DRAGEN reference genome tar file
      --annotation_file          Gene annotation file for DRAGEN
      --dragen_container         Container image for DRAGEN execution
      --qc_container             Container image for QC processing

    Optional DRAGEN arguments:
      --use_direct_capture_mode          Whether to use direct capture mode in DRAGEN (default: ${params.use_direct_capture_mode})
      --scrna_feature_barcode_reference  Feature barcode reference file for DRAGEN (required if fastq_list has RGTY=feature rows)
      --scrna_barcode_sequence_list      Barcode sequence list file for DRAGEN (optional)
      --scrna_cell_hashing_reference     Cell hashing reference file for DRAGEN (required if fastq_list has RGTY=hashing rows)
      --additional_dragen_args           Additional arguments to pass to DRAGEN command line (optional string)

    Fastq_list format:
      CSV file with columns: RGID, RGSM, RGTY, Read1File, Read2File
      - RGSM values represent subsample IDs
      - RGTY indicates readgroup type: ${VALID_RGTY.join(', ')}
      - All rows with the same RGSM belong to the same subsample

    Optional arguments:
      --run_guide_assignment     Whether to run guide assignment (default: ${params.run_guide_assignment}). Toggles both
                                  the CRISPAT and purity-based guide-assignment methods together -- CRISPR feature
                                  extraction and the concatenated supersample AnnData are always produced regardless.
      --guide_assignment_num_processes  Number of processes to use for guide assignment (default: all available cores)
      --outdir                   Output directory (default: ${params.outdir})
      --help                     Show this help message

    Behavior:
      - Runs DRAGEN scRNA for each subsample
      - Concatenates all subsamples into a supersample AnnData (handles single subsample case automatically) -- always runs
      - Runs CRISPAT and purity-based guide assignment on the concatenated CRISPR features -- only if --run_guide_assignment is true
      - Per-subsample QC reports are generated in outdir/<supersample_basename>/<subsample_id>/qc/
      - Concatenated AnnData outputs to outdir/<supersample_basename>/adata/
      - CRISPAT guide assignments are output to outdir/<supersample_basename>/crispat_ga/
      - Purity-based guide assignments are output to outdir/<supersample_basename>/purity_ga/
    """.stripIndent()
}

/*
 * Main workflow
 */
workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate required/typed params against nextflow_schema.json (required fields, types,
    // patterns, min/max). This is the single source of truth for "what's required" --
    // don't add hand-rolled `if (!params.x) exit 1` checks here for anything the schema
    // already declares; update nextflow_schema.json instead.
    validateParameters()
    log.info paramsSummaryLog(workflow)

    log.info "Reading subsamples from fastq_list..."

    // Read fastq_list and parse subsample information
    fastq_list_ch = Channel
        .fromPath(params.fastq_list, checkIfExists: true)
        .splitCsv(header: true)
        .filter { row ->
            // Filter out rows with empty or null required fields
            row.RGID && row.RGSM && row.RGTY && row.Read1File && row.Read2File &&
            row.RGID.trim() && row.RGSM.trim() && row.RGTY.trim() &&
            row.Read1File.trim() && row.Read2File.trim()
        }
        .map { row ->
            [
                RGID: row.RGID,
                RGSM: row.RGSM,  // subsample_id
                RGTY: row.RGTY,  // 'expression' or 'feature' or 'hashing'
                Read1File: file(row.Read1File),
                Read2File: file(row.Read2File)
            ]
        }

    // Group by subsample (RGSM) and collect feature RGIDs
    subsample_info = fastq_list_ch
        .toList()
        .flatMap { rows ->
            // Validate RGTY values here (fastq_list is now fully materialized) so a typo
            // fails immediately with a clear message instead of silently producing an
            // empty feature/hashing group, or failing deep inside a DRAGEN job.
            def invalid_rows = rows.findAll { !VALID_RGTY.contains(it.RGTY) }
            if (invalid_rows) {
                log.error "ERROR: Invalid RGTY value(s) in --fastq_list: " +
                    invalid_rows.collect { "RGID=${it.RGID} RGTY=${it.RGTY}" }.join(', ') +
                    ". RGTY must be one of: ${VALID_RGTY.join(', ')}."
                exit 1
            }

            def has_feature_rows = rows.any { it.RGTY == 'feature' }
            def has_hashing_rows = rows.any { it.RGTY == 'hashing' }
            if (has_feature_rows && !params.scrna_feature_barcode_reference) {
                log.error "ERROR: --fastq_list contains RGTY=feature rows, but --scrna_feature_barcode_reference was not provided."
                exit 1
            }
            if (has_hashing_rows && !params.scrna_cell_hashing_reference) {
                log.error "ERROR: --fastq_list contains RGTY=hashing rows, but --scrna_cell_hashing_reference was not provided."
                exit 1
            }

            // Get unique subsamples
            def subsamples = rows.collect { it.RGSM }.unique()

            // For each subsample, collect feature RGIDs and FASTQ files
            subsamples.collect { rgsm ->
                def feature_rgids = rows
                    .findAll { it.RGSM == rgsm && it.RGTY == 'feature' }
                    .collect { it.RGID }
                    .join(',')

                def hashing_rgids = rows
                    .findAll { it.RGSM == rgsm && it.RGTY == 'hashing' }
                    .collect { it.RGID }
                    .join(',')

                // Collect all unique FASTQ files for this subsample
                def fastq_files = rows
                    .findAll { it.RGSM == rgsm }
                    .collectMany { [it.Read1File, it.Read2File] }
                    .unique()

                [
                    rgsm: rgsm,
                    feature_rgids: feature_rgids,
                    hashing_rgids: hashing_rgids,
                    fastq_files: fastq_files
                ]
            }
        }

    PIPSEQ_CORE(subsample_info, file(params.fastq_list))
}

workflow.onComplete {
    if (workflow.success) writeOutputManifest()
}
