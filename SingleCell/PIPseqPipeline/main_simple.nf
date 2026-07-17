#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Single-Cell QC Metrics Processing Pipeline -- simple single-subsample entrypoint.
 *
 * For one-off runs with a single subsample, where hand-writing a --fastq_list CSV (see
 * main.nf) is unnecessary friction. Takes flat expression/feature/hashing FASTQ lists
 * instead; subsample_id is set to --supersample_id. Both entrypoints share the same
 * engine -- see workflows/pipseq_core.nf.
 */

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { PIPSEQ_CORE; writeOutputManifest } from './workflows/pipseq_core'

// Define parameters
params.num_input_cells = null          // Number of input cells (integer)
params.expression_r1_fastqs = null     // List of expression R1 FASTQ files (required)
params.expression_r2_fastqs = null     // List of expression R2 FASTQ files (required)
params.feature_r1_fastqs = null        // List of feature (CRISPR/feature-barcode) R1 FASTQ files (optional)
params.feature_r2_fastqs = null        // List of feature R2 FASTQ files (optional)
params.hashing_r1_fastqs = null        // List of cell-hashing R1 FASTQ files (optional)
params.hashing_r2_fastqs = null        // List of cell-hashing R2 FASTQ files (optional)
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

// Help message
def helpMessage() {
    log.info"""
    Usage:
      nextflow run main_simple.nf --num_input_cells <int> \\
        --expression_r1_fastqs <r1.fastq.gz,...> --expression_r2_fastqs <r2.fastq.gz,...> \\
        --supersample_id <id> --supersample_basename <name> [options]

    For multi-subsample production runs, see main.nf instead -- it takes a --fastq_list CSV
    describing potentially many subsamples per run.

    This entrypoint always describes exactly one subsample; its subsample_id is set to
    --supersample_id.

    Required arguments:
      --num_input_cells          Number of input cells (integer)
      --expression_r1_fastqs     Expression R1 FASTQ file(s)
      --expression_r2_fastqs     Expression R2 FASTQ file(s), same count as --expression_r1_fastqs
      --supersample_id           Supersample identifier (used as this run's subsample_id too)
      --supersample_basename     Supersample basename for output organization
      --min_valid_guides         Minimum number of valid guides for guide assignment QC (integer)
      --max_valid_guides         Maximum number of valid guides for guide assignment QC (integer)
      --ref_tar                  DRAGEN reference genome tar file
      --annotation_file          Gene annotation file for DRAGEN
      --dragen_container         Container image for DRAGEN execution
      --qc_container             Container image for QC processing

    Optional DRAGEN arguments:
      --feature_r1_fastqs / --feature_r2_fastqs  Feature (CRISPR/feature-barcode) FASTQ file(s)
                                                   (requires --scrna_feature_barcode_reference)
      --hashing_r1_fastqs / --hashing_r2_fastqs   Cell-hashing FASTQ file(s)
                                                   (requires --scrna_cell_hashing_reference)
      --use_direct_capture_mode          Whether to use direct capture mode in DRAGEN (default: ${params.use_direct_capture_mode})
      --scrna_feature_barcode_reference  Feature barcode reference file for DRAGEN (required if feature FASTQs are given)
      --scrna_barcode_sequence_list      Barcode sequence list file for DRAGEN (optional)
      --scrna_cell_hashing_reference     Cell hashing reference file for DRAGEN (required if hashing FASTQs are given)
      --additional_dragen_args           Additional arguments to pass to DRAGEN command line (optional string)

    Optional arguments:
      --run_guide_assignment     Whether to run guide assignment (default: ${params.run_guide_assignment}). Toggles both
                                  the CRISPAT and purity-based guide-assignment methods together -- CRISPR feature
                                  extraction and the concatenated supersample AnnData are always produced regardless.
      --guide_assignment_num_processes  Number of processes to use for guide assignment (default: all available cores)
      --outdir                   Output directory (default: ${params.outdir})
      --help                     Show this help message
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

    // Validate required/typed params against nextflow_schema_simple.json. Same rules as
    // main.nf: don't add hand-rolled `if (!params.x) exit 1` checks for anything the schema
    // can express -- update nextflow_schema_simple.json instead.
    validateParameters(parameters_schema: 'nextflow_schema_simple.json')
    log.info paramsSummaryLog(workflow)

    // --- Business-logic checks that can't be expressed in JSON Schema ---

    if (params.expression_r1_fastqs.size() != params.expression_r2_fastqs.size()) {
        log.error "ERROR: --expression_r1_fastqs and --expression_r2_fastqs must have the same number of files " +
            "(got ${params.expression_r1_fastqs.size()} and ${params.expression_r2_fastqs.size()})."
        exit 1
    }

    def n_feature_r1 = params.feature_r1_fastqs ? params.feature_r1_fastqs.size() : 0
    def n_feature_r2 = params.feature_r2_fastqs ? params.feature_r2_fastqs.size() : 0
    if (n_feature_r1 != n_feature_r2) {
        log.error "ERROR: --feature_r1_fastqs and --feature_r2_fastqs must have the same number of files " +
            "(got ${n_feature_r1} and ${n_feature_r2})."
        exit 1
    }
    if (n_feature_r1 > 0 && !params.scrna_feature_barcode_reference) {
        log.error "ERROR: feature FASTQs were provided, but --scrna_feature_barcode_reference was not."
        exit 1
    }

    def n_hashing_r1 = params.hashing_r1_fastqs ? params.hashing_r1_fastqs.size() : 0
    def n_hashing_r2 = params.hashing_r2_fastqs ? params.hashing_r2_fastqs.size() : 0
    if (n_hashing_r1 != n_hashing_r2) {
        log.error "ERROR: --hashing_r1_fastqs and --hashing_r2_fastqs must have the same number of files " +
            "(got ${n_hashing_r1} and ${n_hashing_r2})."
        exit 1
    }
    if (n_hashing_r1 > 0 && !params.scrna_cell_hashing_reference) {
        log.error "ERROR: hashing FASTQs were provided, but --scrna_cell_hashing_reference was not."
        exit 1
    }

    // Build one synthetic fastq_list row per (R1, R2) pair, across all three categories,
    // all sharing RGSM = supersample_id (this entrypoint always describes one subsample).
    def lane = 0
    def buildRows = { category, r1_files, r2_files ->
        r1_files.indices.collect { i ->
            lane++
            [
                RGID: "${category}_${i + 1}",
                RGSM: params.supersample_id,
                RGLB: 'Lib1',
                Lane: lane,
                Read1File: file(r1_files[i]),
                Read2File: file(r2_files[i]),
                RGTY: category
            ]
        }
    }

    def rows = buildRows('expression', params.expression_r1_fastqs, params.expression_r2_fastqs)
    if (n_feature_r1 > 0) {
        rows += buildRows('feature', params.feature_r1_fastqs, params.feature_r2_fastqs)
    }
    if (n_hashing_r1 > 0) {
        rows += buildRows('hashing', params.hashing_r1_fastqs, params.hashing_r2_fastqs)
    }

    // Synthesize a DRAGEN-compatible fastq-list CSV (DRAGEN's --fastq-list requires RGID,
    // RGSM, RGLB, Lane, Read1File, Read2File; RGTY is a custom RG* tag DRAGEN passes through
    // -- same column set main.nf's --fastq_list CSVs already use) and pass it to PIPSEQ_CORE
    // explicitly (a plain params.fastq_list assignment here would NOT be visible inside
    // workflows/pipseq_core.nf -- see the note in that file).
    def workdir_path = file(workflow.workDir)
    workdir_path.mkdirs()
    def synthesized_fastq_list = workdir_path.resolve("main_simple_fastq_list_${params.supersample_id}.csv")
    synthesized_fastq_list.text = (
        ['RGID,RGSM,RGLB,Lane,Read1File,Read2File,RGTY'] +
        rows.collect { r -> "${r.RGID},${r.RGSM},${r.RGLB},${r.Lane},${r.Read1File},${r.Read2File},${r.RGTY}" }
    ).join('\n') + '\n'

    def feature_rgids = rows.findAll { it.RGTY == 'feature' }.collect { it.RGID }.join(',')
    def hashing_rgids = rows.findAll { it.RGTY == 'hashing' }.collect { it.RGID }.join(',')
    def fastq_files = rows.collectMany { [it.Read1File, it.Read2File] }.unique()

    subsample_info = Channel.of([
        rgsm: params.supersample_id,
        feature_rgids: feature_rgids,
        hashing_rgids: hashing_rgids,
        fastq_files: fastq_files
    ])

    PIPSEQ_CORE(subsample_info, synthesized_fastq_list)
}

workflow.onComplete {
    if (workflow.success) writeOutputManifest()
}
