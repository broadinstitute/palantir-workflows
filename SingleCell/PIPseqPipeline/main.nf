#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Single-Cell QC Metrics Processing Pipeline
 *
 * This pipeline ingests metrics CSVs and data files,
 * processes them with a Python script, and generates output CSVs.
 */

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

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

// concatenate_samples.py hard-caps the number of CRISPR guides it can process for runtime
// reasons. Checked here (fast, before any DRAGEN job runs) and again after concatenation
// in concatenate_samples.py as a safety net, in case the reference and DRAGEN's actual
// feature output ever disagree.
def MAX_CRISPR_GUIDES = 300

// Help message
def helpMessage() {
    log.info"""
    Usage:
      nextflow run main.nf --num_input_cells <int> --fastq_list <fastq_list.csv> --supersample_id <id> --supersample_basename <name> [options]

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
      --run_guide_assignment     Whether to run the CRISPAT guide-assignment step (default: ${params.run_guide_assignment}).
                                  This only toggles the statistical guide-assignment step -- CRISPR feature
                                  extraction and the concatenated supersample AnnData are always produced.
      --guide_assignment_num_processes  Number of processes to use for guide assignment (default: all available cores)
      --outdir                   Output directory (default: ${params.outdir})
      --help                     Show this help message

    Behavior:
      - Runs DRAGEN scRNA for each subsample
      - Concatenates all subsamples into a supersample AnnData (handles single subsample case automatically) -- always runs
      - Runs CRISPAT guide assignment on the concatenated CRISPR features -- only if --run_guide_assignment is true
      - Per-subsample QC reports are generated in outdir/<supersample_basename>/<subsample_id>/qc/
      - Concatenated AnnData outputs to outdir/<supersample_basename>/adata/
      - Guide assignments are output to outdir/<supersample_basename>/crispat_ga/
    """.stripIndent()
}

// Import modules
include { DRAGEN_SCRNA } from './modules/dragen_scrna'
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

    // Validate required/typed params against nextflow_schema.json (required fields, types,
    // patterns, min/max). This is the single source of truth for "what's required" --
    // don't add hand-rolled `if (!params.x) exit 1` checks here for anything the schema
    // already declares; update nextflow_schema.json instead.
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // --- Business-logic checks that can't be expressed in JSON Schema ---

    if (params.min_valid_guides > params.max_valid_guides) {
        log.error "ERROR: --min_valid_guides (${params.min_valid_guides}) must be <= --max_valid_guides (${params.max_valid_guides})"
        exit 1
    }

    if (params.scrna_feature_barcode_reference) {
        def guide_count = file(params.scrna_feature_barcode_reference).readLines().size() - 1 // minus header row
        if (guide_count > MAX_CRISPR_GUIDES) {
            log.error "ERROR: --scrna_feature_barcode_reference contains ${guide_count} guides, but this pipeline cannot process more than ${MAX_CRISPR_GUIDES} guides right now due to runtime constraints."
            exit 1
        }
    }

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

    log.info "Running DRAGEN scRNA for each subsample..."

    // Prepare DRAGEN inputs. Nextflow path inputs can't be cleanly optional, so a 'NO_*'
    // placeholder filename stands in for "not provided" -- DRAGEN_SCRNA and its stub check
    // `.name != 'NO_*'` to tell a real reference apart from this placeholder.
    dragen_input_ch = subsample_info.map { info ->
        tuple(
            [
                subsample_id: info.rgsm,
                feature_barcode_groups: info.feature_rgids,
                hto_barcode_groups: info.hashing_rgids,
                use_direct_capture_mode: params.use_direct_capture_mode,
                additional_dragen_args: params.additional_dragen_args ?: ''
            ],
            file(params.ref_tar),
            file(params.fastq_list),
            file(params.annotation_file),
            params.scrna_feature_barcode_reference ? file(params.scrna_feature_barcode_reference) : file('NO_FEATURE_BARCODE_REF'),
            params.scrna_barcode_sequence_list ? file(params.scrna_barcode_sequence_list) : file('NO_BARCODE_SEQ_LIST'),
            params.scrna_cell_hashing_reference ? file(params.scrna_cell_hashing_reference) : file('NO_CELL_HASHING_REF'),
            info.fastq_files
        )
    }

    // Run DRAGEN
    DRAGEN_SCRNA(dragen_input_ch)

    // Rejoin DRAGEN's named per-file-type outputs by subsample_id. Each of these channels
    // already carries (subsample_id, file) pairs (see modules/dragen_scrna.nf), so this
    // replaces the previous approach of flattening the whole output glob and re-deriving
    // subsample_id/file-type from filenames.
    all_subsamples = DRAGEN_SCRNA.out.metrics
        .join(DRAGEN_SCRNA.out.barcode_summary)
        .join(DRAGEN_SCRNA.out.matrix)
        .join(DRAGEN_SCRNA.out.barcodes)
        .join(DRAGEN_SCRNA.out.features)
        .map { subsample_id, metrics, barcode_summary, matrix, barcodes, features ->
            [
                subsample_id: subsample_id,
                metrics: metrics,
                barcode_summary: barcode_summary,
                matrix: matrix,
                barcodes: barcodes,
                features: features
            ]
        }

    // Generate per-subsample QC reports
    qc_input_ch = all_subsamples.map { s ->
        tuple(
            [
                subsample_id: s.subsample_id,
                supersample_id: params.supersample_id,
                num_input_cells: params.num_input_cells
            ],
            s.metrics,
            s.barcode_summary
        )
    }

    GENERATE_REPORT_DATA(qc_input_ch)

    log.info "Concatenating subsamples into supersample AnnData..."

    // Collect all subsample data for concatenation
    concatenate_input_ch = all_subsamples
        .toList()
        .map { subsamples ->
            tuple(
                subsamples.collect { it.matrix },
                subsamples.collect { it.barcodes },
                subsamples.collect { it.features },
                subsamples.collect { it.subsample_id }
            )
        }

    CONCATENATE(concatenate_input_ch)

    if (params.run_guide_assignment) {
        log.info "Running CRISPR guide assignment..."

        GUIDE_ASSIGNMENT(CONCATENATE.out.concatenated_crispr_adata)

        // Set guide assignments channel
        guide_assignments_ch = GUIDE_ASSIGNMENT.out.guide_assignments
    } else {
        // Use placeholder for guide assignments -- see the NO_* placeholder note above.
        guide_assignments_ch = Channel.of(file('NO_FILE'))
    }

    // Generate supersample QC (always runs)
    supersample_qc_input = GENERATE_REPORT_DATA.out.qc_metrics
        .collect()
        .map { qc_metrics_list -> [qc_metrics_list] }  // Wrap list in tuple to preserve it
        .combine(guide_assignments_ch)
        .map { qc_metrics_list, guide_assignments ->
            // qc_metrics_list is the collected list of qc files
            // guide_assignments is the guide assignments file (or NO_FILE)
            tuple(
                [
                    num_input_cells: params.num_input_cells,
                    supersample_basename: params.supersample_basename,
                    supersample_id: params.supersample_id,
                    min_valid_guides: params.min_valid_guides,
                    max_valid_guides: params.max_valid_guides
                ],
                qc_metrics_list,
                guide_assignments
            )
        }

    GENERATE_SUPERSAMPLE_QC(supersample_qc_input)
}

// Write a short manifest describing the output layout once the run finishes successfully.
workflow.onComplete {
    if (workflow.success) {
        def manifest = file("${params.outdir}/${params.supersample_basename}/README.txt")
        manifest.text = """
            Output layout for supersample '${params.supersample_id}' (${params.supersample_basename}):

              <subsample_id>/dragen_output/   Raw DRAGEN scRNA outputs for that subsample (metrics, barcode
                                               summary, filtered matrix/barcodes/features, and any other files
                                               DRAGEN produced for it)
              <subsample_id>/logs/            DRAGEN logs for that subsample
              <subsample_id>/qc/              Per-subsample QC metrics (qc_metrics.tsv, qc_barcode_metrics.tsv)
              adata/                          Concatenated supersample AnnData (<basename>.h5ad) and CRISPR-
                                               features-only subset (<basename>.crispr.h5ad) -- always produced
              crispat_ga/                     CRISPAT guide assignment output (only if --run_guide_assignment true)
              supersample_qc/                 Final supersample-level QC report, and (if guide assignment ran)
                                               the guide-assignment distribution plot
              pipeline_info/                  Nextflow execution reports (timeline, report, trace, DAG)

            See README.md in the pipeline repository for parameter and output details.
            """.stripIndent()
    }
}
