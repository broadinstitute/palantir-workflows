/*
 * Shared pipeline engine, called by both entrypoints (main.nf for --fastq_list-described
 * runs, main_simple.nf for flat-FASTQ single-subsample runs). Everything here is agnostic
 * to how `subsample_info` was built -- it only cares that it's a list of
 * [rgsm, feature_rgids, hashing_rgids, fastq_files] maps, one per subsample.
 */

include { DRAGEN_SCRNA } from '../modules/dragen_scrna'
include { GENERATE_SUBSAMPLE_QC } from '../modules/generate_subsample_qc'
include { GENERATE_SUPERSAMPLE_QC } from '../modules/generate_supersample_qc'
include { CRISPAT_GUIDE_ASSIGNMENT } from '../modules/crispat_guide_assignment'
include { PURITY_BASED_GUIDE_ASSIGNMENT } from '../modules/purity_based_guide_assignment'
include { CONCATENATE } from '../modules/concatenate'

// Every param this file reads, declared here too so Nextflow doesn't warn about "access to
// undefined parameter" -- the entrypoint (main.nf / main_simple.nf) is what actually sets
// these; these declarations just mirror them for this included module's own scope.
//
// NOTE: params set via CLI/--params-file/nextflow_schema.json at session startup ARE shared
// across included modules like this one, but a plain runtime `params.x = ...` assignment made
// by an entrypoint script does NOT propagate across the include boundary -- Nextflow binds
// each script/module's own `params` view independently after session startup. That's why the
// fastq-list path is threaded through explicitly via `take:` below instead of `params.fastq_list`.
params.num_input_cells = null
params.ref_tar = null
params.annotation_file = null
params.supersample_id = null
params.supersample_basename = null
params.min_valid_guides = null
params.max_valid_guides = null
params.scrna_feature_barcode_reference = null
params.scrna_barcode_sequence_list = null
params.scrna_cell_hashing_reference = null
params.run_guide_assignment = true
params.outdir = "out"
params.use_direct_capture_mode = true
params.additional_dragen_args = null

workflow PIPSEQ_CORE {
    take:
    subsample_info   // channel of maps: [rgsm, feature_rgids, hashing_rgids, fastq_files]
    fastq_list_path  // path to the fastq-list CSV DRAGEN itself reads via --fastq-list

    main:
    // --- Business-logic checks that can't be expressed in JSON Schema ---
    // (Entrypoint-specific checks -- e.g. RGTY validation, or FASTQ-list-length pairing --
    // are the caller's responsibility; these are the ones shared by every entrypoint.)

    if (params.min_valid_guides > params.max_valid_guides) {
        log.error "ERROR: --min_valid_guides (${params.min_valid_guides}) must be <= --max_valid_guides (${params.max_valid_guides})"
        exit 1
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
            file(fastq_list_path),
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

    GENERATE_SUBSAMPLE_QC(qc_input_ch)

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
        log.info "Running CRISPR guide assignment (CRISPAT and purity-based)..."

        // Both methods run independently on the same concatenated CRISPR AnnData.
        CRISPAT_GUIDE_ASSIGNMENT(CONCATENATE.out.concatenated_crispr_adata)
        PURITY_BASED_GUIDE_ASSIGNMENT(CONCATENATE.out.concatenated_crispr_adata)

        // Only CRISPAT's assignments feed into GENERATE_SUPERSAMPLE_QC below; the
        // purity-based assignments are published on their own (see purity_ga/) and
        // aren't otherwise consumed by this pipeline.
        crispat_guide_assignments_ch = CRISPAT_GUIDE_ASSIGNMENT.out.guide_assignments
    } else {
        // Use placeholder for guide assignments -- see the NO_* placeholder note above.
        crispat_guide_assignments_ch = Channel.of(file('NO_FILE'))
    }

    // Generate supersample QC (always runs)
    supersample_qc_input = GENERATE_SUBSAMPLE_QC.out.qc_metrics
        .collect()
        .map { qc_metrics_list -> [qc_metrics_list] }  // Wrap list in tuple to preserve it
        .combine(crispat_guide_assignments_ch)
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

// Write a short manifest describing the output layout. Called from each entrypoint's own
// `workflow.onComplete` block (onComplete is a session-level hook tied to the launched
// script, so it can't live inside this included named workflow itself).
def writeOutputManifest() {
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
          purity_ga/                      Purity-based guide assignment output (only if --run_guide_assignment true)
          supersample_qc/                 Final supersample-level QC report, and (if guide assignment ran)
                                           the guide-assignment distribution plot
          pipeline_info/                  Nextflow execution reports (timeline, report, trace, DAG)

        See README.md in the pipeline repository for parameter and output details.
        """.stripIndent()
}
