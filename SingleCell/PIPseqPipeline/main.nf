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
      --scrna_feature_barcode_reference  Feature barcode reference file for DRAGEN (optional)
      --scrna_barcode_sequence_list      Barcode sequence list file for DRAGEN (optional)
      --scrna_cell_hashing_reference     Cell hashing reference file for DRAGEN (optional)

    Fastq_list format:
      CSV file with columns: RGID, RGSM, RGTY, Read1File, Read2File
      - RGSM values represent subsample IDs
      - RGTY indicates readgroup type ('expression' or 'feature' or 'hashing')
      - All rows with the same RGSM belong to the same subsample

    Optional arguments:
      --run_guide_assignment     Whether to run CRISPR guide assignment (default: ${params.run_guide_assignment})
      --guide_assignment_num_processes  Number of processes to use for guide assignment (default: all available cores)
      --outdir                   Output directory (default: ${params.outdir})
      --help                     Show this help message
    
    Behavior:
      - Runs DRAGEN scRNA for each subsample
      - Concatenates all subsamples (handles single subsample case automatically)
      - Runs GUIDE_ASSIGNMENT on concatenated CRISPR features (if enabled)
      - Per-subsample QC reports are generated in outdir/<supersample_basename>/<subsample_basename>/qc/
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

    // Validate inputs
    if (!params.num_input_cells) {
        log.error "ERROR: --num_input_cells is required"
        helpMessage()
        exit 1
    }

    if (!params.fastq_list) {
        log.error "ERROR: --fastq_list is required"
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

    if (!params.ref_tar) {
        log.error "ERROR: --ref_tar is required"
        helpMessage()
        exit 1
    }

    if (!params.annotation_file) {
        log.error "ERROR: --annotation_file is required"
        helpMessage()
        exit 1
    }

    if (!params.dragen_container) {
        log.error "ERROR: --dragen_container is required"
        helpMessage()
        exit 1
    }

    if (!params.qc_container) {
        log.error "ERROR: --qc_container is required"
        helpMessage()
        exit 1
    }

    if (!params.concatenate_cpus) {
        log.error "ERROR: --concatenate_cpus is required"
        helpMessage()
        exit 1
    }

    if (!params.concatenate_memory_gb) {
        log.error "ERROR: --concatenate_memory_gb is required"
        helpMessage()
        exit 1
    }

    if (!params.dragen_scratch_tb) {
        log.error "ERROR: --dragen_scratch_tb is required"
        helpMessage()
        exit 1
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
    
    // Prepare DRAGEN inputs
    dragen_input_ch = subsample_info.map { info ->
        tuple(
            info.rgsm,
            file(params.ref_tar),
            file(params.fastq_list),
            file(params.annotation_file),
            params.scrna_feature_barcode_reference ? file(params.scrna_feature_barcode_reference) : file('NO_FEATURE_BARCODE_REF'),
            params.scrna_barcode_sequence_list ? file(params.scrna_barcode_sequence_list) : file('NO_BARCODE_SEQ_LIST'),
            params.scrna_cell_hashing_reference ? file(params.scrna_cell_hashing_reference) : file('NO_CELL_HASHING_REF'),
            info.feature_rgids,
            info.hashing_rgids,
            info.fastq_files,
            params.use_direct_capture_mode
        )
    }
    
    // Run DRAGEN
    DRAGEN_SCRNA(dragen_input_ch)
    
    // Extract DRAGEN outputs and create channel for downstream processing
    // DRAGEN outputs are in <sample_id>/ directory with prefix <sample_id>
    // The output channel emits all files (including subdirectories)
    all_subsamples = DRAGEN_SCRNA.out.output
            .flatten()  // Flatten all output items
            .filter { it.isFile() }  // Keep only files, not directories
            .map { output_file ->
                def file_name = output_file.getName()
                
                // Extract subsample_id from filename
                def subsample_id = file_name.tokenize('.')[0]
                
                // Match specific output files
                if (file_name.endsWith('.scRNA_metrics.csv')) {
                    tuple('metrics', subsample_id, output_file)
                } else if (file_name.endsWith('.scRNA.barcodeSummary.tsv')) {
                    tuple('barcode_summary', subsample_id, output_file)
                } else if (file_name.endsWith('.scRNA.filtered.matrix.mtx.gz')) {
                    tuple('matrix', subsample_id, output_file)
                } else if (file_name.endsWith('.scRNA.filtered.barcodes.tsv.gz')) {
                    tuple('barcodes', subsample_id, output_file)
                } else if (file_name.endsWith('.scRNA.filtered.features.tsv.gz')) {
                    tuple('features', subsample_id, output_file)
                } else {
                    null
                }
            }
            .filter { it != null }
            .groupTuple(by: 1)  // Group by subsample_id
            .map { file_types, subsample_id, files ->
                // Reorganize into expected structure
                def file_map = [file_types, files].transpose().collectEntries()
                tuple(
                    file_map['metrics'],
                    file_map['barcode_summary'],
                    file_map['matrix'],
                    file_map['barcodes'],
                    file_map['features'],
                    subsample_id
                )
            }
    
    // Generate per-subsample QC reports
    qc_input_ch = all_subsamples.map { metrics, barcode_summary, _matrix, _barcodes, _features, subsample_id ->
        tuple(
            params.num_input_cells,
            metrics,
            barcode_summary,
            subsample_id,
            params.supersample_id
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
        
        // Set guide assignments channel
        guide_assignments_ch = GUIDE_ASSIGNMENT.out.guide_assignments
    } else {
        // Use placeholder for guide assignments
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
                params.num_input_cells,
                qc_metrics_list,
                params.supersample_basename,
                params.supersample_id,
                guide_assignments,
                params.min_valid_guides,
                params.max_valid_guides
            )
        }
    
    GENERATE_SUPERSAMPLE_QC(supersample_qc_input)
}
