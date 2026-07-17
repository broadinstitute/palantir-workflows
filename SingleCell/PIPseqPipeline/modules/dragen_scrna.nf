// Builds the DRAGEN command-line fragments that are conditional on optional inputs.
// Shared between the script: and stub: blocks below so they can't silently drift apart --
// previously each block duplicated this logic independently.
def buildDragenOptionalArgs(meta, scrna_feature_barcode_reference, scrna_barcode_sequence_list, scrna_cell_hashing_reference) {
    def has_hashing_ref = scrna_cell_hashing_reference.name != 'NO_CELL_HASHING_REF'
    [
        scrna_barcode_sequence_list_arg: scrna_barcode_sequence_list.name != 'NO_BARCODE_SEQ_LIST' ? "--scrna-barcode-sequence-list ${scrna_barcode_sequence_list}" : "",
        scrna_feature_barcode_reference_arg: scrna_feature_barcode_reference.name != 'NO_FEATURE_BARCODE_REF' ? "--scrna-feature-barcode-reference ${scrna_feature_barcode_reference}" : "",
        scrna_cell_hashing_reference_arg: has_hashing_ref ? "--scrna-cell-hashing-reference ${scrna_cell_hashing_reference}" : "",
        scrna_hto_barcode_groups_arg: has_hashing_ref ? "--scrna-hto-barcode-groups ${meta.hto_barcode_groups}" : "",
        scrna_direct_capture_mode_arg: meta.feature_barcode_groups ? "--scrna-enable-direct-capture-mode ${meta.use_direct_capture_mode}" : "",
        scrna_direct_capture_barcode_groups_arg: meta.feature_barcode_groups ? (meta.use_direct_capture_mode ? "--scrna-direct-capture-barcode-groups ${meta.feature_barcode_groups}" : "--scrna-feature-barcode-groups ${meta.feature_barcode_groups}") : ""
    ]
}

process DRAGEN_SCRNA {
    tag "${meta.subsample_id}"

    container "${params.dragen_container}"
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'fpga2-large'
    // add scratch space for intermediate files
    pod annotation: 'volumes.illumina.com/scratchSize', value: "${params.dragen_scratch_tb}TiB"
    // ICA will upload everything in the "out" folder to cloud storage
    publishDir "${params.outdir}/${params.supersample_basename}/${meta.subsample_id}", mode: 'copy'


    input:
        // meta: [subsample_id, feature_barcode_groups, hto_barcode_groups, use_direct_capture_mode, additional_dragen_args]
        tuple val(meta),
              path(ref_tar),
              path(fastq_list),
              path(annotation_file),
              path(scrna_feature_barcode_reference),
              path(scrna_barcode_sequence_list),
              path(scrna_cell_hashing_reference),
              path(fastq_files)
    output:
        stdout emit: result
        // Full glob: publishes EVERY file DRAGEN writes to dragen_output/, not just the five
        // consumed by downstream processes below. Keep this even though it overlaps with the
        // named emits below -- it's what guarantees nothing DRAGEN produces gets silently
        // dropped just because this pipeline doesn't happen to read it.
        path 'dragen_output/*', emit: output
        path 'logs/*', emit: logs
        // Named, order-independent outputs for downstream wiring (see main.nf's use of .join()).
        tuple val(meta.subsample_id), path("dragen_output/${meta.subsample_id}.scRNA_metrics.csv"), emit: metrics
        tuple val(meta.subsample_id), path("dragen_output/${meta.subsample_id}.scRNA.barcodeSummary.tsv"), emit: barcode_summary
        tuple val(meta.subsample_id), path("dragen_output/${meta.subsample_id}.scRNA.filtered.matrix.mtx.gz"), emit: matrix
        tuple val(meta.subsample_id), path("dragen_output/${meta.subsample_id}.scRNA.filtered.barcodes.tsv.gz"), emit: barcodes
        tuple val(meta.subsample_id), path("dragen_output/${meta.subsample_id}.scRNA.filtered.features.tsv.gz"), emit: features

    script:
        def subsample_id = meta.subsample_id
        def fastqList = fastq_files.collect{ it.toString() }
        def dragenArgs = buildDragenOptionalArgs(meta, scrna_feature_barcode_reference, scrna_barcode_sequence_list, scrna_cell_hashing_reference)
        """
        set -ex
        ### create temporary directory for DRAGEN reference to get extracted to
        mkdir -p /scratch/reference
        tar -C /scratch/reference -xf ${ref_tar}

        ### create output folder
        mkdir -p dragen_output

        ### pre-requisite steps to ensure robust DRAGEN pipeline running
        /opt/edico/bin/dragen_reset
        /opt/edico/bin/dragen --partial-reconfig HMM --ignore-version-check true 2>&1

        ### Actual DRAGEN command
        /opt/edico/bin/dragen --lic-instance-id-location /opt/instance-identity \\
            --annotation-file ${annotation_file} \\
            --fastq-list ${fastq_list} \\
            --output-format CRAM \\
            --umi-source read1 \\
            --scrna-barcode-position 0_7+11_16+20_25+31_38 \\
            --scrna-umi-position 39_41 \\
            --rna-library-type SF \\
            --scrna-enable-pipseq-mode true \\
            --scrna-demux-detect-doublets false \\
            --enable-rna true \\
            --enable-single-cell-rna true \\
            --logging-to-output-dir true \\
            --autodetect-reference-validate true \\
            --enable-bam-indexing true \\
            --enable-map-align-output false \\
            --enable-map-align true \\
            --generate-sa-tags true \\
            --dump-map-align-registers true \\
            --rrna-filter-enable true \\
            --enable-variant-caller false \\
            --enable-cnv false \\
            --enable-sv false \\
            --repeat-genotype-enable false \\
            --vc-enable-profile-stats true \\
            --qc-enable-depth-metrics false \\
            --enable-metrics-json true \\
            --output-file-prefix ${subsample_id} \\
            --output-directory dragen_output \\
            --intermediate-results-dir /scratch \\
            --fastq-list-sample-id ${subsample_id} \\
            --ref-dir /scratch/reference \\
            ${dragenArgs.scrna_direct_capture_mode_arg} \\
            ${dragenArgs.scrna_feature_barcode_reference_arg} \\
            ${dragenArgs.scrna_direct_capture_barcode_groups_arg} \\
            ${dragenArgs.scrna_barcode_sequence_list_arg} \\
            ${dragenArgs.scrna_cell_hashing_reference_arg} \\
            ${dragenArgs.scrna_hto_barcode_groups_arg} \\
            --bin_memory 64424509440 \\
            --bin-split-threshold 32212254720 \\
            --force \\
            -v \\
            ${meta.additional_dragen_args}


        # copy logs
        mkdir -p logs
        cp -rvL /var/log/dragen/* logs/
        rm -f logs/dragen_last_good_run.log

        /opt/edico/bin/dragen_reset
        """

    stub:
        def subsample_id = meta.subsample_id
        def dragenArgs = buildDragenOptionalArgs(meta, scrna_feature_barcode_reference, scrna_barcode_sequence_list, scrna_cell_hashing_reference)
    """
    echo "[STUB] Would run DRAGEN with:"
    echo "  Sample ID: ${subsample_id}"
    echo "  Feature barcode groups: ${meta.feature_barcode_groups}"
    echo "  Feature barcode reference: ${scrna_feature_barcode_reference}"
    echo "  Barcode sequence list: ${scrna_barcode_sequence_list.name != 'NO_BARCODE_SEQ_LIST' ? scrna_barcode_sequence_list : 'Not provided'}"
    echo "  Cell hashing reference: ${scrna_cell_hashing_reference.name != 'NO_CELL_HASHING_REF' ? scrna_cell_hashing_reference : 'Not provided'}"
    echo "  Cell hashing barcode groups: ${meta.hto_barcode_groups}"
    echo "  scrna_barcode_sequence_list_arg: ${dragenArgs.scrna_barcode_sequence_list_arg}"
    echo "  scrna_cell_hashing_reference_arg: ${dragenArgs.scrna_cell_hashing_reference_arg}"
    echo "  scrna_hto_barcode_groups_arg: ${dragenArgs.scrna_hto_barcode_groups_arg}"
    echo "  scrna_direct_capture_barcode_groups_arg: ${dragenArgs.scrna_direct_capture_barcode_groups_arg}"
    echo "  scrna_direct_capture_mode_arg: ${dragenArgs.scrna_direct_capture_mode_arg}"
    echo "  Reference: ${ref_tar}"
    echo "  Annotation: ${annotation_file}"
    echo "  FASTQ files: ${fastq_files.join(', ')}"

    mkdir -p dragen_output
    touch dragen_output/${subsample_id}.scRNA_metrics.csv
    touch dragen_output/${subsample_id}.scRNA.barcodeSummary.tsv
    touch dragen_output/${subsample_id}.scRNA.filtered.matrix.mtx.gz
    touch dragen_output/${subsample_id}.scRNA.filtered.barcodes.tsv.gz
    touch dragen_output/${subsample_id}.scRNA.filtered.features.tsv.gz

    mkdir -p logs
    echo "[STUB] DRAGEN completed successfully" > logs/stub.log
    """
}
