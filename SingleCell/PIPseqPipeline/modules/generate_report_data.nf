/*
 * Process module for metrics processing
 */

process GENERATE_REPORT_DATA {
    tag "${meta.subsample_id}"
    publishDir "${params.outdir}/${params.supersample_basename}/${meta.subsample_id}/qc", mode: 'copy'
    container "${params.qc_container}"

    input:
    // meta: [subsample_id, supersample_id, num_input_cells]
    tuple val(meta),
          path(scrna_metrics),
          path(barcode_summary)
    output:
    path "${meta.subsample_id}.qc_metrics.tsv", emit: qc_metrics
    path "${meta.subsample_id}.qc_barcode_metrics.tsv", emit: qc_barcode_metrics

    script:
    """
    set -ex
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}
    export TORCHINDUCTOR_CACHE_DIR=${workflow.launchDir}

    # Run the Python processing script
    # The script should be in the bin/ directory and will be automatically available
    generate_report_data.py \\
        --num-input-cells ${meta.num_input_cells} \\
        --scrna-metrics ${scrna_metrics} \\
        --barcode-summary ${barcode_summary} \\
        --sample-id ${meta.subsample_id} \\
        --supersample-id ${meta.supersample_id}
    """

    stub:
    """
    echo "[STUB] Would generate report data with:"
    echo "  Sample ID: ${meta.subsample_id}"
    echo "  Supersample ID: ${meta.supersample_id}"
    echo "  Num input cells: ${meta.num_input_cells}"

    touch ${meta.subsample_id}.qc_metrics.tsv
    touch ${meta.subsample_id}.qc_barcode_metrics.tsv
    """
}
