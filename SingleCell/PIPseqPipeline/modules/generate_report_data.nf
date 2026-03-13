/*
 * Process module for metrics processing
 */

process GENERATE_REPORT_DATA {
    tag "${subsample_id}"
    publishDir "${params.outdir}/${params.supersample_basename}/${subsample_id}/qc", mode: 'copy'
    container "${params.qc_container}"
    
    input:
    tuple val(num_input_cells), 
          path(scrna_metrics), 
          path(barcode_summary), 
          val(subsample_id),
          val(supersample_id)    
    output:
    path "${subsample_id}.qc_metrics.tsv", emit: qc_metrics
    path "${subsample_id}.qc_barcode_metrics.tsv", emit: qc_barcode_metrics
    
    script:
    """
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}
    export TORCHINDUCTOR_CACHE_DIR=${workflow.launchDir}

    # Run the Python processing script
    # The script should be in the bin/ directory and will be automatically available
    generate_report_data.py \\
        --num-input-cells ${num_input_cells} \\
        --scrna-metrics ${scrna_metrics} \\
        --barcode-summary ${barcode_summary} \\
        --sample-id ${subsample_id} \\
        --supersample-id ${supersample_id}
    """
    
    stub:
    """
    echo "[STUB] Would generate report data with:"
    echo "  Sample ID: ${subsample_id}"
    echo "  Sample basename: ${subsample_id}"
    echo "  Supersample ID: ${supersample_id}"
    echo "  Num input cells: ${num_input_cells}"
    
    touch ${subsample_id}.qc_metrics.tsv
    touch ${subsample_id}.qc_barcode_metrics.tsv
    """
}
