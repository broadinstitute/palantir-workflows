/*
 * Process module for metrics processing
 */

process GENERATE_REPORT_DATA {
    tag "${subsample_basename}"
    publishDir "${params.outdir}/${params.supersample_basename}/${subsample_basename}/qc", mode: 'copy'
    container "${params.qc_container}"
    
    input:
    tuple val(num_input_cells), 
          path(scrna_metrics), 
          path(barcode_summary), 
          val(subsample_id),
          val(subsample_basename)    
    output:
    path "${subsample_basename}.qc_metrics.tsv", emit: qc_metrics
    path "${subsample_basename}.qc_barcode_metrics.tsv", emit: qc_barcode_metrics
    
    script:
    """
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}

    # Run the Python processing script
    # The script should be in the bin/ directory and will be automatically available
    generate_report_data.py \\
        --num-input-cells ${num_input_cells} \\
        --scrna-metrics ${scrna_metrics} \\
        --barcode-summary ${barcode_summary} \\
        --sample-id ${subsample_id} \\
        --sample-basename ${subsample_basename}
    """
    
    stub:
    """
    echo "[STUB] Would generate report data with:"
    echo "  Sample ID: ${subsample_id}"
    echo "  Sample basename: ${subsample_basename}"
    echo "  Num input cells: ${num_input_cells}"
    
    touch ${subsample_basename}.qc_metrics.tsv
    touch ${subsample_basename}.qc_barcode_metrics.tsv
    """
}
