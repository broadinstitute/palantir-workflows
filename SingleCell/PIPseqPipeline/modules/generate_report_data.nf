/*
 * Process module for metrics processing
 */

process GENERATE_REPORT_DATA {
    tag "${subsample_basename}"
    container 'us.gcr.io/broad-dsde-methods/pipseq-qc:latest'
    publishDir "${params.outdir}/${params.supersample_basename}/${subsample_basename}/qc", mode: 'copy'
    
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
    # Run the Python processing script
    # The script should be in the bin/ directory and will be automatically available
    generate_report_data.py \\
        --num-input-cells ${num_input_cells} \\
        --scrna-metrics ${scrna_metrics} \\
        --barcode-summary ${barcode_summary} \\
        --sample-id ${subsample_id} \\
        --sample-basename ${subsample_basename}
    """
}
