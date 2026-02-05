/*
 * Process module for metrics processing
 */

process GENERATE_REPORT_DATA {
    tag "${subsample_basename}"
    publishDir "${params.outdir}/${params.supersample_basename}/${subsample_basename}/qc", mode: 'copy'
    
    input:
    tuple val(num_input_cells), 
          path(scrna_metrics), 
          path(barcode_summary), 
          val(subsample_id),
          val(subsample_basename)    
    output:
    path "${subsample_basename}.qc_*", emit: qc_files
    
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
