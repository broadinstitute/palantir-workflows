/*
 * Process module for metrics processing
 */

process PROCESS_METRICS {
    tag "${scrna_metrics.baseName}"
    publishDir "${params.outdir}/processed", mode: 'copy'
    
    input:
    tuple val(num_input_cells), 
          path(scrna_metrics), 
          path(data_filtered_matrix), 
          path(data_filtered_barcodes), 
          path(data_filtered_features),
          path(assignments)
    
    output:
    path "*.output.csv", emit: results
    path "*.log", emit: logs
    
    script:
    def assignments_arg = assignments.name != 'NO_FILE' ? "--assignments ${assignments}" : ""
    """
    # Run the Python processing script
    # The script should be in the bin/ directory and will be automatically available
    process_metrics.py \\
        --num-input-cells ${num_input_cells} \\
        --scrna-metrics ${scrna_metrics} \\
        --data-filtered-matrix ${data_filtered_matrix} \\
        --data-filtered-barcodes ${data_filtered_barcodes} \\
        --data-filtered-features ${data_filtered_features} \\
        ${assignments_arg} \\
        --output ${scrna_metrics.baseName}.output.csv \\
        2>&1 | tee ${scrna_metrics.baseName}.log
    """
}
