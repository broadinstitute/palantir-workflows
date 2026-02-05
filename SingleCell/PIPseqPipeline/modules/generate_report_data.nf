/*
 * Process module for metrics processing
 */

process GENERATE_REPORT_DATA {
    tag "${scrna_metrics.baseName}"
    publishDir "${params.outdir}/processed", mode: 'copy'
    
    input:
    tuple val(num_input_cells), 
          path(scrna_metrics), 
          path(barcode_summary), 
          val(sample_id),
          val(output_basename),
          path(assignments)
    
    output:
    path "*.output.csv", emit: results
    path "*.log", emit: logs
    
    script:
    def assignments_arg = assignments.name != 'NO_FILE' ? "--assignments ${assignments}" : ""
    """
    # Run the Python processing script
    # The script should be in the bin/ directory and will be automatically available
    generate_report_data.py \\
        --num-input-cells ${num_input_cells} \\
        --scrna-metrics ${scrna_metrics} \\
        --barcode-summary ${barcode_summary} \\
        --sample-id ${sample_id} \\
        ${assignments_arg} \\
        --output-basename ${output_basename} \\
        2>&1 | tee ${output_basename}.generate_report_data.log
    """
}
