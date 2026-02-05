/*
 * Process module for metrics processing
 */

process GENERATE_REPORT_DATA {
    tag "${scrna_metrics.baseName}"
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    tuple val(num_input_cells), 
          path(scrna_metrics), 
          path(barcode_summary), 
          val(sample_id),
          val(output_basename),
          path(guide_assignments)
    
    output:
    path "${output_basename}.qc_*", emit: qc_files
    path "*.log", emit: logs
    
    script:
    def guide_assignments_arg = guide_assignments.name != 'NO_FILE' ? "--guide-assignments ${guide_assignments}" : ""
    """
    # Run the Python processing script
    # The script should be in the bin/ directory and will be automatically available
    generate_report_data.py \\
        --num-input-cells ${num_input_cells} \\
        --scrna-metrics ${scrna_metrics} \\
        --barcode-summary ${barcode_summary} \\
        --sample-id ${sample_id} \\
        ${guide_assignments_arg} \\
        --output-basename ${output_basename}
    """
}
