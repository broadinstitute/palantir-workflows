/*
 * Process module for metrics processing
 */

process GENERATE_REPORT_DATA {
    tag "${sample_basename}"
    publishDir "${params.outdir}/${sample_basename}/qc", mode: 'copy'
    
    input:
    tuple val(num_input_cells), 
          path(scrna_metrics), 
          path(barcode_summary), 
          val(sample_id),
          val(sample_basename),
          path(guide_assignments)
    
    output:
    path "${sample_basename}.qc_*", emit: qc_files
    
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
        --sample-basename ${sample_basename}
    """
}
