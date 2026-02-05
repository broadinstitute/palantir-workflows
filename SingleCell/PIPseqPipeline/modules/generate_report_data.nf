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
          path(guide_assignments)
    
    output:
    path "${output_basename}.qc_*", emit: qc_files
    path "*.log", emit: logs
    
    script:
    def guide_assignments_arg = guide_assignments.name != 'NO_FILE' ? "--guide-assignments ${guide_assignments}" : ""
    """
    set -euo pipefail
    
    # Run the Python processing script
    # The script should be in the bin/ directory and will be automatically available
    generate_report_data.py \\
        --num-input-cells ${num_input_cells} \\
        --scrna-metrics ${scrna_metrics} \\
        --barcode-summary ${barcode_summary} \\
        --sample-id ${sample_id} \\
        ${guide_assignments_arg} \\
        --output-basename ${output_basename} \\
        2>&1 | tee ${output_basename}.generate_report_data.log
    
    exit_code=\${PIPESTATUS[0]}
    if [ \$exit_code -ne 0 ]; then
        echo "ERROR: generate_report_data.py failed with exit code \$exit_code" >&2
        echo "Check ${output_basename}.generate_report_data.log for details" >&2
        exit \$exit_code
    fi
    """
}
