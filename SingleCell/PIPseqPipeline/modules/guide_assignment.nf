/*
 * Process module for CRISPR guide assignment using CRISPAT
 */

process GUIDE_ASSIGNMENT {
    tag "${crispr_h5ad.baseName}"
    publishDir "${params.outdir}/ga_crispat", mode: 'copy'
    
    input:
    tuple path(crispr_h5ad),
          path(output_basename)
    
    output:
    path "${output_basename}.assignments.csv", emit: guide_assignments
    path "*.log", emit: logs
    
    script:
    """
    # Run the CRISPAT guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    run_guide_assignment.py \\
        --crispr-h5ad ${crispr_h5ad} \\
        --output-basename ${output_basename} \\
        2>&1 | tee ${output_basename}.guide_assignment.log
    """
}
