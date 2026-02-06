/*
 * Process module for CRISPR guide assignment using CRISPAT
 */

process GUIDE_ASSIGNMENT {
    tag "guide_assignment"
    publishDir "${params.outdir}/${params.supersample_basename}/crispat_ga", mode: 'copy'
    
    input:
    path(crispr_adata)
    
    output:
    path "assignments.csv", emit: guide_assignments
    
    script:
    """
    # Run the CRISPAT guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    run_guide_assignment.py \\
        --crispr-adata ${crispr_adata}
    """
}
