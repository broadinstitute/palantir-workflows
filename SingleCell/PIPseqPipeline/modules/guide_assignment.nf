/*
 * Process module for CRISPR guide assignment using CRISPAT
 */

process GUIDE_ASSIGNMENT {
    tag "guide_assignment"
    publishDir "${params.outdir}/ga_crispat", mode: 'copy'
    
    input:
        path(crispr_h5ad)
    
    output:
    path "assignments.csv", emit: guide_assignments
    
    script:
    """
    # Run the CRISPAT guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    run_guide_assignment.py \\
        --crispr-h5ad ${crispr_h5ad}
    """
}
