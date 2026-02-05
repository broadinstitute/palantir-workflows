/*
 * Process module for CRISPR guide assignment using CRISPAT
 */

process GUIDE_ASSIGNMENT {
    tag "${crispr_h5ad.baseName}"
    publishDir "${params.outdir}/guide_assignment", mode: 'copy'
    
    input:
    tuple path(crispr_h5ad),
          path(scrna_metrics)
    
    output:
    path "assignments.csv", emit: assignments
    path "*.log", emit: logs
    
    script:
    """
    # Run the CRISPAT guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    run_guide_assignment.py \\
        --crispr-h5ad ${crispr_h5ad} \\
        --scrna-metrics ${scrna_metrics} \\
        --output assignments.csv \\
        2>&1 | tee guide_assignment.log
    """
}
