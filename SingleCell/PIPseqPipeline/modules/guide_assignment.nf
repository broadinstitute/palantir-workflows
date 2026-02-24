/*
 * Process module for CRISPR guide assignment using CRISPAT
 */

process GUIDE_ASSIGNMENT {
    tag "guide_assignment"    
    container 'us.gcr.io/broad-dsde-methods/pipseq-qc:latest'    
    publishDir "${params.outdir}/${params.supersample_basename}", mode: 'copy'
    
    input:
    path(crispr_adata)
    
    output:
    path "crispat_ga/", emit: crispat_ga_dir
    path "crispat_ga/poisson_gauss/assignments.csv", emit: guide_assignments
    
    script:
    """
    # Run the CRISPAT guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    run_guide_assignment.py \\
        --crispr-adata ${crispr_adata}
    """
}
