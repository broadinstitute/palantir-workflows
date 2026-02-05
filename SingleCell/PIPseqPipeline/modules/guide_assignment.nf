/*
 * Process module for CRISPR guide assignment using CRISPAT
 */

process GUIDE_ASSIGNMENT {
    tag "${sample_basename}"
    publishDir "${params.outdir}/${sample_basename}/ga_crispat", mode: 'copy'
    
    input:
    tuple path(crispr_h5ad),
          val(sample_basename)
    
    output:
    path "${sample_basename}.assignments.csv", emit: guide_assignments
    
    script:
    """
    # Run the CRISPAT guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    run_guide_assignment.py \\
        --crispr-h5ad ${crispr_h5ad} \\
        --sample-basename ${sample_basename}
    """
}
