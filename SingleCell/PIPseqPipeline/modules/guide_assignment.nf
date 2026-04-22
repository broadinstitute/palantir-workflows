/*
 * Process module for CRISPR guide assignment using CRISPAT
 */

process GUIDE_ASSIGNMENT {
    tag "guide_assignment"
    publishDir "${params.outdir}/${params.supersample_basename}", mode: 'copy'
    container "${params.qc_container}"
    
    input:
    path(crispr_adata)
    
    output:
    path "crispat_ga/", emit: crispat_ga_dir
    path "crispat_ga/poisson_gauss/assignments.csv", emit: guide_assignments
    
    script:
    def num_processes_arg = params.guide_assignment_num_processes != null ? "--num-processes ${params.guide_assignment_num_processes}" : ""
    """
    set -ex
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}
    export TORCHINDUCTOR_CACHE_DIR=${workflow.launchDir}
    
    # Run the CRISPAT guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    run_guide_assignment.py \\
        --crispr-adata ${crispr_adata} \\
        ${num_processes_arg}
    """
    
    stub:
    """
    echo "[STUB] Would run guide assignment on: ${crispr_adata}"
    
    mkdir -p crispat_ga/poisson_gauss
    touch crispat_ga/poisson_gauss/assignments.csv
    echo "guide_id,cell_id,assignment_score" > crispat_ga/poisson_gauss/assignments.csv
    """
}
