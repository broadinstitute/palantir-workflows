/*
 * Process module for purity-based CRISPR guide assignment
 */

process PURITY_BASED_GUIDE_ASSIGNMENT {
    tag "${subsample_id}"
    publishDir "${params.outdir}/${params.supersample_basename}/${subsample_id}/purity_ga", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_purity_based_guide_assignment
    memory "${params.memory_gb_purity_based_guide_assignment}.GB"

    input:
    tuple val(subsample_id), path(adata)

    output:
    tuple val(subsample_id), path("${subsample_id}.purity_based_guide_assignments.csv"), emit: guide_assignments

    script:
    """
    set -ex
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}
    export TORCHINDUCTOR_CACHE_DIR=${workflow.launchDir}

    # Run the purity-based guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    purity_based_guide_assignment.py \\
        --adata ${adata} \\
        --subsample-id ${subsample_id}
    """

    stub:
    """
    echo "[STUB] Would run purity-based guide assignment on: ${adata} (subsample ${subsample_id})"

    echo "cell,gRNA,purity_1st_vs_2nd,total_count,count_1st,count_2nd" > ${subsample_id}.purity_based_guide_assignments.csv
    """
}
