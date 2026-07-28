/*
 * Process module for purity-based CRISPR guide assignment
 */

process PURITY_BASED_GUIDE_ASSIGNMENT {
    tag "purity_based_guide_assignment"
    publishDir "${params.outdir}/${params.supersample_basename}/purity_ga", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_purity_based_guide_assignment
    memory "${params.memory_gb_purity_based_guide_assignment}.GB"

    input:
    path(crispr_adata)

    output:
    path "${params.supersample_id}.purity_based_guide_assignments.csv", emit: guide_assignments

    script:
    """
    set -ex
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}
    export TORCHINDUCTOR_CACHE_DIR=${workflow.launchDir}

    # Run the purity-based guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    purity_based_guide_assignment.py \\
        --crispr-adata ${crispr_adata} \\
        --supersample-id ${params.supersample_id}
    """

    stub:
    """
    echo "[STUB] Would run purity-based guide assignment on: ${crispr_adata}"

    echo "cell,gRNA,purity,total_count,count_1st,count_2nd" > ${params.supersample_id}.purity_based_guide_assignments.csv
    """
}
