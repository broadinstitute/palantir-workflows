/*
 * Process module for CRISPR guide assignment using CRISPAT
 */

process CRISPAT_GUIDE_ASSIGNMENT {
    tag "${subsample_id}"
    publishDir "${params.outdir}/${params.supersample_basename}/${subsample_id}", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_crispat_guide_assignment
    memory "${params.memory_gb_crispat_guide_assignment}.GB"

    input:
    tuple val(subsample_id), path(adata)

    output:
    path "crispat_ga/", emit: crispat_ga_dir
    tuple val(subsample_id), path("${subsample_id}.crispat_guide_assignments.csv"), emit: guide_assignments

    script:
    def num_processes_arg = params.guide_assignment_num_processes != null ? "--num-processes ${params.guide_assignment_num_processes}" : ""
    """
    set -ex
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}
    export TORCHINDUCTOR_CACHE_DIR=${workflow.launchDir}

    # Run the CRISPAT guide assignment script
    # The script should be in the bin/ directory and will be automatically available
    run_crispat_guide_assignment.py \\
        --adata ${adata} \\
        --subsample-id ${subsample_id} \\
        ${num_processes_arg}
    """

    stub:
    """
    echo "[STUB] Would run CRISPAT guide assignment on: ${adata} (subsample ${subsample_id})"

    mkdir -p crispat_ga/poisson_gauss
    echo "guide_id,cell_id,assignment_score" > crispat_ga/poisson_gauss/assignments.csv
    echo "cell,gRNA" > ${subsample_id}.crispat_guide_assignments.csv
    """
}
