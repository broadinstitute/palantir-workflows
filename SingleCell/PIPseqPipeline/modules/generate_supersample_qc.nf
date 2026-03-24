/*
 * Process module for generating supersample-level QC
 */

process GENERATE_SUPERSAMPLE_QC {
    tag "${params.supersample_basename}"
    publishDir "${params.outdir}/${params.supersample_basename}/supersample_qc", mode: 'copy'
    container "${params.qc_container}"
    
    input:
    tuple val(num_input_cells),
          path(subsample_qc_files),
          val(supersample_basename),
          val(supersample_id),
          path(guide_assignments),
          val(min_valid_guides),
          val(max_valid_guides)
    
    output:
    path "${params.supersample_basename}.guide_assignment_distribution.png", emit: guide_assignment_distribution
    //path "${params.supersample_basename}.sankey.png", emit: sankey
    path "${params.supersample_basename}.supersample_qc_metrics.tsv", emit: supersample_qc_metrics
    
    script:
    def guide_arg = guide_assignments.name != 'NO_FILE' ? "--guide-assignments ${guide_assignments} --min-valid-guides ${min_valid_guides} --max-valid-guides ${max_valid_guides}" : ""
    """
    set -ex
    
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}
    export TORCHINDUCTOR_CACHE_DIR=${workflow.launchDir}
    
    # Run the supersample QC script
    generate_supersample_qc.py \\
        --num-input-cells ${num_input_cells} \\
        --subsample-qc-files ${subsample_qc_files.join(' ')} \\
        --supersample-basename ${supersample_basename} \\
        --supersample-id ${supersample_id} \\
        ${guide_arg}
    """
    
    stub:
    def has_guides = guide_assignments.name != 'NO_FILE'
    """
    echo "[STUB] Would generate supersample QC with:"
    echo "  Supersample ID: ${supersample_id}"
    echo "  Supersample basename: ${supersample_basename}"
    echo "  Num subsamples: ${subsample_qc_files.size()}"
    echo "  Has guide assignments: ${has_guides}"
    
    touch ${params.supersample_basename}.supersample_qc_metrics.tsv
    
    # Create optional outputs if guides are present
    if [ "${has_guides}" = "true" ]; then
        touch ${params.supersample_basename}.guide_assignment_distribution.png
        #touch ${params.supersample_basename}.sankey.png
    fi
    """
}
