/*
 * Process module for generating supersample-level QC
 */

process GENERATE_SUPERSAMPLE_QC {
    tag "${meta.supersample_basename}"
    publishDir "${params.outdir}/${params.supersample_basename}/supersample_qc", mode: 'copy'
    container "${params.qc_container}"

    input:
    // meta: [num_input_cells, supersample_basename, supersample_id, min_valid_guides, max_valid_guides]
    tuple val(meta),
          path(subsample_qc_files),
          path(guide_assignments)

    output:
    path "${meta.supersample_basename}.guide_assignment_distribution.png", emit: guide_assignment_distribution, optional: true
    //path "${meta.supersample_basename}.sankey.png", emit: sankey, optional: true
    path "${meta.supersample_basename}.supersample_qc_metrics.tsv", emit: supersample_qc_metrics

    script:
    def guide_arg = guide_assignments.name != 'NO_FILE' ? "--guide-assignments ${guide_assignments} --min-valid-guides ${meta.min_valid_guides} --max-valid-guides ${meta.max_valid_guides}" : ""
    """
    set -ex

    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}
    export TORCHINDUCTOR_CACHE_DIR=${workflow.launchDir}

    # Run the supersample QC script
    generate_supersample_qc.py \\
        --num-input-cells ${meta.num_input_cells} \\
        --subsample-qc-files ${subsample_qc_files.join(' ')} \\
        --supersample-basename ${meta.supersample_basename} \\
        --supersample-id ${meta.supersample_id} \\
        ${guide_arg}
    """

    stub:
    def guide_arg = guide_assignments.name != 'NO_FILE' ? "--guide-assignments ${guide_assignments} --min-valid-guides ${meta.min_valid_guides} --max-valid-guides ${meta.max_valid_guides}" : ""
    """
    echo "[STUB] Would generate supersample QC with:"
    echo "  Supersample ID: ${meta.supersample_id}"
    echo "  Supersample basename: ${meta.supersample_basename}"
    echo "  Num subsamples: ${subsample_qc_files.size()}"
    echo "  Subsample QC files: ${subsample_qc_files.join(', ')}"
    echo "  Guide arguments: ${guide_arg}"

    touch ${meta.supersample_basename}.supersample_qc_metrics.tsv

    # Create optional outputs if guides are present
    if [ "${guide_arg}" != "" ]; then
        touch ${meta.supersample_basename}.guide_assignment_distribution.png
        #touch ${meta.supersample_basename}.sankey.png
    fi
    """
}
