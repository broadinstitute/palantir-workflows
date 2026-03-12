/*
 * Process module for concatenating multiple subsamples using AnnData
 */

process CONCATENATE {
    tag "Concatenating ${subsample_ids.size()} subsamples"
    publishDir "${params.outdir}/${params.supersample_basename}/adata", mode: 'copy'
    
    input:
    tuple path(matrices),
          path(barcodes),
          path(features),
          val(subsample_ids)
    
    output:
    path "${params.supersample_basename}.h5ad", emit: concatenated_adata
    path "${params.supersample_basename}.crispr.h5ad", emit: concatenated_crispr_adata
    
    script:
    // Convert lists to space-separated strings for passing to Python
    def matrix_files = matrices.join(' ')
    def barcode_files = barcodes.join(' ')
    def feature_files = features.join(' ')
    def subsample_id_list = subsample_ids.join(',')
    """
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}
    # Concatenate multiple subsamples using AnnData
    concatenate_samples.py \\
        --matrices ${matrix_files} \\
        --barcodes ${barcode_files} \\
        --features ${feature_files} \\
        --subsample-ids ${subsample_id_list} \\
        --supersample-basename ${params.supersample_basename}
    """
    
    stub:
    def subsample_id_list = subsample_ids.join(',')
    """
    echo "[STUB] Would concatenate ${subsample_ids.size()} subsamples:"
    echo "  Subsamples: ${subsample_id_list}"
    echo "  Supersample basename: ${params.supersample_basename}"
    
    touch ${params.supersample_basename}.h5ad
    touch ${params.supersample_basename}.crispr.h5ad
    """
}
