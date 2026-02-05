/*
 * Process module for concatenating multiple subsamples using AnnData
 */

process CONCATENATE {
    tag "Concatenating ${subsample_ids.size()} subsamples"
    publishDir "${params.outdir}/${params.supersample_basename}/concatenated", mode: 'copy'
    
    input:
    tuple path(matrices),
          path(barcodes),
          path(features),
          val(subsample_ids)
    
    output:
    path "concatenated.h5ad", emit: concatenated_adata
    
    script:
    // Convert lists to space-separated strings for passing to Python
    def matrix_files = matrices.join(' ')
    def barcode_files = barcodes.join(' ')
    def feature_files = features.join(' ')
    def subsample_id_list = subsample_ids.join(',')
    """
    # Concatenate multiple subsamples using AnnData
    concatenate_samples.py \\
        --matrices ${matrix_files} \\
        --barcodes ${barcode_files} \\
        --features ${feature_files} \\
        --sample-ids ${subsample_id_list} \\
        --output concatenated.h5ad
    """
}
