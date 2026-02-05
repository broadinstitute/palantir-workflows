/*
 * Process module for concatenating multiple samples using AnnData
 */

process CONCATENATE {
    tag "Concatenating ${sample_ids.size()} samples"
    publishDir "${params.outdir}/concatenated", mode: 'copy'
    
    input:
    tuple path(matrices),
          path(barcodes),
          path(features),
          val(sample_ids)
    
    output:
    path "concatenated.h5ad", emit: concatenated_adata
    
    script:
    // Convert lists to space-separated strings for passing to Python
    def matrix_files = matrices.join(' ')
    def barcode_files = barcodes.join(' ')
    def feature_files = features.join(' ')
    def sample_id_list = sample_ids.join(',')
    """
    # Concatenate multiple samples using AnnData
    concatenate_samples.py \\
        --matrices ${matrix_files} \\
        --barcodes ${barcode_files} \\
        --features ${feature_files} \\
        --sample-ids ${sample_id_list} \\
        --output concatenated.h5ad
    """
}
