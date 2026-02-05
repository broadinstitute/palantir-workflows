/*
 * Process module for extracting CRISPR Direct Capture features
 */

process EXTRACT_CRISPR_FEATURES {
    tag "${data_filtered_matrix.baseName}"
    publishDir "${params.outdir}/crispr_adata", mode: 'copy'
    
    input:
    tuple path(data_filtered_matrix),
          path(data_filtered_barcodes),
          path(data_filtered_features),
          val(output_basename)
    
    output:
    path "${output_basename}.crispr.h5ad", emit: crispr_adata
    
    script:
    """
    # Extract CRISPR Direct Capture feature types and write to h5ad
    # The script should be in the bin/ directory and will be automatically available
    extract_crispr_features.py \\
        --data-filtered-matrix ${data_filtered_matrix} \\
        --data-filtered-barcodes ${data_filtered_barcodes} \\
        --data-filtered-features ${data_filtered_features} \\
        --output-basename ${output_basename}
    """
}
