#!/usr/bin/env python3
"""
Extract CRISPR features from 10x-style filtered matrix.
"""

import argparse
import pandas as pd
import sys
from pathlib import Path
from anndata import AnnData
from typing import Literal
import scanpy as sc
import crispat
import tempfile
import os
import shutil

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract CRISPR features from 10x-style filtered matrix"
    )
    parser.add_argument(
        "--data-filtered-matrix",
        type=str,
        required=True,
        help="Path to filtered matrix file"
    )
    parser.add_argument(
        "--data-filtered-barcodes",
        type=str,
        required=True,
        help="Path to filtered barcodes file"
    )
    parser.add_argument(
        "--data-filtered-features",
        type=str,
        required=True,
        help="Path to filtered features file"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to output h5ad file"
    )
    return parser.parse_args()

def read_10x_mtx_feature_types(
        path: sc.readwrite.Path | str,
        feature_types: str,
        *,
        var_names: Literal["gene_symbols", "gene_ids"] = "gene_symbols",
        make_unique: bool = True,
        cache: bool = False,
        cache_compression: Literal["gzip", "lzf"] | None = None,
        prefix: str | None = None,
    ) -> AnnData:
        """Read 10x-Genomics-formatted mtx directory.

        Parameters
        ----------
        path
            Path to directory for `.mtx` and `.tsv` files,
            e.g. './filtered_gene_bc_matrices/hg19/'.
        var_names
            The variables index.
        make_unique
            Whether to make the variables index unique by appending '-1',
            '-2' etc. or not.
        cache
            If `False`, read from source, if `True`, read from fast 'h5ad' cache.
        cache_compression
            See the h5py :ref:`dataset_compression`.
            (Default: `settings.cache_compression`)
        gex_only
            Only keep 'Gene Expression' data and ignore other feature types,
            e.g. 'Antibody Capture', 'CRISPR Guide Capture', or 'Custom'
        prefix
            Any prefix before `matrix.mtx`, `genes.tsv` and `barcodes.tsv`. For instance,
            if the files are named `patientA_matrix.mtx`, `patientA_genes.tsv` and
            `patientA_barcodes.tsv` the prefix is `patientA_`.
            (Default: no prefix)

        Returns
        -------
        An :class:`~anndata.AnnData` object

        """
        path = sc.readwrite.Path(path)
        prefix = "" if prefix is None else prefix
        is_legacy = (path / f"{prefix}genes.tsv").is_file()
        adata = sc.readwrite._read_10x_mtx(
            path,
            var_names=var_names,
            make_unique=make_unique,
            cache=cache,
            cache_compression=cache_compression,
            prefix=prefix,
            is_legacy=is_legacy,
        )
        gex_rows = adata.var["feature_types"] == feature_types
        return adata[:, gex_rows].copy()

def extract_crispr_features(matrix_path, barcodes_path, features_path, output_path):
    with tempfile.TemporaryDirectory() as tmpdirname:
        temp_local_path = f'{tmpdirname}'
        # Copy files to temp directory
        shutil.copy(matrix_path, f'{temp_local_path}/matrix.mtx.gz')
        shutil.copy(barcodes_path, f'{temp_local_path}/barcodes.tsv.gz')
        shutil.copy(features_path, f'{temp_local_path}/features.tsv.gz')

        print('Loading CRISPR Direct Capture data...')
        adata_crispr = read_10x_mtx_feature_types(temp_local_path, feature_types="CRISPR Direct Capture", var_names='gene_symbols', cache=True, prefix=None)
        
        print('Writing CRISPR adata to file...')
        adata_crispr.write(output_path, compression='gzip')
    
    return None

def main():
    """Main execution function."""
    args = parse_args()
    
    try:
        # Load 10x data
        extract_crispr_features(
            args.data_filtered_matrix,
            args.data_filtered_barcodes,
            args.data_filtered_features,
            args.output
        )
        print("\n✓ CRISPR feature extraction completed successfully")
        return 0
        
    except Exception as e:
        print(f"\n✗ Error during CRISPR feature extraction: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
