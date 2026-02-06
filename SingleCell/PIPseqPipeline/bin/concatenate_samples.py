#!/usr/bin/env python3

import argparse
import sys
import scanpy as sc
import anndata as ad
from pathlib import Path

sc.settings.n_jobs = -1

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Concatenate multiple subsamples using AnnData"
    )
    parser.add_argument(
        "--matrices",
        type=str,
        nargs='+',
        required=True,
        help="Paths to filtered matrix files"
    )
    parser.add_argument(
        "--barcodes",
        type=str,
        nargs='+',
        required=True,
        help="Paths to filtered barcodes files"
    )
    parser.add_argument(
        "--features",
        type=str,
        nargs='+',
        required=True,
        help="Paths to filtered features files"
    )
    parser.add_argument(
        "--subsample-ids",
        type=str,
        required=True,
        help="Comma-separated list of subsample IDs"
    )
    parser.add_argument(
        "--supersample-basename",
        type=str,
        required=True,
        help="Basename for output files"
    )
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_args()
    try:
        print(f"Concatenating {len(args.matrices)} subsamples...")
        
        # Parse subsample IDs
        subsample_ids = args.subsample_ids.split(',')
        
        if len(args.matrices) != len(args.barcodes) or len(args.matrices) != len(args.features):
            raise ValueError("Number of matrix, barcode, and feature files must match")
        
        if len(args.matrices) != len(subsample_ids):
            raise ValueError("Number of files must match number of subsample IDs")
        
        # Read all subsamples
        adatas = []
        for i, (matrix, barcode, feature) in enumerate(zip(args.matrices, args.barcodes, args.features)):
            print(f"Reading subsample {subsample_ids[i]} from {matrix}...")
            
            # Create a temporary directory structure for scanpy
            # scanpy expects the files to be in a directory with specific naming
            import tempfile
            import shutil
            import gzip
            
            with tempfile.TemporaryDirectory() as tmpdir:
                # Copy files to temp directory with expected names
                shutil.copy(matrix, f"{tmpdir}/matrix.mtx.gz")
                shutil.copy(barcode, f"{tmpdir}/barcodes.tsv.gz")
                shutil.copy(feature, f"{tmpdir}/features.tsv.gz")
                
                # Read the 10x format data
                adata = sc.read_10x_mtx(tmpdir, var_names='gene_symbols', gex_only=False, cache=False)
                        
            # Filter for CRISPR Direct Capture features
            
            print(f"subsample {subsample_ids[i]}: {adata.n_obs} cells, {adata.n_vars} features")
            adatas.append(adata)
        
        if len(adatas) == 0:
            raise ValueError("No subsamples found")
        
        if len(adatas) == 1:
            print("Only one subsample provided, skipping concatenation.")
            concatenated = adatas[0]
        else:
            # Concatenate all subsamples
            print(f"\nConcatenating {len(adatas)} subsamples...")
            concatenated = ad.concat(adatas, join='outer', merge='same', label='subsample_id', keys=subsample_ids, index_unique='_')

            print(f"Concatenated dataset: {concatenated.n_obs} cells, {concatenated.n_vars} features")

        # Write output
        output_path = f"{args.supersample_basename}.h5ad"
        print(f"Writing data to {output_path}...")
        concatenated.write_h5ad(output_path, compression='gzip')
        
        print("Extracting CRISPR Direct Capture features...")
        print(concatenated.var_keys())
        crispr_mask = concatenated.var['feature_types'] == 'CRISPR Direct Capture'
        concatenated_crispr = concatenated[:, crispr_mask]

        crispr_output_path = f"{args.supersample_basename}.crispr.h5ad"
        print(f"Writing CRISPR data to {crispr_output_path}...")
        concatenated_crispr.write_h5ad(crispr_output_path, compression='gzip')
        
        print("\n✓ Concatenation completed successfully")
        return 0
        
    except Exception as e:
        print(f"\n✗ Error during concatenation: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
