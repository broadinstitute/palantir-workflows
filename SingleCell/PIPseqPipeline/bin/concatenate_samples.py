#!/usr/bin/env python3

import argparse
import os
import sys
import scanpy as sc
import anndata as ad
from pathlib import Path
import tempfile
import shutil

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
        with tempfile.TemporaryDirectory() as tmpdir:

            for i, (matrix, barcode, feature) in enumerate(zip(args.matrices, args.barcodes, args.features)):
                print(f"Reading subsample {subsample_ids[i]} from {matrix}...")
                
                os.mkdir(f"{tmpdir}/{i}")
                # Copy files to temp directory with expected names
                shutil.copy(matrix, f"{tmpdir}/{i}/matrix.mtx.gz")
                shutil.copy(barcode, f"{tmpdir}/{i}/barcodes.tsv.gz")
                shutil.copy(feature, f"{tmpdir}/{i}/features.tsv.gz")
                
                # Read the 10x format data
                adata = sc.read_10x_mtx(f"{tmpdir}/{i}", var_names='gene_symbols', gex_only=False, cache=False)
                adata.write_h5ad(f"{tmpdir}/{i}.h5ad", compression='gzip')
            
                print(f"subsample {subsample_ids[i]}: {adata.n_obs} cells, {adata.n_vars} features")        

            if len(args.matrices) == 1:
                print("Only one subsample provided, skipping concatenation.")
                shutil.move(f"{tmpdir}/0.h5ad", f"{args.supersample_basename}.h5ad")
            else:
                # Concatenate all subsamples
                print(f"\nConcatenating {len(args.matrices)} subsamples...")
                ad.experimental.concat_on_disk(
                    in_files=[f"{tmpdir}/{i}.h5ad" for i in range(len(args.matrices))],
                    out_file=f"{args.supersample_basename}.h5ad",
                    axis='obs',
                    label='subsample_id',
                    keys=subsample_ids,
                    index_unique='_',
                    join='outer',
                    merge='same')

        concatenated = sc.read_h5ad(f"{args.supersample_basename}.h5ad")
        print(f"Concatenated data: {concatenated.n_obs} cells, {concatenated.n_vars} features")
        
        print("Extracting CRISPR Direct Capture features...")
        print(concatenated.var_keys())
        crispr_mask = concatenated.var['feature_types'] == 'CRISPR Direct Capture'
        num_crispr_features = crispr_mask.sum()
        if num_crispr_features > 300:
            raise RuntimeError('We cannot process more than 300 guides right now due to runtime constraints.')
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
