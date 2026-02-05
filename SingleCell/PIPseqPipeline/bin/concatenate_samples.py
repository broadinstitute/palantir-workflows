#!/usr/bin/env python3

import argparse
import sys
import scanpy as sc
import anndata as ad
from pathlib import Path


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Concatenate multiple samples using AnnData"
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
        "--sample-ids",
        type=str,
        required=True,
        help="Comma-separated list of sample IDs"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output path for concatenated h5ad file"
    )
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_args()
    
    try:
        print(f"Concatenating {len(args.matrices)} samples...")
        
        # Parse sample IDs
        sample_ids = args.sample_ids.split(',')
        
        if len(args.matrices) != len(args.barcodes) or len(args.matrices) != len(args.features):
            raise ValueError("Number of matrix, barcode, and feature files must match")
        
        if len(args.matrices) != len(sample_ids):
            raise ValueError("Number of files must match number of sample IDs")
        
        # Read all samples
        adatas = []
        for i, (matrix, barcode, feature) in enumerate(zip(args.matrices, args.barcodes, args.features)):
            print(f"Reading sample {sample_ids[i]} from {matrix}...")
            
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
                adata = sc.read_10x_mtx(tmpdir, var_names='gene_symbols', cache=False)
            
            # Add sample ID to observations
            adata.obs['sample_id'] = sample_ids[i]
            
            # Filter for CRISPR Direct Capture features
            crispr_mask = adata.var['feature_types'] == 'CRISPR Direct Capture'
            if crispr_mask.sum() == 0:
                print(f"Warning: No CRISPR Direct Capture features found in sample {sample_ids[i]}")
            else:
                adata = adata[:, crispr_mask].copy()
                print(f"Sample {sample_ids[i]}: {adata.n_obs} cells, {adata.n_vars} CRISPR features")
                adatas.append(adata)
        
        if len(adatas) == 0:
            raise ValueError("No samples with CRISPR Direct Capture features found")
        
        # Concatenate all samples
        print(f"\nConcatenating {len(adatas)} samples...")
        concatenated = ad.concat(adatas, join='outer', merge='same')
        
        print(f"Concatenated dataset: {concatenated.n_obs} cells, {concatenated.n_vars} features")
        
        # Write output
        print(f"Writing concatenated data to {args.output}...")
        concatenated.write_h5ad(args.output)
        
        print("\n✓ Concatenation completed successfully")
        return 0
        
    except Exception as e:
        print(f"\n✗ Error during concatenation: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
