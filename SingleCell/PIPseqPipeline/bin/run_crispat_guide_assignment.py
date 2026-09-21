#!/usr/bin/env python3
"""
CRISPR guide assignment using CRISPAT.

This script performs guide assignment on single-cell data using CRISPAT's
Gaussian mixture model approach.
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
        description="Perform CRISPR guide assignment using CRISPAT"
    )
    parser.add_argument(
        "--adata",
        type=str,
        required=True,
        help="Path to the subsample's DRAGEN-filtered AnnData (h5ad) file"
    )
    parser.add_argument(
        "--subsample-id",
        type=str,
        required=True,
        help="Subsample identifier, used to name the output file and disambiguate cell barcodes across subsamples"
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        default=None,
        help="Number of processes to use for parallelization (default: all available cores)"
    )
    return parser.parse_args()

def run_crispat_guide_assignment(adata_path, subsample_id, num_processes):
    print(f"Extracting CRISPR Direct Capture features for {subsample_id}...")
    adata = sc.read_h5ad(adata_path)
    crispr_adata = adata[:, adata.var['feature_types'] == 'CRISPR Direct Capture']

    with tempfile.TemporaryDirectory() as tmpdir:
        crispr_adata_path = os.path.join(tmpdir, 'crispr_adata.h5ad')
        crispr_adata.write_h5ad(crispr_adata_path)

        print('Running CRISPAT Gaussian Mixture model...')
        crispat.ga_poisson_gauss(crispr_adata_path, f'crispat_ga/poisson_gauss/', parallelize=True, n_processes=num_processes, report_interval_seconds=30)

    # CRISPAT names its output file identically for every subsample -- rename it and suffix
    # 'cell' with the subsample ID so per-subsample outputs can be safely combined downstream
    # (barcodes can otherwise collide across subsamples), matching the suffixing convention
    # bin/concatenate_samples.py uses (ad.concat(..., index_unique='_')).
    assignments = pd.read_csv('crispat_ga/poisson_gauss/assignments.csv')
    assignments['cell'] = assignments['cell'].astype(str) + '_' + subsample_id
    assignments.to_csv(f'{subsample_id}.crispat_guide_assignments.csv', index=False)

def main():
    """Main execution function."""
    args = parse_args()

    try:
        run_crispat_guide_assignment(
            args.adata,
            args.subsample_id,
            args.num_processes
        )
        print("\n✓ CRISPAT guide assignment completed successfully")
        return 0

    except Exception as e:
        print(f"\n✗ Error during CRISPAT guide assignment: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
