#!/usr/bin/env python3
"""
Purity-based CRISPR guide assignment.

Assigns each cell to its most abundant guide when that guide's share of the
cell's total guide UMI counts (purity) and the cell's total guide UMI count
both clear fixed thresholds; otherwise the cell is left unassigned.
"""

import argparse
import sys
import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc

MIN_TOTAL_COUNT = 10
MIN_PURITY = 0.75


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Perform purity-based CRISPR guide assignment"
    )
    parser.add_argument(
        "--crispr-adata",
        type=str,
        required=True,
        help="Path to CRISPR-features-only AnnData (h5ad) file"
    )
    parser.add_argument(
        "--supersample-id",
        type=str,
        required=True,
        help="Supersample identifier, used to name the output file"
    )
    return parser.parse_args()


def assign_guides(crispr_adata_path):
    adata = sc.read_h5ad(crispr_adata_path)

    counts = adata.X
    if scipy.sparse.issparse(counts):
        counts = counts.toarray()
    counts = np.asarray(counts)

    n_cells, n_guides = counts.shape
    total_count = counts.sum(axis=1)

    # Descending sort per cell so index 0/1 are the top two guides by count.
    sorted_idx = np.argsort(-counts, axis=1)
    top1_idx = sorted_idx[:, 0]
    count_1st = counts[np.arange(n_cells), top1_idx]

    if n_guides >= 2:
        top2_idx = sorted_idx[:, 1]
        count_2nd = counts[np.arange(n_cells), top2_idx]
    else:
        count_2nd = np.zeros(n_cells)

    denom = count_1st + count_2nd
    purity = np.divide(
        count_1st, denom,
        out=np.zeros(n_cells, dtype=float),
        where=denom > 0
    )

    passes = (total_count > MIN_TOTAL_COUNT) & (purity > MIN_PURITY)

    guide_names = np.asarray(adata.var_names)
    top1_guide = guide_names[top1_idx]
    gRNA = np.where(passes, top1_guide, '')

    return pd.DataFrame({
        'cell': adata.obs_names,
        'gRNA': gRNA,
        'purity_1st_vs_2nd': purity,
        'total_count': total_count,
        'count_1st': count_1st,
        'count_2nd': count_2nd,
    })


def main():
    """Main execution function."""
    args = parse_args()

    try:
        print(f"Running purity-based guide assignment for {args.supersample_id}...")
        assignments = assign_guides(args.crispr_adata)

        output_path = f"{args.supersample_id}.purity_based_guide_assignments.csv"
        assignments.to_csv(output_path, index=False)

        print(f"\n✓ Purity-based guide assignment completed successfully ({output_path})")
        return 0

    except Exception as e:
        print(f"\n✗ Error during purity-based guide assignment: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
