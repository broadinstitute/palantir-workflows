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
        "--crispr-adata",
        type=str,
        required=True,
        help="Path to feature barcode reference file"
    )
    parser.add_argument(
        "--output-basename",
        type=str,
        required=True,
        help="Base name for output assignments CSV file"
    )
    return parser.parse_args()

def run_guide_assignment(crispr_adata_path, output_basename):
    
    print('Running CRISPAT Gaussian Mixture model...')
    crispat.ga_poisson_gauss(f'{crispr_adata_path}', f'crispat_ga/poisson_gauss/')

    shutil.copy(f'crispat_ga/poisson_gauss/assignments.csv', f'{output_basename}.assignments.csv')
    
    return None

def main():
    """Main execution function."""
    args = parse_args()
    
    try:
        # Load 10x data
        run_guide_assignment(
            args.crispr_adata,
            args.output_basename
        )
        print("\n✓ Guide assignment completed successfully")
        return 0
        
    except Exception as e:
        print(f"\n✗ Error during guide assignment: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
