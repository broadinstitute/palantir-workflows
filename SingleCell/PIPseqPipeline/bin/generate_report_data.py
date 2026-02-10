#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
matplotlib.rcParams['font.family'] = 'Serif'


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract CRISPR features from 10x-style filtered matrix"
    )
    parser.add_argument(
        "--num-input-cells",
        type=int,
        required=True,
        help="Number of input cells"
    )
    parser.add_argument(
        "--scrna-metrics",
        type=str,
        required=True,
        help="Path to scRNA metrics CSV file"
    )
    parser.add_argument(
        "--barcode-summary",
        type=str,
        required=True,
        help="Path to barcode summary TSV file"
    )
    parser.add_argument(
        "--sample-id",
        type=str,
        required=True,
        help="Sample identifier"
    )
    parser.add_argument(
        "--sample-basename",
        type=str,
        required=True,
        help="Base name for output files"
    )
    return parser.parse_args()

def get_barcode_metrics(barcode_summary_path, sample_id):
    data = pd.read_table(barcode_summary_path, usecols=['Molecules', 'Filter'])
    data = data.sort_values(by='Molecules', ascending=False).reset_index(drop=True)
    data = pd.concat([
        data.drop_duplicates(keep='first'),
        data.drop_duplicates(keep='last')
    ]).sort_values(by=['Molecules'], ascending=[False])
    data = data.loc[~data.index.duplicated(keep='first')]
    data.index = data.index + 1
    data = data.reset_index(names='Rank')
    data['sample_id'] = sample_id
    print('Done')
    return data

def write_rankplot(barcode_summary_path, sample_id, sample_basename):
    print('    Generating barcode rank plot...')
    barcode_metrics = get_barcode_metrics(barcode_summary_path, sample_id)
    barcode_metrics.to_csv(f'{sample_basename}.qc_barcode_metrics.tsv', index=False, sep='\t')

def get_sc_metrics(sc_metrics_path, sample_id):
    metrics = pd.read_csv(sc_metrics_path, names=['metric_type', 'sample_id', 'metric', 'value', 'frac'], dtype={'sample_id': str})
    
    metrics = metrics.pivot(index=['sample_id'], columns='metric', values='value')
    if metrics['sample_id'].values[0] != sample_id:
        raise ValueError(f"Sample ID in sc metrics file ({metrics['sample_id'].values[0]}) does not match expected sample ID ({sample_id})")
    return metrics

def write_metrics(sc_metrics_path, sample_id, sample_basename):
    print('    Writing single-cell RNA QC metrics...')
    metrics = get_sc_metrics(sc_metrics_path, sample_id)
    metrics.to_csv(f'{sample_basename}.qc_metrics.tsv', index=False, sep='\t')

def main():
    """Main execution function."""
    args = parse_args()
    
    try:
        print('Generating report data...')
        # Load inputs and generate plots/metrics
        write_rankplot(
            args.barcode_summary,
            args.sample_id,
            args.sample_basename
        )
        write_metrics(
            args.scrna_metrics,
            args.sample_id,
            args.sample_basename
        )
        print("\n✓ Report data generation completed successfully")
        return 0
        
    except Exception as e:
        print(f"\n✗ Error during report data generation: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())