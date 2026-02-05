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
        "--output-basename",
        type=str,
        required=True,
        help="Base name for output files"
    )
    parser.add_argument(
        "--guide-assignments",
        type=str,
        required=False,
        default=None,
        help="Optional: Path to guide assignments CSV file from CRISPAT"
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

def write_rankplot(barcode_summary_path, sample_id, output_basename):
    print('    Generating barcode rank plot...')
    barcode_metrics = get_barcode_metrics(barcode_summary_path, sample_id)
    barcode_metrics.to_csv(f'{output_basename}.qc_barcode_metrics.tsv', index=False, sep='\t')

def get_sc_metrics(sc_metrics_path, sample_id):
    columns = ['metrics_name', 'value_numeric', 'value_frac']
    metrics = pd.read_csv(sc_metrics_path, names=columns)
    metrics = metrics.rename(columns={
        'metrics_name': 'metric',
        'value_numeric': 'value',})
    
    calculated_rows = list()
    calculated_rows.append({'metric': 'Total Gene Input Reads', 'value': '{:.4f}'.format(metrics.loc[metrics['metric'] == 'Total barcoded reads', 'value'].values[0] + metrics.loc[metrics['metric'] == 'Reads with non-matching barcodes', 'value'].values[0] + metrics.loc[metrics['metric'] == 'Reads missing barcodes', 'value'].values[0])})
    calculated_rows.append({'metric': 'Fraction Reads with Valid Barcodes', 'value_frac': '{:.4f}'.format(metrics.loc[metrics['metric'] == 'Total barcoded reads', 'value'].values[0] / (metrics.loc[metrics['metric'] == 'Total barcoded reads', 'value'].values[0] + metrics.loc[metrics['metric'] == 'Reads with non-matching barcodes', 'value'].values[0] + metrics.loc[metrics['metric'] == 'Reads missing barcodes', 'value'].values[0]))})
    calculated_rows.append({'metric': 'Fraction Reads with valid IMI', 'value_frac': '{:.4f}'.format((metrics.loc[metrics['metric'] == 'Reads with valid molecular identifier sequences', 'value'].values[0] + metrics.loc[metrics['metric'] == 'Reads with corrected molecular identifier sequences', 'value'].values[0]) / metrics.loc[metrics['metric'] == 'Total barcoded reads', 'value'].values[0])})
    calculated_rows.append({'metric': 'Fraction Mapped Reads', 'value_frac': '{:.4f}'.format((metrics.loc[metrics['metric'] == 'Unique exon matching reads', 'value'].values[0] + metrics.loc[metrics['metric'] == 'Unique intron matching reads', 'value'].values[0] + metrics.loc[metrics['metric'] == 'Mitochondrial reads', 'value'].values[0]) / metrics.loc[metrics['metric'] == 'Total barcoded reads', 'value'].values[0])})
    calculated_rows.append({'metric': 'Mean reads per cell', 'value': '{:.4f}'.format(metrics.loc[metrics['metric'] == 'Total gene reads', 'value'].values[0] / metrics.loc[metrics['metric'] == 'Passing cells', 'value'].values[0])})
    calculated_rows.append({'metric': 'Fraction CRISPR Barcoded Reads', 'value_frac': '{:.4f}'.format(metrics.loc[metrics['metric'] == 'Total CRISPR reads matching known barcodes', 'value'].values[0] / metrics.loc[metrics['metric'] == 'Total feature reads', 'value'].values[0])})

    metrics = pd.concat([metrics, pd.DataFrame(calculated_rows)], ignore_index=True)
    metrics['sample_id'] = sample_id
    return metrics

def write_metrics(sc_metrics_path, sample_id, output_basename):
    print('    Writing single-cell RNA QC metrics...')
    metrics = get_sc_metrics(sc_metrics_path, sample_id)
    metrics.to_csv(f'{output_basename}.qc_metrics.tsv', index=False, sep='\t')

from typing import OrderedDict


def guide_qc(guide_assignments_path, metrics_path, sample_id, output_basename):
    guide_assignments = pd.read_csv(guide_assignments_path)
    metrics = pd.read_csv(metrics_path, names=['metric_type', 'sample_id', 'metric', 'value', 'frac'])
    
    guides_per_cell = guide_assignments.groupby('cell')['gRNA'].nunique()
    
    num_passing_cells = metrics.query('metric == "Passing cells"').value.values[0]
    plotdata = guides_per_cell.value_counts()
    plotdata[0] = num_passing_cells - sum(plotdata)
    frac_no_guides = plotdata[0] / num_passing_cells
    frac_1_2_guides = (plotdata[1]+plotdata[2]) / num_passing_cells
    frac_3_plus_guides = 1 - (plotdata[0] + plotdata[1] + plotdata[2]) / num_passing_cells

    plotdata = plotdata[plotdata.index <= 10]
    plotdata = plotdata.sort_index()

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(plotdata.index, plotdata.values, color='C0')
    bars[0].set_color('gray')
    for bar in bars[3:]:
        bar.set_color('C1')
    ax.set_xlabel('Number of guides assigned to cell')
    ax.set_ylabel('Number of cells')
    ax.set_title(f'{sample_id}: Distribution of number of guides assigned to cells')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}'))
    ax.set_xticks(plotdata.index)

    # Legend handles
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(color='gray', label=f'0 guides: {frac_no_guides:.1%}'),
        Patch(color='C0', label=f'1-2 guides: {frac_1_2_guides:.1%}'),
        Patch(color='C1', label=f'3+ guides: {frac_3_plus_guides:.1%}'),
    ]
    ax.legend(handles=legend_handles, loc='upper right')

    plt.tight_layout()
    
    fig.savefig(f'{output_basename}.qc_guide_assignment_distribution.png', dpi=300)

    guide_assignment_stats = OrderedDict()
    guide_assignment_stats['sample_id'] = sample_id
    guide_assignment_stats['num_cells_with_zero_guides'] = plotdata[0]
    guide_assignment_stats['num_cells_with_one_or_two_guides'] = plotdata[1] + plotdata[2]
    guide_assignment_stats['num_cells_with_three_or_more_guides'] = num_passing_cells - (plotdata[0] + plotdata[1] + plotdata[2])
    guide_assignment_stats['frac_cells_with_zero_guides'] = frac_no_guides
    guide_assignment_stats['frac_cells_with_one_or_two_guides'] = frac_1_2_guides
    guide_assignment_stats['frac_cells_with_three_or_more_guides'] = frac_3_plus_guides

    with open(f'{output_basename}.qc_guide_assignment_stats.tsv', 'w') as f:
        f.write('\t'.join(guide_assignment_stats.keys()) + '\n')
        f.write('\t'.join(str(v) for v in guide_assignment_stats.values()) + '\n')

def main():
    """Main execution function."""
    args = parse_args()

    asdfasdfjo;i
    
    try:
        print('Generating report data...')
        # Load inputs and generate plots/metrics
        write_rankplot(
            args.barcode_summary,
            args.sample_id,
            args.output_basename
        )
        write_metrics(
            args.scrna_metrics,
            args.sample_id,
            args.output_basename
        )
        if args.guide_assignments is not None:
            guide_qc(
                args.guide_assignments,
                args.scrna_metrics,
                args.sample_id,
                args.output_basename
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