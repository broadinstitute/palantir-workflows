#!/usr/bin/env python3

import argparse
from collections import OrderedDict
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import subprocess
matplotlib.rcParams['font.family'] = 'Serif'

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate supersample-level QC report"
    )
    parser.add_argument(
        "--num-input-cells",
        type=int,
        required=True,
        help="Number of input cells"
    )
    parser.add_argument(
        "--subsample-qc-files",
        type=str,
        nargs='+',
        required=True,
        help="Paths to subsample QC files"
    )
    parser.add_argument(
        "--supersample-basename",
        type=str,
        required=True,
        help="Basename for output files"
    )
    parser.add_argument(
        "--supersample-id",
        type=str,
        required=True,
        help="Supersample identifier"
    )
    parser.add_argument(
        "--guide-assignments",
        type=str,
        required=False,
        help="Path to guide assignments file (optional)"
    )
    parser.add_argument(
        "--min-valid-guides",
        type=int,
        required=False,
        help="Minimum number of valid guides for guide assignment QC (must be set if --guide-assignments is set)"
    )
    parser.add_argument(
        "--max-valid-guides",
        type=int,
        required=False,
        help="Maximum number of valid guides for guide assignment QC (must be set if --guide-assignments is set)"
    )
    return parser.parse_args()

def generate_supersample_qc(guide_assignments, subsample_qc_files, supersample_basename, supersample_id, num_input_cells, min_valid_guides, max_valid_guides):
    subsample_metrics = pd.concat([pd.read_table(subsample_qc_file) for subsample_qc_file in subsample_qc_files], ignore_index=True)

    supersample_metrics = pd.DataFrame({'sample_id': [supersample_id]})
    supersample_metrics['N0 Input cells'] = num_input_cells
    supersample_metrics['N1 Passing cells'] = subsample_metrics['Passing cells'].sum()
    #supersample_metrics['N2 Guide containing passing cells'] = (subsample_metrics['Fraction passing cells with CRISPR reads'] * subsample_metrics['Passing cells']).sum()
    
    if guide_assignments is not None:
        num_guides_per_cell = guide_assignments.groupby('cell')['gRNA'].nunique().value_counts()
        n3_guide_assignment_passing_cells = num_guides_per_cell[(num_guides_per_cell.index >= min_valid_guides) & (num_guides_per_cell.index <= max_valid_guides)].sum()
        supersample_metrics['N4 Guide assignment passing cells'] = n3_guide_assignment_passing_cells

    return supersample_metrics

def guide_qc(guide_assignments, supersample_metrics, supersample_basename, supersample_id, min_valid_guides, max_valid_guides):
    guides_per_cell = guide_assignments.groupby('cell')['gRNA'].nunique()
    
    num_passing_cells = supersample_metrics.loc[supersample_metrics['sample_id'] == supersample_id, 'N1 Passing cells'].values[0]

    plotdata = guides_per_cell.value_counts()
    plotdata[0] = num_passing_cells - sum(plotdata)
    frac_too_few_guides = plotdata[plotdata.index < min_valid_guides].sum() / num_passing_cells
    frac_too_many_guides = plotdata[plotdata.index > max_valid_guides].sum() / num_passing_cells
    frac_valid_guides = plotdata[(plotdata.index >= min_valid_guides) & (plotdata.index <= max_valid_guides)].sum() / num_passing_cells

    plotdata = plotdata[plotdata.index <= 10]
    plotdata = plotdata.sort_index()

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(plotdata.index, plotdata.values, color='C0')
    bars[0].set_color('gray')
    for bar in bars[3:]:
        bar.set_color('C1')
    ax.set_xlabel('Number of guides assigned to cell')
    ax.set_ylabel('Number of cells')
    ax.set_title(f'{supersample_id}: Distribution of number of guides assigned to cells')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}'))
    ax.set_xticks(plotdata.index)

    # Legend handles
    from matplotlib.patches import Patch
    label_too_few = '0' if min_valid_guides == 1 else f'0-{min_valid_guides-1}'
    label_valid_guides = f'{min_valid_guides}-{max_valid_guides}' if max_valid_guides > min_valid_guides else f'{min_valid_guides}'
    label_too_many = f'>{max_valid_guides+1}+'
    legend_handles = [
        Patch(color='gray', label=f'{label_too_few} guides: {frac_too_few_guides:.1%}'),
        Patch(color='C0', label=f'{label_valid_guides} guides: {frac_valid_guides:.1%}'),
        Patch(color='C1', label=f'{label_too_many} guides: {frac_too_many_guides:.1%}'),
    ]
    ax.legend(handles=legend_handles, loc='upper right')

    plt.tight_layout()
    
    fig.savefig(f'{supersample_basename}.guide_assignment_distribution.png', dpi=300)

def generate_sankey_plot(supersample_metrics, supersample_id, supersample_basename):
    this_supersample_metrics = supersample_metrics[supersample_metrics['sample_id'] == supersample_id].iloc[0]
    sankey_source = f'''
// === Nodes and Flows ===

N0: Input cells [{this_supersample_metrics['N1 Passing cells']:.0f}] N1: GEX passing cells
N0: Input cells [{this_supersample_metrics['N0 Input cells'] - this_supersample_metrics['N1 Passing cells']:.0f}] Empty droplets
N1: GEX passing cells [{this_supersample_metrics['N2 Guide containing passing cells']:.0f}] N2: Guide containing passing cells
N1: GEX passing cells [{this_supersample_metrics['N1 Passing cells'] - this_supersample_metrics['N2 Guide containing passing cells']:.0f}] Cells with no guide reads
N2: Guide containing passing cells [{this_supersample_metrics['N4 Guide assignment passing cells']:.0f}] N4: Guide assignment passing cells
N2: Guide containing passing cells [{this_supersample_metrics['N2 Guide containing passing cells'] - this_supersample_metrics['N4 Guide assignment passing cells']:.0f}] Cells with no valid guide assignment

:N0: Input cells #777777
:Empty droplets #aaaaaa
:Cells with no guide reads #aaaaaa
:Cells with no valid guide assignment #aaaaaa

// === Settings ===

size w 1200
h 600
margin l 12
r 12
t 18
b 20
bg color #ffffff
transparent N
node w 12
h 50
spacing 34.5
border 0
theme none
color #006db6
opacity 1
flow curvature 0.5
inheritfrom target
color #999999
opacity 0.45
layout order exact
justifyorigins N
justifyends N
reversegraph N
attachincompletesto nearest
labels color #000000
hide N
highlight 0.35
fontface sans-serif
linespacing 0.2
relativesize 110
magnify 100
labelname appears Y
size 16
weight 400
labelvalue appears Y
fullprecision Y
position below
weight 400
labelposition autoalign -1
scheme auto
first before
breakpoint 6
value format ',.'
prefix ''
suffix ''
themeoffset a 6
b 0
c 0
d 0
meta mentionsankeymatic N
listimbalances Y'''
    with open(f'{supersample_basename}.sankey_input.txt', 'w') as f:
        f.write(sankey_source)
    
    sankey_process = subprocess.run(
        ['node', '/app/sankeymatic_local/cli.js',
         f'{supersample_basename}.sankey_input.txt', f'{supersample_basename}.sankey.png'],
         capture_output=True)
    if sankey_process.returncode != 0:
        raise RuntimeError(f'''Sankeymatic process failed.
                           stdout:
                           {sankey_process.stdout.decode()}
                           stderr:
                            {sankey_process.stderr.decode()}''')

def main():
    args = parse_args()
    try:
        guide_assignments = None
        if args.guide_assignments is not None:
            if (args.min_valid_guides is None or args.max_valid_guides is None):
                raise ValueError("Both --min-valid-guides and --max-valid-guides must be set if --guide-assignments is provided")
            guide_assignments = pd.read_csv(args.guide_assignments) if args.guide_assignments is not None and args.guide_assignments != 'NO_FILE' else None
        #else:
        #    raise RuntimeError("Right now, we require guide assignments to generate the supersample QC report. Change this if you want to allow generating the report without guide assignments.")

        print(f"Generating supersample QC for {args.supersample_id}...")
        supersample_metrics = generate_supersample_qc(
            guide_assignments,
            args.subsample_qc_files,
            args.supersample_basename,
            args.supersample_id,
            args.num_input_cells,
            min_valid_guides=1,
            max_valid_guides=2
        )
        supersample_metrics.to_csv(f'{args.supersample_basename}.supersample_qc_metrics.tsv', index=False, sep='\t')


        if guide_assignments is not None:
            guide_qc(guide_assignments, supersample_metrics, args.supersample_basename, args.supersample_id, min_valid_guides=args.min_valid_guides, max_valid_guides=args.max_valid_guides)
            #generate_sankey_plot(supersample_metrics, args.supersample_id, args.supersample_basename)

        print("\n✓ Supersample QC generation completed successfully")
        return 0
    except Exception as e:
        print(f"\n✗ Error during supersample QC generation: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    main()
