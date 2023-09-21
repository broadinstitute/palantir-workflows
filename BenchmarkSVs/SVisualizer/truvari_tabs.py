import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import quickboard.base as qbb
import quickboard.plugins as plg

from common_utils import read_and_postprocess, convert_missing_to_pass_filter, add_bbend_stats, sort_svtypes, sort_svlen_bins
from decorators import axes_mode
from plugins import make_type_selector, make_stat_selector, make_length_selector, make_interval_selector, make_axes_mode_selector
from user_config import COVARIATE_X, EXPERIMENT_ORDER, TRUVARI_DUP_TO_INS
from truvari_data import postprocess_truvari_bench, postprocess_truvari_closest
from upset_plot_utils import create_upset, make_disqualified_df


TRUVARI_BENCH_PATH = 'wdl_outputs/truvari_bench_summary.tsv'
TRUVARI_FN_CLOSEST_PATH = 'wdl_outputs/truvari_fn_closest.tsv'
TRUVARI_FP_CLOSEST_PATH = 'wdl_outputs/truvari_fp_closest.tsv'
TRUVARI_FN_INTERVALS_PATH = 'wdl_outputs/truvari_fn_intervals.tsv'
TRUVARI_FP_INTERVALS_PATH = 'wdl_outputs/truvari_fp_intervals.tsv'
TRUVARI_TP_BASE_INTERVALS_PATH = 'wdl_outputs/truvari_tp-base_intervals.tsv'
TRUVARI_TP_COMP_INTERVALS_PATH = 'wdl_outputs/truvari_tp-comp_intervals.tsv'

truvari_bench_df = read_and_postprocess(path=TRUVARI_BENCH_PATH, postprocessor=postprocess_truvari_bench)
truvari_fn_closest_df = read_and_postprocess(path=TRUVARI_FN_CLOSEST_PATH, postprocessor=postprocess_truvari_closest)
truvari_fp_closest_df = read_and_postprocess(path=TRUVARI_FP_CLOSEST_PATH, postprocessor=postprocess_truvari_closest)

if TRUVARI_DUP_TO_INS:
    truvari_fp_closest_df['LTYPE'] = truvari_fp_closest_df['LTYPE'].apply(lambda x: 'INS' if x == 'DUP' else x)
    truvari_fn_closest_df['RTYPE'] = truvari_fn_closest_df['RTYPE'].apply(lambda x: 'INS' if x == 'DUP' else x)

truvari_svlen_bins = sort_svlen_bins(list(truvari_bench_df['SVLEN_Bin'].unique()))
truvari_bench_df['SVLEN_Bin'] = pd.Categorical(truvari_bench_df['SVLEN_Bin'], ordered=True, categories=truvari_svlen_bins)

truvari_svtypes = sort_svtypes(list(truvari_bench_df['SVTYPE'].unique()))
truvari_bench_df['SVTYPE'] = pd.Categorical(truvari_bench_df['SVTYPE'], ordered=True, categories=truvari_svtypes)

truvari_bench_df = truvari_bench_df.sort_values(['SVTYPE', 'SVLEN_Bin'])

## Basic Truvari Bench Tab
@axes_mode
def make_truvari_bench_plot(df, axes_mode, stat, pct_overlap, svtype, svlen_bin, view_mode='Single'):
    if view_mode == 'Single':
        df = df[(df['SVTYPE'] == svtype) & (df['SVLEN_Bin'] == svlen_bin)]

    facet_row = None if view_mode == 'Single' else 'SVTYPE'
    facet_col = None if view_mode == 'Single' else 'SVLEN_Bin'

    if COVARIATE_X is not None and len(df[COVARIATE_X].unique()) > 1:
        x = COVARIATE_X
        y = stat
        title = f'{y} vs {x} Plot'
    else:
        x = 'Recall'
        y = 'Precision'
        title = 'Precision vs Recall Plot'
    
    type_ = df['SVTYPE'].iloc[0] if len(df) > 0 else "Empty df"
    title += f' for SVTYPE {type_}'
    
    interval = df['Interval'].iloc[0] if len(df) > 0 else "Empty df"
    title += f' over {interval}'

    df = df[df['Pct_Overlap'] == pct_overlap]
    title += f' w/ >= {pct_overlap}% overlap'
    
    hover_data = ['Base_Sample_Name', 'tp_base_count', 'fn_count', 'tp_comp_count', 'fp_count']
    category_orders = {'Experiment': EXPERIMENT_ORDER} if EXPERIMENT_ORDER is not None else None
    fig = px.scatter(df, x=x, y=y, color='Experiment', title=title, hover_name='Comp_Sample_Name', hover_data=hover_data, 
                      marginal_x='box', marginal_y='box', category_orders=category_orders, facet_row=facet_row, facet_col=facet_col)

    if view_mode == 'Grid':
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
        # fig.for_each_trace(lambda t: t.update(name=t.name.split("=")[1]))
        fig.update_xaxes(showticklabels=False, title_text='')
    return fig

def make_gt_concordance_plot(df, axes_mode, pct_overlap, svlen_bin, view_mode='Single'):
    if view_mode == 'Single':
        df = df[df['SVLEN_Bin'] == svlen_bin]

    facet_col = None if view_mode == 'Single' else 'SVLEN_Bin'

    category_orders = {'Experiment': EXPERIMENT_ORDER} if EXPERIMENT_ORDER is not None else None
    title = 'GT Concordance'
    if len(df) > 0:
        interval_name = df['Interval'].values[0]
        title += f' over {interval_name}'

        df = df[df['Pct_Overlap'] == pct_overlap]
        title += f' w/ >= {pct_overlap}% overlap'
    
    fig = px.box(df, x='SVTYPE', y='GT_concordance', color='Experiment', 
                 title=title, category_orders=category_orders, facet_col=facet_col)
    if axes_mode == 'Fixed':
        fig.update_layout(yaxis_range=[0, 1])

    if view_mode == 'Grid':
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
        # fig.for_each_trace(lambda t: t.update(name=t.name.split("=")[1]))
    return fig

truvari_bench_plugins = [
    plg.PlotInputRadioButtons(
        header='Variant Type',
        plot_input='svtype',
        data_values=list(truvari_bench_df['SVTYPE'].sort_values().unique())
    )    
]
#[make_type_selector(truvari_bench_df)]
if COVARIATE_X is not None and len(truvari_bench_df[COVARIATE_X].unique()) > 1:
    truvari_bench_plugins += [make_stat_selector(['Precision', 'Recall', 'F1_Score'])]

truvari_bench_plot = qbb.PlotPanel(
    header='Truvari Benchmarks',
    plotter=make_truvari_bench_plot,
    plot_inputs={},
    data_source=truvari_bench_df,
    plugins=truvari_bench_plugins
)

truvari_gt_concordance_plot = qbb.PlotPanel(
    header='GT Concordance',
    plotter=make_gt_concordance_plot,
    plot_inputs={},
    data_source=truvari_bench_df,
    plugins=[
        
    ]
)

truvari_bench_tab = qbb.BaseTab(
    tab_label='Truvari Bench',
    tab_header='',
    content_list=[
        truvari_bench_plot,
        truvari_gt_concordance_plot
    ],
    sidebar_plugins=[
        make_interval_selector(truvari_bench_df),
        plg.PlotInputSlider(
            header='Min Pct Overlap w/ Chosen Region',
            plot_input='pct_overlap',
            slider_min=truvari_bench_df['Pct_Overlap'].min(),
            slider_max=truvari_bench_df['Pct_Overlap'].max(),
            slider_default_value=truvari_bench_df['Pct_Overlap'].min(),
            # slider_step=10,
            slider_marks={
                str(i): {'label': str(i), 'style': {"transform": "rotate(-45deg)"}} for i in truvari_bench_df['Pct_Overlap'].unique()
            },
            updatemode='mouseup'
        ),
        plg.PlotInputRadioButtons(
            header='Select SVLEN Bin',
            plot_input='svlen_bin',
            data_values=truvari_svlen_bins
        ),
        plg.PlotInputRadioButtons(
            header='Select Viewing Mode',
            plot_input='view_mode',
            data_values=['Single', 'Grid']
        ),
        make_axes_mode_selector()
    ]
)

## Truvari Errors Tab
def make_closest_plot(df, sort_by, asc, mode, disq_values):
    title = f'Counts of Disqualified ({disq_values}) Sites (N = {len(df)})'
    asc = asc == 'Ascending'
    disq_df = make_disqualified_df(df, dist_threshold=500, size_ratio_threshold=0.7, color='Experiment')
    fig = create_upset(disq_df, title=title, sort_by=sort_by, asc=asc, mode=mode, color='Experiment', category_orders={'Experiment': EXPERIMENT_ORDER})
    return fig

def make_fp_closest_plot(df, sort_by, asc, mode):
    return make_closest_plot(df, sort_by, asc, mode, disq_values='FP')

def make_fn_closest_plot(df, sort_by, asc, mode):
    return make_closest_plot(df, sort_by, asc, mode, disq_values='FN')

fp_plot = qbb.PlotPanel(
    header='FP Stats',
    plotter=make_fp_closest_plot,
    plot_inputs={
        
    },
    data_source=truvari_fp_closest_df,
    plugins=[
        plg.PlotInputRadioButtons(
            header='Sort By',
            plot_input='sort_by',
            data_values=['Counts', 'Intersections']
        ),
        plg.PlotInputRadioButtons(
            header='Sort Order',
            plot_input='asc',
            data_values=['Descending', 'Ascending']
        ),
        plg.PlotInputRadioButtons(
            header='Bar Mode',
            plot_input='mode',
            data_values=['Percent', 'Counts']
        )
    ],
    plugin_wrap=3
)

fn_plot = qbb.PlotPanel(
    header='FN Stats',
    plotter=make_fn_closest_plot,
    plot_inputs={
        
    },
    data_source=truvari_fn_closest_df,
    plugins=[
        plg.PlotInputRadioButtons(
            header='Sort By',
            plot_input='sort_by',
            data_values=['Counts', 'Intersections']
        ),
        plg.PlotInputRadioButtons(
            header='Sort Order',
            plot_input='asc',
            data_values=['Descending', 'Ascending']
        ),
        plg.PlotInputRadioButtons(
            header='Bar Mode',
            plot_input='mode',
            data_values=['Percent', 'Counts']
        )
    ],
    plugin_wrap=3
)

truvari_errors_tab = qbb.BaseTab(
    tab_label='Truvari Errors',
    tab_header='Truvari Stats on Mismatched Variants',
    content_list=[
        fp_plot,
        fn_plot
    ],
    sidebar_plugins=[
        plg.DataFilterRadioButtons(
            header='kth Closest',
            data_col='kth_closest',
            data_values=[1, 2, 3]
        )
    ]
)