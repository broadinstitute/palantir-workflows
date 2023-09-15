#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go

import quickboard.base as qbb
import quickboard.plugins as plg
from quickboard.app import start_app, deploy_app, app

from upset_plot_utils import create_upset, make_disqualified_df


# In[ ]:





# # User Config
# 
# Use this section to set variables that change the output and settings of the dashboard. Anything after this header's section should be untouched, and modifying the following code is advanced usage.

# In[ ]:


# (Optional) Covariate column names for benchmarking plots
# Should be numeric to plot along x-axis, e.g. average coverage
COVARIATE_X = 'Experiment'

# Option to skip second HWE plot for variant in Comp VCFs
# Useful to skip for single sample vs panel (Base) comparisons
SKIP_COMP_HWE = True

# Coerice missing FILTER to PASS
# Ensure '.' is interpreted as PASS FILTER in case base/comp have different conventions
# (Though this should be done before running through Wittyer...)
MAKE_MISSING_PASS_FILTER = True

# Use to fix explicit order for experiment groups
EXPERIMENT_ORDER = None

# Use if running Truvari w/ dup-to-ins on
TRUVARI_DUP_TO_INS = True

# Toggle which categories of data to plot in Dashboard
INCLUDE_QC = True
INCLUDE_WITTYER = True
INCLUDE_TRUVARI = True


# In[ ]:





# # Import Data

# ## File Paths

# In[ ]:


QC_DATA_PATH = 'wdl_outputs/combined_qc_stats.tsv'

TRUVARI_BENCH_PATH = 'wdl_outputs/truvari_bench_summary.tsv'
TRUVARI_FN_CLOSEST_PATH = 'wdl_outputs/truvari_fn_closest.tsv'
TRUVARI_FP_CLOSEST_PATH = 'wdl_outputs/truvari_fp_closest.tsv'
TRUVARI_FN_INTERVALS_PATH = 'wdl_outputs/truvari_fn_intervals.tsv'
TRUVARI_FP_INTERVALS_PATH = 'wdl_outputs/truvari_fp_intervals.tsv'
TRUVARI_TP_BASE_INTERVALS_PATH = 'wdl_outputs/truvari_tp-base_intervals.tsv'
TRUVARI_TP_COMP_INTERVALS_PATH = 'wdl_outputs/truvari_tp-comp_intervals.tsv'

WITTYER_STATS_PATH = 'wdl_outputs/wittyer_stats-cleaned.tsv'
WITTYER_TRUTH_PATH = 'wdl_outputs/wittyer_truth_intervals.tsv'
WITTYER_QUERY_PATH = 'wdl_outputs/wittyer_query_intervals.tsv'
WITTYER_NOGT_PATH = 'wdl_outputs/wittyer_nogt_intervals.tsv'

# Useful to have full list for ad hoc scripting
# PATHS = [
#     QC_DATA_PATH, TRUVARI_BENCH_PATH, TRUVARI_FN_CLOSEST_PATH, TRUVARI_FP_CLOSEST_PATH, TRUVARI_FN_INTERVALS_PATH, TRUVARI_FP_INTERVALS_PATH,
#     TRUVARI_TP_BASE_INTERVALS_PATH, TRUVARI_TP_COMP_INTERVALS_PATH, WITTYER_STATS_PATH, WITTYER_TRUTH_PATH, WITTYER_QUERY_PATH, WITTYER_NOGT_PATH
# ]


# In[ ]:


def read_and_postprocess(path, toggle, postprocessor):
    # Check file path exists (was computed in WDL) and if not toggled off read and postprocess
    if toggle and os.path.isfile(path):
        df = pd.read_csv(path, sep='\t')
        df = postprocessor(df)
        return df
    else:
        return pd.DataFrame()


# In[ ]:





# ## General Postprocessing Functions

# In[ ]:


# Add BBEND stats
def add_bbend_stats(df):
    interval_names = sorted(list(set([c.split('-')[0] for c in df.columns if '-count' in c])))
    for i in interval_names:
        df[f'{i}-BBEND-overlap'] = (df[f'{i}-LBEND-overlap'] + df[f'{i}-RBEND-overlap'])/2
        df[f'{i}-BBEND-count'] = df[f'{i}-LBEND-count'] + df[f'{i}-RBEND-count']


# In[ ]:


# Convert '.' FILTER to 'PASS'
def convert_missing_to_pass_filter(df, MAKE_MISSING_PASS_FILTER):
    if 'FILTER' in df.columns:
        df['FILTER'] = df['FILTER'].replace('.', 'PASS')


# In[ ]:





# ## QC Data

# ### Postprocessing

# In[ ]:


# Bin lengths
qc_length_bins = [
    '<100bp',
    '100 - 500bp',
    '500bp - 2.5kb',
    '2.5 - 10kb',
    '10 - 50kb',
    '> 50kb'
]

def bin_length(length):
    if length < 100:
        return '<100bp'
    elif length <= 500:
        return '100 - 500bp'
    elif length < 2500:
        return '500bp - 2.5kb'
    elif length < 10_000:
        return '2.5 - 10kb'
    elif length < 50_000:
        return '10 - 50kb'
    else:
        return '> 50kb'

qc_af_bins = [
    'AC=1',
    '< 1%',
    '1-10%',
    '10-50%',
    '> 50%'
]

def bin_af(af):
    if af == 0:
        return '0%'
    elif af < .01:
        return '< 1%'
    elif af < .1:
        return '1-10%'
    elif af < .5:
        return '10-50%'
    else:
        return '> 50%'

def postprocess_qc_df(qc_df):
    # Separate base vs comp experiment groups
    qc_df['Experiment_Suffix'] = qc_df['Experiment'].apply(lambda x: x.split('-')[-1])
    qc_num_exp_groups = len(qc_df['Experiment'].unique())

    # Add AC_Ref counts
    qc_df['AC_Ref'] = qc_df['NS'] - qc_df['AC_Het'] - qc_df['AC_Hom'] / 2

    # Bin by length
    qc_df['SVLEN'] = qc_df['SVLEN'].replace('.', 0).astype(int)
    qc_df['SVLEN_Bin'] = qc_df['SVLEN'].apply(bin_length)

    # Bin AFs
    qc_df['AF'] = qc_df['AF'].replace('.', 0).astype(float)
    qc_df['AF_Bin'] = qc_df['AF'].replace('.', 0).astype(float)

    # Idea: First write AF_Bin column using AF, then overwrite using AC column if AC == 1
    qc_df['AF_Bin'] = qc_df['AF'].apply(bin_af)
    qc_df['AF_Bin'] = qc_df.apply(lambda x: 'AC=1' if x['AC'] == 1 else x['AF_Bin'], axis=1)

    # Make categorical columns in proper order
    qc_df['SVLEN_Bin'] = pd.Categorical(qc_df['SVLEN_Bin'], ordered=True, categories=qc_length_bins)
    qc_df['AF_Bin'] = pd.Categorical(qc_df['AF_Bin'], ordered=True, categories=qc_af_bins)

    add_bbend_stats(qc_df)
    convert_missing_to_pass_filter(qc_df, MAKE_MISSING_PASS_FILTER)

    return qc_df


# ### Read Data

# In[ ]:


qc_df = read_and_postprocess(path=QC_DATA_PATH, toggle=INCLUDE_QC, postprocessor=postprocess_qc_df)


# In[ ]:





# ## Truvari Data

# ### Postprocess

# In[ ]:





# In[ ]:


def postprocess_truvari_bench(df):
    df = df.rename(columns={'SV_Type': 'SVTYPE'})
    df['Experiment'] = df['Experiment'].astype(str)
    convert_missing_to_pass_filter(df, MAKE_MISSING_PASS_FILTER)
    return df

def postprocess_truvari_closest(df):
    df['LLEN'] = df['LLEN'].replace('.', 0).astype(int)
    df['RLEN'] = df['RLEN'].replace('.', 0).astype(int)
    df['SIZE_RATIO'] = df[['LLEN', 'RLEN']].min(axis=1) / df[['LLEN', 'RLEN']].max(axis=1)
    
    convert_missing_to_pass_filter(df, MAKE_MISSING_PASS_FILTER)

    k = 3  # Toggle this to compare other k-closest; should match WDL input
    df['kth_closest'] = df.index % k + 1
    return df


# ### Read Data

# In[ ]:


truvari_bench_df = read_and_postprocess(path=TRUVARI_BENCH_PATH, toggle=INCLUDE_TRUVARI, postprocessor=postprocess_truvari_bench)

truvari_fn_closest_df = read_and_postprocess(path=TRUVARI_FN_CLOSEST_PATH, toggle=INCLUDE_TRUVARI, postprocessor=postprocess_truvari_closest)
truvari_fp_closest_df = read_and_postprocess(path=TRUVARI_FP_CLOSEST_PATH, toggle=INCLUDE_TRUVARI, postprocessor=postprocess_truvari_closest)

truvari_fn_intervals_df = read_and_postprocess(path=TRUVARI_FN_INTERVALS_PATH, toggle=INCLUDE_TRUVARI, postprocessor=lambda df: df)
truvari_fp_intervals_df = read_and_postprocess(path=TRUVARI_FP_INTERVALS_PATH, toggle=INCLUDE_TRUVARI, postprocessor=lambda df: df)
truvari_tpbase_intervals_df = read_and_postprocess(path=TRUVARI_TP_BASE_INTERVALS_PATH, toggle=INCLUDE_TRUVARI, postprocessor=lambda df: df)
truvari_tpcomp_intervals_df = read_and_postprocess(path=TRUVARI_TP_COMP_INTERVALS_PATH, toggle=INCLUDE_TRUVARI, postprocessor=lambda df: df)


# In[ ]:


if TRUVARI_DUP_TO_INS:
    truvari_fp_closest_df['LTYPE'] = truvari_fp_closest_df['LTYPE'].apply(lambda x: 'INS' if x == 'DUP' else x)
    truvari_fn_closest_df['RTYPE'] = truvari_fn_closest_df['RTYPE'].apply(lambda x: 'INS' if x == 'DUP' else x)


# In[ ]:





# In[ ]:





# ## Wittyer Data

# ### Postprocessing

# In[ ]:


def postprocess_adv_wittyer(df):
    add_bbend_stats(df)
    return df


# In[ ]:





# ### Read Data

# In[ ]:


wittyer_stats_df = read_and_postprocess(path=WITTYER_STATS_PATH, toggle=INCLUDE_WITTYER, postprocessor=lambda df: df)
wittyer_truth_df = read_and_postprocess(path=WITTYER_TRUTH_PATH, toggle=INCLUDE_WITTYER, postprocessor=lambda df: df)
wittyer_query_df = read_and_postprocess(path=WITTYER_QUERY_PATH, toggle=INCLUDE_WITTYER, postprocessor=lambda df: df)
wittyer_nogt_df = read_and_postprocess(path=WITTYER_NOGT_PATH, toggle=INCLUDE_WITTYER, postprocessor=lambda df: df)

adv_wittyer_df = pd.concat([wittyer_query_df, wittyer_truth_df])
if INCLUDE_WITTYER:
    adv_wittyer_df = postprocess_adv_wittyer(adv_wittyer_df)
    wittyer_bins = ['All'] + sorted([str(x) for x in wittyer_stats_df['Bin'].unique() if (x != 'All') & (str(x) != 'nan')]) + ['nan']


# In[ ]:





# In[ ]:





# # Quickboard

# ## Plugin Methods

# In[ ]:


# Interval Plugins
def make_interval_selector(df):
    return plg.DataFilterRadioButtons(
        header="Interval List",
        data_col='Interval',
        data_values=list(df['Interval'].unique())
    )

def make_interval_name_selector(df):
    names = sorted(list(set([c.split('-')[0] for c in df.columns if '-count' in c])))
    return plg.PlotInputRadioButtons(
        header="Interval List",
        plot_input='interval_name',
        data_values=names
    )

def make_breakpoint_selector():
    return plg.PlotInputRadioButtons(
        header="Breakpoint Restriction",
        plot_input='breakpoint',
        data_values=['Full', 'Left', 'Right', 'Both']
    )

def make_pct_overlap_slider():
    return plg.PlotInputRangeSlider(
        header='Pct Overlap w/ Chosen Region',
        plot_input='pct_overlap',
        slider_min=0,
        slider_max=100,
        slider_default_values=[0, 100],
        slider_step=10,
        slider_marks={
            str(i): {'label': str(i), 'style': {"transform": "rotate(-45deg)"}} for i in range(0, 101, 10)
        },
        updatemode='mouseup'
    )

def make_interval_plugin_bundle(df):
    return [
        make_interval_name_selector(df),
        make_breakpoint_selector(),
        make_pct_overlap_slider()
    ]

# SV Plugins
def make_type_selector(df):
    return plg.DataFilterRadioButtons(
        header="Variant Type",
        data_col='SVTYPE',
        data_values=list(df['SVTYPE'].unique())
    )

def make_length_selector(df):
    return plg.DataFilterChecklist(
        header="SV Lengths",
        data_col='SVLEN_Bin',
        data_values=qc_length_bins
    )

# Filter Plugins
def make_filter_selector(df):
    filters = df['FILTER'].unique()
    data_values = []
    if 'PASS' in filters:
        data_values += ['PASS']
    if '.' in filters:
        data_values += ['.']
    data_values += sorted([f for f in filters if (f != 'PASS') and (f != '.')])
    
    return plg.DataFilterRadioButtons(
        header='Value for FILTER',
        data_col='FILTER',
        data_values=data_values
    )

# Experiment Plugins
def make_experiment_selector(df):
    return plg.DataFilterChecklist(
        header="Experimental Groups",
        data_col='Experiment',
        data_values=list(df['Experiment'].unique()),
    )

# Plot Format Plugins
def make_axes_mode_selector():
    return plg.PlotInputRadioButtons(
        header="Axis Scaling",
        plot_input='axes_mode',
        data_values=['Dynamic', 'Fixed']
    )

def make_stat_selector(data_values):
    return plg.PlotInputRadioButtons(
        header='Stat to Plot',
        plot_input='stat',
        data_values=data_values
    )


# In[ ]:





# ## Decorators

# In[ ]:


def interval_filter(plotter):
    # Decorator to wrap plotter function to filter df based on interval plugin inputs
    def interval_plotter(df, *args, **kwargs):
        interval_name = kwargs['interval_name'] if 'interval_name' in kwargs else None
        breakpoint = kwargs['breakpoint'] if 'breakpoint' in kwargs else None
        pct_overlap = kwargs['pct_overlap'] if 'pct_overlap' in kwargs else None
        
        # Resolve logic on breakpoint stats
        breakpoint_label = ''
        if breakpoint == 'Left':
            breakpoint_label = '-LBEND'
        elif breakpoint == 'Right':
            breakpoint_label == '-RBEND'
        elif breakpoint == 'Both':
            breakpoint_label = '-BBEND'
        
        # Perform filtering on df using interval name and breakpoint preference
        overlap_label = f'{interval_name}{breakpoint_label}-overlap'
        query = (df[overlap_label] >= pct_overlap[0]/100) & (df[overlap_label] <= pct_overlap[1]/100)
        sub_df = df[query]
        return plotter(sub_df, *args, **kwargs)
    return interval_plotter


# In[ ]:


def axes_mode(plotter):
    # Decorator to wrap plotter to use fixed or dynamic axes
    def axes_plotter(*args, **kwargs):
        fig = plotter(*args, **kwargs)
        axes_mode = kwargs['axes_mode'] if 'axes_mode' in kwargs else None
        if axes_mode == 'Fixed':
            fig.update_layout(yaxis_range=[0, 1])
            if COVARIATE_X is None:
                fig.update_layout(xaxis_range=[0, 1])
        return fig
    return axes_plotter


# In[ ]:


def filter_or_all_factory(col_name):
    # A decorator factory (i.e. decorator w/ args) for filtering by col_name
    def filter_or_all(plotter):
        # A decorator that filters by col_name unless value is 'ALL'
        def filter_or_all_plotter(df, *args, **kwargs):
            if kwargs[col_name] != 'ALL':
                sub_df = df[df[col_name] == kwargs[col_name]]
            else:
                sub_df = df
            return plotter(sub_df, *args, **kwargs)
        return filter_or_all_plotter
    return filter_or_all


# In[ ]:


def gt_match(plotter):
    # A decorator to wrap plotter to force gt_match or not
    # Note: assumes that plotter will only see DataFrames from Wittyer w/ WHY column
    def gt_match_plotter(df, *args, **kwargs):
        if kwargs['gt_match'] == 'True':
            sub_df = df[df['WHY'] != 'GtMismatch']
        else:
            sub_df = df
        return plotter(sub_df, *args, **kwargs)
    return gt_match_plotter


# In[ ]:





# In[ ]:





# ## Tabs

# ### QC Tabs

# #### Counts Tab

# In[ ]:


@interval_filter
def make_bar_counts(df, x, interval_name, breakpoint, pct_overlap):
    category_orders = {'Experiment': [f'{x}-{y}' for x in EXPERIMENT_ORDER for y in ['Base', 'Comp']]} if EXPERIMENT_ORDER is not None else None
    title = f'Mean Count of {x}'
    if pct_overlap[0] > 0 or pct_overlap[1] < 100:
        title += f' for events with {pct_overlap[0]:.0f}%-{pct_overlap[1]:.0f}% overlap with {interval_name}'
    
    plot_args = {
        'x': x,
        'y': 'Mean_Count',
        'title': title,
        'error_y': 'std',
        'color': 'Experiment',
        'barmode': 'group',
        'category_orders': category_orders
    }

    counts_df = df.groupby(['Experiment_Suffix', 'Experiment', 'Sample']).apply(lambda df: df[x].value_counts().reset_index()).reset_index() \
        .groupby(['Experiment', x]).apply(lambda df: df['count'].describe()).reset_index().rename(columns={'mean': 'Mean_Count'})
    
    fig = px.bar(counts_df, **plot_args) # text_auto='.2s'
    fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
    return fig


# In[ ]:


if INCLUDE_QC:
    count_plot = qbb.PlotPanel(
        header='Counts Bar Chart',
        plotter=make_bar_counts,
        plot_inputs={
            # 'x': 'SVTYPE'
        },
        data_source=qc_df,
        plugins=[
            plg.PlotInputRadioButtons(
                header='Value for x',
                plot_input='x',
                data_values=['SVTYPE', 'SVLEN_Bin', 'AF_Bin']
            ),
        ]
    )
    
    counts_tab = qbb.BaseTab(
        tab_label='SV Counts',
        tab_header='Count Plots',
        content_list=[
            count_plot
        ],
        sidebar_plugins=make_interval_plugin_bundle(qc_df) + [
            make_length_selector(qc_df),
            make_filter_selector(qc_df)
        ]
    )


# In[ ]:





# #### Counts Distributions

# In[ ]:


@interval_filter
def make_qc_histogram(df, x, barmode, interval_name, breakpoint, pct_overlap):
    category_orders = {'Experiment': [f'{x}-{y}' for x in EXPERIMENT_ORDER for y in ['Base', 'Comp']]} if EXPERIMENT_ORDER is not None else None
    title = f'Histogram of {x}'
    if pct_overlap[0] > 0 or pct_overlap[1] < 100:
        title += f' for events with {pct_overlap[0]:.0f}%-{pct_overlap[1]:.0f}% overlap with {interval_name}'
    fig = px.histogram(df, x=x, barmode=barmode, color='Experiment', category_orders=category_orders, title=title)
    return fig


# In[ ]:


if INCLUDE_QC:
    histogram_plot = qbb.PlotPanel(
        header='Histogram of Values',
        plotter=make_qc_histogram,
        plot_inputs={},
        data_source=qc_df,
        plugins=[
            plg.PlotInputRadioButtons(
                header='Value for x',
                plot_input='x',
                data_values=['SVLEN', 'QUAL', 'AF']
            ),
            plg.PlotInputRadioButtons(
                header='Bar Mode',
                plot_input='barmode',
                data_values=['group', 'overlay', 'relative']
            )
        ]
    )
    
    histogram_tab = qbb.BaseTab(
        tab_label='QC Histogram',
        tab_header='QC Histogram',
        content_list=[
            histogram_plot
        ],
        sidebar_plugins=make_interval_plugin_bundle(qc_df) + [
            make_length_selector(qc_df),
            make_filter_selector(qc_df)
        ]
    )


# In[ ]:





# In[ ]:





# #### HWE Tab

# In[ ]:


def check_single_sample_hwe(df):
    # Check to see if HWE related stats are likely from single sample, i.e should not waste time creating plot
    pass


# In[ ]:


@interval_filter
def make_hwe_plot(df, interval_name, breakpoint, pct_overlap):
    # Get type info
    type_ = df['SVTYPE'].iloc[0] if len(df['SVTYPE']) > 0 else "Empty df"
    title = f'HWE Plot for SVs of type {type_}'
    
    # Get interval info
    title = title + f' over {interval_name} (w/ {pct_overlap[0]}-{pct_overlap[1]}% overlap)<br></br>'
    
    # Apply filters
    filter_value = df['FILTER'].iloc[0] if len(df['FILTER']) > 0 else "Empty df"
    title = title + f' (FILTER is {filter_value} only)'
    
    if df['HWE'].sum() == len(df) and len(df) > 0:
        fig = go.Figure()
        fig.add_annotation(
            text="All HWE values are 1; <br></br>data likely generated from single sample so skipping HWE plot (All points at corners).",
            xref="paper", yref="paper", x=0,y=1, showarrow=False
        )
        return fig
    
    # Create plot
    fig = px.scatter_ternary(df, 
                   a='AC_Het', b='AC_Ref', c='AC_Hom', color='HWE', color_continuous_scale="haline",
                   width=700, height=700, title=title)

    fig.update_layout(coloraxis_colorbar_x=1.15)
    
    # Choose cut-off for HWE significance
    if len(df) > 0:
        absolute_cutoff = 0.05
        bonferroni_cutoff = absolute_cutoff / len(df)
        
        HWE_below_abs_cutoff = len(df[df['HWE'] < absolute_cutoff])
        HWE_below_bonferroni_cutoff = len(df[df['HWE'] < bonferroni_cutoff])
        # HWE_above_abs_cutoff = len(df[df['HWE'] >= absolute_cutoff])
        # HWE_above_bonferroni_cutoff = len(df[df['HWE'] >= bonferroni_cutoff])
        total_sites = len(df)

        fig.add_annotation(
            text=f"HWE Significant Sites (Bonferroni):<br></br> {HWE_below_bonferroni_cutoff} ({round(HWE_below_bonferroni_cutoff/total_sites*100, 2)}%)",
            xref="paper", yref="paper",
            x=0, y=.9,
            showarrow=False
        )
        
        fig.add_annotation(
            text=f"HWE Significant Sites (Absolute):<br></br> {HWE_below_abs_cutoff} ({round(HWE_below_abs_cutoff/total_sites*100, 2)}%)",
            xref="paper", yref="paper",
            x=0, y=.8,
            showarrow=False
        )
    else:
        fig.add_annotation(
            text=f"No sites fit criteria selected",
            xref="paper", yref="paper",
            x=0, y=.9,
            showarrow=False
        )
    
    return fig


# In[ ]:


if INCLUDE_QC:
    hwe_base_plot = qbb.PlotPanel(
        header='HWE Plot for Base Variants',
        plotter=make_hwe_plot,
        plot_inputs={},
        data_source=qc_df[qc_df['Experiment_Suffix'] == 'Base'],
        plugins=[]
    )
    
    hwe_comp_plot = qbb.PlotPanel(
        header='HWE Plot for Comp Variants',
        plotter=make_hwe_plot,
        plot_inputs={},
        data_source=qc_df[qc_df['Experiment_Suffix'] == 'Comp'],
        plugins=[]
    )
    
    hwe_cg_list = [hwe_base_plot]
    if not SKIP_COMP_HWE:
        hwe_cg_list += [hwe_comp_plot]
    
    hwe_cg = qbb.ContentGrid(
        header='HWE Plots',
        content_list=hwe_cg_list
    )
    
    hwe_tab = qbb.BaseTab(
        tab_label='HWE Plots',
        tab_header='Hardy-Weinberg Equilibrium Plots',
        content_list=[
            hwe_cg
        ],
        sidebar_plugins=make_interval_plugin_bundle(qc_df) + [
            make_type_selector(qc_df),
            make_length_selector(qc_df),
            make_filter_selector(qc_df),
        ]
    )


# In[ ]:





# ### Wittyer Tabs

# #### Basic Wittyer Tab

# In[ ]:


@axes_mode
def make_basic_wittyer_plot(df, axes_mode, stat='Precision'):
    if COVARIATE_X is not None:
        x = COVARIATE_X
        y = stat
        title = f'{y} vs {x} Plot'
    else:
        x = 'Recall'
        y = 'Precision'
        title = 'Precision vs Recall Plot'
    hover_data = ['TruthName', 'TruthTpCount', 'TruthFnCount', 'QueryTpCount', 'QueryFpCount']
    category_orders = {'Experiment': EXPERIMENT_ORDER} if EXPERIMENT_ORDER is not None else None
    return px.scatter(df, x=x, y=y, color='Experiment', title=title, hover_name='QueryName', hover_data=hover_data, 
                      marginal_x='box', marginal_y='box', category_orders=category_orders)


# In[ ]:


plugins = [
    plg.DataFilterRadioButtons(
        header='Resolution of Stats',
        data_col='StatsType',
        data_values=['Event', 'Base']
    )
]
if COVARIATE_X is not None:
    plugins += [make_stat_selector(['Precision', 'Recall', 'Fscore'])]


# In[ ]:


if INCLUDE_WITTYER:
    basic_wittyer_plot = qbb.PlotPanel(
        header='Precision vs Recall Plots by Bin',
        plotter=make_basic_wittyer_plot,
        plot_inputs={},
        data_source=wittyer_stats_df,
        plugins=plugins
    )
    
    sidebar_plugins = [
        make_type_selector(wittyer_stats_df),
        plg.DataFilterRadioButtons(
            header='SV Size Bin',
            data_col='Bin',
            data_values=wittyer_bins
        ),
        make_axes_mode_selector()
    ]
    
    basic_wittyer_tab = qbb.BaseTab(
        tab_label='Basic Wittyer',
        tab_header='Basic Wittyer Stats',
        content_list=[
            basic_wittyer_plot
        ],
        sidebar_plugins=sidebar_plugins
    )


# In[ ]:





# In[ ]:





# #### Adv Wittyer Tab

# In[ ]:


def add_recall(df):
    if len(df[(df['VCF'] == 'truth') & (df['WIT'] == 'TP')]['count']) > 0:
        tp = df[(df['VCF'] == 'truth') & (df['WIT'] == 'TP')]['count'].values[0]
    else:
        tp = 0
    
    if len(df[(df['WIT'] == 'FN')]['count']) > 0:
        fn = df[(df['WIT'] == 'FN')]['count'].values[0]
    else:
        fn = 0
    df['Recall'] = tp / (tp + fn) if tp + fn > 0 else np.nan
    df['TP-Base'] = tp
    df['FN'] = fn
    return df

def add_precision(df):
    if len(df[(df['VCF'] == 'query') & (df['WIT'] == 'TP')]['count']) > 0:
        tp = df[(df['VCF'] == 'query') & (df['WIT'] == 'TP')]['count'].values[0]
    else:
        tp = 0
    
    if len(df[(df['WIT'] == 'FP')]['count']) > 0:
        fp = df[(df['WIT'] == 'FP')]['count'].values[0]
    else:
        fp = 0
    df['Precision'] = tp / (tp + fp) if tp + fp > 0 else np.nan
    df['TP-Comp'] = tp
    df['FP'] = fp
    return df

def add_f1(df):
    df['F1_Score'] = 2 * df['Recall'] * df['Precision'] / (df['Recall'] + df['Precision'])
    return df


# In[ ]:


@axes_mode
@interval_filter
@filter_or_all_factory('SVTYPE')
@gt_match
def make_adv_wittyer_plot(df, interval_name, breakpoint, pct_overlap, gt_match, SVTYPE, axes_mode, stat='Precision'):
    title = 'Advanced Plot of Wittyer Stats'
    title = title + f' over {interval_name} (w/ {int(pct_overlap[0])}-{int(pct_overlap[1])}% overlap)<br></br>'
    y = stat
    
    grouping = ['TruthSample', 'QuerySample', 'Experiment']
    if COVARIATE_X is not None:
        if COVARIATE_X not in grouping:
            grouping += [COVARIATE_X]
        x = COVARIATE_X
    else:
        x = 'Recall'
    counts_df = df.groupby(grouping + ['VCF'])['WIT'].value_counts().reset_index(name='count')
    
    # Add recall/precision stats
    counts_df = counts_df.groupby(grouping).apply(add_recall).reset_index(drop=True)
    counts_df = counts_df.groupby(grouping).apply(add_precision).reset_index(drop=True)
    counts_df = counts_df.groupby(grouping).apply(add_f1).reset_index(drop=True)
    
    # Remove rows redundant for stats
    counts_df = counts_df.groupby(grouping).apply(lambda df: df.iloc[0]).reset_index(drop=True)
    if len(counts_df) == 0:
        # counts_df.columns += ['Precision', 'Recall', 'F1_Score']    # Fix unhelpful error message when filters reduce to nothing
        fig = go.Figure()
        fig.update_layout(title='No variants found with given conditions')
        return fig
    
    hover_data = ['TruthSample', 'TP-Base', 'FN', 'TP-Comp', 'FP']
    category_orders = {'Experiment': EXPERIMENT_ORDER} if EXPERIMENT_ORDER is not None else None
    fig = px.scatter(counts_df, x=x, y=y, title=title, color='Experiment', hover_name='QuerySample', hover_data=hover_data, 
                     marginal_x='box', marginal_y='box', category_orders=category_orders)
    return fig


# In[ ]:


plugins = [
    plg.PlotInputRadioButtons(
    header='Force Match GT',
    plot_input='gt_match',
    data_values=['False', 'True']
)]

if COVARIATE_X is not None:
    plugins += [make_stat_selector(['Precision', 'Recall', 'F1_Score'])]


# In[ ]:


if INCLUDE_WITTYER:
    adv_wittyer_plot = qbb.PlotPanel(
        header='Advanced Wittyer Plot',
        plotter=make_adv_wittyer_plot,
        plot_inputs={},
        data_source=adv_wittyer_df,
        plugins=plugins
    )
    
    sidebar_plugins = make_interval_plugin_bundle(adv_wittyer_df) + [
        plg.PlotInputRadioButtons(
            header="Variant Type",
            plot_input='SVTYPE',
            data_values=['ALL'] + list(adv_wittyer_df['SVTYPE'].unique())
        ),
        make_axes_mode_selector(),
        make_filter_selector(adv_wittyer_df)
    ]
    
    adv_wittyer_tab = qbb.BaseTab(
        tab_label='Adv Wittyer',
        tab_header='Adv Wittyer Plots',
        content_list=[
            adv_wittyer_plot
        ],
        sidebar_plugins=sidebar_plugins
    )


# In[ ]:





# In[ ]:





# ### Truvari Tabs

# #### Truvari Bench Tab

# In[ ]:


@axes_mode
def make_truvari_bench_plot(df, axes_mode, stat='precision'):
    if COVARIATE_X is not None:
        x = COVARIATE_X
        y = stat
        title = f'{y} vs {x} Plot'
    else:
        x = 'recall'
        y = 'precision'
        title = 'Precision vs Recall Plot'
    
    type_ = df['SVTYPE'].iloc[0] if len(df) > 0 else "Empty df"
    title += f' for SVTYPE {type_}'
    
    interval = df['Interval'].iloc[0] if len(df) > 0 else "Empty df"
    title += f' over {interval}'
    
    hover_data = ['Base_Name', 'TP-base', 'FN', 'TP-comp', 'FP']
    category_orders = {'Experiment': EXPERIMENT_ORDER} if EXPERIMENT_ORDER is not None else None
    return px.scatter(df, x=x, y=y, color='Experiment', title=title, hover_name='Comp_Name', hover_data=hover_data, 
                      marginal_x='box', marginal_y='box', category_orders=category_orders)
    # return px.scatter(df, x='recall', y='precision', color='Comp_Name')


# In[ ]:


def make_gt_concordance_plot(df, axes_mode):
    category_orders = {'Experiment': EXPERIMENT_ORDER} if EXPERIMENT_ORDER is not None else None
    title = 'GT Concordance'
    if len(df) > 0:
        interval_name = df['Interval'].values[0]
        title += f' over {interval_name}'
    fig = px.box(df[df['TP-comp_TP-gt'] > 0], x='SVTYPE', y='gt_concordance', color='Experiment', 
                 title=title, category_orders=category_orders)
    if axes_mode == 'Fixed':
        fig.update_layout(yaxis_range=[0, 1])
    return fig


# In[ ]:


if INCLUDE_TRUVARI:
    truvari_bench_plot = qbb.PlotPanel(
        header='Truvari Benchmarks',
        plotter=make_truvari_bench_plot,
        plot_inputs={},
        data_source=truvari_bench_df,
        plugins=[make_type_selector(truvari_bench_df), make_stat_selector(['precision', 'recall', 'f1'])] if COVARIATE_X is not None \
            else [make_type_selector(truvari_bench_df)]
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
            make_axes_mode_selector()
        ]
    )


# In[ ]:





# #### Truvari Errors Tab

# In[ ]:


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


# In[ ]:


if INCLUDE_TRUVARI:
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


# In[ ]:





# ## Main Board

# In[ ]:


tab_list = []
if INCLUDE_QC:
    tab_list += [counts_tab, histogram_tab, hwe_tab]
if INCLUDE_WITTYER:
    tab_list += [basic_wittyer_tab, adv_wittyer_tab]
if INCLUDE_TRUVARI:
    tab_list += [truvari_bench_tab, truvari_errors_tab]

board = qbb.Quickboard(
    tab_list=tab_list
)


# In[ ]:


start_app(board, debug=False, app_title='SVisualizer', mode='external', port=8050)


# In[ ]:





# In[ ]:





# In[ ]:




