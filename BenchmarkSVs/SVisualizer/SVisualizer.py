#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff

import quickboard.base as qbb
import quickboard.plugins as plg
from quickboard.app import start_app, deploy_app, app

from upset_plot_utils import create_upset, make_disqualified_df


# In[ ]:





# # User Config
# 
# Use this section to set variables that change the output and settings of the dashboard. Anything after this header's section should be untouched, and modifying the following code is advanced usage.

# In[ ]:





# In[ ]:





# In[ ]:





# # Import Data

# ## File Paths

# In[2]:


QC_DATA_PATH = 'wdl_outputs/qc_stats.tsv'

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


# In[ ]:





# ## QC Data

# In[3]:


qc_df = pd.read_csv(QC_DATA_PATH, sep='\t')


# In[4]:


# Add AC_Ref counts
qc_df['AC_Ref'] = qc_df['NS'] - qc_df['AC_Het'] - qc_df['AC_Hom'] / 2


# In[5]:


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

qc_df['SVLEN_Bin'] = qc_df['SVLEN'].replace('.', '0').astype(int).apply(bin_length)


# In[6]:


# Bin AFs
qc_df['AF'] = qc_df['AF'].replace('.', 0).astype(float)
qc_df['AF_Bin'] = qc_df['AF'].replace('.', 0).astype(float)

qc_af_bins = [
    'AC=1',
    '< 1%',
    '1-10%',
    '10-50%',
    '> 50%'
]

def bin_af(row):
    if row.AF == 0:
        return row
    elif row.AC == 1:
        row['AF_Bin'] = 'AC=1'
        return row
    elif row.AF < .01:
        row['AF_Bin'] = '< 1%'
        return row
    elif row.AF < .1:
        row['AF_Bin'] = '1-10%'
        return row
    elif row.AF < .5:
        row['AF_Bin'] = '10-50%'
        return row
    else:
        row['AF_Bin'] = '> 50%'
        return row

qc_df = qc_df.apply(bin_af, axis=1)


# In[ ]:





# ## Truvari Data

# ### Bench Data

# In[7]:


truvari_bench_df = pd.read_csv(TRUVARI_BENCH_PATH, sep='\t')

truvari_bench_df = truvari_bench_df.rename(columns={'SV_Type': 'SVTYPE'})


# In[ ]:





# ### Closest Data

# In[8]:


truvari_fn_closest_df = pd.read_csv(TRUVARI_FN_CLOSEST_PATH, sep='\t')
truvari_fp_closest_df = pd.read_csv(TRUVARI_FP_CLOSEST_PATH, sep='\t')


# In[9]:


# Add SIZE_RATIO fields
truvari_fp_closest_df['SIZE_RATIO'] = truvari_fp_closest_df[['LLEN', 'RLEN']].min(axis=1) / truvari_fp_closest_df[['LLEN', 'RLEN']].max(axis=1)
truvari_fn_closest_df['SIZE_RATIO'] = truvari_fn_closest_df[['LLEN', 'RLEN']].min(axis=1) / truvari_fn_closest_df[['LLEN', 'RLEN']].max(axis=1)


# In[10]:


# Add kth closest labels
# TODO: Infer from output dataframe, or include from WDL outputs...
k = 3
truvari_fp_closest_df['kth_closest'] = truvari_fp_closest_df.index % k + 1
truvari_fn_closest_df['kth_closest'] = truvari_fn_closest_df.index % k + 1


# In[ ]:





# ### Intervals Data

# In[11]:


truvari_fn_intervals_df = pd.read_csv(TRUVARI_FN_INTERVALS_PATH, sep='\t')
truvari_fp_intervals_df = pd.read_csv(TRUVARI_FP_INTERVALS_PATH, sep='\t')
truvari_tpbase_intervals_df = pd.read_csv(TRUVARI_TP_BASE_INTERVALS_PATH, sep='\t')
truvari_tpcomp_intervals_df = pd.read_csv(TRUVARI_TP_COMP_INTERVALS_PATH, sep='\t')


# In[ ]:





# ## Wittyer Data

# In[12]:


wittyer_stats_df = pd.read_csv(WITTYER_STATS_PATH, sep='\t')
wittyer_truth_df = pd.read_csv(WITTYER_TRUTH_PATH, sep='\t')
wittyer_query_df = pd.read_csv(WITTYER_QUERY_PATH, sep='\t')
wittyer_nogt_df = pd.read_csv(WITTYER_NOGT_PATH, sep='\t')


# In[ ]:





# In[ ]:





# ## Interval Processing

# In[13]:


# Infer interval labels from Truvari df
# TODO: Make this a bit more robust...
interval_labels = [c.split('-count')[0] for c in truvari_fn_intervals_df.columns if '-count' in c]


# In[14]:


interval_columns = [f'{x}-{s}' for x in interval_labels for s in ('count', 'overlap')]


# In[15]:


# Add interval labels to qc_df
# TODO: Make this still work without Wittyer files...

query_interval_df = pd.concat([
    wittyer_query_df[['CHROM', 'POS', 'END'] + interval_columns],
    wittyer_nogt_df[['CHROM', 'POS', 'END'] + interval_columns]
]).drop_duplicates()

qc_df = pd.merge(qc_df, query_interval_df, on=['CHROM', 'POS', 'END'])


# In[ ]:





# In[ ]:





# # Quickboard

# ## Plugin Methods

# In[16]:


def make_interval_selector(df):
    return plg.DataFilterRadioButtons(
        header="Interval List",
        data_col='Interval',
        data_values=list(df['Interval'].unique())
    )

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

def make_experiment_selector(df):
    return plg.DataFilterChecklist(
        header="Experimental Groups",
        data_col='Experiment',
        data_values=list(df['Experiment'].unique()),
    )

def make_filter_selector(df):
    return plg.DataFilterRadioButtons(
        header='Value for FILTER',
        data_col='FILTER',
        data_values=['PASS'] + sorted([f for f in df['FILTER'].unique() if f != 'PASS'])
    )


# In[ ]:





# ## Tabs

# ### Counts Tab

# In[17]:


# TODO: Make sure x-axis is sorted appropriately for various categories, e.g. AF
def make_bar_counts(df, x):
    counts_df = df[x].value_counts().reset_index(name='count').rename(columns={'index': x})
    fig = px.bar(counts_df, x=x, y='count', barmode='group', title=f'Count of {x}') # text_auto='.2s'
    fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
    return fig


# In[18]:


count_plot = qbb.PlotPanel(
    header='Counts Bar Chart',
    plotter=make_bar_counts,
    plot_inputs={
        'x': 'SVTYPE'
    },
    data_source=qc_df,
    plugins=[
        plg.PlotInputRadioButtons(
            header='Value for x',
            plot_input='x',
            data_values=['SVTYPE', 'SVLEN_Bin', 'AF_Bin']
        ),
        make_filter_selector(qc_df)
    ]
)


# In[19]:


counts_tab = qbb.BaseTab(
    tab_label='SV Counts',
    tab_header='Count Plots',
    content_list=[
        count_plot
    ],
    sidebar_plugins=[
        # make_interval_selector(qc_df),
        make_length_selector(qc_df),
        # make_experiment_selector(qc_df)
    ]
)


# In[ ]:





# ### HWE Tab

# In[20]:


def make_hwe_plot(df, interval, min_pct_overlap):
    # Get type info
    type_ = df['SVTYPE'].iloc[0] if len(df['SVTYPE']) > 0 else "Empty df"
    title = f'HWE Plot for SVs of type {type_}'
    
    # Get interval info
    # interval = df['Interval'].iloc[0] if len(df['Interval']) > 0 else "Empty df"
    df = df[(df[f'{interval}-overlap'] >= min_pct_overlap[0]) & (df[f'{interval}-overlap'] <= min_pct_overlap[1])]
    title = title + f' over {interval} (w/ {100*min_pct_overlap[0]}-{100*min_pct_overlap[1]}% overlap)<br></br>'
        
    # Apply filters
    filter_value = df['FILTER'].iloc[0] if len(df['FILTER']) > 0 else "Empty df"
    title = title + f' (FILTER is {filter_value} only)'
    
    # Create plot
    fig = px.scatter_ternary(df, 
                   a='AC_Het', b='AC_Ref', c='AC_Hom', color='HWE', color_continuous_scale="haline",
                   width=800, height=800, title=title)

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


# In[21]:


hwe_plot = qbb.PlotPanel(
    header='HWE Plot',
    plotter=make_hwe_plot,
    plot_inputs={},
    data_source=qc_df,
    plugins=[
        make_filter_selector(qc_df),
        plg.PlotInputRangeSlider(
            header='Min Pct Overlap w/ Chosen Region',
            plot_input='min_pct_overlap',
            slider_min=0,
            slider_max=1,
            slider_default_values=[0,1],
            slider_step=0.05,
            updatemode='mouseup'
        )
    ]
)


# In[22]:


hwe_tab = qbb.BaseTab(
    tab_label='HWE Plots',
    tab_header='Hardy-Weinberg Equilibrium Plots',
    content_list=[
        hwe_plot
    ],
    sidebar_plugins=[
        plg.PlotInputRadioButtons(
            header='Interval List',
            plot_input='interval',
            data_values=interval_labels
        ),
        # make_interval_selector(qc_df),
        make_type_selector(qc_df),
        make_length_selector(qc_df),
        # make_experiment_selector(qc_df)
    ]
)


# In[ ]:





# ### Basic Wittyer Tab

# In[23]:


wittyer_bins = ['All'] + sorted([str(x) for x in wittyer_stats_df['Bin'].unique() if (x != 'All') & (str(x) != 'nan')]) + ['nan']


# In[24]:


basic_wittyer_plot = qbb.PlotPanel(
    header='Precision vs Recall Plots by Bin',
    plotter=px.scatter,
    plot_inputs={
        'x': 'Recall',
        'y': 'Precision',
        'color': 'Experiment',
        'title': 'Precision vs Recall Plot'
    },
    data_source=wittyer_stats_df,
    plugins=[
        plg.DataFilterRadioButtons(
            header='Resolution of Stats',
            data_col='StatsType',
            data_values=['Event', 'Base']
        )
    ]
)


# In[25]:


basic_wittyer_tab = qbb.BaseTab(
    tab_label='Basic Wittyer',
    tab_header='Basic Wittyer Stats',
    content_list=[
        basic_wittyer_plot
    ],
    sidebar_plugins=[
        make_type_selector(wittyer_stats_df),
        plg.DataFilterRadioButtons(
            header='SV Size Bin',
            data_col='Bin',
            data_values=wittyer_bins
        )
    ]
)


# In[ ]:





# In[ ]:





# ### Adv Wittyer Tab

# In[26]:


wittyer_query_df['VCF'] = 'query'
wittyer_truth_df['VCF'] = 'truth'

adv_wittyer_df = pd.concat([wittyer_query_df, wittyer_truth_df])

# counts_df = adv_wittyer_df.groupby(['TruthSample', 'QuerySample', 'VCF'])['WIT'].value_counts().reset_index(name='count')

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
    return df


# In[27]:


def make_adv_wittyer_plot(df, interval, min_pct_overlap, gt_match, sv_type):
    title = 'Advanced Plot of Wittyer Stats'
    df = df[(df[f'{interval}-overlap'] >= min_pct_overlap[0]) & (df[f'{interval}-overlap'] <= min_pct_overlap[1])]
    title = title + f' over {interval} (w/ {int(100*min_pct_overlap[0])}-{int(100*min_pct_overlap[1])}% overlap)<br></br>'
    
    if gt_match == 'True':
        df = df[df['WHY'] != 'GtMismatch']
        
    if sv_type != 'ALL':
        df = df[df['SVTYPE'] == sv_type]
        
    counts_df = df.groupby(['TruthSample', 'QuerySample', 'VCF'])['WIT'].value_counts().reset_index(name='count')
    
    # Add recall/precision stats
    counts_df = counts_df.groupby(['TruthSample', 'QuerySample']).apply(add_recall).reset_index(drop=True)
    counts_df = counts_df.groupby(['TruthSample', 'QuerySample']).apply(add_precision).reset_index(drop=True)
    
    # Remove rows redundant for stats
    counts_df = counts_df.groupby(['TruthSample', 'QuerySample']).apply(lambda df: df.iloc[0]).reset_index(drop=True)
    
    fig = px.scatter(counts_df, x='Recall', y='Precision', title=title)
    return fig


# In[28]:


adv_wittyer_plot = qbb.PlotPanel(
    header='Advanced Wittyer Plot',
    plotter=make_adv_wittyer_plot,
    plot_inputs={
        
    },
    data_source=adv_wittyer_df,
    plugins=[
        make_filter_selector(adv_wittyer_df),
        plg.PlotInputRangeSlider(
            header='Min Pct Overlap w/ Chosen Region',
            plot_input='min_pct_overlap',
            slider_min=0,
            slider_max=1,
            slider_default_values=[0,1],
            slider_step=0.05,
            updatemode='mouseup'
        )
    ]
)


# In[29]:


adv_wittyer_tab = qbb.BaseTab(
    tab_label='Adv Wittyer',
    tab_header='Adv Wittyer Plots',
    content_list=[
        adv_wittyer_plot
    ],
    sidebar_plugins=[
        plg.PlotInputRadioButtons(
            header='Interval List',
            plot_input='interval',
            data_values=interval_labels
        ),
        plg.PlotInputRadioButtons(
            header="Variant Type",
            plot_input='sv_type',
            data_values=['ALL'] + list(adv_wittyer_df['SVTYPE'].unique())
        ),
        plg.PlotInputRadioButtons(
            header='Force Match GT',
            plot_input='gt_match',
            data_values=['False', 'True']
        )
    ]
)


# In[ ]:





# In[ ]:





# ### Truvari Bench Tab

# In[30]:


truvari_bench_plot = qbb.PlotPanel(
    header='Truvari Benchmarks',
    plotter=px.scatter,
    plot_inputs={
        'x': 'precision',
        'y': 'recall',
        'color': 'Comp_Name'
    },
    data_source=truvari_bench_df,
    plugins=[
        
    ]
)


# In[31]:


truvari_bench_tab = qbb.BaseTab(
    tab_label='Truvari Bench',
    tab_header='',
    content_list=[
        truvari_bench_plot
    ],
    sidebar_plugins=[
        make_interval_selector(truvari_bench_df),
        make_type_selector(truvari_bench_df)
    ]
)


# In[ ]:





# ### Truvari Errors Tab

# In[32]:


def make_closest_plot(df, sort_by, asc, mode):
    title = f'Counts of Disqualified (FP) Sites (N = {len(df)})'
    asc = asc == 'Ascending'
    disq_df = make_disqualified_df(df, dist_threshold=500, size_ratio_threshold=0.7, color='kth_closest')
    fig = create_upset(disq_df, title=title, sort_by=sort_by, asc=asc, mode=mode, color='kth_closest')
    return fig


# In[33]:


fp_plot = qbb.PlotPanel(
    header='FP Stats',
    plotter=make_closest_plot,
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
            data_values=['Counts', 'Percent']
        )
    ],
    plugin_wrap=3
)

fn_plot = qbb.PlotPanel(
    header='FN Stats',
    plotter=make_closest_plot,
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
            data_values=['Counts', 'Percent']
        )
    ],
    plugin_wrap=3
)


# In[34]:


truvari_errors_tab = qbb.BaseTab(
    tab_label='Truvari Errors',
    tab_header='Truvari Stats on Mismatched Variants',
    content_list=[
        fp_plot,
        fn_plot
    ],
    sidebar_plugins=[
        
    ]
)


# In[ ]:





# In[ ]:





# ## Main Board

# In[35]:


board = qbb.Quickboard(
    tab_list=[
        counts_tab,
        hwe_tab,
        basic_wittyer_tab,
        adv_wittyer_tab,
        truvari_bench_tab,
        truvari_errors_tab,
    ]
)


# In[36]:


start_app(board, debug=False, app_title='SVisualizer', mode='external', port=8050)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




