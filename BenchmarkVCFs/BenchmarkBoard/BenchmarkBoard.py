#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import plotly.express as px

import quickboard.base as qbb
import quickboard.plugins as plg
from quickboard.app import start_app


# In[ ]:





# In[ ]:





# ## User Inputs

# In[2]:


# Users must specify column names for quantities to use in making some plots below. See the docs for details.
# Toggle optional values here to affect behavior of generated Dashboard.

# Fill with column name (or list of column names) to plot x-axis vs benchmark stats, e.g. 'Mean_Coverage' 
# Users must edit data files to provide the data, or include in pipeline 'Experiment' or 'Extra_Column' values
COVARIATE_X = None

# Display only usual 6 variant types as options in variant type selection; Toggle to False to include other categories like MNP, MA, etc.
SIMPLE_VARIANTS_ONLY = True

# Manually set order for experiment groups to appear in plots
EXPERIMENT_ORDER = None

# Color map to use for manually coloring experiment categories
EXPERIMENT_COLOR_MAP = None

# Toggle whether to show error bars on multi-sample plots or not
SHOW_ERROR_BARS = True


# In[ ]:





# In[3]:


# DO NOT MODIFY ANY CODE BELOW HERE


# In[ ]:





# ## Load Data

# In[4]:


# These file paths must point to the corresponding output files from the WDL, either saved locally or bucket links

summary_df = pd.read_csv('SimpleSummary.tsv', sep='\t')
roc_df = pd.read_csv('ROCStats.tsv', sep='\t')
st_df = pd.read_csv('SNPSubstitutionStats.tsv', sep='\t')
idd_df = pd.read_csv('IndelDistributionStats.tsv', sep='\t')


# In[ ]:





# ### Fill Null Interval Entries

# In[5]:


for df in [idd_df, st_df, summary_df]:
    strat_values = df['Interval'].unique()
    if 'Whole Genome' not in strat_values:
        df['Interval'] = df['Interval'].fillna('Whole Genome')
    elif 'Whole Genome (default)' in strat_values:
        raise ValueError("Error: Rename your Interval labels to not include either 'Whole Genome' or 'Whole Genome (default)'")
    else:
        df['Interval'] = df['Interval'].fillna('Whole Genome (default)')


# In[ ]:





# ### Fill Experiment Column if not Provided

# In[6]:


for df in [idd_df, roc_df, st_df, summary_df]:
    if 'Experiment' not in df.columns:
        df['Experiment'] = 'No_ExpGroups_Provided'


# In[ ]:





# ### Other Environment Variables

# In[7]:


# Check if only one distinguishable sample in summary_df
# Affects behavior of some plots downstream
SINGLE_SAMPLE_MODE = len(summary_df[['Experiment', 'Query_Name', 'Base_Name']].value_counts()) == 1


# In[8]:


CATEGORY_ORDERS = {
    'Experiment': EXPERIMENT_ORDER
} if EXPERIMENT_ORDER is not None else None


# In[ ]:





# ## Plugin Utilities

# In[9]:


simple_variants = ['SNP', 'HetSNP', 'HomVarSNP', 'INDEL', 'HetINDEL', 'HomVarINDEL']

def make_strat_selector(df):
    return plg.DataFilterRadioButtons(
        header="Interval List",
        data_col='Interval',
        data_values=list(df['Interval'].unique())
    )

def make_type_selector(df):
    if SIMPLE_VARIANTS_ONLY:
        variant_values = [x for x in simple_variants if x in df['Type'].unique()]
    else:
        variant_values = list(df['Type'].unique())
    return plg.DataFilterRadioButtons(
        header="Variant Type",
        data_col='Type',
        data_values=variant_values
    )

def make_sample_selector(df):
    return plg.DataFilterChecklist(
        header="Query Sample",
        data_col='Query_Name',
        data_values=list(df['Query_Name'].unique())
    )

def make_experiment_selector(df):
    return plg.DataFilterChecklist(
        header="Experimental Groups",
        data_col='Experiment',
        data_values=list(df['Experiment'].unique()),
    )

def make_axes_mode_selector():
    return plg.PlotInputRadioButtons(
        header='Axes Mode',
        plot_input='axes_mode',
        data_values=['Dynamic', 'Fixed']
    )


# In[ ]:





# ## Data Utilities

# In[10]:


def make_exp_average(df, group):
    cols = ['TP_Base', 'TP_Query', 'FP', 'FN', 'F1_Score', 'Precision', 'Recall']
    
    # Handle case where no IGN/OUT cols like ROC plots by adding back when present
    if 'IGN' in df.columns:
        cols += ['IGN']
    if 'OUT' in df.columns:
        cols += ['OUT']

    df_means = df.groupby(['Experiment', group])[cols].mean().reset_index()
    # Use naive noise model for error bars
    df_conf = df.groupby(['Experiment', group])[cols].sem().apply(lambda x: 1.96*x).reset_index()
    
    df_means = df_means.round(2)
    df_conf = df_conf.round(4)

    plot_df = df_means.merge(df_conf, on=['Experiment', group], suffixes=('_mean', '_conf'))
    return plot_df


# In[ ]:





# ## Summary Tab

# In[11]:


def make_prec_recall_plot(df, marginal, axes_mode):
    strat = df['Interval'].iloc[0]
    type_ = df['Type'].iloc[0]
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    if not SINGLE_SAMPLE_MODE:
        marginal = marginal.lower() if marginal != 'None' else None
        fig = px.scatter(df, x='Recall', y='Precision', color=color, marginal_x=marginal, marginal_y=marginal,
                          hover_data=['Query_Name'], title=f'Precision vs Recall Plot over {strat} for {type_}', 
                          category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
        if axes_mode == 'Fixed':
            fig.update_layout(xaxis_range=[0, 1.1], yaxis_range=[0, 1.1])
    else:
        melted_df = df.melt(id_vars=['Experiment', 'Query_Name', 'Base_Name', 'Interval', 'Type'], value_vars=['Precision', 'Recall', 'F1_Score'])
        melted_df = melted_df.rename(columns={'variable': 'Stat', 'value': 'Value'})
        fig = px.bar(melted_df, x='Stat', y='Value', title=f'Performance Stats over {strat} for {type_}', 
                     category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
        fig.update_layout(yaxis_range=[0, 1.1])
    
    return fig


def make_stat_covariate_plot(df, covaraite, stat, axes_mode):
    strat = df['Interval'].iloc[0]
    type_ = df['Type'].iloc[0]
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    fig = px.scatter(df, x=stat_corr, y=stat, color=color, hover_data=['Query_Name', 'TP_Base', 'TP_Query', 'FP', 'FN'],
                      title=f'Plot of {stat} by {covariate} over {strat} for {type_}', 
                      category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
    if axes_mode == 'Fixed':
        fig.update_layout(yaxis_range=[0, 1.1])
    return fig


# In[12]:


summary_sidebar_plugins = [
    make_strat_selector(summary_df),
    make_type_selector(summary_df),
    make_experiment_selector(summary_df),
    make_axes_mode_selector()
]

# Prec vs Recall plot 

prec_recall_plot = qbb.PlotPanel(
    header="Precision vs Recall Plot",
    plotter=make_prec_recall_plot,
    plot_inputs={
        'marginal': 'Box'
    },
    data_source=summary_df,
    plugins=[
        plg.PlotInputRadioButtons(
            header='Marginal Plot Type',
            plot_input='marginal',
            data_values=['None', 'Box', 'Violin', 'Histogram', 'Rug']
        ),
    ] if not SINGLE_SAMPLE_MODE else []
)


# Stat correlator plot

stat_covariate_plugins = [
    plg.PlotInputRadioButtons(
        header='Statistic to Plot',
        plot_input='stat',
        data_values=['F1_Score', 'Precision', 'Recall']
    ),
]

if COVARIATE_X is not None:
    if isinstance(COVARIATE_X, str):
        correlators = [COVARIATE_X]
    else:
        try: 
            assert isinstance(COVARIATE_X, list)
            correlators = COVARIATE_X
            stat_corr_plugins += [
                plg.PlotInputRadioButtons(
                    header='x-axis Covariate to Plot',
                    plot_input='covariate',
                    data_values=correlators
                )
            ]
        except AssertionError:
            print('COVARIATE_X must be a string or list of strings!')

    stat_covariate_plot = qbb.PlotPanel(
        header="Stat Correlation Plot",
        plotter=make_stat_covariate_plot,
        plot_inputs={
            'stat_corr': correlators[0],
            'stat': 'F1_Score'
        },
        data_source=summary_df,
        plugins=stat_covariate_plugins
    )

bench_stat_cg = qbb.ContentGrid(
    header='Benchmarking Stat Scatter Plots',
    content_list=[prec_recall_plot] + ([stat_covariate_plot] if COVARIATE_X is not None else [])
)


summary_tab = qbb.BaseTab(
    tab_label="Simple Summary",
    tab_header="Plot for Simple Summary Stats",
    content_list=[
        bench_stat_cg
    ],
    sidebar_plugins=summary_sidebar_plugins
)


# In[ ]:





# ## ROC Tab

# In[13]:


def make_roc_plot(df, roc_mode, error_bars, axes_mode):
    type_ = df['Type'].iloc[0]
    if roc_mode in ['Precision', 'Recall']:
        error_y = None
        if not SINGLE_SAMPLE_MODE:
            plot_df = make_exp_average(df, 'Score')
            error_y = f'{roc_mode}_conf' if SHOW_ERROR_BARS & (error_bars == 'Show') else None
        else:
            plot_df = df
    
        fig = px.line(plot_df, x='Score', y=f'{roc_mode}_mean', error_y=error_y, 
                      hover_data=['Score', 'TP_Base_mean', 'TP_Query_mean', 'FP_mean', 'FN_mean'],
                      title=f'{roc_mode} Plot for {type_} stratified by Score', color='Experiment',
                      category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
    else:
        fig = px.line(df, x='Recall', y='Precision', color='Experiment', line_group='Query_Name', hover_data=['Score'],
                      title=f'ROC Plot for {type_} stratified by Score')
        fig.update_layout(xaxis_range=[0, 1.1])

    if axes_mode == 'Fixed':
        fig.update_layout(yaxis_range=[0, 1.1])
    return fig


# In[14]:


roc_sidebar_plugins = [
    make_type_selector(roc_df),
    make_experiment_selector(roc_df),
    make_axes_mode_selector(),
    make_sample_selector(roc_df),
]

roc_plot_plugins = [
    plg.PlotInputRadioButtons(
        header='ROC Mode',
        plot_input='roc_mode',
        data_values=['Traditional', 'Precision', 'Recall']
    )
]
roc_plot_plugins += [
    plg.PlotInputRadioButtons(
        header='Toggle Error Bar Visibility',
        plot_input='error_bars',
        data_values=['Show', 'Hide']
    )
] if SHOW_ERROR_BARS else []

roc_plot = qbb.PlotPanel(
    header="ROC Plot (TP vs FP by Score Field)",
    plotter=make_roc_plot,
    plot_inputs={},
    data_source=roc_df,
    plugins=roc_plot_plugins
)


# In[15]:


roc_tab = qbb.BaseTab(
    tab_label='ROC Plots',
    tab_header='ROC Plots',
    content_list=[
        roc_plot
    ],
    sidebar_plugins=roc_sidebar_plugins
)


# In[ ]:





# ## SNP Tab

# In[16]:


st_df['Ref_Nucleotide'] = st_df['Substitution'].apply(lambda x: x.split('>')[0])
st_df['Var_Nucleotide'] = st_df['Substitution'].apply(lambda x: x.split('>')[1])


# In[17]:


def make_titv_plot(df, stat, axes_mode):
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    strat = df['Interval'].iloc[0]
    type_ = df['Type'].iloc[0]

    error_y = None
    if not SINGLE_SAMPLE_MODE:
        plot_df = make_exp_average(df, 'Substitution_Type')
        error_y = f'{stat}_conf' if SHOW_ERROR_BARS else None
    else:
        plot_df = df

    plot_df = plot_df.replace('Ti', 'Transitions').replace('Tv', 'Transversions')

    fig = px.bar(plot_df, x='Substitution_Type', y=f'{stat}_mean', error_y=error_y, 
                  title=f'Plot of {stat} mean by Substitution Type on {strat} for {type_}', 
                  hover_data=['TP_Base_mean', 'TP_Query_mean', 'FP_mean', 'FN_mean', 'IGN_mean', 'OUT_mean'], color=color, barmode='group',
                  category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
    if axes_mode == 'Fixed' and stat in ['Precision', 'Recall', 'F1_Score']:
        fig.update_layout(yaxis_range=[0, 1.1])
    return fig

def make_snp_substitution_plot(df, stat, axes_mode):
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    strat = df['Interval'].iloc[0]
    type_ = df['Type'].iloc[0]
    df_means = df.groupby(['Experiment', 'Type', 'Substitution_Type', 
                           'Ref_Nucleotide', 'Var_Nucleotide'])[['TP_Base', 'TP_Query', 
                                                                 'FP', 'FN', 'F1_Score', 'Precision', 'Recall', 'IGN', 'OUT']].mean().reset_index()
    df_conf = df.groupby(['Experiment', 'Type', 'Substitution_Type', 'Ref_Nucleotide', 'Var_Nucleotide'])[[
        'TP_Base', 'TP_Query', 'FP', 'FN', 'F1_Score', 'Precision', 'Recall', 'IGN', 'OUT']].sem().apply(lambda x: 1.96*x).reset_index()
    
    plot_df = df_means.merge(df_conf, on=['Experiment', 'Ref_Nucleotide', 'Var_Nucleotide'], suffixes=('_mean', '_conf'))
    counts = ['TP_Base_mean', 'TP_Query_mean', 'FP_mean', 'FN_mean']
    plot_df[counts] = plot_df[counts].round(2)
    
    category_orders = {
        'Ref_Nucleotide': ['A', 'G', 'C', 'T'],
        'Var_Nucleotide': ['A', 'G', 'C', 'T']
    }
    if EXPERIMENT_ORDER is not None:
        category_orders = {**category_orders, **{'Experiment': EXPERIMENT_ORDER}}

    fig = px.scatter_3d(plot_df, x='Ref_Nucleotide', y='Var_Nucleotide', z=f'{stat}_mean', error_z=f'{stat}_conf', color=color, 
                        hover_data=['TP_Base_mean', 'TP_Query_mean', 'FP_mean', 'FN_mean', 'IGN_mean', 'OUT_mean'],
                        title=f'Plot of {stat} per Substitution Type on {strat} for {type_}', 
                        category_orders=category_orders, symbol='Substitution_Type_mean', color_discrete_map=EXPERIMENT_COLOR_MAP,
                        height=700, width=1000,
                       )
    fig.update_layout(scene_camera=dict(eye=dict(x=1.5, y=1.5, z=0.4)))
    return fig


# In[18]:


snp_sidebar_plugins = [
    make_strat_selector(st_df),
    make_type_selector(st_df),
    make_experiment_selector(st_df),
    make_axes_mode_selector()
]

snp_titv_plot = qbb.PlotPanel(
    header='Distribution of Stats for Ti or Tv',
    plotter=make_titv_plot,
    plot_inputs={
        'stat': 'F1_Score',
    },
    data_source=st_df,
    plugins=[
        plg.PlotInputRadioButtons(
            header='Stat to Plot',
            plot_input='stat',
            data_values=['F1_Score', 'Precision', 'Recall', 'TP_Base', 'TP_Query', 'FP', 'FN', 'IGN', 'OUT']
        ),
    ]
)

snp_substitution_plot = qbb.PlotPanel(
    header='SNP Substitution Plot',
    plotter=make_snp_substitution_plot,
    plot_inputs={
        'stat': 'F1_Score',
    },
    data_source=st_df,
    plugins=[
        plg.PlotInputRadioButtons(
            header='Stat to Plot',
            plot_input='stat',
            data_values=['F1_Score', 'Precision', 'Recall', 'TP_Base', 'TP_Query', 'FP', 'FN', 'IGN', 'OUT']
        )
    ]
)

snp_cg = qbb.ContentGrid(
    content_list=[
        snp_titv_plot,
        snp_substitution_plot
    ],
    col_wrap=1
)

snp_tab = qbb.BaseTab(
    tab_label='SNP Plots',
    tab_header='SNP Plots',
    content_list=[
        snp_cg
    ],
    sidebar_plugins=snp_sidebar_plugins
)


# In[ ]:





# ## INDEL Tab

# In[19]:


def make_idd_plot(df, stat, axes_mode):
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    strat = df['Interval'].iloc[0]
    type_ = df['Type'].iloc[0]

    error_y = None
    if not SINGLE_SAMPLE_MODE:
        plot_df = make_exp_average(df, 'INDEL_Length')
        error_y = f'{stat}_conf' if SHOW_ERROR_BARS else None
    else:
        plot_df = df
    
    fig = px.line(plot_df, x='INDEL_Length', y=f'{stat}_mean', error_y=error_y, title=f'Plot of {stat} mean by INDEL Length on {strat} for {type_}', 
                  hover_data=['TP_Base_mean', 'TP_Query_mean', 'FP_mean', 'FN_mean', 'IGN_mean', 'OUT_mean'], color=color, line_group=color, 
                  category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
    if axes_mode == 'Fixed' and stat in ['Precision', 'Recall', 'F1_Score']:
        fig.update_layout(yaxis_range=[0, 1.1])
    return fig


# In[20]:


indel_sidebar_plugins = [
    make_strat_selector(idd_df),
    make_type_selector(idd_df),
    make_experiment_selector(idd_df),
    make_axes_mode_selector()
]

min_indel_len = min(idd_df['INDEL_Length'])
max_indel_len = max(idd_df['INDEL_Length'])
indel_len_slider_marks = {
    i: str(i) for i in range(min_indel_len, max_indel_len+1, 5)
}

idd_plot = qbb.PlotPanel(
    header='INDEL Distribution Plot',
    plotter=make_idd_plot,
    plot_inputs={
        'stat': 'F1_Score'
    },
    data_source=idd_df,
    plugins=[
        plg.PlotInputRadioButtons(
            header='Stat to Plot',
            plot_input='stat',
            data_values=['F1_Score', 'Precision', 'Recall', 'IGN', 'OUT']
        ),
        plg.DataFilterRangeSlider(
            header='Range for INDEL Length',
            data_col='INDEL_Length',
            slider_min=min_indel_len,
            slider_max=max_indel_len,
            slider_default_values=[max(-15, min_indel_len), min(15, max_indel_len)],
            slider_step=1,
            slider_marks=indel_len_slider_marks
        )
    ]
)

indel_tab = qbb.BaseTab(
    tab_label='INDEL Plots',
    tab_header='INDEL Plots',
    content_list=[
        idd_plot
    ],
    sidebar_plugins=indel_sidebar_plugins
)


# In[ ]:





# ## Main Board

# In[21]:


board = qbb.Quickboard(
    tab_list=[
        summary_tab,
        roc_tab,
        snp_tab,
        indel_tab
    ]
)


# In[22]:


start_app(board, app_title='BenchmarkBoard', mode='external', port=8055)


# In[ ]:





# In[ ]:




