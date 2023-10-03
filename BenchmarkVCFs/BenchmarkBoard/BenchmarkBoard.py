#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import plotly.express as px

import quickboard.base as qbb
import quickboard.plugins as plg
from quickboard.app import start_app


# In[ ]:





# In[ ]:





# ## User Inputs

# In[ ]:


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


# In[ ]:





# In[ ]:


# DO NOT MODIFY ANY CODE BELOW HERE


# In[ ]:





# ## Load Data

# In[ ]:


# These file paths must point to the corresponding output files from the WDL, either saved locally or bucket links

summary_df = pd.read_csv('SimpleSummary.tsv', sep='\t')
roc_df = pd.read_csv('ROCStats.tsv', sep='\t')
st_df = pd.read_csv('SNPSubstitutionStats.tsv', sep='\t')
idd_df = pd.read_csv('IndelDistributionStats.tsv', sep='\t')


# In[ ]:





# ### Fill Null Stratifier Entries

# In[ ]:


for df in [idd_df, st_df, summary_df]:
    strat_values = df['Stratifier'].unique()
    if 'Whole Genome' not in strat_values:
        df['Stratifier'] = df['Stratifier'].fillna('Whole Genome')
    elif 'Whole Genome (default)' in strat_values:
        raise ValueError("Error: Rename your stratifier labels to not include either 'Whole Genome' or 'Whole Genome (default)'")
    else:
        df['Stratifier'] = df['Stratifier'].fillna('Whole Genome (default)')


# In[ ]:





# ### Fill Experiment Column if not Provided

# In[ ]:


for df in [idd_df, roc_df, st_df, summary_df]:
    if 'Experiment' not in df.columns:
        df['Experiment'] = 'No_ExpGroups_Provided'


# In[ ]:





# ### Other Environment Variables

# In[ ]:


# Check if only one distinguishable sample in summary_df
# Affects behavior of some plots downstream
SINGLE_SAMPLE_MODE = len(summary_df[['Experiment', 'Call_Name', 'Base_Name']].value_counts()) == 1


# In[ ]:


CATEGORY_ORDERS = {
    'Experiment': EXPERIMENT_ORDER
} if EXPERIMENT_ORDER is not None else None


# In[ ]:





# ## Plugin Utilities

# In[ ]:


simple_variants = ['SNP', 'HetSNP', 'HomVarSNP', 'INDEL', 'HetINDEL', 'HomVarINDEL']

def make_strat_selector(df):
    return plg.DataFilterRadioButtons(
        header="Interval List",
        data_col='Stratifier',
        data_values=list(df['Stratifier'].unique())
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
        header="Callset Sample",
        data_col='Call_Name',
        data_values=list(df['Call_Name'].unique())
    )

def make_experiment_selector(df):
    return plg.DataFilterChecklist(
        header="Experimental Groups",
        data_col='Experiment',
        data_values=list(df['Experiment'].unique()),
    )


# In[ ]:





# ## Summary Tab

# In[ ]:


def make_prec_recall_plot(df, marginal, axes_mode):
    strat = df['Stratifier'].iloc[0]
    type_ = df['Type'].iloc[0]
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    if not SINGLE_SAMPLE_MODE:
        marginal = marginal.lower() if marginal != 'None' else None
        fig = px.scatter(df, x='Recall', y='Precision', color=color, marginal_x=marginal, marginal_y=marginal,
                          hover_data=['Call_Name'], title=f'Precision vs Recall Plot over {strat} for {type_}', 
                          category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
        if axes_mode == 'Fixed':
            fig.update_layout(xaxis_range=[0, 1.1], yaxis_range=[0, 1.1])
    else:
        melted_df = df.melt(id_vars=['Experiment', 'Call_Name', 'Base_Name', 'Stratifier', 'Type'], value_vars=['Precision', 'Recall', 'F1_Score'])
        melted_df = melted_df.rename(columns={'variable': 'Stat', 'value': 'Value'})
        fig = px.bar(melted_df, x='Stat', y='Value', title=f'Performance Stats over {strat} for {type_}', 
                     category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
        fig.update_layout(yaxis_range=[0, 1.1])
    
    return fig


def make_stat_covariate_plot(df, covaraite, stat):
    strat = df['Stratifier'].iloc[0]
    type_ = df['Type'].iloc[0]
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    return px.scatter(df, x=stat_corr, y=stat, color=color, hover_data=['Call_Name', 'TP_Base', 'TP_Call', 'FP', 'FN'],
                      title=f'Plot of {stat} by {covariate} over {strat} for {type_}', 
                      category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)


# In[ ]:


summary_sidebar_plugins = [
    make_strat_selector(summary_df),
    make_type_selector(summary_df),
    make_experiment_selector(summary_df)
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
        plg.PlotInputRadioButtons(
            header='Axes Mode',
            plot_input='axes_mode',
            data_values=['Dynamic', 'Fixed']
        )
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

# In[ ]:


def make_roc_plot(df, tp):
    tp_value = f'TP_{tp}'
    color = f'{tp}_Name'
    type_ = df['Type'].iloc[0]
    fig = px.line(df, x='Score', y='Precision', hover_data=['Score', 'TP_Base', 'TP_Call', 'FP', 'FN'], 
                  title=f'ROC Plot on for {type_} by Score', color=color, 
                  category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
    return fig


# In[ ]:


roc_sidebar_plugins = [
    make_type_selector(roc_df),
    make_sample_selector(roc_df)
]

roc_plot = qbb.PlotPanel(
    header="ROC Plot (TP vs FP by Score Field)",
    plotter=make_roc_plot,
    plot_inputs={
        'tp': 'Call'
    },
    data_source=roc_df,
    plugins=[]
)


# In[ ]:


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

# In[ ]:


st_df['Ref_Nucleotide'] = st_df['Substitution'].apply(lambda x: x.split('>')[0])
st_df['Var_Nucleotide'] = st_df['Substitution'].apply(lambda x: x.split('>')[1])


# In[ ]:


def make_snp_substitution_plot(df, stat):
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    strat = df['Stratifier'].iloc[0]
    type_ = df['Type'].iloc[0]
    df_means = df.groupby(['Experiment', 'Type', 'Substitution_Type', 'Ref_Nucleotide', 'Var_Nucleotide'])[['TP_Base', 'TP_Call', 
                                                                               'FP', 'FN', 'F1_Score', 'Precision', 'Recall']].mean().reset_index()
    df_conf = df.groupby(['Experiment', 'Type', 'Substitution_Type', 'Ref_Nucleotide', 'Var_Nucleotide'])[[
        'TP_Base', 'TP_Call', 'FP', 'FN', 'F1_Score', 'Precision', 'Recall']].sem().apply(lambda x: 1.96*x).reset_index()
    
    plot_df = df_means.merge(df_conf, on=['Experiment', 'Ref_Nucleotide', 'Var_Nucleotide'], suffixes=('_mean', '_conf'))
    counts = ['TP_Base_mean', 'TP_Call_mean', 'FP_mean', 'FN_mean']
    plot_df[counts] = plot_df[counts].round(2)
    
    category_orders = {
        'Ref_Nucleotide': ['A', 'G', 'C', 'T'],
        'Var_Nucleotide': ['A', 'G', 'C', 'T']
    }
    if EXPERIMENT_ORDER is not None:
        category_orders = {**category_orders, **{'Experiment': EXPERIMENT_ORDER}}

    fig = px.scatter_3d(plot_df, x='Ref_Nucleotide', y='Var_Nucleotide', z=f'{stat}_mean', error_z=f'{stat}_conf', color=color, 
                        hover_data=['TP_Base_mean', 'TP_Call_mean', 'FP_mean', 'FN_mean'],
                        title=f'Plot of {stat} per Substitution Type on {strat} for {type_}', 
                        category_orders=category_orders, symbol='Substitution_Type_mean', color_discrete_map=EXPERIMENT_COLOR_MAP,
                        height=700, width=1000,
                       )
    fig.update_layout(scene_camera=dict(eye=dict(x=1.5, y=1.5, z=0.4)))
    return fig

def make_titv_plot(df, stat, titv):
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    strat = df['Stratifier'].iloc[0]
    type_ = df['Type'].iloc[0]

    df_counts = df.groupby(['Call_Name', 'Substitution_Type', 'Experiment'])[['TP_Base', 'TP_Call', 'FP', 'FN']].sum().reset_index()
    df_counts['Precision'] = df_counts['TP_Call'] / (df_counts['TP_Call'] + df_counts['FP'])
    df_counts['Recall'] = df_counts['TP_Base'] / (df_counts['TP_Base'] + df_counts['FN'])
    df_counts['F1_Score'] = 2 * df_counts['Precision'] * df_counts['Recall'] / (df_counts['Precision'] + df_counts['Recall'])
    
    fig = px.histogram(df_counts, x=stat, color=color, barmode='overlay', title=f'Plot of {stat} for {titv} over {strat} for {type_}',
                       category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP)
    return fig
    


# In[ ]:


snp_sidebar_plugins = [
    make_strat_selector(st_df),
    make_type_selector(st_df),
    make_experiment_selector(st_df)
]

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
            data_values=['F1_Score', 'Precision', 'Recall', 'TP_Base', 'TP_Call', 'FP', 'FN']
        )
    ]
)

snp_titv_plot = qbb.PlotPanel(
    header='Distribution of Stats for Ti or Tv',
    plotter=make_titv_plot,
    plot_inputs={
        'stat': 'F1_Score',
        'titv': 'Ti',
    },
    data_source=st_df,
    plugins=[
        plg.PlotInputRadioButtons(
            header='Stat to Plot',
            plot_input='stat',
            data_values=['F1_Score', 'Precision', 'Recall', 'TP_Base', 'TP_Call', 'FP', 'FN']
        ),
        plg.PlotInputRadioButtons(
            header='Substitution Type',
            plot_input='titv',
            data_values=['Ti', 'Tv']
        )
    ]
)

snp_cg = qbb.ContentGrid(
    content_list=[
        snp_substitution_plot,
        snp_titv_plot
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





# In[ ]:





# ## INDEL Tab

# In[ ]:


def make_idd_plot(df, stat):
    color = None if df['Experiment'].iloc[0] == 'No_ExpGroups_Provided' else 'Experiment'
    strat = df['Stratifier'].iloc[0]
    type_ = df['Type'].iloc[0]

    if not SINGLE_SAMPLE_MODE:
        df_means = df.groupby(['Experiment', 'INDEL_Length'])[['TP_Base', 'TP_Call', 'FP', 'FN', 'F1_Score', 'Precision', 'Recall']].mean().reset_index()
        # Use naive noise model for error bars
        df_conf = df.groupby(['Experiment', 'INDEL_Length'])[['TP_Base', 'TP_Call', 'FP', 'FN', 
                                                              'F1_Score', 'Precision', 'Recall']].sem().apply(lambda x: 1.96*x).reset_index()
        
        df_means = df_means.round(2)
        df_conf = df_conf.round(4)

        error_y = f'{stat}_conf'
    else:
        error_y = None
    
    plot_df = df_means.merge(df_conf, on=['Experiment', 'INDEL_Length'], suffixes=('_mean', '_conf'))
    fig = px.bar(
        plot_df, x='INDEL_Length', y=f'{stat}_mean', error_y=error_y, title=f'Plot of {stat} mean by INDEL Length on {strat} for {type_}',
        hover_data=['TP_Base_mean', 'TP_Call_mean', 'FP_mean', 'FN_mean'], color=color, barmode='group',
        category_orders=CATEGORY_ORDERS, color_discrete_map=EXPERIMENT_COLOR_MAP
    )
    return fig


# In[ ]:


indel_sidebar_plugins = [
    make_strat_selector(idd_df),
    make_type_selector(idd_df),
    make_experiment_selector(idd_df)
]

min_indel_len = min(idd_df['INDEL_Length'])
max_indel_len = max(idd_df['INDEL_Length'])
indel_len_slider_marks = {
    i: str(i) for i in range(min_indel_len, max_indel_len+1, 5)
}

idd_plot = qbb.PlotPanel(
    header='InDel Distribution Plot',
    plotter=make_idd_plot,
    plot_inputs={
        'stat': 'F1_Score'
    },
    data_source=idd_df,
    plugins=[
        plg.PlotInputRadioButtons(
            header='Stat to Plot',
            plot_input='stat',
            data_values=['F1_Score', 'Precision', 'Recall']
        ),
        plg.DataFilterRangeSlider(
            header='Range for INDEL Length',
            data_col='INDEL_Length',
            slider_min=min_indel_len,
            slider_max=max_indel_len,
            slider_default_values=[max(-20, min_indel_len), min(20, max_indel_len)],
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

# In[ ]:


board = qbb.Quickboard(
    tab_list=[
        summary_tab,
        roc_tab,
        snp_tab,
        indel_tab
    ]
)


# In[ ]:


start_app(board, app_title='BenchmarkBoard', mode='external', port=8050)

