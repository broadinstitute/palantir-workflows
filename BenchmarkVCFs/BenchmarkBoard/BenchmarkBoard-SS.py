import pandas as pd
import numpy as np
import plotly.express as px

import quickboard.base as qbb
import quickboard.plugins as plg
from quickboard.app import start_app


## Load Data

# These file paths must point to the corresponding output files from the WDL, either saved locally or bucket links

idd_df = pd.read_csv('Full_IDD.tsv', sep='\t')
roc_df = pd.read_csv('Full_ROC.tsv', sep='\t')
st_df = pd.read_csv('Full_ST.tsv', sep='\t')
summary_df = pd.read_csv('SimpleSummary.tsv', sep='\t')


### Fill Null Stratifier Entries

for df in [idd_df, roc_df, st_df, summary_df]:
    strat_values = df['Stratifier'].unique()
    if 'Whole Genome' not in strat_values:
        df['Stratifier'] = df['Stratifier'].fillna('Whole Genome')
    elif 'Whole Genome (default)' in strat_values:
        raise ValueError("Error: Rename your stratifier labels to not include either 'Whole Genome' or 'Whole Genome (default)'")
    else:
        df['Stratifier'] = df['Stratifier'].fillana('Whole Genome (default)')

## Plugin Utilities

def make_strat_selector(df):
    return plg.DataFilterRadioButtons(
        header="Interval List",
        data_col='Stratifier',
        data_values=list(df['Stratifier'].unique())
    )

def make_type_selector(df):
    return plg.DataFilterRadioButtons(
        header="Variant Type",
        data_col='Type',
        data_values=list(df['Type'].unique())
    )

def make_sample_selector(df):
    return plg.DataFilterRadioButtons(
        header="Callset Sample",
        data_col='Call_Name',
        data_values=list(df['Call_Name'].unique())
    )

## Summary Tab

summary_sidebar_plugins = [
    plg.DataFilterChecklist(
        header="Interval List",
        data_col='Stratifier',
        data_values=list(df['Stratifier'].unique())
    ),
    plg.DataFilterChecklist(
        header="Variant Type",
        data_col='Type',
        data_values=list(df['Type'].unique())
    )
]

summary_body_text = f"""This table contains counts of benchmarking statistics (FP, Precision, etc.) for each 
sample across various interval stratifiers and variant subtypes. 

You can filter the table below by typing into the fields above each column. 

E.g. type 'SNP' in the Type column to see stats on variant types containing the string 'SNP' only, or '=SNP' for just those with value exactly 'SNP'.

You can also sort values in a column alphabetically/ascending using the arrows on the left of the column title. Use syntax like '>0.75' to set 
thresholds to filter by.
"""

# Round entries to shorten table column width
summary_df = summary_df.round(4)

summary_table = qbb.DataPanel(
    header='Summary Table',
    body=summary_body_text,
    data_source=summary_df,
)

summary_tab = qbb.BaseTab(
    tab_label='Summary Stats',
    tab_header='Summary Statistics',
    content_list=[
        summary_table
    ],
    sidebar_plugins=summary_sidebar_plugins
)

## ROC Tab

def make_roc_plot(df, score_field, tp):
    tp_value = f'TP_{tp}'
    color = f'{tp}_Name'
    strat = df['Stratifier'].iloc[0]
    type_ = df['Type'].iloc[0]
    fig = px.line(df, x='Score', y='Precision', hover_data=['Score', 'Score_Field', 'TP_Base', 'TP_Call', 'FP', 'FN'], 
                  title=f'ROC Plot on {strat} for {type_} by {score_field}', color=color)
    return fig

roc_sidebar_plugins = [
    make_strat_selector(roc_df),
    make_type_selector(roc_df),
    plg.DataFilterRadioButtons(
        header="Score Field",
        data_col='Score_Field',
        data_values=list(roc_df['Score_Field'].unique())
    ),
    make_sample_selector(roc_df)
]

roc_plot = qbb.PlotPanel(
    header="ROC Plot (TP vs FP by Score Field)",
    plotter=make_roc_plot,
    plot_inputs={
        'score_field': roc_df['Score_Field'].iloc[0],
        'tp': 'Call'
    },
    data_source=roc_df,
    plugins=[]
)

roc_tab = qbb.BaseTab(
    tab_label='ROC Plots',
    tab_header='ROC Plots',
    content_list=[
        roc_plot
    ],
    sidebar_plugins=roc_sidebar_plugins
)

## SNP Tab

st_df['Ref_Nucleotide'] = st_df['Substitution'].apply(lambda x: x.split('>')[0])
st_df['Var_Nucleotide'] = st_df['Substitution'].apply(lambda x: x.split('>')[1])

def make_snp_plot(df, stat):
    strat = df['Stratifier'].iloc[0]
    type_ = df['Type'].iloc[0]
    category_orders = {
        'Ref_Nucleotide': ['A', 'G', 'C', 'T'],
        'Var_Nucleotide': ['A', 'G', 'C', 'T']
    }
    fig = px.scatter_3d(df, x='Ref_Nucleotide', y='Var_Nucleotide', z=stat, color='Substitution_Type', 
                        hover_data=['TP_Base', 'TP_Call', 'FP', 'FN'],
                        title=f'Plot of {stat} per Substitution Type on {strat} for {type_}', 
                        category_orders=category_orders,
                        height=700, width=1000,
                       )
    fig.update_layout(scene_camera=dict(eye=dict(x=1.5, y=1.5, z=0.4)))
    return fig

snp_sidebar_plugins = [
    make_strat_selector(st_df),
    make_type_selector(st_df),
    make_sample_selector(st_df)
]

snp_plot = qbb.PlotPanel(
    header='SNP Substitution Plot',
    plotter=make_snp_plot,
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

snp_tab = qbb.BaseTab(
    tab_label='SNP Plots',
    tab_header='SNP Plots',
    content_list=[
        snp_plot
    ],
    sidebar_plugins=snp_sidebar_plugins
)

## INDEL Tab

def make_idd_plot(df, stat):
    strat = df['Stratifier'].iloc[0]
    type_ = df['Type'].iloc[0]
    fig = px.bar(df, x='INDEL_Length', y=stat, title=f'Plot of {stat} by INDEL Length on {strat} for {type_}',
                hover_data=['TP_Base', 'TP_Call', 'FP', 'FN'])
    return fig

indel_sidebar_plugins = [
    make_strat_selector(idd_df),
    make_type_selector(idd_df),
    make_sample_selector(idd_df)
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

## Main Board

board = qbb.Quickboard(
    tab_list=[
        summary_tab,
        roc_tab,
        snp_tab,
        indel_tab
    ]
)

start_app(board, app_title='BenchmarkBoard', mode='external', port=8050)
