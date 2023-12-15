import pandas as pd
import plotly.express as px
import quickboard.base as qbb
import quickboard.plugins as plg

from common_utils import read_and_postprocess, add_bbend_stats
from plugins import make_interval_plugin_bundle, make_type_selector, make_filter_selector, make_stat_selector, make_axes_mode_selector
from decorators import axes_mode, filter_or_all_factory, gt_match, interval_filter
from user_config import COVARIATE_X, EXPERIMENT_ORDER


WITTYER_STATS_PATH = 'wdl_outputs/wittyer_stats-cleaned.tsv'
WITTYER_TRUTH_PATH = 'wdl_outputs/wittyer_truth_intervals.tsv'
WITTYER_QUERY_PATH = 'wdl_outputs/wittyer_query_intervals.tsv'
WITTYER_NOGT_PATH = 'wdl_outputs/wittyer_nogt_intervals.tsv'

wittyer_stats_df = read_and_postprocess(path=WITTYER_STATS_PATH, postprocessor=lambda df: df)
wittyer_truth_df = read_and_postprocess(path=WITTYER_TRUTH_PATH, postprocessor=lambda df: df)
wittyer_query_df = read_and_postprocess(path=WITTYER_QUERY_PATH, postprocessor=lambda df: df)
wittyer_nogt_df = read_and_postprocess(path=WITTYER_NOGT_PATH, postprocessor=lambda df: df)

def postprocess_adv_wittyer(df):
    add_bbend_stats(df)
    return df

adv_wittyer_df = pd.concat([wittyer_query_df, wittyer_truth_df])
adv_wittyer_df = postprocess_adv_wittyer(adv_wittyer_df)
wittyer_bins = ['All'] + sorted([str(x) for x in wittyer_stats_df['Bin'].unique() if (x != 'All') & (str(x) != 'nan')]) + ['nan']

## Basic Wittyer Tab
@axes_mode
def make_basic_wittyer_plot(df, axes_mode, stat='Precision'):
    if COVARIATE_X is not None and len(df[COVARIATE_X].unique()) > 1:
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

basic_wittyer_stats_plugins = [
    plg.DataFilterRadioButtons(
        header='Resolution of Stats',
        data_col='StatsType',
        data_values=['Event', 'Base']
    )
]
if COVARIATE_X is not None and len(wittyer_stats_df[COVARIATE_X].unique()) > 1:
    basic_wittyer_stats_plugins += [make_stat_selector(['Precision', 'Recall', 'Fscore'])]

basic_wittyer_plot = qbb.PlotPanel(
    header='Precision vs Recall Plots by Bin',
    plotter=make_basic_wittyer_plot,
    plot_inputs={},
    data_source=wittyer_stats_df,
    plugins=basic_wittyer_stats_plugins
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

## Adv Wittyer Tab
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

adv_wittyer_plugins = [
    plg.PlotInputRadioButtons(
    header='Force Match GT',
    plot_input='gt_match',
    data_values=['False', 'True']
)]

if COVARIATE_X is not None:
    adv_wittyer_plugins += [make_stat_selector(['Precision', 'Recall', 'F1_Score'])]

adv_wittyer_plot = qbb.PlotPanel(
    header='Advanced Wittyer Plot',
    plotter=make_adv_wittyer_plot,
    plot_inputs={},
    data_source=adv_wittyer_df,
    plugins=adv_wittyer_plugins
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

