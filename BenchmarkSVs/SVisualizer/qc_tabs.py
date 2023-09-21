import plotly.express as px
import plotly.graph_objects as go
import quickboard.base as qbb
import quickboard.plugins as plg

from common_utils import read_and_postprocess
from decorators import interval_filter
from plugins import make_interval_plugin_bundle, make_type_selector, make_length_selector, make_filter_selector
from user_config import EXPERIMENT_ORDER, SKIP_COMP_HWE
from qc_data import postprocess_qc_df


QC_DATA_PATH = 'wdl_outputs/combined_qc_stats.tsv'
qc_df = read_and_postprocess(path=QC_DATA_PATH, postprocessor=postprocess_qc_df)

## Counts Tab
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

## Counts Distributions
@interval_filter
def make_qc_histogram(df, x, barmode, interval_name, breakpoint, pct_overlap):
    category_orders = {'Experiment': [f'{x}-{y}' for x in EXPERIMENT_ORDER for y in ['Base', 'Comp']]} if EXPERIMENT_ORDER is not None else None
    title = f'Histogram of {x}'
    if pct_overlap[0] > 0 or pct_overlap[1] < 100:
        title += f' for events with {pct_overlap[0]:.0f}%-{pct_overlap[1]:.0f}% overlap with {interval_name}'
    fig = px.histogram(df, x=x, barmode=barmode, color='Experiment', category_orders=category_orders, title=title)
    return fig

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

## HWE Tab
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

