import quickboard.plugins as plg

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
        data_values=list(df['SVTYPE'].sort_values().unique())
    )

def make_length_selector(df):
    return plg.DataFilterChecklist(
        header="SV Lengths",
        data_col='SVLEN_Bin',
        data_values=list(df['SVLEN_Bin'].sort_values().unique())
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