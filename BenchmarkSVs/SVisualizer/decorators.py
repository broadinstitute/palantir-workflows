from user_config import COVARIATE_X

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

def axes_mode(plotter):
    # Decorator to wrap plotter to use fixed or dynamic axes
    def axes_plotter(df, *args, **kwargs):
        fig = plotter(df, *args, **kwargs)
        axes_mode = kwargs['axes_mode'] if 'axes_mode' in kwargs else None
        if axes_mode == 'Fixed':
            fig.update_layout(yaxis_range=[0, 1])
            if COVARIATE_X is None or len(df[COVARIATE_X].unique()) < 2:
                fig.update_layout(xaxis_range=[0, 1])
        return fig
    return axes_plotter

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

