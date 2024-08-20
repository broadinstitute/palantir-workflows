version 1.0

task F1Evaluation {
    input {
        Array[File] roc_tables
        String tool1_label
        String tool2_label
        String? additional_label
        Boolean signed_difference = false
        Int plot_qual_limit = 30
        Int? mem_gb
        Int? preemptible
    }
    
    Int machine_mem_gb = select_first([mem_gb, 8])

    command <<<
        set -xeuo pipefail
        
        python3 << 'EOF'
        import matplotlib
        import matplotlib.pyplot as plt
        import pandas as pd
        import numpy as np
        import itertools

        # Set plot parameters
        matplotlib.rcParams['text.usetex'] = False
        matplotlib.rcParams['mathtext.default'] = 'regular'
        matplotlib.rcParams['font.family'] = 'serif'

        # Setup user parameters
        VARIANT_TYPES = ["SNP", "INDEL"]
        PLOT_QUAL_LIMIT = ~{plot_qual_limit}
        SIGNED_DIFFERENCE = ~{true="True" false="False" signed_difference}
        tool1_label = "~{tool1_label}"
        tool2_label = "~{tool2_label}"
        additional_label = "~{additional_label}"

        # Import data
        roc_df = pd.concat([pd.read_csv(roc, sep='\t') for roc in ["~{sep="\", \"" roc_tables}"]])

        roc_df = roc_df.rename(columns={
            '#score': 'Score',
            'true_positives_baseline': 'TP_Base',
            'false_positives': 'FP',
            'true_positives_call': 'TP_Query',
            'false_negatives': 'FN',
            'precision': 'Precision',
            'sensitivity': 'Recall',
            'f_measure': 'F1_Score'
        })

        roc_df['F1_Score'] = roc_df['TP_Base'] / (roc_df['TP_Base'] + 0.5 * (roc_df['FP'] + roc_df['FN']))

        ## Set some plotting and data utility functions
        def get_discrete_stats(stats_df, X):
            # For every x in X range, find closest value in stats_df
            indices = []
            stats_df = stats_df.sort_values(by='Score', ascending=True)
            sorted_scores = stats_df['Score']
            new_df = pd.DataFrame()
            rows = []
            for x in X:
                index = sorted_scores.searchsorted(x, side='right')
                if index < len(stats_df):
                    rows += [stats_df.iloc[[index]]]

            new_df = pd.concat(rows)
            new_df['Score'] = X
            return new_df

        def get_intra_difference_data(roc_df, exp):
            # Works on specific Experiment subset to take all cross-differences
            sub_df = roc_df[roc_df['Experiment'] == exp]
            replicate_dict = {
                i: sub_df.loc[sub_df['Replicate'] == i, :] for i in sorted(sub_df['Replicate'].unique())
            }

            full_diff_df = pd.DataFrame()

            # For each pair of distinct replicates, take differences of F1 scores
            # for c in itertools.combinations(replicate_dict.keys(), 2):    # This is probably "better" but not backwards compatible
            for c in itertools.product(replicate_dict.keys(), replicate_dict.keys()):
                if c[0] == c[1]:
                    continue
                df1 = replicate_dict[c[0]]
                df2 = replicate_dict[c[1]]

                # Discretize score before taking diffs
                X = np.arange(0, PLOT_QUAL_LIMIT+1)
                df1 = df1.groupby(['Type', 'Interval', 'Dataset', 'Replicate']).apply(lambda df: get_discrete_stats(df, X), include_groups=False).reset_index()
                df2 = df2.groupby(['Type', 'Interval', 'Dataset', 'Replicate']).apply(lambda df: get_discrete_stats(df, X), include_groups=False).reset_index()

                combined_df = df1.merge(df2, on=['Score', 'Type', 'Interval', 'Dataset'], how='outer').sort_values(by='Score')
                combined_df['F1_Score_diff'] = combined_df['F1_Score_x'].interpolate().ffill().bfill() - combined_df['F1_Score_y'].interpolate().ffill().bfill()

                if not SIGNED_DIFFERENCE:
                    combined_df['F1_Score_diff'] = combined_df['F1_Score_diff'].apply(np.abs)

                full_diff_df = pd.concat([full_diff_df, combined_df])

            return full_diff_df

        def get_inter_difference_data(roc_df, exp1="EvalVsTruthTool1", exp2="EvalVsTruthTool2"):
            df1 = roc_df[roc_df['Experiment'] == exp1]
            df2 = roc_df[roc_df['Experiment'] == exp2]

            # Discretize score before taking diffs
            X = np.arange(0, PLOT_QUAL_LIMIT+1)
            df1 = df1.groupby(['Type', 'Interval', 'Dataset', 'Replicate']).apply(lambda df: get_discrete_stats(df, X), include_groups=False).reset_index()
            df2 = df2.groupby(['Type', 'Interval', 'Dataset', 'Replicate']).apply(lambda df: get_discrete_stats(df, X), include_groups=False).reset_index()

            combined_df = df1.merge(df2, on=['Score', 'Type', 'Interval', 'Dataset', 'Replicate'], how='outer').sort_values(by='Score')
            combined_df['F1_Score_diff'] = combined_df['F1_Score_x'].interpolate().ffill().bfill() - combined_df['F1_Score_y'].interpolate().ffill().bfill()

            if not SIGNED_DIFFERENCE:
                combined_df['F1_Score_diff'] = combined_df['F1_Score_diff'].apply(np.abs)

            return combined_df

        intra_diffs1 = get_intra_difference_data(roc_df, exp='EvalVsTruthTool1')
        intra_diffs2 = get_intra_difference_data(roc_df, exp='EvalVsTruthTool2')
        inter_diff_df = get_inter_difference_data(roc_df, 'EvalVsTruthTool1', 'EvalVsTruthTool2')

        def get_stats(diff_df):
            return diff_df.groupby(['Score', 'Type', 'Interval', 'Dataset'])['F1_Score_diff'].describe().reset_index()

        intra_stats1 = get_stats(intra_diffs1)
        intra_stats2 = get_stats(intra_diffs2)
        inter_stats = get_stats(inter_diff_df)

        def make_single_line_plot(ax, stats_df, min_score, var_type, interval, color):
            if SIGNED_DIFFERENCE:
                ax.set_ylabel(r'$\Delta\;F_1$')
                ax.axhline(y=0, c='grey', linewidth=1)    # If using signed_difference, draw a line to indicate 0 on the vertical axis
            else:
                ax.set_ylabel(r'$|\Delta\;F_1|$')

            sub_df = stats_df[(stats_df['Type'] == var_type) & (stats_df['Interval'] == interval)]
            X = np.arange(0, PLOT_QUAL_LIMIT+1)

            if min_score > PLOT_QUAL_LIMIT:
                min_score_f1_mean = sub_df.loc[sub_df['Score'] == PLOT_QUAL_LIMIT, 'mean'].values[0]
                min_score_f1_std = sub_df.loc[sub_df['Score'] == PLOT_QUAL_LIMIT, 'std'].values[0]
                Y_mean = np.array([min_score_f1_mean] * (PLOT_QUAL_LIMIT+1))
                Y_std = np.array([min_score_f1_std] * (PLOT_QUAL_LIMIT+1))
                linestyle = 'dotted'
            else:
                sub_df = sub_df[sub_df['Score'] <= PLOT_QUAL_LIMIT]    # Apply user input bound
                Y_mean = sub_df['mean'].reset_index(drop=True)
                Y_std = sub_df['std'].reset_index(drop=True)
                linestyle = '-'

            ax.plot(X, Y_mean, c=color, linestyle=linestyle)
            ax.fill_between(X, Y_mean - Y_std, Y_mean + Y_std, color=color, alpha=0.1)
            return Y_mean, Y_std


        min_df = roc_df.groupby(['Experiment', 'Type', 'Interval', 'Dataset'])['Score'].min().reset_index()    # Get min score values for dotted lines

        def draw_line_plot(ax, inter_stats, intra_stats1, intra_stats2, var_type, interval, dataset):
            intra1_min_score = min_df[
                (min_df['Experiment'] == 'EvalVsTruthTool1') & (min_df['Type'] == var_type) & (min_df['Interval'] == interval) & (min_df['Dataset'] == dataset)
            ]['Score'].values[0]
            intra1_diff_means, intra1_diff_std = make_single_line_plot(ax, intra_stats1, intra1_min_score, var_type, interval, color='C0')
            intra2_min_score = min_df[
                (min_df['Experiment'] == 'EvalVsTruthTool2') & (min_df['Type'] == var_type) & (min_df['Interval'] == interval) & (min_df['Dataset'] == dataset)
            ]['Score'].values[0]
            intra2_diff_means, intra2_diff_std = make_single_line_plot(ax, intra_stats2, intra2_min_score, var_type, interval, color='C1')
            inter_min_score = max(intra1_min_score, intra2_min_score)
            inter_diff_means, inter_diff_std = make_single_line_plot(ax, inter_stats, inter_min_score, var_type, interval, color='C2')

            # Check the FE criterion
            fe_status = 0
            if np.any(inter_diff_means > intra1_diff_means) or np.any(inter_diff_means > intra2_diff_means):
                titlecolor = 'orange'
                fe_status = 2
            else:
                titlecolor = 'white'
                fe_status = 0

            ax.set_xlabel('Quality threshold q')
            ax.set_xlim(0, PLOT_QUAL_LIMIT)
            y_min = min(min(inter_diff_means - inter_diff_std), min(intra1_diff_means - intra1_diff_std), min(intra2_diff_means - intra2_diff_std), 0)
            y_max = max(max(inter_diff_means + inter_diff_std), max(intra1_diff_means + intra1_diff_std), max(intra2_diff_means + intra2_diff_std))
            ax.set_ylim(y_min, y_max * 1.1)

            ax.set_title(f'{var_type} {interval}', backgroundcolor=titlecolor, zorder=0)

            return fe_status

        def make_dataset_plot(roc_df, dataset):
            fe_status = 0
            stratifiers = inter_stats['Interval'].unique()
            stratifiers = ['WholeGenome'] + sorted([x for x in stratifiers if x != 'WholeGenome'])

            num_columns = max(len(stratifiers), 3)
            num_rows = len(VARIANT_TYPES)+1

            fig, axes = plt.subplots(num_rows, num_columns, figsize=(3*num_columns, 8))
            for row, var_type in enumerate(VARIANT_TYPES):
                for col, interval in enumerate(stratifiers):
                    column_to_plot = col if len(stratifiers) > 1 else 1
                    ax = axes[row, column_to_plot]
                    fe_status = max(fe_status, draw_line_plot(ax, inter_stats, intra_stats1, intra_stats2, var_type, interval, dataset))

            if SIGNED_DIFFERENCE:
                inter_label = 'Inter: $F_{1, ' + tool1_label + '\\ rep 1} - F_{1, ' + tool2_label + '\\ rep 1}$'
            else:
                inter_label = 'Inter: $|F_{1, ' + tool1_label + '\\ rep 1} - F_{1, ' + tool2_label + '\\ rep 1}|$'

            legend_inter_line = matplotlib.lines.Line2D([], [], color='C2', label=inter_label)
            legend_intra_tool1_line = matplotlib.lines.Line2D([], [], color='C0', label='Intra: $|F_{1, ' + tool1_label + '\\ rep 1} - F_{1, ' + tool1_label + '\\ rep 2}|$')
            legend_intra_tool2_line = matplotlib.lines.Line2D([], [], color='C1', label='Intra: $|F_{1, ' + tool2_label + '\\ rep 1} - F_{1, ' + tool2_label + '\\ rep 2}|$')
            fig.legend(bbox_to_anchor=(0.5, 0.2), loc='center', handles=[legend_inter_line, legend_intra_tool1_line, legend_intra_tool2_line])

            # Clear non-used axes if plotting less than 3 columns
            for c in range(num_columns):
                axes[num_rows-1, c].axis('off')

            fig.suptitle(f'Dataset: {dataset}' + ('' if not additional_label else f', {additional_label}') + '\n' +
                        ('Signed ' if SIGNED_DIFFERENCE else 'Absolute ') +
                        r'$F_1 = \frac{TP}{TP + \frac{1}{2} (FP + FN)}$ score differences for calls with $QUAL \geq q$' +
                        f' (# replicates: {len(roc_df.loc[roc_df["Dataset"] == dataset]["Replicate"].unique())})')

            fig.tight_layout()
            fig.savefig(f'f1_plot_{dataset}.png', dpi=100)
            return fe_status

        fe_status = 0
        for dataset in roc_df['Dataset'].unique():
            fe_status = max(fe_status, make_dataset_plot(roc_df, dataset))

        # Collect stats for final output
        remove_stats = ['count', 'min', '25%', '50%', '75%', 'max']
        data_df = pd.DataFrame()
        data_df = intra_stats1.drop(columns=remove_stats).merge(intra_stats2.drop(columns=remove_stats), on=['Score', 'Type', 'Interval', 'Dataset']).rename(columns={'mean_x': 'intra1_mean', 'std_x': 'intra1_std', 'mean_y': 'intra2_mean', 'std_y': 'intra2_std'})
        data_df = data_df.merge(inter_stats, on=['Score', 'Type', 'Interval', 'Dataset']).drop(columns=remove_stats).rename(columns={'mean': 'inter_mean', 'std': 'inter_std'})
        data_df.to_csv('f1_summary.tsv', sep='\t', index=False)

        # Write final fe_status
        with open('fe_status.txt', 'w') as file:
            file.write(f'{fe_status}')

        EOF
    >>>

    output {
        Array[File] f1_plots = glob("*.png")
        File f1_summary = "f1_summary.tsv"
        Int fe_status = read_int("fe_status.txt")
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim-plots:1.1"
        preemptible: select_first([preemptible, 0])
        memory: machine_mem_gb + " GB"
        disks: "local-disk 20 HDD"
    }
}