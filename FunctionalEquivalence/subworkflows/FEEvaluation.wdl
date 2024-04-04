version 1.0

task FEEvaluation {
    input {
        Array[File] benchmark_summaries
        String tool1_label
        String tool2_label
        String? additional_label
    }

    String title_label = if (defined(additional_label)) then ", " + additional_label else ""

    command <<<
        set -xueo pipefail

        python3 << CODE
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import matplotlib

        ## A few settings for the rest of the script
        VARIANT_TYPES = ['SNP', 'INDEL']
        matplotlib.rcParams['text.usetex'] = False
        matplotlib.rcParams['mathtext.default'] = 'regular'
        matplotlib.rcParams['font.family'] = 'serif'

        full_df = pd.DataFrame()
        for file in ["~{sep="\", \"" benchmark_summaries}"]:
            df = pd.read_csv(file, sep="\t")
            full_df = pd.concat([full_df, df])

        ## Collect interval names and sort so WholeGenome is always first
        intervals = full_df['Interval'].unique()
        intervals = ['WholeGenome'] + sorted([i for i in intervals if i != 'WholeGenome'])
        full_df['Interval'] = pd.Categorical(full_df['Interval'], intervals)
        full_df = full_df.sort_values('Interval')

        ## Subset/augment the full_df a bit
        full_df = full_df[full_df['Type'].isin(VARIANT_TYPES)]
        full_df['Jaccard'] = full_df['TP_Base'] / (full_df['TP_Base'] + full_df['FP'] + full_df['FN'])

        ## Create the plot_df from stats of Jaccard
        mean_df = full_df.groupby(['Interval', 'Type', 'Experiment', 'Dataset'])['Jaccard'].mean().reset_index().rename(columns={'Jaccard': 'Jaccard_mean'})
        print(f'Means of Jaccard are: {mean_df}')
        sem_df = full_df.groupby(['Interval', 'Type', 'Experiment', 'Dataset'])['Jaccard'].sem().reset_index().rename(columns={'Jaccard': 'Jaccard_sem'})
        print(f'Stderr of Jaccard are: {sem_df}')
        plot_df = mean_df.merge(sem_df, on=['Interval', 'Type', 'Experiment', 'Dataset'])

        ## Function to fill in FE plot for given parameters and return fe_status indicator value
        def make_grid_plot(ax, plot_df, var_type, stratifier):
            fe_status = 0
            x = np.array([0, 1, 2])
            xticklabels = ["~{tool1_label}", 'Inter', "~{tool2_label}"]

            sub_df = plot_df[(plot_df['Type'] == var_type) & (plot_df['Interval'] == stratifier)]
            sub_df['Experiment'] = pd.Categorical(sub_df['Experiment'], ['EvalIntraTool1', 'EvalInterTool', 'EvalIntraTool2'])
            sub_df = sub_df.sort_values(by='Experiment')
            print("sub_df Experiment sort DEBUG!")
            print(sub_df)

            y = sub_df['Jaccard_mean']
            y_err = sub_df['Jaccard_sem']

            ax.set_xticks(x)
            ax.set_xticklabels(xticklabels)
            ax.set_ylabel('Jaccard score')
            ax.errorbar(x, y, yerr=y_err, c='k', fmt='o', capsize=15)
            ax.set_xlim(-0.5, 2.5)
            ax.grid(axis='y')

            inter_mean = sub_df[sub_df['Experiment'] == "EvalInterTool"]['Jaccard_mean'].values[0]
            inter_sem = sub_df[sub_df['Experiment'] == "EvalInterTool"]['Jaccard_sem'].values[0]
            tool1_mean = sub_df[sub_df['Experiment'] == "EvalIntraTool1"]['Jaccard_mean'].values[0]
            tool1_sem = sub_df[sub_df['Experiment'] == "EvalIntraTool1"]['Jaccard_sem'].values[0]
            tool2_mean = sub_df[sub_df['Experiment'] == "EvalIntraTool2"]['Jaccard_mean'].values[0]
            tool2_sem = sub_df[sub_df['Experiment'] == "EvalIntraTool2"]['Jaccard_sem'].values[0]

            if inter_mean <= max(tool1_mean, tool2_mean):
                titlecolor = 'orange'
                fe_status = 2
            elif inter_mean - np.nan_to_num(inter_sem) < max(tool1_mean + np.nan_to_num(tool1_sem), tool2_mean + np.nan_to_num(tool2_sem)):
                titlecolor = 'yellow'
                fe_status = 1
            else:
                titlecolor = 'white'
            ax.set_title(f'{var_type} {stratifier}', backgroundcolor=titlecolor, zorder=0)
            return fe_status


        ## Make plots
        for dataset in plot_df['Dataset'].unique():
            num_columns = max(3, len(intervals))   # Why cap at 3?
            fig, axes = plt.subplots(len(VARIANT_TYPES), num_columns, figsize=(3*num_columns , 6))

            fe_status = 0
            dataset_plot_df = plot_df[plot_df['Dataset'] == dataset]
            for row, var_type in enumerate(VARIANT_TYPES):
                for col, stratifier in enumerate(intervals[:3]):
                    fe_status = max(fe_status, make_grid_plot(axes[row, col], dataset_plot_df, var_type, stratifier))

            ## Collect number of replicates per experiment
            sub_df = full_df[(full_df['Type'] == 'SNP') & (full_df['Interval'] == 'WholeGenome') & (full_df['Dataset'] == dataset)]
            tool1_conc_count = len(sub_df[sub_df['Experiment'] == 'EvalIntraTool1'])
            tool2_conc_count = len(sub_df[sub_df['Experiment'] == 'EvalIntraTool2'])
            inter_conc_count = len(sub_df[sub_df['Experiment'] == 'EvalInterTool'])

            concordance_text = f'Concordance (# values ~{tool1_label}: {tool1_conc_count} / ~{tool2_label}: {tool2_conc_count} / inter: {inter_conc_count})'
            title = f'Dataset: {dataset}~{title_label}\n{concordance_text}'
            fig.suptitle(title)
            fig.tight_layout()
            fig.savefig(f'fe_plot_{dataset}.png', dpi=100)

        with open('fe_status.txt', 'w') as file:
            file.write(f'{fe_status}')

        full_df.to_csv('fe_summary.tsv', sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim-plots:1.0"
        preemptible: select_first([2, 0])
        memory: 8 + " GB"
        disks: "local-disk 20 HDD"
    }

    output {
        Array[File] fe_plots = glob("fe_plot_*.png")
        File fe_summary = "fe_summary.tsv"
        Int fe_status = read_int("fe_status.txt")
    }
}