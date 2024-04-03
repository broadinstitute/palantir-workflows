version 1.0

workflow PlotROC {
    input {
        String sample_id
        Array[File] roc_tables
        String tool1_label
        String tool2_label
        String? additional_label
        Int? mem_gb
        Int? preemptible
    }

    call PlotROCTask {
        input:
            sample_id = sample_id,
            roc_tables = roc_tables,
            tool1_label = tool1_label,
            tool2_label = tool2_label,
            additional_label = additional_label,
            mem_gb = mem_gb,
            preemptible = preemptible
    }

    output {
        Array[File] plots = PlotROCTask.plots
    }
}

task PlotROCTask {
    input {
        String sample_id
        Array[File] roc_tables
        String tool1_label
        String tool2_label
        String? additional_label
        Int? mem_gb
        Int? preemptible
    }
    
    Int machine_mem_gb = select_first([mem_gb, 8])

    String additional_label_arg = if defined(additional_label) then "--additional-label \"" + additional_label + "\"" else ""

    command <<<
        set -xeuo pipefail

        cat <<'EOF' > script.py
        import pandas as pd
        import matplotlib
        import matplotlib.pyplot as plt
        import argparse
        import os

        matplotlib.rcParams['text.usetex'] = False
        matplotlib.rcParams['mathtext.default'] = 'regular'
        matplotlib.rcParams['font.family'] = 'serif'

        VARIANT_TYPES = ["SNP", "INDEL"]

        def parse_roc_to_dicts(df):
            data = {}
            max_f1 = 0
            best_qual = None

            # Transform these columns into a list of dicts over the rows with keys column names
            for d in df[['Score', 'FP', 'Recall', 'F1_Score']].to_dict(orient='records'):
                data[d['Score']] = (d['FP'], d['Recall'])
                if d['F1_Score'] > max_f1:
                    best_qual = (d['Score'], d['F1_Score'])

            return data, best_qual


        def read_file(data, best_qual, filename):
            df = pd.read_csv(filename, sep='\t')

            for region in df['Interval-test'].unique():
                for var_type in df['Type'].unique():
                    for tool in df['Tool'].unique():
                        sub_df = df[(df['Interval-test'] == region) & (df['Type'] == var_type) & (df['Tool'] == tool)]
                        data[(region, var_type, tool)], best_qual[(region, var_type, tool)] = parse_roc_to_dicts(sub_df)

        def plot_roc(ax, data, best_qual, region, var_type):
            X1 = [x[0] for x in data[(region, var_type, 'tool1')].values()]
            Y1 = [x[1] for x in data[(region, var_type, 'tool1')].values()]
            X2 = [x[0] for x in data[(region, var_type, 'tool2')].values()]
            Y2 = [x[1] for x in data[(region, var_type, 'tool2')].values()]

            ax.plot(X1, Y1, c='C0', marker='.', markersize=2, linestyle=' ')
            ax.plot(X2, Y2, c='C1', marker='.', markersize=2, linestyle=' ')

            best_qual_tool1 = best_qual[(region, var_type, 'tool1')]
            best_qual_tool2 = best_qual[(region, var_type, 'tool2')]

            best_coordinates_tool1 = data[(region, var_type, 'tool1')][best_qual_tool1[0]]
            best_coordinates_tool2 = data[(region, var_type, 'tool2')][best_qual_tool2[0]]

            ax.plot([best_coordinates_tool1[0]], [best_coordinates_tool1[1]], c='C0', marker='x', markersize=15, linestyle=' ')
            ax.plot([best_coordinates_tool2[0]], [best_coordinates_tool2[1]], c='C1', marker='x', markersize=15, linestyle=' ')

            ax.annotate(r'$F_{1, max}$ = '+'{:.3f} @ Q{:.0f}'.format(best_qual_tool1[1], best_qual_tool1[0]), (0.95, 0.05), xycoords='axes fraction', xytext=(0, 12), textcoords='offset points', c='C0', ha='right', va='bottom')
            ax.annotate(r'$F_{1, max}$ = '+'{:.3f} @ Q{:.0f}'.format(best_qual_tool2[1], best_qual_tool2[0]), (0.95, 0.05), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', c='C1', ha='right', va='bottom')

            ax.set_xlabel(r'FP')
            ax.set_ylabel(r'Sensitivity')

            ax.set_title(r'{} {}'.format(var_type, region), zorder=0)

        def plot_data(data, best_qual, stratifiers, sample_id, tool1_label, tool2_label, additional_label):
            num_columns = max(len(stratifiers), 3)
            fig, axes = plt.subplots(3, num_columns, figsize=(6*num_columns,16))
            for row, var_type in enumerate(VARIANT_TYPES):
                for col, region in enumerate(stratifiers):
                    column_to_plot = col if len(stratifiers) > 1 else 1
                    ax = axes[row, column_to_plot]
                    plot_roc(ax, data, best_qual, region, var_type)

            # Clear axes for legend
            for i in range(num_columns):
                axes[2, i].axis('off')

            # Clear non-used axes if plotting less than 3 columns
            if len(stratifiers) == 1:
                axes[0, 0].axis('off')
                axes[0, 2].axis('off')
                axes[1, 0].axis('off')
                axes[1, 2].axis('off')
            if len(stratifiers) == 2:
                axes[0, 2].axis('off')
                axes[1, 2].axis('off')


            fig.suptitle('Sample: {}'.format(sample_id) + ('' if not additional_label else ', {}'.format(additional_label)) + '\n')

            legend_tool1_line = matplotlib.lines.Line2D([], [], color='C0', label=tool1_label)
            legend_tool2_line = matplotlib.lines.Line2D([], [], color='C1', label=tool2_label)
            fig.legend(bbox_to_anchor=(0.5, 0.2), loc='center', handles=[legend_tool1_line, legend_tool2_line])
            plt.tight_layout()
            fig.savefig('roc_plot_{}.png'.format(sample_id), dpi=100)

        def main(roc_tables, sample_id, tool1_label, tool2_label, additional_label):
            stratifiers = list(pd.read_csv(roc_tables[0], sep='\t')['Interval-test'].unique())    # Grab list of stratifiers used from first file

            data = dict()
            best_qual = dict()
            for filename in roc_tables:
                read_file(data, best_qual, filename)

            plot_data(data, best_qual, stratifiers, sample_id, tool1_label, tool2_label, additional_label)

        if __name__ == '__main__':
            parser = argparse.ArgumentParser(description='Create F1 functional equivalence plots.')
            parser.add_argument('--additional-label', type=str)
            required_named = parser.add_argument_group('Required named arguments')
            required_named.add_argument('--sample-id', type=str)
            required_named.add_argument('--tool1', required=True, type=str)
            required_named.add_argument('--tool2', required=True, type=str)
            required_named.add_argument('--roc-tables', type=str, nargs='+')
            args = parser.parse_args()
            main(args.roc_tables, args.sample_id, args.tool1, args.tool2, args.additional_label)
        EOF

        python script.py --sample-id "~{sample_id}" --tool1 "~{tool1_label}" --tool2 "~{tool2_label}" ~{additional_label_arg} --roc-tables ~{sep=' ' roc_tables}
    >>>

    output {
        Array[File] plots = glob("*.png")
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim-plots:1.0"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk 200 HDD"
    }
}