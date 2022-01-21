version 1.0

workflow F1Evaluation {
    input {
        Array[File] roc_tables
        String tool1_label
        String tool2_label
        String? additional_label
        Boolean signed_difference = false
        Int? mem_gb
        Int? preemptible
    }
    
    call F1EvaluationTask {
        input:
            roc_tables = roc_tables,
            tool1_label = tool1_label,
            tool2_label = tool2_label,
            additional_label = additional_label,
            signed_difference = signed_difference,
            mem_gb = mem_gb,
            preemptible = preemptible
    }

    output {
        Array[File] f1_plots = F1EvaluationTask.f1_plots
        File f1_summary = F1EvaluationTask.f1_summary
        Int fe_status = F1EvaluationTask.fe_status
    }
}


task F1EvaluationTask {
    input {
        Array[File] roc_tables
        String tool1_label
        String tool2_label
        String? additional_label
        Boolean signed_difference = false
        Int? mem_gb
        Int? preemptible
    }
    
    Int machine_mem_gb = select_first([mem_gb, 8])

    String additional_label_arg = if defined(additional_label) then "--additional-label \"" + additional_label + "\"" else ""

    command <<<
        set -xeuo pipefail
        
        source activate fe_evaluation
        
        cat <<'EOF' > script.py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd
import numpy as np

matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.family'] = 'serif'

plot_qual_limit = 30

def f1(tp, fp, fn):
    return tp / (tp + 0.5 * (fp + fn))

def read_file(filename):
    basename = os.path.basename(filename)
    info = basename.split('_')
    evalinfo = info[0].split('.')
    dataset = evalinfo[0]
    replicate = evalinfo[1]
    tool = evalinfo[2]
    
    if info[2] == 'vcfeval':
        region = 'all'
    else:
        region = info[2]
    if info[3 if region == 'all' else 4] == 'non':
        var_type = 'indel'
    else:
        var_type = 'snp'
    
    file_data = pd.read_csv(filename, header=6, delimiter='\t')
    file_data['region'] = region
    file_data['var_type'] = var_type
    file_data['dataset'] = dataset
    file_data['replicate'] = replicate
    file_data['tool'] = tool
    file_data['f1'] = f1(file_data['true_positives_baseline'], file_data['false_positives'], file_data['false_negatives'])
    file_data.rename(columns={'#score': 'score'}, inplace=True)
    return file_data

def read_datasets(roc_tables):
    dataframes = []
    for filename in roc_tables:
        dataframes.append(read_file(filename))
    return pd.concat(dataframes)

def filter_data(data, dataset, var_type, region):
    return data.loc[(data['dataset'] == dataset) & (data['var_type'] == var_type) & (data['region'] == region)]

def query_at_score(data, score):
    # Data is sorted descending, so by supplying a reversed indexing list to searchsorted we achieve an ascending sort.
    # Consequently, we need to return len(data) - index as the actual index.
    index = data['score'].searchsorted(score, side='right', sorter=list(reversed(np.arange(len(data)))))

    # If a QUAL score is requested that is higher than the highest one in the file, return NaN
    if len(data) - index - 1 < 0:
        return np.nan
    
    # Otherwise return the value at the next sampled threshold
    return data.iloc[len(data) - index - 1]['f1']

def compute_differences(region_data, signed_difference):
    inter_difference_means_list = []
    inter_difference_sd_list = []
    intra_difference_means_tool1_list = []
    intra_difference_sd_tool1_list = []
    intra_difference_means_tool2_list = []
    intra_difference_sd_tool2_list = []

    X = np.arange(0, plot_qual_limit + 1)
    for x in X:
        tool1_score_differences = []
        tool2_score_differences = []
        inter_score_differences = []
        replicates = region_data.replicate.unique()
        for replicate_i in replicates:
            replicate_i_tool1_score = query_at_score(region_data.loc[(region_data['replicate'] == replicate_i) & (region_data['tool'] == 'tool1')], x)
            replicate_i_tool2_score = query_at_score(region_data.loc[(region_data['replicate'] == replicate_i) & (region_data['tool'] == 'tool2')], x)
            for replicate_j in replicates:
                if replicate_i == replicate_j:
                    continue
                replicate_j_tool1_score = query_at_score(region_data.loc[(region_data['replicate'] == replicate_j) & (region_data['tool'] == 'tool1')], x)
                replicate_j_tool2_score = query_at_score(region_data.loc[(region_data['replicate'] == replicate_j) & (region_data['tool'] == 'tool2')], x)

                tool1_score_differences.append(abs(replicate_i_tool1_score - replicate_j_tool1_score))
                tool2_score_differences.append(abs(replicate_i_tool2_score - replicate_j_tool2_score))

            if signed_difference:
                inter_score_differences.append(replicate_i_tool1_score - replicate_i_tool2_score)
            else:
                inter_score_differences.append(abs(replicate_i_tool1_score - replicate_i_tool2_score))
        inter_difference_means_list.append(np.nanmean(inter_score_differences))
        inter_difference_sd_list.append(np.nanstd(inter_score_differences, ddof=1))
        intra_difference_means_tool1_list.append(np.nanmean(tool1_score_differences))
        intra_difference_sd_tool1_list.append(np.nanstd(tool1_score_differences, ddof=1))
        intra_difference_means_tool2_list.append(np.nanmean(tool2_score_differences))
        intra_difference_sd_tool2_list.append(np.nanstd(tool2_score_differences, ddof=1))
    
    inter_difference_means = np.array(inter_difference_means_list)
    inter_difference_sd = np.array(inter_difference_sd_list)
    intra_difference_means_tool1 = np.array(intra_difference_means_tool1_list)
    intra_difference_sd_tool1 = np.array(intra_difference_sd_tool1_list)
    intra_difference_means_tool2 = np.array(intra_difference_means_tool2_list)
    intra_difference_sd_tool2 = np.array(intra_difference_sd_tool2_list)

    return X, inter_difference_means, inter_difference_sd, intra_difference_means_tool1, intra_difference_sd_tool1, intra_difference_means_tool2, intra_difference_sd_tool2

def plot_region(ax, region_data, dataset, var_type, region, signed_difference, summary_file):
    X, \
        inter_difference_means, \
        inter_difference_sd, \
        intra_difference_means_tool1, \
        intra_difference_sd_tool1, \
        intra_difference_means_tool2, \
        intra_difference_sd_tool2 \
        = compute_differences(region_data, signed_difference)
    
    min_qual_score_tool1 = int(region_data.loc[region_data['tool'] == 'tool1'].score.min() + 1)
    min_qual_score_tool2 = int(region_data.loc[region_data['tool'] == 'tool2'].score.min() + 1)
    min_inter_qual_score = max(min_qual_score_tool1, min_qual_score_tool2)

    # Draw inter curve
    # If there are no variants with QUAL < plot_qual_limit then draw a dotted line with
    # the F1 score of the lowest QUAL threshold encountered
    if min_inter_qual_score > plot_qual_limit:
        ax.plot(X, [inter_difference_means[plot_qual_limit]] * (plot_qual_limit + 1), c='C2', linestyle='dotted')
        ax.fill_between(X, [inter_difference_means[plot_qual_limit] - inter_difference_sd[plot_qual_limit]] * (plot_qual_limit + 1), [inter_difference_means[plot_qual_limit] + inter_difference_sd[plot_qual_limit]] * (plot_qual_limit + 1), color='C2', alpha=0.1)
    else:
        ax.plot(X[min_inter_qual_score:], inter_difference_means[min_inter_qual_score:], c='C2')
        ax.fill_between(X[min_inter_qual_score:], inter_difference_means[min_inter_qual_score:] - inter_difference_sd[min_inter_qual_score:], inter_difference_means[min_inter_qual_score:] + inter_difference_sd[min_inter_qual_score:], color='C2', alpha=0.1)

    # Draw tool1 curve
    if min_qual_score_tool1 > plot_qual_limit:
        ax.plot(X, [intra_difference_means_tool1[plot_qual_limit]] * (plot_qual_limit + 1), c='C0', linestyle='dotted')
        ax.fill_between(X, [intra_difference_means_tool1[plot_qual_limit] - intra_difference_sd_tool1[plot_qual_limit]] * (plot_qual_limit + 1), [intra_difference_means_tool1[plot_qual_limit] + intra_difference_sd_tool1[plot_qual_limit]] * (plot_qual_limit + 1), color='C0', alpha=0.1)
    else:
        ax.plot(X[min_qual_score_tool1:], intra_difference_means_tool1[min_qual_score_tool1:], c='C0')
        ax.fill_between(X[min_qual_score_tool1:], intra_difference_means_tool1[min_qual_score_tool1:] - intra_difference_sd_tool1[min_qual_score_tool1:], intra_difference_means_tool1[min_qual_score_tool1:] + intra_difference_sd_tool1[min_qual_score_tool1:], color='C0', alpha=0.1)
    
    # Draw tool2 curve
    if min_qual_score_tool2 > plot_qual_limit:
        ax.plot(X, [intra_difference_means_tool2[plot_qual_limit]] * (plot_qual_limit + 1), c='C1', linestyle='dotted')
        ax.fill_between(X, [intra_difference_means_tool2[plot_qual_limit] - intra_difference_sd_tool2[plot_qual_limit]] * (plot_qual_limit + 1), [intra_difference_means_tool2[plot_qual_limit] + intra_difference_sd_tool2[plot_qual_limit]] * (plot_qual_limit + 1), color='C1', alpha=0.1)
    else:
        ax.plot(X[min_qual_score_tool2:], intra_difference_means_tool2[min_qual_score_tool2:], c='C1')
        ax.fill_between(X[min_qual_score_tool2:], intra_difference_means_tool2[min_qual_score_tool2:] - intra_difference_sd_tool2[min_qual_score_tool2:], intra_difference_means_tool2[min_qual_score_tool2:] + intra_difference_sd_tool2[min_qual_score_tool2:], color='C1', alpha=0.1)

    # Find axis limits
    max_y = max(
        np.nanmax(inter_difference_means + inter_difference_sd),
        np.nanmax(intra_difference_means_tool1 + intra_difference_sd_tool1),
        np.nanmax(intra_difference_means_tool2 + intra_difference_sd_tool2))

    if signed_difference:
        min_y = np.nanmin(inter_difference_means - inter_difference_sd)
    else:
        min_y = 0

    if np.isfinite(max_y) and np.isfinite(min_y):
        ax.set_ylim(min(0, min_y * 1.1), max(0, max_y * 1.1))
    else:
        print('Invalid plot data for {} {} {}'.format(dataset, var_type, region))

    ax.set_xlim(0, plot_qual_limit)
    ax.set_xlabel('Quality threshold q')

    if signed_difference:
        ax.set_ylabel(r'$\Delta\;F_1$')
    else:
        ax.set_ylabel(r'$|\Delta\;F_1|$')

    # If using signed_difference, draw a line to indicate 0 on the vertical axis
    if signed_difference:
        ax.axhline(y=0, c='grey', linewidth=1)

    # Check the FE criterion
    if np.any(inter_difference_means > intra_difference_means_tool1) or np.any(inter_difference_means > intra_difference_means_tool2):
        titlecolor = 'orange'
        fe_status = 2
    else:
        titlecolor = 'white'
        fe_status = 0
        
    ax.set_title(r'{} {}'.format(var_type, region), backgroundcolor=titlecolor, zorder=0)

    # Write values to the summary file
    for threshold in range(plot_qual_limit):
        summary_file.write('{}\t{}\t{}\t{}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\n'.format(
            dataset,
            var_type,
            region,
            threshold,
            intra_difference_means_tool1[threshold],
            intra_difference_sd_tool1[threshold],
            intra_difference_means_tool2[threshold],
            intra_difference_sd_tool2[threshold],
            inter_difference_means[threshold],
            inter_difference_sd[threshold]
        ))
    return fe_status

def plot_dataset(data, dataset, tool1_label, tool2_label, additional_label, signed_difference, summary_file):
    fe_status = 0
    stratifiers = data.region.unique()

    num_columns = max(len(stratifiers), 3)
    fig, axes = plt.subplots(3, num_columns, figsize=(3*num_columns,8))
    for row, var_type in enumerate(['snp', 'indel']):
        for col, region in enumerate(stratifiers):
            column_to_plot = col if len(stratifiers) > 1 else 1
            ax = axes[row, column_to_plot]
            fe_status = max(fe_status, plot_region(ax, filter_data(data, dataset, var_type, region), dataset, var_type, region, signed_difference, summary_file))
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
    
            
    fig.suptitle('Dataset: {}'.format(dataset) + ('' if not additional_label else ', {}'.format(additional_label)) + '\n' +
                ('Signed ' if signed_difference else 'Absolute ') +
                r'$F_1 = \frac{TP}{TP + \frac{1}{2} (FP + FN)}$ score differences for calls with $QUAL \geq q$' +
                ' (# replicates: {})'.format(len(data.loc[data['dataset'] == dataset].replicate.unique())))
    if signed_difference:
        inter_label = 'Inter: $F_{1, ' + tool1_label + '\\ rep 1} - F_{1, ' + tool2_label + '\\ rep 1}$'
    else:
        inter_label = 'Inter: $|F_{1, ' + tool1_label + '\\ rep 1} - F_{1, ' + tool2_label + '\\ rep 1}|$'
    legend_inter_line = matplotlib.lines.Line2D([], [], color='C2', label=inter_label)
    legend_intra_tool1_line = matplotlib.lines.Line2D([], [], color='C0', label='Intra: $|F_{1, ' + tool1_label + '\\ rep 1} - F_{1, ' + tool1_label + '\\ rep 2}|$')
    legend_intra_tool2_line = matplotlib.lines.Line2D([], [], color='C1', label='Intra: $|F_{1, ' + tool2_label + '\\ rep 1} - F_{1, ' + tool2_label + '\\ rep 2}|$')
    fig.legend(bbox_to_anchor=(0.5, 0.2), loc='center', handles=[legend_inter_line, legend_intra_tool1_line, legend_intra_tool2_line])
    plt.tight_layout()
    fig.savefig('f1_plot_{}.png'.format(dataset), dpi=100)
    return fe_status

def main(roc_tables, tool1_label, tool2_label, additional_label, signed_difference):
    fe_status = 0

    data = read_datasets(roc_tables)
    datasets = data.dataset.unique()
    with open('f1_summary.tsv', 'w') as summary_file:
        summary_file.write('Dataset\tVar_Type\tRegion\tThreshold\tabs_deltaF1_tool1_mean\tabs_deltaF1_tool1_sd\tabs_deltaF1_tool2_mean\tabs_deltaF1_tool2_sd\tabs_deltaF1_inter_mean\tabs_deltaF1_inter_sd\n')
        for dataset in datasets:
            fe_status = max(fe_status, plot_dataset(data, dataset, tool1_label, tool2_label, additional_label, signed_difference, summary_file))
    with open('fe_status.txt', 'w') as fe_status_file:
        fe_status_file.write(str(fe_status))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create F1 functional equivalence plots.')
    parser.add_argument('--additional-label', type=str)
    parser.add_argument('--signed-difference', action='store_true')
    required_named = parser.add_argument_group('Required named arguments')
    required_named.add_argument('--tool1', required=True, type=str)
    required_named.add_argument('--tool2', required=True, type=str)
    required_named.add_argument('--roc-tables', required=True, type=str, nargs='+')
    args = parser.parse_args()
    main(args.roc_tables, args.tool1, args.tool2, args.additional_label, args.signed_difference)
EOF
        
        python script.py --tool1 "~{tool1_label}" --tool2 "~{tool2_label}" ~{additional_label_arg} ~{true="--signed-difference" false="" signed_difference} --roc-tables ~{sep=' ' roc_tables}
    >>>

    output {
        Array[File] f1_plots = glob("*.png")
        File f1_summary = "f1_summary.tsv"
        Int fe_status = read_int("fe_status.txt")
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/functionalequivalence/fe_evaluation:1.0.0"
        preemptible: select_first([preemptible, 0])
        memory: machine_mem_gb + " GB"
        disks: "local-disk 20 HDD"
    }
}