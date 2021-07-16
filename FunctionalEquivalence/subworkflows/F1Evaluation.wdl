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
        
        cat <<EOF > script.py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import csv
import statistics as stat
import gzip
import argparse

matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.family'] = 'serif'

def f1(tp, fp, fn):
    return tp / (tp + 0.5 * (fp + fn))

def read_roc_to_f1_dict(filename):
    f1_mean_dict = dict()
    f1_list_dict = dict()
    with gzip.open(filename, 'rt') as roc_file:
        for line in roc_file:
            if line.startswith('#score field'):
                break
        for line in csv.DictReader(roc_file, delimiter='\t'):
            f1_score = f1(float(line['true_positives_baseline']), float(line['false_positives']), float(line['false_negatives']))
            qual = round(float(line['#score']))
            if qual not in f1_list_dict:
                f1_list_dict[qual] = [f1_score]
            else:
                f1_list_dict[qual].append(f1_score)
    
    for qual, f1_list in f1_list_dict.items():
        f1_mean_dict[qual] = stat.mean(f1_list)
                                    
    return f1_mean_dict

def read_file(data, filename):
    basename = filename[filename.rfind('/') + 1:]
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
    
    f1_data = read_roc_to_f1_dict(filename)
    
    if dataset not in data:
        data[dataset] = dict()
    if replicate not in data[dataset]:
        data[dataset][replicate] = {'tool1': {'all': dict(), 'HCR': dict(), 'LCR': dict()},
                                    'tool2': {'all': dict(), 'HCR': dict(), 'LCR': dict()}, }
    data[dataset][replicate][tool][region][var_type] = f1_data

def extract_data_from_dict(score_range, data, dataset, tool, region, var_type):
    return [
                np.array([
                    data[dataset][replicate][tool][region][var_type][score] if score in data[dataset][replicate][tool][region][var_type] else np.nan
                    for replicate in data[dataset].keys()
                ])
                for score in score_range
            ]
    


def aggregate_data(data, dataset, signed_difference):
    max_score = 200
    # Don't start the range at anything other than 0 because it will mess up the conversion to numpy arrays
    score_range = range(0,max_score)

    aggregated_data = dict()
    
    for region in ('all', 'HCR', 'LCR'):
        for var_type in ('snp', 'indel'):
            for dataset in (dataset,):
                tool1_data = extract_data_from_dict(score_range, data, dataset, 'tool1', region, var_type)
                tool2_data = extract_data_from_dict(score_range, data, dataset, 'tool2', region, var_type)
                
                means_tool1 = [np.nanmean(tool1_data[score]) if not np.isnan(tool1_data[score]).all() else np.nan for score in score_range]
                sd_tool1 = [np.std(tool1_data[score], ddof=1) for score in score_range]
                means_tool2 = [np.nanmean(tool2_data[score]) if not np.isnan(tool2_data[score]).all() else np.nan for score in score_range]
                sd_tool2 = [np.std(tool2_data[score], ddof=1) for score in score_range]
                
                tool1_differences = dict()
                tool2_differences = dict()
                num_replicates = len(data[dataset].keys())
                
                for score in score_range:
                    tool1_differences[score] = []
                    tool2_differences[score] = []
                    for replicate_index_i in range(num_replicates):
                        for replicate_index_j in range(replicate_index_i + 1, num_replicates):
                            tool1_differences[score].append(abs(tool1_data[score][replicate_index_i] - tool1_data[score][replicate_index_j]))
                            tool2_differences[score].append(abs(tool2_data[score][replicate_index_i] - tool2_data[score][replicate_index_j]))
                
                intra_difference_means_tool1 = [np.nanmean(tool1_differences[score]) for score in score_range]
                if num_replicates > 2:
                    intra_difference_sd_tool1 = [np.std(tool1_differences[score], ddof=1) for score in score_range]
                else:
                    intra_difference_sd_tool1 = [0] * max_score
                
                intra_difference_means_tool2 = [np.nanmean(tool2_differences[score]) for score in score_range]
                if num_replicates > 2:
                    intra_difference_sd_tool2 = [np.std(tool2_differences[score], ddof=1) for score in score_range]
                else:
                    intra_difference_sd_tool2 = [0] * max_score
                
                
                inter_difference_means = [
                    np.nanmean(
                        (tool1_data[score] - tool2_data[score]) if signed_difference else np.abs(tool1_data[score] - tool2_data[score])
                    ) for score in score_range
                ]
                inter_difference_sd = [
                    np.std(
                        (tool1_data[score] - tool2_data[score]) if signed_difference else np.abs(tool1_data[score] - tool2_data[score]), ddof=1
                    ) for score in score_range
                ]
                
                if dataset not in aggregated_data:
                    aggregated_data[dataset] = dict()
                aggregated_data[dataset][(region, var_type)] = {
                    'means_tool1': np.array(means_tool1),
                    'sd_tool1': np.array(sd_tool1),
                    'means_tool2': np.array(means_tool2),
                    'sd_tool2': np.array(sd_tool2),
                    'intra_difference_means_tool1': np.array(intra_difference_means_tool1),
                    'intra_difference_sd_tool1': np.array(intra_difference_sd_tool1),
                    'intra_difference_means_tool2': np.array(intra_difference_means_tool2),
                    'intra_difference_sd_tool2': np.array(intra_difference_sd_tool2),
                    'inter_difference_means': np.array(inter_difference_means),
                    'inter_difference_sd': np.array(inter_difference_sd)
                }
    return aggregated_data
    

def read_datasets(roc_tables):
    data = dict()
    for filename in roc_tables:
        read_file(data, filename)
    return data

def diff_plot(ax, data, dataset, var_type, region, tool1_label, tool2_label, signed_difference, summary_file):
    X = np.arange(0, 200)
    ax.plot(X, data['inter_difference_means'], label='Inter: $|F_{1, ' + tool1_label + ' rep 1} - F_{1, ' + tool2_label + ' rep 1}|$', c='C2')
    ax.fill_between(X, data['inter_difference_means'] - data['inter_difference_sd'], data['inter_difference_means'] + data['inter_difference_sd'], color='C2', alpha=0.1)

    ax.plot(X, data['intra_difference_means_tool1'], label='Intra: $|F_{1, ' + tool1_label + ' rep 1} - F_{1, ' + tool1_label + ' rep 2}|$', c='C0')
    ax.fill_between(X, data['intra_difference_means_tool1'] - data['intra_difference_sd_tool1'], data['intra_difference_means_tool1'] + data['intra_difference_sd_tool1'], color='C0', alpha=0.1)

    ax.plot(X, data['intra_difference_means_tool2'], label='Intra: $|F_{1, ' + tool2_label + ' rep 1} - F_{1, ' + tool2_label + ' rep 2}|$', c='C1')
    ax.fill_between(X, data['intra_difference_means_tool2'] - data['intra_difference_sd_tool2'], data['intra_difference_means_tool2'] + data['intra_difference_sd_tool2'], color='C1', alpha=0.1)

    max_y = max(
        np.nanmax(data['inter_difference_means'][0:30] + data['inter_difference_sd'][0:30]),
        np.nanmax(data['intra_difference_means_tool1'][0:30] + data['intra_difference_sd_tool1'][0:30]),
        np.nanmax(data['intra_difference_means_tool2'][0:30] + data['intra_difference_sd_tool2'][0:30]))

    if signed_difference:
        min_y = np.nanmin(data['inter_difference_means'][0:30] - data['inter_difference_sd'][0:30])
    else:
        min_y = 0

    if np.isfinite(max_y) and np.isfinite(min_y):
        ax.set_ylim(min(0, min_y * 1.1), max(0, max_y * 1.1))
    else:
        print('Invalid plot data for {} {} {}'.format(dataset, var_type, region))

    ax.set_xlim(0, 30)
    ax.set_xlabel('Quality threshold q')

    if signed_difference:
        ax.set_ylabel(r'$\Delta\;F_1$')
    else:
        ax.set_ylabel(r'$|\Delta\;F_1|$')

    if signed_difference:
        ax.axhline(y=0, c='grey', linewidth=1)

    if np.any(data['inter_difference_means'][0:30] > data['intra_difference_means_tool1'][0:30]) or np.any(data['inter_difference_means'][0:30] > data['intra_difference_means_tool2'][0:30]):
        titlecolor = 'orange'
    else:
        titlecolor = 'white'
        
    ax.set_title(r'{} {}'.format(var_type, region), backgroundcolor=titlecolor, zorder=0)

    for threshold in range(30):
        summary_file.write('{}\t{}\t{}\t{}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\n'.format(
            dataset,
            var_type,
            region,
            threshold,
            data['intra_difference_means_tool1'][threshold],
            data['intra_difference_sd_tool1'][threshold],
            data['intra_difference_means_tool2'][threshold],
            data['intra_difference_sd_tool2'][threshold],
            data['inter_difference_means'][threshold],
            data['inter_difference_sd'][threshold]
        ))

def plot_aggregated_data_compact(data, dataset, tool1_label, tool2_label, additional_label, signed_difference, summary_file):
    aggregated_data = aggregate_data(data, dataset, signed_difference)
    fig, axes = plt.subplots(3, 3, figsize=(9,8))
    for row, var_type in enumerate(['snp', 'indel']):
        for col, region in enumerate(['all', 'HCR', 'LCR']):
            plot_data = aggregated_data[dataset][(region, var_type)]
            ax = axes[row, col]
            diff_plot(ax, plot_data, dataset, var_type, region, tool1_label, tool2_label, signed_difference, summary_file)
    
    axes[2, 0].axis('off')
    axes[2, 1].axis('off')
    axes[2, 2].axis('off')
            
    fig.suptitle('Dataset: {}'.format(dataset) + ('' if not additional_label else ', {}'.format(additional_label)) + '\n' +
                ('Signed ' if signed_difference else 'Absolute ') +
                r'$F_1 = \frac{TP}{TP + \frac{1}{2} (FP + FN)}$ score differences for calls with $QUAL \geq q$' +
                ' (# replicates: {})'.format(len(data[dataset])))
    if signed_difference:
        inter_label = 'Inter: $F_{1, ' + tool1_label + ' rep 1} - F_{1, ' + tool2_label + ' rep 1}$'
    else:
        inter_label = 'Inter: $|F_{1, ' + tool1_label + ' rep 1} - F_{1, ' + tool2_label + ' rep 1}|$'
    legend_inter_line = matplotlib.lines.Line2D([], [], color='C2', label=inter_label)
    legend_intra_tool1_line = matplotlib.lines.Line2D([], [], color='C0', label='Intra: $|F_{1, ' + tool1_label + ' rep 1} - F_{1, ' + tool1_label + ' rep 2}|$')
    legend_intra_tool2_line = matplotlib.lines.Line2D([], [], color='C1', label='Intra: $|F_{1, ' + tool2_label + ' rep 1} - F_{1, ' + tool2_label + ' rep 2}|$')
    fig.legend(bbox_to_anchor=(0.5, 0.2), loc='center', handles=[legend_inter_line, legend_intra_tool1_line, legend_intra_tool2_line])
    plt.tight_layout()
    fig.savefig('f1_plot_{}.png'.format(dataset), dpi=100)

def main(roc_tables, tool1_label, tool2_label, additional_label, signed_difference):
    data = read_datasets(roc_tables)
    datasets = data.keys()
    with open('f1_summary.tsv', 'w') as summary_file:
        summary_file.write('Dataset\tVar_Type\tRegion\tThreshold\tabs_deltaF1_tool1_mean\tabs_deltaF1_tool1_sd\tabs_deltaF1_tool2_mean\tabs_deltaF1_tool2_sd\tabs_deltaF1_inter_mean\tabs_deltaF1_inter_sd\n')
        for dataset in datasets:
            plot_aggregated_data_compact(data, dataset, tool1_label, tool2_label, additional_label, signed_difference, summary_file)

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
    }

    # Disable call caching since we fetch the script from a HTTP(S) URL which might have changed
    meta {
        volatile: true
    }

    runtime {
        docker: "michaelgatzen/fe_evaluation"
        preemptible: select_first([preemptible, 0])
        memory: machine_mem_gb + " GB"
        disks: "local-disk 20 HDD"
    }
}