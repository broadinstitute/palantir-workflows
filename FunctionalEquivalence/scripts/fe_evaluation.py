import matplotlib.pyplot as plt
import argparse
import numpy as np
import re
import matplotlib
import csv
import re

matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.family'] = 'serif'

class FEEvaluation:
    def __init__(self):
        self.data = dict()
        self.datasets = set()

    def read_jaccards_from_summary_csv(self, filename):
        jaccard_data = dict()
        with open(filename) as file:
            for line in csv.DictReader(file):
                if line['Type'] in ('SNP', 'INDEL'):
                    var_type = line['Type'].lower()

                    if line['Stratifier'] == 'NA':
                        region = 'all'
                    elif line['Stratifier'] in ('LCR', 'HCR'):
                        region = line['Stratifier']
                    else:
                        raise RuntimeError('Invalid stratifier {} in file {}'.format(line['Stratifier'], filename))
                    
                    tp_base = int(line['TP_Base'])
                    fp = int(line['FP'])
                    fn= int(line['FN'])
                    jaccard = tp_base / (tp_base + fn + fp)

                    jaccard_data[(var_type, region)] = jaccard
        return jaccard_data

    def _add_to_data(self, filename):
        basename = re.sub(r'\.csv$', '', filename[filename.rfind('/') + 1:])

        evaluation = basename[:5]
        if evaluation not in ('tool1', 'tool2', 'inter'):
            raise RuntimeError('Invalid evaluation in file {}'.format(filename))

        # For tool1/tool2, remove last two digits, for inter, remove the last digit.
        # These are the replicate numbers.
        # Modify here to allow more than 9 replicates
        dataset = basename[5:(-1 if evaluation == 'inter' else -2)]

        jaccard_data = self.read_jaccards_from_summary_csv(filename)

        for var_type in ('snp', 'indel'):
            for region in ('all', 'HCR', 'LCR'):
                index = (dataset, var_type, region, evaluation)
                if index not in self.data:
                    self.data[index] = []
                self.data[index].append(jaccard_data[(var_type, region)])

        self.datasets.add(dataset)
    
    def read_data(self, summaries):
        for filename in summaries:
            self._add_to_data(filename)
        
        for dataset in sorted(self.datasets, reverse=True):
            print('{}: # values for tool1: {} / tool2: {} / inter: {}'.format(dataset,
                len(self.data[dataset, 'snp', 'all', 'tool1']) if (dataset, 'snp', 'all', 'tool1') in self.data else 0,
                len(self.data[dataset, 'snp', 'all', 'tool2']) if (dataset, 'snp', 'all', 'tool2') in self.data else 0,
                len(self.data[dataset, 'snp', 'all', 'inter']) if (dataset, 'snp', 'all', 'inter') in self.data else 0))

    def calculate_data(self, summary_file, dataset, var_type, region):
        tool1_data = self.data[dataset, var_type, region, 'tool1']
        tool2_data = self.data[dataset, var_type, region, 'tool2']
        inter_data = self.data[dataset, var_type, region, 'inter']

        tool1_mean = np.mean(tool1_data)
        tool2_mean = np.mean(tool2_data)
        inter_mean = np.mean(inter_data)

        tool1_sd = np.std(tool1_data, ddof=1)
        tool2_sd = np.std(tool2_data, ddof=1)
        inter_sd = np.std(inter_data, ddof=1)

        print('{} {} {} tool1: mean: {} sd: {} data: {}'.format(dataset, var_type, region, tool1_mean, tool1_sd, tool1_data))
        print('{} {} {} tool2: mean: {} sd: {} data: {}'.format(dataset, var_type, region, tool2_mean, tool2_sd, tool2_data))
        print('{} {} {} inter: mean: {} sd: {} data: {}'.format(dataset, var_type, region, inter_mean, inter_sd, inter_data))

        summary_file.write('{}\t{}\t{}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\n'.format(
            dataset,
            var_type,
            region,
            tool1_mean,
            tool1_sd,
            tool2_mean,
            tool2_sd,
            inter_mean,
            inter_sd
        ))
        return tool1_mean, tool1_sd, tool2_mean, tool2_sd, inter_mean, inter_sd

    def plot_data(self, ax, summary_file, dataset, var_type, region, tool1_label, tool2_label):
        tool1_mean, tool1_sd, tool2_mean, tool2_sd, inter_mean, inter_sd = self.calculate_data(summary_file, dataset, var_type, region)

        x = np.array([0, 1, 2])
        xticklabels = [tool2_label, 'Inter', tool1_label]

        y = np.array([tool2_mean, inter_mean, tool1_mean])
        y_err = np.array([tool2_sd, inter_sd, tool1_sd])

        ax.set_xticks(x)
        ax.set_xticklabels(xticklabels)
        ax.set_ylabel('Jaccard score')
        ax.errorbar(x, y, yerr=y_err, c='k', fmt='o', capsize=15)
        ax.set_xlim(-0.5, 2.5)
        ax.grid(axis='y')

        if inter_mean <= max(tool1_mean, tool2_mean):
            titlecolor = 'orange'
        elif inter_mean - np.nan_to_num(inter_sd) < max(tool1_mean + np.nan_to_num(tool1_sd), tool2_mean + np.nan_to_num(tool2_sd)):
            titlecolor = 'yellow'
        else:
            titlecolor = 'white'
        ax.set_title('{} {}'.format(var_type, region), backgroundcolor=titlecolor, zorder=0)

    
    def plot(self, tool1_label, tool2_label, additional_label):
        with open('fe_summary.tsv', 'w') as summary_file:
            summary_file.write('Dataset\tVar_Type\tRegion\tJ_tool1_mean\tJ_tool1_sd\tJ_tool2_mean\tJ_tool2_sd\tJ_inter_mean\tJ_inter_sd\n')
            for dataset in sorted(self.datasets, reverse=True):
                if (dataset, 'snp', 'all', 'tool1') not in self.data:
                    continue
                fig, axes = plt.subplots(2, 3, figsize=(9,6))
                for row, var_type in enumerate(['snp', 'indel']):
                    for col, region in enumerate(['all', 'HCR', 'LCR']):
                        if (dataset, var_type, region, 'tool1') not in self.data:
                            print('No data for {}'.format((dataset, var_type, region, 'tool1')))
                            continue
                        if (dataset, var_type, region, 'tool2') not in self.data:
                            print('No data for {}'.format((dataset, var_type, region, 'tool2')))
                            continue
                        if (dataset, var_type, region, 'inter') not in self.data:
                            print('No data for {}'.format((dataset, var_type, region, 'inter')))
                            continue

                        self.plot_data(axes[row, col], summary_file, dataset, var_type, region, tool1_label, tool2_label)
                        
                fig.suptitle('Dataset: {}'.format(dataset) + ('' if not additional_label else ', {}'.format(additional_label)) + '\n' + 
                r'Concordance (# values {}: {} / {}: {} / inter: {})'.format(
                    tool1_label,
                    len(self.data[dataset, 'snp', 'all', 'tool1']) if (dataset, 'snp', 'all', 'tool1') in self.data else 0,
                    tool2_label,
                    len(self.data[dataset, 'snp', 'all', 'tool2']) if (dataset, 'snp', 'all', 'tool2') in self.data else 0,
                    len(self.data[dataset, 'snp', 'all', 'inter']) if (dataset, 'snp', 'all', 'inter') in self.data else 0))
                plt.tight_layout()
                fig.savefig('fe_plot_{}.png'.format(dataset), dpi=100)

def main(tool1_label, tool2_label, additional_label, summaries):
    fe = FEEvaluation()
    fe.read_data(summaries)
    fe.plot(tool1_label, tool2_label, additional_label)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create functional equivalence plots.')
    parser.add_argument('--additional-label', type=str)
    required_named = parser.add_argument_group('Required named arguments')
    required_named.add_argument('--tool1', required=True, type=str)
    required_named.add_argument('--tool2', required=True, type=str)
    required_named.add_argument('--summaries', required=True, type=str, nargs='+')
    args = parser.parse_args()
    main(args.tool1, args.tool2, args.additional_label, args.summaries)
