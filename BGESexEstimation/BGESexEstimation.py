import argparse
import csv

class ThresholdsAndLabels:
    def __init__(self, lower_bound_xx, upper_bound_xx, lower_bound_xy, upper_bound_xy, label_xx, label_xy, label_inconclusive):
        self.lower_bound_xx = lower_bound_xx
        self.upper_bound_xx = upper_bound_xx
        self.lower_bound_xy = lower_bound_xy
        self.upper_bound_xy = upper_bound_xy
        self.label_xx = label_xx
        self.label_xy = label_xy
        self.label_inconclusive = label_inconclusive

def read_thresholds_and_labels(thresholds_and_labels_path):
    with open(thresholds_and_labels_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        line = next(reader)
        if line is None:
            raise RuntimeError("Thresholds and labels file is empty.")
        if not all(key in line for key in ['lower_bound_xx', 'upper_bound_xx', 'lower_bound_xy', 'upper_bound_xy', 'label_xx', 'label_xy', 'label_inconclusive']):
            raise RuntimeError("Thresholds and labels file is missing required fields.\n" + 
                               "Expected fields: lower_bound_xx, upper_bound_xx, lower_bound_xy, upper_bound_xy, label_xx, label_xy, label_inconclusive\n" + 
                               "Found fields:    " + ', '.join(line.keys()))
        thresholds_and_labels = ThresholdsAndLabels(
            lower_bound_xx=float(line['lower_bound_xx']),
            upper_bound_xx=float(line['upper_bound_xx']),
            lower_bound_xy=float(line['lower_bound_xy']),
            upper_bound_xy=float(line['upper_bound_xy']),
            label_xx=line['label_xx'],
            label_xy=line['label_xy'],
            label_inconclusive=line['label_inconclusive']
        )
        if next(reader, None) is not None:
            raise RuntimeError("Thresholds and labels file contains more than one line.")
        return thresholds_and_labels

def read_coverage_metrics(metrics_file_path):
    covX = None
    covY = None

    with open(metrics_file_path, 'r') as file:
        reader = csv.DictReader(file, fieldnames=['contig', 'num_reads', 'mean_coverage'])
        for row in reader:
            if row['contig'] == 'chrX':
                covX = float(row['mean_coverage'])
            elif row['contig'] == 'chrY':
                covY = float(row['mean_coverage'])
    return covX, covY

def estimate_sex_ploidy(covX, covY, thresholds_and_labels):
    if covX is None or covY is None:
        raise RuntimeError("Coverage metrics for X or Y chromosome not found.")
    
    if covX == 0 or covY == 0:
        return 0, thresholds_and_labels.label_inconclusive
    
    ratio = covY / covX

    estimated_sex_ploidy = thresholds_and_labels.label_xx if thresholds_and_labels.lower_bound_xx <= ratio <= thresholds_and_labels.upper_bound_xx else (thresholds_and_labels.label_xy if thresholds_and_labels.lower_bound_xy <= ratio <= thresholds_and_labels.upper_bound_xy else thresholds_and_labels.label_inconclusive)
    return ratio, estimated_sex_ploidy

def main(metrics_file_path, thresholds_and_labels_path):
    thresholds_and_labels = read_thresholds_and_labels(thresholds_and_labels_path)
    covX, covY = read_coverage_metrics(metrics_file_path)
    ratio, estimated_sex_ploidy = estimate_sex_ploidy(covX, covY, thresholds_and_labels)

    with open('sex_ploidy_ratio.txt', 'w') as sex_ploidy_ratio_file:
        sex_ploidy_ratio_file.write(str(ratio) + '\n')
    with open('sex_ploidy_estimation.txt', 'w') as sex_ploidy_estimation_file:
        sex_ploidy_estimation_file.write(estimated_sex_ploidy + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BGESexEstimation: Estimate sample sex ploidy from exome region coverage metrics file",
    )
    parser.add_argument('--coverage-metrics', type=str, help='Path to the exome region coverage metrics file.', required=True)
    parser.add_argument('--thresholds-and-labels', type=str, help='Path to the thresholds and labels file.', required=True)

    args = parser.parse_args()
    main(args.coverage_metrics, args.thresholds_and_labels)
