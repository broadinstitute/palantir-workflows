import argparse
import csv

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

def estimate_sex_ploidy(covX, covY):
    lower_bound_xx = 0
    upper_bound_xx = 0.05
    lower_bound_xy = 0.25
    upper_bound_xy = 0.8

    label_xx = 'F'
    label_xy = 'M'
    label_inconclusive = 'U'

    if covX is None or covY is None:
        raise RuntimeError("Coverage metrics for X or Y chromosome not found.")
    
    if covX == 0 or covY == 0:
        return 0, label_inconclusive
    
    ratio = covY / covX

    estimated_sex_ploidy = label_xx if lower_bound_xx <= ratio <= upper_bound_xx else (label_xy if lower_bound_xy <= ratio <= upper_bound_xy else label_inconclusive)
    return ratio, estimated_sex_ploidy

def main(metrics_file_path):
    covX, covY = read_coverage_metrics(metrics_file_path)
    ratio, estimated_sex_ploidy = estimate_sex_ploidy(covX, covY)

    with open('sex_ploidy_ratio.txt', 'w') as sex_ploidy_ratio_file:
        sex_ploidy_ratio_file.write(str(ratio) + '\n')
    with open('sex_ploidy_estimation.txt', 'w') as sex_ploidy_estimation_file:
        sex_ploidy_estimation_file.write(estimated_sex_ploidy + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BGESexEstimation: Estimate sample sex ploidy from exome region coverage metrics file",
    )
    parser.add_argument(
        "coverage_metrics",
        type=str,
        help="Path to the exome region coverage metrics file."
    )
    args = parser.parse_args()
    main(args.coverage_metrics)