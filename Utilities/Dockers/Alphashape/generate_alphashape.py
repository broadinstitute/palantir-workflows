import argparse
import pickle
import sys
import numpy as np
import pandas as pd
import alphashape

# Parse command-line arguments
parser = argparse.ArgumentParser(description = "Generate and save alphashape with the given training data and the alpha parameter.")
parser.add_argument("-t", "--training-data", required = True)
parser.add_argument("-a", "--alpha", default = 8.0)
args = parser.parse_args()

# Get training data
# Expected input format:
# Tab-delimited file with at least two columns where the first two columns are PC1 and PC2 (in this order)
# SAMPLE_ID might be added as the first column later based on code review
training_points = []
df_train = pd.read_csv(args.training_data, sep = '\t', header = 0)
for row in df_train.itertuples():
	training_points.append((float(row[1]), float(row[2])))

# Generate alphashape
alpha_shape = alphashape.alphashape(training_points, float(args.alpha))

# Pickle alphashape and save to file
with open("alphashape.pickle", "wb") as outfile:
	pickle.dump(alpha_shape, outfile)
