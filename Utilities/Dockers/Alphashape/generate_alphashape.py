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
with open(args.training_data) as infile:
	header_training = infile.readline()
	for line in infile:
		line = line.rstrip()
		columns = line.split('\t')
		training_points.append((float(columns[0]), float(columns[1])))

# Generate alphashape
alpha = float(args.alpha)
alpha_shape = alphashape.alphashape(training_points, alpha)

# Pickle alphashape and save to file
outfile = open("alphashape.pickle", "wb")
pickle.dump(alpha_shape, outfile)
outfile.close()
