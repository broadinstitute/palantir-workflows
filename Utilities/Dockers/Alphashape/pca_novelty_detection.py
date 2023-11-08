import argparse
import pickle
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from descartes import PolygonPatch
from shapely.geometry import Point
import alphashape

# Parse command-line arguments
parser = argparse.ArgumentParser(description = "Automatically flag novelties in 2D PCA plots using Concave Hulls generated via alphashapes.")
parser.add_argument("-i", "--input", required = True)
parser.add_argument("-a", "--alphashape", required = True)
parser.add_argument("-t", "--training-data", required = True)
parser.add_argument("-d", "--distance-threshold", default = 0.01)
args = parser.parse_args()

# Load alphashape
alphashape_pickle_file = open(args.alphashape, "rb")
alpha_shape = pickle.load(alphashape_pickle_file)
alphashape_pickle_file.close()

# Set the distance threshold
dist_thresh = float(args.distance_threshold)

# Get training data (this is just for visualization in this script)
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

# Get test data
# Expected input format:
# Tab-delimited file with at least three columns where the first three columns are SAMPLE_ID, PC1, and PC2 (in this order)
# We assume that the header exists and is mandatory. This can be changed later if needed
test_points = {}
with open(args.input) as infile:
	header_test = infile.readline()
	for line in infile:
		line = line.rstrip()
		columns = line.split('\t')
		test_point = Point(float(columns[1]), float(columns[2]))
		test_points[columns[0]] = test_point

# Prepare optional output plot for nice visualization
fig, ax = plt.subplots()
ax.scatter(*zip(*training_points), c = 'green', alpha = 1.0, s = 10)
ax.add_patch(PolygonPatch(alpha_shape, alpha = 0.2))

# Open the output file
# Simple name for now, can move to a naming convention later if requested
outfile = open("out.tsv", "w")

# Test and label each sample as a novelty or a regular observation
for sample_id, test_point in test_points.items():
	dist = test_point.distance(alpha_shape)
	if alpha_shape.contains(test_point) or dist < dist_thresh:
		plt.scatter(test_point.x, test_point.y, c = 'blue', alpha = 1.0, s = 10)
		outfile.write(sample_id + "\t" + "PASS" + "\n")
	else:
		plt.scatter(test_point.x, test_point.y, c = 'red', alpha = 1.0, s = 10)
		outfile.write(sample_id + "\t" + "FAIL" + "\n")

# Flush and close the output file
outfile.flush()
outfile.close()

# Optional plotting for nice visualization
plt.title("Alphashape Novelty Detection")
plt.xlabel("PC1")
plt.ylabel("PC2")
colors = ['green', 'blue', 'red']
labels = ['Baseline', 'Pass', 'Fail']
patches = [(plt.Line2D([], [], color = colors[i], label = "{:s}".format(labels[i]), marker = "o", linewidth = 0)) for i in range(len(labels))]
plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), handles = patches)
plt.savefig("out_visualization.png", dpi = 300, bbox_inches = 'tight')
