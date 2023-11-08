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
with open(args.alphashape, "rb") as infile:
	alpha_shape = pickle.load(infile)

# Set the distance threshold
dist_thresh = float(args.distance_threshold)

# Get training data (this is just for visualization in this script)
# Expected input format:
# Tab-delimited file with at least two columns where the first two columns are PC1 and PC2 (in this order)
# SAMPLE_ID might be added as the first column later based on code review
training_points = []
df_train = pd.read_csv(args.training_data, sep = '\t', header = 0)
for row in df_train.itertuples():
	training_points.append((float(row[1]), float(row[2])))

# Get test data
# Expected input format:
# Tab-delimited file with at least three columns where the first three columns are SAMPLE_ID, PC1, and PC2 (in this order)
# We assume that the header exists and is mandatory. This can be changed later if needed
test_points = {}
df_test = pd.read_csv(args.input, sep = '\t', header = 0)
for row in df_test.itertuples():
	test_points[row[1]] = Point(float(row[2]), float(row[3]))

# Prepare output plot for nice visualization
fig, ax = plt.subplots()
ax.scatter(*zip(*training_points), c = 'green', alpha = 1.0, s = 10)
ax.add_patch(PolygonPatch(alpha_shape, alpha = 0.2))

# Simple name for now, can move to a naming convention later if requested
with open("out.tsv", "w") as outfile:
	# Test and label each sample as a novelty or a regular observation
	for sample_id, test_point in test_points.items():
		dist = test_point.distance(alpha_shape)
		if alpha_shape.contains(test_point) or dist < dist_thresh:
			outfile.write(str(sample_id) + "\t" + "PASS" + "\n")
			plt.scatter(test_point.x, test_point.y, c = 'blue', alpha = 1.0, s = 10)
		else:
			outfile.write(str(sample_id) + "\t" + "FAIL" + "\n")
			plt.scatter(test_point.x, test_point.y, c = 'red', alpha = 1.0, s = 10)

# Plotting for nice visualization
plt.title("Alphashape Novelty Detection")
plt.xlabel("PC1")
plt.ylabel("PC2")
colors = ['green', 'blue', 'red']
labels = ['Baseline', 'Pass', 'Fail']
patches = [(plt.Line2D([], [], color = colors[i], label = "{:s}".format(labels[i]), marker = "o", linewidth = 0)) for i in range(len(labels))]
plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), handles = patches)
plt.savefig("out_visualization.png", dpi = 300, bbox_inches = 'tight')
