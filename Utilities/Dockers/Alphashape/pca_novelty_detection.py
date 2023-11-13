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

# Read training data (this is just for visualization in this script)
# Expected input format: tsv with columns PC1 and PC2
df_train = pd.read_csv(args.training_data, sep = '\t', header = 0)
training_points = list(zip(df_train.PC1, df_train.PC2))

# Read test data
# Expected input format: tsv with columns SAMPLE_ID, PC1, and PC2
# Assumption: the header will exist and is considered mandatory for the format
df_test = pd.read_csv(args.input, sep = '\t', header = 0)
test_points = {n: Point(x, y) for n, x, y in zip(df_test.SAMPLE_ID, df_test.PC1, df_test.PC2)}

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
plt.title("Alphashape Automated Novelty Flagging")
plt.xlabel("PC1")
plt.ylabel("PC2")
colors = ['green', 'blue', 'red']
labels = ['Baseline', 'Pass', 'Fail']
patches = [(plt.Line2D([], [], color = colors[i], label = "{:s}".format(labels[i]), marker = "o", linewidth = 0)) for i in range(len(labels))]
plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), handles = patches)
plt.savefig("out_visualization.png", dpi = 300, bbox_inches = 'tight')
