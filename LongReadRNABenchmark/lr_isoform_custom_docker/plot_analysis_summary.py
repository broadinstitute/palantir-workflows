import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys

parser = argparse.ArgumentParser(description = "Generate plot for the isoform reconstruction analysis summary.")
parser.add_argument("-i", "--input", required = True)
parser.add_argument("-d", "--dataset-name", required = True)
parser.add_argument("-t", "--type", required = True)
parser.add_argument("-s", "--save", action = "store_true", required = False)
parser.add_argument("-n", "--no-save", dest = "save", action = "store_false", required = False)
parser.set_defaults(save = True)
args = parser.parse_args()

tools = []
sensitivity_list = []
precision_list = []
f1_score_list = []

with open(args.input) as infile:
	header = infile.readline()
	for line in infile:
		line = line.rstrip()
		tokens = line.split('\t')
		tools.append(tokens[0])
		sensitivity_list.append(float(tokens[1]))
		precision_list.append(float(tokens[2]))
		f1_score_list.append(float(tokens[3]))
		
sns.set_style("darkgrid")
colors = sns.color_palette("bright")

fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (12, 8), tight_layout = True)

fig.suptitle("Analysis Summary: " + args.dataset_name)

ax[0].set_ylim(0.0, 1.0)
ax[1].set_ylim(0.0, 1.0)
ax[2].set_ylim(0.0, 1.0)

ax[0].set_ylabel("Sensitivity")
ax[1].set_ylabel("Precision")
ax[2].set_ylabel("F-1 Score")

ax[0].bar(tools, sensitivity_list, color = colors[1:7])
ax[1].bar(tools, precision_list, color = colors[1:7])
ax[2].bar(tools, f1_score_list, color = colors[1:7])

if args.save == True:
	plt.savefig(args.dataset_name + "_analysis_summary_" + args.type + ".png")
else:
	plt.show()
