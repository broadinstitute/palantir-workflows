import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys

parser = argparse.ArgumentParser(description = "Generate plot for the isoform reconstruction denovo analysis statistics.")
parser.add_argument("-i", "--input", required = True)
parser.add_argument("-d", "--dataset-name", required = True)
parser.add_argument("-t", "--type", choices =["known", "novel"], required = True)
parser.add_argument("-s", "--save", action = "store_true", required = False)
parser.add_argument("-n", "--no-save", dest = "save", action = "store_false", required = False)
parser.set_defaults(save = True)
args = parser.parse_args()

tools = []
reliable = []
almost_reliable = []
mult_pred = []
unique = []
missed = []

with open(args.input) as infile:
	tools = infile.readline().split('\t')[1:]
	reliable = list(map(int, infile.readline().split('\t')[1:]))
	almost_reliable = list(map(int, infile.readline().split('\t')[1:]))
	mult_pred = list(map(int, infile.readline().split('\t')[1:]))
	unique = list(map(int, infile.readline().split('\t')[1:]))
	missed = list(map(int, infile.readline().split('\t')[1:]))
	missed = list(map(lambda x: -x, missed))

transcript_stats = {
	"reliable": reliable,
	"almost_reliable": almost_reliable,
	"mult_pred": mult_pred,
	"unique": unique,
	"missed": missed,
}

sns.set_style("darkgrid")
colors = sns.color_palette("bright")

x = np.arange(len(tools))
width = 0.15
multiplier = 0

#fig, ax = plt.subplots(layout = "constrained")
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (14, 10), layout = "constrained")

for key, value in transcript_stats.items():
	offset = width * multiplier
	rects = ax.bar(x + offset, value, width, label = key)
	ax.bar_label(rects, padding = 3)
	multiplier += 1

ax.set_ylabel("Count")
ax.set_title("Denovo Analysis Summary: " + args.dataset_name)
ax.set_xticks(x + width, tools)
ax.legend(loc = "upper left", ncols = len(tools))

if args.save == True:
	plt.savefig(args.dataset_name + "_analysis_summary_denovo_" + args.type + ".png")
else:
	plt.show()
