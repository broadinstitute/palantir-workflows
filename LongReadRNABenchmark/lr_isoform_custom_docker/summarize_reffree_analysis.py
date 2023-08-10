import argparse
import math
import sys

def get_accuracy_metrics(path):
	infile = open(path)

	lines = infile.read().split('\n')

	# A proper gffcompare <outprefix>.stats file has 31 lines
	# If the isoform detection tool has made no calls at all, you get an "empty" file with 12 lines
	# Since the file format is pretty much hard coded, this breaks the MultiQC-based parser
	# Solution: if one of the output files has 12 lines, just return 0.0 for sensitivity, precision, and the f-1 score
	if len(lines) == 12:
		return (0.0, 0.0, 0.0)

	# Parse the file for accuracy metrics (Sensitivity/Precision)
	gffcompare_data = {}
	gffcompare_data["accuracy"] = {}
	for line in lines[10:16]:
		split = line.replace("|", "").replace("Intron chain", "Intron_chain").split()
		gffcompare_data["accuracy"][split[0]] = {}
		gffcompare_data["accuracy"][split[0]]["sensitivity"] = float(split[2])
		gffcompare_data["accuracy"][split[0]]["precision"] = float(split[3])

	# We just want the accuracy metrics at the Transcript level
	sensitivity = gffcompare_data["accuracy"]["Transcript"]["sensitivity"]
	precision = gffcompare_data["accuracy"]["Transcript"]["precision"]

	# Compute the F1-Score using the precision and sensitivity
	f1_score = 0.0
	if not math.isclose(precision + sensitivity, 0.0):
		f1_score = 2 * precision * sensitivity / (precision + sensitivity)

	# Return the result as a tuple
	return (sensitivity, precision, f1_score)

# Main Program
parser = argparse.ArgumentParser(description = "Parse the main gffcompare output from various tools and summarize isoform reconstruction results.")
parser.add_argument("-i", "--input-list", nargs = '+', help = "List of <outprefix>.stats files that were obtained from the gffcompare output", required = True)
parser.add_argument("-t", "--tool-names", nargs = '+', help = "List of tool names corresponding to the tools used for input-list. Must be in the same order.", required = True)
parser.add_argument("-d", "--dataset-name", required = True)
args = parser.parse_args()

# Only the length match is checked
# It is the user's responsibility to make sure that the order of the input files match the order of the tool names
if len(args.input_list) != len(args.tool_names):
	print("Length of the input list must be equal to the length of tool names list. Order must also be the same, but this is not checked by the program!")
	sys.exit(0)

sensitivity = []
precision = []
f1_score = []
for i in range(len(args.input_list)):
	accuracy_metrics = get_accuracy_metrics(args.input_list[i])
	sensitivity.append(accuracy_metrics[0] / 100.0)
	precision.append(accuracy_metrics[1] / 100.0)
	f1_score.append(accuracy_metrics[2] / 100.0)

# Open the output file to write the analysis summary
outfile = open(args.dataset_name + "_analysis_summary_reffree.tsv", 'w')

# Write the header line
outfile.write("Tool" + "\t" + "Sensitivity" + "\t" + "Precision" + "\t" + "F1-Score\n")

# Write the data rows
for i in range(len(args.tool_names)):
	outfile.write(args.tool_names[i] + "\t" + str(sensitivity[i]) + "\t" + str(precision[i]) + "\t" + str(f1_score[i]) + "\n")

# Don't forget to properly close the output file
outfile.close()
