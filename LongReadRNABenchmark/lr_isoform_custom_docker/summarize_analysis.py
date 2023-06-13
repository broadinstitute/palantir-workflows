import sys
import os
import argparse
import math

def process_tracking(input_tracking_file):
	novel_tp = 0
	novel_fp = 0
	novel_fn = 0
	known_tp = 0
	known_fn = 0
	
	for line in open(input_tracking_file, "r"):
		# column[0]: unique internal id for the transfrag
		# column[1]: unique internal id for the super-locus containing these transcripts across all samples and the reference annotation
		# column[2]: gene name and transcript id of the reference record associated to this transcript
		# column[3]: type of overlap or relationship between the reference transcripts and the transcript structure represented by this row
		# columns[4:]: each following column showns the transcript for each sample/tool
		transcript_columns = line.strip().split()

		if transcript_columns[4] != '-' and transcript_columns[5] == '-' and transcript_columns[6] == '-':
			novel_fn += 1
		elif transcript_columns[4] != '-' and transcript_columns[5] != '-' and transcript_columns[6] == '-':
			known_fn += 1
		elif transcript_columns[4] != '-' and transcript_columns[5] != '-' and transcript_columns[6] != '-':
			known_tp += 1
		elif transcript_columns[4] != '-' and transcript_columns[5] == '-' and transcript_columns[6] != '-':
			novel_tp += 1
		elif transcript_columns[4] == '-' and transcript_columns[5] == '-' and transcript_columns[6] != '-':
			novel_fp += 1
		else:
			print("WARNING: This should not have happened! Current line: " + str(transcript_columns))

	return (novel_tp, novel_fp, novel_fn, known_tp, known_fn)

def compute_novel_accuracy_stats(novel_tp, novel_fp, novel_fn):
	sensitivity = 0.0
	precision = 0.0
	f1_score = 0.0

	if not math.isclose(novel_tp + novel_fn, 0.0):
		sensitivity = novel_tp / (novel_tp + novel_fn)

	if not math.isclose(novel_tp + novel_fp, 0.0):
		precision = novel_tp / (novel_tp + novel_fp)

	if not math.isclose(precision + sensitivity, 0.0):
		f1_score = 2 * precision * sensitivity / (precision + sensitivity)

	return (sensitivity, precision, f1_score)

# Main starts here
parser = argparse.ArgumentParser(description = "Report accuracy statistics for novel isoform detection from long reads.")
parser.add_argument("-r", "--tracking", nargs = '+', help = "List of <outprefix>.tracking files for the tools being benchmarked.", required = True)
parser.add_argument("-t", "--tool-names", nargs = '+', help = "List of tool names. MUST be in the same order as their corresponding tracking files.", required = True)
parser.add_argument("-d", "--dataset-name", help = "Name of the dataset.", required = True)
args = parser.parse_args()

# Only the length match is checked
# It is the user's responsibility to make sure that the order of the input tracking files match the order of the tool names
if len(args.tracking) != len(args.tool_names):
	print("Length of the input tracking file list must be equal to the length of tool name list. Order MUST also be the same, but this is NOT checked by the program!")
	sys.exit(0)

# Parse the .tracking file to get the tracking statistics and use the output to compute the novel detection accuracy statistics for each tool
novel_tps = []
novel_fps = []
novel_fns = []
known_tps = []
known_fns = []
sensitivity = []
precision = []
f1_score = []

for i in range(len(args.tracking)):
	tracking_stats = process_tracking(args.tracking[i])
	accuracy_stats = compute_novel_accuracy_stats(tracking_stats[0], tracking_stats[1], tracking_stats[2])

	novel_tps.append(tracking_stats[0])
	novel_fps.append(tracking_stats[1])
	novel_fns.append(tracking_stats[2])
	known_tps.append(tracking_stats[3])
	known_fns.append(tracking_stats[4])

	sensitivity.append(accuracy_stats[0])
	precision.append(accuracy_stats[1])
	f1_score.append(accuracy_stats[2])

# Open output file
outfile = open(args.dataset_name + "_analysis_summary.tsv", "w")

# Write the header row
outfile.write("Tool" + "\t" + "Sensitivity(Novel)" + "\t" + "Precision(Novel)" + "\t" + "F1-Score(Novel)" + "\t")
outfile.write("TP_Novel" + "\t" + "FP_Novel" + "\t" + "FN_Novel" + "\t" + "TP_Known" + "\t" + "FN_Known" + "\n")

# Write the data rows
for i in range(len(args.tool_names)):
	outfile.write(args.tool_names[i] + "\t" + str(sensitivity[i]) + "\t" + str(precision[i]) + "\t" + str(f1_score[i]) + "\t")
	outfile.write(str(novel_tps[i]) + "\t" + str(novel_fps[i]) + "\t" + str(novel_fns[i]) + "\t" + str(known_tps[i]) + "\t" + str(known_fns[i]) + "\n")

# Close output file
outfile.close()
