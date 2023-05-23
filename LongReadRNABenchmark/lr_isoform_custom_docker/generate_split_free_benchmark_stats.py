import sys
import os
import argparse

def process_tracking(input_tracking_file):
	novel_fn = 0
	known_fn = 0
	known_tp = 0
	novel_tp = 0
	novel_fp = 0
	
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

	return novel_fn, known_fn, known_tp, novel_tp, novel_fp

def print_stats(novel_fn, known_fn, known_tp, novel_tp, novel_fp, out = sys.stdout):
	out.write("novel_fn" + "\t" + str(novel_fn) + "\n")
	out.write("known_fn" + "\t" + str(known_fn) + "\n")
	out.write("known_tp" + "\t" + str(known_tp) + "\n")
	out.write("novel_tp" + "\t" + str(novel_tp) + "\n")
	out.write("novel_fp" + "\t" + str(novel_fp) + "\n")
	out.flush()

# Main starts here
parser = argparse.ArgumentParser(description = "Extract accuracy statistics for long read isoform reconstruction benchmarking.")
parser.add_argument("-r", "--tracking", help = "Tracking file that accompanies the input GTFs.", required = True)
parser.add_argument("-t", "--tool-name", help = "Name of software used for the analysis.", required = True)
parser.add_argument("-d", "--dataset-name", help = "Name of the dataset.", required = True)
args = parser.parse_args()

# Open output file
out_stats = open(args.tool_name + "_" + args.dataset_name + "_accuracy_stats.tsv", "w")

# Parse the .tracking file to get the transcript statistics
novel_fn, known_fn, known_tp, novel_tp, novel_fp = process_tracking(args.tracking)

# Print the transcript statistics
print_stats(novel_fn, known_fn, known_tp, novel_tp, novel_fp, out = out_stats)

# Close output file
out_stats.close()
