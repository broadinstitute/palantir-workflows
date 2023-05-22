import sys
import os
import argparse

# reliable: all methods detected
# almost_reliable: n - 1 methods detected
# mult_pred: multiple (i.e. 2: n - 2) methods detected
# unique: only 1 method predicted
# missed: all others predicted except this one

def process_tracking(input_tracking_file, num_gtfs):
	reliable_transcripts = 0
	almost_reliable_transcripts = [0] * num_gtfs
	unique_transcripts = [0] * num_gtfs
	missed_transcripts = [0] * num_gtfs
	total_transcripts = [0] * num_gtfs

	for line in open(input_tracking_file, "r"):
		# column[0]: unique internal id for the transfrag
		# column[1]: unique internal id for the super-locus containing these transcripts across all samples and the reference annotation
		# column[2]: gene name and transcript id of the reference record associated to this transcript
		# column[3]: type of overlap or relationship between the reference transcripts and the transcript structure represented by this row
		# columns[4:]: each following column showns the transcript for each sample/tool
		vals = line.strip().split()[4:]
		assert len(vals) == num_gtfs

		equality_vector = [0 if vals[i] == '-' else 1 for i in range(num_gtfs)]
		for i, e in enumerate(equality_vector):
			if e == 1:
				total_transcripts[i] += 1

		matches_transcripts = equality_vector.count(1)
		assert matches_transcripts > 0

		if matches_transcripts == num_gtfs:
			reliable_transcripts += 1
		elif matches_transcripts == num_gtfs - 1:
			for i, e in enumerate(equality_vector):
				if e == 1:
					almost_reliable_transcripts[i] += 1
			index = equality_vector.index(0)
			missed_transcripts[index] += 1
		elif matches_transcripts == 1:
			index = equality_vector.index(1)
			unique_transcripts[index] += 1

	return reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts


def print_denovo_stats(tool_names, reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts, out = sys.stdout):
	total_columns = len(total_transcripts)
	out.write("tool_name\t" + "\t".join(tool_names) + "\n")
	out.write("reliable_transcripts\t" + "\t".join([str(reliable_transcripts)] * total_columns) + "\n")
	out.write("almost_reliable_transcripts\t" + "\t".join([str(almost_reliable_transcripts[i]) for i in range(total_columns)]) + "\n")
	out.write("mult_pred\t" + "\t".join([str(total_transcripts[i] - unique_transcripts[i] - almost_reliable_transcripts[i] - reliable_transcripts) for i in range(total_columns)]) + "\n")
	out.write("unique\t" + "\t".join([str(unique_transcripts[i]) for i in range(total_columns)]) + "\n")
	out.write("missed\t" + "\t".join([str(-missed_transcripts[i]) for i in range(total_columns)]) + "\n")
	out.flush()

# Main starts here
parser = argparse.ArgumentParser(description = "Extract and print denovo model statistics for long read isoform reconstruction benchmarking.")
parser.add_argument("-d", "--dataset-name", help = "Name of the dataset.", required = True)
parser.add_argument("-s", "--split-type", help = "Type of the split (full, known, or novel) used to obtain the tracking file.", required = True)
parser.add_argument("-r", "--tracking", help = "Tracking file that accompanies the input GTFs, obtained from a reference-free run of gffcompare with the compared tools.", required = True)
parser.add_argument("-t", "--tool-names", nargs = '+', help = "List of the names of software used for the analysis. Order must match tracking file.", required = True)
parser.add_argument("-n", "--num-tools", type = int, help = "Number of tools (including this one) included in the benchmarking.", required = True)
args = parser.parse_args()

# Open output file
out_stats = open(args.dataset_name + "_" + args.split_type + "_denovo_model_stats.tsv", "w")

# Parse the .tracking file to get the transcript statistics
reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts = process_tracking(args.tracking, args.num_tools)

# Print the transcript statistics
print_denovo_stats(args.tool_names, reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts, out = out_stats)

# Close output file
out_stats.close()
