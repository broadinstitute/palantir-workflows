import sys
import os
import argparse

# reliable: all methods detected
# almost_reliable: n - 1 methods detected
# mult_pred: multiple (i.e. 2: n - 2) methods detected
# unique: only 1 method predicted
# missed: all others predicted except this one

def process_tracking(input_tracking_file, num_gtfs):
	reliable_transcripts = {}
	almost_reliable_transcripts = {}
	unique_transcripts = {}
	missed_transcripts = {}
	total_transcripts = {}

	reliable_transcripts["novel"] = 0
	reliable_transcripts["known"] = 0
	almost_reliable_transcripts["novel"] = [0] * num_gtfs
	almost_reliable_transcripts["known"] = [0] * num_gtfs
	unique_transcripts["novel"] = [0] * num_gtfs
	unique_transcripts["known"] = [0] * num_gtfs
	missed_transcripts["novel"] = [0] * num_gtfs
	missed_transcripts["known"] = [0] * num_gtfs
	total_transcripts["novel"] = [0] * num_gtfs
	total_transcripts["known"] = [0] * num_gtfs

	for line in open(input_tracking_file, "r"):
		# column[0]: unique internal id for the transfrag
		# column[1]: unique internal id for the super-locus containing these transcripts across all samples and the reference annotation
		# column[2]: gene name and transcript id of the reference record associated to this transcript
		# column[3]: type of overlap or relationship between the reference transcripts and the transcript structure represented by this row
		# columns[4:]: each following column showns the transcript for each sample/tool
		# In our case, however, column[4] is used to determine known/novel so tool outputs are reflected in columns [5:]
		prediction_type = ""
		expressed_kept = line.strip().split()[4]
		if expressed_kept == '-':
			prediction_type = "novel"
		else:
			prediction_type = "known"

		vals = line.strip().split()[5:]
		assert len(vals) == num_gtfs

		equality_vector = [0 if vals[i] == '-' else 1 for i in range(num_gtfs)]
		for i, e in enumerate(equality_vector):
			if e == 1:
				total_transcripts[prediction_type][i] += 1

		matches_transcripts = equality_vector.count(1)

		if matches_transcripts == num_gtfs:
			reliable_transcripts[prediction_type] += 1
		elif matches_transcripts == num_gtfs - 1:
			for i, e in enumerate(equality_vector):
				if e == 1:
					almost_reliable_transcripts[prediction_type][i] += 1
			index = equality_vector.index(0)
			missed_transcripts[prediction_type][index] += 1
		elif matches_transcripts == 1:
			index = equality_vector.index(1)
			unique_transcripts[prediction_type][index] += 1

	return reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts


def print_denovo_stats(tool_names, reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts, out = sys.stdout):
	total_columns = len(total_transcripts)
	out.write("Tool_Names\t" + "\t".join(tool_names) + "\n")
	out.write("Reliable_Transcripts\t" + "\t".join([str(reliable_transcripts)] * total_columns) + "\n")
	out.write("Almost_Reliable_Transcripts\t" + "\t".join([str(almost_reliable_transcripts[i]) for i in range(total_columns)]) + "\n")
	out.write("Mult_Pred\t" + "\t".join([str(total_transcripts[i] - unique_transcripts[i] - almost_reliable_transcripts[i] - reliable_transcripts) for i in range(total_columns)]) + "\n")
	out.write("Unique\t" + "\t".join([str(unique_transcripts[i]) for i in range(total_columns)]) + "\n")
	out.write("Missed\t" + "\t".join([str(-missed_transcripts[i]) for i in range(total_columns)]) + "\n")
	out.flush()

# Main starts here
parser = argparse.ArgumentParser(description = "Extract and print denovo model statistics for long read isoform reconstruction benchmarking.")
parser.add_argument("-d", "--dataset-name", help = "Name of the dataset.", required = True)
parser.add_argument("-r", "--tracking", help = "Tracking file that accompanies the input GTFs, obtained from a reference-free run of gffcompare with the compared tools.", required = True)
parser.add_argument("-t", "--tool-names", nargs = '+', help = "List of the names of software used for the analysis. Order must match tracking file.", required = True)
parser.add_argument("-n", "--num-tools", type = int, help = "Number of tools (including this one) included in the benchmarking.", required = True)
args = parser.parse_args()

# Open output files
outfile_novel = open(args.dataset_name + "_denovo_analysis_summary_novel.tsv", "w")
outfile_known = open(args.dataset_name + "_denovo_analysis_summary_known.tsv", "w")

# Parse the .tracking file to get the transcript statistics
reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts = process_tracking(args.tracking, args.num_tools)

# Print the transcript statistics for novel/known
print_denovo_stats(args.tool_names, reliable_transcripts["novel"], almost_reliable_transcripts["novel"], total_transcripts["novel"], unique_transcripts["novel"], missed_transcripts["novel"], out = outfile_novel)
print_denovo_stats(args.tool_names, reliable_transcripts["known"], almost_reliable_transcripts["known"], total_transcripts["known"], unique_transcripts["known"], missed_transcripts["known"], out = outfile_known)

# Close output files
outfile_novel.close()
outfile_known.close()
