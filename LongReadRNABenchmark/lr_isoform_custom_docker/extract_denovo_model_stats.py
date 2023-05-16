import sys
import os
import argparse

def process_tracking(input_tracking_file, num_gtfs):
	reliable_transcripts = 0
	almost_reliable_transcripts = [0] * num_gtfs # supported by (num_gtfs - 2) tools
	unique_transcripts = [0] * num_gtfs
	missed_transcripts = [0] * num_gtfs
	total_transcripts = [0] * num_gtfs

	for line in open(input_tracking_file, "r"):
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


def print_denovo_stats(split_type, reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts, out = sys.stdout):
	out.write(split_type + "\n")
	total_columns = len(total_transcripts)
	out.write("\t".join([str(reliable_transcripts)] * total_columns) + "\n")
	out.write("\t".join([str(almost_reliable_transcripts[i]) for i in range(total_columns)]) + "\n")
	out.write("\t".join([str(total_transcripts[i] - unique_transcripts[i] - almost_reliable_transcripts[i] - reliable_transcripts) for i in range(total_columns)]) + "\n")
	out.write("\t".join([str(unique_transcripts[i]) for i in range(total_columns)]) + "\n")
	out.write("\t".join([str(-missed_transcripts[i]) for i in range(total_columns)]) + "\n")
	out.flush()

# Main starts here
parser = argparse.ArgumentParser(description = "Extract and print denovo model statistics for long read isoform reconstruction benchmarking.")
parser.add_argument("-s", "--split-type", help = "Type of the split (full, known, or novel) used to obtain the tracking file.", required = True)
parser.add_argument("-r", "--tracking", help = "Tracking file that accompanies the input GTFs, obtained from a reference-free run of gffcompare with the compared tools.", required = True)
parser.add_argument("-n", "--num-tools", type = int, help = "Number of tools (including this one) included in the benchmarking.", required = True)
args = parser.parse_args()

out_stats = open(args.split_type + ".denovo_model_stats.tsv", "w")
reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts = process_tracking(args.tracking, args.num_tools)
print_denovo_stats(args.split_type, reliable_transcripts, almost_reliable_transcripts, total_transcripts, unique_transcripts, missed_transcripts, out = out_stats)
out_stats.close()
