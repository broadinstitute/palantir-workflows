import os
import argparse
from collections import defaultdict

def split_isoquant(line):
	if line.find(".nic") != -1 or line.find(".nnic") != -1:
		return "novel"
	elif line.find('transcript_id "SIRV') != -1 or line.find('transcript_id "ENS') != -1:
		return "known"
	return "undefined"

def split_stringtie(line):
	if line.find("reference_id") != -1:
		return "known"
	else:
		return "novel"

def split_transcript_id(line):
	if line.find('transcript_id "SIRV') != -1: # for SIRVs
		return "known"
	elif line.find('transcript_id "ENS') != -1 and line.find("aligned_") == -1 and line.find("PB.") == -1: # for simulated data
		return "known"
	else:
		return "novel"

def initialize_bambu_count_dict(input_counts):
	bambu_count_dict = defaultdict(float)

	for line in open(input_counts):
		if line.startswith("#") or line.startswith("TXNAME"):
			continue
		t = line.strip().split()
		tid = t[0]
		bambu_count_dict[tid] = max(bambu_count_dict[tid], float(t[2]))

	return bambu_count_dict

def split_bambu(line, bambu_count_dict):
	tpos = line.find('transcript_id')
	if tpos == -1:
		return "undefined"

	idpos = tpos + len('transcript_id') + 2
	endpos = line.find(";", idpos)
	if endpos == -1:
		print("Warning, unable to find ;")
		return "undefined"
	tid = line[idpos:endpos-1]

	if tid not in bambu_count_dict or bambu_count_dict[tid] == 0:
		return "undefined"
	elif tid.startswith('ENS') or tid.startswith('SIRV'):
		return "known"
	else:
		return "novel"

def split_flames(line):
	if line.split('\t')[1] == "known":
		return "known"
	else:
		return "novel"

def split_gtf(input_gtf, tool, out_path_full, out_path_known, out_path_novel, bambu_count_dict):
	out_full = open(out_path_full, "w")
	out_known = open(out_path_known, "w")
	out_novel = open(out_path_novel, "w")

	for line in open(input_gtf):
		if line.startswith("#"):
			continue

		transcript_type = ""
		if tool == "isoquant":
			transcript_type = split_isoquant(line)
		elif tool == "stringtie":
			transcript_type = split_stringtie(line)
		elif tool == "flair":
			transcript_type = split_transcript_id(line)
		elif tool == "talon":
			transcript_type = split_transcript_id(line)
		elif tool == "bambu":
			transcript_type = split_bambu(line, bambu_count_dict)
		elif tool == "flames":
			transcript_type = split_flames(line)

		if transcript_type == "undefined":
			continue

		out_full.write(line)

		if transcript_type == "novel":
			out_novel.write(line)
		elif transcript_type == "known":
			out_known.write(line)

	out_full.close()
	out_novel.close()
	out_known.close()

# Main Program
parser = argparse.ArgumentParser(description = "Split the reconstructed isoform output GTF from tools such as IsoQuant and StringTie in to known and novel based on annotations.")
parser.add_argument("-g", "--input-gtf", help = "Input GTF file that contains the reconstructed isoforms to be split into full/known/novel.", required = True)
parser.add_argument("-t", "--tool", choices = ["isoquant", "stringtie", "flair", "talon", "bambu", "flames"], help = "Name of the tool used to obtain the input GTF.", required = True)
parser.add_argument("-c", "--input-bambu-counts", help = "Input .gtf.counts file, which is required only if the tool used for isoform reconstruction is Bambu.", required = False)
args = parser.parse_args()

if args.tool == "bambu" and (args.input_bambu_counts is None):
	parser.error("--tool bambu requires --input-bambu-counts.")

bambu_count_dict = defaultdict(float)
if args.tool == "bambu":
	bambu_count_dict = initialize_bambu_count_dict(args.input_bambu_counts)

base = os.path.basename(args.input_gtf)
root_ext = os.path.splitext(base)
out_path_full = root_ext[0] + ".full.gtf"
out_path_known = root_ext[0] + ".known.gtf"
out_path_novel = root_ext[0] + ".novel.gtf"

split_gtf(args.input_gtf, args.tool, out_path_full, out_path_known, out_path_novel, bambu_count_dict)
