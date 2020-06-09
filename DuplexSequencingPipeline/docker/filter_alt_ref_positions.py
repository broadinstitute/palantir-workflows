import argparse
import gzip
import os
import logging
import subprocess
import sys
import tempfile

import pandas as pd

logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    format="%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("FilterAltRefPositions")

BLASTN_CMD = ('{0} -query {1} -db {2} -perc_identity 90 -qcov_hsp_perc '
              '62.5 -dust no -soft_masking false -outfmt')

BLASTN_FIELDS = "6 std sseq sstrand"

class VCFReader(object):
    """Reader for VCF files."""

    def __init__(self, filename):
        gzipped = filename.endswith(".gz")
        if gzipped:
            self.__reader = gzip.open(filename, mode="rt")
        else:
            self.__reader = open(filename, "rt")

        self.header = ""

        line = next(self.__reader)
        while line.startswith("##"):
            self.header += line
            line = next(self.__reader)

        self.header += line


    def __iter__(self):
        return self

    def __next__(self):
        raw_line = next(self.__reader)
        line = raw_line.strip().split("\t")

        chrom = line[0]
        pos = int(line[1])

        if line[2] == ".":
            ident = None
        else:
            ident = line[2]

        ref = line[3]
        alt = line[4].split(",")

        try:
            qual = int(line[5])
        except ValueError:
            try:
                qual = float(line[5])
            except ValueError:
                qual = None

        if line[6] == ".":
            filt = []
        else:
            filt = line[6].split(";")

        record = VCFRecord(chrom, pos, ident, ref, alt, qual, filt, raw_line)

        return record

    next = __next__  ## Python 2/3 compatibility


class VCFRecord(object):
    """Record from a VCF file"""

    def __init__(self, chrom, pos, ident, ref, alt, qual, filt, raw_line):
        self.CHROM = chrom
        self.POS = pos
        self.ID = ident
        self.REF = ref
        self.ALT = alt
        self.QUAL = qual
        self.FILTER = filt
        self.RAW = raw_line

    def __str__(self):
        self.ID = "." if not self.ID else self.ID
        self.QUAL = "." if not self.QUAL else self.QUAL
        self.FILTER = ["."] if not self.FILTER else self.FILTER

        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(
            self.CHROM, self.POS, self.ID, self.REF, ",".join(self.ALT),
            self.QUAL, ";".join(self.FILTER))


class FastaReader(object):
    """Reference fasta reader"""

    def __init__(self, fasta, fasta_idx=None):
        if not fasta_idx:
            fasta_idx = fasta + ".fai"

        if not os.path.exists(fasta_idx):
            raise FileNotFoundError("Fasta file must be indexed (.fasta.fai)")

        self.__fasta_idx_df = pd.read_table(
            fasta_idx, header=None, index_col="contig", dtype={"contig": str},
            names=["contig", "length", "start",
                   "bases_per_line", "bytes_per_line"])
        self.__fasta_fh = open(fasta, "r")

    def get(self, contig, seq_start, seq_end):
        bases_per_line = self.__fasta_idx_df.ix[contig, "bases_per_line"]
        bytes_per_line = self.__fasta_idx_df.ix[contig, "bytes_per_line"]
        seq_length = seq_end - seq_start
        contig_start = self.__fasta_idx_df.ix[contig, "start"]
        num_newlines = seq_start // bases_per_line

        adj_seq_start = contig_start + seq_start + num_newlines
        max_seq_newlines = (seq_length // bytes_per_line) + 2
        max_seq_length = seq_length + max_seq_newlines

        self.__fasta_fh.seek(adj_seq_start)
        raw_seq = self.__fasta_fh.read(max_seq_length)
        raw_seq = raw_seq.replace("\n", "")

        seq = raw_seq[:seq_length]

        return seq

    def close(self):
        self.__fasta_fh.close()


def parse_cl_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf", help="VCF file", required=True)
    parser.add_argument("-o", "--outfile", help="output filepath", required=True)
    parser.add_argument("-r", "--reference_fasta", help="reference fasta file", required=True)
    parser.add_argument("--blastn", help="path to blastn executable", default="blastn")

    return parser.parse_args()


def check_alt_ref_positions(contig, position, alt_allele, reference, blastn):
    ref_fasta_reader = FastaReader(reference)

    # check for homologous regions that have mutant
    tmp_blast_seq_chrom = str(contig)
    tmp_blast_seq_start = int(position) - 71
    tmp_blast_seq_end = tmp_blast_seq_start + 80

    tmp_blast_seq_dict = {}
    for i in range(1, 8):
        tmp_blast_seq = ref_fasta_reader.get(tmp_blast_seq_chrom,
                                             tmp_blast_seq_start,
                                             tmp_blast_seq_end)

        # switch reference base for alt allele at mutation site
        tmp_blast_seq_list = list(tmp_blast_seq)
        tmp_blast_seq_list[(8-i) * 10] = str(alt_allele).upper()
        tmp_blast_seq = "".join(tmp_blast_seq_list)
        
        tmp_blast_seq_dict[((8-i) * 10) + 1] = tmp_blast_seq

        tmp_blast_seq_start += 10
        tmp_blast_seq_end += 10

    blast_df = blast_nucs(blastn, tmp_blast_seq_dict, reference)

    specificity_flag = 1
    for idx, hit in blast_df.iterrows():
        seq_pos = int(hit["qseqid"])
        if int(hit["qstart"]) <= seq_pos <= int(hit["qend"]):
            alt_dist = seq_pos - int(hit["qstart"])
            ref_base = hit["sseq"][alt_dist].upper()

            if ref_base == str(alt_allele).upper():
                log.info("found potential alt site in reference")
                specificity_flag = 0
                break

    ref_fasta_reader.close()

    return specificity_flag


def blast_nucs(blastn, seq_dict, blast_db, tmp_dir="./"):
    #if not shutil.which("blastn"):
    #    log.critical("BLAST+ software required. Cannot find blastn command.")
    #    sys.exit(1)

    tmp_query_fh = tempfile.NamedTemporaryFile(
        "w", suffix=".fasta", prefix="blast_query_", dir=tmp_dir)
    write_fasta(seq_dict, tmp_query_fh)
    tmp_query_fh.flush()
    os.fsync(tmp_query_fh.fileno())

    try:
        blast_out_str = subprocess.check_output(
            (BLASTN_CMD.format(blastn, tmp_query_fh.name, blast_db).split() 
             + [BLASTN_FIELDS]))
    except Exception as err:
        log.critical(err)
        sys.exit(1)
    finally:
        tmp_query_fh.close()

    blast_results_mat = blast_out_str.strip().split("\n")
    blast_results_mat = [line.split("\t") for line in blast_results_mat]
    
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
               "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
               "sseq", "sstrand"]
    blast_results_df = pd.DataFrame(blast_results_mat, columns=columns)

    return blast_results_df


def write_fasta(seq_dict, output_file, line_length=80):
    if isinstance(output_file, str):
        outfile = open(output_file, "w")
    else:
        outfile = output_file

    for seq_id in seq_dict:
        outfile.write(">{0}\n".format(seq_id))
        seq = seq_dict[seq_id]
        
        # split seq into lines
        seq_lst = [seq[i:i+line_length]
                   for i in range(0, len(seq), line_length)]

        outfile.write("\n".join(seq_lst) + "\n")

    if isinstance(output_file, str):
        outfile.close()


def main():
    args = parse_cl_args()

    params = vars(args).copy()
    params = ["{0}={1}".format(param, params[param]) for param in params]
    log.info("Params: [%s]", ", ".join(params))

    vcf_path = args.vcf
    outfile_name = args.outfile
    ref_fasta_path = args.reference_fasta
    blastn = args.blastn

    vcf_reader = VCFReader(vcf_path)

    if outfile_name.endswith(".gz"):
        outfile = gzip.open(outfile_name, "wt")
    else:
        outfile = open(outfile_name, "w")
    outfile.write(vcf_reader.header)

    num_filtered = 0
    for mut in vcf_reader:
        log.info("Checking {0}_{1}_{2}_{3}".format(mut.CHROM, mut.POS, mut.REF, mut.ALT[0]))
        spec_flag = check_alt_ref_positions(mut.CHROM, mut.POS, mut.ALT[0], 
                                            ref_fasta_path, blastn)

        if spec_flag == 1:
            outfile.write(mut.RAW)
        else:
            num_filtered += 1

    log.info("Number of mutations filtered: {0}".format(num_filtered))

    outfile.close()


if __name__ == "__main__":
    main()