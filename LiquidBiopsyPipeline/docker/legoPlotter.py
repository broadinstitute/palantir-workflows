## Author: Justin Rhoades <rhoades@broadinstitute.org>
# use Python-2.7

# Imports
from __future__ import division
import argparse
import gzip
import logging
import sys
import warnings
warnings.filterwarnings("ignore") # filters warnings from matplotlib plotting

from Bio import SeqIO
import matplotlib
matplotlib.use("Agg")
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Constants
COMP_DICT = {"A": "T", "T": "A", "G": "C", "C": "G"}

DZ_IDX = ["TG_C_T","TA_C_T","TC_C_T","TT_C_T","TG_C_A","TA_C_A","TC_C_A","TT_C_A","TG_C_G","TA_C_G","TC_C_G","TT_C_G",
          "CG_C_T","CA_C_T","CC_C_T","CT_C_T","CG_C_A","CA_C_A","CC_C_A","CT_C_A","CG_C_G","CA_C_G","CC_C_G","CT_C_G",
          "AG_C_T","AA_C_T","AC_C_T","AT_C_T","AG_C_A","AA_C_A","AC_C_A","AT_C_A","AG_C_G","AA_C_G","AC_C_G","AT_C_G",
          "GG_C_T","GA_C_T","GC_C_T","GT_C_T","GG_C_A","GA_C_A","GC_C_A","GT_C_A","GG_C_G","GA_C_G","GC_C_G","GT_C_G",
          "TG_A_G","TA_A_G","TC_A_G","TT_A_G","TG_A_C","TA_A_C","TC_A_C","TT_A_C","TG_A_T","TA_A_T","TC_A_T","TT_A_T",
          "CG_A_G","CA_A_G","CC_A_G","CT_A_G","CG_A_C","CA_A_C","CC_A_C","CT_A_C","CG_A_T","CA_A_T","CC_A_T","CT_A_T",
          "AG_A_G","AA_A_G","AC_A_G","AT_A_G","AG_A_C","AA_A_C","AC_A_C","AT_A_C","AG_A_T","AA_A_T","AC_A_T","AT_A_T",
          "GG_A_G","GA_A_G","GC_A_G","GT_A_G","GG_A_C","GA_A_C","GC_A_C","GT_A_C","GG_A_T","GA_A_T","GC_A_T","GT_A_T"]

DZ_IDX_SET = set(DZ_IDX)

CONT_LIST = ["T_G", "T_A", "T_C", "T_T",
             "C_G", "C_A", "C_C", "C_T",
             "A_G", "A_A", "A_C", "A_T",
             "G_G", "G_A", "G_C", "G_T"]

FILE_COL_NAMES = {"maf": ["ref_context", "Reference_Allele", "Tumor_Seq_Allele2", "i_tumor_f", "i_judgement"],
                  "call_stats": ["context", "ref_allele", "alt_allele", "tumor_f", "judgement"],
                  "vcf": ["context", "ref_allele", "alt_allele", "tumor_f", "judgement"]}


# Classes
class VCFReader(object):
    '''Reader for VCF files.'''
    
    def __init__(self, filename):
        gzipped = filename.endswith(".gz")
        read_mode = "rb" if gzipped else "rt"
        self.__reader = open(filename, read_mode)
        if gzipped:
            self.__reader = gzip.GzipFile(fileobj=self.__reader)
            
        line = next(self.__reader)
        while line.startswith("##"):
            line = next(self.__reader)

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.__reader)
        line = line.strip().split("\t")
        
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
        
        record = VCFRecord(chrom, pos, ident, ref, alt, qual, filt)
        
        return record
        
    next = __next__  ## Python 2/3 compatibility


class VCFRecord(object):
    '''Record from a VCF file'''
    
    def __init__(self, chrom, pos, ident, ref, alt, qual, filt):
        self.CHROM = chrom
        self.POS = pos
        self.ID = ident
        self.REF = ref
        self.ALT = alt
        self.QUAL = qual
        self.FILTER = filt


# Functions
def parse_cl_args():
    parser = argparse.ArgumentParser(description="Create lego plots.")

    parser.add_argument("-m", "--muts", help="mutation file", required=True)
    parser.add_argument("-t", "--fileType", choices=["call_stats", "maf", "vcf"], help="mutation file type", required=True)
    parser.add_argument("-r", "--reference", help="reference genome fasta")
    parser.add_argument("-b", "--basename", help="output base name", default="lego_plot")
    parser.add_argument("--afCut", help="only plot mutations <= allele fraction", type=float, default=1.0)
    mut_excl_grp = parser.add_mutually_exclusive_group()
    mut_excl_grp.add_argument("--passOnly", action="store_true", help="PASS mutations only")
    mut_excl_grp.add_argument("--filterOnly", action="store_true", help="Filtered mutations only")
    parser.add_argument("--logging", choices=["INFO", "DEBUG"], help="set logging level", default="DEBUG")

    return parser.parse_args()


def gen_lego_plot(dz, basename):
    n_muts = sum(dz)
    fig = plt.figure(figsize=(13,8), dpi=300)
    ax = fig.add_subplot(111, projection="3d")

    xpos = []
    ypos = []
    for j in range(1,9):
        for i in range(1,13):
            xpos.append(i)
            ypos.append(j)

    for j in range(9,13):
        for i in range(3,7):
            xpos.append(i)
            ypos.append(j)

    df = pd.DataFrame({"xpos": xpos, "ypos": ypos, "barCol": ["white"]*112})
    df.ix[(df.xpos <= 4) & (df.ypos <= 4), "barCol"] = "#FEFD38" ##yellow
    df.ix[(df.xpos > 4) & (df.xpos <= 8) & (df.ypos <= 4), "barCol"] = "#1CB2B1" ##aqua
    df.ix[(df.xpos > 8) & (df.ypos <= 4), "barCol"] = "#FC0F1C" ##red
    df.ix[(df.xpos <= 4) & (df.ypos > 4) & (df.ypos <= 8), "barCol"] = "#2BCA2E" ##green
    df.ix[(df.xpos > 4) & (df.xpos <= 8) & (df.ypos > 4) & (df.ypos <= 8), "barCol"] = "#0C3BC8" ##blue
    df.ix[(df.xpos > 8) & (df.ypos > 4) & (df.ypos <= 8), "barCol"] = "#7D4DAC" ##purple

    zpos = np.zeros(112)

    dx = np.repeat(0.8, 112)
    dy = np.repeat(0.8, 112)
    dz = dz + list(np.zeros(16))

    ax.view_init(elev=50, azim=55)
    ax.set_zlim(bottom=0, top=max([5,max(dz)]))
    ax.set_xlim(left=1, right=12.5)
    ax.set_ylim(bottom=1, top=12)

    for idx in range(len(xpos)):
        ax.bar3d(xpos[idx], ypos[idx], zpos[idx], dx[idx], dy[idx], dz[idx], linewidth=0.3, color=df.barCol.tolist()[idx], zsort="max")
        ax.collections[idx].set_sort_zpos(idx)
        ax.collections[idx].set_facecolors([df.barCol.tolist()[idx]]*6)  

    # Get rid of the panes
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # Get rid of the spines
    ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    # Get rid of the ticks
    ax.set_xticks([]) 
    ax.set_yticks([])

    # Add Legend
    yellow_proxy = plt.Rectangle((0,0), 1, 1, fc="#FEFD38")
    aqua_proxy = plt.Rectangle((0,0), 1, 1, fc="#1CB2B1")
    red_proxy = plt.Rectangle((0,0), 1, 1, fc="#FC0F1C")
    green_proxy = plt.Rectangle((0,0), 1, 1, fc="#2BCA2E")
    blue_proxy = plt.Rectangle((0,0), 1, 1, fc="#0C3BC8")
    purple_proxy = plt.Rectangle((0,0), 1, 1, fc="#7D4DAC")

    ax.legend([yellow_proxy,aqua_proxy,red_proxy,green_proxy,blue_proxy,purple_proxy],
              ["C>T","C>A","C>G","A>G","A>C","A>T"], loc=(0.8,0.6))

    ax.set_zlabel('SNV Count')

    # Add context text
    for idx, pos in enumerate(range(96,112)):
        ax.text(xpos[pos]+0.45, ypos[pos]+0.45, 0, CONT_LIST[idx], size="small", ha="center", va="center", zorder=112)  

    # Add title and n text
    ax.text2D(0.15, 0.9, basename, transform=ax.transAxes)
    ax.text2D(0.15, 0.87, "n = {0}".format(n_muts), transform=ax.transAxes)

    fig.savefig(basename+".pdf", bbox_inches="tight")


def gen_dz_arr_from_maf(muts_df, cc, rc, ac):
    dz_arr = pd.Series(data=np.zeros(96), index=DZ_IDX, dtype=np.int64)
    muts_df["cont_ext"] = muts_df[cc].apply(extract_context)

    for idx, mut in muts_df.iterrows():
        if not isProperSNV(mut[rc], mut[ac]):
            continue
        cont_mut = mut["cont_ext"] + "_" + mut[rc].upper() + "_" + mut[ac].upper()
        if cont_mut not in DZ_IDX_SET:
            cont_mut = flip_cont_mut(cont_mut)
            if cont_mut not in DZ_IDX_SET:
                logging.error("Cannot find context for {0}".format(cont_mut))
                continue

        dz_arr[cont_mut] += 1

    return dz_arr.tolist()


def extract_context(raw_cont_str):
    raw_cont_str = raw_cont_str.upper()
    cont_len = len(raw_cont_str)
    cont_len_mid = cont_len // 2

    return raw_cont_str[cont_len_mid-1] + raw_cont_str[cont_len_mid+1]


def flip_cont_mut(cont_mut):
    cont_mut = cont_mut.split("_")
    cont_mut[0] = COMP_DICT[cont_mut[0][1]] + COMP_DICT[cont_mut[0][0]]
    cont_mut[1] = COMP_DICT[cont_mut[1]]
    cont_mut[2] = COMP_DICT[cont_mut[2]]

    return "_".join(cont_mut)


def vcf_to_muts_df(muts_file, ref_fasta_file, ref_col, alt_col, context_col, af_col, judge_col):
    logging.info("Reading reference fasta into memory.")
    ref_fasta = SeqIO.to_dict(SeqIO.parse(ref_fasta_file, "fasta"))
    logging.info("Done reading reference.")
    
    in_vcf = VCFReader(muts_file)
    
    muts_mat = []
    for mut in in_vcf:
        for alt in mut.ALT:
            muts_mat.append([str(mut.REF).upper(), str(alt).upper(), str(ref_fasta[mut.CHROM].seq[mut.POS-2:mut.POS+1]), 1.0, ",".join(list(mut.FILTER))])
    
    if not muts_mat:
        logging.warning("No mutations found.")
        muts_mat = None
        
    muts_df = pd.DataFrame(muts_mat, columns=[ref_col, alt_col, context_col, af_col, judge_col])
    
    return muts_df


def isProperSNV(ref, alt):
    return len(ref) == 1 and len(alt) == 1 and ref.upper() in COMP_DICT and alt.upper() in COMP_DICT


# Main
def main():
    args = parse_cl_args()

    muts_file = args.muts
    file_type = args.fileType
    ref_fasta_file = args.reference
    basename = args.basename
    af_cut = args.afCut
    pass_only = args.passOnly
    filter_only = args.filterOnly
    log_level = args.logging
    
    logging.basicConfig(format='%(levelname)s\t%(asctime)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=log_level)
    logging.info("Starting lego plotter.")    
    logging.debug("Input file type: {0}".format(file_type))
    if file_type == "vcf" and af_cut != 1.0:
        logging.warning("AF cutoff not yet supported for VCF input. Using cutoff of 1.0.")
        af_cut = 1.0 
    logging.debug("Allele fraction cutoff: {0}, Pass Only: {1}, Filter Only: {2}".format(af_cut, pass_only, filter_only))

    context_col, ref_col, alt_col, af_col, judge_col = FILE_COL_NAMES[file_type]

    # read in mutation file to data frame
    if file_type == "vcf":
        if not ref_fasta_file:
            logging.critical("Reference fasta required for VCF input.")
            sys.exit(1)
        muts_df = vcf_to_muts_df(muts_file, ref_fasta_file, ref_col, alt_col, context_col, af_col, judge_col)
    else:
        try:
            muts_df = pd.read_table(muts_file, comment="#")
        except UnicodeDecodeError:
            muts_df = pd.read_table(muts_file, comment="#", encoding="latin-1")

    logging.debug("Total mutations: {0}".format(muts_df.shape[0]))
    
    if pass_only:
        # remove mutations that were filtered out
        if file_type == "vcf":
            muts_df = muts_df.ix[muts_df[judge_col] == "PASS",]
        else:
            muts_df = muts_df.ix[muts_df[judge_col] == "KEEP",]
    
    if filter_only:
        # remove mutations marked as KEEP/PASS
        if file_type == "vcf":
            muts_df = muts_df.ix[muts_df[judge_col] != "PASS",]
        else:
            muts_df = muts_df.ix[muts_df[judge_col] == "REJECT",]

    if af_cut < 1.0:
        # filter mutations above allele fraction cutoff
        muts_df = muts_df.ix[muts_df[af_col] <= af_cut,]

    dz = gen_dz_arr_from_maf(muts_df, context_col, ref_col, alt_col)

    logging.debug("SNVs after filtering: {0}".format(sum(dz)))
    logging.info("Plotting {0} mutations.".format(sum(dz)))
    
    gen_lego_plot(dz, basename)
    
    logging.info("Plotting to {0}.pdf complete.".format(basename))


if __name__ == "__main__":
    main()
