truth_path = 'test/CNVControlEventsQC/test_data/inputs/truth.tsv'
intervals_path = 'test/CNVControlEventsQC/test_data/inputs/intervals.tsv'
path = 'test/CNVControlEventsQC/test_data/inputs/eval.vcf'


import pandas as pd
import numpy as np

def add_exon_idxs(df, exons):
        contigs = set(df.contig)
        for contig in contigs:
            df.loc[df.contig==contig,"start_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].end,
                                                                       df.loc[df.contig==contig].start,"left")
            df.loc[df.contig==contig,"end_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].start,
                                                                       df.loc[df.contig==contig].end,"right")
def get_exon_expanded_events(df, exons):
    add_exon_idxs(df, exons)
    df = df.loc[df.start_exon_idx != df.end_exon_idx].reset_index().astype({'start_exon_idx':int,'end_exon_idx':int})
    df_expanded = df.loc[df.index.repeat(df.end_exon_idx-df.start_exon_idx)]
    df_expanded['exon_idx'] = df_expanded.groupby(df_expanded.index).cumcount() + df_expanded.start_exon_idx
    df_expanded = df_expanded.set_index(df_expanded.contig + "_" + df_expanded.exon_idx.astype(str))
    df_expanded = df_expanded.join(exons[['start','end']], rsuffix='_exon')
    df_expanded['event_exon_start']=np.maximum(df_expanded.start, df_expanded.start_exon)
    df_expanded['event_exon_end']=np.minimum(df_expanded.end, df_expanded.end_exon)
    df_expanded = df_expanded.set_index(df_expanded.index + "_" + df_expanded.svtype)
    return df_expanded

def extract_end_from_info(info):
    infos = info.split(';')
    for i in infos:
        if i.startswith('END='):
            return i.replace('END=', '')    

def read_vcf_to_df(vcf_path):
    df = pd.read_csv(vcf_path, sep='\t',comment='#',
                         names=['contig','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'], dtype={'start': int})
    
    df = df.query("ALT != '.' and FILTER=='PASS'").copy()
    df['end'] = df['INFO'].apply(extract_end_from_info).astype(int)
    df["svtype"] = df.ALT.str.replace("<","").str.replace(">","")
    return df[['contig', 'start', 'end', 'ALT', 'svtype', 'ID']]

def read_intervals_to_df(intervals_path):
    intervals = pd.read_csv(intervals_path, sep="\t", comment="@", names = ["contig","start","end","dummy1","dummy2"], dtype={'start': int, 'end': int})
    intervals['contig_idx'] = intervals.groupby('contig').cumcount()
    intervals = intervals.set_index(intervals.contig + "_" + intervals.contig_idx.astype(str))
    return intervals

def read_truth(truth_path):
    truth = pd.read_csv(truth_path, sep="\t")
    truth = truth.reset_index(names='event_index')
    return truth
intervals = read_intervals_to_df(intervals_path)
truth = read_truth(truth_path)
eval_df = read_vcf_to_df(path)
eval_df_expanded = get_exon_expanded_events(eval_df, intervals)

truth_expanded = get_exon_expanded_events(truth, intervals)

truth_with_eval = truth_expanded.merge(eval_df_expanded, how="outer", suffixes=("_truth", "_eval"), indicator=True, left_index=True, right_index=True)

truth_with_eval['overlap_start']=np.maximum(truth_with_eval.event_exon_start_truth, truth_with_eval.event_exon_start_eval)
truth_with_eval['overlap_end']=np.minimum(truth_with_eval.event_exon_end_truth, truth_with_eval.event_exon_end_eval)

truth_with_eval['overlap_length']=np.maximum(truth_with_eval.overlap_end - truth_with_eval.overlap_start,0).fillna(0)

total_truth_length = (truth_expanded.event_exon_end - truth_expanded.event_exon_start).sum()
total_eval_length = (eval_df_expanded.event_exon_end - eval_df_expanded.event_exon_start).sum()

# Sensitivity
sensitivity = truth_with_eval['overlap_length'].sum()/total_truth_length

# Precision
precision = truth_with_eval['overlap_length'].sum()/total_eval_length

# Expected events seen
def calculate_truth_event_overlap(df):
    overlap_length = df.overlap_length.sum()
    truth_dedup = df.drop_duplicates(subset=['exon_idx_truth'])
    truth_dedup_length = (truth_dedup.event_exon_end_truth - truth_dedup.event_exon_start_truth).sum()
    return overlap_length / truth_dedup_length

truth_events_fraction_covered = truth_with_eval.groupby(['event_index']).apply(calculate_truth_event_overlap)
expected_events_seen = truth_events_fraction_covered.sum()

# Expected events not seen
expected_events_not_seen = len(truth) - expected_events_seen

# Unexpected events seen
def calculate_eval_event_overlap(df):
    overlap_length = df.overlap_length.sum()
    eval_dedup = df.drop_duplicates(subset=['exon_idx_eval'])
    eval_dedup_length = (eval_dedup.event_exon_end_eval - eval_dedup.event_exon_start_eval).sum()
    return overlap_length / eval_dedup_length

eval_events_fraction_covered = truth_with_eval.groupby(['ID']).apply(calculate_eval_event_overlap)
unexpected_events_seen = len(eval_df) - eval_events_fraction_covered.sum()

print(f"Total Truth Length: {total_truth_length}")
#eval_only_events.to_csv('unexpected_events.tsv', sep='\t', index=False)