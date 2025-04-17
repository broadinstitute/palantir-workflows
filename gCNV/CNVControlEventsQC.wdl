version 1.0

# Use CNVControlEventsQCGenerateTruth.ipynb to generate the control sample common events file

workflow CNVControlEventsQC {
    input {
        File eval_control_sample
        File control_sample_common_events
        File exon_intervals

        Float threshold_min_sensitivity
        Float threshold_min_precision
        Float threshold_min_expected_events_seen
        Float threshold_max_expected_events_not_seen
        Int threshold_max_unexpected_events_seen
    }

    call CNVControlEventsQCTask {
        input: 
            eval_control_sample = eval_control_sample,
            control_sample_common_events = control_sample_common_events,
            exon_intervals = exon_intervals,
            threshold_max_expected_events_not_seen = threshold_max_expected_events_not_seen,
            threshold_max_unexpected_events_seen = threshold_max_unexpected_events_seen,
            threshold_min_expected_events_seen = threshold_min_expected_events_seen,
            threshold_min_precision = threshold_min_precision,
            threshold_min_sensitivity = threshold_min_sensitivity
    }

    output {
        Boolean qc_passed = CNVControlEventsQCTask.qc_passed
        Float sensitivity = CNVControlEventsQCTask.sensitivity
        Float precision = CNVControlEventsQCTask.precision
        Float expected_events_seen = CNVControlEventsQCTask.expected_events_seen
        Float expected_events_not_seen = CNVControlEventsQCTask.expected_events_not_seen
        Int unexpected_events_seen = CNVControlEventsQCTask.unexpected_events_seen
        File unexpected_events = CNVControlEventsQCTask.unexpected_events
    }
}

task CNVControlEventsQCTask {
    input {
        File eval_control_sample
        File control_sample_common_events
        File exon_intervals

        Float threshold_min_sensitivity
        Float threshold_min_precision
        Float threshold_min_expected_events_seen
        Float threshold_max_expected_events_not_seen
        Int threshold_max_unexpected_events_seen
    }

    command <<<
        python3 << "EOF"
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

intervals = read_intervals_to_df('~{exon_intervals}')
truth = read_truth('~{control_sample_common_events}')
eval_df = read_vcf_to_df('~{eval_control_sample}')
eval_df_expanded = get_exon_expanded_events(eval_df, intervals)

truth_expanded = get_exon_expanded_events(truth, intervals)

truth_with_eval = truth_expanded.merge(eval_df_expanded, how="outer", suffixes=("_truth", "_eval"), indicator=True, left_index=True, right_index=True)

truth_with_eval['overlap_start']=np.maximum(truth_with_eval.event_exon_start_truth, truth_with_eval.event_exon_start_eval)
truth_with_eval['overlap_end']=np.minimum(truth_with_eval.event_exon_end_truth, truth_with_eval.event_exon_end_eval)

truth_with_eval['overlap_length']=np.maximum(truth_with_eval.overlap_end - truth_with_eval.overlap_start,0).fillna(0)

truth_with_eval['truth_length']=truth_with_eval.event_exon_end_truth - truth_with_eval.event_exon_start_truth
truth_with_eval['eval_length']=np.maximum(truth_with_eval.event_exon_end_eval - truth_with_eval.event_exon_start_eval,0).fillna(0)

truth_length_dedup = truth_with_eval.drop_duplicates(subset=['event_index', 'exon_idx_truth']).truth_length.sum()
eval_length_dedup = truth_with_eval.drop_duplicates(subset=['ID', 'exon_idx_eval']).eval_length.sum()

# Sensitivity
sensitivity = truth_with_eval['overlap_length'].sum()/truth_length_dedup

# Precision
precision = truth_with_eval['overlap_length'].sum()/eval_length_dedup

# Expected events seen
def group_truth_events(df):
    overlap_length = df.overlap_length.sum()
    truth_length_dedup = df.drop_duplicates(subset=['exon_idx_truth']).truth_length.sum()
    return overlap_length / truth_length_dedup

truth_events_fraction_covered = truth_with_eval.groupby(['event_index']).apply(group_truth_events)
expected_events_seen = truth_events_fraction_covered.sum()

# Expected events not seen
expected_events_not_seen = len(truth) - expected_events_seen

# Unexpected events seen
event_ids_matching_truth = truth_with_eval.loc[truth_with_eval['_merge'] == 'both', 'ID'].unique()
eval_only_events = eval_df[~eval_df['ID'].isin(event_ids_matching_truth)]
unexpected_events_seen = len(eval_only_events)

# QC pass
qc_passed = ((sensitivity >= ~{threshold_min_sensitivity}) and
    (precision >= ~{threshold_min_precision}) and
    (expected_events_seen >= ~{threshold_min_expected_events_seen}) and 
    (expected_events_not_seen <= ~{threshold_max_expected_events_not_seen}) and 
    (unexpected_events_seen <= ~{threshold_max_unexpected_events_seen}))

# Write results
with open('precision.txt', 'w') as f:
    f.write(f'{precision:.4f}\n')
with open('sensitivity.txt', 'w') as f:
    f.write(f'{sensitivity:.4f}\n')
with open('expected_events_seen.txt', 'w') as f:
    f.write(f'{expected_events_seen:.4f}\n')
with open('expected_events_not_seen.txt', 'w') as f:
    f.write(f'{expected_events_not_seen:.4f}\n')
with open('unexpected_events_seen.txt', 'w') as f:
    f.write(f'{unexpected_events_seen}\n')
eval_only_events.to_csv('unexpected_events.tsv', sep='\t', index=False)
with open('qc_passed.txt', 'w') as f:
    f.write(f'{qc_passed}\n')
EOF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.1"
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "1 GiB"
        preemptible: 1
    }

    output {
        Boolean qc_passed = read_boolean('qc_passed.txt')
        Float sensitivity = read_float('sensitivity.txt')
        Float precision = read_float('precision.txt')
        Float expected_events_seen = read_float('expected_events_seen.txt')
        Float expected_events_not_seen = read_float('expected_events_not_seen.txt')
        Int unexpected_events_seen = read_int('unexpected_events_seen.txt')
        File unexpected_events = 'unexpected_events.tsv'
    }
}
