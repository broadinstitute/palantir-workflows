version 1.0

workflow cnv_control_events_qc {
    input {
        File eval_control_sample
        String eval_control_sample_name

        Float exon_sensitivity_threshold = 0.7
        Float exon_precision_threshold = 0.7

        File control_sample_common_events_na12878
        File control_sample_common_events_na24385
        File exon_intervals
    }

    call cnv_control_events_qc_task {
        input: 
            eval_control_sample = eval_control_sample,
            eval_control_sample_name = eval_control_sample_name,
            exon_sensitivity_threshold = exon_sensitivity_threshold,
            exon_precision_threshold = exon_precision_threshold,
            control_sample_common_events_na12878 = control_sample_common_events_na12878,
            control_sample_common_events_na24385 = control_sample_common_events_na24385,
            exon_intervals = exon_intervals
    }

    output {
        Boolean qc_passed = cnv_control_events_qc_task.qc_passed
        Float sensitivity = cnv_control_events_qc_task.sensitivity
        Float precision = cnv_control_events_qc_task.precision
    }
}

task cnv_control_events_qc_task {
    input {
        File eval_control_sample
        String eval_control_sample_name

        Float exon_sensitivity_threshold
        Float exon_precision_threshold

        File control_sample_common_events_na12878
        File control_sample_common_events_na24385
        File exon_intervals
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
    return df.astype({'start_exon_idx': int, 'end_exon_idx': int})


def get_exon_expanded_events(df, exons):
    add_exon_idxs(df, exons)
    df = df.loc[df.start_exon_idx != df.end_exon_idx].reset_index()
    df_expanded = df.loc[df.index.repeat(df.end_exon_idx-df.start_exon_idx)]
    df_expanded['exon_idx'] = df_expanded.groupby(df_expanded.index).cumcount() + df_expanded.start_exon_idx

    # Add exon start and end positions
    df_expanded = df_expanded.merge(exons[['contig', 'contig_idx', 'start','end']], left_on=['contig', 'exon_idx'], right_on=['contig', 'contig_idx'], suffixes=('', '_exon'))

    # Calculate overlap start and end positions
    df_expanded['overlap_start']=np.maximum(df_expanded.start, df_expanded.start_exon)
    df_expanded['overlap_end']=np.minimum(df_expanded.end, df_expanded.end_exon)
    return df_expanded

def read_vcf_to_df(vcf_path):
    df = pd.read_csv(vcf_path, sep='\t',comment='#',compression='gzip',
                         names=['contig','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'])
    
    # Filter passing events
    df = df.query("ALT != '.' and FILTER=='PASS'").copy()

    # Extract end position
    df['end'] = df['INFO'].str.split(";").apply(lambda x: x[0]).str.split('=').apply(lambda x: x[1]).astype(int)

    return df[['contig', 'start', 'end', 'ALT']]

# Select the right control sample
if '~{eval_control_sample_name}' == 'NA12878':
    control_sample_common_events_path = '~{control_sample_common_events_na12878}'
elif '~{eval_control_sample_name}' == 'NA24385':
    control_sample_common_events_path = '~{control_sample_common_events_na24385}'
else:
    raise RuntimeError('Invalid control sample name. Must be either "NA12878" or "NA24385".')

# Read intervals
intervals = pd.read_csv('~{exon_intervals}', sep="\t", comment="@", names = ["contig","start","end","dummy1","dummy2"])
intervals['contig_idx'] = intervals.groupby('contig').cumcount()
intervals['exon_length'] = intervals['end'] - intervals['start'] + 1

# Read truth and eval data
truth = pd.read_table(control_sample_common_events_path)
eval_control_events = read_vcf_to_df('~{eval_control_sample}')
eval_control_events_expanded = get_exon_expanded_events(add_exon_idxs(eval_control_events, intervals), intervals)

# Full outer join, includes TP, FP, and FN
merged = truth.merge(eval_control_events_expanded, how='outer', on=['contig', 'exon_idx', 'ALT'], indicator=True)

# Add exon_length from intervals
merged = merged.merge(intervals[['contig', 'contig_idx', 'exon_length']], left_on=['contig', 'exon_idx'], right_on=['contig', 'contig_idx'], how='left')

# For exons that overlap, calculate the overlap length. For non-overlapping exons, overlap_length is 0.
merged['overlap_length'] = (merged['overlap_end'] - merged['overlap_start'] + 1).fillna(0)

# overlap_length is only > 0 for overlapping exons, therefore this is the number of TP
TP = merged.overlap_length.sum()

# FP are entire exons that are in eval_control_events_expanded but not in truth
# For these, the overlap_length is 0, so we can sum the exon_length for these rows
FP = merged.query('_merge == "right_only"').exon_length.sum()

# FN are all the bases that are not in eval_control_events_expanded but are in truth
# This is the difference between the exon_length of all exons, and the overlap_length,
# which is only greater than 0 for overlapping exons
FN = (merged.exon_length - merged.overlap_length).sum()

sensitivity = TP / (TP + FN)
precision = TP / (TP + FP)

if sensitivity >= ~{exon_sensitivity_threshold}:
    qc_passed = True
else:
    qc_passed = False

if qc_passed and precision >= ~{exon_precision_threshold}:
    qc_passed = True
else:
    qc_passed = False

with open('qc_passed.txt', 'w') as f:
    f.write(str(qc_passed))

with open('sensitivity.txt', 'w') as f:
    f.write(f'{sensitivity:.4f}')

with open('precision.txt', 'w') as f:
    f.write(f'{precision:.4f}')
EOF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "1 GiB"
        preemptible: 1
    }

    output {
        Boolean qc_passed = read_boolean('qc_passed.txt')
        Float sensitivity = read_float('sensitivity.txt')
        Float precision = read_float('precision.txt')
    }
}