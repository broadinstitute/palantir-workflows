version 1.0

workflow cnv_control_events_qc {
    input {
        File eval_control_sample
        String eval_control_sample_name

        Float exon_sensitivity_threshold = 0.7

        File control_sample_common_events_na12878
        File control_sample_common_events_na24385
    }

    call cnv_control_events_qc_task {
        input: 
            eval_control_sample = eval_control_sample,
            eval_control_sample_name = eval_control_sample_name,
            exon_sensitivity_threshold = exon_sensitivity_threshold,
            control_sample_common_events_na12878 = control_sample_common_events_na12878,
            control_sample_common_events_na24385 = control_sample_common_events_na24385
    }

    output {
        Boolean qc_passed = cnv_control_events_qc_task.qc_passed
        Float sensitivity = cnv_control_events_qc_task.sensitivity
    }
}

task cnv_control_events_qc_task {
    input {
        File eval_control_sample
        String eval_control_sample_name

        Float exon_sensitivity_threshold

        File control_sample_common_events_na12878
        File control_sample_common_events_na24385
    }

    command <<<
        python3 << "EOF"
import pandas as pd

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
    return df_expanded


def read_vcf_to_df(vcf_path, sample_name, sample_alias):
    df = pd.read_csv(vcf_path, sep='\t',comment='#',compression='gzip',
                         names=['contig','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'])
    
    df = df.query("ALT != '.' and FILTER=='PASS'").copy()
    df['end'] = df['INFO'].str.split(";").apply(lambda x: x[0]).str.split('=').apply(lambda x: x[1]).astype(int)
    return df[['contig', 'start', 'end', 'ALT']]

if '~{eval_control_sample_name}' == 'NA12878':
    control_sample_common_events_path = '~{control_sample_common_events_na12878}'
elif '~{eval_control_sample_name}' == 'NA24385':
    control_sample_common_events_path = '~{control_sample_common_events_na24385}'
else:
    raise RuntimeError('Invalid control sample name. Must be either "NA12878" or "NA24385".')

truth = pd.read_table(control_sample_common_events_path)
eval_control_events = read_vcf_to_df(eval_control_sample_path, 'eval', '')

merged = truth.merge(eval_control_events[['contig', 'exon_idx', 'ALT']], how='left', on=['contig', 'exon_idx', 'ALT'], indicator=True)

sensitivity = len(merged[merged['_merge'] == 'both']) / len(truth)
precision = len(merged[merged['_merge'] == 'both']) / len(eval_control_events)

if sensitivity >= ~{exon_sensitivity_threshold}:
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
        Boolean qc_passed = read_string('qc_passed.txt')
        Float sensitivity = read_string('sensitivity.txt')
    }
}