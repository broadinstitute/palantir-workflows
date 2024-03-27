version 1.0

workflow Glimpse2CheckQC {
    input {
        File qc_metrics
        File qc_metrics_thresholds
        String output_basename
        
        Int preemptible = 3
        String docker = "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        Int cpu = 1
        Int mem_gb = 4
    }

    call Glimpse2CheckQCTask {
        input:
            qc_metrics = qc_metrics,
            qc_metrics_thresholds = qc_metrics_thresholds,
            output_basename = output_basename,
            preemptible = preemptible,
            docker = docker,
            cpu = cpu,
            mem_gb = mem_gb
    }

    output {
        Boolean qc_passed = Glimpse2CheckQCTask.qc_passed
        File qc_failures = Glimpse2CheckQCTask.qc_failures
    }
}

task Glimpse2CheckQCTask {
    input {
        File qc_metrics
        File qc_metrics_thresholds
        String output_basename
        
        Int preemptible
        String docker
        Int cpu
        Int mem_gb
    }

    Int disk_size_gb = 10
    
    command <<<
        set -euo pipefail

        cat <<'EOF' > script.py
import pandas as pd

data = pd.read_csv('~{qc_metrics}', sep='\t')
qc_metric_thresholds = pd.read_csv('~{qc_metrics_thresholds}', sep='\t')

data = data.melt(id_vars=['s'], var_name='metric', value_name='value', value_vars=qc_metric_thresholds.metric)
data = data.merge(qc_metric_thresholds, on=['metric'])

samples_out_of_spec = data.loc[(data['value'] < data['min']) | (data['value'] > data['max'])].sort_values(['s', 'metric'])
samples_out_of_spec.rename(columns={'s': 'sample_id'}).to_csv('~{output_basename}.qc_failures.tsv', sep='\t', index=False)

with open('~{output_basename}.qc_passed.txt', 'w') as qc_passed:
    qc_passed.write('true\n' if len(samples_out_of_spec) == 0 else 'false\n')
EOF
        python3 script.py
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        Boolean qc_passed = read_boolean("~{output_basename}.qc_passed.txt")
        File qc_failures = "~{output_basename}.qc_failures.tsv"
    }
}
