version 1.0

workflow Glimpse2ImputationQC {
    input {
        File imputed_vcf
        
        Int preemptible = 1
        String docker = "hailgenetics/hail:0.2.126-py3.11"
        Int cpu = 4
        Int mem_gb = 16
        File? monitoring_script
    }

    call Glimpse2ImputationQCTask {
        input:
            imputed_vcf = imputed_vcf,
            preemptible = preemptible,
            docker = docker,
            cpu = cpu,
            mem_gb = mem_gb,
            monitoring_script = monitoring_script
    }

    output {
        File qc_metrics = Glimpse2ImputationQCTask.qc_metrics
        Boolean qc_passed = Glimpse2ImputationQCTask.qc_passed
        File qc_failures = Glimpse2ImputationQCTask.qc_failures
    }
}

task Glimpse2ImputationQCTask {
    input {
        File imputed_vcf
        
        Int preemptible = 1
        String docker = "hailgenetics/hail:0.2.126-py3.11"
        Int cpu = 4
        Int mem_gb = 16
        File? monitoring_script
    }

    parameter_meta {
        imputed_vcf: {
                        localization_optional: true
                    }
    }

    Int disk_size_gb = 100
    
    command <<<
        set -euo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        cat <<'EOF' > script.py
import hail as hl
import pandas as pd

# Calculate metrics
hl.init(default_reference='GRCh38', idempotent=True)
vcf = hl.import_vcf('~{imputed_vcf}', force_bgz=True)
qc = hl.sample_qc(vcf)
qc.cols().flatten().export('qc.tsv')

# QC check
data = pd.read_csv('qc.tsv', sep='\t')

metrics = ['n_hom_ref', 'n_het', 'n_hom_var', 'n_non_ref', 'n_snp', 'n_insertion', 'n_deletion', 'n_transition',
           'n_transversion', 'n_star', 'r_ti_tv', 'r_het_hom_var', 'r_insertion_deletion']

qc_thresholds = {
    'n_hom_ref':        (80_000_000, 100_000_000),
    'n_het':            ( 1_000_000,   6_000_000),
    'n_hom_var':        (   800_000,   3_000_000),
    'n_non_ref':        ( 2_500_000,   7_500_000),
    'n_snp':            ( 3_000_000,   8_000_000),
    'n_insertion':      (   200_000,   1_000_000),
    'n_deletion':       (   200_000,   1_000_000),
    'n_transition':     ( 2_000_000,   6_000_000),
    'n_transversion':   ( 1_000_000,   3_000_000),
    'n_star':           (0, 0),
    'r_ti_tv':          (2, 2.1),
    'r_het_hom_var':    (0.5, 5),
    'r_insertion_deletion': (0.98, 1.1)
}

with open('qc_failures.txt', 'w') as qc_failures:
    all_ok = True
    for metric in metrics:
        out_of_spec_samples = data.loc[(data['sample_qc.' + metric] < qc_thresholds[metric][0]) | (data['sample_qc.' + metric] > qc_thresholds[metric][1])]
        for index, row in out_of_spec_samples.iterrows():
            all_ok = False
            qc_failures.write(f'Sample {row.s} out of spec for {metric}. Value: {row["sample_qc." + metric]}. Acceptable range: [{qc_thresholds[metric][0]}, {qc_thresholds[metric][1]}].\n')

    with open('qc_passed.txt', 'w') as qc_passed:
        qc_passed.write('true\n' if all_ok else 'false\n')
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
        File qc_metrics = "qc.tsv"
        Boolean qc_passed = read_boolean("qc_passed.txt")
        File qc_failures = "qc_failures.txt"
        File? monitoring = "monitoring.log"
    }
}
