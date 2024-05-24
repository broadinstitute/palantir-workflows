version 1.0

# Simple PRS QC wdl for the PROGRESS VA project
workflow PRSQC {
    input {
        File prs_sample
        File prs_control
        File sample_thresholds
        File control_thresholds
        File output_basename
    }

    Int cpu = 1
    Int mem_gb = 4
    Int disk_size_gb = 32
    String docker = "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"

    call PRSQCTask {
        input:
            prs_sample = prs_sample,
            prs_control = prs_control,
            sample_thresholds = sample_thresholds,
            control_thresholds = control_thresholds,
            output_basename = output_basename,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb,
            docker = docker
    }

    output {
        Boolean qc_passed = PRSQCTask.qc_passed
        File qc_failures = PRSQCTask.qc_failures
    }
}

task PRSQCTask {
    input {
        File prs_sample
        File prs_control
        File sample_thresholds
        File control_thresholds
        File output_basename

        Int cpu
        Int mem_gb
        Int disk_size_gb
        String docker
    }

    command <<<
        set -euo pipefail

        cat <<'EOF' > script.py
        import pandas as pd

        prs_sample = pd.read_csv('~{prs_sample}', sep = '\t', header = None, names = ["prs_score", "pc1", "pc2", "full_model_score"])
        prs_control = pd.read_csv('~{prs_control}', sep = '\t', header = None, names = ["prs_score", "pc1", "pc2", "full_model_score"])

        sample_thresholds = pd.read_csv('~{sample_thresholds}', sep = '\t', header = None, names = ["min", "max"])
        control_thresholds = pd.read_csv('~{control_thresholds}', sep = '\t', header = None, names = ["min", "max"])

        # Check control sample to confirm score and pcs are as expected
        control_pass = True
        if prs_control.iloc[0,0] < control_thresholds.iloc[0,0] or prs_control.iloc[0,0] > control_thresholds.iloc[0,1]:
            control_pass = False

        # Check sample of interest to confirm scores are within some boundaries
        sample_pass = True
        if prs_sample.iloc[0,0] < sample_thresholds.iloc[0,0] or prs_sample.iloc[0,0] > sample_thresholds.iloc[0,1]:
            sample_pass = False

        with open('~{output_basename}.qc_passed.txt', 'w') as qc_passed:
            if control_pass and sample_pass:
                qc_passed.write("true\n")
            else:
                qc_passed.write("false\n")
                with open('~{output_basename}.qc_failures.txt', 'w') as qc_failures:
                    if not control_pass:
                        qc_failures.write("ERROR: PRS score for control outside given thresholds.")
                    if not sample_pass:
                        qc_failures.write("ERROR: PRS score for sample outside given thresholds.")
        EOF
        python3 script.py
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
    }

    output {
        Boolean qc_passed = read_boolean("~{output_basename}.qc_passed.txt")
        File qc_failures = "~{output_basename}.qc_failures.txt"
    }
}