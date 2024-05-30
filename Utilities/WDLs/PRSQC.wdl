version 1.0

# Simple PRS QC wdl for the PROGRESS VA project
workflow PRSQC {
    input {
        String sample_name
        File prs_control
        File prs_sample
        File control_thresholds
        File sample_thresholds
        File alphashape
        #Boolean qc_control_only = false
    }

    Int cpu = 1
    Int mem_gb = 4
    Int disk_size_gb = 32
    String docker = "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"

    call CheckThresholds as CheckControl {
        input:
            scores = prs_control,
            thresholds = control_thresholds,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb,
            docker = docker
    }

    call DetectPCANovelties as DetectPCANoveltiesControl {
        input:
            test_sample = prs_control,
            alphashape = alphashape,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb
    }

    call CheckThresholds as CheckSample {
        input:
            scores = prs_sample,
            thresholds = sample_thresholds,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb,
            docker = docker
    }

    call DetectPCANovelties as DetectPCANoveltiesSample {
        input:
            test_sample = prs_sample,
            alphashape = alphashape,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb
    }

    call FinalizeQCOutputs {
        input:
            sample_name = sample_name,
            qc_passed_control = CheckControl.qc_passed,
            qc_passed_sample = CheckSample.qc_passed,
            pca_qc_passed_control = DetectPCANoveltiesControl.pca_qc_passed,
            pca_qc_passed_sample = DetectPCANoveltiesSample.pca_qc_passed,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb,
            docker = docker
    }

    output {
        Boolean qc_passed = FinalizeQCOutputs.qc_passed
        File qc_failures_control = CheckControl.qc_failures
        File qc_failures_sample = CheckSample.qc_failures
    }
}

task CheckThresholds {
    input {
        File scores
        File thresholds

        Int cpu
        Int mem_gb
        Int disk_size_gb
        String docker
    }

    String output_basename = basename(scores, ".tsv")

    command <<<
        set -euo pipefail

        cat <<'EOF' > script.py
        import pandas as pd

        scores = pd.read_csv('~{scores}', sep = '\t', header = None, names = ["prs_score", "pc1", "pc2", "full_model_score"])
        thresholds = pd.read_csv('~{thresholds}', sep = '\t', header = None, names = ["min", "max"])

        # Confirm scores are within the boundaries
        prs_score_passed = True
        if scores.iloc[0,0] < thresholds.iloc[0,0] or scores.iloc[0,0] > thresholds.iloc[0,1]:
            prs_score_passed = False

        full_model_score_passed = True
        if scores.iloc[0,3] < thresholds.iloc[1,0] or scores.iloc[0,3] > thresholds.iloc[1,1]:
            full_model_score_passed = False

        with open('~{output_basename}.qc_passed.txt', 'w') as qc_passed:
            if prs_score_passed and full_model_score_passed:
                qc_passed.write("true\n")
            else:
                qc_passed.write("false\n")

        with open('~{output_basename}.qc_failures.txt', 'w') as qc_failures:
            if not prs_score_passed:
                prs_score_failure = "PRS_SCORE " + str(scores.iloc[0,0]) + " outside given boundaries with min: " + str(thresholds.iloc[0,0]) + " , max: " + str(thresholds.iloc[0,1])
                qc_failures.write(prs_score_failure + "\n")
            if not full_model_score_passed:
                full_model_score_failure = "FULL_MODEL_SCORE " + str(scores.iloc[0,3]) + " outside given boundaries with min: " + str(thresholds.iloc[1,0]) + " , max: " + str(thresholds.iloc[1,1])
                qc_failures.write(full_model_score_failure + "\n")

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

task DetectPCANovelties{
    input {
        File test_sample
        File alphashape
        Float distance_threshold = 0.01
        Int cpu
        Int mem_gb
        Int disk_size_gb
        String docker = "us.gcr.io/broad-dsde-methods/kockan/alphashape@sha256:96d5a34da2ff6da6e2cd0f85cca2fa2400c2d73e02e1957def111fbf04c2fdda"
    }

    String output_basename = basename(test_sample, ".tsv")

    command <<<
        set -euo pipefail

        cat <<'EOF' > detect_pca_novelties.py
        # Automatically flag novelties in 2D PCA plots using Concave Hulls generated via alphashapes
        import sys
        import pickle
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from descartes import PolygonPatch
        from shapely.geometry import Point
        import alphashape

        # Load alphashape
        with open("~{alphashape}", "rb") as infile:
            alpha_shape = pickle.load(infile)

        # Set the distance threshold
        dist_thresh = ~{distance_threshold}

        # Read test data
        # Expected input format: tsv with columns PRS_SCORE, PC1, PC2, and FULL_MODEL_SCORE
        # There is no header at the moment but it can be added to the format
        df_test = pd.read_csv("~{test_sample}", sep = '\t', header = None, names = ["prs_score", "pc1", "pc2", "full_model_score"])
        test_point = Point(df_test.iloc[0,1], df_test.iloc[0,2])

        # Prepare output plot for nice visualization
        fig, ax = plt.subplots()
        ax.add_patch(PolygonPatch(alpha_shape, alpha = 0.2))

        with open("~{output_basename}.pca_qc.txt", 'w') as outfile:
            # Test and label the sample as a novelty or a regular observation
            dist = test_point.distance(alpha_shape)
            if alpha_shape.contains(test_point) or dist < dist_thresh:
                outfile.write("true" + "\n")
                plt.scatter(test_point.x, test_point.y, c = 'green', alpha = 1.0, s = 10)
            else:
                outfile.write("false" + "\n")
                plt.scatter(test_point.x, test_point.y, c = 'red', alpha = 1.0, s = 10)

        # Plotting for nice visualization
        plt.title("Alphashape Automated Novelty Flagging: ~{output_basename}")
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        colors = ['green', 'red']
        labels = ['Pass', 'Fail']
        patches = [(plt.Line2D([], [], color = colors[i], label = "{:s}".format(labels[i]), marker = "o", linewidth = 0)) for i in range(len(labels))]
        plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), handles = patches)
        plt.savefig("~{output_basename}.pca_qc.png", dpi = 300, bbox_inches = 'tight')

        EOF
        python3 detect_pca_novelties.py
    >>>

    output {
        Boolean pca_qc_passed = read_boolean("~{output_basename}.pca_qc.txt")
        File pca_qc_plot = "~{output_basename}.pca_qc.png"
    }

    runtime {
        cpu: cpu
        memory: "~{mem_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: docker
    }
}

task FinalizeQCOutputs {
    input {
        String sample_name
        Boolean qc_passed_control
        Boolean qc_passed_sample
        Boolean pca_qc_passed_control
        Boolean pca_qc_passed_sample

        Int cpu
        Int mem_gb
        Int disk_size_gb
        String docker
    }

    String qpc = if (qc_passed_control) then "True" else "False"
    String qps = if (qc_passed_sample) then "True" else "False"
    String pqpc = if (pca_qc_passed_control) then "True" else "False"
    String pqps = if (pca_qc_passed_sample) then "True" else "False"

    command <<<
        set -euo pipefail

        cat <<'EOF' > script.py

        with open('~{sample_name}.qc_passed.txt', 'w') as qc_passed:
            if ~{qpc} and ~{qps} and ~{pqpc} and ~{pqps}:
                qc_passed.write("true\n")
            else:
                qc_passed.write("false\n")

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
        Boolean qc_passed = read_boolean("~{sample_name}.qc_passed.txt")
    }
}
