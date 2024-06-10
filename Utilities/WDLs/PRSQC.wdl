version 1.0

# Simple PRS QC wdl for the PROGRESS VA project
workflow PRSQC {
    input {
        String sample_name
        File prs_control
        File prs_sample
        File control_acceptable_range
        File sample_acceptable_range
        File control_expected_pcs
        Float margin_control_pcs
        File alphashape
        Boolean manually_passed_control = false
    }

    Int cpu = 1
    Int mem_gb = 4
    Int disk_size_gb = 32
    String docker = "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"

    if(!manually_passed_control) {
        call CheckScoresAgainstExpectedValues as CheckScoresAgainstExpectedControlValues {
            input:
                scores = prs_control,
                acceptable_range = control_acceptable_range,
                cpu = cpu,
                mem_gb = mem_gb,
                disk_size_gb = disk_size_gb,
                docker = docker
        }

        call CheckControlPCsAgainstExpectedValues {
            input:
                scores = prs_control,
                expected_values = control_expected_pcs,
                margin = margin_control_pcs,
                cpu = cpu,
                mem_gb = mem_gb,
                disk_size_gb = disk_size_gb,
                docker = docker
        }
    }

    call CheckScoresAgainstExpectedValues as CheckScoresAgainstExpectedSampleValues {
        input:
            scores = prs_sample,
            acceptable_range = sample_acceptable_range,
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

    output {
        Boolean qc_passed = CheckScoresAgainstExpectedSampleValues.qc_passed && DetectPCANoveltiesSample.pca_qc_passed && select_first([CheckScoresAgainstExpectedControlValues.qc_passed, true]) && select_first([CheckControlPCsAgainstExpectedValues.pca_qc_passed, true])
        File? qc_failures_control = CheckScoresAgainstExpectedControlValues.qc_failures
        File qc_failures_sample = CheckScoresAgainstExpectedSampleValues.qc_failures
    }
}

task CheckScoresAgainstExpectedValues {
    input {
        File scores
        File acceptable_range

        Int cpu
        Int mem_gb
        Int disk_size_gb
        String docker
    }

    String output_basename = basename(scores, ".tsv")

    command <<<
        set -euo pipefail

        cat <<'EOF' > check_scores_against_expected_values.py
        import pandas as pd

        scores = pd.read_csv('~{scores}', sep = '\t', header = 0)
        acceptable_range = pd.read_csv('~{acceptable_range}', sep = '\t', header = 0)

        # Check whether PRS_SCORE and FULL_MODEL_SCORE are within the acceptable range
        prs_score_passed = scores.loc[0, "PRS_SCORE"] >= acceptable_range.loc[0, "MIN"] and scores.loc[0, "PRS_SCORE"] <= acceptable_range.loc[0, "MAX"]
        full_model_score_passed = scores.loc[0, "FULL_MODEL_SCORE"] >= acceptable_range.loc[1, "MIN"] and scores.loc[0, "FULL_MODEL_SCORE"] <= acceptable_range.loc[1, "MAX"]

        with open('~{output_basename}.qc_passed.txt', 'w') as qc_passed:
            if prs_score_passed and full_model_score_passed:
                qc_passed.write("true\n")
            else:
                qc_passed.write("false\n")

        with open('~{output_basename}.qc_failures.txt', 'w') as qc_failures:
            if not prs_score_passed:
                qc_failures.write("PRS_SCORE " + str(scores.loc[0, "PRS_SCORE"]) + " outside acceptable range with MIN: " + str(acceptable_range.loc[0, "MIN"]) + " , MAX: " + str(acceptable_range.loc[0, "MAX"]) + "\n")
            if not full_model_score_passed:
                qc_failures.write("FULL_MODEL_SCORE " + str(scores.loc[0, "FULL_MODEL_SCORE"]) + " outside acceptable range with MIN: " + str(acceptable_range.loc[1, "MIN"]) + " , MAX: " + str(acceptable_range.loc[1, "MAX"]) + "\n")

        EOF
        python3 check_scores_against_expected_values.py
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

task CheckControlPCsAgainstExpectedValues {
    input {
        File scores
        File expected_values
        Float margin

        Int cpu
        Int mem_gb
        Int disk_size_gb
        String docker
    }

    String output_basename = basename(scores, ".tsv")

    command <<<
        set -euo pipefail

        cat <<'EOF' > check_control_pcs_against_expected_values.py
        import pandas as pd

        scores = pd.read_csv('~{scores}', sep = '\t', header = 0)
        expected_values = pd.read_csv('~{expected_values}', sep = '\t', header = 0)

        # Check whether PC1 and PC2 are within the acceptable range
        pc1_within_range = abs(scores.loc[0, "PC1"] - expected_values.loc[0, "PC1"]) <= ~{margin}
        pc2_within_range = abs(scores.loc[0, "PC2"] - expected_values.loc[0, "PC2"]) <= ~{margin}

        with open('~{output_basename}.pca_qc_passed.txt', 'w') as pca_qc_passed:
            if pc1_within_range and pc2_within_range:
                pca_qc_passed.write("true\n")
            else:
                pca_qc_passed.write("false\n")

        EOF
        python3 check_control_pcs_against_expected_values.py
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
    }

    output {
        Boolean pca_qc_passed = read_boolean("~{output_basename}.pca_qc_passed.txt")
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

        # Read input data
        df_test = pd.read_csv("~{test_sample}", sep = '\t', header = 0)
        test_point = Point(df_test.loc[0, "PC1"], df_test.loc[0, "PC2"])

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
