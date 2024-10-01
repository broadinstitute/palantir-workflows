version 1.0

# Simple PRS QC wdl for the PROGRESS VA project
workflow PRSQC {
    input {
        File prs_full_risk
        File acceptable_range
        String output_basename

        # This alphashape is essentially a bounding polygon containing the set of points in the form (PC1, PC2)
        # from a training set, such as 1kG. It was generated using the library: https://github.com/bellockk/alphashape
        # We have a script that allows generating alphashapes with an arbitrary training set (of valid format),
        # which can be found here: https://github.com/broadinstitute/palantir-workflows/blob/main/Utilities/Dockers/Alphashape/generate_alphashape.py
        # The main (and pretty much only one besides training data) parameter of interest is "alpha", which determines how tight
        # the bounding polygon is. Our default is set at 8.0, which was determined to be a good value experimentally.
        # Users who intend to use a training set different than the default (1kG) should be aware of this.
        File alphashape
    }

    Int cpu = 1
    Int mem_gb = 4
    Int disk_size_gb = 32
    String docker = "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"

    call CheckScores {
        input:
            prs_full_risk = prs_full_risk,
            acceptable_range = acceptable_range,
            output_basename = output_basename,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb,
            docker = docker
    }

    call DetectPCANovelties {
        input:
            prs_full_risk = prs_full_risk,
            output_basename = output_basename,
            alphashape = alphashape,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb
    }

    output {
        Boolean qc_passed = CheckScores.qc_passed && DetectPCANovelties.pcs_within_shape
        Boolean pcs_within_shape = DetectPCANovelties.pcs_within_shape
        File pca_qc_plot = DetectPCANovelties.pca_qc_plot
    }
}

task CheckScores {
    input {
        File prs_full_risk
        File acceptable_range
        String output_basename

        Int cpu
        Int mem_gb
        Int disk_size_gb
        String docker
    }

    command <<<
        set -euo pipefail

        cat <<'EOF' > check_scores.py
        import pandas as pd

        prs_full_risk = pd.read_csv('~{prs_full_risk}', sep = '\t', header = 0, index_col = 0)
        acceptable_range = pd.read_csv('~{acceptable_range}', sep = '\t', header = 0, index_col = 0)

        # Get sample ids as a list; works with both single-sample or multi-sample files
        sample_ids = prs_full_risk["sample_id"].tolist()

        all_metrics_within_range = True
        for sample_id in sample_ids:
            # Check whether prs_score and combined_risk_score are within acceptable range
            prs_score_within_range = prs_full_risk.loc[sample_id, "prs_score"] >= acceptable_range.loc["prs_score", "min"] and prs_full_risk.loc[sample_id, "prs_score"] <= acceptable_range.loc["prs_score", "max"]
            combined_risk_score_within_range = prs_full_risk.loc[sample_id, "combined_risk_score"] >= acceptable_range.loc["combined_risk_score", "min"] and prs_full_risk.loc[sample_id, "combined_risk_score"] <= acceptable_range.loc["combined_risk_score", "max"]

            # Check whether pc1 and pc2 are within acceptable range
            pc1_within_range = prs_full_risk.loc[sample_id, "pc1"] >= acceptable_range.loc["pc1", "min"] and prs_full_risk.loc[sample_id, "pc1"] <= acceptable_range.loc["pc1", "max"]
            pc2_within_range = prs_full_risk.loc[sample_id, "pc2"] >= acceptable_range.loc["pc2", "min"] and prs_full_risk.loc[sample_id, "pc2"] <= acceptable_range.loc["pc2", "max"]

            if not (prs_score_within_range and combined_risk_score_within_range and pc1_within_range and pc2_within_range):
                all_metrics_within_range = False
                break

        with open('~{output_basename}.qc_passed.txt', 'w') as qc_passed:
            if all_metrics_within_range:
                qc_passed.write("true\n")
            else:
                qc_passed.write("false\n")
        EOF
        python3 check_scores.py
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
    }

    output {
        Boolean qc_passed = read_boolean("~{output_basename}.qc_passed.txt")
    }
}

task DetectPCANovelties {
    input {
        File prs_full_risk
        String output_basename
        File alphashape
        Float distance_threshold = 0.01
        Int cpu
        Int mem_gb
        Int disk_size_gb
        String docker = "us.gcr.io/broad-dsde-methods/kockan/alphashape@sha256:96d5a34da2ff6da6e2cd0f85cca2fa2400c2d73e02e1957def111fbf04c2fdda"
    }

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
        prs_full_risk = pd.read_csv("~{prs_full_risk}", sep = '\t', header = 0, index_col = 0)

        # Get sample ids as a list; works with both single-sample or multi-sample files
        sample_ids = prs_full_risk["sample_id"].tolist()

        test_points = []
        for sample_id in sample_ids:
            test_points.append(Point(prs_full_risk.loc[sample_id, "pc1"], prs_full_risk.loc[sample_id, "pc2"]))

        # Prepare output plot for nice visualization
        fig, ax = plt.subplots()
        ax.add_patch(PolygonPatch(alpha_shape, alpha = 0.2))

        # Test and label each sample as a novelty or a regular observation
        all_pcs_within_shape = True
        for test_point in test_points:
            dist = test_point.distance(alpha_shape)
            if alpha_shape.contains(test_point) or dist < dist_thresh:
                plt.scatter(test_point.x, test_point.y, c = 'green', alpha = 1.0, s = 10)
            else:
                plt.scatter(test_point.x, test_point.y, c = 'red', alpha = 1.0, s = 10)
                all_pcs_within_shape = False

        with open("~{output_basename}.pca_qc.txt", 'w') as outfile:
            if all_pcs_within_shape:
                outfile.write("true" + "\n")
            else:
                outfile.write("false" + "\n")

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
        Boolean pcs_within_shape = read_boolean("~{output_basename}.pca_qc.txt")
        File pca_qc_plot = "~{output_basename}.pca_qc.png"
    }

    runtime {
        cpu: cpu
        memory: "~{mem_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: docker
    }
}
