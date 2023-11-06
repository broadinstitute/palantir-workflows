version 1.0

task DetectPCANoveltiesTask{
    input {
        File test
        File training
        Float? alpha
        Float? distanceThreshold
        Int cpu = 1
        Int memoryGB = 8
        Int diskSizeGB = 64
        String docker = "us.gcr.io/broad-dsde-methods/kockan/alphashape@sha256:ea0679e7b22250ddf641772ac12b45c71d950bed690b29edd6fc1d8862e0643b"
    }

    command <<<
        python3 /usr/local/src/pca_novelty_detection.py \
        --input ~{test} \
        --training-data ~{training}
    >>>

    output {
        File testSetPredictions = "out.tsv"
        File runVisualization = "out_visualization.png"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow DetectPCANovelties {
    input {
        File test
        File training
        Float? alpha
        Float? distanceThreshold
    }

    call DetectPCANoveltiesTask {
        input:
            test = test,
            training = training,
            alpha = alpha,
            distanceThreshold = distanceThreshold
    }

    output {
        File testSetPredictions = DetectPCANoveltiesTask.testSetPredictions
        File runVisualization = DetectPCANoveltiesTask.runVisualization
    }
}