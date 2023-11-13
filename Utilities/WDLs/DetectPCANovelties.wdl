version 1.0

task GenerateAlphashape{
    input {
        File training
        Float alpha = 8.0
        Int cpu = 1
        Int memoryGB = 8
        Int diskSizeGB = 64
        String docker = "us.gcr.io/broad-dsde-methods/kockan/alphashape@sha256:96d5a34da2ff6da6e2cd0f85cca2fa2400c2d73e02e1957def111fbf04c2fdda"
    }

    command <<<
        python3 /usr/local/src/generate_alphashape.py \
        --training-data ~{training} \
        --alpha ~{alpha}
    >>>

    output {
        File alphashapePickle = "alphashape.pickle"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task DetectPCANoveltiesTask{
    input {
        File test
        File alphashape
        File training
        Float distanceThreshold = 0.01
        Int cpu = 1
        Int memoryGB = 8
        Int diskSizeGB = 64
        String docker = "us.gcr.io/broad-dsde-methods/kockan/alphashape@sha256:96d5a34da2ff6da6e2cd0f85cca2fa2400c2d73e02e1957def111fbf04c2fdda"
    }

    command <<<
        python3 /usr/local/src/pca_novelty_detection.py \
        --input ~{test} \
        --alphashape ~{alphashape} \
        --training-data ~{training} \
        --distance-threshold ~{distanceThreshold}
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
    }

    call GenerateAlphashape {
        input:
            training = training
    }

    call DetectPCANoveltiesTask {
        input:
            test = test,
            alphashape = GenerateAlphashape.alphashapePickle,
            training = training
    }

    output {
        File testSetPredictions = DetectPCANoveltiesTask.testSetPredictions
        File runVisualization = DetectPCANoveltiesTask.runVisualization
    }
}