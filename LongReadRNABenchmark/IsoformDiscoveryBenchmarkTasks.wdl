version 1.0

task SplitGTF {
    input {
        File inputGTF
        File? inputCounts
        String toolName
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:3b3e48a2360ae12dc3f015cb94932388b62dae2e1098f5c86f29946d0f45b75b"
    }

    String base = if toolName == "flames" then basename(inputGTF, ".gff3") else basename(inputGTF, ".gtf")

    command <<<
        python3 /usr/local/src/split_gtf.py \
        --input-gtf ~{inputGTF} \
        --tool ~{toolName} \
        ~{"--input-bambu-counts " + inputCounts}
    >>>

    output {
        File full = "~{base}.full.gtf"
        File known = "~{base}.known.gtf"
        File novel = "~{base}.novel.gtf"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task ReducedAnnotationAnalysis {
    input {
        File inputFullGTF
        File inputKnownGTF
        File inputNovelGTF
        File expressedGTF
        File expressedKeptGTF
        File excludedGTF
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:b7208e67cb52ef41f0b9f9182414b8f12617a079546bbc2a4dbd826590ec63d2"
    }

    String baseFull = basename(inputFullGTF, ".gtf")
    String baseKnown = basename(inputKnownGTF, ".gtf")
    String baseNovel = basename(inputNovelGTF, ".gtf")

    command <<<
        gffcompare -r ~{expressedGTF} -o ~{baseFull} ~{inputFullGTF}
        gffcompare -r ~{expressedKeptGTF} -o ~{baseKnown} ~{inputKnownGTF}
        gffcompare -r ~{excludedGTF} -o ~{baseNovel} ~{inputNovelGTF}
    >>>

    output {
        File full = "~{baseFull}"
        File known = "~{baseKnown}"
        File novel = "~{baseNovel}"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task ReferenceFreeAnalysis {
    input {
        File inputGTF
        File expressedGTF
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:b7208e67cb52ef41f0b9f9182414b8f12617a079546bbc2a4dbd826590ec63d2"
    }

    String base = basename(inputGTF, ".gtf")

    command <<<
        gffcompare -r ~{expressedGTF} -o ~{base}.denovo ~{inputGTF}
    >>>

    output {
        File stats = "~{base}.denovo"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task DenovoAnalysis {
    input {
        String datasetName
        String splitType
        Array[File]+ gtfList
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:b7208e67cb52ef41f0b9f9182414b8f12617a079546bbc2a4dbd826590ec63d2"
    }

    command <<<
        gffcompare -o ~{datasetName}_~{splitType} ~{sep=" " gtfList}
    >>>

    output {
        File tracking = "~{datasetName}_~{splitType}.tracking"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task DenovoStats {
    input {
        String datasetName
        String splitType
        File trackingFile
        Int numTools = 6
        Array[String]+ toolNames
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:3b3e48a2360ae12dc3f015cb94932388b62dae2e1098f5c86f29946d0f45b75b"
    }

    command <<<
        python3 /usr/local/src/extract_denovo_model_stats.py \
        --dataset-name ~{datasetName} \
        --split-type ~{splitType} \
        --tracking ~{trackingFile} \
        --tool-names ~{sep=" " toolNames} \
        --num-tools ~{numTools}
    >>>

    output {
        File denovoStats = "~{datasetName}_~{splitType}_denovo_model_stats.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task SummarizeAnalysis {
    input {
        Array[File]+ inputList
        Array[String]+ toolNames
        String datasetName
        String analysisType
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:3b3e48a2360ae12dc3f015cb94932388b62dae2e1098f5c86f29946d0f45b75b"
    }

    command <<<
        python3 /usr/local/src/summarize_results.py \
        --input-list ~{sep=" " inputList} \
        --tool-names ~{sep=" " toolNames} \
        --dataset-name ~{datasetName} \
        --analysis-type ~{analysisType}
    >>>

    output {
        File summary = "~{datasetName}_~{analysisType}_analysis_summary.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task PlotAnalysisSummary {
    input {
        File summary
        String datasetName
        String analysisType
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:3b3e48a2360ae12dc3f015cb94932388b62dae2e1098f5c86f29946d0f45b75b"
    }

    command <<<
        python3 /usr/local/src/plot_summary_results.py \
        --input ~{summary} \
        --dataset-name ~{datasetName} \
        --analysis-type ~{analysisType} \
        --save
    >>>

    output {
        File analysisSummaryPlots = "~{datasetName}_~{analysisType}_analysis_summary.png"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task PlotDenovoStats {
    input {
        File stats
        String datasetName
        String splitType
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:3b3e48a2360ae12dc3f015cb94932388b62dae2e1098f5c86f29946d0f45b75b"
    }

    command <<<
        python3 /usr/local/src/plot_denovo_stats.py \
        --input ~{stats} \
        --dataset-name ~{datasetName} \
        --split-type ~{splitType} \
        --save
    >>>

    output {
        File denovoStatsPlot = "~{datasetName}_~{splitType}_denovo.png"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task GenerateSplitFreeTracking {
    input {
        String datasetName
        File toolGTF
        File expressedGTF
        File expressedKeptGTF
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:b7208e67cb52ef41f0b9f9182414b8f12617a079546bbc2a4dbd826590ec63d2"
    }

    command <<<
        gffcompare -o ~{datasetName} ~{expressedGTF} ~{expressedKeptGTF} ~{toolGTF}
    >>>

    output {
        File tracking = "~{datasetName}.tracking"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task GenerateSplitFreeTrackingDenovo {
    input {
        String datasetName
        Array[File]+ toolGTFs
        File expressedKeptGTF
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:b7208e67cb52ef41f0b9f9182414b8f12617a079546bbc2a4dbd826590ec63d2"
    }

    command <<<
        gffcompare -o ~{datasetName}.denovo ~{expressedKeptGTF} ~{sep=" " toolGTFs}
    >>>

    output {
        File tracking = "~{datasetName}.denovo.tracking"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task SplitFreeStats {
    input {
        Array[File]+ trackingFiles
        Array[String]+ toolNames
        String datasetName
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:3b3e48a2360ae12dc3f015cb94932388b62dae2e1098f5c86f29946d0f45b75b"
    }

    command <<<
        python3 /usr/local/src/generate_split_free_benchmark_stats.py \
        --tracking ~{sep=" " trackingFiles} \
        --tool-name ~{sep=" " toolNames} \
        --dataset-name ~{datasetName}
    >>>

    output {
        File splitFreeStats = "~{datasetName}_accuracy_stats.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}