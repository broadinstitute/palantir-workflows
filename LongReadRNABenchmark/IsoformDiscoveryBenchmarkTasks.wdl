version 1.0

task GffCompareTrack {
    input {
        String datasetName
        String toolName
        File toolGTF
        File expressedGTF
        File expressedKeptGTF
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:b7208e67cb52ef41f0b9f9182414b8f12617a079546bbc2a4dbd826590ec63d2"
    }

    command <<<
        gffcompare -o ~{datasetName}.~{toolName} ~{expressedGTF} ~{expressedKeptGTF} ~{toolGTF}
    >>>

    output {
        File tracking = "~{datasetName}.~{toolName}.tracking"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task GffCompareTrackDenovo {
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

task ReferenceFreeAnalysis {
    input {
        File inputGTF
        File expressedGTF
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:b7208e67cb52ef41f0b9f9182414b8f12617a079546bbc2a4dbd826590ec63d2"
    }

    String base = select_first([basename(inputGTF, ".gtf"), basename(inputGTF, ".gff")])

    command <<<
        gffcompare -r ~{expressedGTF} -o ~{base}.reffree ~{inputGTF}
    >>>

    output {
        File stats = "~{base}.reffree"
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
        String datasetName
        Array[File]+ trackingFiles
        Array[String]+ toolNames
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:21c4ba6d275c9ddf18f1250b3aa35066255e5c5574a994f4dfa1b69d20a679d7"
    }

    command <<<
        python3 /usr/local/src/summarize_analysis.py \
        --tracking ~{sep=" " trackingFiles} \
        --tool-names ~{sep=" " toolNames} \
        --dataset-name ~{datasetName}
    >>>

    output {
        File summary = "~{datasetName}_analysis_summary.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task SummarizeReferenceFreeAnalysis {
    input {
        String datasetName
        Array[File]+ inputList
        Array[String]+ toolNames
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:21c4ba6d275c9ddf18f1250b3aa35066255e5c5574a994f4dfa1b69d20a679d7"
    }

    command <<<
        python3 /usr/local/src/summarize_reffree_analysis.py \
        --input-list ~{sep=" " inputList} \
        --tool-names ~{sep=" " toolNames} \
        --dataset-name ~{datasetName}
    >>>

    output {
        File summary = "~{datasetName}_analysis_summary_reffree.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task SummarizeDenovoAnalysis {
    input {
        String datasetName
        File trackingFile
        Array[String]+ toolNames
        Int numTools = 6
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:21c4ba6d275c9ddf18f1250b3aa35066255e5c5574a994f4dfa1b69d20a679d7"
    }

    command <<<
        python3 /usr/local/src/summarize_denovo_analysis.py \
        --dataset-name ~{datasetName} \
        --tracking ~{trackingFile} \
        --tool-names ~{sep=" " toolNames} \
        --num-tools ~{numTools}
    >>>

    output {
        File denovoSummaryNovel = "~{datasetName}_denovo_analysis_summary_novel.tsv"
        File denovoSummaryKnown = "~{datasetName}_denovo_analysis_summary_known.tsv"
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
        String type
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:21c4ba6d275c9ddf18f1250b3aa35066255e5c5574a994f4dfa1b69d20a679d7"
    }

    command <<<
        python3 /usr/local/src/plot_analysis_summary.py \
        --input ~{summary} \
        --dataset-name ~{datasetName} \
        --type ~{type} \
        --save
    >>>

    output {
        File analysisSummaryPlot = "~{datasetName}_analysis_summary_~{type}.png"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task PlotDenovoAnalysisSummary {
    input {
        File denovoSummary
        String datasetName
        String type
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:21c4ba6d275c9ddf18f1250b3aa35066255e5c5574a994f4dfa1b69d20a679d7"
    }

    command <<<
        python3 /usr/local/src/plot_denovo_analysis_summary.py \
        --input ~{denovoSummary} \
        --dataset-name ~{datasetName} \
        --type ~{type} \
        --save
    >>>

    output {
        File denovoAnalysisSummaryPlot = "~{datasetName}_analysis_summary_denovo_~{type}.png"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}