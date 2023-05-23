version 1.0

task IsoQuant {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation
        String datasetName
        String dataType
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant@sha256:9bd8cd8c3a04e02599e10e0b484127fb763a39499302d4c859d230942f9a2d15"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String outputPrefix = if defined(referenceAnnotation) then "IsoQuant_out_~{datasetName}" else "IsoQuant_denovo_out_~{datasetName}"
    String completeGeneDBOption = if defined(referenceAnnotation) then "--complete_genedb" else ""

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/IsoQuant-3.2.0/isoquant.py \
        --reference ~{referenceGenome} \
        ~{completeGeneDBOption} \
        ~{"--genedb " + referenceAnnotation} \
        --bam ~{inputBAM} \
        --data_type ~{dataType} \
        --threads ~{numThreads} \
        --labels ~{datasetName} \
        --output ~{outputPrefix}
    >>>

    output {
        File isoQuantGTF = "~{outputPrefix}/~{datasetName}/~{datasetName}.transcript_models.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task StringTie {
    input {
        File inputBAM
        File? referenceAnnotation
        String datasetName
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie@sha256:ca2a163c7acdcacba741ea98d573080c15f153de18bd1566d24e8d2f1729ce89"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String outputPrefix = if defined(referenceAnnotation) then "StringTie_out_~{datasetName}" else "StringTie_denovo_out_~{datasetName}"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        stringtie \
        -o "~{outputPrefix}.gtf" \
        ~{"-G " + referenceAnnotation} \
        -p ~{numThreads} \
        -L ~{inputBAM}
    >>>

    output {
        File stringTieGTF = "~{outputPrefix}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task Bambu {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
        String dataType
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bambu@sha256:109fcdec65637eaca9f465808f3cc2aba3a9d2a0b1f967b4ed1c87989c3969de"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String bambuOutDir = "Bambu_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir ~{bambuOutDir}

        Rscript -<< "EOF"
        library(bambu)
        fa.file <- "~{referenceGenome}"
        gtf.file <- "~{referenceAnnotation}"
        bambuAnnotations <- prepareAnnotations(gtf.file)
        lr.bam <- "~{inputBAM}"
        lr.se <- bambu(reads = lr.bam, rcOutDir = "~{bambuOutDir}", annotations = bambuAnnotations, genome = fa.file, ncore = ~{numThreads})
        writeBambuOutput(lr.se, path = "~{bambuOutDir}")
        EOF

        awk ' $3 >= 1 ' ~{bambuOutDir}/counts_transcript.txt | sort -k3,3n > ~{bambuOutDir}/expressed_annotations.gtf.counts
        cut -f1 ~{bambuOutDir}/expressed_annotations.gtf.counts > ~{bambuOutDir}/expressed_transcripts.txt
        grep -Ff ~{bambuOutDir}/expressed_transcripts.txt ~{bambuOutDir}/extended_annotations.gtf > ~{bambuOutDir}/Bambu_out_~{datasetName}.gtf
    >>>

    output {
        File bambuGTF = "~{bambuOutDir}/Bambu_out_~{datasetName}.gtf"
        File bambuCounts = "~{bambuOutDir}/expressed_annotations.gtf.counts"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task Flair {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 500
        String docker = "brookslab/flair@sha256:994a5f6dd6bee041c8a2a82e84b77293d9bf5f3a2f172d440a72daee33474043"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String flairPrefix = "Flair_out_~{datasetName}"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        samtools fastq ~{inputBAM} > "~{flairPrefix}_temp.fastq"

        bam2Bed12 -i ~{inputBAM} > "~{flairPrefix}.bed"

        flair correct \
        -q "~{flairPrefix}.bed" \
        -g ~{referenceGenome} \
        -f ~{referenceAnnotation} \
        -o ~{flairPrefix} \
        -t ~{numThreads}

        flair collapse \
        -g ~{referenceGenome} \
        -f ~{referenceAnnotation} \
        -r "~{flairPrefix}_temp.fastq" \
        -q "~{flairPrefix}_all_corrected.bed" \
        -o ~{flairPrefix} \
        -t ~{numThreads}
    >>>

    output {
        File flairGTF = "Flair_out_~{datasetName}.isoforms.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task Talon {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
        String dataType
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/talon@sha256:7b9d92e0fa6e3f83c164a51c2eab9382e80ae963c61ec686ed462d8634bf8009"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String talonPrefix = "Talon_out_~{datasetName}"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        talon_label_reads --f ~{inputBAM} --t 1 --o ~{talonPrefix} --g ~{referenceGenome}
        samtools calmd -@ ~{numThreads} --reference ~{referenceGenome} "~{talonPrefix}_labeled.sam" > "~{talonPrefix}_labeled.md.sam"
        talon_initialize_database --f ~{referenceAnnotation} --g ~{datasetName} --a ~{datasetName} --o ~{datasetName}
        echo ~{datasetName},~{datasetName},~{dataType},"~{talonPrefix}_labeled.md.sam" > "~{talonPrefix}.csv"
        talon --build ~{datasetName} --db "~{datasetName}.db" --o "~{talonPrefix}_raw" --f "~{talonPrefix}.csv" --threads ~{numThreads}
        talon_filter_transcripts --db "~{datasetName}.db" -a ~{datasetName} --datasets ~{datasetName} --o "~{talonPrefix}_filter"
        talon_create_GTF --build ~{datasetName} --db "~{datasetName}.db" -a ~{datasetName} --o ~{talonPrefix} --whitelist "~{talonPrefix}_filter"
    >>>

    output {
        File talonGTF = "Talon_out_~{datasetName}_talon.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task IsoSeq {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        String datasetName
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoseq3@sha256:e715dda61f295d6825c0f4bea5133d583158db5d63c550ff186ee59f1ff10385"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        samtools fastq ~{inputBAM} > temp.fastq

        pbmm2 align --num-threads ~{numThreads} --preset ISOSEQ --sort ~{referenceGenome} temp.fastq  pbmm_realigned.bam

        isoseq3 collapse pbmm_realigned.bam "IsoSeq_out_~{datasetName}.gff"
    >>>

    output {
        File isoSeqGFF = "IsoSeq_out_~{datasetName}.gff"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task Flames {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
        Int cpu = 16
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/flames@sha256:e9b5d5152179e1a820afde3b147586a8ce7440738bf456af74b22ca4cfa0e8cb"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir fq

        samtools fastq ~{inputBAM} > ./fq/temp.fastq

        python3 /usr/local/src/FLAMES/python/bulk_long_pipeline.py \
        --gff3 ~{referenceAnnotation} \
        --genomefa ~{referenceGenome} \
        --fq_dir ./fq \
        --inbam ~{inputBAM} \
        --outdir .
    >>>

    output {
        File flamesGFF = "isoform_annotated.gff3"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task Cupcake {
    input {
        File inputBAM
        File inputBAMIndex
        String datasetName
        Int cpu = 16
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/cdna-cupcake@sha256:fca085dde170c995b5691d07aae2d56ab4426b7651a913957f029e628a0167c2"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String outputPrefix = "Cupcake_out_~{datasetName}"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        samtools fastq ~{inputBAM} > temp.fastq

        python3 /usr/local/src/remove_fastq_duplicates.py temp.fastq

        python3 /usr/local/src/cDNA_Cupcake/cupcake/tofu/collapse_isoforms_by_sam.py \
        --input out.fastq --fq \
        --bam ~{inputBAM} \
        --prefix ~{outputPrefix} \
        --cpus 1
    >>>

    output {
        File cupcakeGFF = "~{outputPrefix}.collapsed.gff"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task SplitGTF {
    input {
        File inputGTF
        File? inputCounts
        String toolName
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:f2c488bc2f75ab9b84934c5937e12dea015bea644fe1e362d0e5d5dec6f21200"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:f2c488bc2f75ab9b84934c5937e12dea015bea644fe1e362d0e5d5dec6f21200"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:f2c488bc2f75ab9b84934c5937e12dea015bea644fe1e362d0e5d5dec6f21200"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:f2c488bc2f75ab9b84934c5937e12dea015bea644fe1e362d0e5d5dec6f21200"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:f2c488bc2f75ab9b84934c5937e12dea015bea644fe1e362d0e5d5dec6f21200"
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

task SplitFreeStats {
    input {
        File trackingFile
        String toolName
        String datasetName
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:f2c488bc2f75ab9b84934c5937e12dea015bea644fe1e362d0e5d5dec6f21200"
    }

    command <<<
        python3 /usr/local/src/generate_split_free_benchmark_stats.py \
        --tracking ~{trackingFile} \
        --tool-name ~{toolName} \
        --dataset-name ~{datasetName}
    >>>

    output {
        File splitFreeStats = "~{toolName}_~{datasetName}_accuracy_stats.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}