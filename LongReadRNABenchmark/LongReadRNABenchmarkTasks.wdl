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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant@sha256:8af2461f7bb9de2137172c8637b35bea1242b28ff367d42b08d62f2a98c3fc8d"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie@sha256:77371c6e81abbc3a0dd169ad21f861903f695ea73ffadca2c84e31651fffb548"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bambu@sha256:330ba8d5e9a70da486dfba3d2271739dc05cf5d157db85d6a1d006de8f1d8953"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/talon@sha256:07f5dda62d29976eded3b8c4afaaff176b370a1237c95253e4a550a4e3a6e629"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/flames@sha256:74bdb17cf2bf092f3358e0ae0edfea33783180d0f3a688e2362b82544049ad63"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/cdna-cupcake@sha256:1d9e8c05fad09223d9f9a0adc2d56b637463365726c3461691ec4aa1f84dcacc"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:775cd637fa3f4befc6b0b7dc1a674a0e64461920df8b10ab1b6e330e439049a3"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:e63da5cf04c83e766ec546bb5f6c0f7b9fef4ec5c5d68c6b4a527ba72012ea3b"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:e63da5cf04c83e766ec546bb5f6c0f7b9fef4ec5c5d68c6b4a527ba72012ea3b"
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
        String splitType
        Array[File]+ gtfList
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare@sha256:e63da5cf04c83e766ec546bb5f6c0f7b9fef4ec5c5d68c6b4a527ba72012ea3b"
    }

    command <<<
        gffcompare -o ~{splitType} ~{sep=" " gtfList}
    >>>

    output {
        File tracking = "~{splitType}.tracking"
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
        String splitType
        File trackingFile
        Int numTools = 6
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:775cd637fa3f4befc6b0b7dc1a674a0e64461920df8b10ab1b6e330e439049a3"
    }

    command <<<
        python3 /usr/local/src/extract_denovo_model_stats.py --split-type ~{splitType} --tracking ~{trackingFile} --num-tools ~{numTools}
    >>>

    output {
        File gffCompareOutput = "~{splitType}.denovo_model_stats.tsv"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:775cd637fa3f4befc6b0b7dc1a674a0e64461920df8b10ab1b6e330e439049a3"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom@sha256:775cd637fa3f4befc6b0b7dc1a674a0e64461920df8b10ab1b6e330e439049a3"
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