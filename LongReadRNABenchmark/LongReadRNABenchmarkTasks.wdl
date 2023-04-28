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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant:latest"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String outputPrefix = if defined(referenceAnnotation) then "IsoQuant_out_~{datasetName}" else "IsoQuant_denovo_out_~{datasetName}"
    String completeGeneDBOption = if defined(referenceAnnotation) then "--complete_genedb" else ""
    String referenceAnnotationBasename = if defined(referenceAnnotation) then basename(select_first([referenceAnnotation]), ".reduced.gtf") else ""

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/IsoQuant-3.1.1/isoquant.py \
        --reference ~{referenceGenome} \
        ~{completeGeneDBOption} \
        ~{"--genedb " + referenceAnnotation} \
        --bam ~{inputBAM} \
        --data_type ~{dataType} \
        --threads ~{numThreads} \
        --labels ~{datasetName} \
        --output ~{outputPrefix}

        mv ~{outputPrefix}/~{datasetName}/~{datasetName}.transcript_models.gtf ~{outputPrefix}.gtf

        tar -zcvf ~{outputPrefix}.tar.gz ~{outputPrefix}/~{datasetName}

        touch ~{outputPrefix}/~{referenceAnnotationBasename}.reduced.db
    >>>

    output {
        File isoQuantGTF = "~{outputPrefix}.gtf"
        File isoQuantOut = "~{outputPrefix}.tar.gz"
        File isoQuantDB = "~{outputPrefix}/~{referenceAnnotationBasename}.reduced.db"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie:latest"
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
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/bambu:latest"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String referenceAnnotationBasename = basename(referenceAnnotation, ".reduced.gtf")
    String bambuOutDir = "Bambu_out"
    String bambuGTFPath = "Bambu_out/Bambu_out_~{datasetName}.gtf"

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
        grep -Ff ~{bambuOutDir}/expressed_transcripts.txt ~{bambuOutDir}/extended_annotations.gtf > ~{bambuGTFPath}
        cp ~{bambuOutDir}/expressed_annotations.gtf.counts "~{bambuGTFPath}.counts"

        tar -zcvf ~{bambuOutDir}_~{datasetName}.tar.gz ~{bambuOutDir}
    >>>

    output {
        File bambuGTF = "Bambu_out/Bambu_out_~{datasetName}.gtf"
        File bambuGTFCounts = "Bambu_out/Bambu_out_~{datasetName}.gtf.counts"
        File bambuOut = "~{bambuOutDir}_~{datasetName}.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "brookslab/flair:latest"
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

        mv Flair_out_~{datasetName}.isoforms.gtf Flair_out_~{datasetName}.gtf
    >>>

    output {
        File flairGTF = "Flair_out_~{datasetName}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/talon:latest"
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

        mv Talon_out_~{datasetName}_talon.gtf Talon_out_~{datasetName}.gtf
    >>>

    output {
        File talonGTF = "Talon_out_~{datasetName}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoseq3:latest"
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
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task Tama {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        String datasetName
        Int cpu = 16
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us.gcr.io/broad-dsde-methods/kockan/tama:latest"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
    }

    String outputPrefix = "TAMA_out_~{datasetName}"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        python2 /usr/local/src/tama/tama_collapse.py -b BAM -s ~{inputBAM} -f ~{referenceGenome} -p ~{outputPrefix} -x capped

        ls -lha

        cat ~{outputPrefix}.bed | perl /usr/local/src/bed12ToGTF.1.pl > ~{outputPrefix}.gtf

        ls -lha
    >>>

    output {
        File tamaGTF = "~{outputPrefix}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/flames:latest"
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

        mv isoform_annotated.gff3 FLAMES_out_~{datasetName}.gff
    >>>

    output {
        File flamesGFF = "FLAMES_out_~{datasetName}.gff"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/cdna-cupcake:latest"
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

        mv ~{outputPrefix}.collapsed.gff ~{outputPrefix}.gff
    >>>

    output {
        File cupcakeGFF = "~{outputPrefix}.gff"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom:latest"
    }

    String base = if toolName == "flames" then basename(inputGTF, ".gff") else basename(inputGTF, ".gtf")

    command <<<
        ls -lha

        python3 /usr/local/src/split_gtf.py \
        --input-gtf ~{inputGTF} \
        --tool ~{toolName} \
        ~{"--input-bambu-counts " + inputCounts}

        ls -lha
    >>>

    output {
        File full = "~{base}.full.gtf"
        File known = "~{base}.known.gtf"
        File novel = "~{base}.novel.gtf"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare:latest"
    }

    String baseFull = basename(inputFullGTF, ".gtf")
    String baseKnown = basename(inputKnownGTF, ".gtf")
    String baseNovel = basename(inputNovelGTF, ".gtf")

    command <<<
        ls -lha

        gffcompare -r ~{expressedGTF} -o ~{baseFull} ~{inputFullGTF}
        gffcompare -r ~{expressedKeptGTF} -o ~{baseKnown} ~{inputKnownGTF}
        gffcompare -r ~{excludedGTF} -o ~{baseNovel} ~{inputNovelGTF}

        ls -lha
    >>>

    output {
        File full = "~{baseFull}"
        File known = "~{baseKnown}"
        File novel = "~{baseNovel}"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task ReferenceFreeAnalysis {
    input {
        File inputGTF
        File expressedGTF
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare:latest"
    }

    String base = basename(inputGTF, ".gtf")

    command <<<
        ls -lha

        gffcompare -r ~{expressedGTF} -o ~{base}.denovo ~{inputGTF}

        ls -lha
    >>>

    output {
        File stats = "~{base}.denovo"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task DenovoAnalysis {
    input {
        String toolName
        Array[File]+ gtfList
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/gffcompare:latest"
    }

    command <<<
        ls -lha

        gffcompare -o ~{toolName} ~{sep=" " gtfList}

        ls -lha
    >>>

    output {
        File tracking = "~{toolName}.tracking"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task DenovoStats {
    input {
        File trackingFile
        String toolName
        Int numTools = 6
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom:latest"
    }

    command <<<
        ls -lha

        python3 /usr/local/src/extract_denovo_model_stats.py --tool ~{toolName} --tracking ~{trackingFile} --num-tools ~{numTools}

        ls -lha
    >>>

    output {
        File gffCompareOutput = "~{toolName}_denovo_model_stats.tsv"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom:latest"
    }

    command <<<
        python3 /usr/local/src/summarize_results.py \
        --input-list ~{sep=" " inputList} \
        --tool-names ~~{sep=" " toolNames} \
        --dataset-name ~{datasetName} \
        --analysis-type ~{analysisType}
    >>>

    output {
        File summary = "~{datasetName}_~{analysisType}_analysis_summary.tsv"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
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
        String docker = "us.gcr.io/broad-dsde-methods/kockan/lr-isoform-reconstruction-benchmarking-custom:latest"
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
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}