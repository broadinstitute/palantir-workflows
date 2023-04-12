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
        Int cpu
        Int numThreads
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
    }

    String outputPrefix = if defined(referenceAnnotation) then "IsoQuant_out_~{datasetName}" else "IsoQuant_denovo_out_~{datasetName}"
    String completeGeneDBOption = if defined(referenceAnnotation) then "--complete_genedb" else ""
    String? referenceAnnotationBasename = basename(referenceAnnotation, ".reduced.gtf")

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
    >>>

    output {
        File? isoQuantGTF = "IsoQuant_out_~{datasetName}.gtf"
        File? isoQuantOut = "IsoQuant_out_~{datasetName}.tar.gz"
        File? isoQuantDB = "IsoQuant_out_~{datasetName}/~{referenceAnnotationBasename}.reduced.db"
        File? isoQuantDenovoGTF = "IsoQuant_denovo_out_~{datasetName}.gtf"
        File? isoQuantDenovoOut = "IsoQuant_denovo_out_~{datasetName}.tar.gz"
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
        Int cpu
        Int numThreads
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
    }

    String referenceAnnotationOption = if defined(referenceAnnotation) then "-G ~{referenceAnnotation}" else ""

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        stringtie \
        -o "StringTie_out_~{datasetName}.gtf" \
        ~{referenceAnnotationOption} \
        -p ~{numThreads} \
        -L ~{inputBAM}
    >>>

    output {
        File? stringTieGTF = "StringTie_out_~{datasetName}.gtf"
        File? stringTieDenovoGTF = "StringTie_denovo_out_~{datasetName}.gtf"
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
        Int cpu
        Int numThreads
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
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
        Int cpu
        Int numThreads
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
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
        Int cpu
        Int numThreads
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
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
        Int cpu
        Int numThreads
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
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
        Int cpu
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
    }

    String outputPrefix = "TAMA_out_~{datasetName}"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        python2 /usr/local/src/tama/tama_collapse.py -b BAM -s ~{inputBAM} -f ~{referenceGenome} -p ~{outputPrefix} -x capped
    >>>

    output {
        File tamaBED = "~{outputPrefix}.bed"
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
        Int cpu
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
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

task ReducedAnnotationGFFCompare {
    input {
        File reducedAnnotationDB
        File expressedGTF
        File expressedKeptGTF
        File excludedGTF
        File inputGTF
        File? counts
        String toolName
        String datasetName
        Int cpu
        Int memoryGB
        Int diskSizeGB
        String docker
    }

    String reducedAnnotationPrefix = basename(reducedAnnotationDB, ".reduced.db")

    command <<<
        mv ~{reducedAnnotationDB} .
        mv ~{excludedGTF} .
        mv ~{expressedGTF} .
        mv ~{expressedKeptGTF} .

        cp ~{counts} .

        /usr/local/src/IsoQuant-3.1.1/misc/reduced_db_gffcompare.py \
        --genedb ~{reducedAnnotationPrefix} \
        --gtf ~{inputGTF} \
        --tool ~{toolName} \
        --output "~{datasetName}_~{toolName}_reduced_db"

        tar -zcvf ~{datasetName}_~{toolName}_reduced_db.tar.gz ~{datasetName}_~{toolName}_reduced_db
    >>>

    output {
        File gffCompareOutput = "~{datasetName}_~{toolName}_reduced_db.tar.gz"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task DenovoAnnotationGFFCompare {
    input {
        File isoQuantGTF
        File stringTieGTF
        File bambuGTF
        File bambuGTFCounts
        File flairGTF
        File talonGTF
        String datasetName
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    Int cpu = 8
    Int memory = 64
    Int diskSizeGB = 300
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    String isoQuantBasename = basename(isoQuantGTF)
    String stringTieBasename = basename(stringTieGTF)
    String bambuBasename = basename(bambuGTF)
    String flairBasename = basename(flairGTF)
    String talonBasename = basename(talonGTF)

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mv ~{isoQuantGTF} .
        mv ~{stringTieGTF} .
        mv ~{bambuGTF} .
        mv ~{bambuGTFCounts} .
        mv ~{flairGTF} .
        mv ~{talonGTF} .

        echo "~{isoQuantBasename} isoquant" > gtfs.list
        echo "~{stringTieBasename} stringtie" >> gtfs.list
        echo "~{bambuBasename} bambu" >> gtfs.list
        echo "~{flairBasename} flair" >> gtfs.list
        echo "~{talonBasename} talon" >> gtfs.list

        /usr/local/src/IsoQuant-3.1.1/misc/denovo_model_stats.py \
        --gtf_list gtfs.list \
        --output "~{datasetName}_denovo_stats"

        tar -zcvf ~{datasetName}_denovo_stats.tar.gz ~{datasetName}_denovo_stats
    >>>

    output {
        File gffCompareOutput = "~{datasetName}_denovo_stats.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task ReferenceFreeGFFCompare {
    input {
        File isoQuantDenovoGTF
        File stringTieDenovoGTF
        File expressedGTF
        String datasetName
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    Int cpu = 8
    Int memory = 64
    Int diskSizeGB = 300
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir ~{datasetName}_reference_free_stats

        gffcompare \
        -r ~{expressedGTF} \
        -o "~{datasetName}_isoquant_gffcompare" ~{isoQuantDenovoGTF}

        gffcompare \
        -r ~{expressedGTF} \
        -o "~{datasetName}_stringtie_gffcompare" ~{stringTieDenovoGTF}

        mv ~{datasetName}* ~{datasetName}_reference_free_stats

        tar -zcvf ~{datasetName}_reference_free_stats.tar.gz ~{datasetName}_reference_free_stats
    >>>

    output {
        File gffCompareOutput = "~{datasetName}_reference_free_stats.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task ReducedAnalysisSummarize {
    input {
        File reducedGffCompareOutIsoQuant
        File reducedGffCompareOutStringTie
        File reducedGffCompareOutBambu
        File reducedGffCompareOutFlair
        File reducedGffCompareOutTalon
        String datasetName
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/kockan-reduced-analysis-summarize:latest"
    Int cpu = 1
    Int memory = 32
    Int diskSizeGB = 50

    command <<<
        cp ~{reducedGffCompareOutIsoQuant} .
        cp ~{reducedGffCompareOutStringTie} .
        cp ~{reducedGffCompareOutBambu} .
        cp ~{reducedGffCompareOutFlair} .
        cp ~{reducedGffCompareOutTalon} .

        tar -xzvf ~{datasetName}_isoquant_reduced_db.tar.gz
        tar -xzvf ~{datasetName}_stringtie_reduced_db.tar.gz
        tar -xzvf ~{datasetName}_bambu_reduced_db.tar.gz
        tar -xzvf ~{datasetName}_flair_reduced_db.tar.gz
        tar -xzvf ~{datasetName}_talon_reduced_db.tar.gz

        python3 /usr/local/src/plot_isoquant_results.py \
            ~{datasetName}_talon_reduced_db/talon.novel.stats,~{datasetName}_flair_reduced_db/flair.novel.stats,~{datasetName}_bambu_reduced_db/bambu.novel.stats,~{datasetName}_stringtie_reduced_db/stringtie.novel.stats,~{datasetName}_isoquant_reduced_db/isoquant.novel.stats \
            talon,flair,bambu,stringtie,isoquant \
            ~{datasetName}
    >>>

    output {
        File reducedAnalysisSummary = "~{datasetName}_reduced_analysis_summary.tsv"
        File reducedAnalysisAccuracyPlots = "~{datasetName}_reduced_analysis_summary.png"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}