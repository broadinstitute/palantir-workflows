version 1.0

task IsoQuant {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
        String dataType
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    String referenceAnnotationBasename = basename(referenceAnnotation, ".reduced.gtf")

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/IsoQuant-3.1.1/isoquant.py \
        --reference ~{referenceGenome} \
        --complete_genedb \
        --genedb ~{referenceAnnotation} \
        --bam ~{inputBAM} \
        --data_type ~{dataType} \
        --threads ~{numThreads} \
        --labels ~{datasetName} \
        --output "IsoQuant_out_~{datasetName}"

        mv \
        IsoQuant_out_~{datasetName}/~{datasetName}/~{datasetName}.transcript_models.gtf \
        IsoQuant_out_~{datasetName}.gtf
    >>>

    output {
        File isoQuantGTF = "IsoQuant_out_~{datasetName}.gtf"
        File isoQuantDB = "IsoQuant_out_~{datasetName}/~{referenceAnnotationBasename}.reduced.db"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task IsoQuantReferenceFree {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        String datasetName
        String dataType
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/IsoQuant-3.1.1/isoquant.py \
        --reference ~{referenceGenome} \
        --bam ~{inputBAM} \
        --data_type ~{dataType} \
        --threads ~{numThreads} \
        --labels ~{datasetName} \
        --output "IsoQuant_denovo_out_~{datasetName}"

        mv \
        IsoQuant_denovo_out_~{datasetName}/~{datasetName}/~{datasetName}.transcript_models.gtf \
        IsoQuant_denovo_out_~{datasetName}.gtf
    >>>

    output {
        File isoQuantDenovoGTF = "IsoQuant_denovo_out_~{datasetName}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task StringTie {
    input {
        File inputBAM
        File referenceAnnotation
        String datasetName
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        stringtie \
        -o "StringTie_out_~{datasetName}.gtf" \
        -G ~{referenceAnnotation} \
        -p ~{numThreads} \
        -L \
        ~{inputBAM}
    >>>

    output {
        File stringTieGTF = "StringTie_out_~{datasetName}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task StringTieReferenceFree {
    input {
        File inputBAM
        String datasetName
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        stringtie \
        -o "StringTie_denovo_out_~{datasetName}.gtf" \
        -p ~{numThreads} \
        -L \
        ~{inputBAM}
    >>>

    output {
        File stringTieDenovoGTF = "StringTie_denovo_out_~{datasetName}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
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
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/bambu:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    String referenceAnnotationBasename = basename(referenceAnnotation, ".reduced.gtf")
    String bambuOutDir = "Bambu_out"
    String bambuGTFPath = "Bambu_out/Bambu_out_~{datasetName}.gtf"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir ~{bambuOutDir}

        echo "library(bambu)" > bambu.R
        echo "fa.file <- \"~{referenceGenome}\"" >> bambu.R
        echo "gtf.file <- \"~{referenceAnnotation}\"" >> bambu.R
        echo "bambuAnnotations <- prepareAnnotations(gtf.file)" >> bambu.R
        echo "lr.bam <- \"~{inputBAM}\"" >> bambu.R
        echo "lr.se <- bambu(reads = lr.bam, rcOutDir =\"~{bambuOutDir}\", annotations = bambuAnnotations, genome = fa.file, ncore = ~{numThreads})" >> bambu.R
        echo "writeBambuOutput(lr.se, path = \"~{bambuOutDir}\")" >> bambu.R

        R < bambu.R --no-save

        awk ' $3 >= 1 ' ~{bambuOutDir}/counts_transcript.txt | sort -k3,3n > ~{bambuOutDir}/expressed_annotations.gtf.counts
        cut -f1 ~{bambuOutDir}/expressed_annotations.gtf.counts > ~{bambuOutDir}/expressed_transcripts.txt
        grep -Ff ~{bambuOutDir}/expressed_transcripts.txt ~{bambuOutDir}/extended_annotations.gtf > ~{bambuGTFPath}
        cp ~{bambuOutDir}/expressed_annotations.gtf.counts "~{bambuGTFPath}.counts"

        ls -lha ~{bambuOutDir}
    >>>

    output {
        File bambuGTF = "Bambu_out/Bambu_out_~{datasetName}.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
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
        Int numThreads
    }

    String docker = "brookslab/flair:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

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
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
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
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/talon:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    String talonPrefix = "Talon_out_~{datasetName}"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        talon_label_reads --f ~{inputBAM} --t ~{numThreads} --o ~{talonPrefix} --g ~{referenceGenome}
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
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task ReducedAnnotationGFFCompare {
    input {
        File reducedAnnotationDB
        File expressedGTF
        File expressedKeptGTF
        File excludedGTF
        File isoQuantGTF
        File stringTieGTF
        String datasetName
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    Int cpu = 8
    Int memory = 64
    Int diskSizeGB = 300
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    String reducedAnnotationPrefix = basename(reducedAnnotationDB, ".reduced.db")

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mv ~{reducedAnnotationDB} .
        mv ~{excludedGTF} .
        mv ~{expressedGTF} .
        mv ~{expressedKeptGTF} .

        /usr/local/src/IsoQuant-3.1.1/misc/reduced_db_gffcompare.py \
        --genedb ~{reducedAnnotationPrefix} \
        --gtf ~{isoQuantGTF} \
        --tool "isoquant" \
        --output "~{datasetName}_isoquant_reduced_db"

        /usr/local/src/IsoQuant-3.1.1/misc/reduced_db_gffcompare.py \
        --genedb ~{reducedAnnotationPrefix} \
        --gtf ~{stringTieGTF} \
        --tool "stringtie" \
        --output "~{datasetName}_stringtie_reduced_db"
    >>>

    output {
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task DenovoAnnotationGFFCompare {
    input {
        File isoQuantGTF
        File stringTieGTF
        File bambuGTF
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
        mv ~{flairGTF} .
        mv ~{talonGTF} .

        echo "~{isoQuantBasename} isoquant" > gtfs.list
        echo "~{stringTieBasename} stringtie" >> gtfs.list
        echo "~{bambuBasename} bambu" >> gtfs.list
        echo "~{flairBasename} flair" >> gtfs.list
        echo "~{talonBasename} talon" >> gtfs.list

        ls -lha

        cat gtfs.list

        /usr/local/src/IsoQuant-3.1.1/misc/denovo_model_stats.py \
        --gtf_list gtfs.list \
        --output "~{datasetName}_denovo_stats"

        ls -lha
    >>>

    output {
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

        gffcompare \
        -r ~{expressedGTF} \
        -o "~{datasetName}_isoquant_gffcompare" ~{isoQuantDenovoGTF}

        gffcompare \
        -r ~{expressedGTF} \
        -o "~{datasetName}_stringtie_gffcompare" ~{stringTieDenovoGTF}
    >>>

    output {
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}