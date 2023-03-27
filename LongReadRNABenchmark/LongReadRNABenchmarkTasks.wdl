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
    >>>

    output {
        File isoQuantGTF = "IsoQuant_out_~{datasetName}/~{datasetName}/~{datasetName}.transcript_models.gtf"
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
    >>>

    output {
        File isoQuantDenovoGTF = "IsoQuant_denovo_out_~{datasetName}/~{datasetName}/~{datasetName}.transcript_models.gtf"
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

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        BAMBU_OUT_DIR="Bambu_out"
        echo "library(bambu)" > bambu.R
        echo "fa.file<-~{referenceGenome}" >> bambu.R
        echo "gtf.file<-~{referenceAnnotation}" >> bambu.R
        echo "bambuAnnotations <- prepareAnnotations(gtf.file)"
        echo "ont.bam<-~{inputBAM}" >> bambu.R
        echo "ont.se <- bambu(reads = ont.bam, rcOutDir = \"$BAMBU_OUT_DIR\", annotations = bambuAnnotations, genome = fa.file, ncore = ~{numThreads}" >> bambu.R
        echo "writeBambuOutput(ont.se, path = \"$BAMBU_OUT_DIR\")" >> bambu.R
        R < bambu.R --no-save

        BAMBU_GTF = $BAMBU_OUT_DIR/"Bambu_out_"$DATASET_NAME".gtf"
        awk ' $3 >= 1 ' $BAMBU_OUT_DIR/counts_transcript.txt | sort -k3,3n > $BAMBU_OUT_DIR/expressed_annotations.gtf.counts
        cut -f1 $BAMBU_OUT_DIR/expressed_annotations.gtf.counts > $BAMBU_OUT_DIR/expressed_transcripts.txt
        grep -Ff $BAMBU_OUT_DIR/expressed_transcripts.txt $BAMBU_OUT_DIR/extended_annotations.gtf > $BAMBU_GTF
        cp $BAMBU_OUT_DIR/expressed_annotations.gtf.counts $BAMBU_GTF".counts"
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
        File inputReads
        File referenceGenome
        File referenceAnnotation
        String datasetName
        Int numThreads
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/flair:latest"
    Int cpu = 16
    Int memory = 256
    Int diskSizeGB = 500
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        FLAIR_PREFIX="Flair_out_~{datasetName}
        FLAIR_GTF=$FLAIR_PREFIX".gtf"
        $FLAIR_PATH"bin/bam2Bed12.py -i ~{inputBAM} > $FLAIR_PREFIX".bed"

        $FLAIR_PATH"flair.py" correct \
        -q $FLAIR_PREFIX".bed" \
        -g ~{referenceGenome} \
        -f ~{referenceAnnotation} \
        -o $FLAIR_PREFIX \
        -t ~{numThreads}

        $FLAIR_PATH"flair.py" collapse \
        -g ~{referenceGenome} \
        -f ~{referenceAnnotation} \
        -r ~{inputReads} \
        -q $FLAIR_PREFIX"_all_corrected.bed" \
        -o $FLAIR_PREFIX \
        -t ~{numThreads}
    >>>

    output {
        File flairGTF = "Flair_out_~{datasetName}.gtf"
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

    # TALON requires config file: dataset name, sample description, platform, sam file (comma-delimited)

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        TALON_PREFIX="Talon_out_~{datasetName}
        TALON_GTF=$TALON_PREFIX".gtf"
        TALON_PATH"talon_label_reads" --f ~{inputBAM} --t=~{numThreads} --o=$TALON_PREFIX --g ~{referenceGenome}
        samtools calmd -@ ~{numThreads} --reference ~{referenceGenome} $TALON_PREFIX"_labeled.sam" > $TALON_PREFIX"_labeled.md.sam"
        $TALON_PATH"talon_initialize_database" --f ~{referenceAnnotation} --g ~{datasetName} --a ~{datasetName} --o ~{datasetName}
        echo ~{datasetName},~{datasetName},ONT,"$TALON_PREFIX"_labeled.md.sam" > $TALON_PREFIX".csv"
        $TALON_PATH"talon" --build ~{datasetName} --db "~{datasetName}.db" --o $TALON_PREFIX"_raw" --f $TALON_PREFIX".csv"
        $TALON_PATH"talon_filter_transcripts" --db "~{datasetName}.db" -a ~{datasetName} --datasets ~{datasetName} --o $TALON_PREFIX"_filter" --f $TALON_CSV
        $TALON_PATH"talon_create_GTF" -- build ~{datasetName} --db "~{datasetName}.db" -a ~{datasetName} --o $TALON_PREFIX --whitelist=$TALON_PREFIX"_filter"
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

        ls -lha "~{datasetName}_isoquant_reduced_db"

        /usr/local/src/IsoQuant-3.1.1/misc/reduced_db_gffcompare.py \
        --genedb ~{reducedAnnotationPrefix} \
        --gtf ~{stringTieGTF} \
        --tool "stringtie" \
        --output "~{datasetName}_stringtie_reduced_db"

        ls -lha "~{datasetName}_stringtie_reduced_db"
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
        #File bambuGTF
        #File flairGTF
        #File talonGTF
        String datasetName
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    Int cpu = 8
    Int memory = 64
    Int diskSizeGB = 300
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        echo ~{isoQuantGTF} > gtfs.list
        echo ~{stringTieGTF} >> gtfs.list

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