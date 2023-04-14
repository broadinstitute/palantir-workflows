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

        cat ~{outputPrefix}.bed | \
        perl -<< "EOF" > ~{outputPrefix}.gtf
        #!/usr/bin/perl
        use strict;
        #use warnings;

        # A script to convert BED12 into a GTF.
        # Input: piped BED12 file.
        # Example: cat clark.gtfs.lift_to_hg38/Clark_DataS3.NoncodingTranscriptsAll.lift_to_hg38.bed | perl /users/rg/rjohnson/science/scripts/bed12ToGTF/bed12ToGTF.1.pl

        my @line=();

        my $trxchr;
        my $trxstart;
        my $trxend;
        my $trxstrand;
        my $trxid;
        my $trxblank;
        my @blocksizes;
        my @blockstarts;

        my $thisstart;
        my $thisend;

        my $exontotal;

        while (<STDIN>) {
            chomp $_;
            @line=split("\t",$_);

            ($trxchr, $trxstart, $trxend, $trxid, $trxblank,$trxstrand)=($line[0],$line[1],$line[2],$line[3],$line[4], $line[5]);

            $trxstart+=1;   # Conversion from BED to GTF!!

            my @trxtoks=split(";",$trxid);

            # Print the Transcript Line
            print "$trxchr\tblank\ttranscript\t$trxstart\t$trxend\t.\t$trxstrand\t.\ttranscript_id \"$trxtoks[0].$trxtoks[1]\"\;\n";

            @blocksizes=split(",", $line[10]);
            @blockstarts=split(",", $line[11]);
            $exontotal=scalar(@blockstarts);

            #print "\n @blocksizes hi @blockstarts";

            my $exon_count=0;
            my $rev_exon_count=$exontotal+1;

            for (my $i=0; $i < $exontotal; $i++) {
                $exon_count++;
                $rev_exon_count--;

                if ($trxstrand eq "+"){
                    $thisstart=$trxstart+$blockstarts[$i];
                    $thisend=$trxstart+$blockstarts[$i]+$blocksizes[$i]-1;    # The -1 added empirically after browser inspection
                    #Print the Exon lines.
                    print "$trxchr\tblank\texon\t$thisstart\t$thisend\t.\t$trxstrand\t.\ttranscript_id \"$trxtoks[0].$trxtoks[1]\"; exon_number $exon_count\;\n";
                }

                elsif ($trxstrand eq "-"){
                    #$thisend=$trxend-$blockstarts[$i];
                    #$thisstart=$thisend-$blocksizes[$i]+1;   # The +1 added empirically after browser inspection
                    $thisstart=$trxstart+$blockstarts[$i];
                    $thisend=$trxstart+$blockstarts[$i]+$blocksizes[$i]-1;    # The -1 added empirically after browser inspection
                    #Print the Exon lines.
                    print "$trxchr\tblank\texon\t$thisstart\t$thisend\t.\t$trxstrand\t.\ttranscript_id \"$trxtoks[0].$trxtoks[1]\"; exon_number $rev_exon_count\;\n";
                }
            }
        }
        EOF
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
        Int cpu = 6
        Int memoryGB = 64
        Int diskSizeGB = 300
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    }

    String reducedAnnotationPrefix = basename(reducedAnnotationDB, ".reduced.db")

    command <<<
        cp ~{reducedAnnotationDB} .
        cp ~{excludedGTF} .
        cp ~{expressedGTF} .
        cp ~{expressedKeptGTF} .

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
        File flamesGFF
        String datasetName
        Int cpu = 8
        Int memoryGB = 64
        Int diskSizeGB = 300
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    }

    String isoQuantBasename = basename(isoQuantGTF)
    String stringTieBasename = basename(stringTieGTF)
    String bambuBasename = basename(bambuGTF)
    String flairBasename = basename(flairGTF)
    String talonBasename = basename(talonGTF)
    String flamesBasename = basename(flamesGFF)

    command <<<
        cp ~{isoQuantGTF} .
        cp ~{stringTieGTF} .
        cp ~{bambuGTF} .
        cp ~{bambuGTFCounts} .
        cp ~{flairGTF} .
        cp ~{talonGTF} .
        cp ~{flamesGFF} .

        echo "~{isoQuantBasename} isoquant" > gtfs.list
        echo "~{stringTieBasename} stringtie" >> gtfs.list
        echo "~{bambuBasename} bambu" >> gtfs.list
        echo "~{flairBasename} flair" >> gtfs.list
        echo "~{talonBasename} talon" >> gtfs.list
        echo "~{flamesBasename} flames" >> gtfs.list

        /usr/local/src/IsoQuant-3.1.1/misc/denovo_model_stats.py \
        --gtf_list gtfs.list \
        --output "~{datasetName}_denovo_stats"

        tar -zcvf ~{datasetName}_denovo_stats.tar.gz ~{datasetName}_denovo_stats
    >>>

    output {
        File gffCompareOutput = "~{datasetName}_denovo_stats.tar.gz"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task ReferenceFreeGFFCompare {
    input {
        File inputGTF
        File expressedGTF
        String toolName
        String datasetName
        Int cpu = 8
        Int memoryGB = 64
        Int diskSizeGB = 300
        String docker = "us.gcr.io/broad-dsde-methods/kockan/isoquant-gffcompare:latest"
    }

    command <<<
        mkdir ~{datasetName}_~{toolName}_reffree

        cp ~{inputGTF} .
        cp ~{expressedGTF} .

        gffcompare -V -r ~{expressedGTF} -o ~{datasetName}_~{toolName}_reffree ~{inputGTF}

        ls -lha

        mv ~{datasetName}_~{toolName}_reffree* ~{datasetName}_~{toolName}_reffree

        ls -lha */*

        tar -zcvf ~{datasetName}_~{toolName}_reffree.tar.gz ~{datasetName}_~{toolName}_reffree
    >>>

    output {
        File gffCompareOutput = "~{datasetName}_~{toolName}_reffree.tar.gz"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
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
        File reducedGffCompareOutFlames
        String datasetName
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/kockan-reduced-analysis-summarize:latest"
    }

    command <<<
        cp ~{reducedGffCompareOutIsoQuant} .
        cp ~{reducedGffCompareOutStringTie} .
        cp ~{reducedGffCompareOutBambu} .
        cp ~{reducedGffCompareOutFlair} .
        cp ~{reducedGffCompareOutTalon} .
        cp ~{reducedGffCompareOutFlames} .

        tar -xzvf ~{datasetName}_isoquant_reduced_db.tar.gz
        tar -xzvf ~{datasetName}_stringtie_reduced_db.tar.gz
        tar -xzvf ~{datasetName}_bambu_reduced_db.tar.gz
        tar -xzvf ~{datasetName}_flair_reduced_db.tar.gz
        tar -xzvf ~{datasetName}_talon_reduced_db.tar.gz
        tar -xzvf ~{datasetName}_flames_reduced_db.tar.gz

        python3 /usr/local/src/plot_isoquant_results.py \
            ~{datasetName}_talon_reduced_db/talon.novel.stats,~{datasetName}_flair_reduced_db/flair.novel.stats,~{datasetName}_bambu_reduced_db/bambu.novel.stats,~{datasetName}_stringtie_reduced_db/stringtie.novel.stats,~{datasetName}_isoquant_reduced_db/isoquant.novel.stats,~{datasetName}_flames_reduced_db/flames.novel.stats \
            talon,flair,bambu,stringtie,isoquant,flames \
            ~{datasetName}
    >>>

    output {
        File reducedAnalysisSummary = "~{datasetName}_reduced_analysis_summary.tsv"
        File reducedAnalysisAccuracyPlots = "~{datasetName}_reduced_analysis_summary.png"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task ReferenceFreeAnalysisSummarize {
    input {
        File referenceFreeGffCompareOutIsoQuant
        File referenceFreeGffCompareOutStringTie
        File referenceFreeGffCompareOutIsoSeq
        File referenceFreeGffCompareOutTama
        File referenceFreeGffCompareOutCupcake
        String datasetName
        Int cpu = 1
        Int memoryGB = 32
        Int diskSizeGB = 100
        String docker = "us.gcr.io/broad-dsde-methods/kockan/kockan-reffree-analysis-summarize:latest"
    }

    command <<<
        cp ~{referenceFreeGffCompareOutIsoQuant} .
        cp ~{referenceFreeGffCompareOutStringTie} .
        cp ~{referenceFreeGffCompareOutIsoSeq} .
        cp ~{referenceFreeGffCompareOutTama} .
        cp ~{referenceFreeGffCompareOutCupcake} .

        tar -xzvf ~{datasetName}_isoquant_reffree.tar.gz
        tar -xzvf ~{datasetName}_stringtie_reffree.tar.gz
        tar -xzvf ~{datasetName}_isoseq_reffree.tar.gz
        tar -xzvf ~{datasetName}_tama_reffree.tar.gz
        tar -xzvf ~{datasetName}_cupcake_reffree.tar.gz

        ls -lha */*

        python3 /usr/local/src/plot_reffree_results.py \
        ~{datasetName}_isoquant_reffree/~{datasetName}_isoquant_reffree.stats,~{datasetName}_stringtie_reffree/~{datasetName}_stringtie_reffree.stats,~{datasetName}_isoseq_reffree/~{datasetName}_isoseq_reffree.stats,~{datasetName}_tama_reffree/~{datasetName}_tama_reffree.stats,~{datasetName}_cupcake_reffree/~{datasetName}_cupcake_reffree.stats \
        isoquant,stringtie,isoseq,tama,cupcake \
        ~{datasetName}
    >>>

    output {
        File referenceFreeAnalysisSummary = "~{datasetName}_reffree_analysis_summary.tsv"
        File referenceFreeAnalysisAccuracyPlots = "~{datasetName}_reffree_analysis_summary.png"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}