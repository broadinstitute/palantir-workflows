version 1.0

task Mapping {
    input {
        String sampleId
        File fq1
        File fq2
        File ref
        File refIdx
        File bwamethIdx
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    String refBasename = basename(ref)
    String bwamethIdxBasename = basename(bwamethIdx)
    String bwamethIdxBaseDir = basename(bwamethIdx, ".tar.gz")

    command <<<
        mv ~{ref} ~{refIdx} .
        mv ~{bwamethIdx} .
        tar -xzvf ~{bwamethIdxBasename}
        mv ~{bwamethIdxBaseDir}/* .

        inst_name=$(zcat -f ~{fq1} | head -n 1 | cut -f 1 -d ':' | sed 's/^@//')
        fastq_barcode=$(zcat -f ~{fq2} | head -n 1 | sed -r 's/.*://')

        if [[ "${inst_name:0:2}" == 'A0' ]] || [[ "${inst_name:0:2}" == 'NS' ]] || \
           [[ "${inst_name:0:2}" == 'NB' ]] || [[ "${inst_name:0:2}" == 'VH' ]] ; then
            trim_polyg='--trim_poly_g'
            echo '2-color instrument: poly-g trim mode on'
        else
            trim_polyg=''
        fi

        fastp --in1 ~{fq1} --out1 ~{sampleId}.1.filtered.fastq.gz --in2 ~{fq2} --out2 ~{sampleId}.2.filtered.fastq.gz -l 2 -Q ${trim_polyg} --overrepresentation_analysis -h ~{sampleId}_fastp.html
        bwameth.py --reference ~{refBasename} --threads ~{numThreads} --read-group "@RG\\tID:~{sampleId}\\tSM:~{sampleId}" ~{sampleId}.1.filtered.fastq.gz ~{sampleId}.2.filtered.fastq.gz > ~{sampleId}.sam
        /usr/local/src/mark-nonconverted-reads-1.1/mark-nonconverted-reads.py --bam ~{sampleId}.sam --out ~{sampleId}.nc_marked.sam 2> ~{sampleId}.nonconverted.tsv
        sambamba view --with-header --sam-input --nthreads 2 --format bam --compression-level 0 --output-filename "~{sampleId}.bam" ~{sampleId}.nc_marked.sam
    >>>

    output {
        File bam = "~{sampleId}.bam"
        File fastpReport = "~{sampleId}_fastp.html"
        File nonconvertedReadCounts = "~{sampleId}.nonconverted.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task MarkDuplicates {
    input {
        String sampleId
        File bam
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    command <<<
        samtools view -h ~{bam} \
        | samblaster 2> ~{sampleId}.log.samblaster \
        | sambamba view -t 2 -l 0 -S -f bam /dev/stdin \
        | sambamba sort --nthreads ~{numThreads} --memory-limit 32GB --out ~{sampleId}.md.bam /dev/stdin
    >>>

    output {
        File mdBam = "~{sampleId}.md.bam"
        File mdBai = "~{sampleId}.md.bam.bai"
        File samblasterLog = "~{sampleId}.log.samblaster"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task MethylDackelMbias {
    input {
        String sampleId
        File bam
        File bai
        File ref
        File refIdx
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    String refBasename = basename(ref)

    command <<<
        mv ~{ref} ~{refIdx} .

        echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > ~{sampleId}_combined_mbias.tsv
        chrs=(`samtools view -H ~{bam} | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

        for chr in ${chrs[*]}; do
            for context in CHH CHG CpG; do
                arg=''
                if [ $context = 'CHH' ]; then
                    arg='--CHH --noCpG'
                elif [ $context = 'CHG' ]; then
                    arg='--CHG --noCpG'
                fi

                # Need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
                # Not sure why we need both --keepDupes and -F, probably a bug in mbias
                join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
                <( \
                    MethylDackel mbias --noSVG $arg -@ ~{numThreads} -r $chr ~{refBasename} ~{bam} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                <( \
                    MethylDackel mbias --noSVG --keepDupes -F 2816 $arg -@ ~{numThreads} -r $chr ~{refBasename} ~{bam} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                | sed "s/^/${chr}\t${context}\t/" \
                >> ~{sampleId}_combined_mbias.tsv
            done
        done

        # Makes the svg files for trimming checks
        MethylDackel mbias -@ !{task.cpus} --noCpG --CHH --CHG -r ${chrs[0]} !{genome} !{md_file} !{library}_chn
        for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CHN <\\/text>/" $f; done;

        MethylDackel mbias -@ !{task.cpus} -r ${chrs[0]} !{genome} !{md_file} !{library}_cpg
        for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CpG<\\/text>/" $f; done;
    >>>

    output {

    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task MethylDackelExtract {
    input {
        String sampleId
        File bam
        File bai
        File ref
        File refIdx
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    String refBasename = basename(ref)

    command <<<
        mv ~{ref} ~{refIdx} .
        MethylDackel extract --methylKit --nOT 0,0,0,5 --nOB 0,0,5,0 -@ ~{numThreads} --CHH --CHG -o ~{sampleId} ~{refBasename} ~{bam}
        pigz -p ~{numThreads} *.methylKit
        ls -lha
    >>>

    output {
        File methylKit = "~{sampleId}.methylKit.gz"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task SelectHumanReads {
    input {
        String sampleId
        File bam
        File bai
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    command <<<
        sambamba view -t 8 -l 0 -f bam ~{bam} \
            chr1 chr2 chr3 chr4 chr5 chr6 \
            chr7 chr8 chr9 chr10 chr11 chr12 \
            chr13 chr14 chr15 chr16 chr17 chr18 \
            chr19 chr20 chr21 chr22 chrX chrY \
        > ~{sampleId}.human.bam
        ls -lha
    >>>

    output {
        File humanBam = "~{sampleId}.human.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task FastQC {
    input {
        String sampleId
        File bam
        File bai
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    command <<<
        fastqc -f bam ~{bam}
        ls -lha
    >>>

    output {
        File fastqcResults = "~{sampleId}_fastqc.zip"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task SamtoolsFlagstats {
    input {
        String sampleId
        File bam
        File bai
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    command <<<
        samtools flagstat -@ ~{numThreads} ~{bam} > ~{sampleId}.flagstat
        samtools idxstats ~{bam} > ~{sampleId}.idxstat
    >>>

    output {
        File flagstat = "~{sampleId}.flagstat"
        File idxstat = "~{sampleId}.idxstat"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task SamtoolsStats {
    input {
        String sampleId
        File bam
        File bai
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    command <<<
        samtools stats -@ ~{numThreads} ~{bam} > ~{sampleId}.samstat
    >>>

    output {
        File samstat = "~{sampleId}.samstat"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task PicardGCBias {
    input {
        String sampleId
        File bam
        File bai
        File ref
        File refIdx
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    String refBasename = basename(ref)

    command <<<
        mv ~{ref} ~{refIdx} .
        picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=LENIENT I=~{bam} O=~{sampleId}.gc_metrics S=~{sampleId}.gc_summary_metrics CHART=~{sampleId}.gc.pdf R=~{refBasename}
    >>>

    output {
        File gcMetrics = "~{sampleId}.gc_metrics"
        File gcSummaryMetrics = "~{sampleId}.gc_summary_metrics"
        File gcChart = "~{sampleId}.gc.pdf"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task PicardStats {
    input {
        String sampleId
        File bam
        File bai
        File ref
        File refIdx
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    String refBasename = basename(ref)

    command <<<
        mv ~{ref} ~{refIdx} .
        picard -Xmx16g CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT  I=~{bam} O=~{sampleId}.insertsize_metrics MINIMUM_PCT=0 HISTOGRAM_FILE=/dev/null
    >>>

    output {
        File insertSizeMetrics = "~{sampleId}.insertsize_metrics"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

task GoLeft {
    input {
        String sampleId
        File bam
        File bai
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    command <<<
        goleft indexcov --directory . ~{bam}
    >>>

    output {
        File goleftPed = "~{sampleId}-indexcov.ped"
        File goleftRoc = "~{sampleId}-indexcov.roc"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}