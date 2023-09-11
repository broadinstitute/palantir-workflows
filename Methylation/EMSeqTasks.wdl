version 1.0
    
task Mapping {
    input {
        String sampleId
        File fq1
        File fq2
        # ile ref
        #File refIdx
        #File bwamethIdx
        Int cpu = 1
        Int numThreads = 8
        Int memoryGB = 16
        Int diskSizeGB = 128
        String docker = "us.gcr.io/broad-dsde-methods/kockan/em-seq:latest"
    }

    #String refBasename = basename(ref)
    #String bwamethIdxBasename = basename(bwamethIdx)
    #String bwamethIdxBaseDir = basename(bwamethIdx, ".tar.gz")

    # command parts that are temporarily commented out
    # mv ~{ref} ~{refIdx} .
    # mv ~{bwamethIdx} .
    # tar -xzvf ~{bwamethIdxBasename}
    # mv ~{bwamethIdxBaseDir}/* .

    # bwameth.py \
        #--reference ~{refBasename} \
        #--threads ~{numThreads} \
        #--read-group '@RG\tID:~{sampleId}\tPL:illumina\tSM:~{sampleId}\tLB:~{sampleId}' \
        #~{fq1} ~{fq2} > ~{sampleId}.sam

    #if [[ "${inst_name:0:2}" == 'A0' ]] || [[ "${inst_name:0:2}" == 'NS' ]] || \
    #        [[ "${inst_name:0:2}" == 'NB' ]] || [[ "${inst_name:0:2}" == 'VH' ]] ;
    #    then
    #        trim_polyg='--trim_poly_g'
    #        echo '2-color instrument: poly-g trim mode on'
    #    else
    #        trim_polyg=''
    #    fi

    #    seqtk mergepe <(zcat -f "!{fq1}") <(zcat -f "!{fq2}") \
    #    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis \
    #        -j "!{fq_set.library}_fastp.json" 2> fastp.stderr \
    #    | bwameth.py -p -t !{task.cpus} --read-group "@RG\\tID:${fastq_barcode}\\tSM:!{fq_set.library}" --reference !{genome} /dev/stdin \
    #             2>  "!{fq_set.library}_${fastq_barcode}!{fq_set.flowcell}_!{fq_set.lane}_!{fq_set.tile}.log.bwamem" \
    #    | mark-nonconverted-reads.py 2> "!{fq_set.library}_${fastq_barcode}_!{fq_set.flowcell}_!{fq_set.lane}_!{fq_set.tile}.nonconverted.tsv" \
    #    | sambamba view -t 2 -S -f bam -o "!{fq_set.library}_${fastq_barcode}_!{fq_set.flowcell}_!{fq_set.lane}_!{fq_set.tile}.aln.bam" /dev/stdin 2> sambamba.stderr;

    command <<<
        inst_name=$(zcat -f ~{fq1} | head -n 1 | cut -f 1 -d ':' | sed 's/^@//')
        echo $inst_name

        fastq_barcode=$(zcat -f ~{fq2} | head -n 1 | sed -r 's/.*://')
        echo $fastq_barcode
    >>>

    output {
        # File sam = "~{sampleId}.sam"
        # set val(fq_set.library), file("*.aln.bam") into aligned_files
        # set val(fq_set.library), file("*.nonconverted.tsv") into nonconverted_counts
        # set val(fq_set.library), file("*_fastp.json") into fastp_log_files
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}