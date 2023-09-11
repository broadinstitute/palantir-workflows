version 1.0

task Mapping {
    input {
        String sampleId
        String flowcell
        String library
        String lane = "all"
        String tile = "all"
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

        seqtk mergepe <(zcat -f ~{fq1}) <(zcat -f ~{fq2}) \
        | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis \
            -j ~{library}_fastp.json 2> fastp.stderr \
        | bwameth.py -p -t ~{numThreads} --read-group "@RG\\tID:${fastq_barcode}\\tSM:~{library}" --reference ~{refBasename} /dev/stdin \
            2>  "~{library}_${fastq_barcode}_~{flowcell}_~{lane}_~{tile}.log.bwamem" \
        | /usr/local/src/mark-nonconverted-reads/mark-nonconverted-reads.py 2> "~{library}_${fastq_barcode}_~{flowcell}_~{lane}_~{tile}.nonconverted.tsv" \
        | sambamba view -t 2 -S -f bam -o "~{library}_${fastq_barcode}_~{flowcell}_~{lane}_~{tile}.aln.bam" /dev/stdin 2> sambamba.stderr;

        ls -lha
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