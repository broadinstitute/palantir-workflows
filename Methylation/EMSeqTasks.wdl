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
        ls -lha

        bwameth.py --reference ~{refBasename} --threads ~{numThreads} --read-group "@RG\\tID:~{sampleId}\\tSM:~{sampleId}" ~{sampleId}.1.filtered.fastq.gz ~{sampleId}.2.filtered.fastq.gz > ~{sampleId}.sam
        ls -lha

        /usr/local/src/mark-nonconverted-reads/mark-nonconverted-reads.py --bam ~{sampleId}.sam --out ~{sampleId}.nc_marked.sam 2> ~{sampleId}.nonconverted.tsv
        ls -lha

        sambamba view --with-header --sam-input --nthreads 2 --format bam --compression-level 0 --output-filename "~{sampleId}.bam" ~{sampleId}.nc_marked.sam
        ls -lha
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