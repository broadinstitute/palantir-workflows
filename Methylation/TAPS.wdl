version 1.0

task Cutadapt {
    input {
        File fastq1
        File fastq2
        String adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        String adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        Int? cpu = 2
        Int? num_cores = 8
        Int? memory_gb = 32
        Int? disk_size_gb = ceil(3 * (size(fastq1, "GiB") + size(fastq2, "GiB")) + 50)
        String? docker = "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/taps@sha256:bbaa4049b71c592eea9ee4634f7b3ac6a3c21529db6ccd853c90fb70db5e7ca7"
    }

    String prefix1 = basename(fastq1, ".fastq.gz")
    String prefix2 = basename(fastq2, ".fastq.gz")

    command <<<
        cutadapt \
        --cores=~{num_cores} \
        -a ~{adapter1} \
        -A ~{adapter2} \
        -o ~{prefix1}.trimmed.fastq.gz \
        -p ~{prefix2}.trimmed.fastq.gz \
        ~{fastq1} \
        ~{fastq2} \
        > cutadapt.log
    >>>

    output {
        File trimmed_fastq1 = "~{prefix1}.trimmed.fastq.gz"
        File trimmed_fastq2 = "~{prefix2}.trimmed.fastq.gz"
        File cutadapt_log = "cutadapt.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: docker
    }
}

task BwaMem {
    input {
        File fastq1
        File fastq2
        File reference
        Int? cpu = 4
        Int? num_threads = 16
        Int? memory_gb = 32
        Int? disk_size_gb = ceil(3 * (size(fastq1, "GiB") + size(fastq2, "GiB")) + size(reference, "GiB") + 50)
        String? docker = "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/taps@sha256:bbaa4049b71c592eea9ee4634f7b3ac6a3c21529db6ccd853c90fb70db5e7ca7"
    }

    String prefix = basename(fastq1, ".trimmed.fastq.gz")

    command <<<
        bwa mem \
        -t ~{num_threads} \
        ~{reference} \
        ~{fastq1} \
        ~{fastq2} \
        | samtools view --threads ~{num_threads} -o ~{prefix}.trimmed.bam -
    >>>

    output {
        File bam = "~{prefix}.trimmed.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: docker
    }
}

task SortBam {
    input {
        File bam
        Int? cpu = 2
        Int? memory_gb = 32
        Int? disk_size_gb = ceil(3 * size(bam, "GiB") + 50)
        String? docker = "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/taps@sha256:bbaa4049b71c592eea9ee4634f7b3ac6a3c21529db6ccd853c90fb70db5e7ca7"
    }

    String prefix = basename(bam, ".trimmed.bam")

    command <<<
        samtools sort \
        -o ~{prefix}.trimmed.sorted.bam \
        -O bam \
        -T $~{prefix}.trimmed.bam.temp \
        ~{bam}
    >>>

    output {
        File sorted_bam = "~{prefix}.trimmed.sorted.bam"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: docker
    }
}

task MarkDuplicates {
    input {
        File bam
        Int? cpu = 2
        Int? memory_gb = 32
        Int? memory_gb_jvm = 28
        Int? disk_size_gb = ceil(3 * size(bam, "GiB") + 50)
        String? docker = "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/taps@sha256:bbaa4049b71c592eea9ee4634f7b3ac6a3c21529db6ccd853c90fb70db5e7ca7"
    }

    String prefix = basename(bam, ".trimmed.sorted.bam")

    command <<<
        python3 /gatk-4.6.2.0/gatk-4.6.2.0/gatk \
        --java-options "-Xmx~{memory_gb_jvm}g" \
        MarkDuplicates \
        --INPUT ~{bam} \
        --OUTPUT ~{prefix}.trimmed.sorted.marked.bam \
        --METRICS_FILE ~{prefix}.trimmed.sorted.marked.metrics \
        --TMP_DIR .
    >>>

    output {
        File marked_bam = "~{prefix}.trimmed.sorted.marked.bam"
        File marked_metrics = "~{prefix}.trimmed.sorted.marked.metrics"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: docker
    }
}

task RemDupsAndIndex {
    input {
        File bam
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil(3 * size(bam, "GiB") + 50)
        String? docker = "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/taps@sha256:bbaa4049b71c592eea9ee4634f7b3ac6a3c21529db6ccd853c90fb70db5e7ca7"
    }

    String prefix = basename(bam, ".trimmed.sorted.marked.bam")

    command <<<
        samtools view \
        -b \
        -f 3 \
        -F 1024 \
        -F 2048 \
        -q 20 \
        -o ~{prefix}.trimmed.sorted.marked.deduped.bam \
        ~{bam}

        samtools index ~{prefix}.trimmed.sorted.marked.deduped.bam
    >>>

    output {
        File deduped_bam = "~{prefix}.trimmed.sorted.marked.deduped.bam"
        File deduped_bam_index = "~{prefix}.trimmed.sorted.marked.deduped.bai"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: docker
    }
}

task RastairMbias {
    input {
        File bam
        File bam_index
        File reference
        Int? cpu = 4
        Int? num_threads = 8
        Int? memory_gb = 32
        Int? disk_size_gb = ceil(3 * size(bam, "GiB") + size(reference, "GiB") + 100)
        String? docker = "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/taps@sha256:bbaa4049b71c592eea9ee4634f7b3ac6a3c21529db6ccd853c90fb70db5e7ca7"
    }

    String prefix = basename(bam, ".trimmed.sorted.marked.deduped.bam")

    command <<<
        rastair mbias \
        --threads ~{num_threads} \
        --fasta-file ~{reference} \
        ~{bam} >
        ~{prefix}.trimmed.sorted.marked.deduped.rastair.mbias.txt

    >>>

    output {
        File rastair_mbias = "~{prefix}.trimmed.sorted.marked.deduped.rastair.mbias.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: docker
    }
}

task RastairCall {
    input {
        File bam
        File bam_index
        File reference
        Int? cpu = 4
        Int? memory_gb = 32
        Int? disk_size_gb = ceil(3 * size(bam, "GiB") + size(reference, "GiB") + 100)
        String? docker = "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/taps@sha256:bbaa4049b71c592eea9ee4634f7b3ac6a3c21529db6ccd853c90fb70db5e7ca7"
    }

    command <<<
        rastair call \
        --nOT 0,0,20,0 \
        --nOB 0,0,20,0 \
        --min-baseq 30 \
        --min-mapq 20 \
        --fasta-file ~{reference} \
        ~{bam} \
        >> rastair_output.mods
    >>>

    output {
        File rastair_mods = "rastair_output.mods"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: docker
    }
}

workflow TAPS {
    input {
        File fastq1
        File fastq2
        File reference
        String adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        String adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    }

    call Cutadapt {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            adapter1 = adapter1,
            adapter2 = adapter2
    }

    call BwaMem {
        input:
            fastq1 = Cutadapt.trimmed_fastq1,
            fastq2 = Cutadapt.trimmed_fastq2,
            reference = reference
    }

    call SortBam {
        input:
            bam = BwaMem.bam
    }

    call MarkDuplicates {
        input:
            bam = SortBam.sorted_bam
    }

    call RemDupsAndIndex {
        input:
            bam = MarkDuplicates.marked_bam
    }

    call RastairMbias {
        input:
            bam = RemDupsAndIndex.deduped_bam,
            bam_index = RemDupsAndIndex.deduped_bam_index,
            reference = reference
    }

    call RastairCall {
        input:
            bam = RemDupsAndIndex.deduped_bam,
            bam_index = RemDupsAndIndex.deduped_bam_index,
            reference = reference
    }

    output {
        File deduped_bam = RemDupsAndIndex.deduped_bam
        File deduped_bam_index = RemDupsAndIndex.deduped_bam_index
        File rastair_mbias = RastairMbias.rastair_mbias
        File rastair_mods = RastairCall.rastair_mods
    }
}