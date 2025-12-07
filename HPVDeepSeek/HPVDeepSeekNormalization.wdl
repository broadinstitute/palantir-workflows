version 1.0

task NormalizeHPVCounts {
    input {
        File bam
        File bai
        String top_hpv_contig
        Float top_hpv_num_reads
        Float ml_plasma
        Float ng_cfdna

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = ceil((2 * size(bam, "GiB")) + 50)
    }

    command <<<
        gapdh_num_reads=$(samtools coverage ~{bam} -r chr12:6534321-6538610 | awk 'NR > 1 {print int($4)}')
        cthpvdna_per_human_genome_equivalents=$(python -c "print(~{top_hpv_num_reads} / ($gapdh_num_reads / 2))")
        cthpvdna_per_ml_plasma=$(python -c "print(~{top_hpv_num_reads} / ~{ml_plasma})")
        cthpvdna_per_ng_cfdna=$(python -c "print(~{top_hpv_num_reads} / ~{ng_cfdna})")
        cthpvdna_per_human_genome_equivalents_per_ml_plasma=$(python -c "print($cthpvdna_per_human_genome_equivalents / ~{ml_plasma})")

        echo $gapdh_num_reads > gapdh_num_reads.txt
        echo $cthpvdna_per_human_genome_equivalents > cthpvdna_per_human_genome_equivalents.txt
        echo $cthpvdna_per_ml_plasma > cthpvdna_per_ml_plasma.txt
        echo $cthpvdna_per_ng_cfdna > cthpvdna_per_ng_cfdna.txt
        echo $cthpvdna_per_human_genome_equivalents_per_ml_plasma > cthpvdna_per_human_genome_equivalents_per_ml_plasma.txt
    >>>

    output {
        Int gapdh_num_reads = read_int("gapdh_num_reads.txt")
        Float cthpvdna_per_human_genome_equivalents = read_float("cthpvdna_per_human_genome_equivalents.txt")
        Float cthpvdna_count_per_ml_plasma = read_float("cthpvdna_per_ml_plasma.txt")
        Float cthpvdna_count_per_ng_cfdna = read_float("cthpvdna_per_ng_cfdna.txt")
        Float cthpvdna_per_human_genome_equivalents_per_ml_plasma = read_float("cthpvdna_per_human_genome_equivalents_per_ml_plasma.txt")
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

workflow HPVDeepSeekNormalization {
    input {
        File bam
        File bai
        String top_hpv_contig
        Float top_hpv_num_reads
        Float ml_plasma
        Float ng_cfdna
    }

    call NormalizeHPVCounts {
        input:
            bam = bam,
            bai = bai,
            top_hpv_contig = top_hpv_contig,
            top_hpv_num_reads = top_hpv_num_reads,
            ml_plasma = ml_plasma,
            ng_cfdna = ng_cfdna
    }

    output {
        Int gapdh_num_reads = NormalizeHPVCounts.gapdh_num_reads
        Float cthpvdna_per_human_genome_equivalents = NormalizeHPVCounts.cthpvdna_per_human_genome_equivalents
        Float cthpvdna_count_per_ml_plasma = NormalizeHPVCounts.cthpvdna_count_per_ml_plasma
        Float cthpvdna_count_per_ng_cfdna = NormalizeHPVCounts.cthpvdna_count_per_ng_cfdna
        Float cthpvdna_per_human_genome_equivalents_per_ml_plasma = NormalizeHPVCounts.cthpvdna_per_human_genome_equivalents_per_ml_plasma
    }
}