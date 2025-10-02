version 1.0

task ExtractFragmentEndMotifs {
    input {
        String sample_id
        File bam
        File bam_index
        File ref_fasta
        File ref_fasta_index
        Int kmer_length
        Int? min_qual

        Int? cpu = 8
        Int? memory_gb = 32
        Int? disk_size_gb = ceil((2.5 * size(bam, "GiB")) + size(ref_fasta, "GiB") + 32)
    }

    command <<<
        set -euxo pipefail

        mkdir results
        mkdir /processing_dir

        tumorbamname=$(basename ~{bam})
        tumorbainame=$(basename ~{bam_index})

        ln -s ~{bam} /processing_dir/$tumorbamname
        ln -s ~{bam_index} /processing_dir/$tumorbamname.bai

        python /deepTools_fragment_motif.py \
        --bam_file /processing_dir/$tumorbamname \
        --sample_name ~{sample_id} \
        --ref_seq ~{ref_fasta} \
        --map_q ~{min_qual} \
        --CPU ~{cpu} \
        --kmer ~{kmer_length} \
        --outfolder ./results/

        mv ./results/~{sample_id}.fragment_end_motif_freq.txt .
    >>>

    output {
        File fragment_end_motif_freq = "~{sample_id}.fragment_end_motif_freq.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "cricker3/fragmentomics:1.0"
    }
}

workflow FragmentEndMotifs {
    input {
        String sample_id
        File bam
        File bam_index
        File ref_fasta
        File ref_fasta_index
        Int kmer_length
        Int min_qual = 20
    }

    call ExtractFragmentEndMotifs {
        input:
            sample_id = sample_id,
            bam = bam,
            bam_index = bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            min_qual = min_qual,
            kmer_length = kmer_length
    }

    output {
        File fragment_end_motif_freq = ExtractFragmentEndMotifs.fragment_end_motif_freq
    }
}