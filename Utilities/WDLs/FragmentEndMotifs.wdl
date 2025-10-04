version 1.0

task CramToBam {
    input {
        String sample_id
        File cram
        File reference_fasta
        File reference_fai
        File reference_dict

        Int? cpu = 4
        Int? mem_gb = 32
        Int? disk_size_gb = ceil((6 * size(cram, "GiB")) + size(reference_fasta, "GiB") + 32)
    }

    command <<<
        set -e
        set -o pipefail

        ln -vs ~{reference_fasta} reference.fasta
        ln -vs ~{reference_fai} reference.fasta.fai
        ln -vs ~{reference_dict} reference.dict

        samtools view -h -T reference.fasta ~{cram} |
        samtools view -b -o ~{sample_id}.bam -
        samtools index -b ~{sample_id}.bam
        mv ~{sample_id}.bam.bai ~{sample_id}.bai
    >>>

    output {
        File output_bam = "~{sample_id}.bam"
        File output_bam_index = "~{sample_id}.bai"
    }

    runtime {
        cpu: cpu
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
    }
}

task ExtractFragmentEndMotifs {
    input {
        String sample_id
        File bam
        File bai
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
        tumorbainame=$(basename ~{bai})

        ln -s ~{bam} /processing_dir/$tumorbamname
        ln -s ~{bai} /processing_dir/$tumorbamname.bai

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
        File? cram
        File? bam
        File? bai
        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict
        Int kmer_length
        Int min_qual = 20
    }

    if(defined(cram)) {
        call CramToBam {
            input:
                sample_id = sample_id,
                cram = select_first([cram]),
                reference_fasta = ref_fasta,
                reference_fai = ref_fasta_index,
                reference_dict = ref_fasta_dict
        }
    }

    call ExtractFragmentEndMotifs {
        input:
            sample_id = sample_id,
            bam = select_first([bam, CramToBam.output_bam]),
            bai = select_first([bai, CramToBam.output_bam_index]),
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            min_qual = min_qual,
            kmer_length = kmer_length
    }

    output {
        File fragment_end_motif_freq = ExtractFragmentEndMotifs.fragment_end_motif_freq
    }
}