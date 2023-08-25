version 1.0

workflow PhaseVCF {
    input {
        File input_vcf
        File input_vcf_index
        String? input_sample_name

        File input_bam
        File input_bam_index

        File reference_fasta
        File reference_fasta_index
    }

    call WhatsHapPhase {
        input:
            input_vcf=input_vcf,
            input_vcf_index=input_vcf_index,
            input_sample_name=input_sample_name,
            input_bam=input_bam,
            input_bam_index=input_bam_index,
            reference_fasta=reference_fasta,
            reference_fasta_index=reference_fasta_index
    }

    output {
        File phased_vcf = WhatsHapPhase.phased_vcf
        File phased_vcf_index = WhatsHapPhase.phased_vcf_index
        File phase_stats = WhatsHapPhase.phase_stats
    }
}

task WhatsHapPhase {
    input {
        File input_vcf
        File input_vcf_index
        String? input_sample_name    # Required if ignore_read_groups false && multi-sample VCF

        File input_bam
        File input_bam_index

        File reference_fasta
        File reference_fasta_index

        # Tool arguments
        Boolean ignore_read_groups = true    # Switch to false to phase by Read Groups (RG)
        Boolean distrust_genotypes = false    # Turn on to allow hets -> hom if it allows optimal phasing solution
        File? pedigree
        String tag = "PS"    # Other option: HP
        Float? recombrate     # Recombination rate to use (1.26 for humans); change for non-human samples

        # Runtime arguments
        Int disk_size = ceil(2 * size(input_vcf, "GB") + size(input_bam, "GB") + size(reference_fasta, "GB")) + 50
        Int cpu = 8
        Int memory_ram = 64
    }

    command <<<
        set -xueo pipefail

        whatshap phase \
            --reference=~{reference_fasta} \
            -o phased.vcf.gz \
            --tag=~{tag} \
            ~{"--ped" + pedigree} \
            ~{true="--ignore-read-groups" false="" ignore_read_groups} \
            ~{true="--distrust-genotypes" false="" distrust_genotypes} \
            ~{"--recombrate" + recombrate} \
            ~{input_vcf} \
            ~{input_bam}

        tabix phased.vcf.gz
        whatshap stats --tsv="phase_stats.tsv" phased.vcf.gz

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/whatshap:v1"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory_ram + " GB"
    }

    output {
        File phased_vcf = "phased.vcf.gz"
        File phased_vcf_index = "phased.vcf.gz.tbi"
        File phase_stats = "phase_stats.tsv"
    }
}
