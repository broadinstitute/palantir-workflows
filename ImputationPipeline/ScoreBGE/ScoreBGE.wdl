version 1.0

workflow ScoreBGE {
    input {
        File exome_gvcf
        File exome_gvcf_index
        File imputed_wgs_vcf
        File imputed_wgs_vcf_index
        String basename
        File weights
        Array[String]? sample_names

        String? score_bge_docker

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int preemptible = 1
    }

    call ScoreGvcfAndVcf {
        input:
            exome_gvcf = exome_gvcf,
            exome_gvcf_index = exome_gvcf_index,
            imputed_wgs_vcf = imputed_wgs_vcf,
            imputed_wgs_vcf_index = imputed_wgs_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            basename = basename,
            weights = weights,
            sample_names = sample_names,
            score_bge_docker = score_bge_docker,
            preemptible = preemptible
    }

    output {
        File exome_gvcf_score = ScoreGvcfAndVcf.exome_gvcf_score
        File imputed_wgs_vcf_score = ScoreGvcfAndVcf.imputed_wgs_vcf_score
        File score = ScoreGvcfAndVcf.score
        File exome_gvcf_sites_scored = ScoreGvcfAndVcf.exome_gvcf_sites_scored
        File imputed_wgs_vcf_sites_scored = ScoreGvcfAndVcf.imputed_wgs_vcf_sites_scored
    }
}

task ScoreGvcfAndVcf {
    input {
        File exome_gvcf
        File exome_gvcf_index
        File imputed_wgs_vcf
        File imputed_wgs_vcf_index
        String basename
        File weights
        Array[String]? sample_names

        String score_bge_docker = "us.gcr.io/broad-dsde-methods/palantir-workflows-score-bge:palantir-workflows_5feb024"

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int preemptible = 1
        Int? disk_gb
        Int mem_gb = 4
        Int cpu = 4
    }

    String sample_names_arg = if defined(sample_names) then "--sample-names" else ""

    Int disk_size_gb = select_first([disk_gb, ceil(size(exome_gvcf, "GiB") + size(imputed_wgs_vcf, "GiB") + size(weights, "GiB") + size(ref_fasta, "GiB") + 50)])

    command <<<
        set -xeuo pipefail
        python3 /ScoreBGE.py --ref-dict ~{ref_dict} --weights ~{weights} --gvcf ~{exome_gvcf} --vcf ~{imputed_wgs_vcf} \
            --basename ~{basename} ~{sample_names_arg} ~{sep=" --sample-names " sample_names}
    >>>

    runtime {
        docker: score_bge_docker
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File exome_gvcf_score = "~{basename}.exome_gvcf.score"
        File imputed_wgs_vcf_score = "~{basename}.imputed_wgs_vcf.score"
        File score = "~{basename}.score"
        File exome_gvcf_sites_scored = "~{basename}.exome_gvcf.sites_scored"
        File imputed_wgs_vcf_sites_scored = "~{basename}.imputed_wgs_vcf.sites_scored"
    }
}
