version 1.0

workflow ScoreBGE {
    input {
        File exome_gvcf
        File exome_gvcf_index
        File imputed_wgs_vcf
        File imputed_wgs_vcf_index
        String basename
        File weights

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
            preemptible = preemptible
    }

    output {
        File exome_gvcf_score = ScoreGvcfAndVcf.exome_gvcf_score
        File imputed_wgs_vcf_score = ScoreGvcfAndVcf.imputed_wgs_vcf_score
        File score = ScoreGvcfAndVcf.score
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

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int preemptible = 1
        Int? disk_gb
        Int mem_gb = 4
        Int cpu = 4
    }

    Int disk_size_gb = select_first([disk_gb, ceil(size(exome_gvcf, "GiB") + size(imputed_wgs_vcf, "GiB") + size(weights, "GiB") + size(ref_fasta, "GiB") + 50)])

    command <<<
        set -xeuo pipefail
        
        cat <<'EOF' > script.py
import sys
sys.path.insert(0, '/ScoreBGE')
import ScoreBGE

bge_scorer = ScoreBGE.BGEScorer('~{ref_dict}', '~{weights}')
bge_scorer.score_wes_gvcf('~{exome_gvcf}')
bge_scorer.score_wgs_vcf('~{imputed_wgs_vcf}')
bge_scorer.write_output('~{basename}')
EOF
            python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/palantir-workflows-score-bge:dev"
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File exome_gvcf_score = "~{basename}.exome_gvcf.score"
        File imputed_wgs_vcf_score = "~{basename}.imputed_wgs_vcf.score"
        File score = "~{basename}.score"
    }
}
