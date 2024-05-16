version 1.0

import "ScoringTasks.wdl"

workflow ScoreVcfWithPreferredGvcf {
    input {
        File exome_gvcf
        File exome_gvcf_index
        File imputed_wgs_vcf
        String basename
        File weights
        File? weights_as_vcf
        Int? base_mem
        String? extra_args
        String chromosome_encoding = "chrMT"
        Boolean use_ref_alt_for_ids = true

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int preemptible = 1
    }

    if (!defined(weights_as_vcf)) {
        call ConvertWeightsTsvToVcf {
            input:
                weights_tsv = weights
        }
    }

    call ExtractSitesFromGvcf {
        input:
            gvcf = exome_gvcf,
            gvcf_index = exome_gvcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            basename = basename,
            sites_to_extract = select_first([weights_as_vcf, ConvertWeightsTsvToVcf.weights_vcf])
    }
  
    call ScoringTasks.ScoreVcf as ExomeGVCFScore {
        input:
            vcf = ExtractSitesFromGvcf.vcf,
            basename = basename,
            weights = weights,
            sites = ExtractSitesFromGvcf.extracted_sites_from_gvcf,
            base_mem = base_mem,
            extra_args = extra_args,
            chromosome_encoding = chromosome_encoding,
            use_ref_alt_for_ids = use_ref_alt_for_ids,
            use_dosage_annotation = false
    }

    call ScoringTasks.ScoreVcf as ImputedWGSVCFScore {
        input:
            vcf = imputed_wgs_vcf,
            basename = basename,
            weights = weights,
            exclude_sites = ExtractSitesFromGvcf.extracted_sites_from_gvcf,
            base_mem = base_mem,
            extra_args = extra_args,
            chromosome_encoding = chromosome_encoding,
            use_ref_alt_for_ids = use_ref_alt_for_ids,
            use_dosage_annotation = true
    }

    call CombinePrimaryAndSecondaryScores {
        input:
            primary_score = ExomeGVCFScore.score,
            secondary_score = ImputedWGSVCFScore.score,
            primary_sites_scored = ExomeGVCFScore.sites_scored,
            secondary_sites_scored = ImputedWGSVCFScore.sites_scored,
            basename = basename,
            preemptible = preemptible
    }

    output {
        File score = CombinePrimaryAndSecondaryScores.score
        File sites_scored = CombinePrimaryAndSecondaryScores.sites_scored
    }
}

task ConvertWeightsTsvToVcf {
    input {
        File weights_tsv
        Int mem_gb = 4
        Int cpu = 4
        Int? disk_gb
        Int preemptible = 1
    }

    Int disk_size_gb = select_first([disk_gb, ceil(size(weights_tsv, "GiB") * 3 + 50)])
    String weights_basename = basename(weights_tsv)

    command <<<
        set -xeuo pipefail

        echo -e '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > ~{weights_basename}.vcf
        tail -n +2 ~{weights_tsv} | awk -F'\t' 'BEGIN {OFS="\t"} {split($1, a, ":"); print a[1], a[2], ".", a[3], a[4], ".", ".", "."}' >> ~{weights_basename}.vcf
    >>>

    runtime {
        docker: "ubuntu:latest"
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File weights_vcf = "~{weights_basename}.vcf"
    }
}

task ExtractSitesFromGvcf {
    input {
        File gvcf
        File gvcf_index
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String basename
        File sites_to_extract
        Int mem_gb = 4
        Int cpu = 4
        Int? disk_gb
        Int preemptible = 1
        String gatk_tag = "4.5.0.0"
    }

    Int disk_size_gb = select_first([disk_gb, ceil(size(gvcf, "GiB") * 2 + size(ref_fasta, "GiB") + 50)])

    command <<<
        set -xeuo pipefail

        gatk SelectVariants \
            -V ~{gvcf} \
            -L ~{sites_to_extract} \
            -O ~{basename}.high_quality.vcf.gz \
            -select "(vc.hasAttribute(\"END\") && vc.getGenotype(\"~{basename}\").getGQ() > 30) || QUAL > 30"
        
        gatk GenotypeGVCFs \
            -R ~{ref_fasta} \
            -V ~{basename}.high_quality.vcf.gz \
            -O ~{basename}.extracted_from_gvcf.vcf.gz \
            --force-output-intervals ~{sites_to_extract}
        
        bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ~{basename}.extracted_from_gvcf.vcf.gz > ~{basename}.extracted_from_gvcf.sites

    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File vcf = "~{basename}.extracted_from_gvcf.vcf.gz"
        File vcf_index = "~{basename}.extracted_from_gvcf.vcf.gz"
        File extracted_sites_from_gvcf = "~{basename}.extracted_from_gvcf.sites"
    }
}

task CombinePrimaryAndSecondaryScores {
    input {
        File primary_score
        File secondary_score
        File primary_sites_scored
        File secondary_sites_scored

        String basename

        Int disk_gb = 50
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 1
    }

    command <<<
        set -xeuo pipefail

        cat ~{primary_sites_scored} ~{secondary_sites_scored} | sort | uniq > ~{basename}.sscore.vars

        cat <<EOF > script.py
import pandas as pd

primary_scores = pd.read_table('~{primary_score}')
secondary_scores = pd.read_table('~{secondary_score}')

scores = primary_scores.merge(secondary_scores, how='inner', on='#IID')
scores['SCORE1_SUM'] = scores['SCORE1_SUM_x'] + scores['SCORE1_SUM_y']
scores['NAMED_ALLELE_DOSAGE_SUM'] = scores['NAMED_ALLELE_DOSAGE_SUM_x'] + scores['NAMED_ALLELE_DOSAGE_SUM_y']
scores['SCORE1_AVG'] = 'NA'
scores[['#IID', 'NAMED_ALLELE_DOSAGE_SUM', 'SCORE1_AVG', 'SCORE1_SUM']].to_csv('~{basename}.sscore', sep='\t', index=None)
EOF
        python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + disk_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File score = "~{basename}.sscore"
        File sites_scored = "~{basename}.sscore.vars"
    }
}