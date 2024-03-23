version 1.0

import "ScoringTasks.wdl"

workflow ScoreVcfWithPreferredGvcf {
    input {
        File preferred_vcf
        File secondary_vcf
        String basename
        File weights
        File? weights_as_vcf
        Int? base_mem
        String? extra_args
        String? chromosome_encoding
        Boolean? use_ref_alt_for_ids

        File ref_fasta
        File ref_fasta_index

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
            gvcf = preferred_vcf,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            basename = basename,
            sites_to_extract = select_first([weights_as_vcf, ConvertWeightsTsvToVcf.weights_vcf])
    }

    call FilterSitesAndGetExtractedSiteIDs {
        input:
            extracted_vcf = ExtractSitesFromGvcf.vcf,
            basename = basename
    }
  
    call ScoringTasks.ScoreVcf as PrimaryScore {
        input:
            vcf = FilterSitesAndGetExtractedSiteIDs.extracted_filtered_vcf,
            basename = basename,
            weights = weights,
            sites = FilterSitesAndGetExtractedSiteIDs.extracted_filtered_sites,
            base_mem = base_mem,
            extra_args = extra_args,
            chromosome_encoding = chromosome_encoding,
            use_ref_alt_for_ids = use_ref_alt_for_ids
    }

    call ScoringTasks.ScoreVcf as SecondaryScore {
        input:
            vcf = secondary_vcf,
            basename = basename,
            weights = weights,
            exclude_sites = FilterSitesAndGetExtractedSiteIDs.extracted_filtered_sites,
            base_mem = base_mem,
            extra_args = extra_args,
            chromosome_encoding = chromosome_encoding,
            use_ref_alt_for_ids = use_ref_alt_for_ids
    }

    call CombinePrimaryAndSecondaryScores {
        input:
            primary_score = PrimaryScore.score,
            secondary_score = SecondaryScore.score,
            primary_sites_scored = PrimaryScore.sites_scored,
            secondary_sites_scored = SecondaryScore.sites_scored,
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
        Int preemptible = 1
    }

    String weights_basename = basename(weights_tsv)

    command <<<
        set -xeuo pipefail

        echo -e '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > ~{weights_basename}.vcf
        awk -F'\t' 'BEGIN {OFS="\t"} {split($1, a, ":"); print a[1], a[2], ".", a[3], a[4], ".", ".", "."}' ~{weights_tsv} >> ~{weights_basename}.vcf
    >>>

    runtime {
        docker: "ubuntu:latest"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File weights_vcf = "~{weights_basename}.vcf"
    }
}

task ExtractSitesFromGvcf {
    input {
        File gvcf
        File ref_fasta
        File ref_fasta_index
        String basename
        File sites_to_extract
        Int mem_gb = 4
        Int cpu = 4
        Int preemptible = 1
        String gatk_tag = "4.5.0.0"
    }

    command <<<
        set -xeuo pipefail

        gatk GenotypeGVCFs \
            -R ~{ref_fasta} \
            -V ~{gvcf} \
            -O ~{basename}.extracted.vcf.gz \
            --force-output-intervals ~{sites_to_extract}
    >>>

    runtime {
        docker: "broadinstitute/gatk:" + gatk_tag
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File vcf = "~{basename}.extracted.vcf.gz"
        File vcf_index = "~{basename}.extracted.vcf.gz"
    }
}

task FilterSitesAndGetExtractedSiteIDs {
    input {
        File extracted_vcf
        String basename
        Int mem_gb = 4
        Int cpu = 4
        Int preemptible = 1
    }

    command <<<
        set -xeuo pipefail

        bcftools view -O z -o ~{basename}.extracted.filtered.vcf.gz -e 'QUAL<=20' ~{extracted_vcf}
        bcftools index ~{basename}.extracted.filtered.vcf.gz

        bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ~{basename}.extracted.filtered.vcf.gz > ~{basename}.extracted.filtered.sites
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File extracted_filtered_vcf = "~{basename}.extracted.filtered.vcf.gz"
        File extracted_filtered_vcf_index = "~{basename}.extracted.filtered.vcf.gz"
        File extracted_filtered_sites = "~{basename}.extracted.filtered.sites"
    }
}

task CombinePrimaryAndSecondaryScores {
    input {
        File primary_score
        File secondary_score
        File primary_sites_scored
        File secondary_sites_scored

        String basename

        Int disk_size_gb = 50
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
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File score = "~{basename}.sscore"
        File sites_scored = "~{basename}.sscore.vars"
    }
}