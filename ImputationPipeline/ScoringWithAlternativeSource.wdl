version 1.0

import "ScoringTasks.wdl"

workflow ScoreVcfWithSecondarySource {
    input {
        File primary_vcf
        File? secondary_vcf
        String basename
        File primary_weights
        File? secondary_weights
        Int? base_mem
        String? extra_args
        File? sites_to_extract_from_secondary_source
        String? chromosome_encoding
        Boolean? use_ref_alt_for_ids

        Int preemptible = 1
    }
  
    call ScoringTasks.ScoreVcf as PrimaryScore {
        input:
            vcf = primary_vcf,
            basename = basename,
            weights = primary_weights,
            exclude_sites = sites_to_extract_from_secondary_source,
            base_mem = base_mem,
            extra_args = extra_args,
            chromosome_encoding = chromosome_encoding,
            use_ref_alt_for_ids = use_ref_alt_for_ids
    }

    # Check for any of these parameters to be set, in order to make sure
    # that this workflow fails if only some of the parameters are set.
    if (defined(secondary_vcf) || defined(secondary_weights) || defined(sites_to_extract_from_secondary_source)) {
        call ScoringTasks.ScoreVcf as SecondaryScore {
            input:
                vcf = select_first([secondary_vcf]),
                basename = basename,
                weights = select_first([secondary_weights]),
                sites = sites_to_extract_from_secondary_source,
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
    }

    output {
        File score = select_first([CombinePrimaryAndSecondaryScores.score, PrimaryScore.score])
        File sites_scored = select_first([CombinePrimaryAndSecondaryScores.sites_scored, PrimaryScore.sites_scored])
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