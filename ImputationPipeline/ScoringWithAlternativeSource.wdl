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

    call UpdateExomeGVCFAltAlleles {
        input:
            exome_gvcf_filtered_vcf = ExtractSitesFromGvcf.vcf,
            exome_gvcf_filtered_vcf_index = ExtractSitesFromGvcf.vcf_index,
            weights = weights,
            ref_dict = ref_dict,
            basename = basename
    }
  
    call ScoringTasks.ScoreVcf as ExomeGVCFScore {
        input:
            vcf = UpdateExomeGVCFAltAlleles.vcf_extracted_from_gvcf,
            basename = basename,
            weights = weights,
            sites = UpdateExomeGVCFAltAlleles.sites_extracted_from_gvcf,
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
            exclude_sites = UpdateExomeGVCFAltAlleles.sites_extracted_from_gvcf,
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
            -O ~{basename}.filtered.vcf.gz \
            -select "(vc.hasAttribute(\"END\") && vc.getGenotype(\"~{basename}\").getGQ() > 30) || QUAL > 30"
        
        gatk GenotypeGVCFs \
            -R ~{ref_fasta} \
            -V ~{basename}.filtered.vcf.gz \
            -O ~{basename}.filtered.genotyped.vcf.gz \
            --force-output-intervals ~{sites_to_extract}

    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File vcf = "~{basename}.filtered.genotyped.vcf.gz"
        File vcf_index = "~{basename}.filtered.genotyped.vcf.gz"
    }
}

task UpdateExomeGVCFAltAlleles {
    input {
        File exome_gvcf_filtered_vcf
        File exome_gvcf_filtered_vcf_index
        File weights
        String basename
        File ref_dict

        Int mem_gb = 4
        Int cpu = 4
        Int? disk_gb
        Int preemptible = 1
    }

    Int disk_size_gb = select_first([disk_gb, ceil(size(exome_gvcf_filtered_vcf, "GiB") * 2 + size(weights, "GiB") + 50)])

    command <<<
        set -xeuo pipefail

        cat <<'EOF' > script.py
import pysam
import pandas as pd
import numpy as np

def read_ref_dict(ref_dict_path):
        with open(ref_dict_path) as ref_dict:
            return [line.split('\t')[1].replace('SN:', '') for line in ref_dict]

def read_prs_weights(ref_dict, prs_weights_path):
    prs_weights = pd.read_table(prs_weights_path)

    prs_weights[['contig', 'position', 'ref', 'alt']] = prs_weights['CHR:BP:REF:ALT'].str.split(':', expand=True)
    
    prs_weights = prs_weights.astype({'effect_allele': str, 'weight': np.float64, 'contig': str, 'position': np.int32, 'ref': str, 'alt': str})

    prs_weights = prs_weights.sort_values(by='position').sort_values(kind='stable', by='contig', key=lambda contig: contig.map(ref_dict.index))
    return prs_weights.set_index(['contig', 'position'])

ref_dict = read_ref_dict('~{ref_dict}')
prs_weights = read_prs_weights(ref_dict, '~{weights}')

gvcf_filtered_genotyped = pysam.VariantFile('~{exome_gvcf_filtered_vcf}', "r")

sites_not_in_gvcf = [(weight_index[0], weight_index[1], weight.ref, weight.alt) for weight_index, weight in prs_weights.iterrows()]
sites_in_gvcf = []

with pysam.VariantFile('~{basename}.extracted_from_gvcf.vcf.gz', 'w', header=gvcf_filtered_genotyped.header) as out:
    for record in gvcf_filtered_genotyped:
        if (record.contig, record.pos) not in prs_weights.index:
            print('Skipping:', record.contig, record.pos)
            continue
        if record.alts is None:
            if len(prs_weights.loc[(record.contig, record.pos)].shape) > 1:
                raise RuntimeError('Multiple indistinguishable weights for the same locus without alt allele:\n' + str(prs_weights.loc[(record.contig, record.pos)]))
            record.alts = [prs_weights.loc[(record.contig, record.pos)].alt]
        elif len(record.alts) != 1:
            raise RuntimeError('Only biallelic variants are supported: ' + str(record))
        
        out.write(record)
        sites_in_gvcf.append((record.contig, record.pos, record.ref, record.alts[0]))
        sites_not_in_gvcf.remove((record.contig, record.pos, record.ref, record.alts[0]))
with open('~{basename}.extracted_from_gvcf.sites', 'w') as out_sites_in_gvcf:
    out_sites_in_gvcf.write('\n'.join([f'{contig}:{position}:{ref}:{alt}' for contig, position, ref, alt in sites_in_gvcf]))
with open('~{basename}.not_extracted_from_gvcf.sites', 'w') as out_sites_not_in_gvcf:
    out_sites_not_in_gvcf.write('\n'.join([f'{contig}:{position}:{ref}:{alt}' for contig, position, ref, alt in sites_not_in_gvcf]))
EOF

        python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim-pysam:v1.0"
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File vcf_extracted_from_gvcf = "~{basename}.extracted_from_gvcf.vcf.gz"
        File sites_extracted_from_gvcf = "~{basename}.extracted_from_gvcf.sites"
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