version 1.0

workflow CollectBGEImputationMetrics {
    input {
        Array[String] ancestries
        Boolean collect_af_0_1_single_number = false

        Array[String] sample_ids
		File eval_vcf
        Array[String] truth_sample_ids
		File truth_vcf
        String configuration_label
		File annotation_vcf
        Map[String, String] ancestry_to_af_annotation_map
        String? intervals

        Int preemptible = 1
	}

    call PearsonCorrelationByAF as PearsonByAF {
        input:
            evalVcf = eval_vcf,
            af_resource = annotation_vcf,
            ancestries = ancestries,
            ancestry_to_af_annotation_map = ancestry_to_af_annotation_map,
            truthVcf = truth_vcf,
            intervals = intervals,
            output_basename = "cohort",
            configuration_label = configuration_label,
            eval_sample_ids = sample_ids,
            truth_sample_ids = truth_sample_ids,
            preemptible = preemptible
    }

    if (collect_af_0_1_single_number) {
        call PearsonCorrelationByAF as PearsonByAF01 {
            input:
                evalVcf = eval_vcf,
                af_resource = annotation_vcf,
                ancestries = ancestries,
                ancestry_to_af_annotation_map = ancestry_to_af_annotation_map,
                truthVcf = truth_vcf,
                min_af_for_accuracy_metrics = 0.1,
                n_bins = 2,
                right_edge_first_bin = 0.1,
                intervals = intervals,
                configuration_label = configuration_label,
                output_basename = "cohort",
                eval_sample_ids = sample_ids,
                truth_sample_ids = truth_sample_ids,
                preemptible = preemptible
        }
    }

    output {
        File correlations_with_info = PearsonByAF.correlations_with_info
    }
}

task PearsonCorrelationByAF {
    input {
        File evalVcf
        File truthVcf
        Array[String] eval_sample_ids
        Array[String] truth_sample_ids
        Map[String, String] ancestry_to_af_annotation_map
        String output_basename
        Array[String] ancestries
        String configuration_label
        File af_resource
        File? sites
        String? intervals
        String? dosage_field
        Int? n_bins
        Float? right_edge_first_bin
        Float? min_af_for_accuracy_metrics
        Int mem_gb = 16
        Int preemptible = 1
    }

    parameter_meta {
        evalVcf : {
            localization_optional : true
        }
        truthVcf : {
            localization_optional : true
        }
        af_resource : {
            localization_optional : true
        }
    }

    command <<<
        set -xeuo pipefail

        python <<'EOF'
eval_sample_ids = ['~{sep="', '" eval_sample_ids}']
truth_sample_ids = ['~{sep="', '" truth_sample_ids}']
ancestries = ['~{sep="', '" ancestries}']

# Write sample-map (eval to truth mapping)
with open('sample_map.list', 'w') as sample_map_file:
    for i in range(len(eval_sample_ids)):
        sample_map_file.write(f'{eval_sample_ids[i]}:{truth_sample_ids[i]}\n')

# Write af_expressions (truth to ancestry-specific af annotation mapping)
ancestry_to_af_annotation_dict = dict()
with open('~{write_map(ancestry_to_af_annotation_map)}') as map_file:
    for line in map_file:
        ancestry, af_annotation = line.strip().split('\t')
        ancestry_to_af_annotation_dict[ancestry] = af_annotation

with open('af_expressions.list', 'w') as af_expressions_file:
    for i in range(len(truth_sample_ids)):
        af_expressions_file.write(f'{truth_sample_ids[i]}:{ancestry_to_af_annotation_dict[ancestries[i]]}\n')
EOF

        gatk --java-options "-Xmx~{mem_gb - 2}G" EvaluateGenotypingPerformance -eval ~{evalVcf} -truth ~{truthVcf} --af-annotations af_expressions.list --resource ~{af_resource} \
        ~{"--ids " + sites} ~{"-L " + intervals} --sample-map sample_map.list ~{"--dosage-field " + dosage_field} -O ~{output_basename}.correlations.tsv \
        -OA ~{output_basename}.accuracy.tsv ~{"-nbins " + n_bins} ~{"-first-bin-right-edge " + right_edge_first_bin} ~{"--min-af-for-accuracy-metrics " + min_af_for_accuracy_metrics} --allow-differing-ploidies

        python <<'EOF'
ancestries = ['~{sep="', '" ancestries}']
eval_sample_ids = ['~{sep="', '" eval_sample_ids}']
ancestry_dict = dict()
for i in range(len(eval_sample_ids)):
    ancestry_dict[eval_sample_ids[i]] = ancestries[i]

df = pd.read_csv('~{output_basename}.correlations.tsv', sep="\t", comment='#')
if df['SNP_SITES'].sum() + df['INDEL_SITES'].sum() == 0:
    raise RuntimeError(f'Correlation file has no sites in it.')
df['ancestry'] = df['SAMPLE'].map(ancestry_dict)
df['configuration'] = '~{configuration_label}'
df.to_csv('~{output_basename}.correlations_with_info.tsv', sep='\t')
EOF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/ckachulis/gatk-array-correlation@sha256:5910defbc137d43145e6cf1f9c3539ae418b6e465250d55465e19e77773445c4"
        disks: "local-disk 100 HDD"
        memory: mem_gb + " GB"
        preemptible: preemptible
    }

    output {
        File correlations = "~{output_basename}.correlations.tsv"
        File correlations_with_info = "~{output_basename}.correlations_with_info.tsv"
        File accuracy = "~{output_basename}.accuracy.tsv"
    }
}