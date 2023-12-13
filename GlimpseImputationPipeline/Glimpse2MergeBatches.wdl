version 1.0

workflow Glimpse2MergeBatches {
    input {
        Array[File] imputed_vcfs
        Array[File] imputed_vcf_indices

        Array[File?]? qc_metrics

        String output_basename

        String docker_extract_annotations = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        String docker_count_samples = "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        String docker_merge = "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
    }

    if (defined(qc_metrics) && (length(imputed_vcfs) != length(select_first([qc_metrics, []])))) {
        call ErrorWithMessage {
            input:
                message = "inputed_vcfs and qc_metrics have different lengths"
        }
    }

    if (length(imputed_vcfs) > 1) {
        scatter(batch_index in range(length(imputed_vcfs))) {
            call ExtractAnnotations {
                input:
                    imputed_vcf = imputed_vcfs[batch_index],
                    imputed_vcf_index = imputed_vcf_indices[batch_index],
                    batch_index = batch_index,
                    docker_extract_annotations = docker_extract_annotations
            }
            call CountSamples {
                input:
                    imputed_vcf = imputed_vcfs[batch_index],
                    imputed_vcf_index = imputed_vcf_indices[batch_index],
                    docker_count_samples = docker_count_samples
            }
        }

        call MergeAndRecomputeAndAnnotate {
            input:
                imputed_vcfs = imputed_vcfs,
                imputed_vcf_indices = imputed_vcf_indices,
                qc_metrics = qc_metrics,
                annotations = ExtractAnnotations.annotations,
                num_samples = CountSamples.num_samples,
                output_basename = output_basename,
                docker_merge = docker_merge
        }
    }

    # This handles the qc_metrics in the case of only one batch. If qc_metrics is not
    # defined then this if statement will be false and qc_metrics_1 will be null.
    # If qc_metrics is defined then set qc_metrics_1 to the first element, which
    # may or may not be null.
    if (length(imputed_vcfs) == 1 && defined(qc_metrics)) {
        File? qc_metrics_1 = select_first([qc_metrics, []])[0]
    }

    output {
        File merged_imputed_vcf = select_first([MergeAndRecomputeAndAnnotate.merged_imputed_vcf, imputed_vcfs[0]])
        File merged_imputed_vcf_index = select_first([MergeAndRecomputeAndAnnotate.merged_imputed_vcf_index, imputed_vcf_indices[0]])

        # If input qc_metrics are defined then we want to return the merged qc_metrics here.
        # MergeAndRecomputeAndAnnotate handles all possible null cases if there are multiple
        # batches. If there is only one batch (and MergeAndRecomputeAndAnnotate is therefore
        # not called), return qc_metrics_1, which implements that logic above.
        File? merged_qc_metrics = if length(imputed_vcfs) > 1 then MergeAndRecomputeAndAnnotate.merged_qc_metrics else qc_metrics_1
    }
}

task ExtractAnnotations {
    input {
        File imputed_vcf
        File imputed_vcf_index
        Int batch_index
        
        String docker_extract_annotations
        Int disk_size_gb = ceil(2 * size(imputed_vcf, "GiB") + 50)
        Int mem_gb = 2
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        gatk VariantsToTable -V ~{imputed_vcf} -O annotations_batch_~{batch_index}.tsv -F CHROM -F POS -F REF -F ALT -F AF -F INFO
        bgzip annotations_batch_~{batch_index}.tsv
    >>>

    runtime {
        docker: docker_extract_annotations
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File annotations = "annotations_batch_~{batch_index}.tsv.gz"
    }
}

task CountSamples { # really?
    input {
        File imputed_vcf
        File imputed_vcf_index
        
        String docker_count_samples
        Int disk_size_gb = ceil(1 * size(imputed_vcf, "GiB") + 10)
        Int mem_gb = 2
        Int cpu = 2
        Int preemptible = 1
    }


    parameter_meta {
        imputed_vcf:
            {
                localization_optional: true
            }
        imputed_vcf_index:
            {
                localization_optional: true
            }
    }

    command <<<
        set -xeuo pipefail
        export GCS_OAUTH_TOKEN=$(/root/google-cloud-sdk/bin/gcloud auth application-default print-access-token)
        bcftools query -l ~{imputed_vcf} | wc -l > "num_samples.txt"
    >>>

    runtime {
        docker: docker_count_samples
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        Int num_samples = read_int("num_samples.txt")
    }
}

task MergeAndRecomputeAndAnnotate {
    input {
        Array[File] imputed_vcfs
        Array[File] imputed_vcf_indices
        Array[File] annotations
        Array[Int] num_samples

        Array[File?]? qc_metrics

        String output_basename

        String docker_merge
        Int disk_size_gb = ceil(2.2 * size(imputed_vcfs, "GiB") + 50)
        Int mem_gb = 2
        Int cpu = 2
        Int preemptible = 1
    }


    parameter_meta {
        imputed_vcfs:
            {
                localization_optional: true
            }
        imputed_vcf_indices:
            {
                localization_optional: true
            }
    }

    # We need this variable because it's tricky to expand a Array[File?]? in the command block.
    # This statement converts the Array[File?]? to an Array[File] according to the following rule:
    # If qc_metrics is either not defined or only contains null values (the conditional statement
    # checks for both) then this returns an empty array, otherwise it returns qc_metrics.
    Array[File] qc_metrics_expanded = if length(select_first([qc_metrics, []])) > 0 then select_all(select_first([qc_metrics, []])) else []

    command <<<
        set -xeuo pipefail
        export GCS_OAUTH_TOKEN=$(/root/google-cloud-sdk/bin/gcloud auth application-default print-access-token)

        bcftools merge -O z -o ~{output_basename}.merged.vcf.gz ~{sep=" " imputed_vcfs}

        cat <<EOF > script.py
import pandas as pd
import functools

input_filenames = ['~{sep="', '" annotations}']
num_samples = [~{sep=", " num_samples}]
if len(num_samples) != len(input_filenames):
    raise RuntimeError('The number of input annotations does not match the number of input number of samples.')
num_batches = len(input_filenames)

def calculate_af(row):
    return sum([row[f'AF_{i}'] * num_samples[i] for i in range(num_batches)]) / sum(num_samples)
def calculate_info(row):
    aggregated_af = row['AF']
    return 1 if aggregated_af == 0 or aggregated_af == 1 else \
                     1 - \
                    (sum([(1 - row[f'INFO_{i}']) * 2 * num_samples[i] * row[f'AF_{i}'] * (1 - row[f'AF_{i}']) for i in range(num_batches)])) / \
                    (2 * sum(num_samples) * aggregated_af * (1 - aggregated_af))

annotation_dfs = [pd.read_csv(input_filename, sep='\t').rename(columns={'AF': f'AF_{i}', 'INFO': f'INFO_{i}'}) for i, input_filename in enumerate(input_filenames)]
annotations_merged = functools.reduce(lambda left, right: pd.merge(left, right, on=['CHROM', 'POS', 'REF', 'ALT'], how='inner', validate='one_to_one'), annotation_dfs)
annotations_merged['AF'] = annotations_merged.apply(lambda row: calculate_af(row), axis=1)
annotations_merged['INFO'] = annotations_merged.apply(lambda row: calculate_info(row), axis=1)
annotations_merged.to_csv('aggregated_annotations.tsv', sep='\t', columns=['CHROM', 'POS', 'REF', 'ALT', 'AF', 'INFO'], header=False, index=False)

if ~{if length(qc_metrics_expanded) > 0 then "True" else "False"}:
    qc_metrics = ['~{sep="', '" qc_metrics_expanded}']
    merged_qc_metrics = pd.concat([pd.read_csv(qc_metric, sep='\t') for qc_metric in qc_metrics])
    merged_qc_metrics.to_csv('~{output_basename}.qc_metrics.tsv', sep='\t')
EOF
        python3 script.py

        bgzip aggregated_annotations.tsv
        tabix -s1 -b2 -e2 aggregated_annotations.tsv.gz

        bcftools annotate -a aggregated_annotations.tsv.gz -c CHROM,POS,REF,ALT,AF,INFO -O z -o ~{output_basename}.vcf.gz ~{output_basename}.merged.vcf.gz
        tabix ~{output_basename}.vcf.gz
    >>>

    runtime {
        docker: docker_merge
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File merged_imputed_vcf = "~{output_basename}.vcf.gz"
        File merged_imputed_vcf_index = "~{output_basename}.vcf.gz.tbi"
        File aggregated_annotations = "aggregated_annotations.tsv.gz"

        File? merged_qc_metrics = "~{output_basename}.qc_metrics.tsv"
    }
}

task ErrorWithMessage{
    input {
        String message
    }
    command <<<
    >&2 echo "Error: ~{message}"
    exit 1
    >>>

    runtime {
        docker: "ubuntu"
    }
}
