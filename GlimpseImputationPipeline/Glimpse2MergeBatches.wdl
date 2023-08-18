version 1.0

workflow Glimpse2MergeBatches {
    input {
        Array[File] imputed_vcfs
        Array[File] imputed_vcf_indices

        String output_basename

        String docker_extract_annotations = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        String docker_count_samples = "us.gcr.io/broad-dsde-methods/bcftools:v1.2"
        String docker_merge = "us.gcr.io/broad-dsde-methods/samtools-suite:v1"
    }

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
            annotations = ExtractAnnotations.annotations,
            num_samples = CountSamples.num_samples,
            output_basename = output_basename,
            docker_merge = docker_merge
    }

    output {
        File merged_imputed_vcf = MergeAndRecomputeAndAnnotate.merged_imputed_vcf
        File merged_imputed_vcf_index = MergeAndRecomputeAndAnnotate.merged_imputed_vcf_index
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

    command <<<
        set -xeuo pipefail
        export GCS_OAUTH_TOKEN=$(/root/google-cloud-sdk/bin/gcloud auth application-default print-access-token)

        bcftools merge -O z -o ~{output_basename}.merged.vcf.gz ~{sep=" " imputed_vcfs}

        cat <<EOF > script.py
import gzip
from contextlib import ExitStack

input_filenames = ['~{sep="', '" annotations}']
num_samples = [~{sep=", " num_samples}]

def check_identical(fields):
    return all(field == fields[0] for field in fields)

if len(num_samples) != len(input_filenames):
    raise RuntimeError('The number of input annotations does not match the number of input number of samples.')
num_batches = len(input_filenames)

with open('aggregated_annotations.tsv', 'w') as output_file:
    with ExitStack() as exit_stack:
        input_files = [exit_stack.enter_context(gzip.open(filename, 'rt')) for filename in input_filenames]

        lines_iterator = iter(zip(*input_files))
        next(lines_iterator)
        for lines in lines_iterator:
            lines_data = [line.strip().split('\t') for line in lines]
            contig = [line_data[0] for line_data in lines_data]
            pos = [line_data[1] for line_data in lines_data]
            ref = [line_data[2] for line_data in lines_data]
            alt = [line_data[3] for line_data in lines_data]
            af = [float(line_data[4]) for line_data in lines_data]
            info = [float(line_data[5]) for line_data in lines_data]

            if not check_identical(contig) or not check_identical(pos) or not check_identical(ref) or not check_identical(alt):
                raise RuntimeError('The site information (CHROM, POS, REF, ALT) is not identical in each line between different batch annotations.')
            
            aggregated_af = sum([af[i] * num_samples[i] for i in range(num_batches)]) / sum(num_samples)
            aggregated_info = 1 if aggregated_af == 0 or aggregated_af == 1 else \
                     1 - \
                    (sum([1 - info[i] * 2 * num_samples[i] * af[i] * (1 - af[i]) for i in range(num_batches)])) / \
                    (2 * sum(num_samples) * aggregated_af * (1 - aggregated_af))
            
            output_file.write(f'{contig}\t{pos}\t{ref}\t{alt}\t{aggregated_af}\t{aggregated_info}\n')
EOF
        python3 script.py

        bgzip aggregated_annotations.tsv
        tabix -s1 -b2 -e2 aggregated_annotations.tsv.gz

        bcftools annotate -a aggregated_annotations.tsv.gz -c CHROM,POS,REF,ALT,AF,INFO -O z -o ~{output_basename}.vcf.gz ~{output_basename}.merged.vcf.gz
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
    }
}