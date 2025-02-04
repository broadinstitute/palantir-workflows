version 1.0

workflow Glimpse2MergeBatches {
    input {
        Array[File] imputed_vcfs
        Array[File] imputed_vcf_indices

        Array[File] qc_metrics

        String output_basename

        String docker_gatk = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        String docker_count_samples = "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        String docker_merge = "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"

        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_metrics_intervals.interval_list"
        Int scatter_count = 100

        Int mem_gb_merge = 16
    }

    if (length(imputed_vcfs) != length(qc_metrics)) {
        call ErrorWithMessage {
            input:
                message = "inputed_vcfs and qc_metrics have different lengths"
        }
    }




    if (length(imputed_vcfs) > 1) {

        call ScatterIntervalList {
            input:
                interval_list = interval_list,
                scatter_count = scatter_count
        }
        call CountVariants as CountVariantsInitial  {
            input:
                vcf = imputed_vcfs[0],
                docker_gatk = docker_gatk
        }
        scatter(batch_index in range(length(imputed_vcfs))) {
            call CountSamples {
                    input:
                        imputed_vcf = imputed_vcfs[batch_index],
                        imputed_vcf_index = imputed_vcf_indices[batch_index],
                        docker_count_samples = docker_count_samples
                }
        }
        
        scatter (intervals in ScatterIntervalList.out) {
            scatter(batch_index in range(length(imputed_vcfs))) {

                call SubsetToIntervals {
                    input:
                        vcf = imputed_vcfs[batch_index],
                        vcf_index = imputed_vcf_indices[batch_index],
                        interval_list = intervals,
                        docker_gatk = docker_gatk
                }
                call ExtractAnnotations {
                    input:
                        imputed_vcf = SubsetToIntervals.subset_vcf,
                        imputed_vcf_index = SubsetToIntervals.subset_vcf_index,
                        batch_index = batch_index,
                        docker_extract_annotations = docker_gatk
                }    
            }

            call MergeVcfs {
                input:
                    imputed_vcfs = SubsetToIntervals.subset_vcf,
                    imputed_vcf_indices = SubsetToIntervals.subset_vcf_index,
                    docker_merge = docker_merge,
                    mem_gb = mem_gb_merge
            }

            call RecomputeAndAnnotate {
                input:
                    merged_vcf = MergeVcfs.merged_vcf,
                    annotations = ExtractAnnotations.annotations,
                    num_samples = CountSamples.num_samples,
                    output_basename = "merged_annotated",
                    mem_gb = mem_gb_merge,
                    docker_merge = docker_merge
            }
        }

        call GatherVcfs {
            input:
                input_vcfs = RecomputeAndAnnotate.merged_imputed_vcf,
                output_vcf_name = output_basename + ".vcf.gz"
        }

        call CountSamples as CountSamplesFinal {
            input:
                imputed_vcf = GatherVcfs.output_vcf,
                imputed_vcf_index = GatherVcfs.output_vcf_index,
                docker_count_samples = docker_count_samples
        }
    
        call MergeQCMetrics {
            input:
                qc_metrics = qc_metrics,
                docker_merge = docker_merge,
                output_basename = output_basename
        }
        
    }

    output {
        File merged_imputed_vcf = select_first([GatherVcfs.output_vcf, imputed_vcfs[0]])
        File merged_imputed_vcf_index = select_first([GatherVcfs.output_vcf_index, imputed_vcf_indices[0]])
        Int? initial_site_count = CountVariantsInitial.count
        Int? final_sample_count = CountSamplesFinal.num_samples

        File merged_qc_metrics = select_first([MergeQCMetrics.merged_qc_metrics, qc_metrics[0]])
    }
}

task GatherVcfs {

  input {
    Array[File] input_vcfs
    String output_vcf_name
    Int disk_size_gb = ceil(1.2 * size(input_vcfs, "GiB") + 50)
    Int machine_mem_mb = 7000
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
  }

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options "-Xms6000m -Xmx6500m" \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input ~{sep=" --input " input_vcfs} \
      --output ~{output_vcf_name}

    tabix ~{output_vcf_name}
  >>>

  runtime {
    memory: "~{machine_mem_mb} MiB"
    cpu: "1"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size_gb + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}


task SubsetToIntervals {
    input {
        File vcf
        File vcf_index
        File interval_list
        Int disk_size_gb = ceil(size(vcf, "GiB") + 50)

        String docker_gatk
    }

    parameter_meta {
        vcf: {
                localization_optional: true
            }

        vcf_index: {
                localization_optional: true
            }
    }

    command <<<
        gatk SelectVariants -V ~{vcf} -L ~{interval_list} -O subset_intervals.vcf.gz
    >>>

    output {
        File subset_vcf = "subset_intervals.vcf.gz"
        File subset_vcf_index = "subset_intervals.vcf.gz.tbi"
    }

    runtime {
        docker: docker_gatk
        memory: "2000 MiB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}

task CountVariants {
    input {
        File vcf
        String docker_gatk
    }

    parameter_meta {
        vcf: {
                localization_optional: true
            }
    }

    command <<<
        gatk CountVariants -V ~{vcf} -O count.txt
    >>>

    output {
        Int count = read_int("count.txt")
    }

    runtime {
        docker: docker_gatk
    }

}

task ScatterIntervalList {
  input {
    File interval_list
    Int scatter_count
  }

  command <<<
    set -e
    mkdir out
    java -Xms1000m -Xmx1500m -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=INTERVAL_SUBDIVISION \
      UNIQUE=true \
      SORT=true \
      INPUT=~{interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1).zfill(4) + filename)
      os.rename(interval, newName)
    with open('num_intervals.txt', 'w') as num_intervals:
      num_intervals.write(f'{len(intervals)}\n')
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int('num_intervals.txt')
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-python:1.0.0-2.26.10-1663951039"
    memory: "2000 MiB"
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

task MergeVcfs {
    input {
        Array[File] imputed_vcfs
        Array[File] imputed_vcf_indices
        String docker_merge
        Int disk_size_gb = ceil(2.2 * size(imputed_vcfs, "GiB") + 50)
        Int mem_gb
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -xeuo pipefail
        bcftools merge -O z -o merged.vcf.gz ~{sep=" " imputed_vcfs}
    >>>

    runtime {
            docker: docker_merge
            disks: "local-disk " + disk_size_gb + " HDD"
            memory: mem_gb + " GiB"
            cpu: cpu
            preemptible: preemptible
        }

    output {
        File merged_vcf = "merged.vcf.gz"
    }
}

task RecomputeAndAnnotate {
    input {
        File merged_vcf
        Array[File] annotations

        Array[Int] num_samples

        String output_basename

        String docker_merge
        Int disk_size_gb = ceil(2.2 * size(merged_vcf, "GiB") + 50)
        Int mem_gb
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
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

EOF
        python3 script.py

        bgzip aggregated_annotations.tsv
        tabix -s1 -b2 -e2 aggregated_annotations.tsv.gz

        bcftools annotate -a aggregated_annotations.tsv.gz -c CHROM,POS,REF,ALT,AF,INFO -O z -o ~{output_basename}.vcf.gz ~{merged_vcf}
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
        File aggregated_annotations = "aggregated_annotations.tsv.gz"
    }
}   

task MergeQCMetrics {
    input {
        Array[File] qc_metrics
        String docker_merge
        String output_basename

        Int disk_size_gb = ceil(2.2 * size(qc_metrics, "GiB") + 50)
        Int mem_gb = 4
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -xeuo pipefail
        

        python3 <<EOF
        import pandas as pd
        qc_metrics = ['~{sep="', '" qc_metrics}']
        merged_qc_metrics = pd.concat([pd.read_csv(qc_metric, sep='\t') for qc_metric in qc_metrics])
        merged_qc_metrics.to_csv('~{output_basename}.qc_metrics.tsv', sep='\t')
        EOF
    >>>

    output {
        File merged_qc_metrics = "~{output_basename}.qc_metrics.tsv"
    }

    runtime {
        docker: docker_merge
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
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
