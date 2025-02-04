version 1.0

workflow MergeSingleSampleMinimacVcfs {
    input {
        Array[File] vcfs
        Int n_per_batch
        String output_basename

        File interval_list
        Int interval_scatter_count = 100
    }

    Int n_batches = ceil(length(vcfs)/n_per_batch)

    scatter (i_batch in range(n_batches)) {
        scatter (i in range(length(vcfs))) {
            Int start_index = i_batch * n_per_batch
            Int end_index = start_index + n_per_batch - 1
            Boolean within_batch = i >= start_index && i <= end_index && i < length(vcfs)
            if (within_batch) {
                File sliced_elements = vcfs[i]
            } 
        }

        Array[File] vcfs_this_batch = select_all(sliced_elements)

        call dfvm_merge {
            input:
                vcfs = vcfs_this_batch
        }

        call CollectQCMetrics {
            input:
                imputed_vcf = dfvm_merge.merged_vcf,
                output_basename = output_basename
        }
    }

    call ScatterIntervalList {
        input:
            interval_list = interval_list,
            scatter_count = interval_scatter_count
    }

    scatter (intervals in ScatterIntervalList.out) {
            scatter(batch_index in range(n_batches)) {

                call SubsetToIntervals {
                    input:
                        vcf = dfvm_merge.merged_vcf[batch_index],
                        vcf_index = dfvm_merge.merged_vcf_index[batch_index],
                        interval_list = intervals
                }

            }

            call dfvm_merge as dfvm_merge_full {
                input:
                    vcfs = SubsetToIntervals.subset_vcf
            }

            call CountSamples {
                input:
                    imputed_vcf = dfvm_merge_full.merged_vcf
            }

            call reannotate_from_dosages {
                input:
                    vcf = dfvm_merge_full.merged_vcf,
                    n_samples = CountSamples.num_samples
            }
        }

    call GatherVcfs {
        input:
            input_vcfs = reannotate_from_dosages.annotated_vcf,
            output_vcf_name = output_basename + ".vcf.gz"
    }

    call MergeQCMetrics {
            input:
                qc_metrics = CollectQCMetrics.qc_metrics,
                output_basename = output_basename
        }

    output {
        File merged_vcf = GatherVcfs.output_vcf
        File merged_imputed_vcf_index = GatherVcfs.output_vcf_index
        File merged_qc_metrics = MergeQCMetrics.merged_qc_metrics
    }
    
}

task CollectQCMetrics {
    input {
        File imputed_vcf
        String output_basename
        
        Int preemptible = 1
        String docker = "hailgenetics/hail:0.2.126-py3.11"
        Int cpu = 4
        Int mem_gb = 16
    }

    parameter_meta {
        imputed_vcf: {
                        localization_optional: true
                    }
    }

    Int disk_size_gb = 100
    
    command <<<
        set -euo pipefail

        cat <<'EOF' > script.py
import hail as hl
import pandas as pd

# Calculate metrics
hl.init(default_reference='GRCh37', idempotent=True)
vcf = hl.import_vcf('~{imputed_vcf}', force_bgz=True)
qc = hl.sample_qc(vcf)
qc_pd = qc.cols().flatten() \
    .rename({'sample_qc.' + col: col for col in list(qc['sample_qc'])}) \
    .rename({'s': 'sample_id'}) \
    .to_pandas()
qc_pd.to_csv('~{output_basename}.qc_metrics.tsv', sep='\t', index=False, float_format='%.4f')
EOF
        python3 script.py
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File qc_metrics = "~{output_basename}.qc_metrics.tsv"
    }
}

task GatherVcfs {

  input {
    Array[File] input_vcfs
    String output_vcf_name
    Int disk_size_gb = ceil(1.2 * size(input_vcfs, "GiB") + 50)
    Int machine_mem_mb = 7000
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int preemptible = 1
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
    preemptible: preemptible
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task MergeQCMetrics {
    input {
        Array[File] qc_metrics
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
        docker: "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

task SubsetToIntervals {
    input {
        File vcf
        File vcf_index
        File interval_list
        Int disk_size_gb = ceil(size(vcf, "GiB") + 50)
        Int preemptible = 3

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
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        memory: "2000 MiB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
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

task dfvm_merge {
    input {
            Array[File] vcfs

             Int disk_size_gb = ceil(2.2 * size(vcfs, "GiB") + 50)
            Int mem_gb = 16
            Int cpu = 2
            Int preemptible = 3
            }
        
        command <<<
            /docker_build/dfvm/dfvm -O merged.vcf.gz -I ~{sep=" -I " vcfs}
            tabix merged.vcf.gz
        >>>

        runtime {
            docker: "us.gcr.io/broad-dsde-methods/dfvm@sha256:387d308771c3073679eea5ecbffc9e629d112ddea5c70869afa456737bcb758e"
            disks: "local-disk " + disk_size_gb + " SSD"
            memory: mem_gb + " GiB"
            cpu: cpu
            preemptible: preemptible
        }

        output {
            File merged_vcf = "merged.vcf.gz"
            File merged_vcf_index = "merged.vcf.gz.tbi"
        }
}

task reannotate_from_dosages {
    input {
        File vcf
        Int n_samples
        Int disk_size_gb = ceil(3.3 * size(vcf, "GiB")) + 50
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
    }

    String output_base = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail

        bcftools query -f '%CHROM,%POS,%REF,%ALT[,%DS]' ~{vcf} | gzip -c > dosage_tbl.csv.gz

        echo "dosages extracted"

        python3 << EOF
        import pandas as pd
        import numpy as np


        # chunksize calculation
        # 2x safety factor
        # bytes per row:
        # 8x(2+n_samples) (index, pos) + 50x3 (chrom, ref, alt)
        bytes_per_row = 8*(2 + ~{n_samples}) + 150
        chunksize = int(~{mem_gb}*1e9/bytes_per_row/2)
        out_annotation_dfs = []
        for chunk in pd.read_csv("dosage_tbl.csv.gz", names=["CHROM","POS","REF","ALT"]+[f'sample_{i}' for i in range(~{n_samples})],dtype={f'sample_{i}':'float' for i in range(306)}, 
                        na_values=".", chunksize = chunksize):
            chunk = chunk.dropna()
            dosages = chunk[[f'sample_{i}' for i in range(~{n_samples})]].to_numpy()
            chunk_annotations = chunk[["CHROM","POS","REF","ALT"]]
            chunk_annotations["AF"] = dosages.mean(axis=1)/2
            chunk_annotations["R2"] = np.where((chunk_annotations["AF"]==0) | (chunk_annotations["AF"]==1),1,dosages.var(axis=1)/(2*chunk_annotations["AF"]*(1-chunk_annotations["AF"])))
            out_annotation_dfs.append(chunk_annotations)
        annotations_df = pd.concat(out_annotation_dfs)
        annotations_df.to_csv("annotations.tsv", sep="\t", index=False, header=False)
        EOF
        
        echo "annotations recomputed"

        bgzip annotations.tsv
        tabix -s1 -b2 -e2 annotations.tsv.gz

        echo "annotating vcf with new annotations"

        bcftools annotate -a annotations.tsv.gz -c CHROM,POS,REF,ALT,AF,R2 -x INFO/MAF -Oz -o ~{output_base}.reannotated.vcf.gz ~{vcf}
        
        

    >>>


    runtime {
        docker: "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File annotated_vcf = "~{output_base}.reannotated.vcf.gz"
    }
    
}

task CountSamples { # really?
    input {
        File imputed_vcf
    }


    parameter_meta {
        imputed_vcf:
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
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        preemptible: 3
    }

    output {
        Int num_samples = read_int("num_samples.txt")
    }
}