version 1.0

workflow Glimpse2Imputation {
    input {
        File sites_vcf
        File sites_vcf_index
        File sites_tsv
        File sites_tsv_index

        Int bcftools_threads
        Boolean use_gatk
        Int calling_batch_size
        Int calling_mem_gb

        Array[File]? crams
        Array[File]? cram_indices
        Array[String] sample_ids
        File? fasta
        File? fasta_index
        String output_basename

        File ref_dict

        Boolean call_indels = false
    }

    call SplitCRAMsIntoBatches {
        input:
            batch_size = calling_batch_size,
            crams = select_first([crams]),
            cram_indices = select_first([cram_indices]),
            sample_ids = sample_ids
    }

    scatter(i in range(length(SplitCRAMsIntoBatches.crams_batches))) {
        if (use_gatk) {
            call GATKCall {
                input:
                    crams = SplitCRAMsIntoBatches.crams_batches[i],
                    cram_indices = SplitCRAMsIntoBatches.cram_indices_batches[i],
                    fasta = select_first([fasta]),
                    fasta_index = select_first([fasta_index]),
                    fasta_dict = ref_dict,
                    sites_tsv = sites_tsv,
                    sites_tsv_index = sites_tsv_index,
                    sample_ids = SplitCRAMsIntoBatches.sample_ids_batches[i],
                    cpu = bcftools_threads,
                    mem_gb = calling_mem_gb
            }
        }
        if (!use_gatk) {
            call BcftoolsCall {
                input:
                    crams = SplitCRAMsIntoBatches.crams_batches[i],
                    cram_indices = SplitCRAMsIntoBatches.cram_indices_batches[i],
                    fasta = select_first([fasta]),
                    fasta_index = select_first([fasta_index]),
                    call_indels = call_indels,
                    sites_vcf = sites_vcf,
                    sites_vcf_index = sites_vcf_index,
                    sites_tsv = sites_tsv,
                    sites_tsv_index = sites_tsv_index,
                    sample_ids = SplitCRAMsIntoBatches.sample_ids_batches[i],
                    cpu = bcftools_threads,
                    mem_gb = calling_mem_gb
            }
        }
        File called_vcf = select_first([BcftoolsCall.output_vcf, GATKCall.output_vcf])
        File called_vcf_index = select_first([BcftoolsCall.output_vcf_index, GATKCall.output_vcf_index])
    }

    if (length(SplitCRAMsIntoBatches.crams_batches) > 1) {
        call BcftoolsMerge {
            input:
                vcfs = called_vcf,
                vcf_indices = called_vcf_index,
                output_basename = output_basename
        }
    }

    output {
        File merged_vcf = select_first([BcftoolsMerge.merged_vcf, called_vcf[0]])
        File merged_vcf_index = select_first([BcftoolsMerge.merged_vcf_index, called_vcf_index[0]])
    }
}
task SplitCRAMsIntoBatches {
    input {
        Int batch_size

        Array[String] crams
        Array[String] cram_indices
        Array[String] sample_ids
    }

    command <<<
        cat <<EOF > script.py
import json

batch_size = ~{batch_size}
crams = ['~{sep="', '" crams}']
cram_indices = ['~{sep="', '" cram_indices}']
sample_ids = ['~{sep="', '" sample_ids}']
        
crams_batches = [crams[i:i + batch_size] for i in range(0, len(crams), batch_size)]
cram_indices_batches = [cram_indices[i:i + batch_size] for i in range(0, len(cram_indices), batch_size)]
sample_ids_batches = [sample_ids[i:i + batch_size] for i in range(0, len(sample_ids), batch_size)]

with open('crams.json', 'w') as json_file:
    json.dump(crams_batches, json_file)
with open('cram_indices.json', 'w') as json_file:
    json.dump(cram_indices_batches, json_file)
with open('sample_ids.json', 'w') as json_file:
    json.dump(sample_ids_batches, json_file)
EOF
        python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "1 GiB"
        preemptible: 3
    }

    output {
        Array[Array[String]] crams_batches = read_json('crams.json')
        Array[Array[String]] cram_indices_batches = read_json('cram_indices.json')
        Array[Array[String]] sample_ids_batches = read_json('sample_ids.json')
    }
}

task SplitArraysIntoBatches {
    input {
        Int batch_size

        Array[String] arr1
        Array[String]? arr2
        Array[String]? arr3
    }

    command <<<
        cat <<EOF > script.py
import json

batch_size = ~{batch_size}

arr = ['~{sep="', '" arr1}']
        
batches = [arr[i:i + batch_size] for i in range(0, len(arr), batch_size)]

with open('batches.json', 'w') as json_file:
    json.dump(batches, json_file)
EOF
        python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "1 GiB"
        preemptible: 3
    }

    output {
        Array[Array[String]] batches = read_json('batches.json')
    }
}

task BcftoolsCall {
    input {
        Array[File] crams
        Array[File] cram_indices
        File fasta
        File fasta_index
        Boolean call_indels
        Array[String] sample_ids

        File sites_vcf
        File sites_vcf_index
        File sites_tsv
        File sites_tsv_index
        File sites_tsv
        File sites_tsv_index
        Int mem_gb = 4
        Int cpu = 2
        Int preemptible = 3
    }

    Int disk_size_gb = ceil(1.5*size(crams, "GiB") + size(fasta, "GiB") + size(sites_tsv, "GiB")) + 10

    String out_basename = "batch"

    command <<<
        set -euo pipefail

        crams=(~{sep=' ' crams})
        sample_ids=(~{sep=' ' sample_ids})

        for i in "${!crams[@]}"; do
            echo "* ${crams[$i]} ${sample_ids[$i]}" >> sample_name_mapping.txt
        done

        bcftools mpileup -f ~{fasta} ~{if !call_indels then "-I" else ""} -G sample_name_mapping.txt -E -a 'FORMAT/DP,FORMAT/AD' -T ~{sites_vcf} -O u ~{sep=" " crams} | \
        bcftools call -Aim -C alleles -T ~{sites_tsv} -O u | \
        bcftools norm -m -both -O z -o ~{out_basename}.vcf.gz
        bcftools index -t ~{out_basename}.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File output_vcf = "~{out_basename}.vcf.gz"
        File output_vcf_index = "~{out_basename}.vcf.gz.tbi"
    }
}

task GATKCall {
    input {
        Array[File] crams
        Array[File] cram_indices
        File fasta
        File fasta_index
        File fasta_dict
        Array[String] sample_ids

        File sites_tsv
        File sites_tsv_index
        Int mem_gb = 4
        Int cpu = 2
        Int preemptible = 1
    }

    Int disk_size_gb = ceil(1.5*size(crams, "GiB") + size(fasta, "GiB") + size(sites_tsv, "GiB")) + 10

    String out_basename = "batch"

    command <<<
        set -xeuo pipefail

        gatk --java-options "-Xmx~{mem_gb}g" HaplotypeCaller \
            -R ~{fasta} \
            -I ~{sep=" -I " crams} \
            -O ~{out_basename}.vcf.gz \
            --alleles ~{sites_tsv} \
            -L ~{sites_tsv}
        bcftools index -t ~{out_basename}.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File output_vcf = "~{out_basename}.vcf.gz"
        File output_vcf_index = "~{out_basename}.vcf.gz.tbi"
    }
}

task BcftoolsMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_indices
        Int mem_gb = 4
        Int cpu = 2
        Int preemptible = 1

        String output_basename
    }

    Int disk_size_gb = ceil(3*size(vcfs, "GiB")) + 10

    command <<<
        set -euo pipefail
        bcftools merge -O z -o ~{output_basename}.bcftools.merged.vcf.gz ~{sep=" " vcfs}
        bcftools index -t ~{output_basename}.bcftools.merged.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File merged_vcf = "~{output_basename}.bcftools.merged.vcf.gz"
        File merged_vcf_index = "~{output_basename}.bcftools.merged.vcf.gz.tbi"
    }
}