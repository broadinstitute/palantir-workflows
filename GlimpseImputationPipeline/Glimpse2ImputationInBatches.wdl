version 1.0

import "Glimpse2Imputation.wdl" as Imputation
import "Glimpse2MergeBatches.wdl" as MergeBatches

workflow Glimpse2ImputationInBatches {
    input {
        Int batch_size

        # List of files, one per line
        File reference_chunks

        Array[File] crams
        Array[File] cram_indices
        Array[String] sample_ids
        File fasta
        File fasta_index
        String output_basename

        File ref_dict

        Boolean? collect_qc_metrics
        Boolean? impute_reference_only_variants
        Boolean? call_indels
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
        
        Int? preemptible
        String? docker
        String? docker_extract_num_sites_from_reference_chunk
        Int? cpu_ligate
        Int? mem_gb_ligate
        Int? mem_gb_merge
        File? monitoring_script

        String? docker_gatk
        String? docker_count_samples
        String? docker_merge

        File interval_list_for_merge
        Int? merge_scatter_count
    }

    call SplitIntoBatches {
        input:
            batch_size = batch_size,
            crams = crams,
            cram_indices = cram_indices,
            sample_ids = sample_ids
    }

    scatter(i in range(length(SplitIntoBatches.crams_batches))) {
        call Imputation.Glimpse2Imputation {
            input:
                reference_chunks = reference_chunks,
                crams = SplitIntoBatches.crams_batches[i],
                cram_indices = SplitIntoBatches.cram_indices_batches[i],
                sample_ids = SplitIntoBatches.sample_ids_batches[i],
                fasta = fasta,
                fasta_index = fasta_index,
                output_basename = output_basename + "_batch_" + i,
                collect_qc_metrics = collect_qc_metrics,
                ref_dict = ref_dict,
                impute_reference_only_variants = impute_reference_only_variants,
                call_indels = call_indels,
                n_burnin = n_burnin,
                n_main = n_main,
                effective_population_size = effective_population_size,
                preemptible = preemptible,
                docker = docker,
                docker_extract_num_sites_from_reference_chunk = docker_extract_num_sites_from_reference_chunk,
                cpu_ligate = cpu_ligate,
                mem_gb_ligate = mem_gb_ligate,
                monitoring_script = monitoring_script
        }
    }

    # We can't access collect_qc_metrics here because we want to inherit the default value
    # from of the Glimpse2Imputation workflow. But we know that if collect_qc_metrics has
    # been set there then the qc_metrics output will be defined, so just check if the output
    # exists for the first batch.
    Boolean merge_qc_metrics = defined(Glimpse2Imputation.qc_metrics[0])

    if (merge_qc_metrics) {
        Array[File] qc_metrics = select_all(Glimpse2Imputation.qc_metrics)
    }

    call MergeBatches.Glimpse2MergeBatches {
        input:
            imputed_vcfs = Glimpse2Imputation.imputed_vcf,
            imputed_vcf_indices = Glimpse2Imputation.imputed_vcf_index,
            output_basename = output_basename,
            qc_metrics = qc_metrics,
            mem_gb_merge = mem_gb_merge,
            docker_gatk = docker_gatk,
            docker_count_samples = docker_count_samples,
            docker_merge = docker_merge,
            interval_list = interval_list_for_merge,
            scatter_count = merge_scatter_count
    }

    output {
        File merged_imputed_vcf = Glimpse2MergeBatches.merged_imputed_vcf
        File merged_imputed_vcf_index = Glimpse2MergeBatches.merged_imputed_vcf_index
        
        # No select_first here because we want this to be null if merge_qc_metrics is false
        File? merged_qc_metrics = Glimpse2MergeBatches.merged_qc_metrics
    }
}

task SplitIntoBatches {
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
