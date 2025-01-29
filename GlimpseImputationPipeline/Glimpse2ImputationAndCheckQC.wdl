version 1.0

import "Glimpse2Imputation.wdl" as Glimpse2Imputation
import "Glimpse2CheckQC.wdl" as Glimpse2CheckQC

workflow Glimpse2ImputationAndCheckQC {
    input {
        # For documentation, please refer to Glimpse2Imputation.wdl
        File reference_chunks

        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        Array[String] sample_ids
        File? fasta
        File? fasta_index
        String output_basename

        File ref_dict

        Boolean? impute_reference_only_variants
        Boolean? call_indels
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
        
        Int? preemptible
        String? docker
        String? docker_extract_num_sites_from_reference_chunk
        Float? mem_scaling_factor_phase
        Int? cpu_ligate
        Int? mem_gb_ligate
        File? monitoring_script

        # For documentation, please refer to CheckQC.wdl
        File qc_metrics_thresholds
        Int? check_qc_preemptible
        String? check_qc_docker
        Int? check_qc_cpu
        Int? check_qc_mem_gb
    }

    call Glimpse2Imputation.Glimpse2Imputation {
        input:
            reference_chunks = reference_chunks,
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            crams = crams,
            cram_indices = cram_indices,
            sample_ids = sample_ids,
            fasta = fasta,
            fasta_index = fasta_index,
            output_basename = output_basename,
            ref_dict = ref_dict,
            impute_reference_only_variants = impute_reference_only_variants,
            call_indels = call_indels,
            n_burnin = n_burnin,
            n_main = n_main,
            effective_population_size = effective_population_size,
            preemptible = preemptible,
            docker = docker,
            docker_extract_num_sites_from_reference_chunk = docker_extract_num_sites_from_reference_chunk,
            mem_scaling_factor_phase = mem_scaling_factor_phase,
            cpu_ligate = cpu_ligate,
            mem_gb_ligate = mem_gb_ligate,
            monitoring_script = monitoring_script
    }

    call Glimpse2CheckQC.Glimpse2CheckQC {
        input:
            qc_metrics = Glimpse2Imputation.qc_metrics,
            qc_metrics_thresholds = qc_metrics_thresholds,
            output_basename = output_basename,
            preemptible = check_qc_preemptible,
            docker = check_qc_docker,
            cpu = check_qc_cpu,
            mem_gb = check_qc_mem_gb
    }

    output {
        File imputed_vcf = Glimpse2Imputation.imputed_vcf
        File imputed_vcf_index = Glimpse2Imputation.imputed_vcf_index
        File imputed_vcf_md5sum = Glimpse2Imputation.imputed_vcf_md5sum
        
        File qc_metrics = Glimpse2Imputation.qc_metrics
        File coverage_metrics = Glimpse2Imputation.coverage_metrics

        Array[File?] glimpse_phase_monitoring = Glimpse2Imputation.glimpse_phase_monitoring
        File? glimpse_ligate_monitoring = Glimpse2Imputation.glimpse_ligate_monitoring

        Boolean qc_passed = Glimpse2CheckQC.qc_passed
        File qc_failures = Glimpse2CheckQC.qc_failures
    }
}
