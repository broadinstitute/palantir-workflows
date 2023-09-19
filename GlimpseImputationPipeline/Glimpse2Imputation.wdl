version 1.0

workflow Glimpse2Imputation {
    input {
        # List of files, one per line
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

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
        
        Int preemptible = 1
        String docker = "us.gcr.io/broad-dsde-methods/glimpse:palantir-workflows_20c9de0"
        Int cpu_phase = 4
        Int mem_gb_phase = 8
        Int cpu_ligate = 4
        Int mem_gb_ligate = 4
        File? monitoring_script
    }

    scatter (reference_chunk in read_lines(reference_chunks)) {
        call GlimpsePhase {
            input:
                reference_chunk = reference_chunk,
                input_vcf = input_vcf,
                input_vcf_index = input_vcf_index,
                impute_reference_only_variants = impute_reference_only_variants,
                n_burnin = n_burnin,
                n_main = n_main,
                effective_population_size = effective_population_size,
                call_indels = call_indels,
                crams = crams,
                cram_indices = cram_indices,
                sample_ids = sample_ids,
                fasta = fasta,
                fasta_index = fasta_index,
                preemptible = preemptible,
                docker = docker,
                cpu = cpu_phase,
                mem_gb = mem_gb_phase,
                monitoring_script = monitoring_script
        }
    }

    call GlimpseLigate {
        input:
            imputed_chunks = GlimpsePhase.imputed_vcf,
            imputed_chunks_indices = GlimpsePhase.imputed_vcf_index,
            output_basename = output_basename,
            ref_dict = ref_dict,
            preemptible = preemptible,
            docker = docker,
            cpu = cpu_ligate,
            mem_gb = mem_gb_ligate,
            monitoring_script = monitoring_script
    }

    output {
        File imputed_vcf = GlimpseLigate.imputed_vcf
        File imputed_vcf_index = GlimpseLigate.imputed_vcf_index
        Array[File?] glimpse_phase_monitoring = GlimpsePhase.monitoring
        File? glimpse_ligate_monitoring = GlimpseLigate.monitoring
    }
}

task GlimpsePhase {
    input {
        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        Array[String] sample_ids
        File? fasta
        File? fasta_index
        File reference_chunk

        Boolean impute_reference_only_variants
        Boolean call_indels
        Int? n_burnin
        Int? n_main
        Int? effective_population_size

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(2.2 * size(input_vcf, "GiB") + size(reference_chunk, "GiB") + 100)
        Int preemptible = 1
        Int max_retries = 3
        String docker
        File? monitoring_script
    }

    parameter_meta {
        crams: {
                        localization_optional: true
                    }
        cram_indices: {
                        localization_optional: true
                    }
    }

    String bam_file_list_input = if defined(crams) then "--bam-list crams.list" else ""
    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(/root/google-cloud-sdk/bin/gcloud auth application-default print-access-token)

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        cram_paths=( ~{sep=" " crams} )
        sample_ids=( ~{sep=" " sample_ids} )

        for i in "${!cram_paths[@]}" ; do
            echo -e "${cram_paths[$i]} ${sample_ids[$i]}" >> crams.list
        done

        /bin/GLIMPSE2_phase \
        ~{"--input-gl " + input_vcf} \
        --reference ~{reference_chunk} \
        --output phase_output.bcf \
        --threads ~{cpu} \
        ~{if impute_reference_only_variants then "--impute-reference-only-variants" else ""} ~{if call_indels then "--call-indels" else ""} \
        ~{"--burnin " + n_burnin} ~{"--main " + n_main} \
        ~{"--ne " + effective_population_size} \
        ~{bam_file_list_input} \
        ~{"--fasta " + fasta}
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        File imputed_vcf = "phase_output.bcf"
        File imputed_vcf_index = "phase_output.bcf.csi"
        File? monitoring = "monitoring.log"
    }
}

task GlimpseLigate {
    input {
        Array[File] imputed_chunks
        Array[File] imputed_chunks_indices
        String output_basename
        File ref_dict

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(2.2 * size(imputed_chunks, "GiB") + 100)
        Int preemptible = 1
        Int max_retries = 3
        String docker
        File? monitoring_script
    }

    command <<<
        set -xeuo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        NPROC=$(nproc)
        echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."
        
        /bin/GLIMPSE2_ligate --input ~{write_lines(imputed_chunks)} --output ligated.vcf.gz --threads ${NPROC}

        # Set correct reference dictionary
        bcftools view -h --no-version ligated.vcf.gz > old_header.vcf        
        java -jar /picard.jar UpdateVcfSequenceDictionary -I old_header.vcf --SD ~{ref_dict} -O new_header.vcf        
        bcftools reheader -h new_header.vcf -o ~{output_basename}.imputed.vcf.gz ligated.vcf.gz
        tabix ~{output_basename}.imputed.vcf.gz
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        File imputed_vcf = "~{output_basename}.imputed.vcf.gz"
        File imputed_vcf_index = "~{output_basename}.imputed.vcf.gz.tbi"
        File? monitoring = "monitoring.log"
    }
}
