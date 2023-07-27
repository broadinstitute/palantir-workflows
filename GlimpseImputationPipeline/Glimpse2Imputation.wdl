version 1.0

workflow Glimpse2Imputation {
    input {
        # List of files, one per line
        File reference_chunks

        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        Array[String]? sample_ids
        File? fasta
        File? fasta_index
        String output_basename

        File ref_dict

        Boolean impute_reference_only_variants
        Boolean call_indels
        Int? n_burnin
        Int? n_main
        
        Int preemptible = 1
        String docker = "us.gcr.io/broad-dsde-methods/glimpse:palantir-workflows_94a0cdd"
        Int cpu_phase = 4
        Int mem_gb_phase = 64
        Int cpu_ligate = 4
        Int mem_gb_ligate = 4
        File? monitoring_script
    }

    scatter (reference_chunk in read_lines(reference_chunks)) {
        call GetNumberOfSitesInChunk {
            input:
                reference_chunk = reference_chunk
        }

        Int n_rare = GetNumberOfSitesInChunk.n_rare
        Int n_common = GetNumberOfSitesInChunk.n_common
        Int n_sites = n_rare + n_common

        if (defined(input_vcf)) {
            call CountSamples {
                input:
                    vcf = select_first([input_vcf])
            }
        }

        Int n_samples = select_first([CountSamples.nSamples, length(select_first([crams]))])

        call SelectResourceParameters {
            input:
                n_rare = n_rare,
                n_common = n_common,
                n_samples = n_samples
        }

        if (SelectResourceParameters.memory_gb > 256 || SelectResourceParameters.request_n_cpus > 32) {
            # force failure if we're accidently goind to request too much resources and spend too much money
            Int safety_check_memory_gb = -1
            Int safety_check_n_cpu = -1
        }

        call GlimpsePhase {
            input:
                reference_chunk = reference_chunk,
                input_vcf = input_vcf,
                input_vcf_index = input_vcf_index,
                impute_reference_only_variants = impute_reference_only_variants,
                n_burnin = n_burnin,
                n_main = n_main,
                call_indels = call_indels,
                crams = crams,
                cram_indices = cram_indices,
                sample_ids = sample_ids,
                fasta = fasta,
                fasta_index = fasta_index,
                preemptible = preemptible,
                docker = docker,
                cpu = select_first([safety_check_memory_gb, SelectResourceParameters.request_n_cpus]),
                mem_gb = select_first([safety_check_n_cpu, SelectResourceParameters.memory_gb]),
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
        Array[String]? sample_ids
        File? fasta
        File? fasta_index
        File reference_chunk

        Boolean impute_reference_only_variants
        Boolean call_indels
        Int? n_burnin
        Int? n_main

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

        function check_for_oom {
            if [ $? -eq 137]; then
                echo 'OutOfMemory' 1>&2
            fi
        }

        trap check_for_oom EXIT

        export GCS_OAUTH_TOKEN=$(/root/google-cloud-sdk/bin/gcloud auth application-default print-access-token)

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        cram_paths=( ~{sep=" " crams} )
        sample_ids=( ~{sep=" " sample_ids} )

        for i in "${!cram_paths[@]}" ; do
            echo -e "${cram_paths[$i]} ${sample_ids[$i]}" >> crams.list
        done

        /GLIMPSE/GLIMPSE2_phase \
        ~{"--input-gl " + input_vcf} \
        --reference ~{reference_chunk} \
        --output phase_output.bcf \
        --threads ~{cpu} \
        ~{if impute_reference_only_variants then "--impute-reference-only-variants" else ""} ~{if call_indels then "--call-indels" else ""} \
        ~{"--burnin " + n_burnin} ~{"--main " + n_main} \
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
        File? picard_jar_override
        File? monitoring_script
    }

    command <<<
        set -xeuo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        NPROC=$(nproc)
        echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."
        
        /GLIMPSE/GLIMPSE2_ligate --input ~{write_lines(imputed_chunks)} --output ligated.vcf.gz --threads ${NPROC}

        # Set correct reference dictionary
        ~{"mv " + picard_jar_override + " /picard.jar"}
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

task GetNumberOfSitesInChunk {
    input {
        File reference_chunk

        String docker = "us.gcr.io/broad-dsde-methods/glimpse:extract_num_sites_from_reference_chunk"
        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(size(reference_chunk, "GiB") + 100)
        Int preemptible = 1
        Int max_retries = 3
    }

    command <<<
        set -xeuo pipefail
        /GLIMPSE/GLIMPSE2_extract_num_sites_from_reference_chunk ~{reference_chunk} > n_sites.txt
        grep "Lrare" n_sites.txt | sed 's/Lrare=//' > n_rare.txt
        grep "Lcommon" n_sites.txt | sed 's/Lcommon=//' > n_common.txt
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
        Int n_rare = read_int("n_rare.txt")
        Int n_common = read_int("n_common.txt")
    }
}

task CountSamples {
  input {
    File vcf

    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 3000
    Int disk_size_gb = 100 + ceil(size(vcf, "GiB"))
  }

  command <<<
    bcftools query -l ~{vcf} | wc -l
  >>>

  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
  }
  output {
    Int nSamples = read_int(stdout())
  }
}

task SelectResourceParameters {
    input {
        Int n_rare
        Int n_common
        Int n_samples
    }

    command <<<
        python3 << EOF
        import math
        n_rare = ~{n_rare}
        n_common = ~{n_common}
        n_samples = ~{n_samples}
        n_sites = n_common + n_rare

        # try to keep expected runtime under 4 hours, but don't ask for more than 32 cpus, or 256 GB memory
        estimated_needed_threads = min(math.ceil(5e-6*n_sites*n_samples/240), 32)
        estimated_needed_memory_gb = min(math.ceil(800e-3 + 0.97e-6 * n_rare * estimated_needed_threads + 14.6e-6 * n_common * estimated_needed_threads + 6.5e-9 * (n_rare + n_common) * n_samples + 13.7e-3 * n_samples + 1.8e-6*(n_rare + n_common)*math.log(n_samples)), 256)
        # add 20% buffer
        estimated_needed_memory_gb = math.ceil(1.2 * estimated_needed_memory_gb)

        with open("n_cpus_request.txt", "w") as f_cpus_request:
            f_cpus_request.write(f'{int(estimated_needed_threads)}')

        with open("memory_gb.txt", "w") as f_mem:
            f_mem.write(f'{int(estimated_needed_memory_gb)}')
        EOF
    >>>

    runtime {
        docker : "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
    }

    output {
        Int memory_gb = read_int("memory_gb.txt")
        Int request_n_cpus = read_int("n_cpus_request.txt")
    }
}