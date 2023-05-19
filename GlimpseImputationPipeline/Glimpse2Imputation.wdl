version 1.0 

workflow Glimpse2Imputation {
    input {
        # List of files, one per line
        File reference_chunks

        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        File? fasta
        File? fasta_index
        String output_basename

        File ref_dict
        
        Int preemptible = 1
        String docker = "us.gcr.io/broad-dsde-methods/glimpse:2.0.0"
        Int cpu_phase = 4
        Int mem_gb_phase = 64
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
                crams = crams,
                cram_indices = cram_indices,
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
        File? fasta
        File? fasta_index
        File reference_chunk

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(2.2 * size(input_vcf, "GiB") + 1.5 * size(select_first([crams, []]), "GiB") + size(reference_chunk, "GiB") + 100)
        Int preemptible = 1
        Int max_retries = 3
        String docker
        File? monitoring_script
    }

    String bam_file_list_input = if defined(crams) then "--bam-list crams.list" else ""
    command <<<
        set -xeuo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}


        #NPROC=$(nproc)
        #echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."

        echo -e "~{sep="\n" crams}" > crams.list

        cram_paths=( ~{sep=" " crams} )
        cram_index_paths=( ~{sep=" " cram_indices} )

        for i in "${!cram_paths[@]}" ; do
            if [[ "${cram_index_paths[$i]}" != "${cram_paths[$i]}.crai" ]]; then
                mv "${cram_index_paths[$i]}" "${cram_paths[$i]}.crai"
            fi
        done

        /GLIMPSE/GLIMPSE2_phase \
        ~{"--input-gl " + input_vcf} \
        --reference ~{reference_chunk} \
        --output phase_output.bcf \
        --threads ~{cpu} \
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
        java -jar /picard.jar UpdateVcfSequenceDictionary -I ligated.vcf.gz --SD ~{ref_dict} -O ~{output_basename}.imputed.vcf.gz
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