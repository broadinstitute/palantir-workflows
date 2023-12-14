version 1.0

workflow GlimpseImputation {
    input {
        Array[File] chunks
        Array[String] reference_panel_contig_names
        Array[String] genetic_map_contig_names
        String reference_panel_prefix
        String reference_panel_suffix
        String reference_panel_index_suffix
        String genetic_map_path_prefix
        String genetic_map_path_suffix
        File input_vcf
        File input_vcf_index
        File ref_dict
        
        Int mem_gb = 4
        Int cpu = 4
        Int preemptible = 1
        File? monitoring_script
    }

    scatter (chunks_and_contig in zip(chunks, zip(reference_panel_contig_names, genetic_map_contig_names))) {
        String reference_filename = reference_panel_prefix + chunks_and_contig.right.left + reference_panel_suffix
        String genetic_map_filename = genetic_map_path_prefix + chunks_and_contig.right.right + genetic_map_path_suffix

        call Glimpse {
            input:
                reference_panel = reference_filename,
                reference_panel_index = reference_filename + reference_panel_index_suffix,
                contig_name = chunks_and_contig.right.right,
                chunks = chunks_and_contig.left,
                ref_dict = ref_dict,
                input_vcf = input_vcf,
                input_vcf_index = input_vcf_index,
                genetic_map = genetic_map_filename,
                mem_gb = mem_gb,
                cpu = cpu,
                preemptible = preemptible,
                monitoring_script = monitoring_script
        }
    }

    call GatherVcfs {
        input:
            input_vcfs = Glimpse.imputed_vcf,
            input_vcf_indices = Glimpse.imputed_vcf_index,
            output_vcf_basename = basename(input_vcf, ".vcf.gz") + ".imputed",
            preemptible = preemptible,
            monitoring_script = monitoring_script
    }

    output {
        File imputed_vcf = GatherVcfs.output_vcf
        File imputed_vcf_index = GatherVcfs.output_vcf_index
        Array[File?] glimpse_monitoring = Glimpse.monitoring
        File? gather_monitoring = GatherVcfs.monitoring
    }
}

task Glimpse {
    input {
        String contig_name
        File input_vcf
        File input_vcf_index
        File reference_panel
        File reference_panel_index
        File ref_dict
        File genetic_map

        File chunks

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(2 * size(input_vcf, "GiB") + size(reference_panel, "GiB") + size(genetic_map, "GiB") + 100)
        Int preemptible = 1

        File? monitoring_script
    }

    command <<<
        set -xeuo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}
        
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            PHASE_OUT=phased.chunk${ID}.vcf.gz
            /glimpse/phase/bin/GLIMPSE_phase --input ~{input_vcf} --reference ~{reference_panel} --map ~{genetic_map} --input-region ${IRG} --output-region ${ORG} --output ${PHASE_OUT} --thread ~{cpu}
            tabix ${PHASE_OUT}
        done < ~{chunks}

        LST=list.txt
        ls phased.chunk*.vcf.gz > ${LST}
        LIGATE_OUT=ligated.~{contig_name}.vcf.gz
        /glimpse/ligate/bin/GLIMPSE_ligate --input ${LST} --output $LIGATE_OUT --thread ~{cpu}
        tabix ${LIGATE_OUT}

        SAMPLE_OUT=ligated.sampled.~{contig_name}.vcf.gz
        /glimpse/sample/bin/GLIMPSE_sample --input ${LIGATE_OUT} --solve --output ${SAMPLE_OUT} --thread ~{cpu}
        tabix ${SAMPLE_OUT}

        UPDATE_OUT=ligated.sampled.dict_updated.~{contig_name}.vcf.gz
        java -jar /picard.jar UpdateVcfSequenceDictionary -I ${SAMPLE_OUT} --SD ~{ref_dict} -O ${UPDATE_OUT}
        tabix ${UPDATE_OUT}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/glimpse:1.1.1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File imputed_vcf = "ligated.sampled.dict_updated." + contig_name + ".vcf.gz"
        File imputed_vcf_index = "ligated.sampled.dict_updated." + contig_name + ".vcf.gz.tbi"
        File? monitoring = "monitoring.log"
    }
}

task GatherVcfs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_vcf_basename

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
        Int cpu = 1
        Int memory_mb = 16000
        Int disk_size_gb = ceil(3*size(input_vcfs, "GiB"))
        Int preemptible = 1

        File? monitoring_script
    }
    Int command_mem = memory_mb - 1000
    Int max_heap = memory_mb - 500

    command <<<
        set -xeuo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        GatherVcfs \
        -I ~{sep=' -I ' input_vcfs} \
        -O ~{output_vcf_basename}.vcf.gz

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        IndexFeatureFile -I ~{output_vcf_basename}.vcf.gz

    >>>
    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: preemptible
    }
    output {
        File output_vcf = "~{output_vcf_basename}.vcf.gz"
        File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
        File? monitoring = "monitoring.log"
    }
}
