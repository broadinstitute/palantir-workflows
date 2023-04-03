version 1.0

import "Glimpse2Imputation.wdl" as Glimpse2ImputationImport
import "CollectBGEImputationMetricsCohort.wdl" as CollectBGEImputationMetricsCohort

workflow ReferencePanelDev {
    input {
        File reference_panel_chr20_vcf
        File reference_panel_chr20_vcf_index
        File genetic_map_chr20

        File subset_intervals_bed

        File input_vcf
        File input_vcf_index

        File ref_dict

        String glimpse_docker
    }

    call BcftoolsTask {
        input:
            input_vcf = reference_panel_chr20_vcf,
            input_vcf_index = reference_panel_chr20_vcf_index,
            subset_intervals_bed = subset_intervals_bed
    }

    call GlimpseSplitReferenceTask {
        input:
            contig = "chr20",
            i_contig = 0,
            reference_panel = BcftoolsTask.output_vcf,
            reference_panel_index = BcftoolsTask.output_vcf_index,
            genetic_map = genetic_map_chr20,
            min_window_cm = 20,
            docker = glimpse_docker
    }

    call Glimpse2ImputationImport.Glimpse2Imputation as Imputation {
        input:
            reference_chunks = GlimpseSplitReferenceTask.split_reference_chunks_file,
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            docker = glimpse_docker,
            ref_dict = ref_dict
    }

    output {
        File imputed_vcf = Imputation.imputed_vcf
        File imputed_vcf_index = Imputation.imputed_vcf_index
    }
}

task BcftoolsTask {
    input {
        File input_vcf
        File input_vcf_index
        
        File subset_intervals_bed

        Int mem_gb = 4
        Int preemptible = 1
    }

    String output_basename = basename(input_vcf, ".vcf.gz")

    command <<<
        bcftools view -R ~{subset_intervals_bed} -O z -o ~{output_basename}.subset.vcf.gz ~{input_vcf}
        bcftools index -t ~{output_basename}.subset.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.0"
        preemptible: preemptible
        cpu: 2
		disks: "local-disk 100 HDD"
		memory: mem_gb + " GB"
    }

    output {
        File output_vcf = "~{output_basename}.subset.vcf.gz"
        File output_vcf_index = "~{output_basename}.subset.vcf.gz.tbi"
    }
}

task GlimpseSplitReferenceTask {
    input {
        String contig
        Int i_contig
        File reference_panel
        File reference_panel_index
        File genetic_map

        Int? seed
        Int? min_window_cm

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(2.2 * size(reference_panel, "GiB") + size(genetic_map, "GiB") + 100)
        Int preemptible = 1
        String docker
        File? monitoring_script
    }

    String reference_output_dir = "reference_output_dir"

    command <<<
        set -xeuo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        NPROC=$(nproc)
        echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."

        # Print chunk index to variable
        CONTIGINDEX=$(printf "%04d" ~{i_contig})

        /GLIMPSE/GLIMPSE2_chunk --input ~{reference_panel} --region ~{contig} --map ~{genetic_map} --sequential --threads ${NPROC} --output chunks_contigindex_${CONTIGINDEX}.txt ~{"--seed "+seed} ~{"--window-cm "+min_window_cm}

        mkdir -p ~{reference_output_dir}

        I_CHUNK=0
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            # Extract coordinates from chunks.txt file
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)

            # Print chunk index to variable
            CHUNKINDEX=$(printf "%04d" $I_CHUNK)

            /GLIMPSE/GLIMPSE2_split_reference --threads ${NPROC} --reference ~{reference_panel} --map ~{genetic_map} --input-region ${IRG} --output-region ${ORG} --output ~{reference_output_dir}/reference_panel_contigindex_${CONTIGINDEX}_chunkindex_${CHUNKINDEX} ~{"--seed "+seed}

            # Increase i (and make sure the exit code is zero)
            (( I_CHUNK++ )) || true
        done < chunks_contigindex_${CONTIGINDEX}.txt
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        Array[File] split_reference_chunks = glob(reference_output_dir + "/*.bin")
        File split_reference_chunks_file = write_lines(glob(reference_output_dir + "/*.bin"))

        # We don't know the exact filename of the chunks.txt file since we need to add leading zeros to the contigindex. Since WDL doesn't
        # have a built-in way to do that, we have to rely on the command section to do that. However, we don't have access to that bash
        # variable in the output section, so we have to use glob here and return the first (and only) result.
        File chunks = glob("chunks_contigindex_*.txt")[0]

        File? monitoring = "monitoring.log"
    }
}