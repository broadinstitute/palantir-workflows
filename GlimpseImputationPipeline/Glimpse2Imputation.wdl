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

        Boolean collect_qc_metrics = true
        
        Int preemptible = 9
        String docker = "us.gcr.io/broad-dsde-methods/ckachulis/glimpse_for_wdl_pipeline:checkpointing_and_extract_num_sites"
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
            # force failure if we're accidently going to request too much resources and spend too much money
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
                effective_population_size = effective_population_size,
                call_indels = call_indels,
                crams = crams,
                cram_indices = cram_indices,
                sample_ids = sample_ids,
                fasta = fasta,
                fasta_index = fasta_index,
                preemptible = preemptible,
                docker = docker,
                cpu = select_first([safety_check_n_cpu, SelectResourceParameters.request_n_cpus]),
                mem_gb = select_first([safety_check_memory_gb, SelectResourceParameters.memory_gb]),
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

    if (collect_qc_metrics) {
        call CollectQCMetrics {
            input:
                imputed_vcf = GlimpseLigate.imputed_vcf,
                output_basename = output_basename,
                monitoring_script = monitoring_script
        }
    }

    output {
        File imputed_vcf = GlimpseLigate.imputed_vcf
        File imputed_vcf_index = GlimpseLigate.imputed_vcf_index
        
        File? qc_metrics = CollectQCMetrics.qc_metrics

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
        Int disk_size_gb = ceil(2.2 * size(input_vcf, "GiB") + size(reference_chunk, "GiB") + 10)
        Int preemptible = 9
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

        duplicate_cram_filenames=$(printf "%s\n" "${cram_paths[@]}" | xargs -I {} basename {} | uniq -d)
        if [ ! -z "$duplicate_cram_filenames" ]; then
            echo "The input CRAMs contain multiple files with the same basename, which leads to an error due to the way that htslib is implemented. Duplicate filenames:"
            echo $duplicate_cram_filenames
            exit 1
        fi

        for i in "${!cram_paths[@]}" ; do
            echo -e "${cram_paths[$i]} ${sample_ids[$i]}" >> crams.list
        done

        cmd="/bin/GLIMPSE2_phase \
        ~{"--input-gl " + input_vcf} \
        --reference ~{reference_chunk} \
        --output phase_output.bcf \
        --threads ~{cpu} \
        ~{if impute_reference_only_variants then "--impute-reference-only-variants" else ""} ~{if call_indels then "--call-indels" else ""} \
        ~{"--burnin " + n_burnin} ~{"--main " + n_main} \
        ~{"--ne " + effective_population_size} \
        ~{bam_file_list_input} \
        ~{"--fasta " + fasta} \
        --checkpoint-file-out checkpoint.bin"

        if [ -s "checkpoint.bin" ]; then
            cmd="$cmd --checkpoint-file-in checkpoint.bin" 
        fi

        eval $cmd
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        checkpointFile: "checkpoint.bin"
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

task CollectQCMetrics {
    input {
        File imputed_vcf
        String output_basename
        
        Int preemptible = 1
        String docker = "hailgenetics/hail:0.2.126-py3.11"
        Int cpu = 4
        Int mem_gb = 16
        File? monitoring_script
    }

    parameter_meta {
        imputed_vcf: {
                        localization_optional: true
                    }
    }

    Int disk_size_gb = 100
    
    command <<<
        set -euo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        cat <<'EOF' > script.py
import hail as hl
import pandas as pd

# Calculate metrics
hl.init(default_reference='GRCh38', idempotent=True)
vcf = hl.import_vcf('~{imputed_vcf}', force_bgz=True)
qc = hl.sample_qc(vcf)
qc.cols().flatten().export('~{output_basename}.hail_qc_metrics.tsv')

qc_pandas = pd.read_csv('~{output_basename}.hail_qc_metrics.tsv', sep='\t')
qc_pandas_renamed = qc_pandas.rename(columns={col: col.replace('sample_qc.', '') for col in qc_pandas.columns if 'sample_qc.' in col})
qc_pandas_renamed.to_csv('~{output_basename}.qc_metrics.tsv', sep='\t', index=False)
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
        Boolean qc_passed = read_boolean("~{output_basename}.qc_passed.txt")
        File qc_failures = "~{output_basename}.qc_failures.txt"
        File? monitoring = "monitoring.log"
    }
}

task GetNumberOfSitesInChunk {
    input {
        File reference_chunk

        String docker = "us.gcr.io/broad-dsde-methods/ckachulis/glimpse_for_wdl_pipeline:checkpointing_and_extract_num_sites "
        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(size(reference_chunk, "GiB") + 10)
        Int preemptible = 1
        Int max_retries = 3
    }

    command <<<
        set -xeuo pipefail
        /bin/GLIMPSE2_extract_num_sites_from_reference_chunk ~{reference_chunk} > n_sites.txt
        cat n_sites.txt
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
    Int disk_size_gb = 10 + ceil(size(vcf, "GiB"))
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
        estimated_needed_memory_gb = min(math.ceil((800e-3 + 0.97e-6 * n_rare * estimated_needed_threads + 14.6e-6 * n_common * estimated_needed_threads + 6.5e-9 * (n_rare + n_common) * n_samples + 13.7e-3 * n_samples + 1.8e-6*(n_rare + n_common)*math.log(n_samples))), 256)
        # recalc allowable threads, may be some additional threads available due to rounding memory up
        threads_to_use = max(math.floor((estimated_needed_memory_gb - (800e-3 + 6.5e-9 * (n_rare + n_common) * n_samples + 13.7e-3 * n_samples + 1.8e-6*(n_rare + n_common)*math.log(n_samples)))/(0.97e-6 * n_rare + 14.6e-6 * n_common)), 1) 
        #estimated_needed_memory_gb = math.ceil(1.2 * estimated_needed_memory_gb)

        with open("n_cpus_request.txt", "w") as f_cpus_request:
            f_cpus_request.write(f'{int(threads_to_use)}')

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
