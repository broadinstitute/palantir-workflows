version 1.0
import "Glimpse2Imputation.wdl" as GLIMPSE
import "CollectBGEImputationMetrics.wdl" as Eval

workflow GlimpseCheckpointTest {
    input {
        File reference_chunks

        Array[File] crams
        Array[File] cram_indices
        Array[String] sample_ids
        File fasta
        File fasta_index

        File ref_dict
        String docker = "us.gcr.io/broad-dsde-methods/ckachulis/glimpse_for_wdl_pipeline:checkpointing_tests"
        Int cpu_phase = 4
        Int mem_gb_phase = 4
        Int cpu_ligate = 4
        Int mem_gb_ligate = 4

        File truth_vcf
        Array[String] truth_sample_ids
        String? intervals
    }

    Array[Int] n_burnin_array = [0,1,3,5]
    Array[Int] n_main_array = [1,5,10]

    scatter (id in sample_ids) {
        String ancestries = "unknown"
    }

    Map[String,String] ancestry_annotation_map = {"unknown" : "RAF"}

    scatter (reference_chunk in read_lines(reference_chunks))
    {

        call GlimpsePhase {
                input:
                    reference_chunk = reference_chunk,
                    crams = crams,
                    cram_indices = cram_indices,
                    sample_ids = sample_ids,
                    fasta = fasta,
                    fasta_index = fasta_index,
                    preemptible = 0,
                    docker = docker,
                    cpu = cpu_phase,
                    mem_gb = mem_gb_phase
            }
        
    }

    call GLIMPSE.GlimpseLigate{
        input:
            imputed_chunks = GlimpsePhase.imputed_vcf,
            imputed_chunks_indices = GlimpsePhase.imputed_vcf_index,
            output_basename = "no_checkpoint",
            ref_dict = ref_dict,
            preemptible = 0,
            docker = docker,
            cpu = cpu_ligate,
            mem_gb = mem_gb_ligate
    }

    scatter (n_burnin in n_burnin_array) {
        scatter (reference_chunk in read_lines(reference_chunks))
            {
                call GlimpsePhase as GlimpseIntoCheckpointBurnin {
                        input:
                            reference_chunk = reference_chunk,
                            n_burnin = n_burnin,
                            n_main = 0,
                            crams = crams,
                            cram_indices = cram_indices,
                            sample_ids = sample_ids,
                            fasta = fasta,
                            fasta_index = fasta_index,
                            preemptible = 0,
                            docker = docker,
                            cpu = cpu_phase,
                            mem_gb = mem_gb_phase
                    }

                call GlimpsePhase as GlimpseFromCheckpointBurnin {
                        input:
                            reference_chunk = reference_chunk,
                            crams = crams,
                            cram_indices = cram_indices,
                            sample_ids = sample_ids,
                            fasta = fasta,
                            fasta_index = fasta_index,
                            preemptible = 0,
                            docker = docker,
                            cpu = cpu_phase,
                            mem_gb = mem_gb_phase,
                            checkpointFileIn = GlimpseIntoCheckpointBurnin.checkpointFileOut
                    }
                
            }

        call GLIMPSE.GlimpseLigate as GlimpseLigateBurnin {
            input:
                imputed_chunks = GlimpseFromCheckpointBurnin.imputed_vcf,
                imputed_chunks_indices = GlimpseFromCheckpointBurnin.imputed_vcf_index,
                output_basename = "from_checkpoint_n_burnin_"+n_burnin,
                ref_dict = ref_dict,
                preemptible = 0,
                docker = docker,
                cpu = cpu_ligate,
                mem_gb = mem_gb_ligate
        }

        call Eval.PearsonCorrelationByAF as CorrelationBurnin{
            input:
                evalVcf = GlimpseLigateBurnin.imputed_vcf,
                truthVcf = truth_vcf,
                eval_sample_ids = sample_ids,
                truth_sample_ids = truth_sample_ids,
                ancestry_to_af_annotation_map = ancestry_annotation_map,
                ancestries = ancestries,
                af_resource = GlimpseLigateBurnin.imputed_vcf,
                intervals = intervals,
                output_basename = "correlation_from_checkpoint_n_burnin_"+n_burnin
        }

        call Eval.PearsonCorrelationByAF as CorrelationBurninVsNoCheckpoint {
            input:
                evalVcf = GlimpseLigateBurnin.imputed_vcf,
                truthVcf = GlimpseLigate.imputed_vcf,
                eval_sample_ids = sample_ids,
                truth_sample_ids = sample_ids,
                ancestry_to_af_annotation_map = ancestry_annotation_map,
                ancestries = ancestries,
                af_resource = GlimpseLigateBurnin.imputed_vcf,
                intervals = intervals,
                output_basename = "correlation_vs_no_checkpoin_from_checkpoint_n_burnin_"+n_burnin
        }

        call AnnotateCorrelation as AnnotateCorrelationBurnin {
            input:
                correlation = CorrelationBurnin.correlations,
                annotation = "from_checkpoint_burnin_"+n_burnin
        }
    }

    scatter (n_main in n_main_array) {
        scatter (reference_chunk in read_lines(reference_chunks))
            {
                call GlimpsePhase as GlimpseIntoCheckpointMain {
                        input:
                            reference_chunk = reference_chunk,
                            n_main = n_main,
                            crams = crams,
                            cram_indices = cram_indices,
                            sample_ids = sample_ids,
                            fasta = fasta,
                            fasta_index = fasta_index,
                            preemptible = 0,
                            docker = docker,
                            cpu = cpu_phase,
                            mem_gb = mem_gb_phase
                    }

                call GlimpsePhase as GlimpseFromCheckpointMain {
                        input:
                            reference_chunk = reference_chunk,
                            crams = crams,
                            cram_indices = cram_indices,
                            sample_ids = sample_ids,
                            fasta = fasta,
                            fasta_index = fasta_index,
                            preemptible = 0,
                            docker = docker,
                            cpu = cpu_phase,
                            mem_gb = mem_gb_phase,
                            checkpointFileIn = GlimpseIntoCheckpointMain.checkpointFileOut
                    }
                
            }

        call GLIMPSE.GlimpseLigate as GlimpseLigateMain {
            input:
                imputed_chunks = GlimpseFromCheckpointMain.imputed_vcf,
                imputed_chunks_indices = GlimpseFromCheckpointMain.imputed_vcf_index,
                output_basename = "from_checkpoint_n_main_"+n_main,
                ref_dict = ref_dict,
                preemptible = 0,
                docker = docker,
                cpu = cpu_ligate,
                mem_gb = mem_gb_ligate
        }

        call Eval.PearsonCorrelationByAF as CorrelationMain{
            input:
                evalVcf = GlimpseLigateMain.imputed_vcf,
                truthVcf = truth_vcf,
                eval_sample_ids = sample_ids,
                truth_sample_ids = truth_sample_ids,
                ancestry_to_af_annotation_map = ancestry_annotation_map,
                ancestries = ancestries,
                af_resource = GlimpseLigateMain.imputed_vcf,
                intervals = intervals,
                output_basename = "correlation_from_checkpoint_n_main_"+n_main
        }

        call Eval.PearsonCorrelationByAF as CorrelationMainVsNoCheckpoint {
            input:
                evalVcf = GlimpseLigateMain.imputed_vcf,
                truthVcf = GlimpseLigate.imputed_vcf,
                eval_sample_ids = sample_ids,
                truth_sample_ids = sample_ids,
                ancestry_to_af_annotation_map = ancestry_annotation_map,
                ancestries = ancestries,
                af_resource = GlimpseLigateMain.imputed_vcf,
                intervals = intervals,
                output_basename = "correlation_vs_no_checkpoint_from_checkpoint_n_burnin_"+n_main
        }

        call AnnotateCorrelation as AnnotateCorrelationMain {
            input:
                correlation = CorrelationMain.correlations,
                annotation = "from_checkpoint_main_"+n_main
        }
    }

    call Eval.PearsonCorrelationByAF {
        input:
            evalVcf = GlimpseLigate.imputed_vcf,
            truthVcf = truth_vcf,
            eval_sample_ids = sample_ids,
            truth_sample_ids = truth_sample_ids,
            ancestry_to_af_annotation_map = ancestry_annotation_map,
            ancestries = ancestries,
            af_resource = GlimpseLigate.imputed_vcf,
            intervals = intervals,
            output_basename = "correlation_no_checkpoint"
    }

    call AnnotateCorrelation {
            input:
                correlation = PearsonCorrelationByAF.correlations,
                annotation = "no_checkpoint"
    }

    call CombineCorrelations {
        input:
            correlations = flatten([[AnnotateCorrelation.annotated_correlation], AnnotateCorrelationMain.annotated_correlation, AnnotateCorrelationBurnin.annotated_correlation])
    }

    output {
        File combined_correlations = CombineCorrelations.combined_correlations
        File correlations_plot = CombineCorrelations.correlations_plot
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

        Int? n_burnin
        Int? n_main
        Int? effective_population_size

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(2.2 * size(input_vcf, "GiB") + size(reference_chunk, "GiB") + 100)
        Int preemptible = 9
        Int max_retries = 3
        String docker
        File? monitoring_script
        File? checkpointFileIn
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

        touch checkpoint.bin

        cmd="/bin/GLIMPSE2_phase \
        ~{"--input-gl " + input_vcf} \
        --reference ~{reference_chunk} \
        --output phase_output.bcf \
        --threads ~{cpu} \
        ~{"--burnin " + n_burnin} ~{"--main " + n_main} \
        ~{"--ne " + effective_population_size} \
        ~{bam_file_list_input} \
        ~{"--fasta " + fasta} \
        --checkpoint-file-out checkpoint.bin"

        ~{"cp " + checkpointFileIn + " checkpoint.bin"}
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
    }

    output {
        File imputed_vcf = "phase_output.bcf"
        File imputed_vcf_index = "phase_output.bcf.csi"
        File? monitoring = "monitoring.log"
        File checkpointFileOut = "checkpoint.bin"
    }
}

task AnnotateCorrelation {
    input {
        File correlation
        String annotation
    }

    String basename = basename(correlation, ".tsv")

    command <<<
        Rscript - <<- "EOF"
        library(dplyr)
        library(readr)

        cor <- read_tsv("~{correlation}", skip=6) %>% mutate(annotation="~{annotation}")

        cor %>% write_tsv("~{basename}.anotate.tsv")

        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    }

    output {
        File annotated_correlation = "~{basename}.anotate.tsv"
    }
}

task CombineCorrelations {
    input {
        Array[File] correlations
    }

    command <<<
        Rscript - <<- "EOF"
        library(dplyr)
        library(readr)
        library(purrr)
        library(tidyr)
        library(ggplot2)

        cor <- c("~{sep='","' correlations}") %>% map(read_tsv) %>% reduce(bind_rows)
        cor %>% write_tsv("correlations.tsv")

        cor %>% pivot_longer(cols=c("SNP_CORRELATION", "INDEL_CORRELATION"),
                     names_sep="_",
                     names_to=c("variant_type", NA),
                     values_to = "correlation") %>%
                filter(!is.na(correlation)) %>%
                ggplot(aes(x=BIN_CENTER, y=correlation^2)) +
                    geom_line(aes(color=annotation)) +
                    facet_grid(SAMPLE~variant_type) +
                    scale_x_log10() +
                    theme_bw()

        ggsave(filename = "checkpoint_test_correlation_plot.png", dpi=300, width = 6, height = 6)

        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse@sha256:aaace6c41a258e13da76881f0b282932377680618fcd5d121583f9455305e727"
    }

    output {
        File combined_correlations = "correlations.tsv"
        File correlations_plot = "checkpoint_test_correlation_plot.png"
    }
}

