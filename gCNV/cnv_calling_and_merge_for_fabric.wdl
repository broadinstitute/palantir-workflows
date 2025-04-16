version 1.0

import "single_sample_cnv_germline_case_filter_workflow.wdl" as cnv_case_and_filter
workflow CNVCallingAndMergeForFabric {
    input {
        File normal_bam
        File normal_bai
        File short_variant_vcf

        File contig_ploidy_model_tar
        File preprocessed_intervals
        File gcnv_model_tar
        Array[File]+ gcnv_panel_genotyped_segments
        Array[File]+ gcnv_panel_copy_ratios
        Array[File]+ gcnv_panel_read_counts

        Float overlap_thresh = 0.5

        String gatk_docker

        Int maximum_number_events_per_sample
        Int maximum_number_pass_events_per_sample
        Int ref_copy_number_autosomal_contigs

        File ref_fasta
        File ref_fasta_fai
        File ref_fasta_dict

        Array[String] allosomal_contigs
    }

    call cnv_case_and_filter.SingleSampleGCNVAndFilterVCFs {
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            contig_ploidy_model_tar = contig_ploidy_model_tar,
            preprocessed_intervals = preprocessed_intervals,
            gcnv_model_tar = gcnv_model_tar,
            pon_genotyped_segments_vcfs = gcnv_panel_genotyped_segments,
            gatk_docker = gatk_docker,
            maximum_number_events_per_sample = maximum_number_events_per_sample,
            maximum_number_pass_events_per_sample = maximum_number_pass_events_per_sample,
            ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            allosomal_contigs = allosomal_contigs,
            overlap_thresh = overlap_thresh
    }

    call ReformatAndMergeForFabric {
        input:
            cnv_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf,
            short_variant_vcf = short_variant_vcf,
            gatk_docker = gatk_docker
    }

    call GCNVVisualzation {
        input:
            filtered_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf,
            case_copy_ratios = SingleSampleGCNVAndFilterVCFs.denoised_copy_ratios,
            case_read_counts = SingleSampleGCNVAndFilterVCFs.read_counts,
            panel_copy_ratios = gcnv_panel_copy_ratios,
            panel_read_counts = gcnv_panel_read_counts,
            interval_list = SingleSampleGCNVAndFilterVCFs.interval_list

    }

    output {
        File filtered_cnv_genotyped_segments_vcf = SingleSampleGCNVAndFilterVCFs.filtered_vcf
        File filtered_cnv_genotyped_segments_vcf_index = SingleSampleGCNVAndFilterVCFs.filtered_vcf_index
        File filtered_cnv_genotyped_segments_vcf_md5sum = SingleSampleGCNVAndFilterVCFs.filtered_vcf_md5sum

        File merged_vcf = ReformatAndMergeForFabric.merged_vcf
        File merged_vcf_index = ReformatAndMergeForFabric.merged_vcf_index
        File merged_vcf_md5sum = ReformatAndMergeForFabric.merged_vcf_md5sum

        Boolean qc_passed = SingleSampleGCNVAndFilterVCFs.qc_passed
        File cnv_metrics = SingleSampleGCNVAndFilterVCFs.cnv_metrics
        File cnv_event_report = GCNVVisualzation.cnv_event_report

    }
}

#Fabric doesn't seem to like ./. genotypes on
task ReformatGCNVForFabric {
    input {
        File cnv_vcf
        Int disk_size_gb = 20
        Int mem_gb = 4
    }

    String output_basename = basename(cnv_vcf, ".filtered.genotyped-segments.vcf.gz")

    command <<<
        set -euo pipefail

        python << CODE
        from pysam import VariantFile

        with VariantFile("~{cnv_vcf}") as cnv_vcf_in:
            header_out = cnv_vcf_in.header
            header_out.info.add("CN", "A", "Integer", "Copy number associated with <CNV> alleles")
            with VariantFile("~{output_basename}.reformatted_for_fabric.vcf.gz",'w', header = header_out) as cnv_vcf_out:
                for rec in cnv_vcf_in.fetch():
                    if 'PASS' in rec.filter:
                        if rec.alts and rec.alts[0] == "<DUP>":
                            for rec_sample in rec.samples.values():
                                ploidy = len(rec_sample.alleles)
                                rec_sample.allele_indices = (None,)*(ploidy - 1) + (1,)
                        cnv_vcf_out.write(rec)
        CODE
    >>>

    runtime {
            docker: "us.gcr.io/broad-dsde-methods/pysam:v1.1"
            preemptible: 3
            cpu: 2
            disks: "local-disk " + disk_size_gb + " HDD"
            memory: mem_gb + " GB"
        }

    output {
        File reformatted_vcf = "~{output_basename}.reformatted_for_fabric.vcf.gz"
    }
}

task MergeVcfs {
    input {
        File cnv_vcf
        File short_variant_vcf

        String gatk_docker
        Int mem_gb=4
        Int disk_size_gb = 100
    }

    String output_basename = basename(short_variant_vcf, ".hard-filtered.vcf.gz")
    String output_cnv_basename = basename(cnv_vcf, ".reformatted_for_fabric.vcf.gz")
    command <<<
        set -euo pipefail

        if [ "~{output_basename}" != "~{output_cnv_basename}" ]; then
            echo "input vcf names do not agree"
            exit 1
        fi

        gatk --java-options "-Dsamjdk.create_md5=true" MergeVcfs -I ~{short_variant_vcf} -I ~{cnv_vcf} -O ~{output_basename}.merged.vcf.gz

        mv ~{output_basename}.merged.vcf.gz.md5 ~{output_basename}.merged.vcf.gz.md5sum
    >>>

    output {
        File merged_vcf = "~{output_basename}.merged.vcf.gz"
        File merged_vcf_index = "~{output_basename}.merged.vcf.gz.tbi"
        File merged_vcf_md5sum = "~{output_basename}.merged.vcf.gz.md5sum"
    }

     runtime {
        docker: gatk_docker
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 5
    }
}

task ReformatAndMergeForFabric {
    input {
            File cnv_vcf
            File short_variant_vcf


            String gatk_docker
            Int mem_gb=4
            Int disk_size_gb = 100
        }

    String output_basename = basename(cnv_vcf, ".filtered.genotyped-segments.vcf.gz")
    String output_short_variant_basename = basename(short_variant_vcf, ".hard-filtered.vcf.gz")

    command <<<
        set -euo pipefail

        if [ "~{output_basename}" != "~{output_short_variant_basename}" ]; then
            echo "input vcf names do not agree"
            exit 1
        fi

        python << CODE
        from pysam import VariantFile

        with VariantFile("~{cnv_vcf}") as cnv_vcf_in:
            header_out = cnv_vcf_in.header
            header_out.info.add("CN", "A", "Integer", "Copy number associated with <CNV> alleles")
            with VariantFile("~{output_basename}.reformatted_for_fabric.vcf.gz",'w', header = header_out) as cnv_vcf_out:
                for rec in cnv_vcf_in.fetch():
                    if 'PASS' in rec.filter:
                        if rec.alts and rec.alts[0] == "<DUP>":
                            for rec_sample in rec.samples.values():
                                ploidy = len(rec_sample.alleles)
                                rec_sample.allele_indices = (None,)*(ploidy - 1) + (1,)
                        cnv_vcf_out.write(rec)
        CODE

        gatk --java-options "-Dsamjdk.create_md5=true" MergeVcfs -I ~{short_variant_vcf} -I ~{output_basename}.reformatted_for_fabric.vcf.gz -O ~{output_basename}.merged.vcf.gz

        mv ~{output_basename}.merged.vcf.gz.md5 ~{output_basename}.merged.vcf.gz.md5sum
    >>>

    output {
        File merged_vcf = "~{output_basename}.merged.vcf.gz"
        File merged_vcf_index = "~{output_basename}.merged.vcf.gz.tbi"
        File merged_vcf_md5sum = "~{output_basename}.merged.vcf.gz.md5sum"
    }

    runtime {
        docker: gatk_docker
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 5
    }
}

task GCNVVisualzation {
    input {
        File filtered_vcf
        File case_copy_ratios
        File case_read_counts
        Array[File]+ panel_copy_ratios
        Array[File]+ panel_read_counts
        File interval_list
        Int mem_gb=4
    }

    String output_prefix = basename(filtered_vcf, ".filtered.genotyped-segments.vcf.gz")

    command <<<
       set -eo pipefail
       gunzip -c ~{filtered_vcf} > ~{output_prefix}.filtered.genotyped-segments.vcf
       cat << EOF > ~{output_prefix}_cnv_event_report.Rmd

       ---
       title: '~{output_prefix}'
       output: html_document
       ---

       \`\`\`{css, echo=FALSE}
       .container-fluid.main-container {
         max-width: 85%;
       }

       .nav {
           position: sticky;
           background-color: #f8f9fa;
           border-bottom: 2px solid #007bff; /* Blue underline */
           top:0;
           z-index: 1;
       }

       .custom-tabset .nav-tabs {
         position: fixed
         max-width: 50%;
         background-color: #f8f9fa;  /* Light gray background */
         border-bottom: 2px solid #007bff; /* Blue underline */
         z-index: 0;
       }

       .tab-content {
         margin-top: 0;
         padding-top: 10px
       }
       \`\`\`

       \`\`\`{r setup, include=FALSE, message=FALSE, warning=FALSE}
       #knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
       library(ggplot2)
       library(dplyr)
       library(readr)
       library(purrr)
       library(scales)
       library(stringr)
       library(tidyr)
       library(rhdf5)
       \`\`\`

       \`\`\`{r load data, include=FALSE, message=FALSE, warning=FALSE}

       panel_tsvs = c("~{sep='","' panel_copy_ratios}")

       read_copy_ratios <- function(file_name) {
         t <- read_tsv(file_name, comment="@")

         lines <- read_lines(file_name)
         sample_name <- lines %>%
            str_subset("^@RG") %>%
            str_extract("SM:([^\\\t]+)")
         t <- t %>% mutate(sample_name=sample_name)
       }
       cr_panel <- panel_tsvs %>% map(read_copy_ratios) %>% reduce(bind_rows)

       panel_extremes <- cr_panel %>% group_by(CONTIG, START, END) %>% summarise(max=max(LINEAR_COPY_RATIO),
                                                               min=min(LINEAR_COPY_RATIO),
                                                               second_max=sort(LINEAR_COPY_RATIO,decreasing=TRUE)[2],
                                                               second_min=sort(LINEAR_COPY_RATIO,decreasing=FALSE)[2])

       cr_case <- read_tsv("~{case_copy_ratios}", comment="@")
       header_lines <- read_lines("~{case_copy_ratios}")
       chromosomes <- header_lines %>%
         str_subset("^@SQ") %>%           # Keep only lines starting with @SQ
         str_extract("(?<=SN:)[^\\\t]+")   # Extract the chromosome name

       cr_case<-cr_case %>% mutate(CONTIG = factor(CONTIG,chromosomes))


       vcf_df <- read_tsv("~{output_prefix}.filtered.genotyped-segments.vcf", comment="#",
                          col_names = c("CONTIG","START","ID","REF","ALT","QUAL","FILTER","INFO")) %>%
       mutate(
         END = as.numeric(str_extract(INFO, "(?<=END=)\\\d+")),
         PANEL_FREQ = as.numeric(str_extract(INFO, "(?<=PANEL_FREQ=)\\\d+(\\\.\\\d+)?")),
         PANEL_COUNT = as.numeric(str_extract(INFO, "(?<=PANEL_COUNT=)\\\d+"))
       ) %>% separate(X10, into=c(NA,"CN"))

       passing_events <- vcf_df %>% filter(ALT != ".", FILTER=="PASS")
       n_passing_events <- nrow(passing_events)

       contigs <- passing_events %>% pull(CONTIG)
       starts <- passing_events %>% pull(START)
       ends <- passing_events %>% pull(END)
       alts <- passing_events %>% pull(ALT)
       cns <- passing_events %>% pull(CN)
       quals <- passing_events %>% pull(QUAL)

       read_counts <- function(file) {
         h5f <- H5Fopen(file, flags = "H5F_ACC_RDONLY")
         contigs <- h5read(h5f, 'intervals/indexed_contig_names') %>% as_tibble()
         colnames(contigs) <- "CONTIG"
         contigs <- contigs %>% mutate(contig_idx = seq(0, nrow(contigs)-1))
         counts <- h5read(h5f, 'intervals/transposed_index_start_end') %>% as_tibble(.name_repair = "minimal")
         colnames(counts) <- c("contig_idx","START","END")
         counts <- counts %>% inner_join(contigs)
         count_values <- h5read(h5f, 'counts/values') %>% as_tibble()
         colnames(count_values) <- "counts"
         counts <- counts %>% bind_cols(count_values)
         sample_names <- h5read(h5f, 'sample_metadata/sample_name') %>% as_tibble()
         colnames(sample_names) <- "sample_name"
         counts <- counts %>% bind_cols(sample_names)
         mean_counts <- counts %>% pull(counts) %>% mean()
         counts <- counts %>% mutate(norm_counts = counts/mean_counts)
         H5Fclose(h5f)
         return(counts %>% select(-contig_idx))
       }

       panel_read_count_hdf5s <- c("~{sep='","' panel_read_counts}")
       panel_read_counts <- panel_read_count_hdf5s %>% map(read_counts) %>% reduce(bind_rows)
       H5garbage_collect()
       panel_mean_counts <- panel_read_counts %>% group_by(CONTIG, START,END) %>% summarise(panel_mean_norm_counts = mean(norm_counts))
       panel_read_counts <- panel_read_counts %>% inner_join(panel_mean_counts)
       panel_read_counts <- panel_read_counts %>% mutate(adjusted_counts = 2*norm_counts/panel_mean_norm_counts) %>% drop_na()
       panel_read_counts <- panel_read_counts %>% semi_join(cr_case)

       case_read_counts <- read_counts("~{case_read_counts}")
       case_read_counts <- case_read_counts %>% inner_join(panel_mean_counts)
       case_read_counts <- case_read_counts %>% mutate(adjusted_counts = 2*norm_counts/panel_mean_norm_counts) %>% drop_na()
       case_read_counts <- case_read_counts %>% semi_join(cr_case)
       intervals <- read_tsv("~{interval_list}",
                             comment="@")
       case_read_counts <- case_read_counts %>% left_join(intervals)

       \`\`\`

       # {.tabset .tabset-pills}

       ## Overview {.tabset .custom-tabset}

       ### By Chromosome
       \`\`\`{r overview-chr, fig.align='center', echo=FALSE, message=FALSE, warning=FALSE}
       ggplot(cr_case, aes(x=START, y=LINEAR_COPY_RATIO)) +
         geom_point(alpha=0.2, size=0.2) +
         facet_wrap(~CONTIG, scales="free_x") +
         theme_bw() + ylim(0,5) +
         theme( axis.text.x=element_blank()) +
         xlab("Position") + ylab("Denoised Linear Copy Ratio")

       \`\`\`

       ### By GC Content
       \`\`\`{r overview-gc, fig.align='center', echo=FALSE, message=FALSE, warning=FALSE}
       ggplot(case_read_counts, aes(x=GC_CONTENT, y=adjusted_counts)) +
         geom_point(size=0.2, alpha=0.2) +
         xlim(0,1) + theme_bw() + ylim(0,5) +
         xlab("GC Content") +
         ylab("Adjusted Counts")
       \`\`\`

       \`\`\`{r events, results='asis', fig.align='center', echo=FALSE, message=FALSE, warning=FALSE}
       if (n_passing_events > 0) {
           for (i in 1:n_passing_events) {
             contig<- contigs[[i]]
             start<- starts[[i]]
             end<- ends[[i]]
             alt <- alts[[i]]
             cn <- cns[[i]]
             qual <- quals[[i]]
             width <- end - start

             cat(paste0("## ", contig,":",comma(start),"-",comma(end), " {.tabset .custom-tabset}"))


             points_before <- cr_case %>% filter(CONTIG==contig, START<start) %>% pull(START)
             points_after <- cr_case %>% filter(CONTIG==contig, START>end) %>% pull(START)
             total_before <- length(points_before)
             total_after <- length(points_after)

             n_before <- 0
             n_after <- 0
             min_non_event_points <- 10
             max_width_factor <- 10
             min_width_factor <- 2.5
             get_point_before <- function(n_step) {
               if (n_step>0) {
                 return(points_before[[total_before-n_step+1]])
               } else {
                return(start)
               }
             }

             get_point_after <- function(n_step) {
               if (n_step>0) {
                 return(points_after[[n_step]])
               } else {
                 return(end)
               }
             }

             while ((n_before < total_before ||
                    n_after < total_after) &&
                    (get_point_after(n_after) - get_point_before(n_before) < min_width_factor * width ||
                     n_after + n_before < min_non_event_points) &&
                    get_point_after(n_after) - get_point_before(n_before) < max_width_factor * width) {
               prev_n_before <- n_before
               prev_n_after <- n_after
               #decide which direction to increase
               # if already to one edge, increase the other
               if (n_before == total_before) {
                 n_after <- n_after +1
               } else if (n_after == total_after) {
                 n_before <- n_before + 1
               } else {
                 # increment whichever direction is a smaller step
                 step_before <- get_point_before(n_before) - get_point_before(n_before+1)
                 step_after <- get_point_after(n_after+1) - get_point_after(n_after)
                 if (step_before < step_after) {
                   n_before <- n_before + 1
                 } else {
                   n_after <- n_after + 1
                 }
               }
             }
             if (get_point_after(n_after) - get_point_before(n_before) > max_width_factor * width) {
               n_after <- prev_n_after
               n_before <- prev_n_before
             }

             window_start <- get_point_before(n_before)
             window_end <- get_point_after(n_after)

             sorted_sample_names <- cr_panel %>% filter(CONTIG==contig, START>=start, START<=end) %>% group_by(sample_name) %>%
                                summarise(mean_cr = mean(LINEAR_COPY_RATIO)) %>% arrange(mean_cr) %>% pull(sample_name)
             extreme_samples <- if (alt=="<DEL>") {
               sorted_sample_names[1:2]
             } else {
               sorted_sample_names[(length(sorted_sample_names)-1):length(sorted_sample_names)]
             }

             cat("\n")
             cat("### Denoised Copy Ratio")
             cat("\n")

             p_cr <- ggplot(cr_panel %>% filter(CONTIG==contig, START>=window_start, START<=window_end), aes(x=START, y=LINEAR_COPY_RATIO)) +
               geom_point(alpha=0.2, aes(color='Panel')) + theme_bw() + ylim(0,7) +
               geom_point(data=cr_case %>% filter(CONTIG==contig, START>=window_start, START<=window_end), aes(color="Case")) +
               geom_line(data=cr_case %>% filter(CONTIG==contig, START>=window_start, START<=window_end),aes(color="Case")) +
               annotate("rect",xmin=start, xmax=end, ymin=0, ymax=Inf, alpha=0.2) +
               guides(color = guide_legend(override.aes = list(alpha = 1))) +
               scale_color_manual(name="",values=c("Panel"="blue",
                                                   "Case"="black"))+
               geom_hline(yintercept=1, linetype="dashed", alpha=0.8) +
               geom_hline(yintercept=2, linetype="dashed", alpha=0.8) +
               geom_hline(yintercept=3, linetype="dashed", alpha=0.8) +
               scale_x_continuous(labels=comma) +
               xlab("Position") + ylab("Denoised Linear Copy Ratio")
             print(p_cr)

             cat("\n\n")
             cat("### Adjusted Read Counts")
             cat("\n")

             p_counts <- ggplot(panel_read_counts %>% filter(CONTIG==contig, START>=window_start, START<=window_end), aes(x=START, y=adjusted_counts)) +
               geom_point(alpha=0.2, aes(color='Panel')) + theme_bw() + ylim(0,7) +
               geom_point(data=case_read_counts %>% filter(CONTIG==contig, START>=window_start, START<=window_end), aes(color="Case")) +
               geom_line(data=case_read_counts %>% filter(CONTIG==contig, START>=window_start, START<=window_end),aes(color="Case")) +
               annotate("rect",xmin=start, xmax=end, ymin=0, ymax=Inf, alpha=0.2) +
               guides(color = guide_legend(override.aes = list(alpha = 1))) +
               scale_color_manual(name="",values=c("Panel"="blue",
                                                   "Panel Second Most Extreme Values"="blue",
                                                   "Case"="black"))+
               geom_hline(yintercept=1, linetype="dashed", alpha=0.8) +
               geom_hline(yintercept=2, linetype="dashed", alpha=0.8) +
               geom_hline(yintercept=3, linetype="dashed", alpha=0.8) +
               scale_x_continuous(labels=comma) +
               xlab("Position") + ylab("Adjusted Read Counts")
             print(p_counts)

             t <- tibble(x=c(1,2), y=c(1,2))

             cat("\n\n")
             cat("\n\n")
           }
       }
       \`\`\`
       EOF

       Rscript -e "library(rmarkdown); rmarkdown::render('~{output_prefix}_cnv_event_report.Rmd', 'html_document')"
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/r-gcnv-viz@sha256:c92a9a26cba4ab10b579776cd4dd70d5fca485d4a7216d0ab9b477662e3d9228"
        disks: "local-disk 100 HDD"
        memory: mem_gb + " GB"
      }

      output {
        File cnv_event_report = "~{output_prefix}_cnv_event_report.html"
      }
}