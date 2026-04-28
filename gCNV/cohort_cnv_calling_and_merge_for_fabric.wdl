version 1.0

import "cnv_germline_cohort_workflow.wdl" as cnv_germline_cohort
import "single_sample_cnv_germline_case_filter_workflow.wdl" as cnv_germline_case_filter
import "cnv_calling_and_merge_for_fabric.wdl" as cnv_calling_and_merge_for_fabric
workflow CohortCNVCallingAndMergeForFabric {
    input {
        ##################################
      #### required basic arguments ####
      ##################################
      File intervals
      Array[File]+ normal_bams
      Array[File]+ normal_bais
      String cohort_entity_id
      File contig_ploidy_priors
      Int num_intervals_per_scatter
      File ref_fasta_dict
      File ref_fasta_fai
      File ref_fasta
      String gatk_docker
      Array[File] short_variant_vcfs

      ####################################################
      #### optional arguments for PreprocessIntervals ####
      ####################################################
      Int? padding
      Int? bin_length

      Int? mem_gb_for_germline_cnv_caller
      Int? cpu_for_germline_cnv_caller

      Int? mem_gb_for_postprocess_germline_cnv_calls
      Int? disk_space_gb_for_postprocess_germline_cnv_calls

      Int ref_copy_number_autosomal_contigs
      

      ##########################
      #### arguments for QC ####
      ##########################
      Int maximum_number_events_per_sample
      Int maximum_number_pass_events_per_sample


      ##################################################################
      #### optional arguments for ExtractPoNFreqAnnotateFilterAndQC ####
      ##################################################################
      Array[String] filter_expressions = ['(GT=="alt" | GT=="mis") & ((FMT/CN>1 & QUAL<50) | (FMT/CN==1 & QUAL<100 ) | (FMT/CN==0 & QUAL<400))','(GT=="alt" | GT=="mis") & (INFO/PANEL_FREQ>0.02)']
      Array[String] filter_names = ['LowQual','PanelOverlap']
      Int? mem_gb_for_extract_pon_freq
      Int? disk_for_extract_pon_freq
      Float? overlap_thresh
    }

    call cnv_germline_cohort.CNVGermlineCohortWorkflow {
        input:
            intervals = intervals,
            normal_bams = normal_bams,
            normal_bais = normal_bais,
            cohort_entity_id = cohort_entity_id,
            contig_ploidy_priors = contig_ploidy_priors,
            num_intervals_per_scatter = num_intervals_per_scatter,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta = ref_fasta,
            gatk_docker = gatk_docker,
            padding = padding,
            bin_length = bin_length,
            mem_gb_for_postprocess_germline_cnv_calls = mem_gb_for_postprocess_germline_cnv_calls,
            disk_space_gb_for_postprocess_germline_cnv_calls = disk_space_gb_for_postprocess_germline_cnv_calls,
            maximum_number_events_per_sample = maximum_number_events_per_sample,
            maximum_number_pass_events_per_sample = maximum_number_pass_events_per_sample,
            mem_gb_for_germline_cnv_caller = mem_gb_for_germline_cnv_caller,
            cpu_for_germline_cnv_caller = cpu_for_germline_cnv_caller,
            ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs
    }

    scatter (i in range(length(CNVGermlineCohortWorkflow.genotyped_segments_vcfs))) {
        call cnv_germline_case_filter.ExtractPoNFreqAnnotateFilterAndQC {
            input:
                vcf = CNVGermlineCohortWorkflow.genotyped_segments_vcfs[i],
                vcf_idx = CNVGermlineCohortWorkflow.genotyped_segments_vcf_indexes[i],
                panel_vcfs = CNVGermlineCohortWorkflow.genotyped_segments_vcfs, #!NonemptyCoercion
                intervals = CNVGermlineCohortWorkflow.preprocessed_intervals,
                overlap_thresh = overlap_thresh,
                filter_expressions = filter_expressions,
                filter_names = filter_names,
                max_events = maximum_number_events_per_sample,
                max_pass_events = maximum_number_pass_events_per_sample,
                sample_name = CNVGermlineCohortWorkflow.read_counts_entity_ids[i],
                mem_gb = mem_gb_for_extract_pon_freq,
                disk_size_gb = disk_for_extract_pon_freq
        }

        Boolean qc_passed_scatter = (ExtractPoNFreqAnnotateFilterAndQC.qc_status == "PASS")

        call cnv_calling_and_merge_for_fabric.ReformatAndMergeForFabric {
            input:
                cnv_vcf = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf,
                short_variant_vcf = short_variant_vcfs[i],
                gatk_docker = gatk_docker
        }

        call cnv_calling_and_merge_for_fabric.GCNVVisualzation {
            input:
            filtered_vcf = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf,
            case_copy_ratios = CNVGermlineCohortWorkflow.denoised_copy_ratios[i],
            case_read_counts = CNVGermlineCohortWorkflow.read_counts[i],
            panel_copy_ratios = CNVGermlineCohortWorkflow.denoised_copy_ratios, #!NonemptyCoercion
            panel_read_counts = CNVGermlineCohortWorkflow.read_counts, #!NonemptyCoercion
            interval_lists = CNVGermlineCohortWorkflow.sharded_interval_lists #!NonemptyCoercion
        }

        call LowGCDropoutQC {
            input:
                counts_hdf5 = CNVGermlineCohortWorkflow.read_counts[i],
                annotated_intervals = select_first([CNVGermlineCohortWorkflow.annotated_intervals]),
        }
    }


    output {
        Array[File] filtered_cnv_genotyped_segments_vcf = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf
        Array[File] filtered_cnv_genotyped_segments_vcf_index = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf_index
        Array[File] filtered_cnv_genotyped_segments_vcf_md5sum = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf_md5sum

        Array[File] merged_vcf = ReformatAndMergeForFabric.merged_vcf
        Array[File] merged_vcf_index = ReformatAndMergeForFabric.merged_vcf_index
        Array[File] merged_vcf_md5sum = ReformatAndMergeForFabric.merged_vcf_md5sum

        Array[Boolean] qc_passed = qc_passed_scatter
        Array[File] cnv_metrics = ExtractPoNFreqAnnotateFilterAndQC.cnv_metrics
        Array[File] cnv_event_report = GCNVVisualzation.cnv_event_report
        Array[File] low_gc_dropout_metric = LowGCDropoutQC.low_gc_drouput_tsv

    }

}

task LowGCDropoutQC {
    input {
        File counts_hdf5
        File annotated_intervals
        Int mem_gb = 4
    }

    String basename = basename(counts_hdf5, ".hdf5")
    command <<<
        set -euo pipefail

        python3 << "EOF"
        import h5py
        import pandas as pd
        import numpy as np

        def read_counts(file_path):
            with h5py.File(file_path, 'r') as f:
                # 1. Read contig names and create index
                # .asstr() handles the common HDF5 byte-string issue in Python
                contigs = pd.DataFrame({
                    'CONTIG': f['intervals/indexed_contig_names'][:].astype(str)
                })
                contigs['contig_idx'] = np.arange(len(contigs))

                # 2. Read interval coordinates (START/END)
                # HDF5 'transposed_index_start_end' is typically a 2D array
                coords = f['intervals/transposed_index_start_end'][:].T
                counts_df = pd.DataFrame(coords, columns=['contig_idx', 'START', 'END'])

                # 3. Join contig names to the coordinates
                counts_df = counts_df.merge(contigs, on='contig_idx')

                # 4. Read the raw count values
                counts_df['counts'] = f['counts/values'][:].flatten()

                # 5. Read Sample Name
                # In GATK HDF5, this is usually a single value or one per bin
                sample_name = f['sample_metadata/sample_name'][0].decode('utf-8')
                counts_df['sample_name'] = sample_name

                # 6. Normalize by the mean
                mean_counts = counts_df['counts'].mean()
                counts_df['norm_counts'] = counts_df['counts'] / mean_counts

                # 7. Final cleanup (remove the internal index)
                return counts_df.drop(columns=['contig_idx']), sample_name

        counts_df, sample_name = read_counts("~{counts_hdf5}")
        intervals = pd.read_csv("~{annotated_intervals}", sep='\t', comment="@")
        counts_df = counts_df.merge(intervals)

        low_gc_bin_counts = counts_df[(counts_df['GC_CONTENT'] > 0.25) && (counts_df['GC_CONTENT'] < 0.3)]['counts']
        low_gc_dropout_frac = (low_gc_bin_counts < 0.25 * low_gc_bin_counts.mean()).mean()


        with open(f'~{basename}_low_gc_dropout_qc.tsv', 'w') as f:
            f.write(f"sample\tlow_gc_dropout_frac\n")
            f.write(f"{sample_name}\t{low_gc_dropout_frac:.4f}\n")        
        
        EOF


    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim"
        preemptible: 3
        cpu: 1
        disks: "local-disk 50 HDD"
        memory: mem_gb + " GB"
    }

    output {
        File low_gc_drouput_tsv = "~{basename}_low_gc_dropout_qc.tsv"
    }
}