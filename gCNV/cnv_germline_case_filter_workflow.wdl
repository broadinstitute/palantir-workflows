version 1.0

import "cnv_germline_case_workflow.wdl" as GermlineCNVCaseWorkflow

workflow SingleSampleGCNVAndFilterVCFs {
    input {
        Array[String]? filter_expressions #= "QUAL < 50"  # Example filter criteria
        Array[String]? filter_names

        ##################################
        ##for case mode
        #### required basic arguments ####
        ##################################
        File intervals
        File? blacklist_intervals
        File filtered_intervals
        String normal_bam
        String normal_bai
        File contig_ploidy_model_tar
        Array[File]+ gcnv_model_tars
        Array[File]+ pon_genotyped_segments_vcfs
        Int num_intervals_per_scatter
        File ref_fasta_dict
        File ref_fasta_fai
        File ref_fasta
        String gatk_docker
        Float overlap_thresh = 0.5

        ##################################
        #### optional basic arguments ####
        ##################################
        File? gatk4_jar_override
        Int? preemptible_attempts

        # Required if BAM/CRAM is in a requester pays bucket
        String? gcs_project_for_requester_pays

        ####################################################
        #### optional arguments for PreprocessIntervals ####
        ####################################################
        Int? padding
        Int? bin_length

        ##############################################
        #### optional arguments for CollectCounts ####
        ##############################################
        Array[String]? disabled_read_filters_for_collect_counts
        String? collect_counts_format
        Boolean? collect_counts_enable_indexing
        Int? mem_gb_for_collect_counts

        ######################################################################
        #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
        ######################################################################
        Float? ploidy_mapping_error_rate
        Float? ploidy_sample_psi_scale
        Int? mem_gb_for_determine_germline_contig_ploidy
        Int? cpu_for_determine_germline_contig_ploidy
        Int? disk_for_determine_germline_contig_ploidy

        ##########################################################
        #### optional arguments for GermlineCNVCallerCaseMode ####
        ##########################################################
        Float? gcnv_p_alt
        Float? gcnv_cnv_coherence_length
        Int? gcnv_max_copy_number
        Int? mem_gb_for_germline_cnv_caller
        Int? cpu_for_germline_cnv_caller
        Int? disk_for_germline_cnv_caller

        # optional arguments for germline CNV denoising model
        Float? gcnv_mapping_error_rate
        Float? gcnv_sample_psi_scale
        Float? gcnv_depth_correction_tau
        String? gcnv_copy_number_posterior_expectation_mode
        Int? gcnv_active_class_padding_hybrid_mode

        # optional arguments for Hybrid ADVI
        Float? gcnv_learning_rate
        Float? gcnv_adamax_beta_1
        Float? gcnv_adamax_beta_2
        Int? gcnv_log_emission_samples_per_round
        Float? gcnv_log_emission_sampling_median_rel_error
        Int? gcnv_log_emission_sampling_rounds
        Int? gcnv_max_advi_iter_first_epoch
        Int? gcnv_max_advi_iter_subsequent_epochs
        Int? gcnv_min_training_epochs
        Int? gcnv_max_training_epochs
        Float? gcnv_initial_temperature
        Int? gcnv_num_thermal_advi_iters
        Int? gcnv_convergence_snr_averaging_window
        Float? gcnv_convergence_snr_trigger_threshold
        Int? gcnv_convergence_snr_countdown_window
        Int? gcnv_max_calling_iters
        Float? gcnv_caller_update_convergence_threshold
        Float? gcnv_caller_internal_admixing_rate
        Float? gcnv_caller_external_admixing_rate
        Boolean? gcnv_disable_annealing

        ###################################################
        #### arguments for PostprocessGermlineCNVCalls ####
        ###################################################
        Int ref_copy_number_autosomal_contigs
        Array[String]? allosomal_contigs

        ##########################
        #### arguments for QC ####
        ##########################
        Int maximum_number_events_per_sample
        Int maximum_number_pass_events_per_sample

    }


    call GermlineCNVCaseWorkflow.CNVGermlineCaseWorkflow as CNVGermlineCaseWorkflow {
        input:
            intervals = intervals,
            blacklist_intervals=blacklist_intervals,
            filtered_intervals=filtered_intervals,
            normal_bams=[normal_bam],
            normal_bais=[normal_bai],
            contig_ploidy_model_tar=contig_ploidy_model_tar,
            gcnv_model_tars=gcnv_model_tars,
            num_intervals_per_scatter=num_intervals_per_scatter,
            ref_fasta_dict=ref_fasta_dict,
            ref_fasta_fai=ref_fasta_fai,
            ref_fasta=ref_fasta,
            gatk_docker=gatk_docker,

            ##################################
            #### optional basic arguments ####
            ##################################
            gatk4_jar_override=gatk4_jar_override,
            preemptible_attempts=preemptible_attempts,

            # Required if BAM/CRAM is in a requester pays bucket
            gcs_project_for_requester_pays=gcs_project_for_requester_pays,

            ####################################################
            #### optional arguments for PreprocessIntervals ####
            ####################################################
            padding=padding,
            bin_length=bin_length,

            ##############################################
            #### optional arguments for CollectCounts ####
            ##############################################
            disabled_read_filters_for_collect_counts=disabled_read_filters_for_collect_counts,
            collect_counts_format=collect_counts_format,
            collect_counts_enable_indexing=collect_counts_enable_indexing,
            mem_gb_for_collect_counts=mem_gb_for_collect_counts,

            ######################################################################
            #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
            ######################################################################
            ploidy_mapping_error_rate=ploidy_mapping_error_rate,
            ploidy_sample_psi_scale=ploidy_sample_psi_scale,
            mem_gb_for_determine_germline_contig_ploidy=mem_gb_for_determine_germline_contig_ploidy,
            cpu_for_determine_germline_contig_ploidy=cpu_for_determine_germline_contig_ploidy,
            disk_for_determine_germline_contig_ploidy=disk_for_determine_germline_contig_ploidy,

            ##########################################################
            #### optional arguments for GermlineCNVCallerCaseMode ####
            ##########################################################
            gcnv_p_alt=gcnv_p_alt,
            gcnv_cnv_coherence_length=gcnv_cnv_coherence_length,
            gcnv_max_copy_number=gcnv_max_copy_number,
            mem_gb_for_germline_cnv_caller=mem_gb_for_germline_cnv_caller,
            cpu_for_germline_cnv_caller=cpu_for_germline_cnv_caller,
            disk_for_germline_cnv_caller=disk_for_germline_cnv_caller,

            # optional arguments for germline CNV denoising model
            gcnv_mapping_error_rate=gcnv_mapping_error_rate,
            gcnv_sample_psi_scale=gcnv_sample_psi_scale,
            gcnv_depth_correction_tau=gcnv_depth_correction_tau,
            gcnv_copy_number_posterior_expectation_mode=gcnv_copy_number_posterior_expectation_mode,
            gcnv_active_class_padding_hybrid_mode=gcnv_active_class_padding_hybrid_mode,

            # optional arguments for Hybrid ADVI
            gcnv_learning_rate=gcnv_learning_rate,
            gcnv_adamax_beta_1=gcnv_adamax_beta_1,
            gcnv_adamax_beta_2=gcnv_adamax_beta_2,
            gcnv_log_emission_samples_per_round=gcnv_log_emission_samples_per_round,
            gcnv_log_emission_sampling_median_rel_error=gcnv_log_emission_sampling_median_rel_error,
            gcnv_log_emission_sampling_rounds=gcnv_log_emission_sampling_rounds,
            gcnv_max_advi_iter_first_epoch=gcnv_max_advi_iter_first_epoch,
            gcnv_max_advi_iter_subsequent_epochs=gcnv_max_advi_iter_subsequent_epochs,
            gcnv_min_training_epochs=gcnv_min_training_epochs,
            gcnv_max_training_epochs=gcnv_max_training_epochs,
            gcnv_initial_temperature=gcnv_initial_temperature,
            gcnv_num_thermal_advi_iters=gcnv_num_thermal_advi_iters,
            gcnv_convergence_snr_averaging_window=gcnv_convergence_snr_averaging_window,
            gcnv_convergence_snr_trigger_threshold=gcnv_convergence_snr_trigger_threshold,
            gcnv_convergence_snr_countdown_window=gcnv_convergence_snr_countdown_window,
            gcnv_max_calling_iters=gcnv_max_calling_iters,
            gcnv_caller_update_convergence_threshold=gcnv_caller_update_convergence_threshold,
            gcnv_caller_internal_admixing_rate=gcnv_caller_internal_admixing_rate,
            gcnv_caller_external_admixing_rate=gcnv_caller_external_admixing_rate,
            gcnv_disable_annealing=gcnv_disable_annealing,

            ###################################################
            #### arguments for PostprocessGermlineCNVCalls ####
            ###################################################
            ref_copy_number_autosomal_contigs=ref_copy_number_autosomal_contigs,
            allosomal_contigs=allosomal_contigs,

            ##########################
            #### arguments for QC ####
            ##########################
            maximum_number_events_per_sample=maximum_number_events_per_sample,
            maximum_number_pass_events_per_sample=maximum_number_pass_events_per_sample
    }

    File vcf = CNVGermlineCaseWorkflow.genotyped_segments_vcfs[0]
    File vcf_index = CNVGermlineCaseWorkflow.genotyped_intervals_vcf_indexes[0]
    String samplename = CNVGermlineCaseWorkflow.entity_id[0]

    call ExtractPoNFreqAnnotateFilterAndQC {
        input:
            vcf = vcf,
            vcf_idx = vcf_index,
            panel_vcfs = pon_genotyped_segments_vcfs,
            intervals = intervals,
            overlap_thresh = overlap_thresh,
            filter_expressions = filter_expressions,
            filter_names = filter_names,
            max_events = maximum_number_events_per_sample,
            max_pass_events = maximum_number_pass_events_per_sample
    }

    output {
        File filtered_vcf = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf
        File filtered_vcf_index = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf_index
        File filtered_vcf_md5sum = ExtractPoNFreqAnnotateFilterAndQC.filtered_vcf_md5sum

        File preprocessed_intervals = CNVGermlineCaseWorkflow.preprocessed_intervals
        File read_counts_entity_id = CNVGermlineCaseWorkflow.read_counts_entity_id[0]
        File read_counts = CNVGermlineCaseWorkflow.read_counts[0]
        File sample_contig_ploidy_calls_tars = CNVGermlineCaseWorkflow.sample_contig_ploidy_calls_tars[0]
        Array[File] gcnv_calls_tars = CNVGermlineCaseWorkflow.gcnv_calls_tars[0]
        File gcnv_tracking_tars = CNVGermlineCaseWorkflow.gcnv_tracking_tars[0]
        File genotyped_intervals_vcfs = CNVGermlineCaseWorkflow.genotyped_intervals_vcfs[0]
        File genotyped_intervals_vcf_indexes = CNVGermlineCaseWorkflow.genotyped_intervals_vcf_indexes[0]
        File genotyped_segments_vcfs = CNVGermlineCaseWorkflow.genotyped_segments_vcfs[0]
        File genotyped_segments_vcf_indexes = CNVGermlineCaseWorkflow.genotyped_segments_vcf_indexes[0]
        String qc_status_string = ExtractPoNFreqAnnotateFilterAndQC.qc_status
        Boolean qc_passed = (ExtractPoNFreqAnnotateFilterAndQC.qc_status == "PASS")
        File cnv_metrics = ExtractPoNFreqAnnotateFilterAndQC.cnv_metrics
        File denoised_copy_ratios = CNVGermlineCaseWorkflow.denoised_copy_ratios[0]
    }
}


task ExtractSamplename {
    input {
        String vcf_filename
        String gatk_docker
    }
    command <<<
        set -euo pipefail
        echo ~{vcf_filename} | sed -E 's/genotyped-segments-(.*)\.cram\.vcf\.gz/\1/' > output.txt
        >>>
    output {
        String samplename = read_string("output.txt")
    }
    runtime {
        docker: gatk_docker
        memory: "2G"
        cpu: 2
    }
}

task ExtractPoNFreq {
    input {
        File vcf
        Array[File] panel_vcfs
        File intervals
        Float overlap_thresh = 0.5
        Int mem_gb = 16
        Int disk_size_gb = 100
    }

    String basename_out = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        python << "EOF"
        import pandas as pd
        import numpy as np
        import gzip

        intervals = pd.read_csv("~{intervals}", sep="\t", comment="@", names = ["contig","start","end","dummy1","dummy2"])
        intervals['contig_idx'] = intervals.groupby('contig').cumcount()
        intervals = intervals.set_index(intervals.contig + "_" + intervals.contig_idx.astype(str))

        def add_exon_idxs(df, exons):
                contigs = set(df.contig)
                for contig in contigs:
                    df.loc[df.contig==contig,"start_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].end,
                                                                               df.loc[df.contig==contig].start,"left")
                    df.loc[df.contig==contig,"end_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].start,
                                                                               df.loc[df.contig==contig].end,"right")

        def get_exon_expanded_events(df, exons):
            add_exon_idxs(df, exons)
            df = df.loc[df.start_exon_idx != df.end_exon_idx].reset_index().astype({'start_exon_idx':int,'end_exon_idx':int})
            df_expanded = df.loc[df.index.repeat(df.end_exon_idx-df.start_exon_idx)]
            df_expanded['exon_idx'] = df_expanded.groupby(df_expanded.index).cumcount() + df_expanded.start_exon_idx
            df_expanded = df_expanded.set_index(df_expanded.contig + "_" + df_expanded.exon_idx.astype(str))
            df_expanded = df_expanded.join(exons[['start','end']], rsuffix='_exon')
            df_expanded['event_exon_start']=np.maximum(df_expanded.start, df_expanded.start_exon)
            df_expanded['event_exon_end']=np.minimum(df_expanded.end, df_expanded.end_exon)
            df_expanded = df_expanded.set_index(df_expanded.index + "_" + df_expanded.svtype)
            return df_expanded

        def standardize_gt_vcf(df):
            df['end'] = df['INFO'].str.replace('END=','').astype(int) #split('=').apply(lambda x: x[1])
            df["svtype"] = df.ALT.str.replace("<","").str.replace(">","")
            for num,f in enumerate(df.loc[0,'FORMAT'].split(':')):
                df[f] = df['SAMPLE'].str.split(':').apply(lambda x: x[num])
            return df

        def get_sample_id(vcf):
            with gzip.open(vcf,"r") as f_vcf:
                for line in f_vcf:
                    if line.startswith(b'#CHROM'):
                        return line.split(b'\t')[-1].decode()

        def read_vcf_to_df(vcf):
            df = pd.read_csv(vcf, sep='\t',comment='#',compression='gzip',
                                 names=['contig','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'])
            df = standardize_gt_vcf(df)
            df['sample_name']=get_sample_id(vcf)
            for c in ['NP','CN','QA','QS']:
                df[c]=df[c].astype(int)
            return df
            
        df = read_vcf_to_df("~{vcf}")
        df_expanded = get_exon_expanded_events(df, intervals)

        vcfs = ["~{sep='","' panel_vcfs}"]
        df_panel = pd.concat([read_vcf_to_df(vcf) for vcf in vcfs])
        df_panel_expanded = get_exon_expanded_events(df_panel, intervals)

        df_expanded_with_panel = df_expanded.join(df_panel_expanded, how="left", lsuffix="_sample", rsuffix="_panel")

        df_expanded_with_panel['overlapping_panel_exon_start']=np.maximum(df_expanded_with_panel.event_exon_start_sample, df_expanded_with_panel.event_exon_start_panel)
        df_expanded_with_panel['overlapping_panel_exon_end']=np.minimum(df_expanded_with_panel.event_exon_end_sample, df_expanded_with_panel.event_exon_end_panel)

        df_expanded['event_exon_length']=df_expanded.event_exon_end - df_expanded.event_exon_start
        df_expanded_with_panel['overlapping_panel_exon_length']=np.maximum(df_expanded_with_panel.overlapping_panel_exon_end - df_expanded_with_panel.overlapping_panel_exon_start,0).fillna(0)

        df_panel_counts = df_expanded_with_panel.groupby(['ID_sample','sample_name_panel']).agg(
            {'overlapping_panel_exon_length':'sum'}
            ).reset_index('sample_name_panel').join(df_expanded.groupby('ID').agg({'event_exon_length':'sum'})
        )

        df_panel_counts['PANEL_COUNT'] = np.where(df_panel_counts.overlapping_panel_exon_length/df_panel_counts.event_exon_length>~{overlap_thresh}, 1, 0)
        df_panel_counts = df_panel_counts.groupby(df_panel_counts.index).agg({'PANEL_COUNT':'sum'})

        df = df.set_index('ID').join(df_panel_counts).fillna({'PANEL_COUNT':0})

        n_panel_samples = len(set(df_panel.sample_name))
        df['PANEL_FREQ'] = df['PANEL_COUNT']/n_panel_samples

        df_annotations = df[['contig','start','PANEL_FREQ','PANEL_COUNT']].copy()
        df_annotations = df_annotations.rename({"contig":"CHROM","start":"POS"}, axis=1)
        df_annotations = df_annotations.astype({"PANEL_COUNT":int})
        df_annotations.to_csv("~{basename_out}.annotations.tsv", index = False, sep="\t", float_format="%g")


        EOF

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim"
        preemptible: 3
        cpu: 2
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GB"
    }

    output {
        File annotations = "~{basename_out}.annotations.tsv"
    }
}

task ExtractPoNFreqAnnotateFilterAndQC {
    input {
        File vcf
        File vcf_idx
        Array[File] panel_vcfs
        File intervals
        Float overlap_thresh = 0.5
        Array[String] filter_expressions = ['(GT=="alt" | GT=="mis") & ((FMT/CN>1 & QUAL<50) | (FMT/CN==1 & QUAL<100 ) | (FMT/CN==0 & QUAL<400))','(GT=="alt" | GT=="mis") & (INFO/PANEL_COUNT>1)']
        Array[String] filter_names = ['LowQual','PanelOverlap']
        Int max_events
        Int max_pass_events
        Int mem_gb = 16
        Int disk_size_gb = 100
    }

    String output_basename = basename(vcf)
    command <<<
        set -euo pipefail

        python << "EOF"
        import pandas as pd
        import numpy as np
        import gzip

        intervals = pd.read_csv("~{intervals}", sep="\t", comment="@", names = ["contig","start","end","dummy1","dummy2"])
        intervals['contig_idx'] = intervals.groupby('contig').cumcount()
        intervals = intervals.set_index(intervals.contig + "_" + intervals.contig_idx.astype(str))

        def add_exon_idxs(df, exons):
                contigs = set(df.contig)
                for contig in contigs:
                    df.loc[df.contig==contig,"start_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].end,
                                                                               df.loc[df.contig==contig].start,"left")
                    df.loc[df.contig==contig,"end_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].start,
                                                                               df.loc[df.contig==contig].end,"right")

        def get_exon_expanded_events(df, exons):
            add_exon_idxs(df, exons)
            df = df.loc[df.start_exon_idx != df.end_exon_idx].reset_index().astype({'start_exon_idx':int,'end_exon_idx':int})
            df_expanded = df.loc[df.index.repeat(df.end_exon_idx-df.start_exon_idx)]
            df_expanded['exon_idx'] = df_expanded.groupby(df_expanded.index).cumcount() + df_expanded.start_exon_idx
            df_expanded = df_expanded.set_index(df_expanded.contig + "_" + df_expanded.exon_idx.astype(str))
            df_expanded = df_expanded.join(exons[['start','end']], rsuffix='_exon')
            df_expanded['event_exon_start']=np.maximum(df_expanded.start, df_expanded.start_exon)
            df_expanded['event_exon_end']=np.minimum(df_expanded.end, df_expanded.end_exon)
            df_expanded = df_expanded.set_index(df_expanded.index + "_" + df_expanded.svtype)
            return df_expanded

        def standardize_gt_vcf(df):
            df['end'] = df['INFO'].str.replace('END=','').astype(int) #split('=').apply(lambda x: x[1])
            df["svtype"] = df.ALT.str.replace("<","").str.replace(">","")
            for num,f in enumerate(df.loc[0,'FORMAT'].split(':')):
                df[f] = df['SAMPLE'].str.split(':').apply(lambda x: x[num])
            return df

        def get_sample_id(vcf):
            with gzip.open(vcf,"r") as f_vcf:
                for line in f_vcf:
                    if line.startswith(b'#CHROM'):
                        return line.split(b'\t')[-1].decode()

        def read_vcf_to_df(vcf):
            df = pd.read_csv(vcf, sep='\t',comment='#',compression='gzip',
                                 names=['contig','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'])
            df = standardize_gt_vcf(df)
            df['sample_name']=get_sample_id(vcf)
            for c in ['NP','CN','QA','QS']:
                df[c]=df[c].astype(int)
            return df

        df = read_vcf_to_df("~{vcf}")
        df_expanded = get_exon_expanded_events(df, intervals)

        vcfs = ["~{sep='","' panel_vcfs}"]
        df_panel = pd.concat([read_vcf_to_df(vcf) for vcf in vcfs])
        df_panel_expanded = get_exon_expanded_events(df_panel, intervals)

        df_expanded_with_panel = df_expanded.join(df_panel_expanded, how="left", lsuffix="_sample", rsuffix="_panel")

        df_expanded_with_panel['overlapping_panel_exon_start']=np.maximum(df_expanded_with_panel.event_exon_start_sample, df_expanded_with_panel.event_exon_start_panel)
        df_expanded_with_panel['overlapping_panel_exon_end']=np.minimum(df_expanded_with_panel.event_exon_end_sample, df_expanded_with_panel.event_exon_end_panel)

        df_expanded['event_exon_length']=df_expanded.event_exon_end - df_expanded.event_exon_start
        df_expanded_with_panel['overlapping_panel_exon_length']=np.maximum(df_expanded_with_panel.overlapping_panel_exon_end - df_expanded_with_panel.overlapping_panel_exon_start,0).fillna(0)

        df_panel_counts = df_expanded_with_panel.groupby(['ID_sample','sample_name_panel']).agg(
            {'overlapping_panel_exon_length':'sum'}
            ).reset_index('sample_name_panel').join(df_expanded.groupby('ID').agg({'event_exon_length':'sum'})
        )

        df_panel_counts['PANEL_COUNT'] = np.where(df_panel_counts.overlapping_panel_exon_length/df_panel_counts.event_exon_length>~{overlap_thresh}, 1, 0)
        df_panel_counts = df_panel_counts.groupby(df_panel_counts.index).agg({'PANEL_COUNT':'sum'})

        df = df.set_index('ID').join(df_panel_counts).fillna({'PANEL_COUNT':0})

        n_panel_samples = len(set(df_panel.sample_name))
        df['PANEL_FREQ'] = df['PANEL_COUNT']/n_panel_samples

        df_annotations = df[['contig','start','PANEL_FREQ','PANEL_COUNT']].copy()
        df_annotations = df_annotations.rename({"contig":"CHROM","start":"POS"}, axis=1)
        df_annotations = df_annotations.astype({"PANEL_COUNT":int})
        df_annotations.to_csv("~{output_basename}.annotations.tsv", index = False, sep="\t", float_format="%g")
        EOF

        bgzip ~{output_basename}.annotations.tsv
        tabix -s1 -b2 -e2 --skip-lines 1 ~{output_basename}.annotations.tsv.gz
        echo '##INFO=<ID=PANEL_FREQ,Number=1,Type=Float,Description="Frequency in panel">' > header_lines.txt
        echo '##INFO=<ID=PANEL_COUNT,Number=1,Type=Float,Description="Count in panel">' >> header_lines.txt
        bcftools annotate --no-version -a ~{output_basename}.annotations.tsv.gz -c CHROM,POS,PANEL_FREQ,PANEL_COUNT -h header_lines.txt -o  ~{output_basename}.vcf.gz ~{vcf}

        cp ~{output_basename}.vcf.gz tmp.vcf
        filters=('~{sep="' '" filter_expressions}')
        filter_names=(~{sep=" " filter_names})

        for i in ${!filters[@]}
        do
            eval bcftools filter --no-version -m + -e \'${filters[$i]}\' --soft-filter ${filter_names[$i]} -o tmp_out.vcf.gz tmp.vcf.gz
            mv tmp_out.vcf.gz tmp.vcf.gz
        done

        mv tmp.vcf.gz ~{output_basename}.filtered.genotyped-segments.vcf.gz

        bcftools index -t ~{output_basename}.filtered.genotyped-segments.vcf.gz

        md5sum ~{output_basename}.filtered.genotyped-segments.vcf.gz | awk '{ print $1 }' > ~{output_basename}.filtered.genotyped-segments.vcf.gz.md5sum

        n_total_events=$(bcftools view --no-header -e '(GT=="ref")' ~{output_basename}.filtered.genotyped-segments.vcf.gz | wc -l)
        n_pass_events=$(bcftools view --no-header -e '(GT=="ref") || (FILTER!~"PASS")' ~{output_basename}.filtered.genotyped-segments.vcf.gz | wc -l)

        if [ $n_total_events -le ~{max_events} ]; then
            if [ $n_pass_events -le ~{max_pass_events} ]; then
                echo "PASS" > qc_status.txt
            else
                echo "EXCESSIVE_NUMBER_OF_PASS_EVENTS" > qc_status.txt
            fi
        else
            echo "EXCESSIVE_NUMBER_OF_EVENTS" > qc_status.txt
        fi

        echo $n_total_events > n_total_events.txt
        echo $n_pass_events > n_pass_events.txt

        echo -e "sample\ttotal_cnv_events\tpassing_cnv_events" > ~{output_basename}.cnv_metrics.tsv
        echo -e "~{output_basename}\t${n_total_events}\t${n_pass_events}" >> ~{output_basename}.cnv_metrics.tsv
    >>>

    runtime {
            docker: "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
            disks: "local-disk " + disk_size_gb + " HDD"
            memory: mem_gb + " GiB"
    }

    output {
        File output_vcf = "~{output_basename}.filtered.genotyped-segments.vcf.gz"
        File output_vcf_index = "~{output_basename}.filtered.genotyped-segments.vcf.gz.tbi"
        File output_vcf_md5sum = "~{output_basename}.filtered.genotyped-segments.vcf.gz.md5sum"
        String qc_status = read_string("qc_status.txt")
        Int n_total_events = read_int("n_total_events.txt")
        Int n_pass_events = read_int("n_pass_events.txt")
        File cnv_metrics = "~{output_basename}.cnv_metrics.tsv"
    }

}

task AnnotateWithPoNFreq {
    input {
        File vcf
        File vcf_idx
        File annotations
        Int mem_gb=4
        Int disk_size_gb = 100
    }

    String output_basename = basename(vcf)
    command <<<
        set -euo pipefail
        bgzip ~{annotations}
        tabix -s1 -b2 -e2 --skip-lines 1 ~{annotations}.gz
        echo '##INFO=<ID=PANEL_FREQ,Number=1,Type=Float,Description="Frequency in panel">' > header_lines.txt
        echo '##INFO=<ID=PANEL_COUNT,Number=1,Type=Float,Description="Count in panel">' >> header_lines.txt
        bcftools annotate --no-version -a ~{annotations}.gz -c CHROM,POS,PANEL_FREQ,PANEL_COUNT -h header_lines.txt -o  ~{output_basename}.vcf.gz ~{vcf}

    >>>

    runtime {
            docker: "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
            disks: "local-disk " + disk_size_gb + " HDD"
            memory: mem_gb + " GiB"
    }

    output {
        File output_vcf = "~{output_basename}.vcf.gz"
    }
}


task FilterVCF {
    input {
        String samplename
        File vcf_file
        Array[String] filter_expressions = ['(GT=="alt" | GT=="mis") & ((FMT/CN>1 & QUAL<50) | (FMT/CN==1 & QUAL<100 ) | (FMT/CN==0 & QUAL<400))','(GT=="alt" | GT=="mis") & (INFO/PANEL_COUNT>1)']
        String gatk_docker
        Array[String] filter_names = ["LowQual","PanelCount"]
    }

    command <<<
        set -euo pipefail
        cp ~{vcf_file} tmp.vcf.gz
        filters=('~{sep="' '" filter_expressions}')
        filter_names=(~{sep=" " filter_names})

        for i in ${!filters[@]}
        do
            eval bcftools filter --no-version -m + -e \'${filters[$i]}\' --soft-filter ${filter_names[$i]} -o tmp_out.vcf.gz tmp.vcf.gz
            mv tmp_out.vcf.gz tmp.vcf.gz
        done

        mv tmp.vcf.gz ~{samplename}.filtered.genotyped-segments.vcf.gz

        bcftools index -t ~{samplename}.filtered.genotyped-segments.vcf.gz

        md5sum ~{samplename}.filtered.genotyped-segments.vcf.gz | awk '{ print $1 }' > ~{samplename}.filtered.genotyped-segments.vcf.gz.md5sum
    >>>

    output {
        File filtered_vcf = samplename + ".filtered.genotyped-segments.vcf.gz"
        File filtered_vcf_index = samplename + ".filtered.genotyped-segments.vcf.gz.tbi"
        File filtered_vcf_md5sum = samplename + ".filtered.genotyped-segments.vcf.gz.md5sum"
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        memory: "8G"
        cpu: 2
    }
}

task QCFilteredVCF {
    input {
        File filtered_vcf
        String samplename
        Int max_events
        Int max_pass_events
    }

    command <<<
        set -euo pipefail

        n_total_events=$(bcftools view --no-header -e '(GT=="ref")' ~{filtered_vcf} | wc -l)
        n_pass_events=$(bcftools view --no-header -e '(GT=="ref") || (FILTER!~"PASS")' ~{filtered_vcf} | wc -l)

        if [ $n_total_events -le ~{max_events} ]; then
            if [ $n_pass_events -le ~{max_pass_events} ]; then
                echo "PASS" > qc_status.txt
            else
                echo "EXCESSIVE_NUMBER_OF_PASS_EVENTS" > qc_status.txt
            fi
        else
            echo "EXCESSIVE_NUMBER_OF_EVENTS" > qc_status.txt
        fi

        echo $n_total_events > n_total_events.txt
        echo $n_pass_events > n_pass_events.txt

        echo -e "sample\ttotal_cnv_events\tpassing_cnv_events" > ~{samplename}.cnv_metrics.tsv
        echo -e "~{samplename}\t${n_total_events}\t${n_pass_events}" >> ~{samplename}.cnv_metrics.tsv
    >>>

    output {
        String qc_status = read_string("qc_status.txt")
        Int n_total_events = read_int("n_total_events.txt")
        Int n_pass_events = read_int("n_pass_events.txt")
        File cnv_metrics = "~{samplename}.cnv_metrics.tsv"
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        memory: "8G"
        cpu: 2
    }
}