version 1.0

import "cnv_common_tasks.wdl" as CNVTasks

workflow SingleSampleGCNVAndFilterVCFs {
    input {
        Array[String]? filter_expressions #= "QUAL < 50"  # Example filter criteria
        Array[String]? filter_names

        ##################################
        ##for case mode
        #### required basic arguments ####
        ##################################
        File preprocessed_intervals
        File filtered_intervals
        String normal_bam
        String normal_bai
        File contig_ploidy_model_tar
        File gcnv_model_tar
        Array[File]+ pon_genotyped_segments_vcfs
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

    call CNVTasks.CollectCounts {
        input:
        intervals = preprocessed_intervals,
        bam = normal_bam,
        bam_idx = normal_bai,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_fasta_dict = ref_fasta_dict,
        format = collect_counts_format,
        enable_indexing = collect_counts_enable_indexing,
        disabled_read_filters = disabled_read_filters_for_collect_counts,
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        mem_gb = mem_gb_for_collect_counts,
        preemptible_attempts = preemptible_attempts,
        gcs_project_for_requester_pays = gcs_project_for_requester_pays
    }

    call DetermineGermlineContigPloidyCaseMode {
        input:
            read_counts_file = CollectCounts.counts,
            contig_ploidy_model_tar = contig_ploidy_model_tar,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_determine_germline_contig_ploidy,
            cpu = cpu_for_determine_germline_contig_ploidy,
            disk_space_gb = disk_for_determine_germline_contig_ploidy,
            mapping_error_rate = ploidy_mapping_error_rate,
            sample_psi_scale = ploidy_sample_psi_scale,
            preemptible_attempts = preemptible_attempts
    }

    call GermlineCNVCallerCaseMode {
        input:
            read_counts_file = CollectCounts.counts,
            contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
            gcnv_model_tar = gcnv_model_tar,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_germline_cnv_caller,
            cpu = cpu_for_germline_cnv_caller,
            p_alt = gcnv_p_alt,
            cnv_coherence_length = gcnv_cnv_coherence_length,
            max_copy_number = gcnv_max_copy_number,
            mapping_error_rate = gcnv_mapping_error_rate,
            sample_psi_scale = gcnv_sample_psi_scale,
            depth_correction_tau = gcnv_depth_correction_tau,
            copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
            active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
            learning_rate = gcnv_learning_rate,
            adamax_beta_1 = gcnv_adamax_beta_1,
            adamax_beta_2 = gcnv_adamax_beta_2,
            log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
            log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
            log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
            max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
            max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
            min_training_epochs = gcnv_min_training_epochs,
            max_training_epochs = gcnv_max_training_epochs,
            initial_temperature = gcnv_initial_temperature,
            num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
            convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
            convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
            convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
            max_calling_iters = gcnv_max_calling_iters,
            caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
            caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
            caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
            disable_annealing = gcnv_disable_annealing,
            preemptible_attempts = preemptible_attempts
    }

     call CNVTasks.PostprocessGermlineCNVCalls {
        input:
            entity_id = CollectCounts.entity_id,
            gcnv_calls_tars = [GermlineCNVCallerCaseMode.gcnv_call_tar],
            gcnv_model_tars = [gcnv_model_tar],
            calling_configs = [GermlineCNVCallerCaseMode.calling_config_json],
            denoising_configs = [GermlineCNVCallerCaseMode.denoising_config_json],
            gcnvkernel_version = [GermlineCNVCallerCaseMode.gcnvkernel_version_json],
            sharded_interval_lists = [GermlineCNVCallerCaseMode.sharded_interval_list],
            allosomal_contigs = allosomal_contigs,
            ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
            contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
            sample_index = 0,
            maximum_number_events = maximum_number_events_per_sample,
            maximum_number_pass_events = maximum_number_pass_events_per_sample,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

#    File vcf = CNVGermlineCaseWorkflow.genotyped_segments_vcfs[0]
#    File vcf_index = CNVGermlineCaseWorkflow.genotyped_intervals_vcf_indexes[0]
#    String samplename = CNVGermlineCaseWorkflow.entity_id[0]

    call ExtractPoNFreqAnnotateFilterAndQC {
        input:
            vcf = PostprocessGermlineCNVCalls.genotyped_segments_vcf,
            vcf_idx = PostprocessGermlineCNVCalls.genotyped_segments_vcf_index,
            panel_vcfs = pon_genotyped_segments_vcfs,
            intervals = preprocessed_intervals,
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

        File read_counts_entity_id = CollectCounts.entity_id
        File read_counts = CollectCounts.counts
        File contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar
        File gcnv_call_tar = GermlineCNVCallerCaseMode.gcnv_call_tar
        File gcnv_tracking_tar = GermlineCNVCallerCaseMode.gcnv_tracking_tar
        File genotyped_intervals_vcf = PostprocessGermlineCNVCalls.genotyped_intervals_vcf
        File genotyped_intervals_vcf_index = PostprocessGermlineCNVCalls.genotyped_intervals_vcf_index
        File genotyped_segments_vcf = PostprocessGermlineCNVCalls.genotyped_segments_vcf
        File genotyped_segments_vcf_index = PostprocessGermlineCNVCalls.genotyped_segments_vcf_index
        String qc_status_string = ExtractPoNFreqAnnotateFilterAndQC.qc_status
        Boolean qc_passed = (ExtractPoNFreqAnnotateFilterAndQC.qc_status == "PASS")
        File cnv_metrics = ExtractPoNFreqAnnotateFilterAndQC.cnv_metrics
        File denoised_copy_ratios = PostprocessGermlineCNVCalls.denoised_copy_ratios
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
        File filtered_vcf = "~{output_basename}.filtered.genotyped-segments.vcf.gz"
        File filtered_vcf_index = "~{output_basename}.filtered.genotyped-segments.vcf.gz.tbi"
        File filtered_vcf_md5sum = "~{output_basename}.filtered.genotyped-segments.vcf.gz.md5sum"
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

task DetermineGermlineContigPloidyCaseMode {
    input {
      File read_counts_file
      File contig_ploidy_model_tar
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts

      # Model parameters
      Float? mapping_error_rate
      Float? sample_psi_scale
    }

    # We do not expose Hybrid ADVI parameters -- the default values are decent

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{default=8 cpu}
        export OMP_NUM_THREADS=~{default=8 cpu}

        mkdir contig-ploidy-model
        tar xzf ~{contig_ploidy_model_tar} -C contig-ploidy-model

        gatk --java-options "-Xmx~{command_mem_mb}m" DetermineGermlineContigPloidy \
            --input ~{read_counts_file} \
            --model contig-ploidy-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --mapping-error-rate ~{default="0.3" mapping_error_rate} \
            --sample-psi-scale ~{default="0.0001" sample_psi_scale}

        tar c -C ~{output_dir_}/case-calls . | gzip -1 > case-contig-ploidy-calls.tar.gz

        rm -rf contig-ploidy-model
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File contig_ploidy_calls_tar = "case-contig-ploidy-calls.tar.gz"
    }
}

task GermlineCNVCallerCaseMode {
    input {
      File read_counts_file
      File contig_ploidy_calls_tar
      File gcnv_model_tar
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts

      # Caller parameters
      Float? p_alt
      Float? cnv_coherence_length
      Int? max_copy_number

      # Denoising model parameters
      Float? mapping_error_rate
      Float? sample_psi_scale
      Float? depth_correction_tau
      String? copy_number_posterior_expectation_mode
      Int? active_class_padding_hybrid_mode

      # Hybrid ADVI parameters
      Float? learning_rate
      Float? adamax_beta_1
      Float? adamax_beta_2
      Int? log_emission_samples_per_round
      Float? log_emission_sampling_median_rel_error
      Int? log_emission_sampling_rounds
      Int? max_advi_iter_first_epoch
      Int? max_advi_iter_subsequent_epochs
      Int? min_training_epochs
      Int? max_training_epochs
      Float? initial_temperature
      Int? num_thermal_advi_iters
      Int? convergence_snr_averaging_window
      Float? convergence_snr_trigger_threshold
      Int? convergence_snr_countdown_window
      Int? max_calling_iters
      Float? caller_update_convergence_threshold
      Float? caller_internal_admixing_rate
      Float? caller_external_admixing_rate
      Boolean? disable_annealing
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{default=8 cpu}
        export OMP_NUM_THREADS=~{default=8 cpu}

        mkdir contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls

        mkdir gcnv-model
        tar xzf ~{gcnv_model_tar} -C gcnv-model

        gatk --java-options "-Xmx~{command_mem_mb}m"  GermlineCNVCaller \
            --run-mode CASE \
            --input ~{read_counts_file} \
            --contig-ploidy-calls contig-ploidy-calls \
            --model gcnv-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --p-alt ~{default="5e-4" p_alt} \
            --cnv-coherence-length ~{default="10000.0" cnv_coherence_length} \
            --max-copy-number ~{default="5" max_copy_number} \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --sample-psi-scale ~{default="0.01" sample_psi_scale} \
            --depth-correction-tau ~{default="10000.0" depth_correction_tau} \
            --copy-number-posterior-expectation-mode ~{default="HYBRID" copy_number_posterior_expectation_mode} \
            --active-class-padding-hybrid-mode ~{default="50000" active_class_padding_hybrid_mode} \
            --learning-rate ~{default="0.05" learning_rate} \
            --adamax-beta-1 ~{default="0.9" adamax_beta_1} \
            --adamax-beta-2 ~{default="0.99" adamax_beta_2} \
            --log-emission-samples-per-round ~{default="50" log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error ~{default="0.005" log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds ~{default="10" log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch ~{default="5000" max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs ~{default="100" max_advi_iter_subsequent_epochs} \
            --min-training-epochs ~{default="10" min_training_epochs} \
            --max-training-epochs ~{default="100" max_training_epochs} \
            --initial-temperature ~{default="2.0" initial_temperature} \
            --num-thermal-advi-iters ~{default="2500" num_thermal_advi_iters} \
            --convergence-snr-averaging-window ~{default="500" convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold ~{default="0.1" convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window ~{default="10" convergence_snr_countdown_window} \
            --max-calling-iters ~{default="10" max_calling_iters} \
            --caller-update-convergence-threshold ~{default="0.001" caller_update_convergence_threshold} \
            --caller-internal-admixing-rate ~{default="0.75" caller_internal_admixing_rate} \
            --caller-external-admixing-rate ~{default="1.00" caller_external_admixing_rate} \
            --disable-annealing ~{default="false" disable_annealing}

        tar czf case-gcnv-tracking.tar.gz -C ~{output_dir_}/case-tracking .

        tar czf case-gcnv-calls.tar.gz -C ~{output_dir_}/case-calls/SAMPLE_0 .

        rm -rf contig-ploidy-calls
        rm -rf gcnv-model
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File gcnv_call_tar = "case-gcnv-calls.tar.gz"
        File gcnv_tracking_tar = "case-gcnv-tracking.tar.gz"
        File calling_config_json = "~{output_dir_}/case-calls/calling_config.json"
        File denoising_config_json = "~{output_dir_}/case-calls/denoising_config.json"
        File gcnvkernel_version_json = "~{output_dir_}/case-calls/gcnvkernel_version.json"
        File sharded_interval_list = "~{output_dir_}/case-calls/interval_list.tsv"
    }
}

task GermlineCNVCallerCaseModeAndPostProcess {
    input {
      File read_counts_file
      File contig_ploidy_calls_tar
      File gcnv_model_tar
      Array[String] allosomal_contigs
      Int ref_copy_number_autosomal_contigs
      String entity_id
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts

      # Caller parameters
      Float? p_alt
      Float? cnv_coherence_length
      Int? max_copy_number

      # Denoising model parameters
      Float? mapping_error_rate
      Float? sample_psi_scale
      Float? depth_correction_tau
      String? copy_number_posterior_expectation_mode
      Int? active_class_padding_hybrid_mode

      # Hybrid ADVI parameters
      Float? learning_rate
      Float? adamax_beta_1
      Float? adamax_beta_2
      Int? log_emission_samples_per_round
      Float? log_emission_sampling_median_rel_error
      Int? log_emission_sampling_rounds
      Int? max_advi_iter_first_epoch
      Int? max_advi_iter_subsequent_epochs
      Int? min_training_epochs
      Int? max_training_epochs
      Float? initial_temperature
      Int? num_thermal_advi_iters
      Int? convergence_snr_averaging_window
      Float? convergence_snr_trigger_threshold
      Int? convergence_snr_countdown_window
      Int? max_calling_iters
      Float? caller_update_convergence_threshold
      Float? caller_internal_admixing_rate
      Float? caller_external_admixing_rate
      Boolean? disable_annealing
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    Array[String] allosomal_contigs_args = if defined(allosomal_contigs) then prefix("--allosomal-contig ", select_first([allosomal_contigs])) else []
    String genotyped_intervals_vcf_filename = "genotyped-intervals-~{entity_id}.vcf.gz"
    String genotyped_segments_vcf_filename = "genotyped-segments-~{entity_id}.vcf.gz"
    String denoised_copy_ratios_filename = "denoised_copy_ratios-~{entity_id}.tsv"

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{default=8 cpu}
        export OMP_NUM_THREADS=~{default=8 cpu}

        mkdir contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls

        mkdir gcnv-model
        tar xzf ~{gcnv_model_tar} -C gcnv-model

        gatk --java-options "-Xmx~{command_mem_mb}m"  GermlineCNVCaller \
            --run-mode CASE \
            --input ~{read_counts_file} \
            --contig-ploidy-calls contig-ploidy-calls \
            --model gcnv-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --p-alt ~{default="5e-4" p_alt} \
            --cnv-coherence-length ~{default="10000.0" cnv_coherence_length} \
            --max-copy-number ~{default="5" max_copy_number} \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --sample-psi-scale ~{default="0.01" sample_psi_scale} \
            --depth-correction-tau ~{default="10000.0" depth_correction_tau} \
            --copy-number-posterior-expectation-mode ~{default="HYBRID" copy_number_posterior_expectation_mode} \
            --active-class-padding-hybrid-mode ~{default="50000" active_class_padding_hybrid_mode} \
            --learning-rate ~{default="0.05" learning_rate} \
            --adamax-beta-1 ~{default="0.9" adamax_beta_1} \
            --adamax-beta-2 ~{default="0.99" adamax_beta_2} \
            --log-emission-samples-per-round ~{default="50" log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error ~{default="0.005" log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds ~{default="10" log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch ~{default="5000" max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs ~{default="100" max_advi_iter_subsequent_epochs} \
            --min-training-epochs ~{default="10" min_training_epochs} \
            --max-training-epochs ~{default="100" max_training_epochs} \
            --initial-temperature ~{default="2.0" initial_temperature} \
            --num-thermal-advi-iters ~{default="2500" num_thermal_advi_iters} \
            --convergence-snr-averaging-window ~{default="500" convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold ~{default="0.1" convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window ~{default="10" convergence_snr_countdown_window} \
            --max-calling-iters ~{default="10" max_calling_iters} \
            --caller-update-convergence-threshold ~{default="0.001" caller_update_convergence_threshold} \
            --caller-internal-admixing-rate ~{default="0.75" caller_internal_admixing_rate} \
            --caller-external-admixing-rate ~{default="1.00" caller_external_admixing_rate} \
            --disable-annealing ~{default="false" disable_annealing}

        gatk --java-options "-Xmx~{command_mem_mb}m" PostprocessGermlineCNVCalls \
                    --calls-shard-path ~{output_dir_}/case-calls \
                    --model-shard-path gcnv-model \
                    ~{sep=" " allosomal_contigs_args} \
                    --autosomal-ref-copy-number ~{ref_copy_number_autosomal_contigs} \
                    --contig-ploidy-calls contig-ploidy-calls \
                    --sample-index 0 \
                    --output-genotyped-intervals genotyped-intervals-~{entity_id}.vcf.gz \
                    --output-genotyped-segments genotyped-segments-~{entity_id}.vcf.gz \
                    --output-denoised-copy-ratios "denoised_copy_ratios-~{entity_id}.tsv"


        tar czf case-gcnv-tracking.tar.gz -C ~{output_dir_}/case-tracking .

        tar czf case-gcnv-calls.tar.gz -C ~{output_dir_}/case-calls/SAMPLE_0 .

    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File gcnv_call_tar = "case-gcnv-calls.tar.gz"
        File gcnv_tracking_tar = "case-gcnv-tracking.tar.gz"
        File calling_config_json = "~{output_dir_}/case-calls/calling_config.json"
        File denoising_config_json = "~{output_dir_}/case-calls/denoising_config.json"
        File gcnvkernel_version_json = "~{output_dir_}/case-calls/gcnvkernel_version.json"
        File sharded_interval_list = "~{output_dir_}/case-calls/interval_list.tsv"

        File genotyped_intervals_vcf = "genotyped-intervals-~{entity_id}.vcf.gz"
        File genotyped_intervals_vcf_index = "genotyped-intervals-~{entity_id}.vcf.gz.tbi"
        File genotyped_segments_vcf = "genotyped-segments-~{entity_id}.vcf.gz"
        File genotyped_segments_vcf_index = "genotyped-segments-~{entity_id}.vcf.gz.tbi"
        File denoised_copy_ratios = "denoised_copy_ratios-~{entity_id}.tsv"
    }
}


