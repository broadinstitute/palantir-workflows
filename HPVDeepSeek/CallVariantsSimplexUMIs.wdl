version 1.0

task CollectSequencingArtifactMetrics {
    input {
        String output_basename
        File bam
        File bai
        File reference
        File reference_fai

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 128)
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        gatk --java-options "-Xmx14g" \
        CollectSequencingArtifactMetrics \
        --INPUT ~{bam} \
        --OUTPUT ~{output_basename} \
        --REFERENCE_SEQUENCE ~{reference} \
        --VALIDATION_STRINGENCY LENIENT
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.6.2.0"
    }

    output {
        File pre_adapter_metrics = "~{output_basename}.pre_adapter_detail_metrics"
    }
}

task Mutect2 {
    input {
        String output_basename
        File tumor_bam
        File tumor_bai
        File reference
        File reference_fai
        File reference_dict
        File gnomad
        File gnomad_idx
        File pon
        File pon_idx
        File intervals
        File variants_for_contamination
        File variants_for_contamination_idx

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(tumor_bam, "GiB")) + 100)
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        gatk --java-options "-Xms8g -Xmx14g" \
        Mutect2 \
        --input ~{tumor_bam} \
        --output ~{output_basename}.vcf.gz \
        --reference ~{reference} \
        --germline-resource ~{gnomad} \
        --panel-of-normals ~{pon} \
        --intervals ~{intervals} \
        --f1r2-tar-gz ~{output_basename}.f1r2.tar.gz \
        --read-filter NotSupplementaryAlignmentReadFilter \
        --dont-use-soft-clipped-bases true \
        --af-of-alleles-not-in-resource 0.001 \
        --tumor-lod-to-emit 0 \
        --max-reads-per-alignment-start 0 \
        --initial-tumor-lod 0 \
        --pcr-snv-qual 70 \
        --pcr-indel-qual 70

        mutect2_exit_code=$?

        set +e

        gatk --java-options "-Xms8g -Xmx14g" \
        GetPileupSummaries \
        --input ~{tumor_bam} \
        --output ~{output_basename}.tumor-pileups.table \
        --reference ~{reference} \
        --interval-set-rule INTERSECTION \
        --intervals ~{intervals} \
        --variant ~{variants_for_contamination} \
        --intervals ~{variants_for_contamination}

        exit $mutect2_exit_code
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.6.2.0"
    }

    output {
        File unfiltered_vcf = "~{output_basename}.vcf.gz"
        File unfiltered_vcf_idx = "~{output_basename}.vcf.gz.tbi"
        File mutect2_stats = "~{output_basename}.vcf.gz.stats"
        File f1r2_counts = "~{output_basename}.f1r2.tar.gz"
        File tumor_pileups = "~{output_basename}.tumor-pileups.table"
    }
}

task LearnReadOrientationModel {
    input {
        String output_basename
        File f1r2_tar_gz

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = 128
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        gatk --java-options "-Xms8g -Xmx14g" \
        LearnReadOrientationModel \
        --input ~{f1r2_tar_gz} \
        --output ~{output_basename}.artifact-priors.tar.gz
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.6.2.0"
    }

    output {
        File artifact_prior_table = "~{output_basename}.artifact-priors.tar.gz"
    }
}

task CalculateContamination {
    input {
        String output_basename
        File tumor_pileups

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = 128
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        gatk ---java-options "-Xms8g -Xmx14g" \
        CalculateContamination \
        --input ~{tumor_pileups} \
        --output ~{output_basename}.contamination.table \
        --tumor-segmentation ~{output_basename}.segments.table
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.6.2.0"
    }

    output {
        File contamination_table = "~{output_basename}.contamination.table"
        File maf_segments = "~{output_basename}.segments.table"
    }
}

task FilterMutectCalls {
    input {
        String output_basename
        File unfiltered_vcf
        File unfiltered_vcf_idx
        File reference
        File reference_fai
        File reference_dict
        File mutect2_stats
        File artifact_priors_tar_gz
        File contamination_table
        File maf_segments
        Boolean compress = true

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = 128
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        gatk --java-options "-Xms8g -Xmx14g" \
        FilterMutectCalls \
        --variant ~{unfiltered_vcf} \
        --reference ~{reference} \
        --output ~{output_basename}.vcf.gz \
        --contamination-table ~{contamination_table} \
        --tumor-segmentation ~{maf_segments} \
        --ob-priors ~{artifact_priors_tar_gz} \
        --stats ~{mutect2_stats} \
        --filtering-stats ~{output_basename}.filtering.stats
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.6.2.0"
    }

    output {
        File filtered_vcf = "~{output_basename}.vcf.gz"
        File filtered_vcf_idx = "~{output_basename}.vcf.gz.tbi"
        File filtering_stats = "~{output_basename}.filtering.stats"
    }
}

task FilterAlignmentArtifacts {
    input {
        String output_basename
        File bam
        File bai
        File input_vcf
        File input_vcf_idx
        File reference
        File reference_fai
        File reference_dict
        File realignment_index_bundle
        Boolean compress = true

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = 128
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        gatk --java-options "-Xms8g -Xmx14g" \
        FilterAlignmentArtifacts \
        --input ~{bam} \
        --output ~{output_basename}.vcf.gz
        --variant ~{input_vcf} \
        --reference ~{reference} \
        --bwa-mem-index-image ~{realignment_index_bundle}
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.6.2.0"
    }

    output {
        File filtered_vcf = "~{output_basename}.vcf.gz"
        File filtered_vcf_idx = "~{output_basename}.vcf.gz.tbi"
    }
}

task RunMappingFilter {
    input {
        String output_basename
        File vcf
        File vcf_idx
        File reference
        File reference_fai
        File reference_dict
        File blastdb_nhr
        File blastdb_nin
        File blastdb_nsq
        String blastn_path
        String mapping_filter_python_script

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = 512
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    command <<<
        python2.7 ~{mapping_filter_python_script} \
        --vcf ~{vcf} \
        --outfile ~{output_basename}.filtered.vcf \
        --reference_fasta ~{reference} \
        --blastn ~{blastn_path}

        bgzip ~{output_basename}.filtered.vcf
        tabix ~{output_basename}.filtered.vcf.gz
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.6.2.0"
    }

    output {
        File map_filtered_vcf = "~{output_basename}.filtered.vcf.gz"
        File map_filtered_vcf_idx = "~{output_basename}.filtered.vcf.gz.tbi"
    }
}

task VariantFiltration {
    input {
        String output_basename
        File mutect_filtered_vcf
        File mutect_filtered_vcf_idx
        File filter_vcf
        File filter_vcf_idx
        File reference
        File reference_fai
        File reference_dict
        String mapping_filter_name
        Boolean filter_not_in_mask

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = 512
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    # Filter out variant with POPAF < 3.0 (this is phred scaled), which is
    # equivalent to filtering population allele frequency > 0.001.
    command <<<
        gatk VariantFiltration \
        --reference ~{reference} \
        --variant ~{mutect_filtered_vcf} \
        --output popaf.filtered.vcf.gz \
        --filter-expression "POPAF < 3.0" \
        --filter-name germline

        gatk VariantFiltration \
        --reference ~{reference} \
        --variant popaf.filtered.vcf.gz \
        --output ~{output_basename}.~{mapping_filter_name}.vcf.gz \
        --mask ~{filter_vcf} \
        --filter-not-in-mask ~{filter_not_in_mask} \
        --mask-name ~{mapping_filter_name}
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.6.2.0"
    }

    output {
        File output_vcf = "~{output_basename}.~{mapping_filter_name}.vcf.gz"
        File output_vcf_idx = "~{output_basename}.~{mapping_filter_name}.vcf.gz.tbi"
    }
}

task Funcotate {
    input {
        String output_basename
        File reference
        File reference_fai
        File reference_dict
        File vcf
        File vcf_idx
        File funcotator_data_source
        String reference_version
        String output_format
        Boolean compress
        Boolean use_gnomad
        String? control_id
        String? case_id
        String? sequencing_center
        String? sequence_source
        String? transcript_selection_mode
        File? transcript_selection_list
        Array[String]? annotation_defaults
        Array[String]? annotation_overrides
        Array[String]? funcotator_excluded_fields
        Boolean? filter_funcotations
        File? interval_list
        String? extra_args

        Int? cpu = 1
        Int? memory_gb = 16
        Int? disk_size_gb = 512
        Int? min_ssd_size_gb = 512
        Boolean use_ssd = true
    }

    String output_maf = output_basename + ".maf"
    String output_maf_index = output_maf + ".idx"
    String output_vcf = output_basename + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"
    String output_file = if output_format == "MAF" then output_maf else output_vcf
    String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx
    String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
    String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
    String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
    String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
    String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
    String interval_list_arg = if defined(interval_list) then " -L " else ""
    String extra_args_arg = select_first([extra_args, ""])
    String dollar = "$"

    command <<<
        echo "Extracting data sources zip file..."
        mkdir datasources_dir
        tar xzvf ~{funcotator_data_source} -C datasources_dir --strip-components 1
        DATA_SOURCES_FOLDER="$PWD/datasources_dir"

        if ~{use_gnomad} ; then
            echo "Enabling gnomAD..."
            for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                if [[ -f ~{dollar}{DATA_SOURCES_FOLDER}/~{dollar}{potential_gnomad_gz} ]] ; then
                    cd ~{dollar}{DATA_SOURCES_FOLDER}
                    tar -zvxf ~{dollar}{potential_gnomad_gz}
                    cd -
                else
                    echo "ERROR: Cannot find gnomAD folder: ~{dollar}{potential_gnomad_gz}" 1>&2
                    false
                fi
            done
        fi

        gatk --java-options "-Xmx14g" \
        Funcotator \
        --data-sources-path $DATA_SOURCES_FOLDER \
        --ref-version ~{reference_version} \
        --output-file-format ~{output_format} \
        --reference ~{reference} \
        --variant ~{vcf} \
        --output ~{output_file} \
        ~{interval_list_arg} ~{default="" interval_list} \
        --annotation-default normal_barcode:~{default="Unknown" control_id} \
        --annotation-default tumor_barcode:~{default="Unknown" case_id} \
        --annotation-default Center:~{default="Unknown" sequencing_center} \
        --annotation-default source:~{default="Unknown" sequence_source} \
        ~{"--transcript-selection-mode " + transcript_selection_mode} \
        ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
        ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \
        ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \
        ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
        ~{filter_funcotations_args} \
        --prefer-mane-transcripts \
        ~{extra_args_arg}

        if [[ "~{output_format}" == "MAF" ]] ; then
            touch ~{output_maf_index}
        fi
    >>>

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk" + if use_ssd then " ~{min_ssd_size_gb} SSD" else " ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/gatk-custom:v4.6.2.0a"
    }

    output {
        File funcotated_output_file = "~{output_file}"
        File funcotated_output_file_index = "~{output_file_index}"
    }
}

workflow CallVariantsSimplexUMIs {
    input {
        String output_basename
        File tumor_bam
        File tumor_bai
        File target_intervals
        File reference
        File reference_fai
        File reference_dict
        File gnomad
        File gnomad_idx
        File pon
        File pon_idx
        File variants_for_contamination
        File variants_for_contamination_idx
        #File realignment_index_bundle
        #String mapping_filter_python_script
        #File blastdb_nhr
        #File blastdb_nin
        #File blastdb_nsq
        #String blastn_path
        File funcotator_data_source
    }

    call CollectSequencingArtifactMetrics {
        input:
            bam = tumor_bam,
            bai = tumor_bai,
            reference = reference,
            reference_fai = reference_fai,
            output_basename = output_basename
    }

    call Mutect2 {
        input:
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            gnomad = gnomad,
            gnomad_idx = gnomad_idx,
            pon = pon,
            pon_idx = pon_idx,
            intervals = target_intervals,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_idx = variants_for_contamination_idx,
            output_basename = output_basename
    }

    call LearnReadOrientationModel {
        input:
            f1r2_tar_gz = Mutect2.f1r2_counts,
            output_basename = output_basename
    }

    call CalculateContamination {
        input:
            tumor_pileups = Mutect2.tumor_pileups,
            output_basename = output_basename
    }

    call FilterMutectCalls {
        input:
            unfiltered_vcf = Mutect2.unfiltered_vcf,
            unfiltered_vcf_idx = Mutect2.unfiltered_vcf_idx,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            mutect2_stats = Mutect2.mutect2_stats,
            contamination_table = CalculateContamination.contamination_table,
            maf_segments = CalculateContamination.maf_segments,
            artifact_priors_tar_gz = LearnReadOrientationModel.artifact_prior_table,
            output_basename = output_basename
    }

# Not sure if this filter is even needed, but the realignment index bundle can be generated with gatk BwaMemIndexImageCreator if needed
#    call FilterAlignmentArtifacts {
#        input:
#            bam = tumor_bam,
#            bai = tumor_bai,
#            input_vcf = FilterMutectCalls.filtered_vcf,
#            input_vcf_idx = FilterMutectCalls.filtered_vcf_idx,
#            reference = reference,
#            reference_fai = reference_fai,
#            reference_dict = reference_dict,
#            realignment_index_bundle = realignment_index_bundle,
#            output_basename = output_basename
#    }

#    call RunMappingFilter {
#        input:
#            vcf = FilterMutectCalls.filtered_vcf,
#            vcf_idx = FilterMutectCalls.filtered_vcf_idx,
#            reference = reference,
#            reference_fai = reference_fai,
#            reference_dict = reference_dict,
#            blastdb_nhr = blastdb_nhr,
#            blastdb_nin = blastdb_nin,
#            blastdb_nsq = blastdb_nsq,
#            blastn_path = blastn_path,
#            mapping_filter_python_script = mapping_filter_python_script,
#            output_basename = output_basename
#    }

#    call VariantFiltration {
#        input:
#            mutect_filtered_vcf = FilterMutectCalls.filtered_vcf,
#            mutect_filtered_vcf_idx = FilterMutectCalls.filtered_vcf_idx,
#            filter_vcf = RunMappingFilter.map_filtered_vcf,
#            filter_vcf_idx = RunMappingFilter.map_filtered_vcf_idx,
#            reference = reference,
#            reference_fai = reference_fai,
#            reference_dict = reference_dict,
#            mapping_filter_name = "mapping_filter",
#            filter_not_in_mask = true,
#            output_basename = output_basename
#    }

#    call Funcotate {
#        input:
#            vcf = VariantFiltration.output_vcf,
#            vcf_idx = VariantFiltration.output_vcf_idx,
#            reference = reference,
#            reference_fai = reference_fai,
#            reference_dict = reference_dict,
#            funcotator_data_source = funcotator_data_source,
#            reference_version = "hg38",
#            filter_funcotations = true,
#            use_gnomad = false,
#            output_format = "MAF",
#            compress = true,
#            transcript_selection_mode = "BEST_EFFECT",
#            output_basename = output_basename
#    }

    call Funcotate {
        input:
            vcf = FilterMutectCalls.filtered_vcf,
            vcf_idx = FilterMutectCalls.filtered_vcf_idx,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            funcotator_data_source = funcotator_data_source,
            reference_version = "hg38",
            filter_funcotations = true,
            use_gnomad = false,
            output_format = "MAF",
            compress = true,
            transcript_selection_mode = "BEST_EFFECT",
            output_basename = output_basename
    }

    output {
        File unfiltered_vcf = Mutect2.unfiltered_vcf
        File unfiltered_vcf_idx = Mutect2.unfiltered_vcf_idx
        File mutect2_stats = Mutect2.mutect2_stats
        File filtered_vcf = FilterMutectCalls.filtered_vcf
        File filtered_vcf_idx = FilterMutectCalls.filtered_vcf_idx
        File filtering_stats = FilterMutectCalls.filtering_stats
        #File filtered_vcf = VariantFiltration.output_vcf
        #File filtered_vcf_idx = VariantFiltration.output_vcf_idx
        File funcotated_maf = Funcotate.funcotated_output_file
    }
}