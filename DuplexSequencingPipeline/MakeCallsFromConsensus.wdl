import "https://api.firecloud.org/ga4gh/v1/tools/LiquidBiopsyDevelopment:Mutect2/versions/6/plain-WDL/descriptor" as m2

workflow MakeCallsFromConsensus {

   # dockers and override jars
   String bloodbiopsydocker 
   String reference_version
   Int? preemptible_attempts

   # Mutect2 specific arguments
   String basename
   File tumor_bam
   File tumor_bam_idx
   String tumor_sample_name
   File? normal_bam
   File? normal_bam_idx
   String? normal_sample_name
   File? raw_tumor_bam
   File? raw_tumor_bam_idx
   File? raw_normal_bam 
   File? raw_normal_bam_idx

   File target_intervals
   File reference 
   File reference_idx 
   File reference_dict 
   File? pon
   File? pon_idx
   File? variants_for_contamination
   File? variants_for_contamination_idx
   File? gnomad
   File? gnomad_idx
   Boolean? run_orientation_bias_mixture_model_filter
   Boolean run_ob_filter = select_first([run_orientation_bias_mixture_model_filter, false])
   Int scatter_count
   String? m2_extra_args
   String? m2_extra_filtering_args
   String m2_extra_filtering_args_or_default = select_first([m2_extra_filtering_args, ""])
   Boolean? compress_vcfs
   Boolean? make_bamout
   Boolean make_bamout_or_default = select_first([make_bamout, false])
   String mapping_filter_python_script 
   File blastdb_nhr 
   File blastdb_nin 
   File blastdb_nsq 
   String blastn_path 
   
   # Fingerprinting haplotype database
   File haplotype_map
   Boolean? fingerprint_tumor_normal

   ## Use as a last resort to increase the disk given to every task in case of ill behaving data
   Int? emergency_extra_disk
   # This is added to every task as padding, should increase if systematically you need more disk for every call
   Int disk_pad = 10 + select_first([emergency_extra_disk,0])
   #if set to true, use the strand bias filter on the raw annotation
   #if set to false (default), uses the strand bias filter on the duplex annotations 

   # Funcotator inputs
   Boolean filter_funcotations 
   File? data_sources_tar_gz
   String? funcotator_extra_args
   
   if(select_first([fingerprint_tumor_normal, true])) {
      call CrosscheckFingerprints {
         input:
            bloodbiopsydocker = bloodbiopsydocker,
            preemptible_attempts = preemptible_attempts,
            disk_pad = disk_pad,
            input_bam = tumor_bam,
            second_input_bam = normal_bam,
            haplotype_map = haplotype_map
      }
   }

   # Collect Sequencing Artifact Metrics after deduplication by start and stop
   # positions (but not including UMIs).
   call CollectSequencingArtifactMetrics as ConsensusArtifactMetricsTumor {
      input:
         bam_file = tumor_bam,
         bam_idx = tumor_bam_idx,
         basename = tumor_sample_name,
         reference = reference,
         reference_idx = reference_idx,
         bloodbiopsydocker = bloodbiopsydocker, 
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   if (defined(normal_bam)) {
      call CollectSequencingArtifactMetrics as ConsensusArtifactMetricsNormal {
         input:
            bam_file =  normal_bam,
            bam_idx =  normal_bam_idx,
            basename = normal_sample_name,
            reference = reference,
            reference_idx = reference_idx,
            bloodbiopsydocker = bloodbiopsydocker, 
            preemptible_attempts = preemptible_attempts,
            disk_pad = disk_pad
      }
   }

   call m2.Mutect2 as M2Duplex {
      input: 
         intervals = target_intervals,
         ref_fasta = reference,
         ref_fai = reference_idx,
         ref_dict = reference_dict,
         tumor_reads = tumor_bam, 
         tumor_reads_index = tumor_bam_idx, 
         normal_reads = normal_bam, 
         normal_reads_index = normal_bam_idx, 
         pon = pon, 
         pon_idx = pon_idx, 
         scatter_count = scatter_count, 
         gnomad = gnomad, 
         gnomad_idx = gnomad_idx,
         variants_for_contamination = variants_for_contamination, 
         variants_for_contamination_idx = variants_for_contamination_idx, 
         run_orientation_bias_mixture_model_filter = run_ob_filter, 
         m2_extra_args = m2_extra_args,
         m2_extra_filtering_args = m2_extra_filtering_args_or_default, 
         gatk_docker = bloodbiopsydocker,
         preemptible_attempts = preemptible_attempts, 
         compress_vcfs = compress_vcfs
   } 

   # Apply MRD mapping filter to duplex
   call RunMappingFilter {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         blastn_path = blastn_path,
         mapping_filter_python_script = mapping_filter_python_script,
         vcf_file = M2Duplex.filtered_vcf,
         vcf_idx = M2Duplex.filtered_vcf_idx,
         reference = reference,
         reference_idx = reference_idx,
         reference_dict = reference_dict, 
         blastdb_nhr = blastdb_nhr,
         blastdb_nin = blastdb_nin,
         blastdb_nsq = blastdb_nsq, 
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }


   #if not in runmapping, mark as filtered by mapping filter
   call VariantFiltration as VariantFiltration_Mapping {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         filter_name = "mapping_filter",
         duplex_vcf_file = M2Duplex.filtered_vcf,
         duplex_vcf_idx = M2Duplex.filtered_vcf_idx,
         filter_vcf_file = RunMappingFilter.map_filtered_vcf,
         filter_vcf_idx = RunMappingFilter.map_filtered_vcf_idx,
         reference_fasta = reference,
         reference_fasta_idx = reference_idx,
         reference_dict = reference_dict, 
         preemptible_attempts = preemptible_attempts, 
         filter_not_in_mask = true,
         disk_pad = disk_pad
   }
   
   call RemoveIrrelevantFilters {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         vcf_file = VariantFiltration_Mapping.output_vcf,
         vcf_idx = VariantFiltration_Mapping.output_vcf_idx,
         disk_pad = disk_pad
   }

   call SplitVCFs {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         vcf_file = RemoveIrrelevantFilters.irrelevant_filters_removed_vcf,
         vcf_idx = RemoveIrrelevantFilters.irrelevant_filters_removed_vcf_idx,
         reference_fasta = reference,
         reference_fasta_idx = reference_idx,
         reference_dict = reference_dict, 
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   }

   call VariantsToTable {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         vcf = SplitVCFs.output_snp_vcf, 
         vcf_idx = SplitVCFs.output_snp_vcf_idx,
         output_name = basename, 
         preemptible_attempts = preemptible_attempts,
         disk_pad = disk_pad
   } 

   call CalculateMeanAlleleFraction {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         output_table = VariantsToTable.output_table,
         preemptible_attempts  = preemptible_attempts,
         tumor_sample_name = tumor_sample_name,
         disk_pad = disk_pad
   }

   call m2.Funcotate as FuncotateMafSnps {
      input:
         gatk_docker = bloodbiopsydocker,
         ref_fasta = reference,
         ref_fai = reference_idx,
         ref_dict = reference_dict,
         input_vcf = SplitVCFs.output_snp_vcf,
         input_vcf_idx = SplitVCFs.output_snp_vcf_idx,
         case_id = tumor_sample_name,
         control_id = normal_sample_name,
         reference_version = reference_version,
         data_sources_tar_gz = data_sources_tar_gz,
         filter_funcotations = filter_funcotations,
         extra_args = funcotator_extra_args,
         use_gnomad = false, 
         output_format = "MAF",
         output_file_base_name = basename, 
         compress = true, 
         preemptible_attempts = preemptible_attempts
   }

   call m2.Funcotate as FuncotateMafIndels {
      input:
         gatk_docker = bloodbiopsydocker,
         ref_fasta = reference,
         ref_fai = reference_idx,
         ref_dict = reference_dict,
         input_vcf = SplitVCFs.output_indel_vcf,
         input_vcf_idx = SplitVCFs.output_indel_vcf_idx,
         case_id = tumor_sample_name,
         control_id = normal_sample_name,
         reference_version = reference_version,
         data_sources_tar_gz = data_sources_tar_gz,
         filter_funcotations = filter_funcotations,
         extra_args = funcotator_extra_args,
         use_gnomad = false, 
         output_format = "MAF",
         output_file_base_name = basename, 
         compress = true, 
         preemptible_attempts = preemptible_attempts
   }

   Array[String] input_files = select_all([SplitVCFs.output_snp_vcf, SplitVCFs.output_indel_vcf, tumor_bam, normal_bam])

   call GenerateIGVSession {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         preemptible_attempts = preemptible_attempts, 
         input_files = input_files,
         file_name =  basename, 
         reference_version = reference_version
   }
   
   output {

      File igv_session = GenerateIGVSession.igv_session

      File filtered_vcf = RemoveIrrelevantFilters.irrelevant_filters_removed_vcf
      File filtered_vcf_idx = RemoveIrrelevantFilters.irrelevant_filters_removed_vcf_idx

      File filtered_snp_vcf = SplitVCFs.output_snp_vcf
      File filtered_snp_vcf_idx = SplitVCFs.output_snp_vcf_idx

      File filtered_indel_vcf = SplitVCFs.output_indel_vcf
      File filtered_indel_vcf_idx = SplitVCFs.output_indel_vcf_idx

      Int n_passing_snps = SplitVCFs.passing_SNP
      Int n_filtered_snps = SplitVCFs.filtered_SNP
      
      Int n_passing_indels = SplitVCFs.passing_INDEL
      Int n_filtered_indels = SplitVCFs.filtered_INDEL

      File? contamination_table = M2Duplex.contamination_table

      File funcotated_snp_maf = FuncotateMafSnps.funcotated_output_file
      File funcotated_indel_maf = FuncotateMafIndels.funcotated_output_file

      File variant_table = VariantsToTable.output_table

      Float mean_af = CalculateMeanAlleleFraction.mean_af

   }
}


task CrosscheckFingerprints {
   String bloodbiopsydocker
   Int? preemptible_attempts
   Int disk_pad
   String input_bam
   String second_input_bam
   File haplotype_map
   Int? memory
   
   Int disk_size = 20 + disk_pad 
   Int mem = select_first([memory, 5])   
   
   command {
      set -e
      
      gatk CrosscheckFingerprints \
         --I ${input_bam} \
         --I ${second_input_bam} \
         --EXIT_CODE_WHEN_MISMATCH 1 \
         --CROSSCHECK_BY SAMPLE \
         --EXPECT_ALL_GROUPS_TO_MATCH true \
         --HAPLOTYPE_MAP ${haplotype_map} \
         --OUTPUT fingerprint_metrics.txt
    }
    runtime {
       docker: bloodbiopsydocker
       memory: mem + " GB"
       disks: "local-disk " + disk_size + " HDD"
       preemptible: select_first([preemptible_attempts, 2])
    }
    output {
       File fingerprint_metrics = "fingerprint_metrics.txt"
    }
}


# Collect Sequencing Artifact Metrics
# This is generally a good set of metrics to collect, but
# we will also need this later for the variant filtration
# step in FilterByOrientationBias
task CollectSequencingArtifactMetrics {
   String bloodbiopsydocker
   File bam_file
   File bam_idx
   File reference
   File reference_idx
   String? basename
   Int disk_pad
   Int? preemptible_attempts
   Int ref_size = ceil(size(reference, "GB") + size(reference_idx, "GB"))
   Int? memory
   Int disk_size = 500
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500

   command {
      set -e

      gatk --java-options "-Xmx4G" \
         CollectSequencingArtifactMetrics \
         --I ${bam_file} \
         --O "${basename}.gatk" \
         --R ${reference} \
         --VALIDATION_STRINGENCY LENIENT
   }
   runtime {
      docker: bloodbiopsydocker
      memory: mem + " GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_attempts, 0])
   }
   output {
      File pre_adapter_metrics = "${basename}.gatk.pre_adapter_detail_metrics"
   }
}


task RunMappingFilter {

   String bloodbiopsydocker
   String basename
   File vcf_file
   File vcf_idx
   File reference 
   File reference_idx 
   File reference_dict
   String mapping_filter_python_script
   String blastn_path
   File blastdb_nhr
   File blastdb_nin
   File blastdb_nsq
   Int? preemptible_attempts
   Int disk_pad
   Int? memory
   Int ref_size = ceil(size(reference, "GB") + size(reference_idx, "GB") + size(reference_dict, "GB"))
   Int blast_size = ceil(size(blastdb_nhr, "GB") + size(blastdb_nin, "GB") + size(blastdb_nsq, "GB"))
   Int disk_size = ceil(size(vcf_file, "GB") * 2) + ref_size + blast_size + disk_pad 
   Int mem = select_first([memory, 10])

   command {

      set -e

      python2.7 ${mapping_filter_python_script} \
         --vcf ${vcf_file} \
         --outfile ${basename}.filtered.vcf \
         --reference_fasta ${reference} \
         --blastn ${blastn_path}

      bgzip "${basename}.filtered.vcf"
      tabix "${basename}.filtered.vcf.gz"

   }
   runtime {
      docker: bloodbiopsydocker
      preemptible: select_first([preemptible_attempts, 0])
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
   }
   output {
      File map_filtered_vcf = "${basename}.filtered.vcf.gz"
      File map_filtered_vcf_idx = "${basename}.filtered.vcf.gz.tbi"
   }
}


task RemoveIrrelevantFilters {
   String bloodbiopsydocker
   String basename
   File vcf_file
   File vcf_idx
   Int? preemptible_attempts
   Int disk_pad
   Int? memory
   Int disk_size = ceil(size(vcf_file, "GB") * 2) + disk_pad 
   Int mem = select_first([memory, 10])

   command {
      set -e
      bcftools annotate -x FILTER/strand_bias -o ${basename}.filters_removed.vcf ${vcf_file}
      gatk IndexFeatureFile -I ${basename}.filters_removed.vcf
   }
   runtime {
      docker: bloodbiopsydocker
      preemptible: select_first([preemptible_attempts, 0])
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
   }
   output {
      File irrelevant_filters_removed_vcf = "${basename}.filters_removed.vcf"
      File irrelevant_filters_removed_vcf_idx = "${basename}.filters_removed.vcf.idx"
   }
}



# This task filters out anything not included in the filter_vcf_file
task VariantFiltration {

   String bloodbiopsydocker
   String basename
   File? gatk_override
   String filter_name
   File duplex_vcf_file
   File duplex_vcf_idx
   File filter_vcf_file
   File filter_vcf_idx
   File reference_fasta
   File reference_fasta_idx
   File reference_dict
   Int? preemptible_attempts
   Boolean filter_not_in_mask
   Int? memory
   Int disk_pad
   Int mem = select_first([memory, 5])
   Int ref_size = ceil(size(reference_fasta, "GB") + size(reference_fasta_idx, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(duplex_vcf_file, "GB") * 1.25) + ceil(size(filter_vcf_file, "GB") * 1.25) + ref_size + disk_pad 


   command {
      set -e

      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      gatk VariantFiltration \
         -R ${reference_fasta} \
         -V ${duplex_vcf_file} \
         -O ${basename}.${filter_name}.vcf.gz \
         --mask ${filter_vcf_file} \
         --filter-not-in-mask ${filter_not_in_mask} \
         --mask-name ${filter_name}
         
   }

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      preemptible: select_first([preemptible_attempts, 0])
   }

   output {
      File output_vcf = "${basename}.${filter_name}.vcf.gz"
      File output_vcf_idx = "${basename}.${filter_name}.vcf.gz.tbi"
   }

}


# select variants that PASS a given filter
task SelectStrictStrandBias {

   String bloodbiopsydocker
   File? gatk_override
   String basename
   File vcf_file
   File vcf_idx
   File reference_fasta
   File reference_fasta_idx
   File reference_dict
   Int? preemptible_attempts
   String extra_filter_arg
   Int? memory
   Int disk_pad
   Int mem = select_first([memory, 10])
   Int compute_mem = mem * 1000 - 500
   Int ref_size = ceil(size(reference_fasta, "GB") + size(reference_fasta_idx, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(vcf_file, "GB") * 2) + ref_size + disk_pad 


   command {

   set -e

   export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

   gatk --java-options -Xmx${compute_mem}m SelectVariants \
      -V ${vcf_file} \
      -O "${basename}.passing.vcf" \
      -R ${reference_fasta} \
      ${extra_filter_arg}

      
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      preemptible: select_first([preemptible_attempts, 0])
         
   }
   output {
      File output_vcf = "${basename}.passing.vcf"
      File output_vcf_idx = "${basename}.passing.vcf.idx"
   }
}

task SplitVCFs {

   String bloodbiopsydocker
   String basename
   File? gatk_override
   File vcf_file
   File vcf_idx
   File reference_fasta
   File reference_fasta_idx
   File reference_dict
   Int? preemptible_attempts
   Int? memory
   Int disk_pad
   Int ref_size = ceil(size(reference_fasta, "GB") + size(reference_fasta_idx, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(vcf_file, "GB") * 2) + ceil(size(vcf_idx)) + ref_size + disk_pad 
   Int mem = select_first([memory, 16])

   command <<<

      set -e

      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      gatk --java-options "-Xmx15G" SelectVariants \
         -V ${vcf_file} \
         -O "${basename}.snp.vcf.gz" \
         -R ${reference_fasta} \
         -select-type SNP

      #passing only 
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -O "${basename}.snp.passing.vcf.gz" \
         -R ${reference_fasta} \
         --exclude-filtered

      passing_SNP="$(gatk CountVariants -V ${basename}.snp.passing.vcf.gz | tail -1)"
      echo "$passing_SNP" > passing_snp.txt

      #filtered only 
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -XL "${basename}.snp.passing.vcf.gz" \
         -O "${basename}.snp.filtered.vcf.gz" \
         -R ${reference_fasta} 

      filtered_SNP="$(gatk CountVariants -V ${basename}.snp.filtered.vcf.gz | tail -1)"
      echo "$filtered_SNP" > filtered_snp.txt


      gatk --java-options "-Xmx15G" SelectVariants \
         -V ${vcf_file} \
         -O "${basename}.indel.vcf.gz" \
         -R ${reference_fasta} \
         -select-type INDEL

      #passing only 
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.indel.vcf.gz" \
         -O "${basename}.indel.passing.vcf.gz" \
         -R ${reference_fasta} \
         --exclude-filtered

      passing_INDEL="$(gatk CountVariants -V ${basename}.indel.passing.vcf.gz | tail -1)"
      echo "$passing_INDEL" > passing_indel.txt

      #filtered only 
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.indel.vcf.gz" \
         -XL "${basename}.indel.passing.vcf.gz" \
         -O "${basename}.indel.filtered.vcf.gz" \
         -R ${reference_fasta} 

      filtered_INDEL="$(gatk CountVariants -V ${basename}.indel.filtered.vcf.gz | tail -1)"
      echo "$filtered_INDEL" > filtered_indel.txt

   >>>
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      preemptible: select_first([preemptible_attempts, 0])
   }
   output {

      File output_snp_vcf = "${basename}.snp.vcf.gz"
      File output_snp_vcf_idx = "${basename}.snp.vcf.gz.tbi"

      Int passing_SNP = read_int("passing_snp.txt")
      Int filtered_SNP = read_int("filtered_snp.txt")

      File output_indel_vcf = "${basename}.indel.vcf.gz"
      File output_indel_vcf_idx = "${basename}.indel.vcf.gz.tbi"

      Int passing_INDEL = read_int("passing_indel.txt")
      Int filtered_INDEL = read_int("filtered_indel.txt")
      
   }   
}

#creates an IGV session
#given a list of IGV compatible files (as strings)
#reference is either "hg19" or "hg38"

task GenerateIGVSession {
   String bloodbiopsydocker
   Array[String] input_files
   String reference_version
   String file_name
   Array[String]? input_names
   Array[String] input_names_prefix = if defined(input_names) then prefix('-n ', select_first([input_names])) else []
   Int? preemptible_attempts
    
   command {
      bash /usr/writeIGV.sh ${reference_version} ${sep=" " input_files} ${sep=" " input_names_prefix}  > "${file_name}.xml"
   }
   runtime {
      docker: bloodbiopsydocker
      preemptible: select_first([preemptible_attempts, 0])
   }
   output {
      File igv_session = "${file_name}.xml"
   }
}


task VariantsToTable {
  
   String bloodbiopsydocker
   File vcf
   File vcf_idx
   String output_name
   File? gatk_override
   Int? memory
   Int mem = if defined(memory) then memory * 1000 else 3500
   Int compute_mem = mem - 500
   Int disk_pad
   Int disk_size = ceil(size(vcf, "GB") * 1.5 ) + ceil(size(vcf_idx, "GB")) + disk_pad 
   Int? preemptible_attempts 
  
   command {

      set -e

      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      gatk --java-options "-Xmx${compute_mem}m" VariantsToTable \
         -V ${vcf} \
         -F CHROM -F POS -F ID -F REF -F ALT -F QUAL \
         -F FILTER -ASGF AD -ASGF AF -GF NCount -ASGF SB \
         -SMA \
         -O "${output_name}.table"

   }

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + 500 + " HDD"
      memory: mem + " MB"
      preemptible: select_first([preemptible_attempts, 0])
   }
   output {
      File output_table = "${output_name}.table"
   }
}

task CalculateMeanAlleleFraction {

   String bloodbiopsydocker
   File output_table
   Int? preemptible_attempts 
   String tumor_sample_name
   Int? memory
   Int mem = select_first([memory, 5])
   Int disk_pad
   Int disk_size = ceil(size(output_table, "GB") * 1.5) + disk_pad

   command <<<

   python2.7 <<CODE

   import pandas as pd
   import numpy as np

   df = pd.read_csv("${output_table}", delim_whitespace=True)
   mean_af = round(df["${tumor_sample_name}"+".AF"].mean(),3)

   if np.isnan(mean_af): 
      print(0)
   else:
      print(mean_af)

   CODE
   >>>
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem +" GB"
      preemptible: select_first([preemptible_attempts, 0])
   }
   output {
      Float mean_af = read_float(stdout())
   }
}
