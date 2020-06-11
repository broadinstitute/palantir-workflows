# Hi this is Mark and Chris
workflow GenerateTruth {

  String ChrisIsCool
  String bloodbiopsydocker

  String base_name 
  String tumor_sample_name
  File duplex_tumor_bam
  File duplex_tumor_bam_index 
  String normal_sample_name 
  File duplex_normal_bam
  File duplex_normal_bam_index

  File reference
  File reference_index
  File reference_dict

  File target_intervals 
  File bait_intervals 


  call HaplotypeCaller as HaplotypeCallerTumor {
     input:
        reference = reference,
        reference_index = reference_index,
        reference_dict = reference_dict,
        target_intervals = target_intervals,
        bam_file = duplex_tumor_bam,
        bam_index = duplex_tumor_bam_index,
        tumor_normal = "tumor",
        base_name = tumor_sample_name,
        bloodbiopsydocker = bloodbiopsydocker
  }
  call HaplotypeCaller as HaplotypeCallerNormal {
     input:
        reference = reference,
        reference_index = reference_index,
        reference_dict = reference_dict,
        target_intervals = target_intervals,
        bam_file = duplex_normal_bam,
        bam_index = duplex_normal_bam_index,
        tumor_normal = "normal",
        base_name = normal_sample_name,
        bloodbiopsydocker = bloodbiopsydocker
  }

  call Mutect2Gatk4 as BlacklistedTumor {
     input:
        ref_fasta = reference,
        ref_fasta_index = reference_index,
        ref_dict = reference_dict,
        tumor_bam = duplex_tumor_bam,
        tumor_bam_index = duplex_tumor_bam_index,
        intervals = target_intervals,
        output_vcf_name = tumor_sample_name,
        tumor_normal = "tumor", 
        bloodbiopsydocker = bloodbiopsydocker
  }

  call Mutect2Gatk4 as BlacklistedNormal {
     input:
        ref_fasta = reference,
        ref_fasta_index = reference_index,
        ref_dict = reference_dict,
        tumor_bam = duplex_normal_bam,
        tumor_bam_index = duplex_normal_bam_index,
        intervals = target_intervals,
        output_vcf_name = normal_sample_name,
        tumor_normal = "normal", 
        bloodbiopsydocker = bloodbiopsydocker
  }

  call GenerateTruthVcfs {
     input:
        tumor_vcf = HaplotypeCallerTumor.output_vcf,
        tumor_vcf_index = HaplotypeCallerTumor.output_vcf_index,
        normal_vcf = HaplotypeCallerNormal.output_vcf,
        normal_vcf_index = HaplotypeCallerNormal.output_vcf_index,
        reference = reference,
        reference_index = reference_index,
        reference_dict = reference_dict,
        target_intervals = target_intervals,
        base_name = base_name,
        normal_sample_name = normal_sample_name,
        tumor_sample_name = tumor_sample_name,
        bloodbiopsydocker = bloodbiopsydocker
  }

  call MakeTruthSet {
     input: 
        bloodbiopsydocker = bloodbiopsydocker, 
        reference = reference, 
        reference_index = reference_index,
        reference_dict = reference_dict,
        input_vcf = GenerateTruthVcfs.confident_truth_vcf,
        input_vcf_index = GenerateTruthVcfs.confident_truth_vcf_index,
        confidence_interval = GenerateTruthVcfs.confident_regions,
        base_name = base_name,
        blacklisted_normal = BlacklistedNormal.output_vcf, 
        blacklisted_normal_index = BlacklistedNormal.output_vcf_index,
        blacklisted_tumor = BlacklistedTumor.output_vcf, 
        blacklisted_tumor_index = BlacklistedTumor.output_vcf_index,
  }
  output {
    File high_confidence_vcf = MakeTruthSet.high_confidence_vcf
    File high_confidence_vcf_index = MakeTruthSet.high_confidence_vcf_index
    File high_confidence_interval = MakeTruthSet.high_confidence_interval
    File high_confidence_bed = MakeTruthSet.high_confidence_bed
  }
}

# Call Mutect2 using GATK4.  At the moment GATK4 does not
# support tagged file inputs.  The way around this is to
# pass the tumor_sample_name, and normal_sample_name so that
# Mutect can distinguish the bams.
task Mutect2Gatk4 {
   String bloodbiopsydocker
   File ref_fasta
   File ref_fasta_index
   File ref_dict
   File tumor_bam
   File tumor_bam_index
   File intervals
   String output_vcf_name
   String tumor_normal

   Int disk_pad = 20
   Int disk_size = ceil(size(ref_fasta, "GB") + 
                   size(ref_fasta_index, "GB") +
                   size(ref_dict, "GB") + 
                   size(tumor_bam, "GB") +
                   size(tumor_bam_index, "GB") +
                   size(intervals, "GB")) * 2 + 
                   disk_pad

   command {
      gatk Mutect2 \
         -R ${ref_fasta} \
         ${"-I " + tumor_bam} \
         ${"-L " + intervals} \
         --dont-use-soft-clipped-bases \
         -O "${output_vcf_name}.${tumor_normal}.vcf.gz" \
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: "16 GB"
      preemptible: 0
   }
   output {
      File output_vcf = "${output_vcf_name}.${tumor_normal}.vcf.gz"
      File output_vcf_index = "${output_vcf_name}.${tumor_normal}.vcf.gz.tbi"
   }
}

task HaplotypeCaller {
   String bloodbiopsydocker
   File bam_file
   File bam_index
   File target_intervals
   File reference
   File reference_index
   File reference_dict
   String base_name
   String tumor_normal

   Int disk_pad = 20
   Int disk_size = ceil(size(reference, "GB") + 
                   size(reference_index, "GB") +
                   size(reference_dict, "GB") + 
                   size(bam_file, "GB") +
                   size(bam_index, "GB") +
                   size(target_intervals, "GB")) * 2 +  
                   disk_pad
   
   command <<<

      set -e

      gatk HaplotypeCaller \
         -R ${reference} \
         -I ${bam_file} \
         -L ${target_intervals} \
         -O ${base_name}.${tumor_normal}.vcf.gz
   >>>
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: "16 GB"
      preemptible: 0
   }  
   output {
      File output_vcf = "${base_name}.${tumor_normal}.vcf.gz"
      File output_vcf_index = "${base_name}.${tumor_normal}.vcf.gz.tbi"
   }  
}


# Generate Truth sets from pure "tumor" and normal vcfs.
# This creates a truthset that only includes het sites from
# the spiked in "tumor".
task GenerateTruthVcfs {

   String bloodbiopsydocker
   File tumor_vcf
   File tumor_vcf_index
   File normal_vcf
   File normal_vcf_index
   File reference
   File reference_index
   File reference_dict
   File target_intervals
   String base_name
   String normal_sample_name
   String tumor_sample_name

   Int disk_pad = 20
   Int disk_size = ceil(size(reference, "GB") +
                   size(reference_index, "GB") +
                   size(reference_dict, "GB") +
                   size(tumor_vcf, "GB") +
                   size(tumor_vcf_index, "GB") +
                   size(normal_vcf, "GB") +
                   size(normal_vcf_index, "GB") +
                   size(target_intervals, "GB")) * 2 +
                   disk_pad
    
   command <<<
      set -e

      # tumor_vcf is the vcf of the spike in, control_vcf is the vcf called on the sample being spiked into
      # These vcfs are generated using UnifiedGenotyper
      java -jar /usr/GenomeAnalysisTK.jar -T CombineVariants -R ${reference} -V ${tumor_vcf} -V ${normal_vcf} -o tumorandnormal.vcf.gz

      # TumorTruth.vcf.gz contains all the heterozygous variants in the tumor that are not in the control sample (either het or homvar).
      gatk SelectVariants -R ${reference} -V tumorandnormal.vcf.gz -select '!vc.getGenotype("${normal_sample_name}").isHomVar() && !vc.getGenotype("${normal_sample_name}").isHet() && (vc.getGenotype("${tumor_sample_name}").isHet())' -O TumorTruth.vcf.gz -L ${target_intervals}
      gatk SelectVariants -R ${reference} -V tumorandnormal.vcf.gz -select '!vc.getGenotype("${normal_sample_name}").isHomVar() && !vc.getGenotype("${normal_sample_name}").isHet() && (vc.getGenotype("${tumor_sample_name}").isHomVar())' -O TumorTruthHomVar.vcf.gz -L ${target_intervals}

      gatk SelectVariants -R ${reference} -V TumorTruth.vcf.gz -sn "${normal_sample_name}" -O HetTruthSites.vcf.gz
      gatk SelectVariants -R ${reference} -V ${tumor_vcf} -L HetTruthSites.vcf.gz -O truth.vcf.gz 

      # Generate Region of Confidence
      gatk IntervalListTools --I ${normal_vcf} --I ${tumor_vcf} --PADDING 50 --O TumorAndControl.padded.interval_list --ACTION UNION

      grep ^@ TumorAndControl.padded.interval_list | grep -v PG > nonConfidentRegions.interval_list
      grep interval.*interval TumorAndControl.padded.interval_list >> nonConfidentRegions.interval_list

      # Remove HomVar sites from confident region
      gatk IntervalListTools --I nonConfidentRegions.interval_list --I TumorTruthHomVar.vcf.gz --ACTION UNION --PADDING 15 --O nonConfidentRegions.interval_list 

      # Subtract non confident regions from target intervals
      gatk IntervalListTools --I ${target_intervals} --SI nonConfidentRegions.interval_list --O ${base_name}.ConfidentRegions.interval_list --ACTION SUBTRACT
      gatk IntervalListToBed --I ${base_name}.ConfidentRegions.interval_list --O ${base_name}.ConfidentRegions.bed

      # Subset the truth.vcf.gz to the confident regions to make confident.truth.vcf.gz
      gatk SelectVariants -V truth.vcf.gz -L ${base_name}.ConfidentRegions.interval_list -O ${base_name}.confident.truth.vcf.gz -R ${reference}

   >>>
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: "4 GB"
      preemptible: 0
   }
   output {
      File confident_regions = "${base_name}.ConfidentRegions.interval_list"
      File confident_regions_bed = "${base_name}.ConfidentRegions.bed"
      File truth_vcf = "truth.vcf.gz"
      File truth_vcf_tbi = "truth.vcf.gz.tbi"
      File confident_truth_vcf = "${base_name}.confident.truth.vcf.gz"
      File confident_truth_vcf_index = "${base_name}.confident.truth.vcf.gz.tbi"
   }
}

task MakeTruthSet {

   String bloodbiopsydocker
   File reference
   File reference_index
   File reference_dict

   File input_vcf
   File input_vcf_index
   File confidence_interval
   String base_name

   File blacklisted_normal
   File blacklisted_normal_index
   File blacklisted_tumor
   File blacklisted_tumor_index

   Int disk_pad = 20
   Int disk_size = ceil(size(reference, "GB") +
                   size(reference_index, "GB") +
                   size(reference_dict, "GB") +
                   size(input_vcf, "GB") +
                   size(input_vcf_index, "GB") +
                   size(confidence_interval, "GB") +
                   size(blacklisted_normal, "GB") +
                   size(blacklisted_normal_index, "GB") +
                   size(blacklisted_tumor, "GB") +
                   size(blacklisted_tumor_index, "GB") * 2) +
                   disk_pad

   command <<<

      set -e

      # Merge all the blacklisted sites
      java -jar /usr/GenomeAnalysisTK.jar -T CombineVariants -R ${reference} -V ${blacklisted_normal} -V ${blacklisted_tumor} -genotypeMergeOptions UNIQUIFY -o blacklist.vcf.gz
      
      # Find the sites that were called by mutect but not HC
      gatk SelectVariants -R ${reference} -V blacklist.vcf.gz -XL ${input_vcf} -O blacklist.only.vcf.gz

      # Subtract out the blacklisted sites
      gatk IntervalListTools --ACTION SUBTRACT --I ${confidence_interval} --SI blacklist.only.vcf.gz --PADDING 0 --O ${base_name}.confident.interval_list

      gatk IntervalListToBed --I ${confidence_interval} --O ${base_name}.confident.bed

      # Identify sites below 30x
      gatk SelectVariants \
         -R ${reference} \
         -V ${input_vcf}  \
         -OVI true \
         -select "DP > 30" \
         -O ${base_name}.above30x.confident.truth.vcf.gz

      # Subtract out the truth sites that are below 30 depth 
      gatk SelectVariants \
         -R ${reference} \
         -V ${input_vcf}  \
         -OVI true \
         -XL ${base_name}.above30x.confident.truth.vcf.gz \
         -O ${base_name}.below30x.confident.truth.vcf.gz

      # Subtract out the low coverage truth from confident interval list
      gatk IntervalListTools --ACTION SUBTRACT --I ${base_name}.confident.interval_list --SI ${base_name}.below30x.confident.truth.vcf.gz --PADDING 0 --O ${base_name}.confident.interval_list
      gatk IntervalListToBed --I ${base_name}.confident.interval_list --O ${base_name}.confident.bed
   >>>

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: "4 GB"
      preemptible: 0
   }
   output {

      File high_confidence_vcf = "${base_name}.above30x.confident.truth.vcf.gz"
      File high_confidence_vcf_index = "${base_name}.above30x.confident.truth.vcf.gz.tbi"
      File high_confidence_interval = "${base_name}.confident.interval_list"
      File high_confidence_bed = "${base_name}.confident.bed"
   }
}

