version 1.0

workflow ImputationPipeline {
  input {
    Int chunkLength = 25000000
    File multi_sample_vcf
    File multi_sample_vcf_index
    File samples_list = "gs://broad-dsde-methods-skwalker/polygenic_risk_scores/partners_subset_samples.args"
    Boolean perform_qc_steps
    File ref_dict = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict"
    String genetic_maps_eagle = "/genetic_map_hg19_withX.txt.gz"
    String output_callset_name = "broad_imputation"
    String path_to_reference_panel = "gs://broad-dsde-methods-skwalker/polygenic_risk_scores/eagle_reference_panels/"
    String path_to_m3vcf = "gs://broad-dsde-methods-skwalker/polygenic_risk_scores/minimac3_files/"
  }

  Array[String] samples = read_lines(samples_list)

  ## TODO: need to make sure each sample is valid 
 

  scatter (chrom in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]) { 
    call CalculateChromsomeLength {
      input:
        ref_dict = ref_dict,
        chrom = chrom
    }
    
    Float chunkLengthFloat = chunkLength
    Int num_chunks = ceil(CalculateChromsomeLength.chrom_length / chunkLengthFloat)

    scatter (i in range(num_chunks)) {

      call GenerateChunk {
        input:
          samples_list = samples_list,
          vcf = multi_sample_vcf,
          vcf_index = multi_sample_vcf_index,
          start = (i * chunkLength) + 1,
          end = if (CalculateChromsomeLength.chrom_length < ((i + 1) * chunkLength)) then CalculateChromsomeLength.chrom_length else ((i + 1) * chunkLength), 
          chrom = chrom,
          basename = "chrom_" + chrom + "_chunk_" + i
      }


       if (perform_qc_steps) {
        call QCSites {
          input:
            input_vcf = GenerateChunk.output_vcf,
            input_vcf_index = GenerateChunk.output_vcf_index,
            output_vcf_basename =  "chrom_" + chrom + "_chunk_" + i

          }
        }

      call CheckChunkValid {
        input: 
          vcf = select_first([QCSites.output_vcf,  GenerateChunk.output_vcf]),
          vcf_index = select_first([QCSites.output_vcf_index, GenerateChunk.output_vcf_index]),
          panel_vcf = path_to_reference_panel + "ALL.chr" + chrom + ".phase3_integrated.20130502.genotypes.vcf.gz",
          panel_vcf_index = path_to_reference_panel + "ALL.chr" + chrom + ".phase3_integrated.20130502.genotypes.vcf.gz.tbi"
      }

      if (CheckChunkValid.valid) {

      call PrePhaseVariantsEagle {
        input:
          dataset_bcf = CheckChunkValid.valid_chunk_bcf,
          dataset_bcf_index = CheckChunkValid.valid_chunk_bcf_index,
          reference_panel_bcf = path_to_reference_panel + "ALL.chr" + chrom + ".phase3_integrated.20130502.genotypes.bcf",
          reference_panel_bcf_index = path_to_reference_panel + "ALL.chr" + chrom + ".phase3_integrated.20130502.genotypes.bcf.csi",
          chrom = chrom,
          genetic_map_file = genetic_maps_eagle,
          start = (i * chunkLength) + 1,
          end = (i + 1) * chunkLength + 5000000 # they do an overlap of 5,000,000 bases in michigan pipeline
      }

        call minimac4 {
          input:
            ref_panel = path_to_m3vcf + chrom + ".1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz",
            phased_vcf = PrePhaseVariantsEagle.dataset_prephased_vcf,
            prefix = "chrom" + "_chunk_" + i +"_imputed",
            chrom = chrom,
            start = (i * chunkLength) + 1,
            end = (i + 1) * chunkLength
        }
      }
     }
     Array[File] chromsome_vcfs = select_all(minimac4.vcf)
  }

  Array[File] phased_vcfs = flatten(chromsome_vcfs)

  call MergeVCFs {
    input:
      input_vcfs = phased_vcfs,
      output_vcf_basename = output_callset_name
  }

  call UpdateHeader {
    input:
      vcf = MergeVCFs.output_vcf,
      vcf_index = MergeVCFs.output_vcf_index,
      ref_dict = ref_dict
  }

  call RemoveSymbolicAlleles {
    input:
      original_vcf = UpdateHeader.output_vcf,
      original_vcf_index = UpdateHeader.output_vcf_index
  }

  call SeparateMultiallelics {
    input:
      original_vcf = RemoveSymbolicAlleles.output_vcf,
      original_vcf_index = RemoveSymbolicAlleles.output_vcf_index
  }

  scatter (i in range(length(samples))) {

    call SplitSample {
      input: 
        sample = samples[i],
        vcf = SeparateMultiallelics.output_vcf,
        vcf_index = SeparateMultiallelics.output_vcf_index
    }
  }
  output {
    Array[File] imputed_output_vcf = SplitSample.output_gzipped_vcf
    Array[File] imputed_output_vcf_index = SplitSample.output_gzipped_vcf_index
  }
}

task RemoveDuplicates {
    input {
      File input_vcfs
      File input_vcf_indices
      String output_basename
    }
     Int disk_size = size([input_vcfs, input_vcf_indices], "GB")
  command <<<
    bcftools norm -d both ~{input_vcfs} | bgzip -c > ~{output_basename}.vcf.gz
    bcftools index -t ~{output_basename}.vcf.gz
  >>>
  runtime {
    docker: "farjoun/impute:0.0.3-1504715575"
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}


task CalculateChromsomeLength {
  input {
    File ref_dict
    Int chrom
  }
  command {
    grep -P "SN:~{chrom}\t" ~{ref_dict} | sed 's/.*LN://' | sed 's/\t.*//'
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    disks: "local-disk 100 HDD"
    memory: "2 GB"
  }
  output {
    Int chrom_length = read_int(stdout())
  }
}

task GenerateChunk {
  input {
    Int start
    Int end
    String chrom
    String basename
    String vcf
    String vcf_index
    File samples_list
    Int disk_size = 400 # not sure how big the disk size needs to be since it should be using
    # 2*size([vcf, vcf_index], "GB")
  }
  command {

    gatk SelectVariants -V ~{vcf} --select-type-to-include SNP --max-nocall-fraction 0.1 \
    --restrict-alleles-to BIALLELIC  -xl-select-type SYMBOLIC --select-type-to-exclude MIXED \
    --exclude-non-variants true -L ~{chrom}:~{start}-~{end} -O ~{basename}.vcf.gz \
    -sn ~{samples_list} --exclude-filtered true
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "8 GB"
  }

  parameter_meta {
    vcf: {
      description: "vcf",
      localization_optional: true
    } vcf_index: {
      description: "vcf index",
      localization_optional: true
    }
  }

  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi" 
  }
}  

task CheckChunkValid {
  input {
    File vcf
    File vcf_index
    File panel_vcf
    File panel_vcf_index
    Int disk_size = size([vcf, vcf_index, panel_vcf, panel_vcf_index], "GB")
  }
  command <<<

    ### TODO: make sure ref and alts match between array and panel or are a simple mismatch
 
    var_in_original=$(gatk CountVariants -V ~{vcf})
    var_in_reference=$(gatk CountVariants -V ~{vcf} -L ~{panel_vcf})

    echo ${var_in_reference} " * 2 - " ${var_in_original} "should be greater than 0 AND " ${var_in_reference} "should be greater than 3"
    if [ $(( ${var_in_reference} * 2 - ${var_in_original})) -gt 0 ] && [ ${var_in_reference} -gt 3 ]; then
      echo true > valid_file.txt
      bcftools convert -Ob ~{vcf} > valid_variants.bcf
      bcftools index -f valid_variants.bcf 
    else
      echo false > valid_file.txt
    fi
  >>>
  output {
    File? valid_chunk_bcf ="valid_variants.bcf"
    File? valid_chunk_bcf_index = "valid_variants.bcf.csi"
    Boolean valid = read_boolean("valid_file.txt")
  }
  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4 GB"
  }
}

task PrePhaseVariantsEagle {
  input {
    File? dataset_bcf
    File? dataset_bcf_index
    File reference_panel_bcf
    File reference_panel_bcf_index
    String chrom
    String genetic_map_file
    Int start
    Int end
  }  
  Int disk_size = 1.5 * size([dataset_bcf, reference_panel_bcf, dataset_bcf_index, reference_panel_bcf_index], "GB")
  command <<<
      /eagle  \
             --vcfTarget ~{dataset_bcf}  \
             --vcfRef ~{reference_panel_bcf} \
             --geneticMapFile ~{genetic_map_file} \
             --outPrefix pre_phased_~{chrom} \
             --vcfOutFormat z \
             --bpStart ~{start} --bpEnd ~{end} --allowRefAltSwap 
  >>>
  output {
    File dataset_prephased_vcf="pre_phased_~{chrom}.vcf.gz"
  }
  runtime {
    docker: "skwalker/imputation:test"
    memory: "32 GB"
    cpu: "8"
    disks: "local-disk " + disk_size + " HDD"
  }
}

task minimac4 {
  input {
    File ref_panel
    File phased_vcf
    String prefix
    String chrom
    Int start
    Int end
  }
  command <<<
    /Minimac4 --refHaps ~{ref_panel} --haps ~{phased_vcf} --start ~{start} --end ~{end} --window 500000 \
      --chr ~{chrom} --noPhoneHome --format GT,DS,GP --allTypedSites --prefix ~{prefix} --minRatio 0.00001 
  >>>
  output {
    File vcf = "~{prefix}.dose.vcf.gz"
    File info = "~{prefix}.info"
  }
  runtime {
    docker: "skwalker/imputation:test"
    memory: "4 GB"
    cpu: "1"
    disks: "local-disk 100 HDD"
  }
}

task MergeVCFs {
  input {
    Array[File] input_vcfs # these are all valid
    String output_vcf_basename
  }
  
  Int disk_size = 1.5*size(input_vcfs, "GB") 
  
  command <<<
    bcftools concat ~{input_vcfs} -Oz -o ~{output_vcf_basename}.vcf.gz
    bcftools index -t ~{output_vcf_basename}.vcf.gz

  >>>
  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "~{output_vcf_basename}.vcf.gz"
    File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}

task UpdateHeader {
  input {
    File vcf
    File vcf_index
    File ref_dict
    Int disk_size = 2*(size(vcf, "GB") + size(vcf_index, "GB"))
  }
  command <<<
    
    ## update the header of the merged vcf
    gatk UpdateVCFSequenceDictionary --source-dictionary ~{ref_dict} --output final_call.vcf.gz \
     --replace -V ~{vcf} --disable-sequence-dictionary-validation
  >>>
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "8 GiB"
  }
  output {
    File output_vcf = "final_call.vcf.gz"
    File output_vcf_index = "final_call.vcf.gz.tbi"   
  }
}

task RemoveSymbolicAlleles {
  input {
    File original_vcf
    File original_vcf_index 
    String output_basename = "allchroms.no_symbolic"
    Int disk_size = 2*(size(original_vcf, "GB") + size(original_vcf_index, "GB"))
  }
  command {
    gatk SelectVariants -V ~{original_vcf} -xl-select-type SYMBOLIC --select-type-to-exclude MIXED \
    --exclude-non-variants TRUE --remove-unused-alternates TRUE -O ~{output_basename}.vcf.gz
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4 GB"
  }
}

task SeparateMultiallelics {
  input {
    File original_vcf
    File original_vcf_index 
    String output_basename = "allchroms.no_multi_symbolic"
    Int disk_size =  2*(size(original_vcf, "GB") + size(original_vcf_index, "GB"))
  }
  command {
    bcftools norm -m - ~{original_vcf} | bgzip -c > ~{output_basename}.vcf.gz
    bcftools index -t ~{output_basename}.vcf.gz
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4 GB"
  }
}

task SplitSample {
  input {
    File vcf
    File vcf_index
    String sample
    Int disk_size = 2*(size(vcf, "GB") + size(vcf_index, "GB"))
  }
  command {
    gatk SelectVariants -V ~{vcf} -sn ~{sample} -O ~{sample}.vcf.gz
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "9 GB"
  }
  output {
    File output_gzipped_vcf = "~{sample}.vcf.gz"
    File output_gzipped_vcf_index = "~{sample}.vcf.gz.tbi"
  }
}

task QCSites {
  input {
    File input_vcf
    File input_vcf_index
    String output_vcf_basename
   }
    Int disk_size = 2*(size(input_vcf, "GB") + size(input_vcf_index, "GB"))

  command <<<
    # site missing rate < 5% ; hwe p > 1e-6
    vcftools --gzvcf ~{input_vcf}  --max-missing 0.05 --hwe 0.000001 --recode -c | bgzip -c > ~{output_vcf_basename}.vcf.gz
    bcftools index -t ~{output_vcf_basename}.vcf.gz 
  >>>

  runtime {
    docker: "skwalker/imputation:with_vcftools" # TODO: use a public one
    memory: "16 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "~{output_vcf_basename}.vcf.gz"
    File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}
