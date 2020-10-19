version 1.0

workflow ImputationPipeline {
  input {
    Int chunkLength = 25000000

    # You can either input a multisample VCF or an array of single sample VCFs
    # The pipeline will just merge the single sample VCFs into one multisample VCF
    # and then impute the multisample VCF
    # If you want to run a single sample VCF, set the multi_sample_vcf input to the
    # single sample VCF
    File? multi_sample_vcf
    File? multi_sample_vcf_index
    Array[File]? single_sample_vcfs
    Array[File]? single_sample_vcf_indices
    
    Boolean perform_extra_qc_steps # these are optional additional extra QC steps from Amit's group that should only be 
    # run for large sample sets, especially a diverse set of samples (it's further limiting called at sites to 95% and by HWE)
    File ref_dict = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict" # for reheadering / adding contig lengths in the header of the ouptut VCF
    String genetic_maps_eagle = "/genetic_map_hg19_withX.txt.gz" # this is for Eagle, it is in the docker image 
    String output_callset_name = "broad_imputation" # the output callset name
    String path_to_reference_panel = "gs://fc-6413177b-e99c-4476-b085-3da80d320081/eagle_panels/" # from the "Imputation and Polygenic Risk Score Files" workspace in Terra
    String path_to_m3vcf = "gs://fc-6413177b-e99c-4476-b085-3da80d320081/minimac3_files/" # from the same workspace ^
    Int par1_end = 2699620
    Int par2_start = 154931044
  }

  if (defined(single_sample_vcfs)) {
    call MergeSingleSampleVcfs {
      input:
        input_vcfs = select_first([single_sample_vcfs, []]),
        input_vcf_indices = select_first([single_sample_vcf_indices, []]),
        output_vcf_basename = "merged_input_samples"
    }
  }

  File vcf_to_impute = select_first([multi_sample_vcf, MergeSingleSampleVcfs.output_vcf])
  File vcf_index_to_impute = select_first([multi_sample_vcf_index, MergeSingleSampleVcfs.output_vcf_index])

  Float chunkLengthFloat = chunkLength

  #impute chromosomeX
  #first split into non-par, par1, par2
  call SplitX {
  	input:
  		input_vcf = vcf_to_impute,
        par1_interval = "X:1-" + par1_end,
        par2_interval = "X:" + par2_start + "-" + CalculateChromosomeLengthX.chrom_length,
        non_par_interval = "X:" + (par1_end + 1) + "-" + (par2_start - 1)
  }

  #par1 and par2 are short, most of X is non_par region
  call CalculateChromosomeLength as CalculateChromosomeLengthX {
  	input:
  		ref_dict = ref_dict,
  		chrom = "X"
  }

  #par1 and par2 can be treated the same as the autosome
  scatter (X_par in zip(zip([SplitX.par1_vcf, SplitX.par2_vcf], [SplitX.par1_vcf_index, SplitX.par2_vcf_index]), [1,2])) {
  	File vcf = X_par.left.left
  	File vcf_index = X_par.left.right
  	Int iPAR = X_par.right

  	call GenerateChunk as GenerateChunkPAR {
		input:
		  vcf = vcf,
		  vcf_index = vcf_index,
		  start = if iPAR == 1 then 1 else par2_start,
		  end = if iPAR ==1 then par1_end else CalculateChromosomeLengthX.chrom_length,
		  chrom = "X",
		  basename = "chrom_X_par_" + iPAR
	  }

	  if (perform_extra_qc_steps) {
		  call OptionalQCSites as OptionalQCSitesPAR {
			input:
			  input_vcf = GenerateChunkPAR.output_vcf,
			  input_vcf_index = GenerateChunkPAR.output_vcf_index,
			  output_vcf_basename =  "chrom_X_par_" + iPAR
			}
		  }

		call CheckChunkValid as CheckChunkValidPAR{
		  input:
			vcf = select_first([OptionalQCSitesPAR.output_vcf,  GenerateChunkPAR.output_vcf]),
			vcf_index = select_first([OptionalQCSitesPAR.output_vcf_index, GenerateChunkPAR.output_vcf_index]),
			panel_vcf = path_to_reference_panel + "ALL.chrX.phase3_integrated.20130502.genotypes.vcf.gz",
			panel_vcf_index = path_to_reference_panel + "ALL.chrX.phase3_integrated.20130502.genotypes.vcf.gz.tbi"
		}

		if (CheckChunkValidPAR.valid) {

		call PrePhaseVariantsEagle as PrePhaseVariantsEaglePAR {
		  input:
			dataset_bcf = CheckChunkValidPAR.valid_chunk_bcf,
			dataset_bcf_index = CheckChunkValidPAR.valid_chunk_bcf_index,
			reference_panel_bcf = path_to_reference_panel + "ALL.chrX.phase3_integrated.20130502.genotypes.bcf",
			reference_panel_bcf_index = path_to_reference_panel + "ALL.chrX.phase3_integrated.20130502.genotypes.bcf.csi",
			chrom = "X",
			genetic_map_file = genetic_maps_eagle,
			start = if iPAR == 1 then 1 else par2_start,
			end = if iPAR ==1 then par1_end else CalculateChromosomeLengthX.chrom_length # they do an overlap of 5,000,000 bases in michigan pipeline
		}

		  call minimac4 as minimac4PAR{
			input:
			  ref_panel = path_to_m3vcf + "X.Pseudo.Auto.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz",
			  phased_vcf = PrePhaseVariantsEaglePAR.dataset_prephased_vcf,
			  prefix = "chrom_X_par_" + iPAR +"_imputed",
			  chrom = "X",
			  start = if iPAR == 1 then 1 else par2_start,
			  end = if iPAR ==1 then par1_end else CalculateChromosomeLengthX.chrom_length
		  }

		  call UpdateHeader as UpdateHeaderPAR {
			input:
			  vcf = minimac4PAR.vcf,
			  vcf_index = minimac4PAR.vcf_index,
			  ref_dict = ref_dict,
			  basename = "chrom_X_par_" + iPAR +"_imputed",
		  }

		  call SeparateMultiallelics as SeparateMultiallelicsPAR {
			input:
			  original_vcf = UpdateHeaderPAR.output_vcf,
			  original_vcf_index = UpdateHeaderPAR.output_vcf_index,
			  output_basename = "chrom_X_par_" + iPAR +"_imputed",
		  }

		  call RemoveSymbolicAlleles as RemoveSymbolicAllelesPAR {
			  input:
				  original_vcf = SeparateMultiallelicsPAR.output_vcf,
				  original_vcf_index = SeparateMultiallelicsPAR.output_vcf_index,
				  output_basename = "chrom_X_par_" + iPAR +"_imputed",
		  }

		  call SortIds as SortIdsPAR {
			input:
			  vcf = RemoveSymbolicAllelesPAR.output_vcf,
			  basename = "chrom_X_par_" + iPAR +"_imputed"
		  }
		}
  	}

	# for non-par region, we must impute sex, and turn males in haploid
	call ImputeSexAndMakeMalesHaploid {
		input:
			vcf = SplitX.non_par_vcf,
			basename = "non_par.X"
	}

  Int num_chunksX = ceil(CalculateChromosomeLengthX.chrom_length / chunkLengthFloat)



  scatter (iX in range(num_chunksX)) {
  	call GenerateChunk as GenerateChunkX {
		input:
		  vcf = ImputeSexAndMakeMalesHaploid.vcf_out,
		  vcf_index = ImputeSexAndMakeMalesHaploid.vcf_out_index,
		  start = (iX * chunkLength) + 1,
		  end = if (CalculateChromosomeLengthX.chrom_length < ((iX + 1) * chunkLength)) then CalculateChromosomeLengthX.chrom_length else ((iX + 1) * chunkLength),
		  chrom = "X",
		  basename = "chrom_X_chunk_" + iX
	  }

	  if (perform_extra_qc_steps) {
		  call OptionalQCSites as OptionalQCSitesX {
			input:
			  input_vcf = GenerateChunkX.output_vcf,
			  input_vcf_index = GenerateChunkX.output_vcf_index,
			  output_vcf_basename =  "chrom_X_chunk_" + iX
			}
		  }

		call CheckChunkValid as CheckChunkValidX {
		  input:
			vcf = select_first([OptionalQCSitesX.output_vcf,  GenerateChunkX.output_vcf]),
			vcf_index = select_first([OptionalQCSitesX.output_vcf_index, GenerateChunkX.output_vcf_index]),
			panel_vcf = path_to_reference_panel + "ALL.chrX.phase3_integrated.20130502.genotypes.vcf.gz",
			panel_vcf_index = path_to_reference_panel + "ALL.chrX.phase3_integrated.20130502.genotypes.vcf.gz.tbi"
		}

		if (CheckChunkValidX.valid) {

		call PrePhaseVariantsEagle as PrePhaseVariantsEagleX {
		  input:
			dataset_bcf = CheckChunkValidX.valid_chunk_bcf,
			dataset_bcf_index = CheckChunkValidX.valid_chunk_bcf_index,
			reference_panel_bcf = path_to_reference_panel + "ALL.chrX.phase3_integrated.20130502.genotypes.bcf",
			reference_panel_bcf_index = path_to_reference_panel + "ALL.chrX.phase3_integrated.20130502.genotypes.bcf.csi",
			chrom = "X",
			genetic_map_file = genetic_maps_eagle,
			start = (iX * chunkLength) + 1,
			end = (iX + 1) * chunkLength + 5000000 # they do an overlap of 5,000,000 bases in michigan pipeline
		}

		  call minimac4 as minimac4X {
			input:
			  ref_panel = path_to_m3vcf + "X.Non.Pseudo.Auto.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz",
			  phased_vcf = PrePhaseVariantsEagleX.dataset_prephased_vcf,
			  prefix = "chromX" + "_chunk_" + iX +"_imputed",
			  chrom = "X",
			  start = (iX * chunkLength) + 1,
			  end = (iX + 1) * chunkLength
		  }

		  call UpdateHeader as UpdateHeaderX {
			input:
			  vcf = minimac4X.vcf,
			  vcf_index = minimac4X.vcf_index,
			  ref_dict = ref_dict,
			  basename = "chromX" + "_chunk_" + iX +"_imputed",
		  }

		  call SeparateMultiallelics as SeparateMultiallelicsX {
			input:
			  original_vcf = UpdateHeaderX.output_vcf,
			  original_vcf_index = UpdateHeaderX.output_vcf_index,
			  output_basename = "chromX_chunk_" + iX +"_imputed",
		  }

		  call RemoveSymbolicAlleles as RemoveSymbolicAllelesX {
			  input:
				  original_vcf = SeparateMultiallelicsX.output_vcf,
				  original_vcf_index = SeparateMultiallelicsX.output_vcf_index,
				  output_basename = "chrom_X_chunk_" + iX +"_imputed",
		  }

		  call SortIds as SortIdsX {
			input:
			  vcf = RemoveSymbolicAllelesX.output_vcf,
			  basename = "chrom_X_chunk_" + iX +"_imputed"
		  }
		}
	}


  scatter (chrom in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]) { 
    call CalculateChromosomeLength {
      input:
        ref_dict = ref_dict,
        chrom = chrom
    }
    

    Int num_chunks = ceil(CalculateChromosomeLength.chrom_length / chunkLengthFloat)

    scatter (i in range(num_chunks)) {

      call GenerateChunk {
        input:
          vcf = vcf_to_impute,
          vcf_index = vcf_index_to_impute,
          start = (i * chunkLength) + 1,
          end = if (CalculateChromosomeLength.chrom_length < ((i + 1) * chunkLength)) then CalculateChromosomeLength.chrom_length else ((i + 1) * chunkLength),
          chrom = chrom,
          basename = "chrom_" + chrom + "_chunk_" + i
      }

       if (perform_extra_qc_steps) {
        call OptionalQCSites {
          input:
            input_vcf = GenerateChunk.output_vcf,
            input_vcf_index = GenerateChunk.output_vcf_index,
            output_vcf_basename =  "chrom_" + chrom + "_chunk_" + i
          }
        }

      call CheckChunkValid {
        input: 
          vcf = select_first([OptionalQCSites.output_vcf,  GenerateChunk.output_vcf]),
          vcf_index = select_first([OptionalQCSites.output_vcf_index, GenerateChunk.output_vcf_index]),
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
        
        call UpdateHeader {
          input:
            vcf = minimac4.vcf,
            vcf_index = minimac4.vcf_index,
            ref_dict = ref_dict,
            basename = "chrom" + "_chunk_" + i +"_imputed",
        }
        
        call SeparateMultiallelics {
          input:
            original_vcf = UpdateHeader.output_vcf,
            original_vcf_index = UpdateHeader.output_vcf_index,
            output_basename = "chrom" + chrom + "_chunk_" + i +"_imputed",
        }

        call RemoveSymbolicAlleles {
            input:
                original_vcf = SeparateMultiallelics.output_vcf,
                original_vcf_index = SeparateMultiallelics.output_vcf_index,
                output_basename = "chrom" + chrom + "_chunk_" + i +"_imputed",
        }
        
        call SortIds {
          input:
            vcf = RemoveSymbolicAlleles.output_vcf,
            basename = "chrom" + chrom + "_chunk_" + i +"_imputed"
        }
      }
    }
    Array[File] chromosome_vcfs = select_all(SortIds.output_vcf)
    Array[File] chromosome_vcf_indices = select_all(SortIds.output_vcf_index)
  }

  Array[File] phased_vcfs_autosome = flatten(chromosome_vcfs)
  Array[File] phased_vcf_indices_autosome = flatten(chromosome_vcf_indices)

  Array[File] par_vcfs = select_all(SortIdsPAR.output_vcf)
  Array[File] par_vcf_indices = select_all(SortIdsPAR.output_vcf_index)

  Array[File] non_par_vcfs = select_all(SortIdsX.output_vcf)
  Array[File] non_par_vcf_indices = select_all(SortIdsX.output_vcf_index)

  Array[File] phased_vcfs = flatten([phased_vcfs_autosome, [par_vcfs[0]], non_par_vcfs, [par_vcfs[1]]])
  Array[File] phased_vcf_indices = flatten([phased_vcf_indices_autosome, [par_vcf_indices[0]], non_par_vcf_indices, [par_vcf_indices[1]]])

  call GatherVcfs {
    input: 
      input_vcfs = phased_vcfs,
      input_vcf_indices = phased_vcf_indices,
      output_vcf_basename = output_callset_name
  }

  output {
    File imputed_multisample_vcf = GatherVcfs.output_vcf
    File imputed_multisample_vcf_index = GatherVcfs.output_vcf_index
  }
}

task SplitX {
	input {
		File input_vcf
		String par1_interval
		String par2_interval
		String non_par_interval
	}

	Int disk_size = ceil(size(input_vcf, "GB")) + 50

	parameter_meta {
		input_vcf : {
			localization_optional : true
		}
	}

	command <<<
		set -xeuo pipefail
		gatk SelectVariants -V ~{input_vcf} -L ~{par1_interval} -O par1.X.vcf.gz
		gatk SelectVariants -V ~{input_vcf} -L ~{par2_interval} -O par2.X.vcf.gz
		gatk SelectVariants -V ~{input_vcf} -L ~{non_par_interval} -O non.par.X.vcf.gz
	>>>

	runtime {
		docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
		disks: "local-disk " + disk_size + "  HDD"
		memory: "16 GB"
	}

	output {
		File par1_vcf = "par1.X.vcf.gz"
		File par1_vcf_index = "par1.X.vcf.gz.tbi"
		File par2_vcf = "par2.X.vcf.gz"
		File par2_vcf_index = "par2.X.vcf.gz.tbi"
		File non_par_vcf = "non.par.X.vcf.gz"
		File non_par_vcf_index = "non.par.X.vcf.gz.tbi"
	}
}

task ImputeSexAndMakeMalesHaploid {
	input {
		File vcf
		String basename
	}

	Int disk_size = 4 * ceil(size(vcf, "GB")) + 50

	command <<<
		plink --vcf ~{vcf} --impute-sex --const-fid --make-bed --out sex_imputed
		plink2 --bfile sex_imputed --export vcf-4.2 bgz --out ~{basename}.sex_imputed
		bcftools index -t ~{basename}.sex_imputed.vcf.gz
	>>>

	runtime {
		docker: "quay.io/ckachuli/plink@sha256:b5ac6b3b4d28e3fd6f9ca09a8d759a0229899ae9626d0b869d3b3f111e83255f"
		disks: "local-disk " + disk_size + " HDD"
		memory: "16 GB"
  	}

  	output {
  		File vcf_out = "~{basename}.sex_imputed.vcf.gz"
  		File vcf_out_index = "~{basename}.sex_imputed.vcf.gz.tbi"
  	}
}


task RemoveDuplicates {
    input {
      File input_vcfs
      File input_vcf_indices
      String output_basename
    }
     Int disk_size = ceil(size([input_vcfs, input_vcf_indices], "GB")) + 50
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

task CalculateChromosomeLength {
  input {
    File ref_dict
    String chrom
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
    Int disk_size = 400 # not sure how big the disk size needs to be since we aren't downloading the entire VCF here 
  }
  command {
    gatk SelectVariants -V ~{vcf} --select-type-to-include SNP --max-nocall-fraction 0.1 -xl-select-type SYMBOLIC \
     --select-type-to-exclude MIXED --restrict-alleles-to BIALLELIC -L ~{chrom}:~{start}-~{end} -O ~{basename}.vcf.gz --exclude-filtered true
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
    Int disk_size =ceil(2*size([vcf, vcf_index, panel_vcf, panel_vcf_index], "GB"))
    File gatk_jar =  "gs://fc-6413177b-e99c-4476-b085-3da80d320081/gatk_jar_for_forcing_genotypegvcf_to_run.jar" # you can use any public gatk jar on the cloud (but I dont think we publish those anywhere? )
  }
  command <<<
 
    var_in_original=$(java -jar ~{gatk_jar} CountVariants -V ~{vcf}  | sed 's/Tool returned://')
    var_in_reference=$(java -jar ~{gatk_jar}  CountVariants -V ~{vcf} -L ~{panel_vcf}  | sed 's/Tool returned://')

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
    docker: "farjoun/impute:0.0.4-1506086533" # need to use this one because you need to be able to run java and bcftools (bcftools one doesn't let you)
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
  Int disk_size = ceil(3 * size([dataset_bcf, reference_panel_bcf, dataset_bcf_index, reference_panel_bcf_index], "GB"))
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
    docker: "skwalker/imputation:test" # this has the exact version of minimac and eagle we want to perfectly match Michigan Server
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
    if [ ! -f ~{prefix}.dose.vcf.gz.tbi ]
    then
      bcftools index -t ~{prefix}.dose.vcf.gz
    fi
  >>>
  output {
    File vcf = "~{prefix}.dose.vcf.gz"
    File vcf_index = "~{prefix}.dose.vcf.gz.tbi"
    File info = "~{prefix}.info"
  }
  runtime {
    docker: "skwalker/imputation:test" # this has the exact version of minimac and eagle we want to perfectly match Michigan Server
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
  
  Int disk_size = ceil(3*size(input_vcfs, "GB"))
  
  command <<<
    bcftools concat ~{sep=' ' input_vcfs} -Oz -o ~{output_vcf_basename}.vcf.gz
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

task GatherVcfs {
  input {
    Array[File] input_vcfs # these are all valid
    Array[File] input_vcf_indices # these are all valid    
    String output_vcf_basename
  }
  
  Int disk_size = ceil(3*size(input_vcfs, "GB"))
  
  command <<<
    gatk GatherVcfs -I ~{sep=' -I=' input_vcfs} -O ~{output_vcf_basename}.vcf.gz # I don't think this creates an output index file
    if [ ! -f ~{output_vcf_basename}.vcf.gz.tbi ]
    then
      gatk IndexFeatureFile -F ~{output_vcf_basename}.vcf.gz
    fi
  >>>
   runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "16 GB"
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
    String basename
    Int disk_size = ceil(4*(size(vcf, "GB") + size(vcf_index, "GB"))) + 20 
  }
  command <<<
    
    ## update the header of the merged vcf
    gatk UpdateVCFSequenceDictionary --source-dictionary ~{ref_dict} --output ~{basename}.vcf.gz \
     --replace -V ~{vcf} --disable-sequence-dictionary-validation
  >>>
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "8 GiB"
  }
  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"   
  }
}

task RemoveSymbolicAlleles {
  input {
    File original_vcf
    File original_vcf_index 
    String output_basename
    Int disk_size = ceil(3*(size(original_vcf, "GB") + size(original_vcf_index, "GB")))
  }
  command {
    gatk SelectVariants -V ~{original_vcf} -xl-select-type SYMBOLIC -O ~{output_basename}.vcf.gz
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
    String output_basename
    Int disk_size =  ceil(2*(size(original_vcf, "GB") + size(original_vcf_index, "GB")))
  }
  command {
    bcftools norm -m - ~{original_vcf} -Oz -o ~{output_basename}.vcf.gz
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
    Int disk_size = ceil(2*(size(vcf, "GB") + size(vcf_index, "GB")))
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

task OptionalQCSites {
  input {
    File input_vcf
    File input_vcf_index
    String output_vcf_basename
   }
    Int disk_size = ceil(2*(size(input_vcf, "GB") + size(input_vcf_index, "GB")))

  command <<<
    # site missing rate < 5% ; hwe p > 1e-6
    vcftools --gzvcf ~{input_vcf}  --max-missing 0.05 --hwe 0.000001 --recode -c | bgzip -c > ~{output_vcf_basename}.vcf.gz
    bcftools index -t ~{output_vcf_basename}.vcf.gz 
  >>>

  runtime {
    docker: "skwalker/imputation:with_vcftools" # TODO: use a public one (not suse bcftools biocontainers also has vcftools )
    memory: "16 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "~{output_vcf_basename}.vcf.gz"
    File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}

task SortIds {
    input {
        File vcf
        String basename
    Int disk_space =  3*ceil(size(vcf, "GB"))
    }

    command <<<
    # what better way is there to do this I really don't know
    zcat ~{vcf} | awk -v OFS='\t' '{split($3, n, ":"); if ( !($1 ~ /^"#"/) && n[4] < n[3])  $3=n[1]":"n[2]":"n[4]":"n[3]; print $0}' | bgzip -c > ~{basename}.vcf.gz
    bcftools index -t ~{basename}.vcf.gz

    >>>

    output {
        File output_vcf = "~{basename}.vcf.gz"
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
        
    }

    runtime {
        docker: "skwalker/imputation:with_vcftools"
        disks: "local-disk " + disk_space + " HDD"
        memory: "16 GB"
    }
}

task MergeSingleSampleVcfs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcf_indices
    String output_vcf_basename
   }

   Int disk_size = 3* ceil(size(input_vcfs, "GB") + size(input_vcf_indices, "GB")) + 20

  command <<<
    bcftools merge ~{sep=' ' input_vcfs} -Oz -o ~{output_vcf_basename}.vcf.gz 
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
