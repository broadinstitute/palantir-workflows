version 1.0

# When you use a new population dataset 

workflow PerformPopulationPCA {
  input {
    File population_vcf # Like Thousand Genomes
    File population_vcf_index # Like Thousand Genomes
    String basename # what the outputs will be named
    Array[File] imputed_array_vcfs # limit to TYPED and TYPED_ONLY sites before LD pruning.  Also will limit population to sites in imputed vcf for scoring correction
    Array[File] original_array_vcfs
    Array[File]? subset_to_sites
  }
 
  # this task seaparates multiallelics and changes variant IDs to chr:pos:ref:alt1 (bc there are no multiallelics now, alt1=alt)
  call SeparateMultiallelics {
    input:
      original_vcf = population_vcf,
      original_vcf_index = population_vcf_index,
      output_basename = basename + ".no_multiallelics"
  }

    #  we use sorted variant IDs so this step makes sure the variant IDs are in the format of chr:pos:allele1:allele2 where allele1
    # and allele2 are sorted
    call SortVariantIds {
      input:
        vcf = SeparateMultiallelics.output_vcf,
        basename = basename + ".sorted_ids"
    }

  scatter (imputed_array_vcf in imputed_array_vcfs) {
	  call UpdateVariantIds as UpdateVariantIdsImputed {
		input:
			vcf = imputed_array_vcf,
			basename = basename + ".imputed_array.updated_ids."
	  }



	  call SortVariantIds as SortVariantIdsImputedArray {
		input:
			vcf = UpdateVariantIdsImputed.output_vcf,
			basename = basename + ".imputed_array.sorted_ids"
	  }

	  call SelectTypedSites {
		input:
			vcf = SortVariantIdsImputedArray.output_vcf,
			basename = basename
	  }

	  call ExtractIDs as ExtractIDsTyped {
		  input:
			  vcf = SelectTypedSites.output_vcf,
			  output_basename = basename
		}
	}

  call SubsetToArrayVCF {
    input:
        vcf = SortVariantIds.output_vcf,
        vcf_index = SortVariantIds.output_vcf_index,
        intervals = SelectTypedSites.output_vcf,
        intervals_index = SelectTypedSites.output_vcf_index,
        basename = basename + ".sorted_ids.subsetted"
  }

	scatter (original_array_vcf in original_array_vcfs) {
	  call SelectSitesOriginalArray {
		input:
			vcf = original_array_vcf,
			basename = basename
	  }
  }
 
  # this performs some basic QC steps (filtering by MAF, HWE, etc.), as well as 
  # generating a plink-style bim,bed,fam format that has been limited to LD pruned
  # sites. alternatively, if you already have a list of LDPruned sites you want to use,
  # you can run the LDPruneToSites task that is at the bottom of this wdl
  call LDPruning {
    input:
      vcf = SubsetToArrayVCF.output_vcf,
      basename = basename,
      imputed_typed_sites = ExtractIDsTyped.ids,
      original_array_sites = SelectSitesOriginalArray.ids,
      selected_sites = subset_to_sites
  }
  
  # perform PCA using flashPCA
  call PerformPCA {
    input:
      bim = LDPruning.bim,
      bed = LDPruning.bed,
      fam = LDPruning.fam,
      basename = basename
  }

  # see how well your PCA performed: according to the flashPCA website, you want your
  # mean squared error to be low (<1e-8). You will have to read the log file to check
  # this value 
  call CheckPCA {
    input:
      bim = LDPruning.bim,
      bed = LDPruning.bed,
      fam = LDPruning.fam,
      eigenvectors = PerformPCA.eigenvectors,
      eigenvalues = PerformPCA.eigenvalues,
      basename = basename
  }
  # these are the files you need to perform the adjusted scoring
  output {
    File population_loadings = PerformPCA.pc_loadings
    File population_meansd = PerformPCA.mean_sd
    File population_pcs = PerformPCA.pcs 
    File pruning_sites_for_pca = LDPruning.prune_in 
    File sorted_variant_id_dataset = SortVariantIds.output_vcf # this is what you should use as your population dataset for the 
    # ScoringPart, since all the IDs will be matching 
    File sorted_variant_id_dataset_index = SortVariantIds.output_vcf_index
  }
}

task SelectSitesOriginalArray {
	input {
		File vcf
		String basename
		Int mem = 8
	}

	Int disk_size =  ceil(size(vcf, "GB")) + 50

	command <<<
		/plink2 --vcf ~{vcf} \
		--set-all-var-ids @:#:\$1:\$2 \
		--rm-dup force-first \
		--geno 0.001 \
		--maf 0.01 \
		--snps-only \
		--write-snplist \
		--out ~{basename}_selected
	>>>

 	runtime {
		docker: "skwalker/plink2:first"
		disks: "local-disk 400 HDD"
		memory: mem + " GB"
  }

  output {
  	File ids = "~{basename}_selected.snplist"
  }
}

task SelectTypedSites {
	input {
		File vcf
		String basename
	}

	Int disk_size =  ceil(size(vcf, "GB")) + 50

	parameter_meta {
		vcf: {
			localization_optional : true
		}
	}

	command <<<
		gatk SelectVariants -V ~{vcf} -select "TYPED || TYPED_ONLY" -O ~{basename}.vcf.gz
	>>>

	runtime {
	docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
	disks: "local-disk " + disk_size + " HDD"
	memory: "16 GB"
  }

  output {
	File output_vcf = "~{basename}.vcf.gz"
	File output_vcf_index = "~{basename}.vcf.gz.tbi"
  }
}

# Note, we exclude sites on chromosome X because we currently do not impute chromosome X
task LDPruning {
  input {
    File vcf
    Array[File] original_array_sites
    Array[File] imputed_typed_sites
    Array[File]? selected_sites
    Int mem = 8
    String basename
  }
   
  # all these numbers are from Wallace Wang
  command {
  
    /plink2 --vcf ~{vcf} \
    --rm-dup force-first \
    --geno 0.05 \
    --hwe 1e-10 \
    --extract-intersect ~{sep=" " original_array_sites} ~{sep = " " imputed_typed_sites} ~{sep=" " selected_sites} \
    --indep-pairwise 1000 50 0.2 \
    --maf 0.01 \
    --allow-extra-chr \
    --not-chr X \
    --out ~{basename} 

    /plink2 --vcf ~{vcf} \
    --rm-dup force-first \
    --keep-allele-order \
    --extract ~{basename}.prune.in \
    --make-bed \
    --allow-extra-chr \
    --not-chr X \
    --out ~{basename}

  }

  output {
    File prune_in = "~{basename}.prune.in"
    File prune_out = "~{basename}.prune.out"
    File bed = "~{basename}.bed"
    File bim = "~{basename}.bim"
    File fam = "~{basename}.fam"
  }

  runtime {
    docker: "skwalker/plink2:first"
    disks: "local-disk 400 HDD"
    memory: mem + " GB"
  }
} 

task PerformPCA {
  input {
    File bim
    File bed
    File fam
    String basename
    Int mem = 8
  }

  # again, based on Wallace commands
  command {
    cp ~{bim} ~{basename}.bim
    cp ~{bed} ~{basename}.bed
    cp ~{fam} ~{basename}.fam
  
    ~/flashpca/flashpca --bfile ~{basename} -n 16 -d 20 --outpc ${basename}.pc \
    --outpve ${basename}.pc.variance --outload ${basename}.pc.loadings \
    --outmeansd ${basename}.pc.meansd
  }

  output {
    File pcs = "${basename}.pc"
    File pc_variance = "${basename}.pc.variance"
    File pc_loadings = "${basename}.pc.loadings"
    File mean_sd = "${basename}.pc.meansd"
    File eigenvectors = "eigenvectors.txt"
    File eigenvalues = "eigenvalues.txt"
  }

  runtime {
    docker: "skwalker/flashpca:v1"
    disks: "local-disk 400 HDD"
    memory: mem + " GB"
  }
} 

task SeparateMultiallelics {
  input {
    File original_vcf
    File original_vcf_index
    String output_basename
    Int disk_size =  2*ceil(size(original_vcf, "GB"))
  }
  command {
    bcftools norm -m - ~{original_vcf} -Ou | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o ~{output_basename}.vcf.gz
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
  }
  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4 GB"
  }
}

task ExtractIDs {
    input {
        File vcf
        String output_basename
        Int disk_size = 2*ceil(size(vcf, "GB")) + 100
    }

    command <<<
        bcftools query -f "%ID\n" ~{vcf} -o ~{output_basename}.original_array.ids
    >>>
    output {
        File ids = "~{output_basename}.original_array.ids"
    }
    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        disks: "local-disk " + disk_size + " HDD"
        memory: "4 GB"
    }
}

task CheckPCA {
  input {
    File bim
    File bed
    File fam
    File eigenvectors
    File eigenvalues
    String basename
    Int mem = 8
  }

  command {
    cp ~{bim} ~{basename}.bim
    cp ~{bed} ~{basename}.bed
    cp ~{fam} ~{basename}.fam
    
    cp ~{eigenvectors} eigenvectors.txt
    cp ~{eigenvalues} eigenvalues.txt
  
    ~/flashpca/flashpca --bfile ~{basename} --check --verbose \
    --outvec eigenvectors.txt --outval eigenvalues.txt 
  }

  output {
    File output_eigenvectors = "eigenvectors.txt"
    File output_eigenvalues = "eigenvalues.txt"    
  }

  runtime {
    docker: "skwalker/flashpca:v1"
    disks: "local-disk 400 HDD"
    memory: mem + " GB"
  }
}

# if choose to run this task, make sure that the IDs in both your population vcf and
# your list of sites to prune to are exactly the same. make sure that the output bed file
# has the same number of sites are your input ld pruning file !
task LDPruneToSites {
  input {
    File vcf
    File pruning_sites
    Int mem = 8
    String basename
  }
  
  command {
    /plink2 --vcf ~{vcf} \
    --rm-dup force-first \
    --keep-allele-order \
    --extract ~{pruning_sites} \
    --make-bed \
    --allow-extra-chr \
    --out ~{basename}
  }

  output {
    File bed = "~{basename}.bed"
    File bim = "~{basename}.bim"
    File fam = "~{basename}.fam"
  }

  runtime {
    docker: "skwalker/plink2:first"
    disks: "local-disk 400 HDD"
    memory: mem + " GB"
  }
} 

task SortVariantIds {
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

task UpdateVariantIds {

  input {
    File vcf
    String basename
    Int disk_space =  3*ceil(size(vcf, "GB"))
  }

  command <<<
    bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' ~{vcf} -O z -o ~{basename}.vcf.gz
  >>>

  output {
    File output_vcf = "~{basename}.vcf.gz"
  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}

task SubsetToArrayVCF {
  input {
    File vcf
    File vcf_index
    Array[File] intervals
    Array[File] intervals_index
    String basename 
    Int disk_size = 3*ceil(size([vcf, intervals, vcf_index], "GB")) + 20
  }

    command {
   gatk SelectVariants -V ~{vcf} -L ~{sep =" -L " intervals} --interval-set-rule INTERSECTION -O ~{basename}.vcf.gz
   }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4.5 GB"
  }

  output {
    File output_vcf = "~{basename}.vcf.gz"
  }

}
