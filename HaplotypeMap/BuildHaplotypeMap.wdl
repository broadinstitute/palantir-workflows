# BuildHaplotypeMap
#
# This workflow creates a haplotype map file containing a set of low LD sites for fingerprinting
#
# Inputs:
#
# 1. Genotyped multisample GVCF
# 2. Header file. Standard SAM header (should match header of files to be fingerprinted)
# 3. Pruning parameters, see https://www.cog-genomics.org/plink/2.0/ld for more info. 
#    In general, higher values for prune window and prune slide result in more pruning, while lower values for prune cutoff and maf result in less pruning
#  
#   
#

task vcftools {
    String intermediates_prefix
    File input_vcf
    
    command <<<
  	vcftools --vcf ${input_vcf} --max-missing-count 0 --remove-indels --recode --recode-INFO-all --out ${intermediates_prefix}
    >>>
	
    runtime {
    disks: "local-disk 1000 HDD"
    memory: "3500 MB"
    docker: "jeltje/vcftools"
  }

	output {
    File recode_vcf = "${intermediates_prefix}.recode.vcf"
  }
}

task plink {
  File input_vcf
  String intermediates_prefix
  String output_prefix
  Int prune_window
  Int prune_slide
  Float prune_cutoff
  Float min_maf
  
  String set_missing = "@:#"

  command <<<
  	plink --vcf ${input_vcf} --snps-only --maf ${min_maf} --biallelic-only --set-missing-var-ids ${set_missing} --indep-pairwise ${prune_window} ${prune_slide} ${prune_cutoff} --out ${intermediates_prefix}
    >>>

  runtime {
    disks: "local-disk 1000 HDD"
    memory: "3500 MB"
    docker: "jrose77/plinkdocker"
  }

  output {
    File prune_in = "${intermediates_prefix}.prune.in"
  }
}

task reformat{
	String output_prefix
	File prune_in
    File recode_vcf
	File header_file

	command <<<
python <<CODE
import pandas as pd
import vcf


vcf_file = "${recode_vcf}"
prune_in = "${prune_in}"
header = open("${header_file}")
new_map = open("${output_prefix}.hapmap.txt", 'w')
vcf_reader = vcf.Reader(open(vcf_file, 'r', encoding='utf-8'))
vcf_chr = False
first = next(vcf_reader)
if str(first.CHROM)[0:3] == 'chr':
	vcf_chr = True
sites_to_include = open(prune_in)
sites_list = []
for line in sites_to_include:
	if vcf_chr and (line[0:3] != 'chr'):
		sites_list.append('chr'+line)
	else:
		sites_list.append(line)
for line in header:
	new_map.write(line)
new_map.write("\n")
header.close()
for record in vcf_reader:
	if record.ID == None:
		if str(record.CHROM) + ":" + str(record.POS) + "\n" not in sites_list:
			continue
	elif record.ID not in sites_list:
		continue
	ref = record.REF[0]
	alt = record.ALT[0]
	if len(alt) > 1:
		continue
	if alt == '*':
		continue
	maf = record.INFO['AF'][0]
	newline = str(record.CHROM) + '\t' \
              + str(record.POS) + '\t' \
              + str(record.CHROM) + ":" + str(record.POS) + '\t' \
              + str(ref) + '\t' \
              + str(alt) + '\t' \
              + str(maf) + '\t' \
              + '\t' + "" + '\t' + '\n'
	new_map.write(newline)
new_map.close()
CODE
  >>>
  
    runtime {
      disks: "local-disk 1000 HDD"
      memory: "3500 MB"
      docker: "docker.io/mollysacks/python:hapmap_builder"
    }
    
    output {
    File map = "${output_prefix}.hapmap.txt"
  }
}

workflow BuildHapMap {
  File input_vcf
  File header_file
  String intermediates_prefix
  String output_prefix
  Int prune_window
  Int prune_slide
  Float prune_cutoff
  Float min_maf
  
  call vcftools {
  	input:
    	input_vcf = input_vcf,
        intermediates_prefix = intermediates_prefix
  
  }

  call plink {
    input:
      input_vcf=vcftools.recode_vcf,
      intermediates_prefix = intermediates_prefix,
      output_prefix = output_prefix,
      prune_window = prune_window,
      prune_slide = prune_slide,
      prune_cutoff = prune_cutoff,
      min_maf = min_maf
  }

  call reformat {
  	input:
    	output_prefix = output_prefix,
  		recode_vcf = vcftools.recode_vcf,
        prune_in = plink.prune_in,
  		header_file = header_file
  }
  
  output {
    File map = reformat.map
  }
}