# BuildHaplotypeMap
#
# This workflow creates a haplotype map file containing a set of low LD sites for fingerprinting
#
# Inputs:
#
# 1. Genotyped multisample VCF or GVCF
# 2. Sequence dictionary for hap map header (should match header of files to be fingerprinted)
# 3. Pruning parameters, see https://www.cog-genomics.org/plink/2.0/ld for more info. 
#    In general, higher values for prune window and prune slide result in more pruning, while lower values for prune cutoff and maf result in less pruning
#  
#   
#

# Filter out indels and variants that weren't called in every sample

task vcftools {
    String? intermediates_prefix
    File input_vcf
    Int? max_missing
    
    Int maxmissing = select_first([max_missing, 0])
    String int_prefix = select_first([intermediates_prefix, "int"])
    
    command <<<
    vcftools --vcf ${input_vcf} --max-missing-count ${maxmissing} --remove-indels --recode --recode-INFO-all --out ${int_prefix}
    >>>
  
    runtime {
    disks: "local-disk 1000 HDD"
    memory: "3500 MB"
    docker: "biocontainers/vcftools:v0.1.16-1-deb_cv1"
  }

  output {
    File recode_vcf = "${int_prefix}.recode.vcf"
  }
}

# Set missing variant IDs to chr:pos

task bcftools {
    String? intermediates_prefix
    File input_vcf
    
    String int_prefix = select_first([intermediates_prefix, "int"])
    
    command <<<
    bcftools annotate ${input_vcf} --set-id '%CHROM\:%POS' -o ${int_prefix}.annotated.vcf
    >>>
  
    runtime {
    disks: "local-disk 1000 HDD"
    memory: "3500 MB"
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
  }

  output {
    File annotated_vcf = "${int_prefix}.annotated.vcf"
  }
}

# Prune variants that are in high LD, as well as variants with a low MAF

task plink {
  File input_vcf
  String? intermediates_prefix
  String output_prefix
  Int? prune_window
  Int? prune_slide
  Float? prune_cutoff
  Float? min_maf
  
  String int_prefix = select_first([intermediates_prefix, "int"])
  Float maf = select_first([min_maf, 0.4])
  Int window = select_first([prune_window, 50])
  Int slide = select_first([prune_slide, 5])
  Float cutoff = select_first([prune_cutoff, 0.5])

  command <<<
    plink1.9 --vcf ${input_vcf} --snps-only --maf ${maf} --biallelic-only --indep-pairwise ${window} ${slide} ${cutoff} --out ${int_prefix}
    >>>

  runtime {
    disks: "local-disk 1000 HDD"
    memory: "3500 MB"
    docker: "biocontainers/plink1.9:v1.90b3.45-170113-1-deb_cv1"
  }

  output {
    File prune_in = "${int_prefix}.prune.in"
  }
}

# Format remaining variants into haplotype map format readable by Picard tools

task reformat{
  String output_prefix
  File prune_in
    File input_vcf
  File sequence_dict

  command <<<
python <<CODE
import pandas as pd
import vcf

with open("${output_prefix}" + ".hapmap.txt", 'w') as new_map:
  sites_set = set()
  with open("${prune_in}") as p:
    for line in p:
      sites_set.add(line[:-1])
  with open("${sequence_dict}") as h:
    for line in h:
      new_map.write(line)
  column_headers= "#CHROM POS ID  REF ALT INFO" + "\n"
  new_map.write(column_headers)
  vcf_reader = vcf.Reader(open("${input_vcf}"), 'r', encoding='utf-8')
  for record in vcf_reader:
    if record.ID == None:
      if str(record.CHROM) + ":" + str(record.POS) + "\n" not in sites_set:
        continue
    if record.ID not in sites_set:
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
                + str(maf) + '\n'
    new_map.write(newline)
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

task Picard{
    File hapmap
    File reference
    File reference_index
    File reference_dict
    File picard_jar
    String output_prefix

  command <<<
    java -jar ${picard_jar} ConvertHaplotypeDatabaseToVcf \
      OUTPUT=${output_prefix}.hapmap.vcf \
      INPUT=${hapmap} \
      R=${reference}
  >>>
  runtime {
    memory: "32 GB"
    disks: "local-disk 1500 HDD"
    docker: "broadinstitute/picard"
  }
  output {
    File map_vcf = "${output_prefix}.hapmap.vcf"    
  }
}

workflow BuildHapMap {
  File input_vcf
  File sequence_dict
  File reference
  File reference_index
  File reference_dict
  File picard_jar
  String? intermediates_prefix
  String output_prefix
  Int? prune_window
  Int? prune_slide
  Float? prune_cutoff
  Float? min_maf
  Int? max_missing
  
  call vcftools {
    input:
      input_vcf = input_vcf,
        intermediates_prefix = intermediates_prefix
  
  }
  
  call bcftools {
    input:
      input_vcf = vcftools.recode_vcf,
        intermediates_prefix = intermediates_prefix
  
  }

  call plink {
    input:
      input_vcf=bcftools.annotated_vcf,
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
      input_vcf = bcftools.annotated_vcf,
        prune_in = plink.prune_in,
      sequence_dict = sequence_dict
  }
  
  call Picard {
    input:
      output_prefix = output_prefix,
        hapmap = reformat.map,
        reference = reference,
        reference_index = reference_index,
        reference_dict = reference_dict,
        picard_jar = picard_jar
  }
  output {
    File map = reformat.map
    File map_vcf = Picard.map_vcf
  }
}