version 1.0

workflow SelectVariantsRemoveInterval {
  input {
    File vcf
    File vcf_index
    File bed
  }

  call SelectVariants {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      bed = bed
  }

  output {
    File vcf_out = SelectVariants.vcf_out
    File vcf_index_out = SelectVariants.vcf_index_out
  }
}

task SelectVariants {
  input {
    File vcf
    File vcf_index
    File bed
  }

  String basename = basename(vcf, ".vcf.gz")
  command <<<
    gatk SelectVariants -V ~{vcf} -XL ~{bed} -O ~{basename}.subset.vcf.gz
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.0.0"
    disks: "local-disk 50 HDD"
  }

  output {
    File vcf_out = "~{basename}.subset.vcf.gz"
    File vcf_index_out = "~{basename}.subset.vcf.gz.tbi"
  }
}