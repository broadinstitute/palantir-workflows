version 1.0

workflow ListMaleAndFemaleSamples {
  input {
    File vcf

    String output_name
  }

  call LDPruneVcf {
    input:
      vcf = vcf,
      basename = output_name
  }

  call ImputeSex {
    input:
      bed = LDPruneVcf.bed,
      bim = LDPruneVcf.bim,
      fam = LDPruneVcf.fam,
      output_name = output_name
  }

  call ListSamplesBySex as ListFemaleSamples {
    input:
      fam = ImputeSex.output_fam,
      sex_code = 2,
      output_name = output_name + ".female.samples"
  }

  call ListSamplesBySex as ListMaleSamples {
    input:
      fam = ImputeSex.output_fam,
      sex_code = 1,
      output_name = output_name + ".male.samples"
  }

  output {
    File female_samples = ListFemaleSamples.sample_list
    File male_samples = ListMaleSamples.sample_list
  }
}

task LDPruneVcf {
  input {
    File vcf
    String basename
    Int mem = 16
  }

  command <<<
    plink --vcf ~{vcf} --chr X --const-fid --indep-pairwise 1000 50 0.2 --make-bed --out ~{basename}
  >>>

  runtime {
    docker: "quay.io/ckachuli/plink@sha256:b5ac6b3b4d28e3fd6f9ca09a8d759a0229899ae9626d0b869d3b3f111e83255f"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB"
  }

  Int disk_size = ceil(2.5*(size(vcf, "GB"))) + 20

  output {
    File bed = "~{basename}.bed"
    File bim = "~{basename}.bim"
    File fam = "~{basename}.fam"
  }
}


task ImputeSex {
  input {
    File bed
    File bim
    File fam

    String output_name
    Int mem = 16
  }

  Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB"))) + 20

  command <<<
    ln -s ~{bed} input.bed
    ln -s ~{bim} input.bim
    ln -s ~{fam} input.fam

    plink --bfile input --chr X --const-fid --impute-sex --make-bed --out sex_imputed
  >>>

  runtime {
    docker: "quay.io/ckachuli/plink@sha256:b5ac6b3b4d28e3fd6f9ca09a8d759a0229899ae9626d0b869d3b3f111e83255f"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB"
  }

  output {
    File output_fam = "sex_imputed.fam"
    File sex_check = "sex_imputed.sexcheck"
  }
}

task ListSamplesBySex {
  input {
    File fam
    Int sex_code # ('1' = male, '2' = female, '0' = unknown)
    String output_name
  }

  Int disk_size =  ceil(2*size(fam, "GB")) + 20

  command <<<
    awk '$5==~{sex_code}{print $2}' ~{fam} > ~{output_name}
  >>>

  runtime {
    docker: "ubuntu:20.04"
    disks: "local-disk " + disk_size + " SSD"
  }

  output {
    File sample_list = "~{output_name}"
  }
}