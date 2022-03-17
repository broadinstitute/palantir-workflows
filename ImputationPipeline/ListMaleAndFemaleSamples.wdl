version 1.0


workflow ListMaleAndFemaleSamples {
  input {
    File vcf

    String output_name
  }

  call ImputeSex {
    input:
      vcf = vcf,
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


task ImputeSex {
  input {
    File vcf

    String output_name
    Int mem = 16
  }

  Int disk_size = ceil(2.5*(size(vcf, "GB"))) + 20

  command <<<

    plink --vcf ~{vcf} --chr X --const-fid --impute-sex --make-bed --out sex_imputed
  >>>

  runtime {
    docker: "quay.io/ckachuli/plink@sha256:b5ac6b3b4d28e3fd6f9ca09a8d759a0229899ae9626d0b869d3b3f111e83255f"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB"
  }

  output {
    File output_fam = "sex_imputed.fam"
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