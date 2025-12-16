version 1.0

  workflow VariantsToTable {
    input {
      File vcf
      File vcf_index
    }

    call VariantsToTable_t {
      input:
        vcf = vcf,
        vcf_index = vcf_index
    }

    output {
      File table = VariantsToTable_t.table
    }
  }

  task VariantsToTable_t {
    input {
      File vcf
      File vcf_index
    }

    String basename = basename(vcf, ".vcf.gz")
    command <<<
    gatk VariantsToTable -V ~{vcf} -F CHROM -F POS -F REF -F ALT -F AF -F AC -O ~{basename}.tsv
  >>>

    runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.0.0"
      disks: "local-disk 50 HDD"
    }

    output {
      File table = "~{basename}.tsv"
    }
  }