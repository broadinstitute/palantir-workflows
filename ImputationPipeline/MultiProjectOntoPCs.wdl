version 1.0

import "ProjectOntoPCs.wdl" as ProjectOntoPCs

workflow MultiProjectOntoPCs {
  input {
    Array[File] vcfs

    File pruning_sites
    File pc_loadings
    File pc_meansd
    String basename
  }

  scatter(vcf in vcfs) {
    call ProjectOntoPCs.ProjectOntoPCs {
      input:
        vcf = vcf,
        pruning_sites = pruning_sites,
        pc_loadings = pc_loadings,
        pc_meansd = pc_meansd,
        basename = basename(vcf, ".vcf.gz")
    }
  }

  call CombineTables {
    input:
      tables = ProjectOntoPCs.projections,
      basename = basename
  }

  output {
    File projections = CombineTables.combined_table
  }
}

task CombineTables {
  input {
    Array[File] tables
    String basename
  }

  command <<<
    Rscript -<< "EOF"
    library(dplyr)
    library(readr)
    library(purrr)

    table <- list("~{sep='","' tables}") %>% map(read_tsv) %>% reduce(bind_rows)

    write_tsv(table, "~{basename}.tsv")

    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse:4.1.0"
    disks: "local-disk 100 HDD"
  }

  output {
    File combined_table = "~{basename}.tsv"
  }
}