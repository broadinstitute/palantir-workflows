version 1.0

import "../Structs.wdl"

workflow SubsetWeightSet {
  input {
    WeightSet weight_set
    File sites_to_subset_to
  }

  call SubsetSitesBasedFile as SubsetLinearWeights {
    input:
      sites_to_subset_to = sites_to_subset_to,
      sites_based_file = weight_set.linear_weights,
      cols_to_subset_on = ["1"]
  }

  if (defined(weight_set.interaction_weights)) {
    call SubsetSitesBasedFile as SubsetInteractionWeights {
      input:
        sites_to_subset_to = sites_to_subset_to,
        sites_based_file = select_first([weight_set.interaction_weights]),
        cols_to_subset_on = ["id_1", "id_2"]
    }
  }

  if (defined(weight_set.interaction_self_exclusive_sites)) {
    SelfExclusiveSites input_self_exclusive_sites = select_first([weight_set.interaction_self_exclusive_sites])
    call SubsetSitesBasedFile as SubsetSelfExclusiveSites {
      input:
        sites_to_subset_to = sites_to_subset_to,
        sites_based_file = input_self_exclusive_sites.sites,
        cols_to_subset_on = ["id"]
    }

    SelfExclusiveSites output_self_exclusive_sites = object{sites : SubsetSelfExclusiveSites.subset_file, maxAllowed : input_self_exclusive_sites.maxAllowed}
  }

  output {
    WeightSet subsetted_weight_set = object{linear_weights : SubsetLinearWeights.subset_file, interaction_weights : SubsetInteractionWeights.subset_file, self_exclusive_sites : output_self_exclusive_sites}
    Int n_original_weights = SubsetLinearWeights.n_original_sites + select_first([SubsetInteractionWeights.n_original_sites, 0])
    Int n_subset_weights = SubsetLinearWeights.n_subset_sites + select_first([SubsetInteractionWeights.n_subset_sites, 0])

  }
}

task SubsetSitesBasedFile {
  input {
    File sites_to_subset_to
    File sites_based_file
    Array[String] cols_to_subset_on
  }

  String output_base = basename(sites_to_subset_to)
  command <<<

    Rscript -<< "EOF"
    library(readr)
    library(dplyr)

    sites_based_file <- read_tsv("~{sites_based_file}")
    sites_to_subset_to <- read_lines("~{sites_to_subset_to}")

    sites_based_file_subset <- sites_based_file %>% filter(if_all(c(~{sep=","cols_to_subset_on}), ~ . %in% sites_to_subset_to))
    write_tsv(sites_based_file_subset, "subset_~{output_base}")
    write_tsv(sites_based_file_subset %>% count(), "n_original.txt", col_names = FALSE)
    write_tsv(sites_based_file_subset %>% count(), "n_subset.txt", col_names = FALSE)

    EOF

  >>>

  runtime {
    docker: "rocker/tidyverse"
    disks: "local-disk 100 HDD"
    memory: "16 GB"
  }

  output {
    File subset_file = "subset_~{output_base}"
    Int n_original_sites = read_int("n_original.txt")
    Int n_subset_sites = read_int("n_subset.txt")
  }
}