version 1.0

import "PCATasks.wdl" as PCATasks

workflow PopulationPCAFromImputedCohort {
  input {
    File population_imputed_vcf
    String output_basename
    Int? n_samples_thin
    File? excluded_samples
    Array[String] contigs
  }

  call PCATasks.ArrayVcfToPlinkDataset {
    input:
      vcf = population_imputed_vcf,
      basename = output_basename,
      additional_arguments = ["--require-info TYPED"]
  }

  scatter (contig in contigs) {
    call LDPrune {
      input:
        bim = ArrayVcfToPlinkDataset.bim,
        bed = ArrayVcfToPlinkDataset.bed,
        fam = ArrayVcfToPlinkDataset.fam,
        contig = contig,
        n_samples_thin = n_samples_thin,
        output_basename = output_basename + "." + contig
    }
  }

  call ConcatenateLists {
    input:
      lists = LDPrune.prune_in,
      output_name = output_basename + ".prune.in"
  }

  call PCATasks.PrunePopulation {
    input:
      bim = ArrayVcfToPlinkDataset.bim,
      bed = ArrayVcfToPlinkDataset.bed,
      fam = ArrayVcfToPlinkDataset.fam,
      excluded_samples = excluded_samples,
      output_basename = output_basename + ".pruned",
      pruning_sites = ConcatenateLists.concatenated_lists
  }

  call PCATasks.PerformPCA {
    input:
      bim = PrunePopulation.out_bim,
      bed = PrunePopulation.out_bed,
      fam = PrunePopulation.out_fam,
      basename = output_basename
  }

  call PCATasks.CheckPCA {
    input:
      bim = PrunePopulation.out_bim,
      bed = PrunePopulation.out_bed,
      fam = PrunePopulation.out_fam,
      basename = output_basename,
      eigenvectors = PerformPCA.eigenvectors,
      eigenvalues = PerformPCA.eigenvalues,
  }

  output {
    File population_loadings = PerformPCA.pc_loadings
    File population_meansd = PerformPCA.mean_sd
    File population_pcs = PerformPCA.pcs
    File pruning_sites_for_pca = ConcatenateLists.concatenated_lists
  }
}


task LDPrune {
  input {
    File bim
    File bed
    File fam
    String contig
    Int? n_samples_thin
    String output_basename
    Int mem = 8
  }

  Int disk_space =  2*ceil(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")) + 20

  command <<<
    ln -s ~{bim} input.bim
    ln -s ~{bed} input.bed
    ln -s ~{fam} input.fam

    /plink2 --bfile input \
    --rm-dup exclude-all \
    --geno 0.001 \
    --hwe 1e-10 \
    --snps-only \
    --chr ~{contig} \
    --maf 0.01 ~{"--thin-indiv-count " + n_samples_thin} \
    --indep-pairwise 1000 50 0.2 \
    --out ~{output_basename}
  >>>

  output {
    File prune_in = "~{output_basename}.prune.in"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    disks: "local-disk " + disk_space + " HDD"
    memory: mem + " GB"
    maxRetries: 3
  }
}

task ConcatenateLists {
  input {
    Array[File] lists
    String output_name
  }

  Int disk_size =  ceil(2.2*size(lists, "GB")) + 20

  command <<<
    cat ~{sep = " " lists} > ~{output_name}
  >>>

  runtime {
    docker: "ubuntu:20.04"
    disks: "local-disk " + disk_size + " SSD"
  }

  output {
    File concatenated_lists = "~{output_name}"
  }
}



