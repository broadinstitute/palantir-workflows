version 1.0

import "PCATasks.wdl" as PCATasks

workflow PopulationPCAFromImputedCohort {
  input {
    Array[File] population_imputed_vcfs
    String input_basename
    String output_basename
    Int? n_samples_thin
    File? excluded_samples
    File? keep_samples
    File? pruning_sites
    Array[String] contigs
  }

  scatter (population_imputed_vcf in population_imputed_vcfs) {
    call PCATasks.ArrayVcfToPlinkDataset {
      input:
        vcf = population_imputed_vcf,
        basename = input_basename,
        additional_arguments = ["--require-info TYPED"]
    }
  }

  if (length(ArrayVcfToPlinkDataset.bed) > 1) {
    call MergePlinkFilesets {
      input:
        beds = ArrayVcfToPlinkDataset.bed,
        bims = ArrayVcfToPlinkDataset.bim,
        fams = ArrayVcfToPlinkDataset.fam,
        output_basename = output_basename
    }
  }

  if (!defined(pruning_sites)) {
    scatter (contig in contigs) {
      call LDPrune {
        input:
          bim = select_first([MergePlinkFilesets.bim, ArrayVcfToPlinkDataset.bim[0]]),
          bed = select_first([MergePlinkFilesets.bed, ArrayVcfToPlinkDataset.bed[0]]),
          fam = select_first([MergePlinkFilesets.fam, ArrayVcfToPlinkDataset.fam[0]]),
          contig = contig,
          keep_samples = keep_samples,
          n_samples_thin = n_samples_thin,
          output_basename = output_basename + "." + contig
      }
    }

    call ConcatenateLists {
      input:
        lists = LDPrune.prune_in,
        output_name = output_basename + ".prune.in"
    }
  }

  call PCATasks.PrunePopulation {
    input:
      bim = select_first([MergePlinkFilesets.bim, ArrayVcfToPlinkDataset.bim[0]]),
      bed = select_first([MergePlinkFilesets.bed, ArrayVcfToPlinkDataset.bed[0]]),
      fam = select_first([MergePlinkFilesets.fam, ArrayVcfToPlinkDataset.fam[0]]),
      excluded_samples = excluded_samples,
      keep_samples = keep_samples,
      output_basename = output_basename + ".pruned",
      pruning_sites = select_first([pruning_sites, ConcatenateLists.concatenated_lists])
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
    File? pruning_sites_for_pca = ConcatenateLists.concatenated_lists
  }
}


task LDPrune {
  input {
    File bim
    File bed
    File fam
    File? keep_samples
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
    --chr ~{contig} ~{"--keep " + keep_samples} \
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

task MergePlinkFilesets {
  input {
    Array[File] beds
    Array[File] bims
    Array[File] fams

    String output_basename
    Int mem = 8
  }

  Int disk_size =  ceil(2.2*(size(beds, "GB") + size(bims, "GB") + size(fams, "GB"))) + 100

  command <<<
    paste ~{write_lines(beds)} ~{write_lines(bims)} ~{write_lines(fams)} > merge.list

    /plink2 --pmerge-list merge.list bfile --make-bed --out ~{output_basename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB"
  }

  output {
    File bed = "~{output_basename}.bed"
    File bim = "~{output_basename}.bim"
    File fam = "~{output_basename}.fam"
  }
}



