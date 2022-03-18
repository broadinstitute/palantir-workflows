version 1.0

import "PCATasks.wdl" as PCATasks

workflow ProjectOntoPCs {
  input {
    File vcf

    File pruning_sites
    File pc_loadings
    File pc_meansd
    String basename
  }

  call PCATasks.ArrayVcfToPlinkDataset {
    input:
      vcf = vcf,
      pruning_sites = pruning_sites,
      basename = basename
  }

  call PCATasks.ProjectArray {
    input:
      bim = ArrayVcfToPlinkDataset.bim,
      bed = ArrayVcfToPlinkDataset.bed,
      fam = ArrayVcfToPlinkDataset.fam,
      pc_loadings = pc_loadings,
      pc_meansd = pc_meansd,
      basename = basename
  }

  output {
    File projections = ProjectArray.projections
  }
}