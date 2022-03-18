version 1.0

workflow Admixture {
  input {
    File bed
    File bim
    File fam
    Int n_ancestral_populations
    Int? thin_count
  }

  if(defined(thin_count)) {
    call ThinVariants {
      input:
        bed = bed,
        bim = bim,
        fam = fam,
        thin_count = select_first([thin_count])
    }
  }

  call Admixture_t {
    input:
      bed = select_first([ThinVariants.thinned_bed, bed]),
      bim = select_first([ThinVariants.thinned_bim, bim]),
      fam = select_first([ThinVariants.thinned_fam, fam]),
      n_ancestral_populations = n_ancestral_populations
  }

}

task Admixture_t {
  input {
    File bed
    File bim
    File fam
    Int n_ancestral_populations
    Boolean cv = false
    Int mem = 16
    Int n_cpus = 4
  }

  Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
  String basename = basename(bed, ".bed")

  command <<<

    /admixture_linux-1.3.0/admixture ~{if (cv) then "--cv" else ""} ~{bed} ~{n_ancestral_populations} -j~{n_cpus}
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/admixture_docker:v1.0.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB"
    cpu: n_cpus
  }

  output {
    File ancestry_fractions = "~{basename}.~{n_ancestral_populations}.Q"
    File allele_frequencies = "~{basename}.~{n_ancestral_populations}.P"
  }
}

task ThinVariants {
  input {
    File bed
    File bim
    File fam
    Int thin_count
    Int mem = 8
  }

  Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
  String basename = basename(bed, ".bed")

  command <<<
    ln -s ~{bim} input.bim
    ln -s ~{bed} input.bed
    ln -s ~{fam} input.fam

    /plink2 --bfile input \
    --thin-count ~{thin_count} \
    --make-bed \
    --out ~{basename}.thined

  >>>

  output {
    File thinned_bed = "~{basename}.thined.bed"
    File thinned_bim = "~{basename}.thined.bim"
    File thinned_fam = "~{basename}.thined.fam"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB"
  }
}