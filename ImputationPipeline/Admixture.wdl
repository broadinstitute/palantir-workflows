version 1.0

workflow Admixture {
  input {
    File bed
    File bim
    File fam
    Int n_ancestral_populations
  }

  call Admixture_t {
    input:
      bed = bed,
      bim = bim,
      fam = fam,
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