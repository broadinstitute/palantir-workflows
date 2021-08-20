version 1.0

workflow APOL1CKDRiskAdjustment {
  input {
    File adjustedScores
    File vcf
    File risk_alleles
  }

  call GenotypeG1G2RiskAlleles {
    input:
      vcf = vcf,
      risk_alleles = risk_alleles
  }

}

task GenotypeG1G2RiskAlleles {
  input {
    File vcf
    File risk_alleles
    Int mem = 8
  }

  Int disk_size = ceil(1.2*size(vcf, "GB")) + 50

  command <<<
    /plink2 --vcf ~{vcf} --export A --export-allele ~{risk_alleles} --out apol1_g1_g2_counts
  >>>

  runtime {
     docker: "skwalker/plink2:first"
     disks: "local-disk " + disk_size + " HDD"
     memory: mem + " GB"
  }

  output {
    File risk_allele_counts = "apol1_g1_g2_counts.raw"
  }
}

task AdjustRisk {
  input {
    File risk_allele_counts
    File adjusted_scores
  }

  Int disk_size = ceil(1.2*(size(risk_allele_counts, "GB") + size(adjusted_scores, "GB"))) + 50

  command <<<
    Rscript -<< "EOF"
      library(dplyr)
      library(readr)

      risk_allele_counts <- read_tsv("~{risk_allele_counts}")
      adjusted_scores <- read_tsv("~{adjusted_scores}")

      risk_allele_counts <- risk_allele_counts %>% transmute(IID, apol1_high_risk = ifelse(.[[7]] + .[[8]] + .[[9]] >= 3, 1, 0))
      adjusted_scores <- inner_join(adjusted_scores, risk_allele_counts)
      adjusted_scores <- adjusted_score %>% mutate(adjusted_score_with_apol1 = adjusted_score + apol1_high_risk)

      write_tsv(adjusted_scores, "adjusted_scores_with_apol1.tsv")
    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse:4.1.0"
    disks : "local-disk " + disk_size + " HDD"
  }

  output {
    File adjusted_scores_with_apol1 = "adjusted_scores_with_apol1.tsv"
  }
}