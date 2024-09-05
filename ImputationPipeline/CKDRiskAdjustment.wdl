version 1.0
import "ScoringTasks.wdl" as ScoringTasks

workflow CKDRiskAdjustment {
  input {
    File adjustedScores
    File vcf
    File risk_alleles
  }

  call ScoringTasks.DetermineChromosomeEncoding {
		input:
			weights = named_weight_set.weight_set.linear_weights
	}

  call GenotypeG1G2RiskAlleles {
    input:
      vcf = vcf,
      risk_alleles = risk_alleles,
      chromosome_encoding = DetermineChromosomeEncoding.chromosome_encoding
  }

  call AdjustRisk {
    input:
      risk_allele_counts = GenotypeG1G2RiskAlleles.risk_allele_counts,
      adjusted_scores = adjustedScores
  }

  output {
    File adjusted_scores_with_apol1 = AdjustRisk.adjusted_scores_with_apol1
  }
}

task GenotypeG1G2RiskAlleles {
  input {
    File vcf
    File risk_alleles
    Int mem = 8
    String? chromosome_encoding
  }

  Int disk_size = ceil(1.2*size(vcf, "GB")) + 50
  Int plink_mem = ceil(mem * 0.75 * 1000)

  command <<<
    /plink2 --vcf ~{vcf} --extract ~{risk_alleles} --export A --export-allele ~{risk_alleles} \
    --allow-extra-chr --set-all-var-ids @:#:\$1:\$2 --new-id-max-allele-len 1000 missing \
    --memory ~{plink_mem} --out apol1_g1_g2_counts ~{"--output-chr " + chromosome_encoding}
  >>>

  runtime {
     docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
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
      library(tidyr)

      risk_allele_counts <- read_tsv("~{risk_allele_counts}")
      adjusted_scores <- read_tsv("~{adjusted_scores}")

      risk_allele_counts <- risk_allele_counts %>% mutate(G1_1 = replace_na(.[[7]], 0), G1_2 = replace_na(.[[8]], 0), G2 = replace_na(.[[9]], 0))
      risk_allele_counts <- risk_allele_counts %>% transmute(IID, apol1_high_risk = ifelse(pmax(G1_1, G1_2) + G2 >= 2, 1, 0))
      adjusted_scores <- inner_join(adjusted_scores, risk_allele_counts)
      adjusted_scores <- adjusted_scores %>% mutate(adjusted_score = adjusted_score + apol1_high_risk)
      adjusted_scores <- adjusted_scores %>% mutate(percentile = pnorm(adjusted_score))

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