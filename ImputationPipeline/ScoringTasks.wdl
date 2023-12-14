version 1.0

import "Structs.wdl"

# score with plink2
task ScoreVcf {
  input {
    File vcf
    String basename
    File weights
    Int base_mem = 8
    String? extra_args
    File? sites
    String? chromosome_encoding
    Boolean use_ref_alt_for_ids = false
  }

  Int runtime_mem = base_mem + 2
  Int plink_mem = ceil(base_mem * 0.75 * 1000)
  Int disk_space =  3*ceil(size(vcf, "GB")) + 20
  String var_ids_string = "@:#:" + if use_ref_alt_for_ids then "\\$r:\\$a" else "\\$1:\\$2"

  command <<<
    /plink2 --score ~{weights} header ignore-dup-ids list-variants no-mean-imputation \
    cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums --set-all-var-ids ~{var_ids_string} --allow-extra-chr ~{extra_args} -vcf ~{vcf} dosage=DS \
    --new-id-max-allele-len 1000 missing ~{"--extract " + sites} --out ~{basename} --memory ~{plink_mem} ~{"--output-chr " + chromosome_encoding}
  >>>

  output {
    File score = "~{basename}.sscore"
    File log = "~{basename}.log"
    File sites_scored = "~{basename}.sscore.vars"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    disks: "local-disk " + disk_space + " HDD"
    memory: runtime_mem + " GB"
  }
}


task AddInteractionTermsToScore {
  input {
    File vcf
    File interaction_weights
    File scores
    File? sites
    String basename
    SelfExclusiveSites? self_exclusive_sites # The interaction term will only be added in no more than selfExclusiveSites.maxAllowed of the
    # effect alleles listed in SelfExclusizeSites.sites is observed

    Float mem = 8
    Int block_buffer=10000000
  }

  Int disk_space =  3*ceil(size(vcf, "GB")) + 20

  command <<<

    tabix ~{vcf}

    python3 << "EOF"
    from cyvcf2 import VCF
    import pandas as pd
    import csv

    vcf = VCF("~{vcf}", lazy=True)
    samples = vcf.samples

    def add_allele_to_count(site, allele, dictionary):
      if site in dictionary:
        dictionary[site][allele]=[0]*len(samples)
      else:
        dictionary[site]={allele:[0]*len(samples)}

    interactions_allele_counts = dict()
    interactions_dict = dict()
    positions = set()
    if ~{if defined(sites) then "True" else "False"}:
      with open("~{sites}") as f_sites:
        sites = {s.strip() for s in f_sites}
    else:
      sites = {}
    with open("~{interaction_weights}") as f:
      for line in csv.DictReader(f, delimiter='\t'):
        site_1 = line['id_1']
        site_2 = line['id_2']
        if len(sites) == 0 or site_1 in sites and site_2 in sites:
          allele_1 = line['allele_1']
          allele_2 = line['allele_2']
          chrom_1 = line['chrom_1']
          chrom_2 = line['chrom_2']
          pos_1 = int(line['pos_1'])
          pos_2 = int(line['pos_2'])
          weight = float(line['weight'])

          add_allele_to_count(site_1, allele_1, interactions_allele_counts)
          add_allele_to_count(site_2, allele_2, interactions_allele_counts)
          interactions_dict[(site_1, allele_1, site_2, allele_2)] = weight
          positions.add((chrom_1, pos_1))
          positions.add((chrom_2, pos_2))

    def add_self_exclusive_site(site, allele, dictionary):
      if site in dictionary:
        dictionary[site].add(allele)
      else:
        dictionary[site]={allele}

    self_exclusive_sites = dict()
    max_self_exclusive_sites = ~{if (defined(self_exclusive_sites)) then select_first([self_exclusive_sites]).maxAllowed else 0}
    self_exclusive_sites_counts = [0]*len(samples)
    if ~{if (defined(self_exclusive_sites)) then "True" else "False"}:
      with open("~{if (defined(self_exclusive_sites)) then select_first([self_exclusive_sites]).sites else ''}") as f_self_exclusive_sites:
        for line in csv.DictReader(f_self_exclusive_sites, delimiter='\t'):
          id = line['id']
          if len(sites) == 0 or id in sites:
            chrom = line['chrom']
            pos = int(line['pos'])
            allele = line['allele']
            add_self_exclusive_site(id, allele, self_exclusive_sites)
            positions.add((chrom, pos))

    #select blocks to read
    positions = sorted(positions)
    current_chrom=positions[0][0]
    current_start=positions[0][1]
    current_end = current_start+1
    buffer=~{block_buffer}

    blocks_to_read=[]
    for site in positions:
      if site[0] != current_chrom or site[1] - current_end > buffer:
        blocks_to_read.append(current_chrom + ":" + str(current_start) + "-" + str(current_end))
        current_chrom=site[0]
        current_start=site[1]
        current_end = current_start+1
      else:
        current_end = site[1] + 1

    #last block
    blocks_to_read.append(current_chrom + ":" + str(current_start) + "-" + str(current_end))

    #count interaction alleles for each sample
    sites_used_in_score = set()
    for block in blocks_to_read:
      for variant in vcf(block):
        alleles = [a for a_l in [[variant.REF], variant.ALT] for a in a_l]
        vid=":".join(s for s_l in [[variant.CHROM], [str(variant.POS)], sorted(alleles)] for s in s_l)
        if vid in interactions_allele_counts:
          sites_used_in_score.add(vid)
          for sample_i,gt in enumerate(variant.genotypes):
            for gt_allele in gt[:-1]:
              allele = alleles[gt_allele]
              if allele in interactions_allele_counts[vid]:
                interactions_allele_counts[vid][allele][sample_i] += 1
        if vid in self_exclusive_sites:
          sites_used_in_score.add(vid)
          for sample_i,gt in enumerate(variant.genotypes):
            for gt_allele in gt[:-1]:
              allele = alleles[gt_allele]
              if allele in self_exclusive_sites[vid]:
                self_exclusive_sites_counts[sample_i] += 1

    #calculate interaction scores for each sample
    interaction_scores = [0] * len(samples)

    def get_interaction_count(site_and_allele_1, site_and_allele_2, sample_i):
      if site_and_allele_1 == site_and_allele_2:
        return interactions_allele_counts[site_and_allele_1[0]][site_and_allele_1[1]][sample_i]//2
      else:
        return min(interactions_allele_counts[site_and_allele_1[0]][site_and_allele_1[1]][sample_i], interactions_allele_counts[site_and_allele_2[0]][site_and_allele_2[1]][sample_i])

    for interaction in interactions_dict:
      for sample_i in range(len(samples)):
        if self_exclusive_sites_counts[sample_i] <= max_self_exclusive_sites:
          site_and_allele_1 = (interaction[0], interaction[1])
          site_and_allele_2 = (interaction[2], interaction[3])
          interaction_scores[sample_i]+=get_interaction_count(site_and_allele_1, site_and_allele_2, sample_i) * interactions_dict[interaction]

    #add interaction scores to linear scores
    df_interaction_score = pd.DataFrame({"sample_id":samples, "interaction_score":interaction_scores}).set_index("sample_id")
    df_scores=pd.read_csv("~{scores}", sep="\t").astype({'#IID':'string'}).set_index("#IID")
    df_scores = df_scores.join(df_interaction_score)
    df_scores['SCORE1_SUM'] = df_scores['SCORE1_SUM'] + df_scores['interaction_score']
    df_scores.to_csv("~{basename}_scores_with_interactions.tsv", sep="\t")
    with open("~{basename}_sites_used_in_interaction_score.ids", "w") as f_sites_used:
      for site in sites_used_in_score:
        f_sites_used.write("%s\n" % site)
    EOF
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/imputation_interaction_python@sha256:40a8fb88fe287c3e3a11022ff63dae1ad5375f439066ae23fe089b2b61d3222e"
    disks: "local-disk " + disk_space + " HDD"
    memory: mem + " GB"
  }

  output {
    File scores_with_interactions = basename + "_scores_with_interactions.tsv"
    File sites_used_in_interaction_score = basename + "_sites_used_in_interaction_score.ids"
  }
}

task CheckWeightsCoverSitesUsedInTraining {
  input {
    File sites_used_in_training
    WeightSet weight_set
  }

  command <<<
    python3 << "EOF"
    import csv
    import sys

    with open("~{sites_used_in_training}") as f_sites_used_in_training:
      sites_used_in_training = {s.strip() for s in f_sites_used_in_training}

    sites_in_weight_set = set()
    with open("~{weight_set.linear_weights}") as f_linear_weights:
      linear_weights_reader = csv.reader(f_linear_weights, delimiter='\t')
      next(linear_weights_reader)
      for line in linear_weights_reader:
        sites_in_weight_set.add(line[0])

    if ~{if (defined(weight_set.interaction_weights)) then "True" else "False"}:
      with open("~{weight_set.interaction_weights}") as f_interaction_weights:
        for line in csv.DictReader(f_interaction_weights, delimiter='\t'):
          sites_in_weight_set.add(line['id_1'])
          sites_in_weight_set.add(line['id_2'])

    sites_missing_from_weight_set = sites_used_in_training - sites_in_weight_set

    if len(sites_missing_from_weight_set) > 0:
      sys.exit(f"Error: {len(sites_missing_from_weight_set)} sites used in model training are missing from weights files.")
    EOF


  >>>

  runtime {
    docker : "python:3.9.10"
  }
}

task CompareScoredSitesToSitesUsedInTraining {
  input {
    File sites_used_in_training
    File sites_used_in_scoring
    WeightSet weight_set
  }

  command <<<
    python3 << "EOF"
    import csv

    with open("~{sites_used_in_training}") as f_sites_used_in_training:
      sites_used_in_training = {s.strip() for s in f_sites_used_in_training}

    with open("~{sites_used_in_scoring}") as f_sites_used_in_scoring:
      sites_used_in_scoring = {s.strip() for s in f_sites_used_in_scoring}

    missing_sites = sites_used_in_training - sites_used_in_scoring

    with open("missing_sites.txt", "w") as f_missing_sites:
      for site in missing_sites:
        f_missing_sites.write(site)

    with open("n_missing_sites.txt", "w") as f_n_missing_sites:
      f_n_missing_sites.write(f"{len(missing_sites)}")

    max_error_up = 0
    max_error_down = 0

    with open("~{weight_set.linear_weights}") as f_weights:
      weights_reader = csv.reader(f_weights, delimiter = "\t")
      next(weights_reader)
      for line in weights_reader:
        id = line[0]
        if id in missing_sites:
          missing_site_weight = float(line[2])
          if missing_site_weight > 0:
            max_error_up += 2 * missing_site_weight
          else:
            max_error_down += 2 * missing_site_weight


    if ~{if defined(weight_set.interaction_weights) then "True" else "False"}:
      with open("~{weight_set.interaction_weights}") as f_interaction_weights:
        for line in csv.DictReader(f_interaction_weights, delimiter='\t'):
          id_1 = line['id_1']
          id_2 = line['id_2']
          if id_1 in missing_sites or id_2 in missing_sites:
            weight_multiplier = 2 if id_1 != id_2 else 1
            missing_site_weight = float(line['weight'])
            if missing_site_weight > 0:
              max_error_up += weight_multiplier * missing_site_weight
            else:
              max_error_down += weight_multiplier * missing_site_weight

    with open("max_error_up.txt", "w") as f_max_error_up:
      f_max_error_up.write(f"{max_error_up}")

    with open("max_error_down.txt", "w") as f_max_error_down:
      f_max_error_down.write(f"{max_error_down}")

    EOF
  >>>

  runtime {
    docker : "python:3.9.10"
  }

  output {
    File missing_sites = "missing_sites.txt"
    Int n_missing_sites = read_int("n_missing_sites.txt")
    Float max_error_up = read_float("max_error_up.txt")
    Float max_error_down = read_float("max_error_down.txt")
  }
}

task CombineScoringSites {
  input {
    File sites_used_linear_score
    File sites_used_interaction_score
    String basename
  }

  Int disk_size = ceil(size(sites_used_linear_score, "GB") + 2*size(sites_used_interaction_score, "GB")) + 50
  command <<<
    cat ~{sites_used_linear_score} ~{sites_used_interaction_score} | sort | uniq > ~{basename}_sites_used_in_score.ids
  >>>

  runtime {
    docker: "ubuntu:20.04"
    disks: "local-disk " + disk_size + " SSD"
  }

  output {
    File combined_scoring_sites = "~{basename}_sites_used_in_score.ids"
  }
}

task AddShiftToRawScores {
  input {
    File raw_scores
    Float shift
    String basename
  }

  command <<<
    Rscript -<< "EOF"
      library(dplyr)
      library(readr)

      scores <- read_tsv("~{raw_scores}")
      shifted_scores <- scores %>% mutate(SCORE1_SUM = SCORE1_SUM + ~{shift})

      write_tsv(shifted_scores, "~{basename}.tsv")
    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse:4.1.0"
    disks: "local-disk 100 HDD"
  }

  output {
    File shifted_scores = "~{basename}.tsv"
  }
}

task CombineMissingSitesAdjustedScores {
  input {
    File adjusted_scores_shifted_up
    File adjusted_scores_shifted_down
    File adjusted_scores
    Int n_missing_sites
    String condition_name
  }

  command <<<
    Rscript -<< "EOF"
    library(dplyr)
    library(readr)

    adjusted_scores <- read_tsv("~{adjusted_scores}") %>% transmute(IID, condition = "~{condition_name}", n_missing_sites = ~{n_missing_sites}, adjusted_score, percentile)
    adjusted_scores_shifted_up <- read_tsv("~{adjusted_scores_shifted_up}") %>% transmute(IID, potential_high_adjusted_score = adjusted_score, potential_high_percentile = percentile)
    adjusted_scores_shifted_down <- read_tsv("~{adjusted_scores_shifted_down}") %>% transmute(IID, potential_low_adjusted_score = adjusted_score, potential_low_percentile = percentile)

    adjusted_scores_shifts <- inner_join(inner_join(adjusted_scores, adjusted_scores_shifted_up), adjusted_scores_shifted_down)
    write_tsv(adjusted_scores_shifts, "missing_sites_shifted_scores.tsv")
    EOF
  >>>

  runtime {
    docker: "rocker/tidyverse:4.1.0"
    disks: "local-disk 100 HDD"
  }

  output {
    File missing_sites_shifted_scores = "missing_sites_shifted_scores.tsv"
  }
}

task TrainAncestryModel {
  input {
    File population_pcs
    File population_scores
    String output_basename
    Int mem = 2
  }

  command <<<
    Rscript -<< "EOF"
      library(dplyr)
      library(readr)
      library(tibble)
      population_pcs = read_tsv("~{population_pcs}")
      population_scores = read_tsv("~{population_scores}")

      population_data = inner_join(population_pcs, population_scores, by=c("IID" = "#IID"))

      # generate the linear model from the population data using the first 4 PCs
      population_model = glm(SCORE1_SUM ~ PC1 + PC2 + PC3 + PC4, data = population_data, family = "gaussian")

      population_data <- population_data %>% mutate(residual_score2 = resid(population_model)^2)

      # generate the linear model for the variance of the score using the first 4 PCs
      population_var_model <- glm(residual_score2 ~ PC1 + PC2 + PC3 + PC4, data = population_data, family = Gamma(link = "log"))

      # use linear model to fit full likelihood model

      # linear transformation to predict variance
      f_sigma2 <- function(t, theta) {
      PC1 = t %>% pull(PC1)
      PC2 = t %>% pull(PC2)
      PC3 = t %>% pull(PC3)
      PC4 = t %>% pull(PC4)
      PC5 = t %>% pull(PC5)
      sigma2 <- exp(theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4)
      }


      # linear transformation to predict mean
      f_mu <- function(t, theta) {
      PC1 = t %>% pull(PC1)
      PC2 = t %>% pull(PC2)
      PC3 = t %>% pull(PC3)
      PC4 = t %>% pull(PC4)
      PC5 = t %>% pull(PC5)
      mu <- theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4
      }


      # negative log likelihood
      nLL_mu_and_var <- function(theta) {
      theta_mu = theta[1:5]
      theta_var = theta[6:10]
      x = population_data %>% pull(SCORE1_SUM)
      sum(log(sqrt(f_sigma2(population_data, theta_var))) + (1/2)*(x-f_mu(population_data, theta_mu))^2/f_sigma2(population_data, theta_var))
      }


      # gradient of negative log likelihood function
      grr <- function(theta) {
      theta_mu = theta[1:5]
      theta_var = theta[6:10]
      d_mu_1 <- 1
      d_mu_2 <- population_data %>% pull(PC1)
      d_mu_3 <- population_data %>% pull(PC2)
      d_mu_4 <- population_data %>% pull(PC3)
      d_mu_5 <- population_data %>% pull(PC4)
      d_sig_7 <- 1 * f_sigma2(population_data, theta_var)
      d_sig_8 <- population_data %>% pull(PC1) * f_sigma2(population_data, theta_var)
      d_sig_9 <- population_data %>% pull(PC2) * f_sigma2(population_data, theta_var)
      d_sig_10 <- population_data %>% pull(PC3) * f_sigma2(population_data, theta_var)
      d_sig_11 <- population_data %>% pull(PC4) * f_sigma2(population_data, theta_var)

      x <- population_data %>% pull(SCORE1_SUM)
      mu_coeff <- -(x - f_mu(population_data, theta_mu))/f_sigma2(population_data, theta_var)
      sig_coeff <- 1/(2*f_sigma2(population_data, theta_var)) -(1/2)*(x - f_mu(population_data, theta_mu))^2/(f_sigma2(population_data, theta_var)^2)


      grad <- c(sum(mu_coeff*d_mu_1),
      sum(mu_coeff*d_mu_2),
      sum(mu_coeff*d_mu_3),
      sum(mu_coeff*d_mu_4),
      sum(mu_coeff*d_mu_5),
      sum(sig_coeff*d_sig_7),
      sum(sig_coeff*d_sig_8),
      sum(sig_coeff*d_sig_9),
      sum(sig_coeff*d_sig_10),
      sum(sig_coeff*d_sig_11)
      )
      }

      # use linear model fits as initial parameters for full likelihood fit
      initial_pars <- c(population_model$coefficients, population_var_model$coefficients)
      initial_pars <- setNames(initial_pars, c("Beta0_mu", "Beta1_mu", "Beta2_mu", "Beta3_mu", "Beta4_mu",
                               "Beta0_var", "Beta1_var", "Beta2_var", "Beta3_var", "Beta4_var"))
      fit_mu_and_var <- optim(nLL_mu_and_var, par = initial_pars, gr = grr, method = "BFGS")

      write(ifelse(fit_mu_and_var$convergence == 0, "true", "false"), "fit_converged.txt")

      write_tsv(enframe(fit_mu_and_var$par), "~{output_basename}_fitted_model_params.tsv")

      population_adjusted <- population_data %>% select(-residual_score2) %>% mutate(adjusted_score =
                                                                        (SCORE1_SUM - f_mu(population_data, fit_mu_and_var$par[1:5]))/
                                                                          sqrt(f_sigma2(population_data, fit_mu_and_var$par[6:10])))
      population_adjusted <- population_adjusted %>% mutate(percentile=pnorm(adjusted_score,0))

      write_tsv(population_adjusted, "population_adjusted_scores.tsv")
    EOF
  >>>

  output {
    File fitted_params = "~{output_basename}_fitted_model_params.tsv"
    File adjusted_population_scores = "population_adjusted_scores.tsv"
    Boolean fit_converged = read_boolean("fit_converged.txt")
  }

  runtime {
    docker: "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    disks: "local-disk 100 HDD"
    memory: mem + " GB"
  }
}

task AdjustScores {
  input {
    File fitted_model_params
    File pcs
    File scores
    Int mem = 2
  }

  command <<<
    Rscript -<< "EOF"
      library(dplyr)
      library(readr)

      # read in model params
      params_tibble <- read_tsv("~{fitted_model_params}")
      params <- params_tibble %>% pull(value)

      # linear transformation to predict variance
      f_sigma2 <- function(t, theta) {
        PC1 = t %>% pull(PC1)
        PC2 = t %>% pull(PC2)
        PC3 = t %>% pull(PC3)
        PC4 = t %>% pull(PC4)
        PC5 = t %>% pull(PC5)
        sigma2 <- exp(theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4)
      }


      # linear transformation to predict mean
      f_mu <- function(t, theta) {
        PC1 = t %>% pull(PC1)
        PC2 = t %>% pull(PC2)
        PC3 = t %>% pull(PC3)
        PC4 = t %>% pull(PC4)
        PC5 = t %>% pull(PC5)
        mu <- theta[[1]] + theta[[2]] * PC1 + theta[[3]] * PC2 + theta[[4]] * PC3 + theta[[5]] * PC4
      }

      scores = inner_join(read_tsv("~{pcs}"),
                                read_tsv("~{scores}"), by=c("IID" = "#IID"))

      adjusted_scores <- scores %>% mutate(adjusted_score =
                                                        (SCORE1_SUM - f_mu(scores, params[1:5]))/
                                                        sqrt(f_sigma2(scores, params[6:10]))
                                                      )
      adjusted_scores <- adjusted_scores %>% mutate(percentile=pnorm(adjusted_score,0))

      # return array scores
      write_tsv(adjusted_scores, "adjusted_scores.tsv")
    EOF
  >>>

  output {
    File adjusted_scores = "adjusted_scores.tsv"
  }

  runtime {
    docker: "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    disks: "local-disk 100 HDD"
    memory: mem + " GB"
  }
}

task MakePCAPlot {
  input {
    File population_pcs
    File target_pcs
  }

  command <<<
    Rscript -<< "EOF"
      library(dplyr)
      library(readr)
      library(ggplot2)

      population_pcs <- read_tsv("~{population_pcs}")
      target_pcs <- read_tsv("~{target_pcs}")

      ggplot(population_pcs, aes(x=PC1, y=PC2, color="Population")) +
        geom_point(size=0.1, alpha=0.1) +
        geom_point(data = target_pcs, aes(x=PC1, y=PC2, color="Target")) +
        labs(x="PC1", y="PC2") +
        theme_bw()

      ggsave(filename = "PCA_plot.png", dpi=300, width = 6, height = 6)

    EOF
  >>>

  output {
    File pca_plot = "PCA_plot.png"
  }

  runtime {
    docker: "rocker/tidyverse@sha256:0adaf2b74b0aa79dada2e829481fa63207d15cd73fc1d8afc37e36b03778f7e1"
    disks: "local-disk 100 HDD"
  }
}


task ExtractIDsPlink {
  input {
    File vcf
    Int disk_size = 2 * ceil(size(vcf, "GB")) + 100
    Int mem = 8
  }

  Int plink_mem = ceil(mem * 0.75 * 1000)

  command <<<
    /plink2 \
      --vcf ~{vcf} \
      --set-all-var-ids @:#:\$1:\$2 \
      --new-id-max-allele-len 1000 missing \
      --rm-dup exclude-all \
      --allow-extra-chr \
      --write-snplist allow-dups \
      --memory ~{plink_mem}
  >>>

  output {
    File ids = "plink2.snplist"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB"
  }
}

#plink chromosome encoding rules: https://www.cog-genomics.org/plink/2.0/data#irreg_output
task DetermineChromosomeEncoding {
  input {
    File weights
  }

  command <<<
    python3 << "EOF"
    code = 'MT'
    with open("~{weights}") as weights_file:
      chroms = {s.split("\t")[0].split(":")[0] for i, s in enumerate(weights_file) if i > 0}
      if any('chr' in c for c in chroms):
          if 'chrM' in chroms:
              code = 'chrM'
          else:
              code = 'chrMT'
      elif 'M' in chroms:
          code = 'M'

    with open('chr_encode_out.txt', 'w') as write_code_file:
        write_code_file.write(f'{code}\n')
    EOF
  >>>

  runtime {
    docker : "python:3.9.10"
  }

  output {
    String chromosome_encoding = read_string("chr_encode_out.txt")
  }
}