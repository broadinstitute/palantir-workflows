version 1.0
## This is assuming that your variant IDs are in the format chr:positionr:allele1:allele2 and just alphabetizes allele1 and
## allele 2. To get variant IDs in this format, you can run `bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ~{vcf}`
task SortIds {
  input {
    File vcf
    String basename
    Int disk_space =  3*ceil(size(vcf, "GB"))
  }

  command <<<
    zcat ~{vcf} | awk -v OFS='\t' '{split($3, n, ":"); if ( !($1 ~ /^"#"/) && n[4] < n[3])  $3=n[1]":"n[2]":"n[4]":"n[3]; print $0}' | bgzip -c > ~{basename}.vcf.gz
  >>>

  output {
    File output_vcf = "~{basename}.vcf.gz"
  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}

# score with plink2
task ScoreVcf {
  input {
    File vcf
    String basename
    File weights
    Int base_mem = 8
    String? extra_args
    File? sites
  }

  Int runtime_mem = base_mem + 2
  Int plink_mem = base_mem * 1000
  Int disk_space =  3*ceil(size(vcf, "GB")) + 20

  command {
    /plink2 --score ~{weights} header ignore-dup-ids list-variants-zs no-mean-imputation \
    cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,scoresums --allow-extra-chr ~{extra_args} -vcf ~{vcf} dosage=DS \
    ~{"--extract " + sites} --out ~{basename} --memory ~{plink_mem}
  }

  output {
    File score = "~{basename}.sscore"
    File log = "~{basename}.log"
    File sites_scored = "~{basename}.sscore.vars.zst"
  }

  runtime {
    docker: "skwalker/plink2:first"
    disks: "local-disk " + disk_space + " HDD"
    memory: runtime_mem + " GB"
  }
}

# This just turns the array vcf into bim/bed/fam format and extracts only the sites
# Used for the original PCA steps
task ArrayVcfToPlinkDataset {
  input {
    File vcf
    File pruning_sites
    File? subset_to_sites
    String basename
    Int mem = 8
  }

  Int disk_space =  3*ceil(size(vcf, "GB")) + 20

  command {

    /plink2 --vcf ~{vcf} --extract-intersect ~{pruning_sites} ~{subset_to_sites} --allow-extra-chr \
    --out ~{basename} --make-bed --rm-dup force-first
  }

  output {
    File bed = "~{basename}.bed"
    File bim = "~{basename}.bim"
    File fam = "~{basename}.fam"
  }

  runtime {
    docker: "skwalker/plink2:first"
    disks: "local-disk " + disk_space + " HDD"
    memory: mem + " GB"
  }
}

# This projects the array dataset using the previously generated PCs, using flashPCA
task ProjectArray {
  input {
    File bim
    File bed
    File fam
    File pc_loadings
    File pc_meansd
    String basename
    Int mem = 8
  }

  command {

    cp ~{bim} ~{basename}.bim
    cp ~{bed} ~{basename}.bed
    cp ~{fam} ~{basename}.fam

    cp ~{pc_loadings} loadings.txt
    cp ~{pc_meansd} meansd.txt

    ~/flashpca/flashpca --bfile ~{basename} --project --inmeansd meansd.txt \
    --outproj projections.txt --inload loadings.txt -v
  }

  output {
    File projections = "projections.txt"
  }

  runtime {
    docker: "quay.io/ckachuli/flashpca@sha256:85e9ee91bc552e46a0d69cc851b893419c8de6588c696458fc770eee526e381d" # a special version of flashpca which allows to project a single sample without erroring out at an unnecessary check
    disks: "local-disk 400 HDD"
    memory: mem + " GB"
  }
}

# This does the scoring adjustment
task AdjustScores {
  input {
    File population_pcs
    File population_scores
    File array_pcs
    File array_scores
    Int mem = 2
  }

  command <<<
    Rscript - <<- "EOF"
    library(ggplot2)
    library(dplyr)

    population_pcs = read.csv("~{population_pcs}", sep="\t", header = T)
    population_scores = read.csv("~{population_scores}", sep="\t", header = T)

    population_data = merge(population_pcs, population_scores, by.x="IID", by.y="X.IID")

    # generate the linear model from the population data using the first 5 PCs
    population_model = glm(SCORE1_SUM ~ PC1 + PC2 + PC3 + PC4 + PC5, data = population_data, family = "gaussian")

    population_data$residual_score = resid(population_model)
    population_resid_mean = mean(population_data$residual_score)
    population_resid_sd = sd(population_data$residual_score)

    # calculate adjusted score on population data,  make sure it's standardized to N(0, 1)
    population_data$adjusted_score = (population_data$residual_score - population_resid_mean)/population_resid_sd

    # this calculates the adjusted score for the new data
    generate_adjusted_scores = function(new_data) {
    subset_data_for_model = new_data %>% transmute(raw_score = SCORE1_SUM, PC1=PC1, PC2=PC2, PC3=PC3, PC4=PC4, PC5=PC5)
    new_data$residual_score = subset_data_for_model$raw_score - predict(population_model, subset_data_for_model) # calculate the residual score from the model
    new_data$adjusted_score = (new_data$residual_score - population_resid_mean)/population_resid_sd # again adjust compared to population
    new_data %>% rowwise() %>% mutate(percentile=pnorm(adjusted_score, round(mean(population_data$adjusted_score),5),
    round(sd(population_data$adjusted_score), 5) == 1))
    }

    array_scores = merge(read.csv("~{array_pcs}",  sep = "\t", header = T),
    read.csv("~{array_scores}",  sep = "\t", header = T), by.x="IID", by.y="X.IID")

    adjusted_array_scores = generate_adjusted_scores(array_scores)

    # make sure the PCs fit well between the array and the population data
    ggplot(population_data, aes(x=PC1, y=PC2, color="Population Data")) + geom_point() + geom_point() +
    geom_point(data = array_scores, aes(x=PC1, y=PC2, color="Array Data")) + labs(x="PC1", y="PC2") + theme_bw()
    ggsave(filename = "PCA_plot.png", dpi=300, width = 6, height = 6)

    # return population scores
    write.table(population_data %>% subset(select = -residual_score), file = "population_data_scores.tsv", sep="\t", row.names=F, quote = F)

    # return array scores
    write.table(adjusted_array_scores %>% subset(select = -residual_score), file = "array_data_scores.tsv", sep="\t", row.names=F, quote = F)
    EOF
  >>>

  output {
    File pca_plot = "PCA_plot.png"
    File adjusted_population_scores = "population_data_scores.tsv"
    File adjusted_array_scores = "array_data_scores.tsv"
  }

  runtime {
    docker: "skwalker/rscripting:with_rutils"
    disks: "local-disk 400 HDD"
    memory: mem + " GB"
  }
}


# This task allows you to sort each variant ID in your weights file. It already assumes they are in the format chr:position:a1:a2
# and just sorts a1, a2. You will have to perform other awk magic to get it into this format otherwise.

task SortWeights {
  input {
    File weights_file
    Int disk_space = 50
    Int id_column # the column # of the variant IDs
    String basename # what you wanted the new weights file to be called
  }

  command <<<

    awk -v id_col="~{id_column}" -v OFS='\t' '{split($id_col, n, ":"); if ( n[4] < n[3])  $id_col=n[1]":"n[2]":"n[4]":"n[3]; print $0}' ~{weights_file} > ~{basename}.txt

  >>>

  output {
    File sorted_weights = "~{basename}.txt"
  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}

task UpdateVariantIds {
  input {
    File? vcf
    String basename
    Int disk_space =  3*ceil(size(vcf, "GB"))
  }

  command <<<
    bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' ~{vcf} -O z -o ~{basename}.vcf.gz
  >>>

  output {
    File output_vcf = "~{basename}.vcf.gz"
  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}

task ExtractIDs {
  input {
    File vcf
    String output_basename
    Int disk_size = 2*ceil(size(vcf, "GB")) + 100
  }

  command <<<
    bcftools query -f "%ID\n" ~{vcf} -o ~{output_basename}.original_array.ids
  >>>
  output {
    File ids = "~{output_basename}.original_array.ids"
  }
  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4 GB"
  }
}

task PerformPCA {
  input {
    File bim
    File bed
    File fam
    String basename
    Int mem = 8
  }

  # again, based on Wallace commands
  command {
    cp ~{bim} ~{basename}.bim
    cp ~{bed} ~{basename}.bed
    cp ~{fam} ~{basename}.fam

    ~/flashpca/flashpca --bfile ~{basename} -n 16 -d 20 --outpc ${basename}.pc \
    --outpve ${basename}.pc.variance --outload ${basename}.pc.loadings \
    --outmeansd ${basename}.pc.meansd
  }

  output {
    File pcs = "${basename}.pc"
    File pc_variance = "${basename}.pc.variance"
    File pc_loadings = "${basename}.pc.loadings"
    File mean_sd = "${basename}.pc.meansd"
    File eigenvectors = "eigenvectors.txt"
    File eigenvalues = "eigenvalues.txt"
  }

  runtime {
    docker: "skwalker/flashpca:v1"
    disks: "local-disk 400 HDD"
    memory: mem + " GB"
  }
}

#Print given message to stderr and return an error
task ErrorWithMessage{
  input {
    String message
  }
  command <<<
    >&2 echo "Error: ~{message}"
    exit 1
  >>>

  runtime {
    docker: "ubuntu"
  }
}
