version 1.0

workflow Glimpse2QCReport {
    input {
        String cohort_name
        File metrics
        Array[File] coverage_metrics
        File ancestries
        File predicted_sex
        File info_score_qc
        File info_mean_quantile
        File info_score_sampling
    }
    call Glimpse2QCReport_t {
        input:
            cohort_name = cohort_name,
            metrics = metrics,
            coverage_metrics = coverage_metrics,
            ancestries = ancestries,
            predicted_sex = predicted_sex,
            info_score_qc = info_score_qc,
            info_mean_quantile = info_mean_quantile,
            info_score_sampling = info_score_sampling
    }


    output {
        File qc_report = Glimpse2QCReport_t.qc_report
    }
}

task Glimpse2QCReport_t {
    input {
        String cohort_name
        File metrics
        Array[File] coverage_metrics
        File ancestries
        File predicted_sex
        File info_score_qc
        File info_mean_quantile
        File info_score_sampling
        Int mem_gb=4
    }

    command <<<
        set -eo pipefail

        cat << EOF > ~{cohort_name}_qc_report.Rmd
        ---
        title: GLIMPSE Imputation QC Report
        subtitle: ~{cohort_name}
        author: "Broad Clinical Labs"
        date: "\`r Sys.Date()\`"
        output: pdf_document
        ---

        \`\`\`{r setup, include=FALSE}
        library(dplyr)
        library(readr)
        library(ggplot2)
        library(knitr)
        library(scales)
        \`\`\`

        \`\`\`{r load data, include=FALSE,echo=FALSE, message=FALSE, warning=FALSE}
        qc_metrics = read_tsv("~{metrics}")
        ancestries = read_csv("~{ancestries}")
        predicted_sex = read_tsv("~{predicted_sex}")

        qc_metrics = qc_metrics %>% inner_join(ancestries, by=c("sample_id"="s")) %>% inner_join(predicted_sex)

        ancestry_counts <- ancestries %>% group_by(ancestry) %>% count() %>% arrange(-n)

        coverage_metrics<-read_tsv(c("~{sep='","' coverage_metrics}"))
        coverage_metrics_per_sample <- coverage_metrics %>% group_by(Sample) %>% summarise(mean_cov=mean(\`Cov.\`),
                                                            mean_frac_sites = mean(1-\`No data pct\`/100))

        info_score_count <- read_tsv("~{info_score_qc}") %>% arrange(-total_sites)

        info_mean_quantile <- read_tsv("~{info_mean_quantile}")

        info_sample <- read_tsv("~{info_score_sampling}")

        \`\`\`
        # Site QC
        ## Site Coverage QC
        During GLIMPSE imputation, the genome is split into 837 chunks, with the majority of the algorithm being applied to each chunk separately.  All the chunks are then ligated back together at the end of the process.  For each chunk, the mean coverage of sites being genotyped is calculated for each sample, along with the fractions of sites covered by at least one read.

        In the histograms below, we show the distributions of mean coverage and fractions of sites covered over all chunk/sample combinations.


        \`\`\`{r mean_cov, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(coverage_metrics, aes(x=\`Cov.\`)) +
        geom_histogram() + theme_bw() + xlim(0,20) +
        xlab("Mean Coverage (per sample x chunk)")
        \`\`\`

        \`\`\`{r frac_covered, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(coverage_metrics, aes(x=1-\`No data pct\`/100)) +
        geom_histogram() + theme_bw() + xlim(0,1) +
        xlab("Fraction of Sites Covered by >=1 Read (per sample x chunk)")
        \`\`\`

        Below we show the same plots, but with the values for each sample averaged across all chunks.

        \`\`\`{r mean_cov_sample, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(coverage_metrics_per_sample, aes(x=mean_cov)) +
        geom_histogram() + theme_bw() + xlim(0,20) +
        xlab("Mean Coverage (per sample)")
        \`\`\`

        \`\`\`{r frac_covered_sample, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(coverage_metrics_per_sample, aes(x=mean_frac_sites)) +
        geom_histogram() + theme_bw() + xlim(0,1) +
        xlab("Fraction of Sites Covered by >=1 Read (per sample)")
        \`\`\`


        ## Info Score QC
        For each genotyped site, GLIMPSE calculates an INFO score, which is it\`s estimate of how well it has imputed the site across all the samples in the cohort.  It is standard practice during downstream analyses using imputed datasets to filter sites according the an INFO score threshold, often with the INFO score required to be greater than either 0.6 or 0.8.  

        The total number and fraction of sites passing passing these thresholds is show in the table below.  The rows in the table indicate various reference panel allele frequency thresholds, since more common sites are generally imputed much more accurately.

        \`\`\`{r info_score_table, echo=FALSE, message=FALSE, warning=FALSE}
        kable(info_score_count, format.args = list(big.mark = ','), digits = 2, col.names = c("","Total Sites", "INFO > 0.6", "INFO > 0.8", 
                                            "Fraction with INFO > 0.6", "Fraction with INFO > 0.8"))
        \`\`\`

        Below we show the relationship between reference panel allele frequency and INFO score.  The line and shaded regions represent the mean and 5%-95% quantiles for various reference panel allele frequency bins.  The points are a random sampling of 10,000 sites.

        \`\`\`{r info_score_plot, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(info_sample, aes(x=RAF,y=INFO)) +
        geom_point(alpha=0.2) +
        geom_line(data=info_mean_quantile, aes(x=raf, y=mean), color='red') +
        geom_ribbon(data=info_mean_quantile, aes(x=raf, y=mean, ymin=fifth_pctile, ymax=ninety_fifth_pctile), fill='red', alpha=0.3) +
        scale_x_log10() + theme_bw() + xlab("Reference Panel Allele Frequency")
        \`\`\`

        # Sample QC

        ## Genetic Ancestry Estimation
        As part of sample QC, genetic ancestry was estimated for all samples, in order to confirm that any stratification of QC metrics are expected based on genetic ancestry.  Genetic ancestry was estimated  by PCA using resources provided by [gnomad](https://gnomad.broadinstitute.org/downloads#v3-ancestry-classification).  Note that these ancestry estimates are used solely for QC, and are not intended to be used for any other purpose.

        \`\`\`{r pca, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(ancestries, aes(x=PC1, y=PC2)) +
        geom_point(aes(color=ancestry)) + theme_bw()
        \`\`\`

        \`\`\`{r ancestry_table, echo=FALSE, message=FALSE, warning=FALSE}
        kable(ancestry_counts) 
        \`\`\`

        ## QC Metrics
        Sample QC Metrics distributions are shown below.  Note that some metrics are stratified by reported sex in addition to by ancestry due to expected differences in the non-PAR region of Chromosome X.  These metrics are calculated using [hail sample_qc](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc).  Note that for genotype counts, each site is counted at most once per sample.  However, metrics like \`Number of SNP Alt Alleles\` count the number of alternate alleles observed in a sample, in which case a heterozygous genotype would count as one, and a homozygous alternate genotype as two. 

        ### Number of Homozygous Reference Genotypes
        \`\`\`{r n_hom_ref, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(qc_metrics, aes(x=n_hom_ref, fill=ancestry)) +
        geom_histogram() + theme_bw() +
        xlab("Number of Homozygous Reference Genotypes") +
        scale_x_continuous(labels=comma)
        \`\`\`

        ### Number of Heterozygous Genotypes
        \`\`\`{r n_het, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(qc_metrics %>% filter(predicted_sex %in% c("F","M")), aes(x=n_het, fill=ancestry)) +
        geom_histogram() + theme_bw() + facet_grid(~predicted_sex) +
        xlab("Number of Heterozygous Genotypes") +
        scale_x_continuous(labels=comma) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        \`\`\`


        ### Number of Homozygous Variant Genotypes
        \`\`\`{r n_hom_var, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(qc_metrics %>% filter(predicted_sex %in% c("F","M")), aes(x=n_hom_var, fill=ancestry)) +
        geom_histogram() + theme_bw() + facet_grid(~predicted_sex)  +
        xlab("Number of Homozygous Variant Genotypes") +
        scale_x_continuous(labels=comma) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
        \`\`\`

        ### Ratio of Heterozygous to Homozygous Variant Genotypes
        \`\`\`{r r_het_hom_var, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(qc_metrics %>% filter(predicted_sex %in% c("F","M")), aes(x=r_het_hom_var, fill=ancestry)) +
        geom_histogram() + theme_bw() + facet_grid(~predicted_sex) +
        xlab("Ratio of Heterozygous to Homozygous Variant Genotypes") +
        scale_x_continuous(labels=comma)
        \`\`\`

        ### Number of Non-Reference Genotypes
        \`\`\`{r n_non_ref, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(qc_metrics, aes(x=n_non_ref, fill=ancestry)) +
        geom_histogram() + theme_bw() +
        xlab("Number of Non-Reference Genotypes") +
        scale_x_continuous(labels=comma)
        \`\`\`

        ### Number of SNP Alt Alleles 
        \`\`\`{r n_snps, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(qc_metrics, aes(x=n_snp, fill=ancestry)) +
        geom_histogram() + theme_bw() +
        xlab("Number of SNP Alt Alleles") +
        scale_x_continuous(labels=comma)
        \`\`\`

        ### Number of Insertion Alt Alleles
        \`\`\`{r n_ins, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(qc_metrics, aes(x=n_insertion, fill=ancestry)) +
        geom_histogram() + theme_bw() +
        xlab("Number of Insertion Alt Alleles") +
        scale_x_continuous(labels=comma)
        \`\`\`

        ### Number of Deletion Alt Alleles
        \`\`\`{r n_del, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(qc_metrics, aes(x=n_deletion, fill=ancestry)) +
        geom_histogram() + theme_bw() +
        xlab("Number of Deletion Alt Alleles") +
        scale_x_continuous(labels=comma)
        \`\`\`

        ### ti/tv Ratio
        \`\`\`{r r_ti_tv, echo=FALSE, message=FALSE, warning=FALSE}
        ggplot(qc_metrics, aes(x=r_ti_tv, fill=ancestry)) +
        geom_histogram() + theme_bw() +
        xlab("ti/tv Ratio") +
        scale_x_continuous(labels=comma)
        \`\`\`

        EOF

        Rscript -e "library(rmarkdown); rmarkdown::render('~{cohort_name}_qc_report.Rmd', 'pdf_document')"

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/r-gcnv-viz@sha256:c92a9a26cba4ab10b579776cd4dd70d5fca485d4a7216d0ab9b477662e3d9228"
        disks: "local-disk 100 HDD"
        memory: mem_gb + " GB"
      }

    output {
        File qc_report = "~{cohort_name}_qc_report.pdf"
    }
}