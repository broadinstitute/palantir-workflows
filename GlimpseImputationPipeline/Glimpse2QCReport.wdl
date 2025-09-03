version 1.0

workflow Glimpse2QCReport {
    input {
        String cohort_name
        File metrics
        File vcf
        File vcf_index
        File pca_loadings
        File onnx_model
        File coverage_metrics
        File predicted_sex
    }

    call EstimateAncestry {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            pca_loadings = pca_loadings,
            onnx_model = onnx_model
    }

    call ExtractInfoScores {
        input:
            vcf = vcf
    }

    call Glimpse2QCReport_t {
        input:
            cohort_name = cohort_name,
            metrics = metrics,
            coverage_metrics = coverage_metrics,
            ancestries = EstimateAncestry.ancestries,
            predicted_sex = predicted_sex,
            info_scores = ExtractInfoScores.info_scores
    }


    output {
        File qc_report = Glimpse2QCReport_t.qc_report
    }
}

task ExtractInfoScores {
    input {
        String cohort_name
        File vcf
        Int disk_size_gb=100
        Int mem_gb=4
        Int cpu=1
    }

    command <<<
        set -euo pipefail
        echo "RAF\tINFO" | gzip > ~{cohort_name}_info_scores.tsv.gz
        bcftools query -f '%RAF\t%INFO/INFO\n' ~{vcf} | gzip >> ~{cohort_name}_info_scores.tsv.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.4"
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
      }

    output {
        File info_scores = "~{cohort_name}_info_scores.tsv.gz"   
    }
}

task EstimateAncestry {
    input {
        File vcf
        File vcf_index
        String cohort_name
        File pca_loadings
        File onnx_model
        Int mem_gb=32
        Int cpu=8
        Int n_splits=8
        Int disk_size_gb=100
    }

    command <<<
        set -euo pipefail

        # split bed into appropriate number of chunks
        awk -v OFS='\t' 'NR > 1 {print $1, $2-1, $2}' ~{pca_loadings} > loadings.bed
        n_line=$(wc -l < loadings.bed)
        line_per_split=$((n_line/~{n_splits}))
        split -d -l $line_per_split loadings.bed loadings_chunk_

        for f in loadings_chunk_*
        do 
            mv "$f" "$f".bed
            printf "contig\tpos\tref\talt\t%s\n" "$(bcftools query -l ~{vcf} | head -c -1 | tr "\n" "\t")" | gzip > dosages_"$f".tsv.gz
            bcftools +dosage -R "$f".bed ~{vcf} -- -t GT | tail -n +2 | gzip >> dosages_"$f".tsv.gz &
        done

        wait

        cat << EOF > predict_ancestry.py
        import pandas as pd
        import numpy as np
        import onnx
        import onnxruntime as ort
        import glob

        loadings = pd.read_csv('~{pca_loadings}', sep='\t')
        dosage_chunk_paths = glob.glob('dosages_loadings_chunk_*.tsv.gz')
        pcs = np.zeros((1,1))
        for path in dosage_chunk_paths:
            dosages = pd.read_csv(path, sep='\t')
            dosages = dosages.merge(loadings, on = ["contig","pos","ref","alt"], how="inner")
            pc_cols = [c for c in dosages.columns if c.startswith("pc") and c != "pca_af"]
            dosage_cols = [c for c in dosages.columns if c not in ['contig','pos','ref','alt','pca_af'] + pc_cols]
            pca_af_np = dosages['pca_af'].to_numpy().reshape(-1,1)
            n_var = loadings.shape[0]
            dosages_np = (dosages[dosage_cols].to_numpy() - 2*pca_af_np)/np.sqrt(n_var*2*pca_af_np*(1-pca_af_np))
            loadings_np = dosages[pc_cols].to_numpy()
            pcs = pcs + np.matmul(dosages_np.transpose(), loadings_np)
        
        ancestry_df = pd.DataFrame(pcs)
        ancestry_df.columns = pc_cols
        ancestry_df['sample'] = dosage_cols
        ancestry_df = ancestry_df[['sample'] + pc_cols].copy()

        with open('~{onnx_model}', 'rb') as f_onnx:
            onnx_model = onnx.load(f_onnx)
        
        sess = ort.InferenceSession(onnx_model.SerializeToString(), providers=["CPUExecutionProvider"])
        input_name = sess.get_inputs()[0].name
        label_name = sess.get_outputs()[0].name
        prob_name = sess.get_outputs()[1].name

        ancestry, _ = sess.run([label_name,prob_name], {input_name:pcs.astype(np.float32)})
        ancestry_df['ancestry'] = ancestry

        ancestry_df.to_csv('~{cohort_name}_ancestry.tsv', index=False, sep='\t')
        EOF

        python3 predict_ancestry.py
    >>>

    runtime {
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/bcftools-onnx:0.1"
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " HDD"
      }

    output {
        File ancestries = "~{cohort_name}_ancestry.tsv"   
    }
}

task Glimpse2QCReport_t {
    input {
        String cohort_name
        File metrics
        File coverage_metrics
        File ancestries
        File predicted_sex
        File info_scores
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
        ancestries = read_tsv("~{ancestries}")
        predicted_sex = read_tsv("~{predicted_sex}")

        qc_metrics = qc_metrics %>% inner_join(ancestries, by=c("sample_id"="sample")) %>% inner_join(predicted_sex)

        ancestry_counts <- ancestries %>% group_by(ancestry) %>% count() %>% arrange(-n)

        coverage_metrics<-read_tsv(c("~{coverage_metrics}"))
        coverage_metrics_per_sample <- coverage_metrics %>% group_by(Sample) %>% summarise(mean_cov=mean(\`Cov.\`),
                                                            mean_frac_sites = mean(1-\`No data pct\`/100))

        info_scores <- read_tsv("~{info_scores}")
        info_scores <- info_scores %>% mutate(raf_bin = 10**(ceiling((log10(RAF)- -4)/(4/10))*(4/10)+-4))

        info_mean_quantile <- info_scores %>% group_by(raf_bin) %>% summarise(q5 = quantile(INFO,0.05),
                                              q95 = quantile(INFO, 0.95),
                                              mean = mean(INFO))

        info_score_count <- bind_rows(info_scores %>% summarise(n=n(),
                         n_gt_0_6 = sum(INFO>0.6),
                         n_gt_0_8 = sum(INFO>0.8)) %>% mutate(raf_thresh = "all"),
                info_scores %>% filter(RAF>0.01) %>% summarise(n=n(),
                                                  n_gt_0_6 = sum(INFO>0.6),
                                                  n_gt_0_8 = sum(INFO>0.8)) %>%
                                mutate(raf_thresh = "raf > 1%"),
                info_scores %>% filter(RAF>0.001) %>% summarise(n=n(),
                                                     n_gt_0_6 = sum(INFO>0.6),
                                                     n_gt_0_8 = sum(INFO>0.8)) %>%
                                mutate(raf_thresh = "raf > 0.1%")   
            ) %>% relocate(raf_thresh)

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
        ggplot(info_scores %>% sample_n(10000), aes(x=RAF,y=INFO)) +
        geom_point(alpha=0.2) +
        geom_line(data=info_mean_quantile, aes(x=raf_bin, y=mean), color='red') +
        geom_ribbon(data=info_mean_quantile, aes(x=raf_bin, y=mean, ymin=q5, ymax=q95), fill='red', alpha=0.3) +
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
        docker: "us.gcr.io/broad-dsde-methods/r-markdown-pdf:0.1"
        disks: "local-disk 100 HDD"
        memory: mem_gb + " GB"
      }

    output {
        File qc_report = "~{cohort_name}_qc_report.pdf"
    }
}