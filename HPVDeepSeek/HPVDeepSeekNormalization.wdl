version 1.0

task NormalizeHPV {
    input {
        String sample_id
        File simplex_bam
        File simplex_bam_index
        File hpv_status
        File fp_intervals
        Float ul_plasma
        Float ng_cfdna

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = 128
    }

    command <<<
        set -e
        python3 <<CODE

        import pysam
        import pandas as pd
        from collections import Counter

        df = pd.read_csv("~{fp_intervals}", sep = '\t', header = None, names = ["chromosome", "start", "end", "info"])
        df_detected_hpv_genotypes = pd.read_csv("~{hpv_status}", sep = '\t')

        with pysam.AlignmentFile("~{simplex_bam}", "rb") as infile_simplex:
            chroms_and_lengths = dict(zip(infile_simplex.references, infile_simplex.lengths))
            chroms_and_lengths_hpv = {k: v for k, v in chroms_and_lengths.items() if k.startswith("HPV")}

            new_rows = []
            for key, value in chroms_and_lengths_hpv.items():
                new_rows.append({"chromosome": key, "start": 0, "end": value, "info": "N/A"})
            df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index = True)

            for idx, row in df.iterrows():
                total_depth = 0
                for pileupcolumn in infile_simplex.pileup(row.chromosome, row.start, row.end, stepper = "all", truncate = True, max_depth = 1000000, ignore_overlaps = True):
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.alignment.get_tag("cD") >= 5:
                            total_depth += 1

                num_positions = row.end - row.start
                mean_depth = 0.0
                if num_positions > 0:
                    mean_depth = total_depth / num_positions
                df.loc[idx, "mean_depth"] = mean_depth

        hg38_median_depth = df.loc[~df["chromosome"].str.startswith("HPV") & ~df["chromosome"].str.startswith("chrX") & ~df["chromosome"].str.startswith("chrY"), "mean_depth"].median()

        df = df[df["chromosome"].isin(df_detected_hpv_genotypes["HPV_Genotype"].tolist())]
        df = df.rename(columns = {"mean_depth": "HPV_Mean_Depth"})
        df['hg38_median_depth'] = hg38_median_depth
        df["HPV_Mean_Depth_Over_hg38_Median_Depth"] = df["HPV_Mean_Depth"] / hg38_median_depth
        df["ng_cfDNA"] = ~{ng_cfdna}
        df["mL_Plasma"] = ~{ul_plasma} / 1000.0
        df["HPV_Quantity"] = df["HPV_Mean_Depth_Over_hg38_Median_Depth"] * ((df["ng_cfDNA"] / 0.0033) / df["mL_Plasma"])

        df = df.rename(columns = {"chromosome": "HPV_Genotype"})
        df = df[["HPV_Genotype", "HPV_Mean_Depth", "hg38_median_depth", "HPV_Mean_Depth_Over_hg38_Median_Depth", "ng_cfDNA", "mL_Plasma", "HPV_Quantity"]]
        df.to_csv("~{sample_id}.normalized_hpv.tsv", sep = '\t', index = False)

        CODE
    >>>

    output {
        File normalized_hpv = "~{sample_id}.normalized_hpv.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:f7f71cf1996056c32a4c7ad5ef6a855093383ab5c919a3539874e56a85539256"
    }
}

workflow HPVDeepSeekNormalization {
    input {
        String sample_id
        File simplex_bam
        File simplex_bam_index
        File hpv_status
        File fp_intervals
        Float ul_plasma
        Float ng_cfdna
    }

    call NormalizeHPV {
        input:
            sample_id = sample_id,
            simplex_bam = simplex_bam,
            simplex_bam_index = simplex_bam_index,
            hpv_status = hpv_status,
            fp_intervals = fp_intervals,
            ul_plasma = ul_plasma,
            ng_cfdna = ng_cfdna
    }

    output {
        File normalized_hpv = NormalizeHPV.normalized_hpv
    }
}