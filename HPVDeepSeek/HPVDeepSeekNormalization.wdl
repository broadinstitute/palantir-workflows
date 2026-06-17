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
        python <<CODE

        import pysam
        import pandas as pd
        from collections import Counter

        infile_simplex = pysam.AlignmentFile("~{simplex_bam}", "rb")
        df = pd.read_csv("~{fp_intervals}", sep = '\t', header = None, names = ["chromosome", "start", "end", "info"])
        df_detected_hpv_genotypes = pd.read_csv("~{hpv_status}", sep = '\t')

        chroms_and_lengths = dict(zip(infile_simplex.references, infile_simplex.lengths))
        chroms_and_lengths_hpv = {k: v for k, v in chrom_lengths.items() if k.startswith("HPV")}

        new_rows = []
        for key, value in chroms_and_lengths_hpv.items():
            new_rows.append({"chromosome": key, "start": 0, "end": value, "info": "N/A"})
        df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index = True)

        for idx, row in df.iterrows():
            total_depth = 0
            num_positions = 0
            for pileupcolumn in infile_simplex.pileup(row.chromosome, row.start, row.end, stepper = "all", truncate = False, max_depth = 1000000, ignore_overlaps = True):
                for pileupread in pileupcolumn.pileups:
                    if pileupread.alignment.get_tag("cD") >= 5:
                        total_depth += 1
                num_positions += 1

            mean_depth = 0.0
            if num_positions > 0:
                mean_depth = total_depth / num_positions
            df.loc[idx, "mean_depth"] = mean_depth

        hg38_median_depth = df.loc[~df["chromosome"].str.startswith("HPV") & ~df["chromosome"].str.startswith("chrX") & ~df["chromosome"].str.startswith("chrY"), "mean_depth"].median()

        mask = df["chromosome"].isin(df_detected_hpv_genotypes["HPV_Genotype"].tolist())
        df.loc[mask, "r"] = df.loc[mask, "mean_depth"].apply(lambda x: x / hg38_median_depth)
        df.loc[mask, "hpv_quantity"] = df.loc[mask, "r"].apply(lambda x: x * ~{ng_cfdna} / 0.0033 / ~{ul_plasma} / 1000.0)

        outfile.write("HPV_Genotype" + "\t" + "HPV_Mean_Depth_Over_hg38_Median_Depth" + "\t" + "ng_cfDNA" + "\t" + "mL_Plasma" + "\t" + "HPV_Quantity" + "\n")
        for row in df.loc[mask].itertuples():
            outfile.write(row.chromosome + "\t" + str(row.r) + "\t" + str(ng_cfdna) + "\t"+ str(ul_plasma / 1000.0) + "\t" + str(row.hpv_quantity) + "\n")

        infile_simplex.close()
        outfile.close()

        CODE
    >>>

    output {
        File normalized_hpv = "~{sample_id}.normalized_hpv.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:1abe088592a1c82c93e6899e03f5aecf1651fc7ecfaaf947417b6ab9b4706884"
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