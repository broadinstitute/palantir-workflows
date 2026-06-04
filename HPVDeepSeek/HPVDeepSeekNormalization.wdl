version 1.0

task CalculateMeanDepthsSimplex {
    input {
        String sample_id
        File simplex_bam
        File simplex_bam_index
        File target_intervals

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    command <<<
        set -e
        python3 <<CODE

        import pysam
        from collections import Counter

        infile_simplex = pysam.AlignmentFile("~{simplex_bam}", "rb")

        target_interval_list = []
        with open("~{target_intervals}", 'r') as f:
            target_interval_list = [line.strip() for line in f]

        hg38_target_mean_depths = open("~{sample_id}.hg38_target_mean_depths.tsv", 'w')
        hpv_target_mean_depths_fs_geq_5 = open("~{sample_id}.hpv_target_mean_depths_fs_geq_5.tsv", 'w')

        for target in target_interval_list:
            tokens = target.split('\t')
            chromosome = tokens[0]
            start = int(tokens[1])
            end = int(tokens[2])

            total_depth = 0
            num_positions = 0

            for pileupcolumn in infile_simplex.pileup(chromosome, start, end, stepper = "all", truncate = False, max_depth = 1000000, ignore_overlaps = True):
                if chromosome.startswith("HPV"):
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.alignment.get_tag("cD") >= 5:
                            total_depth += 1
                else:
                    total_depth += pileupcolumn.nsegments
                num_positions += 1

            mean_depth = 0.0
            if num_positions > 0:
                mean_depth = total_depth / num_positions

            if chromosome.startswith("HPV"):
                hpv_target_mean_depths_fs_geq_5.write(chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + str(mean_depth) + "\n")
            else:
                hg38_target_mean_depths.write(chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + str(mean_depth) + "\n")

        hpv_target_mean_depths_fs_geq_5.close()
        hg38_target_mean_depths.close()

        CODE
    >>>

    output {
        File hg38_target_mean_depths = "~{sample_id}.hg38_target_mean_depths.tsv"
        File hpv_target_mean_depths_fs_geq_5 = "~{sample_id}.hpv_target_mean_depths_fs_geq_5.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:a302f9efe0bf1d4f9998ee1e9dda406223454ccaea0b5619046742221c1d2a74"
    }
}

task GetMedianOfHg38MeanDepthsSimplex {
    input {
        String sample_id
        File hg38_target_mean_depths
        File fp_intervals

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    command <<<
        set -e
        python3 <<CODE

        import statistics

        fp_interval_list = []
        with open("~{fp_intervals}", 'r') as f:
            fp_interval_list = [line.strip() for line in f]

        hg38_non_xy_fp_mean_depths = []

        with open("~{hg38_target_mean_depths}", 'r') as f:
            for line in f:
                line = line.rstrip()
                columns = line.split('\t')

                chromosome = tokens[0]
                start = int(tokens[1])
                end = int(tokens[2])
                mean_depth = float(tokens[3])

                key = chromosome + ":" + str(start) + "-" + str(end)

                if (key in fp_interval_list) and (not key.startswith("chrX")) and (not key.startswith("chrY")):
                    hg38_non_xy_simplex_mean_depths.append(mean_depth)

        median_hg38 = statistics.median(hg38_non_xy_fp_mean_depths)

        outfile = open("~{sample_id}.median_hg38.tsv", 'w')
        outfile.write(str(median_hg38) + "\n")
        outfile.close()

        CODE
    >>>

    output {
        File median_hg38 = "~{sample_id}.median_hg38.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:a302f9efe0bf1d4f9998ee1e9dda406223454ccaea0b5619046742221c1d2a74"
    }
}

task NormalizeHPV {
    input {
        String sample_id
        File hpv_status
        File hpv_target_mean_depths_fs_geq_5
        File median_hg38
        Float ml_plasma
        Float ng_cfdna

        Int cpu = 2
        Int memory_gb = 16
        Int disk_size_gb = 128
    }

    command <<<
        set -e
        python3 <<CODE

        detected_hpv_genotypes = []
        with open("~{hpv_status}", 'r') as f:
            header = f.readline()
            for line in f:
                line = line.rstrip()
                columns = line.split('\t')

                is_detected = columns[3].strip().lower() == "true"
                if is_detected:
                    detected_hpv_genotypes.append(columns[0])

        median_hg38_val = 0.0
        with open(~{median_hg38}", 'r') as f:
            median_hg38_val = float(f.readline().strip())

        outfile = open("~{sample_id}.normalized_hpv.tsv", 'w')

        with open("~{hpv_target_mean_depths_fs_geq_5}", 'r') as f:
            for line in f:
                line = line.rstrip()
                columns = line.split('\t')

                if columns[0] not in detected_hpv_genotypes:
                    continue

                normalized_hpv_quantity = float(columns[3]) / median_hg38_val
                normalized_hpv_quantity_ml_plasma = normalized_hpv_quantity / ~{ml_plasma}
                normalized_hpv_quantity_ng_cfdna = normalized_hpv_quantity / ~{ng_cfdna}

                outfile.write(columns[0] + "\t" + str(normalized_hpv_quantity) + "\t" + str(normalized_hpv_quantity_ml_plasma) + "\t" + str(normalized_hpv_quantity_ng_cfdna) + "\n")

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
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

workflow HPVDeepSeekNormalization {
    input {
        String sample_id
        File simplex_bam
        File simplex_bam_index
        File hpv_status
        File target_intervals
        File fp_intervals
        Float ml_plasma
        Float ng_cfdna
    }

    call CalculateMeanDepthsSimplex {
        input:
            sample_id = sample_id,
            simplex_bam = simplex_bam,
            simplex_bam_index = simplex_bam_index,
            target_intervals = target_intervals
    }

    call GetMedianOfHg38MeanDepthsSimplex {
        input:
            sample_id = sample_id,
            hg38_target_mean_depths = CalculateMeanDepthsSimplex.hg38_target_mean_depths,
            fp_intervals = fp_intervals

    }

    call NormalizeHPV {
        input:
            sample_id = sample_id,
            hpv_status = hpv_status,
            hpv_target_mean_depths_fs_geq_5 = CalculateMeanDepthsSimplex.hpv_target_mean_depths_fs_geq_5,
            median_hg38 = GetMedianOfHg38MeanDepthsSimplex.median_hg38,
            ml_plasma = ml_plasma,
            ng_cfdna = ng_cfdna
    }

    output {
        File normalized_hpv = NormalizeHPV.normalized_hpv
    }
}