version 1.0

workflow BenchmarkPhasing {
    input {
        File base_vcf
        File base_vcf_index
        String? base_sample_name

        File call_vcf
        File call_vcf_index
        String? call_sample_name

        String? experiment

        Array[File]? interval_beds
        Array[String]? interval_bed_labels
    }

    call WhatsHapCompare {
        input:
            base_vcf=base_vcf,
            base_vcf_index=base_vcf_index,
            base_sample_name=base_sample_name,
            call_vcf=call_vcf,
            call_vcf_index=call_vcf_index,
            call_sample_name=call_sample_name,
            experiment=experiment
    }

    if (defined(interval_beds)) {
        call EvaluateRegions {
            input:
                switch_errors_bed=WhatsHapCompare.switch_errors_bed,
                call_sample_name=call_sample_name,
                experiment=experiment,
                interval_beds=select_first([interval_beds]),
                interval_bed_labels=select_first([interval_bed_labels])
        }
    }

    output {
        File whatshap_eval = WhatsHapCompare.eval
        File switch_errors_bed = WhatsHapCompare.switch_errors_bed
        File longest_blocks = WhatsHapCompare.longest_blocks
        File jaccard_table = select_first([EvaluateRegions.jaccard_table, ""])
    }
}

task WhatsHapCompare {
    input {
        File base_vcf
        File base_vcf_index
        String? base_sample_name

        File call_vcf
        File call_vcf_index
        String? call_sample_name

        String? experiment

        Boolean only_snvs = false    # Toggle to only consider phasing for SNVs

        # Runtime arguments
        Int disk_size = 2 * ceil(size(base_vcf, "GB") + size(call_vcf, "GB")) + 50
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xue

        whatshap compare \
            --names ~{default="truth" base_sample_name},~{default="callset" call_sample_name} \
            --tsv-pairwise eval.tsv \
            --ignore-sample-name \
            ~{true="--only-snvs" false="" only_snvs} \
            --switch-error-bed switch_errors.bed \
            --longest-block-tsv longest_blocks.tsv \
            ~{base_vcf} ~{call_vcf}

        echo "Experiment" > exp_column.txt
        yes ~{experiment} | head -n $(cat eval.tsv | wc -l) | sed -e '1d' >> exp_column.txt
        paste exp_column.txt eval.tsv > final-eval.tsv

        echo "Experiment" > exp_column.txt
        yes ~{experiment} | head -n $(cat longest_blocks.tsv | wc -l) | sed -e '1d' >> exp_column.txt
        paste exp_column.txt longest_blocks.tsv > final-longest_blocks.tsv
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/whatshap:v1"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory_ram + " GB"
    }

    output {
        File eval = "final-eval.tsv"
        File switch_errors_bed = "switch_errors.bed"
        File longest_blocks = "final-longest_blocks.tsv"
    }
}

task EvaluateRegions {
    input {
        File switch_errors_bed

        String? call_sample_name
        String? experiment

        Array[File] interval_beds
        Array[String] interval_bed_labels

        # Runtime arguments
        Int disk_size = 2 * ceil(size(switch_errors_bed, "GB") + size(interval_beds, "GB")) + 50
        Int cpu = 4
        Int memory_ram = 16
    }

    command <<<
        set -xue

        INTERVAL_BEDS=(~{sep=' ' interval_beds})
        INTERVAL_LABELS=(~{sep=' ' interval_bed_labels})

        # For each interval bed, we want to compute the proportion of switch errors occuring in that region.
        # This is achieved by first intersecting the two bed files, and then using bedtools to compute the
        # Jaccard index of the two (compute the size of the intersection over the full switch error file).
        # Finally add a column with the interval label in it for downstream data analysis.
        for i in {0..~{length(interval_beds)-1}}
        do
            bedtools intersect -a "${INTERVAL_BEDS[$i]}" -b ~{switch_errors_bed} > intersect.bed
            bedtools sort -i intersect.bed > intersect-sorted.bed
            bedtools jaccard -a intersect-sorted.bed -b ~{switch_errors_bed} > jaccard.tsv
            tail -n 1 jaccard.tsv > contents.tsv

            echo "${INTERVAL_LABELS[$i]}" > names.txt
            paste names.txt contents.tsv >> full-jaccard.tsv
        done

        # Add header back
        echo -e "Interval_Name\tintersection\tunion\tjaccard\tn_intersections" > header.txt
        cat header.txt full-jaccard.tsv > header-jaccard.tsv

        # Add constant columns
        echo "Experiment" > exp_column.txt
        yes ~{experiment} | head -n $(cat header-jaccard.tsv | wc -l) | sed -e '1d' >> exp_column.txt
        paste exp_column.txt header-jaccard.tsv > jaccard-no_sample.tsv

        echo "Sample" > sample_column.txt
        yes ~{call_sample_name} | head -n $(cat jaccard-no_sample.tsv | wc -l) | sed -e '1d' >> sample_column.txt
        paste sample_column.txt jaccard-no_sample.tsv > final-jaccard.tsv
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bedtools:v1.0"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory_ram + " GB"
    }

    output {
        File jaccard_table = "final-jaccard.tsv"
    }
}
