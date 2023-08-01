version 1.0

workflow DownsampleAndCollectCoverage {
    input {
        File input_cram
        File input_cram_index
        File ref_fasta
        File ref_fasta_index

        Float? target_coverage
        Float? downsample_probability
        Float? fail_if_below_coverage

        File? coverage_intervals
        String downsample_strategy = "ConstantMemory"
        Int read_length = 150
        Boolean use_fast_algorithm = true
        Boolean output_bam_instead_of_cram = false

        String docker = "us.gcr.io/broad-gatk/gatk:4.4.0.0"
        File? picard_jar_override

        Int preemptible = 1
    }

    # If downsample_probability is defined, use that. Otherwise, use target_coverage.
    # If target_coverage isn't defined either, skip this task too. The Downsample task
    # contains the input validation and will output an error if neither are defined.
    if (!defined(downsample_probability) && defined(target_coverage)) {
        call CollectWgsMetrics as CollectOriginalCoverage {
            input:
                input_cram = input_cram,
                input_cram_index = input_cram_index,
                coverage_intervals = coverage_intervals,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                read_length = read_length,
                use_fast_algorithm = use_fast_algorithm,
                docker = docker,
                picard_jar_override = picard_jar_override,
                preemptible = preemptible
        }
    }

    call Downsample {
        input:
            input_cram = input_cram,
            downsample_probability = downsample_probability,
            target_coverage = target_coverage,
            original_coverage = CollectOriginalCoverage.mean_coverage,
            downsample_strategy = downsample_strategy,
            output_bam_instead_of_cram = output_bam_instead_of_cram,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            docker = docker,
            picard_jar_override = picard_jar_override,
            preemptible = preemptible
    }

    call CollectWgsMetrics as CollectDownsampledCoverage {
        input:
            input_cram = Downsample.downsampled_cram,
            input_cram_index = Downsample.downsampled_cram_index,
            coverage_intervals = coverage_intervals,
            fail_if_below_coverage = fail_if_below_coverage,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            read_length = read_length,
            use_fast_algorithm = use_fast_algorithm,
            docker = docker,
            picard_jar_override = picard_jar_override,
            preemptible = preemptible
    }

    output {
        File downsampled_cram = Downsample.downsampled_cram
        File downsampled_cram_index = Downsample.downsampled_cram_index
        Float downsampled_mean_coverage = CollectDownsampledCoverage.mean_coverage
        File downsampled_wgs_metrics = CollectDownsampledCoverage.wgs_metrics
        Float? original_mean_coverage = CollectOriginalCoverage.mean_coverage
    }
}

task CollectWgsMetrics {
    input {
        File input_cram
        File input_cram_index
        File? coverage_intervals
        File ref_fasta
        File ref_fasta_index
        Int read_length
        Boolean use_fast_algorithm
        Float? fail_if_below_coverage

        String docker
        File? picard_jar_override

        Int preemptible
    }

    Int disk_size = ceil(size(input_cram, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")) + 50
    String output_basename = sub(sub(basename(input_cram), "\\.bam$", ""), "\\.cram$", "")

    command <<<
        set -e -o pipefail

        ~{if defined(picard_jar_override) then "java -Xms2000m -Xmx2500m -jar " + picard_jar_override else 'gatk --java-options "-Xms2000m -Xmx2500m"'} \
        CollectWgsMetrics \
        -INPUT ~{input_cram} \
        -VALIDATION_STRINGENCY SILENT \
        -REFERENCE_SEQUENCE ~{ref_fasta} \
        ~{"-INTERVALS " + coverage_intervals} \
        -OUTPUT ~{output_basename}.wgs_metrics \
        -READ_LENGTH ~{read_length} \
        -USE_FAST_ALGORITHM ~{use_fast_algorithm} \
        -COUNT_UNPAIRED true

        grep -v -e '^#' -e "^$" "~{output_basename}.wgs_metrics" | head -2 > wgs.tsv
        COL_NUM=$(sed 's/\t/\n/g' wgs.tsv | grep -n 'MEAN_COVERAGE' | cut -d':' -f1)
        awk -v col=$COL_NUM ' { print $col }' wgs.tsv | tail -1 > "~{output_basename}.mean_coverage"

        mean_coverage=$(cat ~{output_basename}.mean_coverage)
        ~{"if (( $(echo '$mean_coverage < " + fail_if_below_coverage + "' |bc -l) )); then echo -e '\nERROR: Downsampled coverage below minimum threshold.\n'; exit 1; fi"}
    >>>

    runtime {
        docker: docker
        preemptible: preemptible
        memory: "3000 MiB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        Float mean_coverage = read_float("~{output_basename}.mean_coverage")
        File wgs_metrics = "~{output_basename}.wgs_metrics"
    }
}

task Downsample {
    input{
        File input_cram

        Float? downsample_probability
        Float? target_coverage
        Float? original_coverage

        String downsample_strategy
        Boolean output_bam_instead_of_cram
        File ref_fasta
        File ref_fasta_index

        String docker
        File? picard_jar_override
        Int preemptible
        Int additional_disk_gb = 200
    }

    Int disk_size = ceil(3.5 * size(input_cram, "GiB") + 20 + additional_disk_gb)
    String output_basename = sub(sub(basename(input_cram), "\\.bam$", ""), "\\.cram$", "")
    String output_extension = if output_bam_instead_of_cram then "bam" else "cram"
    String output_index_extension = if output_bam_instead_of_cram then "bam.bai" else "cram.crai"

    command <<<
        set -e -o pipefail

        ~{if !defined(downsample_probability) && !defined(target_coverage) then "echo -e '\nERROR: Must define either downsample_probability or target_coverage.\n'; exit 1" else ""}

        PROBABILITY=~{if defined(downsample_probability) then downsample_probability else (if select_first([target_coverage, 0]) > select_first([original_coverage, 0]) then "1" else "$(bc -l <<< 'scale=2; " + target_coverage + "/" + original_coverage + "')")}
        
        ~{if defined(picard_jar_override) then "java -Xms2000m -Xmx2500m -jar " + picard_jar_override else 'gatk --java-options "-Xms2000m -Xmx2500m"'} \
            DownsampleSam \
            -I ~{input_cram} \
            -O ~{output_basename}.downsampled.~{output_extension} \
            -STRATEGY ~{downsample_strategy} \
            -P $PROBABILITY \
            -CREATE_INDEX false \
            -REFERENCE_SEQUENCE ~{ref_fasta}
        
        samtools index ~{output_basename}.downsampled.~{output_extension}
    >>>

    runtime {
        docker: docker
        preemptible: preemptible
        memory: "14 GiB"
        cpu: "16"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File downsampled_cram = output_basename + ".downsampled." +output_extension
        File downsampled_cram_index = output_basename + ".downsampled." + output_index_extension
    }
}
