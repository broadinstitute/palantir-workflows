version 1.0

workflow DownsampleAndCollectCoverage {
    input {
        File input_cram
        File input_cram_index
        File ref_fasta
        File ref_fasta_index

        Float? target_coverage
        Float? downsample_probability

        File? coverage_intervals
        String downsample_strategy = "ConstantMemory"
        Int read_length = 150
        Boolean use_fast_algorithm = true

        String docker = "us.gcr.io/broad-gatk/gatk:4.4.0.0"
        File? picard_jar_override

        Int preemptible = 1
    }

    if (!defined(downsample_probability)) {
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

        String docker
        File? picard_jar_override

        Int preemptible
    }

    Int disk_size = ceil(size(input_cram, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")) + 50

    String output_basename = sub(sub(basename(input_cram), "\\.bam$", ""), "\\.cram$", "")
    String command_line = if defined(picard_jar_override) then "java -Xms2000m -Xmx2500m -jar " + picard_jar_override else 'gatk --java-options "-Xms2000m -Xmx2500m"'

    command <<<
        set -e -o pipefail

        ~{command_line} \
        CollectWgsMetrics \
        -INPUT ~{input_cram} \
        -VALIDATION_STRINGENCY SILENT \
        -REFERENCE_SEQUENCE ~{ref_fasta} \
        ~{"-INTERVALS " + coverage_intervals} \
        -OUTPUT ~{output_basename}.wgs_metrics \
        -READ_LENGTH ~{read_length} \
        -USE_FAST_ALGORITHM ~{use_fast_algorithm}

        cat ~{output_basename}.wgs_metrics | grep -v -e "^#" -e "^$" -e "^GENOME_TERRITORY" | head -n 1 | awk '{ print $2 }' > ~{output_basename}.mean_coverage
    >>>
    runtime {
        docker: docker
        preemptible: preemptible
        memory: "3000 MiB"
        disks: "local-disk " + disk_size + " HDD"
    }
    output {
        Float mean_coverage = read_float("~{output_basename}.mean_coverage")
    }
}

task Downsample {
    input{
        File input_cram

        Float? downsample_probability
        Float? target_coverage
        Float? original_coverage

        String downsample_strategy
        File ref_fasta
        File ref_fasta_index

        String docker
        File? picard_jar_override
        Int preemptible
        Int additional_disk_gb = 200
    }
    Int disk_size = ceil(3.5 * size(input_cram, "GiB") + 20 + additional_disk_gb)

    String output_basename = sub(sub(basename(input_cram), "\\.bam$", ""), "\\.cram$", "")
    String command_line = if defined(picard_jar_override) then "java -Xms10000m -Xmx10000m -jar " + picard_jar_override else 'gatk --java-options "-Xms10000m -Xmx10000m"'


    command {
        set -e -o pipefail

        PROBABILITY=~{if defined(downsample_probability) then downsample_probability else "$(bc -l <<< 'scale=2; " + target_coverage + "/" + original_coverage + "')"}
        
        ~{command_line} \
            DownsampleSam \
            -I ~{input_cram} \
            -O ~{output_basename}.downsampled.cram \
            -STRATEGY ~{downsample_strategy} \
            -P $PROBABILITY \
            -CREATE_INDEX true \
            -REFERENCE_SEQUENCE ~{ref_fasta}
            
    }

    runtime {
        docker: docker # TODO: update docker to use the new Picard options
        preemptible: preemptible
        memory: "14 GiB"
        cpu: "16"
        disks: "local-disk " + disk_size + " HDD"
    }
    output {
        File downsampled_cram = output_basename + ".downsampled.cram"
        File downsampled_cram_index = output_basename + ".downsampled.crai"
    }
}
