version 1.0

workflow CrosscheckVcfs {
    input {
        File vcf1
        File vcf1_index
        File vcf2
        File vcf2_index
        File haplotype_database
        String output_prefix = "crosscheck"
        Int memory_gb = 8
        Int disk_size_gb = 100
    }

    call CrosscheckFingerprints {
        input:
            vcf1 = vcf1,
            vcf1_index = vcf1_index,
            vcf2 = vcf2,
            vcf2_index = vcf2_index,
            haplotype_database = haplotype_database,
            output_prefix = output_prefix,
            memory_gb = memory_gb,
            disk_size_gb = disk_size_gb
    }

    output {
        File crosscheck_metrics = CrosscheckFingerprints.crosscheck_metrics
    }
}

task CrosscheckFingerprints {
    input {
        File vcf1
        File vcf1_index
        File vcf2
        File vcf2_index
        File haplotype_database
        String output_prefix
        Int memory_gb = 8
        Int disk_size_gb = 100
        Int preemptible = 1
    }

    Int command_mem_mb = memory_gb * 1000 - 500

    command <<<
        set -e

        # Run Picard CrosscheckFingerprints
        java -Xmx~{command_mem_mb}m -jar /usr/gitc/picard.jar CrosscheckFingerprints \
            INPUT=~{vcf1} \
            SECOND_INPUT=~{vcf2} \
            HAPLOTYPE_MAP=~{haplotype_database} \
            OUTPUT=~{output_prefix}.crosscheck_metrics.txt \
            CROSSCHECK_MODE = CHECK_SAME_SAMPLE \
            OUTPUT_ERRORS_ONLY=true \
            NUM_THREADS=4 

        
    >>>

    output {
        File crosscheck_metrics = "~{output_prefix}.crosscheck_metrics.txt"
        Boolean samples_match = read_boolean("match_status.txt")
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
        cpu: 4
    }
}
