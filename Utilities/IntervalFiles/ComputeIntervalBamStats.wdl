version 1.0

# For computing stats/distributions over various interval files for BAM files
workflow ComputeIntervalBamStats {
    input {
        String input_name
        File input_bam
        File input_bam_index
        String experiment = ""

        File ref_fasta
        File ref_index
        File ref_dict

        Array[File] interval_files
        Array[String] interval_names
    }

    scatter(interval in zip(interval_names, interval_files)) {
        call SubsetBam {
            input:
                input_bam=input_bam,
                input_bam_index=input_bam_index,
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                ref_dict=ref_dict,
                interval_file=interval.right,
                interval_name=interval.left
        }

        call ComputeWgsMetrics {
            input:
                input_bam=SubsetBam.bam_subset,
                input_bam_index=SubsetBam.bam_subset_index,
                ref_fasta=ref_fasta,
                ref_index=ref_index,
                ref_dict=ref_dict,
                interval_file=interval.right,
                interval_name=interval.left
        }

        call ComputeMAPQDistribution {
            input:
                input_bam=SubsetBam.bam_subset,
                input_bam_index=SubsetBam.bam_subset_index,
                interval_file=interval.right,
                interval_name=interval.left
        }
    }

    call CollectCoverageResults {
        input:
            input_name=input_name,
            experiment=experiment,
            wgs_files=ComputeWgsMetrics.wgs_file
    }

    call CollectMAPQResults {
        input:
            input_name=input_name,
            experiment=experiment,
            mapq_files=ComputeMAPQDistribution.mapq_dist
    }

    output {
        File wgs_summary = CollectCoverageResults.wgs_summary
        File cov_summary = CollectCoverageResults.cov_summary
        File mapq_summary = CollectMAPQResults.mapq_summary
    }
}

task SubsetBam {
    input {
        File input_bam
        File input_bam_index

        File ref_fasta
        File ref_index
        File ref_dict

        File interval_file
        String interval_name

        String gatk_tag = "4.3.0.0"
        Int disk_size = 2 * ceil(size(input_bam, "GB")) + 25
        Int cpu = 8
        Int memory = 32
    }

    parameter_meta {
        input_bam: {
           localization_optional: true
       }
    }

    command <<<
        set -xueo pipefail

        gatk PrintReads \
            -I ~{input_bam} \
            -L ~{interval_file} \
            -O "~{interval_name}.bam"

    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + " GB"
    }

    output {
        File bam_subset = "~{interval_name}.bam"
        File bam_subset_index = "~{interval_name}.bai"
    }
}

# For coverage-related stats and distribution
task ComputeWgsMetrics {
    input {
        File input_bam
        File input_bam_index

        File ref_fasta
        File ref_index
        File ref_dict

        File interval_file
        String interval_name

        # Tool arguments
        Int min_BQ = 20
        Int min_MQ = 20
        Int coverage_cap = 250

        # Runtime
        String gatk_tag = "4.3.0.0"
        Int disk_size = 2 * ceil(size(input_bam, "GB")) + 25
        Int cpu = 8
        Int memory = 32
    }

    command <<<
        set -xueo pipefail

        gatk CollectWgsMetrics \
            -I ~{input_bam} \
            -R ~{ref_fasta} \
            -O "~{interval_name}.wgs_metrics" \
            --INTERVALS ~{interval_file} \
            ~{"--MINIMUM_BASE_QUALITY " + min_BQ} \
            ~{"--MINIMUM_MAPPING_QUALITY " + min_MQ} \
            ~{"--COVERAGE_CAP " + coverage_cap}

    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:" + gatk_tag
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + " GB"
    }

    output {
        File wgs_file = "~{interval_name}.wgs_metrics"
    }
}

task ComputeMAPQDistribution {
    input {
        File input_bam
        File input_bam_index

        File interval_file
        String interval_name

        # Runtime
        Int cpu = 8
        Int memory = 16
        Int extra_disk = 25

    }

    Int disk_size = 2*ceil(size(input_bam, "GB")) + extra_disk

    command <<<
        set -xueo pipefail

        samtools view ~{input_bam} | awk '{print $5}' | sort -n | uniq -c > "~{interval_name}_mapq_dist.txt"

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/samtools:v1"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + " GB"
    }

    output {
        File mapq_dist = "~{interval_name}_mapq_dist.txt"
    }
}

task CollectMAPQResults {
    input {
        String input_name
        String experiment
        Array[File] mapq_files

        Int disk_size = 50
        Int cpu = 2
        Int memory = 8
    }

    command <<<
        python3 << CODE

        import os
        import pandas as pd

        file_paths = ["~{default="" sep="\", \"" mapq_files}"]
        full_mapq_summary = pd.DataFrame()

        for file in file_paths:
            interval_name = os.path.basename(file).split('_mapq_dist.txt')[0]

            mapq_summary = pd.read_csv(file, sep='\s+', names=['Count', 'MAPQ'])
            mapq_summary['Interval_List'] = interval_name
            mapq_summary['Sample'] = "~{input_name}"
            mapq_summary['Experiment'] = "~{experiment}"
            full_mapq_summary = pd.concat([full_mapq_summary, mapq_summary])

        full_mapq_summary.to_csv("mapq_summary.tsv", sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + "GB"
    }

    output {
        File mapq_summary = "mapq_summary.tsv"
    }
}


task CollectCoverageResults {
    input {
        String input_name
        String experiment
        Array[File] wgs_files

        Int disk_size = 50
        Int cpu = 2
        Int memory = 8
    }

    command <<<
        python3 << CODE

        import os
        import pandas as pd

        file_paths = ["~{default="" sep="\", \"" wgs_files}"]
        full_wgs_summary = pd.DataFrame()
        full_cov_summary = pd.DataFrame()

        for file in file_paths:
            interval_name = os.path.basename(file).split('.wgs_metrics')[0]

            wgs_summary = pd.read_csv(file, sep='\t', comment='#', nrows=1)
            wgs_summary['Interval_List'] = interval_name
            wgs_summary['Sample'] = "~{input_name}"
            wgs_summary['Experiment'] = "~{experiment}"
            full_wgs_summary = pd.concat([full_wgs_summary, wgs_summary])

            cov_dist = pd.read_csv(file, sep='\t', skiprows=10)
            cov_dist['Interval_List'] = interval_name
            cov_dist['Sample'] = "~{input_name}"
            cov_dist['Experiment'] = "~{experiment}"
            full_cov_summary = pd.concat([full_cov_summary, cov_dist])

        full_wgs_summary.to_csv("wgs_summary.tsv", sep='\t', index=False)
        full_cov_summary.to_csv("cov_summary.tsv", sep='\t', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        memory: memory + "GB"
    }

    output {
        File wgs_summary = "wgs_summary.tsv"
        File cov_summary = "cov_summary.tsv"
    }
}
