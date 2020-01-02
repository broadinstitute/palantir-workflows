version 1.0

import "../../CoverageAnalysis/GetUndercoveredRegions.wdl"
workflow testGetUndercoveredRegions {

    input {
        Array[File] sam_files
        File expected_interval_list
    }

    scatter (sam in sam_files) {
        call SamToBam{
            input:
                sam = sam
        }
    }

    call GetUndercoveredRegions.GetUndercoveredRegions {
        input:
            bam_files = SamToBam.bam,
            bam_indicies = SamToBam.bam_index
    }

    call AssertPassed {
        input:
            expected_interval_list = expected_interval_list,
            observed_interval_list = GetUndercoveredRegions.interval_list_output
    }
}

task SamToBam {
    input {
        File sam
    }

    String bamName = basename(sam, ".sam") + ".bam"
    command <<<
        samtools view -h ~{sam} -o ~{bamName}
        samtools index ~{bamName}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
    }

    output {
        File bam = "~{bamName}"
        File bam_index = "~{bamName}.bai"
    }
}

task AssertPassed {
    input {
        File expected_interval_list
        File observed_interval_list
    }

    command <<<
        bases_expected=$(java -jar /usr/gitc/picard.jar IntervalListTools I=~{expected_interval_list} OUTPUT_VALUE=BASES)
        bases_observed=$(java -jar /usr/gitc/picard.jar IntervalListTools I=~{observed_interval_list} OUTPUT_VALUE=BASES)
        bases_intersect=$(java -jar /usr/gitc/picard.jar IntervalListTools ACTION=INTERSECT I=~{observed_interval_list} I=~{expected_interval_list} OUTPUT_VALUE=BASES)

        if [[ $bases_expected != $bases_observed || $bases_expected != $bases_intersect ]]; then
            >&2 echo "observed interval list ~{observed_interval_list} differs from expected interval list ~{expected_interval_list}"
            exit 1
        fi
    >>>

    runtime {
    	docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
    }
}