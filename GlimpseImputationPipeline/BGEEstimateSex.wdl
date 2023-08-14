version 1.0

workflow BGEEstimateSex {
    input {
        File input_cram
        File input_cram_index
        File ref_fasta
        String docker = "us.gcr.io/broad-dsde-methods/samtools:v1"
        Int preemptible = 1
    }

    call BGEEstimateSexTask {
        input:
            input_cram = input_cram,
            input_cram_index = input_cram_index,
            ref_fasta = ref_fasta,
            docker = docker,
            preemptible = preemptible
    }

    output {
        Float sex_ratio = BGEEstimateSexTask.sex_ratio
    }
}

task BGEEstimateSexTask {
    input {
        File input_cram
        File input_cram_index
        File ref_fasta

        String docker
        Int preemptible
    }

    parameter_meta {
        input_cram: {
            localization_optional: true
        }
        input_cram_index: {
            localization_optional: true
        }
        ref_fasta: {
            localization_optional: true
        }
    }

    Int disk_size = 20

    command <<<
        set -xe -o pipefail

        export GCS_OAUTH_TOKEN=$(/root/google-cloud-sdk/bin/gcloud auth application-default print-access-token)
        
        reads_x=$(samtools view -h -T ~{ref_fasta} -F 0x0400 -q 20 -X ~{input_cram} ~{input_cram_index} chrX | samtools idxstats - | grep "chrX" | head -n 1 | awk '{ print $3 }')
        reads_y=$(samtools view -h -T ~{ref_fasta} -F 0x0400 -q 20 -X ~{input_cram} ~{input_cram_index} chrY | samtools idxstats - | grep "chrY" | head -n 1 | awk '{ print $3 }')

        python3 -c "print(${reads_y}/${reads_x})" > sex_ratio.txt
    >>>

    runtime {
        docker: docker
        preemptible: preemptible
        memory: "4 GiB"
        cpu: "2"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        Float sex_ratio = read_float("sex_ratio.txt")
    }
}