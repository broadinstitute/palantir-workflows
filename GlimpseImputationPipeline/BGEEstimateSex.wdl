version 1.0

workflow BGEEstimateSex {
    input {
        File input_cram
        String docker = "us.gcr.io/broad-dsde-methods/samtools:v1"
        Int preemptible = 1
    }

    call BGEEstimateSexTask {
        input:
            input_cram = input_cram,
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

        String docker
        Int preemptible
    }

    parameter_meta {
        input_cram: {
            localization_optional: true
        }
    }

    Int disk_size = 20

    command <<<
        set -e -o pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        
        cov_x=$(samtools view -h ~{input_cram} chrX | samtools idxstats - | grep "chrX" | head -n 1 | awk '{printf "%.5f\n", $3/$2}')
        cov_y=$(samtools view -h ~{input_cram} chrY | samtools idxstats - | grep "chrY" | head -n 1 | awk '{printf "%.5f\n", $3/$2}')

        python3 -c "print(${cov_y}/${cov_x})" > sex_ratio.txt
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