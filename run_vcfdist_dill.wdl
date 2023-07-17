version 1.0

workflow vcfdist_eval {
    input {
        File truth_vcf
        File eval_vcf
        File bed_file
        File fasta_file
        String docker = "us.gcr.io/broad-dsde-methods/vcfdist:v0.1"
    }

    call run_vcfdist_task{
        input:
            truth_vcf = truth_vcf,
            eval_vcf = eval_vcf,
            bed_file = bed_file,
            fasta_file = fasta_file,
            docker = docker
    }

    output {
        File prs_tsv = run_vcfdist_task.prs_tsv
    }
}

task run_vcfdist_task{
    input {
        File truth_vcf
        File eval_vcf
        File bed_file
        File fasta_file
        
        String docker
        Int disk_size_gb = ceil(1 * size(truth_vcf, "GiB") + 10)
        Int mem_gb = 2
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<

        docker build -t "ubuntu:vcfdist_dill" .
        docker run -i -t ubuntu:vcfdist_dill /bin/bash
        cd htslib-1.17/vcfdist/src
        ./vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{fasta_file} \
            -b ~{bed_file} \
            -v 0
        exit
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File prs_tsv = "htslib-1.17/vcfdist/src/precision-recall-summary.tsv"
    }
}