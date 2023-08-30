
version 1.0

workflow vcfdist_eval {
    input {
        File truth_vcf
        File eval_vcf
        File bed_file
        File fasta_file
        String docker = "us.gcr.io/broad-dsde-methods/vcfdist:v0.1"
    }

    call run_vcfdist_task {
        input:
            truth_vcf = truth_vcf,
            eval_vcf = eval_vcf,
            bed_file = bed_file,
            fasta_file = fasta_file,
            docker = docker
    }

    output {
        File prs_tsv = run_vcfdist_task.prs_tsv
        File vcfdistsummary = run_vcfdist_task.summary
        File precrec_tsv = run_vcfdist_task.precrec
        File query_tsv = run_vcfdist_task.querytsv
        File truth_tsv = run_vcfdist_task.truthtsv
        
    }
}

task run_vcfdist_task {
    input {
        File truth_vcf
        File eval_vcf
        File bed_file
        File fasta_file
        
        String docker
        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        /htslib-1.17/vcfdist/src/vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{fasta_file} \
            -b ~{bed_file} \
            -v 0
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File prs_tsv = "precision-recall-summary.tsv"
        File summary = "summary.vcf"
        File precrec = "precision-recall.tsv"
        File querytsv = "query.tsv"
        File truthtsv = "truth.tsv"
    }
}