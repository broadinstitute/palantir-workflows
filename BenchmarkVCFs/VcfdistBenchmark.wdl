version 1.0

workflow VcfdistEval {
    input {
        File truth_vcf
        File eval_vcf
        File bed_file
        File fasta_file
        String docker = "us.gcr.io/broad-dsde-methods/vcfdist:v0.1"
        Int? verbosity
    }

    call RunVcfdistTask {
        input:
            truth_vcf = truth_vcf,
            eval_vcf = eval_vcf,
            bed_file = bed_file,
            fasta_file = fasta_file,
            docker = docker,
            verbosity = select_first([0,verbosity])
    }

    output {
        File prs_tsv = RunVcfdistTask.prs_tsv
        File vcfdistsummary = RunVcfdistTask.summary
        File precrec_tsv = RunVcfdistTask.precrec
        File query_tsv = RunVcfdistTask.query_tsv
        File truth_tsv = RunVcfdistTask.truth_tsv
        
    }
}

task RunVcfdistTask {
    input {
        File truth_vcf
        File eval_vcf
        File bed_file
        File fasta_file
        Int verbosity
        
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
            -v ~{verbosity}
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
        File query_tsv = "query.tsv"
        File truth_tsv = "truth.tsv"
    }
}