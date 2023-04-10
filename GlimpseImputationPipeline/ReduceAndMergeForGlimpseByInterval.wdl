version 1.0

workflow ReduceAndMergeForGlimpseByInterval {
    input {
        Array[File] vcfs
        Array[File] vcf_indices
        File ref_fasta
        File ref_fasta_fai

        File subset_intervals_bed

        String output_basename_for_merged_vcf

        Int preemptible = 1
    }

    scatter(vcf_and_index in zip(vcfs, vcf_indices)) {
        call ReduceAnnotations {
            input:
                vcf = vcf_and_index.left,
                vcf_index = vcf_and_index.right,
                subset_intervals_bed = subset_intervals_bed,
                preemptible = preemptible
        }
    }

    call Merge {
        input:
            input_vcfs = ReduceAnnotations.reduced_vcf,
            input_vcf_indices = ReduceAnnotations.reduced_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            output_basename = output_basename_for_merged_vcf,
            preemptible = preemptible
    }

    output {
        File reduced_and_merged_vcf = Merge.merged_vcf
        File reduced_and_merged_vcf_index = Merge.merged_vcf_index
    }
}

task ReduceAnnotations {
    input {
        File vcf
        File vcf_index

        File subset_intervals_bed

        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int mem_gb = 8
        Int preemptible = 1
    }

    Int disk_size_gb = 3 * ceil(size(vcf, "GiB")) + 20

    String basename = basename(vcf, ".vcf.gz")

    command <<<
        set -e -o pipefail

        # Remove all annotations except GT and PL
        bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PL -O z -o ~{basename}.input_for_glimpse_temp.vcf.gz --threads 4 ~{vcf}
        bcftools index -t ~{basename}.input_for_glimpse_temp.vcf.gz
        
        # Intersect the output above with the reference panel sites to exclude any sites that are not in the panel
        bcftools view -R ~{subset_intervals_bed} -O z -o ~{basename}.input_for_glimpse.vcf.gz ~{basename}.input_for_glimpse_temp.vcf.gz
        bcftools index -t ~{basename}.input_for_glimpse.vcf.gz

    >>>

    runtime {
        docker: bcftools_docker
        cpu: 4
        memory: "~{mem_gb} GiB"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: preemptible
    }

    output {
        File reduced_vcf = "~{basename}.input_for_glimpse.vcf.gz"
        File reduced_vcf_index = "~{basename}.input_for_glimpse.vcf.gz.tbi"
    }
}

task Merge {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        File ref_fasta
        File ref_fasta_fai
        String output_basename
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int mem_gb = 8
        Int preemptible = 1
    }

    Int disk_size_gb = ceil(2.2 * size(input_vcfs, "GiB") + 20)

    command {
        set -e
        bcftools merge -o ~{output_basename}.unnormed.vcf.gz -O z --threads 4 ~{sep=' ' input_vcfs}
        bcftools index -t ~{output_basename}.unnormed.vcf.gz
        bcftools norm -O z -o ~{output_basename}.vcf.gz -f ~{ref_fasta} -m - ~{output_basename}.unnormed.vcf.gz
        bcftools index -t ~{output_basename}.vcf.gz
    }

    runtime {
        docker: bcftools_docker
        cpu: 2
        memory: "~{mem_gb} GiB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File merged_vcf = output_basename + ".vcf.gz"
        File merged_vcf_index = output_basename + ".vcf.gz.tbi"
    }
}
