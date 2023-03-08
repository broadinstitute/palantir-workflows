version 1.0

workflow ReduceAndMergeForGlimpse {
    input {
        Array[File] vcfs
        Array[File] vcf_indices

        File? reference_panel_sites_vcf
        File? reference_panel_sites_vcf_index

        String output_basename_for_merged_vcf

        Int preemptible = 1
    }

    scatter(vcf_and_index in zip(vcfs, vcf_indices)) {
        call ReduceAnnotations {
            input:
                vcf = vcf_and_index.left,
                vcf_index = vcf_and_index.right,
                reference_panel_sites_vcf = reference_panel_sites_vcf,
                reference_panel_sites_vcf_index = reference_panel_sites_vcf_index,
                preemptible = preemptible
        }
    }

    call Merge {
        input:
            input_vcfs = ReduceAnnotations.reduced_vcf,
            input_vcf_indices = ReduceAnnotations.reduced_vcf_index,
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

        File? reference_panel_sites_vcf
        File? reference_panel_sites_vcf_index

        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int mem_gb = 8
        Int preemptible = 1
    }

    Int disk_size_gb = 3 * ceil(size(vcf, "GiB")) + 20

    String basename = basename(vcf, ".vcf.gz")

    command <<<
        set -e -o pipefail

        bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PL -O z -o ~{basename}.input_for_glimpse.vcf.gz --threads 4 ~{vcf}
        bcftools index -t ~{basename}.input_for_glimpse.vcf.gz
        # This line will only be executed if reference_panel_sites_vcf is defined
        ~{"bcftools isec -p isec_output -n =2 -w 1 -O z " + basename + ".input_for_glimpse.vcf.gz " + reference_panel_sites_vcf + " && mv isec_output/0000.vcf.gz " + basename + ".input_for_glimpse.vcf.gz && mv isec_output/0000.vcf.gz.tbi " + basename + ".input_for_glimpse.vcf.gz.tbi"}
    >>>

    runtime {
        docker: bcftools_docker
        cpu: 1
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
        String output_basename
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int mem_gb = 8
        Int preemptible = 1
    }

    Int disk_size_gb = ceil(2.2 * size(input_vcfs, "GiB") + 20)

    command {
        set -e
        bcftools merge -o ~{output_basename}.vcf.gz -O z --threads 4 ~{sep=' ' input_vcfs}
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
