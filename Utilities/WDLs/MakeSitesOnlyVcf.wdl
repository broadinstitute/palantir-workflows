version 1.0

  workflow MakeSitesOnlyVcfAndLiftover {
    input {
      File vcf
      File chain
      File ref_fasta
      File ref_dict
      File ref_index
    }

    call MakeSitesOnlyVcf {
      input:
        vcf = vcf
    }

    call LiftoverVcf {
      input:
        vcf = MakeSitesOnlyVcf.sites_only_vcf,
        vcf_index = MakeSitesOnlyVcf.sites_only_vcf_index,
        chain = chain,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict
    }
  }

  task MakeSitesOnlyVcf {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf.gz")

    Int disk_size = 100 + 2*ceil(size(vcf, "GB"))


    command <<<
      gatk MakeSitesOnlyVcf -I ~{vcf} -O ~{basename}.sites_only.vcf.gz
  >>>

    runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.0.0"
      disks: "local-disk " + disk_size + " HDD"
      bootDiskSizeGb: "16"
    }

    output {
      File sites_only_vcf = "~{basename}.sites_only.vcf.gz"
      File sites_only_vcf_index = "~{basename}.sites_only.vcf.gz.tbi"
    }
  }

  task LiftoverVcf {
    input {
      File vcf
      File vcf_index
      File chain
      File ref_fasta
      File ref_dict
      File ref_index

      Int memory_jvm=12
    }

    Int memory_ram = memory_jvm + 4

    String basename = basename(vcf, ".vcf.gz")
    Int disk_size = 100 + ceil(2*size(vcf, "GB")+size(ref_fasta, "GB"))

    command <<<
      gatk --java-options "-Xmx~{memory_jvm}G" LiftoverVcf -I ~{vcf} -O ~{basename}.lifted_over.vcf.gz -C ~{chain} --REJECT ~{basename}.reject.vcf.gz -R ~{ref_fasta}
  >>>

    runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
      disks: "local-disk " + disk_size + " HDD"
      memory: memory_ram + " GB"
    }

    output {
      File lifted_over_vcf = "~{basename}.lifted_over.vcf.gz"
      File lifter_over_vcf_index = "~{basename}.lifted_over.vcf.gz.tbi"
    }
  }