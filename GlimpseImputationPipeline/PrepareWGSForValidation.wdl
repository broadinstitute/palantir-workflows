version 1.0

workflow PrepareWGSForValidation {
  input {
    Array[File] gvcfs
    Array[File] gvcf_indices

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Int scatter_count

    String output_basename
    File reference_panel
    File reference_panel_index
  }

  call SplitIntervalList {
    input:
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = 100
  }
  scatter (i in range(length(gvcfs))) {
    String gvcf = gvcfs[i]
    String gvcf_index = gvcf_indices[i]
    String gvcf_basename = basename(gvcf, ".g.vcf.gz")
    call Reblock {
      input:
        gvcf = gvcf,
        gvcf_index = gvcf_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        output_vcf_filename = gvcf_basename + ".reblocked.g.vcf.gz",
        docker_image = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    }

    call RemoveAFAnnotation_t {
      input:
        vcf = Reblock.output_vcf,
        basename = basename(Reblock.output_vcf,".g.vzf.gz") + ".AF.removed.g.vcf.gz"
    }
  }

  scatter(interval in SplitIntervalList.output_intervals) {
    call ImportGVCFs {
      input:
        vcfs = RemoveAFAnnotation_t.output_vcf,
        workspace_dir_name = "genomicsdb",
        interval = interval
    }

    call GnarlyGenotyper {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        output_vcf_filename = output_basename + ".g.vcf.gz",
        interval = interval
    }

    call UpdateAlternateAlleles {
      input:
        input_vcf = GnarlyGenotyper.output_vcf,
        input_vcf_index = GnarlyGenotyper.output_vcf_index,
        basename = output_basename + "_" + basename(interval, ".interval_list"),
        reference_panel = reference_panel,
        reference_panel_index = reference_panel_index
    }

    call RemoveAnnotations {
      input:
        vcf = UpdateAlternateAlleles.output_vcf,
        basename = output_basename + "_" + basename(interval, ".interval_list")
    }

    call SeparateMultiallelics {
      input:
        original_vcf = RemoveAnnotations.output_vcf,
        original_vcf_index = RemoveAnnotations.output_vcf_index,
        output_basename = output_basename + "_" + basename(interval, ".interval_list"),
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict
    }
  }

  call GatherVcfs {
    input:
      input_vcfs = SeparateMultiallelics.output_vcf,
      output_vcf_name = output_basename + ".vcf.gz",
      disk_size = 100 + ceil(2.2*size(SeparateMultiallelics.output_vcf, "GB"))
  }

  output {
    File output_vcf = GatherVcfs.output_vcf
    File output_vcf_index = GatherVcfs.output_vcf_index
  }

}

task SplitIntervalList {

  input {
    Int scatter_count
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  command <<<
    gatk --java-options -Xms3g SplitIntervals \
    -O  scatterDir -scatter ~{scatter_count} -R ~{ref_fasta}
  >>>

  runtime {
    memory: "3.75 GiB"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
  }
}


task ImportGVCFs {

  input {
    Array[File] vcfs
    File interval

    String workspace_dir_name

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    vcfs : {
             localization_optional: true
           }
  }

  command <<<
    set -euo pipefail

    rm -rf ~{workspace_dir_name}

    # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg
    # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    gatk --java-options -Xms8g \
    GenomicsDBImport \
    --genomicsdb-workspace-path ~{workspace_dir_name} \
    -V ~{sep=" -V " vcfs} \
    --merge-input-intervals \
    --consolidate \
    -L ~{interval} \


    tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}
  >>>

  runtime {
    memory: "26 GiB"
    cpu: 4
    bootDiskSizeGb: 15
    disks: "local-disk 200 HDD"
    docker: gatk_docker
  }

  output {
    File output_genomicsdb = "~{workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {

  input {
    File workspace_tar
    File interval

    String output_vcf_filename

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # This is needed for gVCFs generated with GATK3 HaplotypeCaller
    Boolean allow_old_rms_mapping_quality_annotation_data = false
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
  }

  Int disk_size = ceil(size(workspace_tar, "GiB") + size(ref_fasta, "GiB")) * 3 + 100

  parameter_meta {
    interval: {
                localization_optional: true
              }
  }

  command <<<
    set -euo pipefail

    tar -xf ~{workspace_tar}
    WORKSPACE=$(basename ~{workspace_tar} .tar)

    gatk --java-options -Xms8g \
    GenotypeGVCFs \
    -R ~{ref_fasta} \
    -O ~{output_vcf_filename} \
    --only-output-calls-starting-in-intervals \
    -V gendb://$WORKSPACE \
    -L ~{interval} \
    --merge-input-intervals \
    --include-non-variant-sites
  >>>

  runtime {
    memory: "26 GiB"
    cpu: 2
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

task GnarlyGenotyper {

  input {
    File workspace_tar
    String interval
    String output_vcf_filename
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    interval: {
                localization_optional: true
              }
  }

  Int disk_size = ceil(size(workspace_tar, "GiB") + size(ref_fasta, "GiB")) * 3 + 100

  command <<<
    set -e

    tar -xf ~{workspace_tar}
    WORKSPACE=$( basename ~{workspace_tar} .tar)

    gatk --java-options -Xms8g \
    GnarlyGenotyper \
    -R ~{ref_fasta} \
    -O ~{output_vcf_filename} \
    --only-output-calls-starting-in-intervals \
    -V gendb://$WORKSPACE \
    -stand-call-conf 10 \
    -L ~{interval} \
    --merge-input-intervals \
    --keep-all-sites
  >>>

  runtime {
    memory: "26 GiB"
    cpu: 2
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

task GatherVcfs {

  input {
    Array[File] input_vcfs
    String output_vcf_name
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    input_vcfs: {
                  localization_optional: true
                }
  }

  command <<<
    set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \
    GatherVcfsCloud \
    --ignore-safety-checks \
    --gather-type BLOCK \
    --input ~{sep=" --input " input_vcfs} \
    --output ~{output_vcf_name}

    tabix ~{output_vcf_name}
  >>>

  runtime {
    memory: "7 GiB"
    cpu: "1"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task UpdateAlternateAlleles {

  input {
    File input_vcf
    File input_vcf_index
    String basename
    File reference_panel
    File reference_panel_index

    File gatk_jar = "gs://broad-dsde-methods-ckachulis/jars/gatk_UpdateAltAlleles.jar"
  }

  Int disk_size = 3*ceil(size(input_vcf, "GB") + size(reference_panel, "GB")) + 20

  command <<<

    java -jar ~{gatk_jar} UpdateAltAlleles -V ~{input_vcf} -RP ~{reference_panel} -O ~{basename}.vcf.gz

  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/ckachulis/gatk-update-alt-alleles@sha256:1853bdf9cbed0f34e21ad5eddf0f4dd2e311aedb832835b08764a7861bec23a1"
    disks: "local-disk " + disk_size + " HDD"
    memory: "8 GB"

  }

  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"

  }
}

task SeparateMultiallelics {
  input {
    File original_vcf
    File original_vcf_index
    String output_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Int disk_size =  ceil(2*(size(original_vcf, "GB") + size(original_vcf_index, "GB")) + size(ref_fasta, "GB")) + 100
  }
  command {
    bcftools norm --fasta-ref ~{ref_fasta} -m - ~{original_vcf} -Oz -o ~{output_basename}.vcf.gz
    bcftools index -t ~{output_basename}.vcf.gz
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4 GB"
  }
}

task UpdateVariantIds {

  input {
    File vcf
    String basename
    Int disk_space =  3*ceil(size(vcf, "GB"))
  }

  command <<<
    bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' ~{vcf} -O z -o ~{basename}.vcf.gz
  >>>

  output {
    File output_vcf = "~{basename}.vcf.gz"
  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}

task SortVariantIds {
  input {
    File vcf
    String basename
    Int disk_space =  3*ceil(size(vcf, "GB"))
  }

  command <<<
    # what better way is there to do this I really don't know
    zcat ~{vcf} | awk -v OFS='\t' '{split($3, n, ":"); if ( !($1 ~ /^"#"/) && n[4] < n[3])  $3=n[1]":"n[2]":"n[4]":"n[3]; print $0}' | bgzip -c > ~{basename}.vcf.gz
    bcftools index -t ~{basename}.vcf.gz
  >>>

  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"

  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}

task RemoveAnnotations {
  input {
    File vcf
    String basename
  }

  Int disk_size = ceil(2.2*size(vcf, "GB")) + 100

  command <<<
    bcftools annotate ~{vcf} -x FORMAT,INFO -Oz -o ~{basename}.vcf.gz
    bcftools index -t ~{basename}.vcf.gz
  >>>

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"
  }
}

task SelectUnfilteredVariants {
  input {
    File vcf
    String basename
  }

  Int disk_size =  ceil(size(vcf, "GB")) + 50

  parameter_meta {
    vcf: {
           localization_optional : true
         }
  }

  command <<<
    gatk SelectVariants -V ~{vcf} --exclude-filtered -O ~{basename}.vcf.gz
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "16 GB"
  }

  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"
  }
}

task RemoveBadAnnotation {
  input {
    File vcf
    String basename
  }

  Int disk_size = ceil(2.2*size(vcf, "GB")) + 100

  command <<<
    bcftools annotate ~{vcf} -x FORMAT/AF,FORMAT/DP,FORMAT/F1R2,FORMAT/F2R1,FORMAT/ICNT,FORMAT/MB,FORMAT/PRI,FORMAT/SB,FORMAT/SPL -Oz -o ~{basename}.vcf.gz
    bcftools index -t ~{basename}.vcf.gz
  >>>

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"
  }
}

task Reblock {

  input {
    File gvcf
    File gvcf_index
    String output_vcf_filename
    String docker_image
    File ref_fasta
    File ref_fasta_index
    File ref_dict
  }

  Int disk_size = ceil(size(gvcf, "GiB")) * 2

  command {
    gatk --java-options "-Xms3g -Xmx3g" \
    ReblockGVCF \
    -V ~{gvcf} \
    --reference ~{ref_fasta} \
    -drop-low-quals \
    -do-qual-approx \
    --floor-blocks -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 \
    -O ~{output_vcf_filename}
  }

  runtime {
    memory: "3.75 GB"
    bootDiskSizeGb: "15"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
    docker: docker_image
  }

  output {
    File output_vcf = output_vcf_filename
    File output_vcf_index = output_vcf_filename + ".tbi"
  }
}

task RemoveAFAnnotation_t {
  input {
    File vcf
    String basename
  }

  Int disk_size = ceil(2.2*size(vcf, "GB")) + 100

  command <<<
    bcftools annotate ~{vcf} -x FORMAT/AF,FORMAT/DP,FORMAT/F1R2,FORMAT/F2R1,FORMAT/ICNT,FORMAT/MB,FORMAT/PRI,FORMAT/SB,FORMAT/SPL -Oz -o ~{basename}.vcf.gz
    bcftools index -t ~{basename}.vcf.gz
  >>>

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"
  }
}



