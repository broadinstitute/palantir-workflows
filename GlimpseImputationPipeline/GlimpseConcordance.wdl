version 1.0

workflow GlimpseConcordance {
  input {
    File eval
    File truth
    File truth_index
    File af_resource
    File af_resource_index
    String af_tag
    File sample_map
  }

  call RenameSamples {
    input:
      vcf = eval,
      sample_map = sample_map
  }

  call GlimpseConcordance_t {
    input:
      eval = RenameSamples.renamed_vcf,
      eval_index = RenameSamples.renamed_vcf_index,
      truth = truth,
      truth_index = truth_index,
      af_resource = af_resource,
      af_resource_index = af_resource_index,
      af_tag = af_tag
  }

  output {
    File concordance_per_sample = GlimpseConcordance_t.concordance_per_sample
    File concordance_by_frequency = GlimpseConcordance_t.concordance_by_frequency
    File concordance_per_cal_bin = GlimpseConcordance_t.concordance_per_cal_bin
    File r2_by_frequency = GlimpseConcordance_t.r2_by_frequency
    File r2_by_sample = GlimpseConcordance_t.r2_by_sample
  }
}


task GlimpseConcordance_t {
  input {
    File eval
    File eval_index
    File truth
    File truth_index
    File af_resource
    File af_resource_index
    String af_tag

    Int mem_gb = 4
    Int cpu = 4
    Int disk_size_gb = ceil(size(eval, "GiB") + size(truth, "GiB") + size(af_resource, "GB") + 100)
    Int preemptible = 1
    Int max_retries = 3
    String docker = "us.gcr.io/broad-dsde-methods/glimpse:2.0.0"
  }

  command <<<

    NPROC=$(nproc)
    echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."

    echo $(seq -f 'chr%g' -s , 22)'chrX,chrY ~{af_resource} ~{truth} ~{eval}' > concordance_input.txt
    /GLIMPSE/GLIMPSE2_concordance --af-tag ~{af_tag} --use-alt-af --bins 0.0005 0.0008972126509288622 0.0016099810819735928 0.0028889907890056895 0.005184078168625719 0.009302441032590253 0.01669253555791943 0.029953508157290396 0.053749332916843626 0.09644916294395851 0.1730708183296376 0.3105626554239233 0.5572814867048095 1.0 \
    -O concordance_glimpse --threads ${NPROC} --input concordance_input.txt --gt-val

  >>>

  output {
    File concordance_per_sample = 'concordance_glimpse.error.spl.txt.gz'
    File concordance_by_frequency = 'concordance_glimpse.error.grp.txt.gz'
    File concordance_per_cal_bin = 'concordance_glimpse.error.cal.txt.gz'
    File r2_by_frequency = 'concordance_glimpse.rsquare.grp.txt.gz'
    File r2_by_sample = 'concordance_glimpse.rsquare.spl.txt.gz'
  }

  runtime {
    docker: docker
    disks: "local-disk " + disk_size_gb + " HDD"
    memory: mem_gb + " GiB"
    cpu: cpu
    preemptible: preemptible
    maxRetries: max_retries
  }
}

task RenameSamples {
  input {
    File vcf
    File sample_map
  }

  String vcf_base = basename(vcf, ".vcf.gz")
  command <<<

    bcftools reheader ~{vcf} -s ~{sample_map} -o ~{vcf_base}_renamed_samples.vcf.gz
    bcftools index -t ~{vcf_base}_renamed_samples.vcf.gz
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    cpu: 4
    disks: "local-disk 100 HDD"
    preemptible: 0
  }

  output {
    File renamed_vcf = "~{vcf_base}_renamed_samples.vcf.gz"
    File renamed_vcf_index = "~{vcf_base}_renamed_samples.vcf.gz.tbi"
  }
}