version 1.0

workflow GlimpseConcordance {
  input {
    File eval
    File eval_index
    File truth
    File truth_index
    File af_resource
    File af_resource_index
    String af_tag
  }

  call GlimpseConcordance_t {
    input:
      eval = eval,
      eval_index = eval_index,
      truth = truth,
      truth_index = truth_index,
      af_resource = af_resource,
      af_resource_index = af_resource_index,
      af_tag = af_tag
  }

  output {
    File concordance = GlimpseConcordance_t.concordance
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
    Int disk_size_gb = ceil(size(eval, "GiB") + size(truth, "GiB") + size(af_resource) + 100)
    Int preemptible = 1
    Int max_retries = 3
    String docker = "us.gcr.io/broad-dsde-methods/glimpse:2.0.0"
  }

  command <<<

    #NPROC=$(nproc)
    #echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."

    echo $(seq -f 'chr%g' -s , 22)'chrX,chrY ~{af_resource} ~{truth} ~{eval}' > concordance_input.txt
    /GLIMPSE/GLIMPSE2_concordance --af-tag ~{af_tag} --use-alt-af --bins 0.0005 0.0008972126509288622 0.0016099810819735928 0.0028889907890056895 0.005184078168625719 0.009302441032590253 0.01669253555791943 0.029953508157290396 0.053749332916843626 0.09644916294395851 0.1730708183296376 0.3105626554239233 0.5572814867048095 1.0 \
    -O concordance_glimpse -T ${NPROC}

  >>>

  output {
    File concordance = 'concordance_glimpse}'
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