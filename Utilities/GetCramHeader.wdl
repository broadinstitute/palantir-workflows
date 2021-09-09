version 1.0

workflow GetCramHeader{
  input {
    Array[File] input_files
  }

  File input_file = select_first(input_files)

  String input_ending = if sub(input_file, ".*\\.cram$", "is_cram") == "is_cram" then ".cram" else ".bam"


  call getHeader{
    input:
      file = input_file,
      extension = input_ending
  }

  output {
    File header = getHeader.header
  }
}

task getHeader{
  input {
    File file
    String extension
    File basename = basename(file, extension)
  }
  command <<<
    REF_PATH='.' samtools view -H ~{file} > ~{basename}.header.txt
  >>>
  output {
    File header = "~{basename}.header.txt"
  }
  runtime {
    preemptible: 3
    memory: "2 GB"
    cpu: "1"
    disks: "local-disk " + ceil(size(file,'GiB') + 30) +" LOCAL"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"
  }
}
