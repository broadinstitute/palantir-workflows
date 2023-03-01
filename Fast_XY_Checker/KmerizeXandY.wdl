version 1.0

workflow KmerizeXandY {
  input {
    File ref_fasta
    File ref_index

    String x_chr_name
    String y_chr_name

    Int? kmer_size
    String? kmer_mask
    Float? bfpp
  }

  call ExtractContig as ExtractContigX {
    input:
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      contig = x_chr_name
  }

  call ExtractContig as ExtractContigY {
    input:
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      contig = y_chr_name
  }

  call KmerizeContig as KmerizeX {
    input:
      contig_fasta = ExtractContigX.contig_fasta,
      contig_index = ExtractContigX.contig_index,
      contig_dict = ExtractContigX.contig_dict,
      kmer_size = kmer_size,
      kmer_mask = kmer_mask,
      bfpp = bfpp
  }

  call KmerizeContig as KmerizeY {
    input:
      contig_fasta = ExtractContigY.contig_fasta,
      contig_index = ExtractContigY.contig_index,
      contig_dict = ExtractContigY.contig_dict,
      kmer_size = kmer_size,
      kmer_mask = kmer_mask,
      bfpp = bfpp
  }

  output {
    File x_kmer_file = KmerizeX.kmerized_contig
    File y_kmers_file = KmerizeY.kmerized_contig
  }

}

task ExtractContig {
  input {
    File ref_fasta
    File ref_index

    String contig
  }

  Int disk_size_gb = 2*ceil(size(ref_fasta, "GB")) + 20
  String basename = basename(ref_fasta,".fasta") + "_" + contig

  command <<<
    samtools faidx ~{ref_fasta} ~{contig} -o ~{basename}.fasta
    samtools faidx ~{basename}.fasta
    samtools dict ~{basename}.fasta -o ~{basename}.dict
  >>>

  output {
    File contig_fasta = "~{basename}.fasta"
    File contig_index = "~{basename}.fasta.fai"
    File contig_dict = "~{basename}.dict"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    disks: "local-disk ~{disk_size_gb} HDD"
  }
}

task KmerizeContig {
  input {
    File contig_fasta
    File contig_index
    File contig_dict

    Int kmer_size = 31
    String kmer_mask = "16"
    Float bfpp = 0.001
  }

  String basename = basename(contig_fasta, ".fasta")

  command <<<
    gatk --java-options "-Xmx12g" PathSeqBuildKmers --reference ~{contig_fasta} \
    --output ~{basename}_kmer_size_~{kmer_size}_kmer_mask_~{kmer_mask}_bfpp_~{bfpp}.bfi \
    --kmer-size ~{kmer_size} ~{"--kmer-mask " + kmer_mask} --bloom-false-positive-probability ~{bfpp}

  >>>

  runtime {
    docker : "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    disks: "local-disk 100 HDD"
    memory: "16 GB"
  }

  output {
    File kmerized_contig = "~{basename}_kmer_size_~{kmer_size}_kmer_mask_~{kmer_mask}_bfpp_~{bfpp}.bfi"
  }
}