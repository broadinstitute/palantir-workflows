version 1.0

struct FastQSet {
    File fastq_r1
    File fastq_r2
    String prefix
}

struct ReferenceFasta {
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File ref_sa
  File ref_amb
  File ref_bwt
  File ref_ann
  File ref_pac
}

struct Sample {
    String name
    String library
    String barcode
}