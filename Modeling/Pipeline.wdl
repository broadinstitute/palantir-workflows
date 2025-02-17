version 1.0

workflow VariantAnnotationWorkflow {
  input {
    File input_vcf
    File input_vcf_index
    File reference_fasta
    File truth_vcf
    File truth_vcf_index
    String output_dir
    File no_inconsistent_bed
    File no_inconsistent_bed_index
    Array[File] bed_files
    Array[String] labels
  }

  call PrepareVCF {
    input:
      input_vcf = input_vcf,
      reference_fasta = reference_fasta
  }

  scatter (i in range(length(bed_files))) {
    call AnnotateVCF {
      input:
        input_vcf = PrepareVCF.output_vcf,
        bed_file = bed_files[i],
        label = labels[i]
    }
  }

  scatter (i in range(length(AnnotateVCF.output_vcf))) {
    Boolean do_remove_gt = i > 0
    call RemoveGT {
      input:
        annotated_vcf = AnnotateVCF.output_vcf[i],
        do_remove_gt = do_remove_gt
    }
  }

  Array[File] all_vcfs = RemoveGT.output_vcf
  Array[File] all_vcf_indices = RemoveGT.output_vcf_index

  call MergeVCFs {
    input:
      annotated_vcfs = all_vcfs,
      annotated_vcfs_indices = all_vcf_indices,
      reference_fasta = reference_fasta
  }

  call FinalizeVCF {
    input:
      input_vcf = MergeVCFs.merged_vcf,
      input_vcf_index = MergeVCFs.merged_vcf_index,
      reference_fasta = reference_fasta,
      truth_vcf = truth_vcf,
      truth_vcf_index = truth_vcf_index,
      output_dir = output_dir,
      no_inconsistent_bed = no_inconsistent_bed,
      no_inconsistent_bed_index = no_inconsistent_bed_index,
      labels = labels
  }

  output {
    File annotated_vcf = FinalizeVCF.annotated_vcf
    File annotated_vcf_index = FinalizeVCF.annotated_vcf_index
    File dataset = FinalizeVCF.dataset
    File vcfeval_output_zip = FinalizeVCF.vcfeval_output_zip
  }
}

task PrepareVCF {
  input {
    File input_vcf
    File reference_fasta
  }

  command <<<
    set -xueo pipefail

    # Decompress the input VCF file
    cp ~{input_vcf} annotated_tmp.vcf.gz
    gunzip annotated_tmp.vcf.gz
    if [ $? -ne 0 ]; then
      echo "Failed to decompress the input VCF file."
      exit 1
    fi

    # Verify the decompressed VCF file
    if [ ! -f annotated_tmp.vcf ]; then
      echo "Decompressed VCF file not found: annotated_tmp.vcf"
      exit 1
    fi

    mv annotated_tmp.vcf prepared_vcf.vcf
  >>>

  output {
    File output_vcf = "prepared_vcf.vcf"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
    memory: "16G"
    cpu: "4"
  }
}

task AnnotateVCF {
  input {
    File input_vcf
    File bed_file
    String label
  }

  command <<<
    set -xueo pipefail

    # Check if the input VCF is compressed
    if [[ "~{input_vcf}" == *.gz ]]; then
      cp ~{input_vcf} annotated_tmp.vcf.gz
      gunzip annotated_tmp.vcf.gz
    else
      cp ~{input_vcf} annotated_tmp.vcf
    fi

    # Verify the decompressed VCF file
    if [ ! -f annotated_tmp.vcf ]; then
      echo "Decompressed VCF file not found: annotated_tmp.vcf"
      exit 1
    fi

    # Create the new BED file with an additional column
    awk '{print $0 "\t1"}' ~{bed_file} > ~{bed_file}_with_flag.bed
    if [ $? -ne 0 ]; then
      echo "Failed to create BED file with flag."
      exit 1
    fi

    # Compress and index the new BED file
    bgzip -c ~{bed_file}_with_flag.bed > ~{bed_file}_with_flag.bed.gz
    if [ $? -ne 0 ]; then
      echo "Failed to compress BED file."
      exit 1
    fi
    tabix -p bed ~{bed_file}_with_flag.bed.gz
    if [ $? -ne 0 ]; then
      echo "Failed to index BED file."
      exit 1
    fi

    # Add the header for this label
    HEADER_FILE="header.txt"
    echo "##INFO=<ID=~{label},Number=1,Type=Integer,Description=\"Overlap with ~{label} regions\">" > ${HEADER_FILE}
    if [ $? -ne 0 ]; then
      echo "Failed to add header for ~{label}."
      exit 1
    fi

    # Annotate the VCF file
    bcftools annotate -a ~{bed_file}_with_flag.bed.gz -h ${HEADER_FILE} -c CHROM,FROM,TO,INFO/~{label} annotated_tmp.vcf -O v -o annotated_vcf.vcf
    if [ $? -ne 0 ]; then
      echo "Failed to annotate the VCF file with ~{label}."
      exit 1
    fi

    # Compress and index the annotated VCF file
    bcftools view annotated_vcf.vcf -o annotated_vcf.vcf.gz
    if [ $? -ne 0 ]; then
      echo "Failed to compress the annotated VCF file."
      exit 1
    fi
    bcftools index -t annotated_vcf.vcf.gz
    if [ $? -ne 0 ]; then
      echo "Failed to index the annotated VCF file."
      exit 1
    fi
  >>>

  output {
    File output_vcf = "annotated_vcf.vcf.gz"
    File output_vcf_index = "annotated_vcf.vcf.gz.tbi"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
    memory: "16G"
    cpu: "4"
  }
}

task RemoveGT {
  input {
    File annotated_vcf
    Boolean do_remove_gt
  }

  command <<<
    set -xueo pipefail

    if [ "~{do_remove_gt}" == "true" ]; then
      # Remove GT column from the VCF
      bcftools view -G ~{annotated_vcf} -O z -o annotated_no_gt.vcf.gz
      if [ $? -ne 0 ]; then
        echo "Failed to remove GT column from VCF."
        exit 1
      fi
    else
      cp ~{annotated_vcf} annotated_no_gt.vcf.gz
    fi

    # Index the modified VCF file
    bcftools index -t annotated_no_gt.vcf.gz
    if [ $? -ne 0 ]; then
      echo "Failed to index the modified VCF file."
      exit 1
    fi
  >>>

  output {
    File output_vcf = "annotated_no_gt.vcf.gz"
    File output_vcf_index = "annotated_no_gt.vcf.gz.tbi"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
    memory: "16G"
    cpu: "4"
  }
}

task MergeVCFs {
  input {
    Array[File] annotated_vcfs
    Array[File] annotated_vcfs_indices
    File reference_fasta
  }

  command <<<
    set -xueo pipefail

    # Merge the VCF files with --force-samples option to handle duplicate sample names
    bcftools merge --force-samples -m all -O z -o merged_vcf.vcf.gz ~{sep=' ' annotated_vcfs}
    if [ $? -ne 0 ]; then
      echo "Failed to merge the VCF files."
      exit 1
    fi

    # Index the merged VCF file
    bcftools index -t merged_vcf.vcf.gz
    if [ $? -ne 0 ]; then
      echo "Failed to index the merged VCF file."
      exit 1
    fi
  >>>

  output {
    File merged_vcf = "merged_vcf.vcf.gz"
    File merged_vcf_index = "merged_vcf.vcf.gz.tbi"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
    memory: "16G"
    cpu: "4"
  }
}

task FinalizeVCF {
  input {
    File input_vcf
    File input_vcf_index
    File reference_fasta
    File truth_vcf
    File truth_vcf_index
    String output_dir
    File no_inconsistent_bed
    File no_inconsistent_bed_index
    Array[String] labels
  }

  command <<<
    set -xueo pipefail

    # Verify input VCF and index
    if [ ! -f ~{input_vcf} ]; then
      echo "Input VCF file not found: ~{input_vcf}"
      exit 1
    fi

    if [ ! -f ~{input_vcf_index} ]; then
      echo "Input VCF index file not found: ~{input_vcf_index}"
      exit 1
    fi

    # Verify truth VCF and index
    if [ ! -f ~{truth_vcf} ]; then
      echo "Truth VCF file not found: ~{truth_vcf}"
      exit 1
    fi

    if [ ! -f ~{truth_vcf_index} ]; then
      echo "Truth VCF index file not found: ~{truth_vcf_index}"
      exit 1
    fi

    # Verify no_inconsistent_bed and index
    if [ ! -f ~{no_inconsistent_bed} ]; then
      echo "No inconsistent BED file not found: ~{no_inconsistent_bed}"
      exit 1
    fi

    if [ ! -f ~{no_inconsistent_bed_index} ]; then
      echo "No inconsistent BED index file not found: ~{no_inconsistent_bed_index}"
      exit 1
    fi

    # Create SDF from the FASTA reference genome
    rtg format -o reference.sdf ~{reference_fasta}
    if [ $? -ne 0 ]; then
      echo "Failed to create reference SDF from the FASTA file."
      exit 1
    fi

    # Compress and index the final annotated VCF file using bcftools
    bcftools view ~{input_vcf} -o annotated_final.vcf.gz
    if [ $? -ne 0 ]; then
      echo "Failed to compress the final annotated VCF file."
      exit 1
    fi

    if [ ! -f annotated_final.vcf.gz ]; then
      echo "Compressed VCF file not found: annotated_final.vcf.gz"
      exit 1
    fi

    bcftools index -t -f annotated_final.vcf.gz
    if [ $? -ne 0 ]; then
      echo "Failed to index the final annotated VCF file."
      exit 1
    fi

    # Run vcfeval
    rtg vcfeval -b ~{truth_vcf} -c annotated_final.vcf.gz -t reference.sdf -o ~{output_dir} --output-mode=annotate -e ~{no_inconsistent_bed}
    if [ $? -ne 0 ]; then
      echo "Failed to run rtg vcfeval."
      exit 1
    fi

    if [ ! -f ~{output_dir}/calls.vcf.gz ]; then
      echo "vcfeval output file does not exist."
      exit 1
    fi

    # Check if the output directory is not empty before zipping
    if [ "$(ls -A ~{output_dir})" ]; then
      # Zip the vcfeval output directory
      tar -cvzf vcfeval_output.zip -C ~{output_dir} .
      if [ $? -ne 0 ]; then
        echo "Failed to zip vcfeval output directory."
        exit 1
      fi
    else
      echo "Output directory is empty, not creating archive."
      exit 1
    fi

    # Create dataset with chosen columns from vcfeval output calls.vcf.gz
    query_format="%CHROM\t%POS\t%REF\t%ALT\t%FILTER"
    for label in ~{sep=' ' labels}; do
      query_format+="\t%INFO/${label}"
    done
    query_format+="\t%INFO/CALL\n"
    bcftools query -f "${query_format}" ~{output_dir}/calls.vcf.gz > dataset.txt
    if [ $? -ne 0 ]; then
      echo "Failed to create dataset from vcfeval output."
      exit 1
    fi

    echo "Annotation and evaluation complete. Output written to annotated_final.vcf.gz and vcfeval results in ~{output_dir}"
    echo "Dataset with chosen columns written to dataset.txt"
  >>>

  output {
    File annotated_vcf = "annotated_final.vcf.gz"
    File annotated_vcf_index = "annotated_final.vcf.gz.tbi"
    File dataset = "dataset.txt"
    File vcfeval_output_zip = "vcfeval_output.zip"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
    memory: "32G"
    cpu: "8"
  }
}