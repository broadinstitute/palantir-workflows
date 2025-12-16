version 1.0

workflow Rename_File {
    input {
        File file_rename_map  # File containing <vcf_path>\t<new_name> on each line
        String output_bucket  # GCS bucket path where renamed VCFs will be stored
        Int chunk_size = 100  # Number of VCFs to process per chunk
    }

    call SplitMapFile {
        input:
            file_rename_map = file_rename_map,
            chunk_size = chunk_size
    }

    scatter (chunk_file in SplitMapFile.chunk_files) {
        call ProcessVcfs {
            input:
                file_rename_map = chunk_file,
                output_bucket = output_bucket
        }
    }
}

task SplitMapFile {
    input {
        File file_rename_map
        Int chunk_size
    }

    command <<<
        split -l ~{chunk_size} -d ~{file_rename_map} chunk_
        for f in chunk_*; do
            mv $f $f.txt
        done
    >>>

    output {
        Array[File] chunk_files = glob("chunk_*.txt")
    }

    runtime {
        docker: "ubuntu:latest"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 HDD"
    }
}

task ProcessVcfs {
    input {
        File file_rename_map
        String output_bucket
    }

    command <<<
        set -eou pipefail
        
        while IFS=$'\t' read -r old_path new_name; do
            /root/google-cloud-sdk/bin/gcloud storage cp "$old_path" ~{output_bucket}/$new_name
        done < ~{file_rename_map}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 50 HDD"
    }
}
