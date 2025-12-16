version 1.0

workflow Rename_Invitae_files {
    input {
        File id_map  # File containing old_ids on each line
        String input_bucket_path  # GCS bucket path where original files are stored
        String output_bucket_path  # GCS bucket path where renamed files will be stored
        Int chunk_size = 100  # Number of files to process per chunk
    }

    call SplitMapFile {
        input:
            id_map = id_map,
            chunk_size = chunk_size
    }

    scatter (chunk_file in SplitMapFile.chunk_files) {
        call ProcessFiles {
            input:
                id_map = chunk_file,
                input_bucket_path = input_bucket_path,
                output_bucket_path = output_bucket_path
        }
    }
}

task SplitMapFile {
    input {
        File id_map
        Int chunk_size
    }

    command <<<
        split -l ~{chunk_size} -d ~{id_map} chunk_
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

task ProcessFiles {
    input {
        File id_map
        String input_bucket_path  # GCS bucket path where original files are stored
        String output_bucket_path
    }

    command <<<
        set -eou pipefail
        
        while IFS=$'\t' read -r old_id new_id; do
            # rename and copy files
            /root/google-cloud-sdk/bin/gcloud storage cp ~{input_bucket_path}/"$old_id"/aligned_requisitioned_cleaned.bam ~{output_bucket_path}/"$new_id".bam
            /root/google-cloud-sdk/bin/gcloud storage cp ~{input_bucket_path}/"$old_id"/aligned_requisitioned_cleaned.bam.bai ~{output_bucket_path}/"$new_id".bam.bai
            /root/google-cloud-sdk/bin/gcloud storage cp ~{input_bucket_path}/"$old_id"/aligned_requisitioned_cleaned.bam.md5 ~{output_bucket_path}/"$new_id".bam.md5
        done < ~{id_map}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 50 HDD"
    }
}
