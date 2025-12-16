version 1.0

workflow Rename_Reheader_vcfs {
    input {
        File vcf_map_file  # File containing <vcf_path>\t<new_name> on each line
        String output_bucket  # GCS bucket path where renamed VCFs will be stored
        Int chunk_size = 100  # Number of VCFs to process per chunk
    }

    call SplitMapFile {
        input:
            vcf_map_file = vcf_map_file,
            chunk_size = chunk_size
    }

    scatter (chunk_file in SplitMapFile.chunk_files) {
        call ProcessVcfs {
            input:
                vcf_map_file = chunk_file,
                output_bucket = output_bucket
        }
    }
}

task SplitMapFile {
    input {
        File vcf_map_file
        Int chunk_size
    }

    command <<<
        split -l ~{chunk_size} -d ~{vcf_map_file} chunk_
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
        File vcf_map_file
        String output_bucket
    }

    command <<<
        set -eou pipefail
        
        while IFS=$'\t' read -r vcf_path new_name; do
            # Copy VCF locally
            /root/google-cloud-sdk/bin/gcloud storage cp "$vcf_path" original.vcf.gz
            
            echo "${new_name}" > new_sample_name.txt
            
            # Reheader the VCF and save with new name
            bcftools reheader -s new_sample_name.txt -o "${new_name}.vcf.gz" original.vcf.gz
            bcftools index -t "${new_name}.vcf.gz"
            
            # Copy to destination bucket
            /root/google-cloud-sdk/bin/gcloud storage cp "${new_name}.vcf.gz" "${new_name}.vcf.gz.tbi" ~{output_bucket}
            
            # Clean up local files
            rm original.vcf.gz new_sample_name.txt "${new_name}.vcf.gz"
        done < ~{vcf_map_file}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 50 HDD"
    }
}
