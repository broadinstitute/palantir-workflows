version 1.0

workflow FindSamplesToCompare{
	
	input {
		File input_reads
		File input_reads_index

		Array[File] ground_truth_files
		Array[File] ground_truth_indexes
		Array[File] ground_truth_intervals

		File haplotype_database 
		File picard_cloud_jar


	}

	Int VCF_disk_size = 50

	File monitoring_script="gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"

	 call CrosscheckFingerprints {
         input:
           input_bam = deduplicated_bam,
           input_bam_index = deduplicated_bam_index,
           ground_truth_files = ground_truth_files,
           haplotype_database = haplotype_database,
           disk_size = VCF_disk_size,
           preemptible_tries = 3,
           docker = jukebox_vc_docker,
           monitoring_script = monitoring_script,
           picard_jar = picard_cloud_jar,
      }


   output {
   		File crosscheck = CrosscheckFingerprints.crosscheck
   }

}



task CrosscheckFingerprints {
  input {
    File monitoring_script
    File input_reads
    File ground_truth_files
    File haplotype_database
    Int disk_size
    Int preemptible_tries
    String docker
    File picard_jar
  }
  parameter_meta {
    input_reads: {
      localization_optional: true
    }
   
    ground_truth_files: {
      localization_optional: true
    }
  }  
  
  String tsv_out="~{final_vcf_base_name}.fingerprints_metrics"

  command <<<
    bash ~{monitoring_script} > /cromwell_root/monitoring.log &

    source ~/.bashrc
    conda activate genomics.py3

    java -jar ~{picard_jar} CrosscheckFingerprints \
      INPUT=~{input_reads} \
      SECOND_INPUT=~{sep=" SECOND_INPUT=" ground_truth_files}
      HAPLOTYPE_MAP=~{haplotype_database} \
      LOD_THRESHOLD=-5 \
      OUTPUT=~{tsv_out} \
      CROSSCHECK_BY=SAMPLE \
      CROSSCHECK_MODE=CHECK_ALL_OTHERS \
      TEST_INPUT_READABILITY=false
      
   >>>
   output {
   	File crosscheck = ~{tsv_out}
   }

   runtime {
    preemptible: preemptible_tries
    memory: "8 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: false
    maxRetries: 1
  }
}


