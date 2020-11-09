version 1.0

workflow FindSamplesToCompare {
	
	input {
		Array[File] input_callset
		
		Array[File] ground_truth_files
		Array[File] ground_truth_intervals

		File haplotype_database 
		File picard_cloud_jar

		String docker 

	}

	Int VCF_disk_size = 50

	File monitoring_script="gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"

	 call CrosscheckFingerprints {
         input:
           input_data = input_callset,
           metrics_basename = "crosscheck",
           ground_truth_files = ground_truth_files,
           haplotype_database = haplotype_database,
           disk_size = VCF_disk_size,
           preemptible_tries = 3,
           docker = docker,
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
    Array[File] input_data
    String metrics_basename
    Array[File] ground_truth_files
    File haplotype_database
    Int disk_size
    Int preemptible_tries
    String docker
    File picard_jar
  }
  parameter_meta {
    input_data: {
      localization_optional: true
    }
   
    ground_truth_files: {
      localization_optional: true
    }
  }  
  
  String tsv_out="~{metrics_basename}.fingerprints_metrics"

  command <<<
    bash ~{monitoring_script} > /cromwell_root/monitoring.log &

    
    java -jar ~{picard_jar} CrosscheckFingerprints \
      INPUT=~{sep=" INPUT=" input_data} \
      SECOND_INPUT=~{sep=" SECOND_INPUT=" ground_truth_files} \
      HAPLOTYPE_MAP=~{haplotype_database} \
      LOD_THRESHOLD=-5 \
      OUTPUT=~{tsv_out} \
      CROSSCHECK_BY=FILE \
      CROSSCHECK_MODE=CHECK_ALL_OTHERS \
      TEST_INPUT_READABILITY=false \
      EXIT_CODE_WHEN_MISMATCH=5
      
   >>>
   output {
   	File crosscheck = tsv_out
   }

   runtime {
    preemptible: preemptible_tries
    memory: "8 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    noAddress: false
    maxRetries: 1
    continueOnReturnCode: [0,5]
  }
}


