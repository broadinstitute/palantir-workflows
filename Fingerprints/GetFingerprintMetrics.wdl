version 1.0

workflow GetFingerprintMetrics {
	input {
		File fingerprint_vcf_or_bam
		File? fingerprint_vcf_or_bam_index
		String output_prefix
		File haplotype_database
	}

	call GetFingerprintMetricsTask {
		input:
			fingerprint_vcf_or_bam = fingerprint_vcf_or_bam,
			fingerprint_vcf_or_bam_index = fingerprint_vcf_or_bam_index,
			output_prefix = output_prefix,
			haplotype_database = haplotype_database
	}

	output {
		File fingerprint_metrics = GetFingerprintMetricsTask.fingerprint_metrics
	}
}

task GetFingerprintMetricsTask{
	input {
		File fingerprint_vcf_or_bam
		File? fingerprint_vcf_or_bam_index
		File haplotype_database
		String output_prefix
	}
	parameter_meta {
	    fingerprint_vcf_or_bam: {
	      description: "A file containing a fingerprint (can be sam/bam/cram or vcf). Should be indexed.",
	      localization_optional: true
	    }
	    fingerprint_vcf_or_bam_index: {
	      description: "The index for the file containing a fingerprint.",
	      localization_optional: true
	    }
  	}
	command <<<
		java -jar /usr/picard/picard.jar CalculateFingerprintMetrics \
			INPUT=~{fingerprint_vcf_or_bam} \
			OUTPUT=~{output_prefix}.fingerprint_metrics \
			HAPLOTYPE_MAP=~{haplotype_database}
	>>>
	output {
		File fingerprint_metrics="~{output_prefix}.fingerprint_metrics"
	}
	runtime {
		preemptible: 3
    	memory: "2 GB"
    	cpu: "1"
    	disks: "local-disk 100 LOCAL"
    	docker: "docker pull broadinstitute/picard:2.25.6"
	}
}
