version 1.0

workflow GetFingerprintMetrics{
	input {
		File fingerprint_vcf
		File picard_jar
		String output_prefix
		File haplotype_database
	}

	call GetFingerprintMetricsTask{
		input:
			fingerprint_vcf = fingerprint_vcf,
			picard_jar = picard_jar,
			output_prefix = output_prefix,
			haplotype_database = haplotype_database

	}

	output {
		File fingerprint_metrics = GetFingerprintMetricsTask.fingerprint_metrics

	}
}

task GetFingerprintMetricsTask{
	input {
		File fingerprint_vcf
		File picard_jar
		File haplotype_database
		String output_prefix
	}
	command <<<
		java -jar ~{picard_jar} CalculateFingerprintMetrics \
			INPUT=~{fingerprint_vcf} \
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
    	docker: "openjdk:8"
	}
}
