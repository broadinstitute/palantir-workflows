version 1.0

workflow ExtractSingleSampleVCFFromCallset{
	input {
		File vcf
		String basename
		String sampleName

	}

	call ExtractSingleSampleVCFFromCallsetTask{
		input:
			vcf_callset=vcf,
			basename=basename,
			sampleName=sampleName
	}

	output {
		File output_vcf=ExtractSingleSampleVCFFromCallsetTask.vcf
		File output_vcf_index = ExtractSingleSampleVCFFromCallsetTask.vcf_index
	}
}


task ExtractSingleSampleVCFFromCallsetTask{
	input {
		File vcf_callset
		String basename
		String sampleName
		Int disk_size=ceil(size(vcf_callset,"GB")+20)
	}

	command <<<
		gatk SelectVariants -V ~{vcf_callset} --sample-name ~{sampleName} -O ~{basename}.vcf.gz
	>>>

	output {
		File vcf="~{basename}.vcf.gz"
		File vcf_index="~{basename}.vcf.gz.tbi"
	}
	runtime {
		docker:"us.gcr.io/broad-dsde-methods/imputation_bcftools_vcftools_docker:v1.0.0"
		memory:"2 GB"
		disk: "local-disk " + disk_size + " HDD"
		noAddress: false
        maxRetries: 1
	}

}