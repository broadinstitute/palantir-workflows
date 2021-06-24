version 1.0

workflow IndexCramOrBam{
	input {
		File input_file
	}

	String input_index_ending = if sub(input_file, ".*\\.cram$", "is_cram") == "is_cram" then ".crai" else ".bai"


	call indexFile{
		input:
			file = input_file,
			file_index_name = basename(input_file) + input_index_ending
	}

	output {
		File file_index = indexFile.index
	}
}

task indexFile{
	input {
		File file
		String file_index_name
	}
	command <<<
		REF_PATH='.' samtools index ~{file} ~{file_index_name}
	>>>
	output {
		File index = file_index_name
	}
	runtime {
		preemptible: 3
		memory: "2 GB"
		cpu: "1"
		disks: "local-disk " + ceil(size(file,'GiB') + 30) +" LOCAL"
		docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"
	}
}
