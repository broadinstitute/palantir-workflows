version 1.0

import "../../BenchmarkVCFs/BenchmarkVCFs.wdl"

workflow testBenchmarkVCFs {
	input {
		Float expected_snpPrecision
		Float expected_indelPrecision
		Float expected_snpRecall
		Float expected_indelRecall

		File evalVcf
		File truthVcf
	}

	call BgzipAndIndex as BgzipAndIndexEval {
		input:
			vcf = evalVcf
	}

	call BgzipAndIndex as BgzipAndIndexTruth {
		input:
			vcf = truthVcf
	}

	call BenchmarkVCFs.Benchmark {
		input:
			evalVcf = BgzipAndIndexEval.bgzipped_vcf,
			evalVcfIndex = BgzipAndIndexEval.vcf_index,
			truthVcf = BgzipAndIndexTruth.bgzipped_vcf,
			truthVcfIndex = BgzipAndIndexTruth.vcf_index
	}

	call AssertPassed {
		input:
			expected_snpPrecision = expected_snpPrecision,
			expected_indelPrecision = expected_indelPrecision,
			expected_snpRecall = expected_snpRecall,
			expected_indelRecall = expected_indelRecall,
			observed_snpPrecision = Benchmark.snpPrecision,
			observed_indelPrecision = Benchmark.indelPrecision,
			observed_snpRecall = Benchmark.snpRecall,
			observed_indelRecall = Benchmark.indelRecall,

	}
}

task BgzipAndIndex {
	input {
		File vcf
	}

	command <<<
		set -xeuo pipefail

		ln -s ~{vcf} .
		bgzip ~{basename(vcf)}
		tabix ~{basename(vcf) + ".gz"}
	>>>

	runtime {
		docker: "biocontainers/tabix@sha256:7e093436d00c01cf6ad7b285680bf1657f9fcb692cc083c972e5df5a7e951f49"
	}

	output {
		File bgzipped_vcf = "~{basename(vcf) + '.gz'}"
		File vcf_index = "~{basename(vcf) + '.gz.tbi'}"
	}
}

task AssertPassed {
	input {
		Float expected_snpPrecision
		Float expected_indelPrecision
		Float expected_snpRecall
		Float expected_indelRecall
		
		Float observed_snpPrecision
		Float observed_indelPrecision
		Float observed_snpRecall
		Float observed_indelRecall
	}

	command <<<
		set -euo pipefail
		
		assert_eq() {
			local variable="$1"
			local expected="$2"
			local observed="$3"
			
			if [[ $expected != $observed ]]; then
				>&2 echo $variable expected to be $expected, observed as $observed
				exit 1;
			fi
		}
		
		assert_eq snpPrecision ~{expected_snpPrecision} ~{observed_snpPrecision}
		assert_eq indelPrecision ~{expected_indelPrecision} ~{observed_indelPrecision}
		assert_eq snpRecall ~{expected_snpRecall} ~{observed_snpRecall}
		assert_eq indelRecall ~{expected_indelRecall} ~{observed_indelRecall}
	>>>

    runtime {
        docker: "ubuntu@sha256:134c7fe821b9d359490cd009ce7ca322453f4f2d018623f849e580a89a685e5d"
    }
}