version 1.0

import "../../BenchmarkVCFs/BenchmarkVCFs.wdl"

workflow testBenchmarkVCFs {
	input {
		Float expected_snpPrecision
		Float expected_indelPrecision
		Float expected_snpRecall
		Float expected_indelRecall
		Float expected_snpF1Score
		Float expected_indelF1Score
	}

	call BenchmarkVCFs.Benchmark

	call AssertPassed {
		input:
			expected_snpPrecision = expected_snpPrecision,
			expected_indelPrecision = expected_indelPrecision,
			expected_snpRecall = expected_snpRecall,
			expected_indelRecall = expected_indelRecall,
			expected_snpF1Score = expected_snpF1Score,
			expected_indelF1Score = expected_indelF1Score,
			observed_snpPrecision = Benchmark.snpPrecision,
			observed_indelPrecision = Benchmark.indelPrecision,
			observed_snpRecall = Benchmark.snpRecall,
			observed_indelRecall = Benchmark.indelRecall,
			observed_snpF1Score = Benchmark.snpF1Score,
			observed_indelF1Score = Benchmark.indelF1Score
	}

}
task AssertPassed {
	input {
		Float expected_snpPrecision
		Float expected_indelPrecision
		Float expected_snpRecall
		Float expected_indelRecall
		Float expected_snpF1Score
		Float expected_indelF1Score
		
		Float observed_snpPrecision
		Float observed_indelPrecision
		Float observed_snpRecall
		Float observed_indelRecall
		Float observed_snpF1Score
		Float observed_indelF1Score
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
		assert_eq snpF1Score ~{expected_snpF1Score} ~{observed_snpF1Score}
		assert_eq indelF1Score ~{expected_indelF1Score} ~{observed_indelF1Score}
	>>>

    runtime {
        docker: "ubuntu@sha256:134c7fe821b9d359490cd009ce7ca322453f4f2d018623f849e580a89a685e5d"
    }
}