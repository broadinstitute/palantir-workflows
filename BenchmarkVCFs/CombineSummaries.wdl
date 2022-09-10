version 1.0

import "BenchmarkVCFs.wdl" as BM
workflow CombineSummaries{
	input{
		Array[File] summaries

	}
	call BM.CombineSummaries as CS{
		input :
			summaries=summaries,
			preemptible=3
		
	}

	output{
		File summary=CS.summaryOut
	}
}