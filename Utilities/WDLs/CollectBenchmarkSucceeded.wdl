version 1.0

workflow CollectBenchmarkSucceeded {
    input {
        String workspace_bucket_name
        String submission_id
    }

    call CombineSucceededVCFs {
        input:
            ws_name = workspace_bucket_name,
            sub_id = submission_id
    }

    output {
        File combined_summaries = CombineSucceededVCFs.combined_summaries
    }
}


task CombineSucceededVCFs {
    input {
        String ws_name
        String sub_id
    }

    command <<<
        i=0

        # Loop over all shards for all benchmarking tasks
        for WORKFLOW in $(gsutil ls gs://~{ws_name}/submissions/~{sub_id}/FindSamplesAndBenchmark/)
        do
            for SHARD in $(gsutil ls "${WORKFLOW}call-BenchmarkVCF")
            do
               benchmark_url=$(gsutil ls "${SHARD}Benchmark")
               summary_url=$(gsutil ls "${benchmark_url}call-CombineSummaries/summary.csv")

                # Check if summary exists, i.e. shard was successful
                return_code=$(echo $?)

                 # If successful, copy benchmark output to ith shard
                 if [ $return_code -eq 0 ]; then
                    gsutil cp $summary_url "shard-${i}.csv"
                    ((++i))
                 fi
            done
        done

        # Create summary file
        > CombinedBenchmarkSummaries.csv

        # Copy header from first shard
        head -1 shard-0.csv >> CombinedBenchmarkSummaries.csv

        # Combine all summaries into one and clean
        awk FNR-1 shard-*.csv >> CombinedBenchmarkSummaries.csv
        rm shard-*.csv
    >>>

    runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:alpine"
    }

    output {
        File combined_summaries = "CombinedBenchmarkSummaries.csv"
    }
}
