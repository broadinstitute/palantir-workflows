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
        # Keep track of shard number in temp file to handle loss of scope in loop
        echo 0 > number.txt

        # Double loop allows compatibility with one or multiple workflow sample pairs
        gsutil ls gs://~{ws_name}/~{sub_id}/FindSamplesAndBenchmark/ | while read -r workflow ; do
            gsutil ls "${workflow}call-BenchmarkVCF" | while read -r shard ; do
                benchmark_url=$(gsutil ls "${shard}Benchmark")
                summary_url=$(gsutil ls "${benchmark_url}call-CombineSummaries/summary.csv")

                # Check if summary exists, i.e. shard was successful
                return_code=$(echo $?)
                if [ $return_code -eq 0 ]; then
                    i=$(cat number.txt)
                    gsutil cp $summary_url "shard-${i}.csv"
                    ((++i))
                    echo $i > number.txt
                fi
            done
        done

        # Create summary file and combine all csv's into one
        > CombinedBenchmarkSummaries.csv
        head -1 shard-0.csv >> CombinedBenchmarkSummaries.csv
        awk FNR-1 shard-*.csv >> CombinedBenchmarkSummaries.csv
        rm shard-*.csv number.txt
    >>>

    runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:alpine"
    }

    output {
        File combined_summaries = "CombinedBenchmarkSummaries.csv"
    }
}