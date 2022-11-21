version 1.0

import "../../BenchmarkVCFs/CompareBenchmarks.wdl" as CompareBenchmarks

workflow testCompareBenchmarks {
    input {
        Array[String] sample_ids
        Array[String] configurations
        Array[File] benchmark_summaries
        Array[String] stratifiers

        Array[String] order_of_samples
        Array[String] order_of_configurations
        Array[Int] deltas

        File expected_csv
        Array[File] expected_gc_plots
    }

    call CompareBenchmarks.CompareBenchmarks as CompareBenchmarks {
        input:
            sample_ids = sample_ids,
            configurations = configurations,
            benchmark_summaries = benchmark_summaries,
            stratifiers = stratifiers,
            include_counts = true,
            order_of_samples = order_of_samples,
            order_of_configurations = order_of_configurations,
            deltas = deltas,
    }

    call AssertPassed {
        input:
            expected_csv = expected_csv,
            expected_gc_plots = expected_gc_plots,
            observed_csv = CompareBenchmarks.comparison_csv,
            observed_gc_plots = select_first([CompareBenchmarks.gc_plots, []])

    }
}

task AssertPassed {
    input {
        File expected_csv
        File observed_csv
        Array[File] expected_gc_plots
        Array[File] observed_gc_plots
    }

    command <<<
        set -euxo pipefail
        
        if ! cmp ~{expected_csv} ~{observed_csv} ; then
            echo "=================================== FAILURE ==================================="
            echo "Expected (~{expected_csv}) and observed (~{observed_csv}) CSV tables differ."
            echo "=================================== FAILURE ==================================="
            exit 1
        fi
        
        declare -a expected_gc=("~{sep='" "' expected_gc_plots}")
        declare -a observed_gc=("~{sep='" "' observed_gc_plots}")

        for i in "${!expected_gc[@]}" ; do
            if ! cmp "${expected_gc[$i]}" "${observed_gc[$i]}" ; then
                echo "=================================== FAILURE ==================================="
                echo "Expected (${expected_gc[$i]}) and observed (${observed_gc[$i]}) GC plots differ."
                echo "=================================== FAILURE ==================================="
                exit 1
            fi
        done
    >>>

    runtime {
        docker: "ubuntu@sha256:134c7fe821b9d359490cd009ce7ca322453f4f2d018623f849e580a89a685e5d"
    }
}