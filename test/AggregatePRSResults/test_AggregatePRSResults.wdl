version 1.0

import "../../ImputationPipeline/AggregatePRSResults.wdl"

workflow test_AggregatePRSResults {
    input {
        Array[File] results
        Array[File] target_pc_projections
        Array[File] missing_sites_shifts
        File high_risk_thresholds
        File population_pc_projections
        String population_name = "Reference Population"
        File expected_control_results
        String lab_batch
        Int group_n

        File expected_batch_all_results
        File expected_batch_control_results
        File expected_batch_summarised_results
        File expected_batch_pcs
    }
    
    call AggregatePRSResults.AggregatePRSResults {
        input:
            results = results,
            target_pc_projections = target_pc_projections,
            missing_sites_shifts = missing_sites_shifts,
            high_risk_thresholds = high_risk_thresholds,
            population_pc_projections = population_pc_projections,
            population_name = population_name,
            expected_control_results = expected_control_results,
            lab_batch = lab_batch,
            group_n = group_n
    }
    
    call CompareTextFiles {
        input:
            test_text_files = [AggregatePRSResults.batch_all_results,
                              AggregatePRSResults.batch_control_results,
                              AggregatePRSResults.batch_summarised_results,
                              AggregatePRSResults.batch_pcs],
            truth_text_files = [expected_batch_all_results,
                               expected_batch_control_results,
                               expected_batch_summarised_results,
                               expected_batch_pcs]
    }
}

    
#copied from WARP
task CompareTextFiles {
    input {
    Array[File] test_text_files
    Array[File] truth_text_files
    }

    command {
        exit_code=0

        test_files_length=~{length(test_text_files)}
        truth_files_length=~{length(truth_text_files)}
        if [ $test_files_length -ne $truth_files_length ]; then
        exit_code=1
        echo "Error: Different number of input files ($test_files_length vs. $truth_files_length).  This is really not OK"
        fi

        while read -r a && read -r b <&3;
        do
        echo "Comparing File $a with $b"
        diff $a $b > diffs.txt
        if [ $? -ne 0 ];
        then
        exit_code=1
        echo "Error: Files $a and $b differ" >&2
        cat diffs.txt >&2
        fi
        # catting the diffs.txt on STDOUT as that's what's expected.
        cat diffs.txt
        done < ~{write_lines(test_text_files)} 3<~{write_lines(truth_text_files)}

        exit $exit_code
    }

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
        disks: "local-disk 10 HDD"
        memory: "2 GiB"
        preemptible: 3
    }
}