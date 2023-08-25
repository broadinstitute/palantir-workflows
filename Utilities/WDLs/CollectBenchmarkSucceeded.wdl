version 1.0

workflow CollectBenchmarkSucceeded {
    input {
        String namespace
        String workspace_name
        String submission_id
    }

    call CombineSucceeded {
        input:
            namespace = namespace,
            workspace_name = workspace_name,
            submission_id = submission_id
    }

    output {
        File combined_summaries = CombineSucceeded.combined_summaries
    }
}


task CombineSucceeded {
    input {
        String namespace
        String workspace_name
        String submission_id
    }

    command <<<
        python <<CODE
        import pandas as pd
        import firecloud.api as fapi


        namespace = "~{namespace}"
        workspace_name = "~{workspace_name}"
        submission_id = "~{submission_id}"

        # Get all workflows associated with given submission id
        benchmark_submission = fapi.get_submission(namespace, workspace_name, submission_id).json()
        benchmark_df = pd.DataFrame()

        # Loop over all workflows to get outputs from those that succeeded
        for wf in benchmark_submission['workflows']:
            if wf['status'] == 'Succeeded':
                wf_meta = fapi.get_workflow_metadata(namespace, workspace_name, submission_id, wf['workflowId']).json()
                summary = wf_meta['outputs']['FindSamplesAndBenchmark.benchmark_vcf_summary']
                new_df = pd.read_csv(summary)
                benchmark_df = pd.concat([benchmark_df, new_df])    # Combine into one df

        benchmark_df.to_csv('CombinedBenchmarkSummaries.csv', index=False)

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim"
    }

    output {
        File combined_summaries = "CombinedBenchmarkSummaries.csv"
    }
}
