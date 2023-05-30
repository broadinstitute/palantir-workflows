import sys, os
import json
import firecloud.api as fapi


# Parse cmd args
NAMESPACE = sys.argv[1]
WORKSPACE = sys.argv[2]
SUBMISSION_ID = sys.argv[3]
WORKFLOW_ID = sys.argv[4]

# Query for workflow output data
response = fapi.get_workflow_outputs(
    namespace=NAMESPACE,
    workspace=WORKSPACE,
    submission_id=SUBMISSION_ID,
    workflow_id=WORKFLOW_ID
)

# Make output dir location
os.mkdir('./wdl_outputs')

# Parse and download output files
output_json = json.loads(response.content)
output_file_dict = output_json['tasks']['BenchmarkSVs']['outputs']
for k in output_file_dict:
    os.system(f'gsutil cp {output_file_dict[k]} wdl_outputs')