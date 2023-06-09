import sys, os, shutil
import json
import tarfile
import firecloud.api as fapi
import pandas as pd


# Parse cmd args
NAMESPACE = sys.argv[1]
WORKSPACE = sys.argv[2]
SUBMISSION_ID = sys.argv[3]

# Query for workflow output data
response = fapi.get_submission(
    namespace=NAMESPACE,
    workspace=WORKSPACE,
    submission_id=SUBMISSION_ID,
)

# Make output dir location
try:
    os.mkdir('./wdl_outputs')
except FileExistsError:
    print("ERROR: The directory wdl_outputs already exists. Please delete or archive it before rerunning.")
    sys.exit(1)
    
file_df_dict = {}
file_basenames = {}

os.mkdir('./wdl_outputs/workflows/')

# Parse and download output files
output_json = json.loads(response.content)
for workflow in output_json['workflows']:
    wf_response = fapi.get_workflow_outputs(
        namespace=NAMESPACE,
        workspace=WORKSPACE,
        submission_id=SUBMISSION_ID,
        workflow_id=workflow['workflowId']
    )
    
    wf_json = json.loads(wf_response.content)
    
    # Get dict of output names -> output files from workflow
    if 'BenchmarkSVs' in wf_json['tasks']:
        wf_path = f'./wdl_outputs/workflows/{workflow["workflowId"]}/'
        gs_path = wf_json['tasks']['BenchmarkSVs']['outputs']['BenchmarkSVs.combined_files']
        
        os.mkdir(wf_path)
        os.system(f'gsutil cp {gs_path} {wf_path}')
        with tarfile.open(f'{wf_path}/benchmark_outputs.tar.gz') as tar:
            tar.extractall(path=wf_path)
        os.remove(f'{wf_path}/benchmark_outputs.tar.gz')
        # output_file_dict = wf_json['tasks']['BenchmarkSVs']['outputs']
    
        # Save each output file locally
        # for k in output_file_dict:
        #     basename = output_file_dict[k].split('/')[-1]
        #     file_basenames[k] = basename
        #     os.system(f'gsutil cp {output_file_dict[k]} wdl_outputs/{basename}')
        #     df = pd.read_csv(f'wdl_outputs/{basename}', sep='\t')
        #     file_df_dict[k] = pd.concat([file_df_dict[k], df]) if k in file_df_dict.keys() else df
    else:
        print(f"WARNING: Workflow {wf_json['workflowId']} seems to have failed... Skipping data collection for this run.")

# for k in file_df_dict:
#     file_df_dict[k].to_csv(f'wdl_outputs/{file_basenames[k]}', sep='\t', index=False)

print('Combining files across workflow runs...')

# Get list of file names across the different stat categories
workflows = os.listdir('./wdl_outputs/workflows/')
dir_names = os.listdir(f"./wdl_outputs/workflows/{workflows[0]}/benchmark_outputs/")
file_names = []
for d in dir_names:
    first_workflow = os.listdir('./wdl_outputs/workflows/')[0]
    file_names += os.listdir(f"./wdl_outputs/workflows/{first_workflow}/benchmark_outputs/{d}")

# Concatenate all file outputs across all workflow runs
for file_name in file_names:
    full_df = pd.DataFrame()
    for wf in workflows:
        for d in dir_names:
            if file_name in os.listdir(f"./wdl_outputs/workflows/{wf}/benchmark_outputs/{d}/"):
                df = pd.read_csv(f"./wdl_outputs/workflows/{wf}/benchmark_outputs/{d}/{file_name}", sep='\t', low_memory=False)
                full_df = pd.concat([full_df, df])
    full_df.to_csv(f"./wdl_outputs/{file_name}", sep='\t', index=False)

shutil.rmtree('./wdl_outputs/workflows/')