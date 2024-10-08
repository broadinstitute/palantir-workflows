import sys
import os
import shutil
import json
import tarfile
import firecloud.api as fapi
import pandas as pd


# Parse cmd args
NAMESPACE = sys.argv[1]
WORKSPACE = sys.argv[2]
SUBMISSION_ID = sys.argv[3]


# Make output dir location
print("Creating directory for WDL outputs...")
try:
    os.mkdir('./wdl_outputs')
except FileExistsError:
    print("ERROR: The directory wdl_outputs already exists. Please delete or archive it before rerunning.")
    sys.exit(1)

# Query for workflow output data
print("Querying Terra API for workflow outputs...")
response = fapi.get_submission(
    namespace=NAMESPACE,
    workspace=WORKSPACE,
    submission_id=SUBMISSION_ID,
)
    
file_df_dict = {}
file_basenames = {}

os.mkdir('./wdl_outputs/workflows/')

# Parse and download output files
print("Parsing outputs and organizing data...")
output_json = json.loads(response.content)
for i, workflow in enumerate(output_json['workflows']):
    num_workflows = len(output_json['workflows'])
    print(f"Copying file {i+1} of {num_workflows}")
    if 'workflowId' in workflow:
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
            print(f'Extracting files from archive [{i+1}/{num_workflows}]')
            with tarfile.open(f'{wf_path}/benchmark_outputs.tar.gz') as tar:
                tar.extractall(path=wf_path)
            os.remove(f'{wf_path}/benchmark_outputs.tar.gz')
        else:
            print(f"WARNING: Workflow {wf_json['workflowId']} seems to have failed... Skipping data collection for this run.")
    else:
        print("WARNING: Workflow seems to have failed to launch... Skipping data collection for this run.")

print('Consolidating files across workflow runs...')
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
    for i, wf in enumerate(workflows):
        print(f'Loading workflow {i} of {len(workflows)}...')
        for d in dir_names:
            if file_name in os.listdir(f"./wdl_outputs/workflows/{wf}/benchmark_outputs/{d}/"):
                df = pd.read_csv(f"./wdl_outputs/workflows/{wf}/benchmark_outputs/{d}/{file_name}", sep='\t', low_memory=False)
                df['Terra_workflow_id'] = wf
                full_df = pd.concat([full_df, df])
    full_df.to_csv(f"./wdl_outputs/{file_name}", sep='\t', index=False)

shutil.rmtree('./wdl_outputs/workflows/')

# Create README file
print("Creating README file...")

USER_COMMENT = output_json['userComment']

with open('./wdl_outputs/README.txt', 'w') as file:
    lines = []
    lines += ['Files in this directory were created using the gather_terra_data.py script provided with the SVisualizer script.\n']
    lines += ['Files copied on: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}\n']
    lines += ['Taken from:\n']
    lines += [f'\tNamespace: {NAMESPACE}\n']
    lines += [f'\tWorkspace: {WORKSPACE}\n']
    lines += [f'\tSubmission ID: {SUBMISSION_ID}\n']
    lines += [f'User Comment: {USER_COMMENT}']
    file.writelines(lines)

print("Finished!")