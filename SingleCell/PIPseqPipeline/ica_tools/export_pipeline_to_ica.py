import requests
import subprocess
import json
import datetime
import os
import time

from ica_common import API_URL, PROJECT_NAMES_AND_IDS, load_api_key, prompt_choice

api_url = API_URL
script_dir = os.path.dirname(os.path.abspath(__file__))

ica_api_key = load_api_key()

current_git_commit_id = subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip()
current_git_commit_id_short = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode().strip()
current_git_commit_message = subprocess.check_output(['git', 'log', '-1', '--pretty=%B']).decode().strip().split('\n')[0].strip()

entrypoints = {
    'Full pipeline (main.nf, --fastq_list)': {
        'name_suffix': '',
        'main_file_path': 'SingleCell/PIPseqPipeline/main.nf',
        'input_form_dir': 'main',
    },
    'Simple single-subsample entrypoint (main_simple.nf, flat FASTQ params)': {
        'name_suffix': '_Simple',
        'main_file_path': 'SingleCell/PIPseqPipeline/main_simple.nf',
        'input_form_dir': 'main_simple',
    },
}

entrypoint_names = list(entrypoints.keys())
entrypoint_choice = prompt_choice('Which entrypoint do you want to export?', entrypoint_names)
if entrypoint_choice is None:
    exit(0)
entrypoint = entrypoints[entrypoint_names[entrypoint_choice]]

input_form_path = os.path.join(script_dir, 'inputforms', entrypoint['input_form_dir'], 'inputForm.json')
if not os.path.isfile(input_form_path):
    print(f'ERROR: input form file not found at {input_form_path}')
    exit(1)

pipeline_name = f'PIPseq_BCL{entrypoint["name_suffix"]}_{current_git_commit_id_short}'

repository_url = 'https://github.com/broadinstitute/palantir-workflows'
main_file_path = entrypoint['main_file_path']
# Both entrypoints share the same process/resource config regardless of which is exported.
nextflow_config_path = 'SingleCell/PIPseqPipeline/nextflow.config'

git_credential_uuid = '5a2282d8-61a7-4222-8969-bfefbbe4f949'

print(f'Exporting pipeline with the following data:')
print(f'  Pipeline name: {pipeline_name}')
print(f'  Pipeline version: {current_git_commit_id_short}')
print(f'  Commit: {current_git_commit_id}')
print(f'  Repository URL: {repository_url}')
print(f'  Main file path: {main_file_path}')
print(f'  Nextflow config path: {nextflow_config_path}')
print(f'  Git credential UUID: {git_credential_uuid}')
print('')
project_names = list(PROJECT_NAMES_AND_IDS.keys())
project_choice = prompt_choice('Which project do you want to export to?', project_names)
if project_choice is None:
    exit(0)
project_id = PROJECT_NAMES_AND_IDS[project_names[project_choice]]

headers = {
    'X-API-Key': ica_api_key,
    'Accept': 'application/vnd.illumina.v4+json',
}

# For multipart/form-data, use the files parameter with (None, value) tuples
# This forces requests to send as multipart/form-data instead of application/x-www-form-urlencoded
files = {
    'language': (None, 'NEXTFLOW'),
    'code': (None, pipeline_name),
    'description': (None, f'Pipeline exported on {datetime.date.today().isoformat()}: {current_git_commit_message}'),
    'defaultStorageType': (None, 'Small'),
    'proprietary': (None, 'false'),
    'version': (None, current_git_commit_id_short),
    'gitCredentialId': (None, git_credential_uuid),
    'commitId': (None, current_git_commit_id),
    'repositoryUrl': (None, repository_url),
    'mainFilePath': (None, main_file_path),
    'configFilePath': (None, nextflow_config_path),
}

response = requests.post(f'{api_url}/projects/{project_id}/pipelines:importGitPipeline', headers=headers, files=files)
print(f'API response status code: {response.status_code}')
print(f'API response body:')
print(json.dumps(response.json(), indent=2))
if not response.ok:
    exit(1)
pipeline_id = response.json()['id']

# The git pipeline import runs asynchronously (status starts as 'Importing'). The input form
# can only be uploaded once ICA has finished parsing the repo and the pipeline reaches 'Draft'.
poll_interval_s = 5
max_attempts = 120  # 10 minutes
FAILURE_STATUSES = {'Import Failed', 'Import Incomplete', 'Import Cancelling'}

print('')
print(f'Waiting for pipeline {pipeline_id} to reach Draft status...')
status = None
for attempt in range(1, max_attempts + 1):
    poll_response = requests.get(f'{api_url}/projects/{project_id}/pipelines/{pipeline_id}', headers=headers, timeout=30)
    poll_response.raise_for_status()
    status = poll_response.json()['pipeline']['statusAsString']
    print(f'  [{attempt}/{max_attempts}] status: {status}')
    if status == 'Draft':
        break
    if status in FAILURE_STATUSES:
        print(f'ERROR: pipeline import ended in status "{status}"; not uploading input form.')
        exit(1)
    time.sleep(poll_interval_s)
else:
    print(f'ERROR: pipeline did not reach Draft status within {max_attempts * poll_interval_s}s (last status: {status}); not uploading input form.')
    exit(1)

print('')
print(f'Uploading input form from {input_form_path}...')
# The inputForm/inputFormFile endpoint only recognizes the v3 API version -- sending the v4
# Accept header used for the other endpoints above gets rejected with "Invalid Accept Header".
input_form_headers = {**headers, 'Accept': 'application/vnd.illumina.v3+json'}
with open(input_form_path, 'rb') as input_form_file:
    form_response = requests.put(
        f'{api_url}/projects/{project_id}/pipelines/{pipeline_id}/inputForm/inputFormFile',
        headers=input_form_headers,
        files={'content': ('inputForm.json', input_form_file, 'application/json')},
    )
print(f'Input form upload status code: {form_response.status_code}')
if form_response.content:
    print(json.dumps(form_response.json(), indent=2))
if not form_response.ok:
    exit(1)
print('Done.')
