import requests
import subprocess
import json

api_url = 'https://ica.illumina.com/ica/rest/api'

ica_api_key_filename = '/Users/mgatzen/.icav2/api_key.txt'
ica_api_key = open(ica_api_key_filename).read().strip()

project_names_and_ids = {
    'BCL Shared Development': '4ab33fe6-c169-4dc9-928d-6ce7a8062d34',
    'MDL Single Cell Dev': 'de460091-a137-4848-8e23-c59b851d7425',
}

current_git_commit_id = subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip()
current_git_commit_id_short = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode().strip()

pipeline_name = f'PIPseq_BCL_{current_git_commit_id_short}'

repository_url = 'https://github.com/broadinstitute/palantir-workflows'
main_file_path = 'SingleCell/PIPseqPipeline/main.nf'
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
print('Which project do you want to export to?')

project_names = list(project_names_and_ids.keys())
for i, project_name in enumerate(project_names, start=1):
    print(f'  {i}. {project_name}')
project_choice = input('Enter the number of the project (or anything else to abort): ')
if not project_choice.isdigit() or int(project_choice) < 1 or int(project_choice) > len(project_names):
    exit(0)
project_id = project_names_and_ids[project_names[int(project_choice) - 1]]

headers = {
    'X-API-Key': ica_api_key,
    'Accept': 'application/vnd.illumina.v4+json',
}

# For multipart/form-data, use the files parameter with (None, value) tuples
# This forces requests to send as multipart/form-data instead of application/x-www-form-urlencoded
files = {
    'language': (None, 'NEXTFLOW'),
    'code': (None, pipeline_name),
    'description': (None, f'PIPseq pipeline exported from git commit {current_git_commit_id} of repository {repository_url}'),
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
