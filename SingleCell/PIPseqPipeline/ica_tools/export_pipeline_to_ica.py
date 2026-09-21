#!/usr/bin/env python3
"""
Imports the current git commit of this pipeline into an ICA project as a git-backed Nextflow
pipeline, then uploads the matching hand-maintained input form.

Interactively prompts for which entrypoint (main.nf or main_simple.nf) and which ICA project.
Also importable as a module (see export_pipeline() below) -- run_tests.py uses this to export a
pipeline on demand when no already-imported one matches the current commit.

Usage:
    python3 export_pipeline_to_ica.py
"""

import requests
import subprocess
import datetime
import os
import sys
import time

from ica_common import API_URL, PROJECT_NAMES_AND_IDS, load_api_key, prompt_choice

script_dir = os.path.dirname(os.path.abspath(__file__))
pipeline_root = os.path.dirname(script_dir)

ENTRYPOINTS = {
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

REPOSITORY_URL = 'https://github.com/broadinstitute/palantir-workflows'
# Both entrypoints share the same process/resource config regardless of which is exported.
NEXTFLOW_CONFIG_PATH = 'SingleCell/PIPseqPipeline/nextflow.config'
GIT_CREDENTIAL_UUID = '5a2282d8-61a7-4222-8969-bfefbbe4f949'

# The git pipeline import runs asynchronously (status starts as 'Importing'). The input form
# can only be uploaded once ICA has finished parsing the repo and the pipeline reaches 'Draft'.
POLL_INTERVAL_S = 5
MAX_POLL_ATTEMPTS = 120  # 10 minutes
FAILURE_STATUSES = {'Import Failed', 'Import Incomplete', 'Import Cancelling'}


def get_current_git_commit_info():
    commit_id = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=pipeline_root).decode().strip()
    commit_id_short = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=pipeline_root).decode().strip()
    commit_message = subprocess.check_output(['git', 'log', '-1', '--pretty=%B'], cwd=pipeline_root).decode().strip().split('\n')[0].strip()
    return commit_id, commit_id_short, commit_message


def pipeline_code_for(entrypoint, commit_id_short):
    return f'PIPseq_BCL{entrypoint["name_suffix"]}_{commit_id_short}'


def import_git_pipeline(api_key, project_id, entrypoint, commit_id, commit_id_short, commit_message):
    pipeline_name = pipeline_code_for(entrypoint, commit_id_short)

    print(f'Exporting pipeline with the following data:')
    print(f'  Pipeline name: {pipeline_name}')
    print(f'  Pipeline version: {commit_id_short}')
    print(f'  Commit: {commit_id}')
    print(f'  Repository URL: {REPOSITORY_URL}')
    print(f'  Main file path: {entrypoint["main_file_path"]}')
    print(f'  Nextflow config path: {NEXTFLOW_CONFIG_PATH}')
    print(f'  Git credential UUID: {GIT_CREDENTIAL_UUID}')
    print('')

    headers = {
        'X-API-Key': api_key,
        'Accept': 'application/vnd.illumina.v4+json',
    }
    # For multipart/form-data, use the files parameter with (None, value) tuples
    # This forces requests to send as multipart/form-data instead of application/x-www-form-urlencoded
    files = {
        'language': (None, 'NEXTFLOW'),
        'code': (None, pipeline_name),
        'description': (None, f'Pipeline exported on {datetime.date.today().isoformat()}: {commit_message}'),
        'defaultStorageType': (None, 'Small'),
        'proprietary': (None, 'false'),
        'version': (None, commit_id_short),
        'gitCredentialId': (None, GIT_CREDENTIAL_UUID),
        'commitId': (None, commit_id),
        'repositoryUrl': (None, REPOSITORY_URL),
        'mainFilePath': (None, entrypoint['main_file_path']),
        'configFilePath': (None, NEXTFLOW_CONFIG_PATH),
    }

    response = requests.post(f'{API_URL}/projects/{project_id}/pipelines:importGitPipeline', headers=headers, files=files)
    if not response.ok:
        raise RuntimeError(f'import failed (HTTP {response.status_code}): {response.json().get("detail", response.text)}')
    pipeline_id = response.json()['id']
    print(f'Import scheduled successfully (HTTP {response.status_code}). Pipeline ID: {pipeline_id}')
    return pipeline_id


def wait_for_draft(api_key, project_id, pipeline_id):
    headers = {
        'X-API-Key': api_key,
        'Accept': 'application/vnd.illumina.v4+json',
    }

    print('')
    print(f'Waiting for pipeline {pipeline_id} to reach Draft status...')
    status = None
    for attempt in range(1, MAX_POLL_ATTEMPTS + 1):
        poll_response = requests.get(f'{API_URL}/projects/{project_id}/pipelines/{pipeline_id}', headers=headers, timeout=30)
        poll_response.raise_for_status()
        status = poll_response.json()['pipeline']['statusAsString']
        print(f'  [{attempt}/{MAX_POLL_ATTEMPTS}] status: {status}')
        if status == 'Draft':
            return
        if status in FAILURE_STATUSES:
            raise RuntimeError(f'pipeline import ended in status "{status}"; not uploading input form.')
        time.sleep(POLL_INTERVAL_S)
    raise RuntimeError(f'pipeline did not reach Draft status within {MAX_POLL_ATTEMPTS * POLL_INTERVAL_S}s (last status: {status}); not uploading input form.')


def upload_input_form(api_key, project_id, pipeline_id, input_form_path):
    print('')
    print(f'Uploading input form from {input_form_path}...')
    # The inputForm/inputFormFile endpoint only recognizes the v3 API version -- sending the v4
    # Accept header used for the other endpoints above gets rejected with "Invalid Accept Header".
    headers = {
        'X-API-Key': api_key,
        'Accept': 'application/vnd.illumina.v3+json',
    }
    with open(input_form_path, 'rb') as input_form_file:
        response = requests.put(
            f'{API_URL}/projects/{project_id}/pipelines/{pipeline_id}/inputForm/inputFormFile',
            headers=headers,
            files={'content': ('inputForm.json', input_form_file, 'application/json')},
        )
    if not response.ok:
        raise RuntimeError(f'input form upload failed (HTTP {response.status_code}): {response.json().get("detail", response.text)}')
    print(f'Input form uploaded successfully (HTTP {response.status_code}).')


def export_pipeline(api_key, project_id, entrypoint_key, commit_id, commit_id_short, commit_message):
    """Imports entrypoint_key's given git commit into project_id, waits for it to reach Draft,
    then uploads its input form. Returns the new pipeline_id."""
    entrypoint = ENTRYPOINTS[entrypoint_key]
    input_form_path = os.path.join(script_dir, 'inputforms', entrypoint['input_form_dir'], 'inputForm.json')
    if not os.path.isfile(input_form_path):
        raise RuntimeError(f'input form file not found at {input_form_path}')

    pipeline_id = import_git_pipeline(api_key, project_id, entrypoint, commit_id, commit_id_short, commit_message)
    wait_for_draft(api_key, project_id, pipeline_id)
    upload_input_form(api_key, project_id, pipeline_id, input_form_path)
    print('Done.')
    return pipeline_id


def main():
    api_key = load_api_key()
    commit_id, commit_id_short, commit_message = get_current_git_commit_info()

    entrypoint_names = list(ENTRYPOINTS.keys())
    entrypoint_choice = prompt_choice('Which entrypoint do you want to export?', entrypoint_names)
    if entrypoint_choice is None:
        sys.exit(0)
    entrypoint_key = entrypoint_names[entrypoint_choice]

    project_names = list(PROJECT_NAMES_AND_IDS.keys())
    project_choice = prompt_choice('Which project do you want to export to?', project_names)
    if project_choice is None:
        sys.exit(0)
    project_id = PROJECT_NAMES_AND_IDS[project_names[project_choice]]

    try:
        export_pipeline(api_key, project_id, entrypoint_key, commit_id, commit_id_short, commit_message)
    except RuntimeError as e:
        print(f'ERROR: {e}')
        sys.exit(1)


if __name__ == '__main__':
    main()
