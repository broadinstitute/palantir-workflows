#!/usr/bin/env python3
"""
Interactively kicks off an ICA analysis run for a PIPseq pipeline already imported via
export_pipeline_to_ica.py.

Prompts for which ICA project to run in, which already-imported PIPseq pipeline to run
(flagging whichever one's commit matches the current git HEAD), and which of the
test/test_inputs_main.json / test/test_inputs_main_simple.json input files to submit.

File/folder inputs in those JSON files are given as project-relative paths (matching how
they appear in the ICA project's data tree) rather than ICA data IDs -- this script resolves
each one to its `fil.<hash>`/`fol.<hash>` ID by fetching the chosen pipeline's live input form
from ICA (which is what already declares which fields are data fields) before submitting.

Usage:
    python3 start_analysis.py [--dry-run]
"""

import argparse
import json
import os
import re
import subprocess
import sys

import requests

from ica_common import API_URL, PROJECT_NAMES_AND_IDS, load_api_key, prompt_choice

script_dir = os.path.dirname(os.path.abspath(__file__))
pipeline_root = os.path.dirname(script_dir)

DATA_ID_PATTERN = re.compile(r'^(fil|fol)\.[0-9a-f]+$')
PIPELINE_NAME_PREFIX = 'PIPseq_BCL'

TEST_INPUT_FILES = {
    'main.nf (test/test_inputs_main.json)': os.path.join(pipeline_root, 'test', 'test_inputs_main.json'),
    'main_simple.nf (test/test_inputs_main_simple.json)': os.path.join(pipeline_root, 'test', 'test_inputs_main_simple.json'),
}


def get_current_commit_id():
    return subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=pipeline_root).decode().strip()


def list_pipseq_pipelines(api_key, project_id):
    headers = {
        'X-API-Key': api_key,
        'Accept': 'application/vnd.illumina.v4+json',
    }
    response = requests.get(f'{API_URL}/projects/{project_id}/pipelines', headers=headers)
    response.raise_for_status()
    pipelines = [item['pipeline'] for item in response.json()['items'] if item['pipeline']['code'].startswith(PIPELINE_NAME_PREFIX)]
    pipelines.sort(key=lambda p: p['timeCreated'], reverse=True)
    return pipelines


def load_input_form_field_types(api_key, project_id, pipeline_id):
    headers = {
        'X-API-Key': api_key,
        'Accept': 'application/octet-stream',
    }
    response = requests.get(
        f'{API_URL}/projects/{project_id}/pipelines/{pipeline_id}/inputForm/inputFormFile',
        headers=headers,
    )
    if not response.ok:
        raise ValueError(f'Failed to fetch input form for pipeline {pipeline_id} in project {project_id}: {response.status_code} {response.text}')
    fields = response.json()['fields']
    return {field['id']: field['type'] for field in fields if field['type'] != 'section'}


def resolve_data_id(api_key, project_id, path):
    if DATA_ID_PATTERN.match(path):
        return path
    if not path.startswith('/'):
        path = '/' + path
    headers = {
        'X-API-Key': api_key,
        'Accept': 'application/vnd.illumina.v3+json',
    }
    response = requests.get(
        f'{API_URL}/projects/{project_id}/data',
        headers=headers,
        params={
            'filePath': path,
            'filePathMatchMode': 'FULL_CASE_INSENSITIVE',
            'type': 'FILE',
        },
    )
    if not response.ok:
        raise ValueError(f'Failed to look up path "{path}" in project {project_id}: {response.status_code} {response.text}')
    items = response.json()['items']
    if len(items) == 0:
        raise ValueError(f'No file found at path "{path}" in project {project_id}')
    if len(items) > 1:
        candidates = ', '.join(item['data']['id'] for item in items)
        raise ValueError(f'Path "{path}" matched multiple files in project {project_id}: {candidates}')
    return items[0]['data']['id']


def to_value_list(value):
    return [str(v).lower() if isinstance(v, bool) else str(v) for v in (value if isinstance(value, list) else [value])]


def build_input_form_fields(api_key, project_id, field_types, inputs):
    unknown_fields = set(inputs) - set(field_types)
    if unknown_fields:
        raise ValueError(f'Unknown input field(s) not present in the input form: {sorted(unknown_fields)}')

    fields = []
    for field_id, value in inputs.items():
        if field_types[field_id] == 'data':
            paths = value if isinstance(value, list) else [value]
            data_ids = [resolve_data_id(api_key, project_id, path) for path in paths]
            fields.append({
                'id': field_id,
                'dataValues': [{'dataId': data_id} for data_id in data_ids],
            })
        else:
            fields.append({
                'id': field_id,
                'values': to_value_list(value),
            })
    return fields


def start_analysis(api_key, project_id, pipeline_id, user_reference, fields):
    headers = {
        'X-API-Key': api_key,
        'Content-Type': 'application/vnd.illumina.v4+json',
        'Accept': 'application/vnd.illumina.v4+json',
    }
    body = {
        'userReference': user_reference,
        'pipelineId': pipeline_id,
        'inputFormValues': {
            'fields': fields,
        },
    }
    return requests.post(f'{API_URL}/projects/{project_id}/analysis:nextflowJson', headers=headers, json=body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true', help='Resolve data paths to IDs and print the request payload without starting the analysis')
    args = parser.parse_args()

    api_key = load_api_key()
    current_commit_id = get_current_commit_id()

    project_names = list(PROJECT_NAMES_AND_IDS.keys())
    project_choice = prompt_choice('Which ICA project do you want to start the test run in?', project_names)
    if project_choice is None:
        sys.exit(0)
    project_id = PROJECT_NAMES_AND_IDS[project_names[project_choice]]

    pipelines = list_pipseq_pipelines(api_key, project_id)
    if not pipelines:
        print(f'No PIPseq pipelines found in project "{project_names[project_choice]}" -- export one first with export_pipeline_to_ica.py.')
        sys.exit(1)
    pipeline_labels = [
        f"{p['code']} [{p['statusAsString']}]" + (' <- matches current HEAD commit' if p.get('gitPipelineImportDto', {}).get('commitId') == current_commit_id else '')
        for p in pipelines
    ]
    pipeline_choice = prompt_choice('Which pipeline do you want to run?', pipeline_labels)
    if pipeline_choice is None:
        sys.exit(0)
    pipeline = pipelines[pipeline_choice]
    pipeline_id = pipeline['id']

    input_file_names = list(TEST_INPUT_FILES.keys())
    input_file_choice = prompt_choice('Which input file do you want to use?', input_file_names)
    if input_file_choice is None:
        sys.exit(0)
    with open(TEST_INPUT_FILES[input_file_names[input_file_choice]]) as f:
        inputs = json.load(f)

    print('')
    try:
        field_types = load_input_form_field_types(api_key, project_id, pipeline_id)
        print(f'Resolving {len(inputs)} input field(s) against project {project_id}...')
        fields = build_input_form_fields(api_key, project_id, field_types, inputs)
    except ValueError as e:
        print(f'ERROR: {e}')
        sys.exit(1)

    user_reference = f"pipseq_test_run_{current_commit_id[:7]}"

    if args.dry_run:
        print('Dry run -- would submit the following analysis:')
        print(json.dumps({
            'projectId': project_id,
            'pipelineId': pipeline_id,
            'userReference': user_reference,
            'inputFormValues': {'fields': fields},
        }, indent=2))
        return

    print(f'Starting analysis "{user_reference}" for pipeline {pipeline["code"]} in project {project_names[project_choice]}...')
    response = start_analysis(api_key, project_id, pipeline_id, user_reference, fields)
    if not response.ok:
        print(f'ERROR: failed to start analysis (HTTP {response.status_code}): {response.json().get("detail", response.text)}')
        sys.exit(1)
    print(f'Analysis started successfully (HTTP {response.status_code}). Analysis ID: {response.json()["id"]}')


if __name__ == '__main__':
    main()
