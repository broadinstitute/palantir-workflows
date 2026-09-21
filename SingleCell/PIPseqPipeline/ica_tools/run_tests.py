#!/usr/bin/env python3
"""
Runs main.nf's test/test_inputs_main.json against the dedicated "BCL Pipeline Testing" ICA
project, exporting a new pipeline for the current git commit first if one isn't already
imported there.

Looks up already-imported pipelines by the exact code export_pipeline_to_ica.py would have
given the current commit (e.g. "PIPseq_BCL_abc1234") and reuses one if found; otherwise imports
one via export_pipeline_to_ica.export_pipeline(). Since ICA's git-backed pipeline import fetches
code from GitHub at a specific commit id -- never from the local working tree -- this warns and
asks for confirmation before exporting if there are uncommitted local changes, since those
changes will *not* be part of whatever gets tested.

Like start_analysis.py, submits the run and exits immediately (fire-and-forget); it doesn't wait
for the analysis to finish.

Usage:
    python3 run_tests.py [--dry-run]
"""

import argparse
import json
import os
import subprocess
import sys

from ica_common import PROJECT_NAMES_AND_IDS, load_api_key
from export_pipeline_to_ica import ENTRYPOINTS, export_pipeline, get_current_git_commit_info, pipeline_code_for
from start_analysis import build_input_form_fields, list_pipseq_pipelines, load_input_form_field_types, start_analysis

script_dir = os.path.dirname(os.path.abspath(__file__))
pipeline_root = os.path.dirname(script_dir)

TEST_PROJECT_NAME = 'BCL Pipeline Testing'
ENTRYPOINT_KEY = 'Full pipeline (main.nf, --fastq_list)'
TEST_INPUT_FILE = os.path.join(pipeline_root, 'test', 'test_inputs_main.json')


def has_uncommitted_changes():
    status = subprocess.check_output(['git', 'status', '--porcelain'], cwd=pipeline_root).decode()
    return bool(status.strip())


def confirm(prompt):
    return input(f'{prompt} [y/N]: ').strip().lower() == 'y'


def find_existing_pipeline(api_key, project_id, expected_code):
    for pipeline in list_pipseq_pipelines(api_key, project_id):
        if pipeline['code'] == expected_code:
            return pipeline
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dry-run', action='store_true', help='Resolve inputs and print the request payload without exporting a pipeline or starting the analysis')
    args = parser.parse_args()

    project_id = PROJECT_NAMES_AND_IDS[TEST_PROJECT_NAME]
    api_key = load_api_key()
    commit_id, commit_id_short, commit_message = get_current_git_commit_info()
    expected_code = pipeline_code_for(ENTRYPOINTS[ENTRYPOINT_KEY], commit_id_short)

    print(f'Looking for pipeline "{expected_code}" in project "{TEST_PROJECT_NAME}"...')
    pipeline = find_existing_pipeline(api_key, project_id, expected_code)

    if pipeline is not None:
        pipeline_id = pipeline['id']
        print(f'Found existing pipeline {pipeline["code"]} ({pipeline_id}).')
        if args.dry_run:
            print(f'Dry run -- would run pipeline {pipeline["code"]} ({pipeline_id}) with inputs from {TEST_INPUT_FILE}.')
            return
    else:
        print(f'No existing pipeline found for commit {commit_id_short}.')
        if has_uncommitted_changes():
            print('')
            print('WARNING: you have uncommitted local changes. ICA\'s git-backed pipeline import')
            print('fetches code from GitHub at a specific commit -- it will NOT see your uncommitted')
            print(f'changes, only what was actually pushed as commit {commit_id_short}.')
            print('')
            if not args.dry_run and not confirm(f'Export and test commit {commit_id_short} anyway (ignoring uncommitted changes)?'):
                print('Aborted.')
                sys.exit(0)

        if args.dry_run:
            print(f'Dry run -- would export "{ENTRYPOINT_KEY}" at commit {commit_id_short} to project '
                  f'"{TEST_PROJECT_NAME}" ({project_id}), then run it with inputs from {TEST_INPUT_FILE}.')
            return

        print(f'Exporting pipeline for commit {commit_id_short} to project "{TEST_PROJECT_NAME}"...')
        try:
            pipeline_id = export_pipeline(api_key, project_id, ENTRYPOINT_KEY, commit_id, commit_id_short, commit_message)
        except RuntimeError as e:
            print(f'ERROR: {e}')
            sys.exit(1)

    with open(TEST_INPUT_FILE) as f:
        inputs = json.load(f)

    print('')
    try:
        field_types = load_input_form_field_types(api_key, project_id, pipeline_id)
        print(f'Resolving {len(inputs)} input field(s) against project {project_id}...')
        fields = build_input_form_fields(api_key, project_id, field_types, inputs)
    except ValueError as e:
        print(f'ERROR: {e}')
        sys.exit(1)

    user_reference = f'pipseq_test_run_{commit_id_short}'
    print(f'Starting analysis "{user_reference}" in project "{TEST_PROJECT_NAME}"...')
    response = start_analysis(api_key, project_id, pipeline_id, user_reference, fields)
    if not response.ok:
        print(f'ERROR: failed to start analysis (HTTP {response.status_code}): {response.json().get("detail", response.text)}')
        sys.exit(1)
    print(f'Analysis started successfully (HTTP {response.status_code}). Analysis ID: {response.json()["id"]}')


if __name__ == '__main__':
    main()
