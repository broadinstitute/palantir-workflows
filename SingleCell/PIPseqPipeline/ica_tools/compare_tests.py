#!/usr/bin/env python3
"""
Downloads the output files matching test/compare/file_patterns.txt from two completed ICA
analyses in the "BCL Pipeline Testing" project and diffs them, to check whether a pipeline
change had any unintended effect on pipeline outputs.

Usage:
    python3 compare_tests.py <base_analysis_id> <eval_analysis_id>

Writes one <relative-path>.diff file per differing output file (mirroring that file's path
under its analysis's output folder) plus a summary.txt, to:
    test/compare/diffs/{base_pipeline_version}_vs_{eval_pipeline_version}/
Exits non-zero if any matched file differs or is missing from one side.
"""

import argparse
import difflib
import fnmatch
import os
import shutil
import sys
import tempfile

import requests

from ica_common import API_URL, PROJECT_NAMES_AND_IDS, load_api_key

script_dir = os.path.dirname(os.path.abspath(__file__))
pipeline_root = os.path.dirname(script_dir)

PROJECT_NAME = 'BCL Pipeline Testing'
PATTERNS_FILE = os.path.join(pipeline_root, 'test', 'compare', 'file_patterns.txt')
DIFFS_DIR = os.path.join(pipeline_root, 'test', 'compare', 'diffs')


def load_patterns():
    with open(PATTERNS_FILE) as f:
        return [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]


def matches_any_pattern(filename, patterns):
    return any(fnmatch.fnmatch(filename, pattern) for pattern in patterns)


def get_analysis(api_key, project_id, analysis_id):
    headers = {'X-API-Key': api_key, 'Accept': 'application/vnd.illumina.v4+json'}
    response = requests.get(f'{API_URL}/projects/{project_id}/analyses/{analysis_id}', headers=headers)
    if not response.ok:
        raise RuntimeError(f'failed to fetch analysis {analysis_id} (HTTP {response.status_code}): {response.text}')
    return response.json()


def get_output_root_paths(api_key, project_id, analysis_id):
    """Returns the project-relative output folder path(s) an analysis published its results under."""
    headers = {'X-API-Key': api_key, 'Accept': 'application/vnd.illumina.v3+json'}
    response = requests.get(f'{API_URL}/projects/{project_id}/analyses/{analysis_id}/outputs', headers=headers)
    if not response.ok:
        raise RuntimeError(f'failed to fetch outputs of analysis {analysis_id} (HTTP {response.status_code}): {response.text}')
    data_ids = [entry['dataId'] for group in response.json()['items'] for entry in group['data']]

    paths = []
    for data_id in data_ids:
        response = requests.get(f'{API_URL}/projects/{project_id}/data/{data_id}', headers=headers)
        response.raise_for_status()
        paths.append(response.json()['data']['details']['path'])
    return paths


def list_files_under(api_key, project_id, root_path):
    """Recursively lists every FILE (not folder) whose path starts with root_path, as (dataId, path) pairs."""
    headers = {'X-API-Key': api_key, 'Accept': 'application/vnd.illumina.v3+json'}
    files = []
    page_token = None
    while True:
        params = {'filePath': root_path, 'filePathMatchMode': 'STARTS_WITH_CASE_INSENSITIVE', 'pageSize': 1000}
        if page_token:
            params['pageToken'] = page_token
        response = requests.get(f'{API_URL}/projects/{project_id}/data', headers=headers, params=params)
        response.raise_for_status()
        page = response.json()
        for item in page['items']:
            d = item['data']
            if d['details']['dataType'] == 'FILE':
                files.append((d['id'], d['details']['path']))
        if page.get('remainingRecords', 0) <= 0:
            break
        page_token = page['nextPageToken']
    return files


def download_file(api_key, project_id, data_id, dest_path):
    headers = {'X-API-Key': api_key, 'Accept': 'application/vnd.illumina.v3+json'}
    response = requests.post(f'{API_URL}/projects/{project_id}/data/{data_id}:createDownloadUrl', headers=headers)
    response.raise_for_status()
    # The presigned S3 URL requires no auth headers of its own.
    download_response = requests.get(response.json()['url'])
    download_response.raise_for_status()

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with open(dest_path, 'wb') as f:
        f.write(download_response.content)


def fetch_matching_outputs(api_key, project_id, analysis_id, patterns, download_dir):
    """Downloads every output file of analysis_id matching patterns into download_dir.

    Returns {path_relative_to_output_root: local_path}.
    """
    downloaded = {}
    for root_path in get_output_root_paths(api_key, project_id, analysis_id):
        for data_id, file_path in list_files_under(api_key, project_id, root_path):
            if not matches_any_pattern(os.path.basename(file_path), patterns):
                continue
            relative_path = os.path.relpath(file_path, root_path)
            local_path = os.path.join(download_dir, relative_path)
            print(f'  downloading {file_path}')
            download_file(api_key, project_id, data_id, local_path)
            downloaded[relative_path] = local_path
    return downloaded


def diff_files(base_path, eval_path):
    """Returns a unified diff string, or None if the two files are identical."""
    try:
        with open(base_path, encoding='utf-8') as f:
            base_lines = f.readlines()
        with open(eval_path, encoding='utf-8') as f:
            eval_lines = f.readlines()
    except UnicodeDecodeError:
        with open(base_path, 'rb') as f:
            base_bytes = f.read()
        with open(eval_path, 'rb') as f:
            eval_bytes = f.read()
        return None if base_bytes == eval_bytes else '(binary files differ)\n'

    diff = list(difflib.unified_diff(base_lines, eval_lines, fromfile='base', tofile='eval'))
    return ''.join(diff) if diff else None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('base', help='ICA analysis ID to treat as the baseline')
    parser.add_argument('eval', help='ICA analysis ID to compare against the baseline')
    args = parser.parse_args()

    project_id = PROJECT_NAMES_AND_IDS[PROJECT_NAME]
    api_key = load_api_key()
    patterns = load_patterns()

    try:
        base_analysis = get_analysis(api_key, project_id, args.base)
        eval_analysis = get_analysis(api_key, project_id, args.eval)
    except RuntimeError as e:
        print(f'ERROR: {e}')
        sys.exit(1)

    # export_pipeline_to_ica.py stores the short commit hash it exported as the pipeline's version.
    base_version = base_analysis['pipeline']['version']
    eval_version = eval_analysis['pipeline']['version']

    print(f'Base analysis {args.base}: pipeline version {base_version}, status {base_analysis["status"]}')
    print(f'Eval analysis {args.eval}: pipeline version {eval_version}, status {eval_analysis["status"]}')
    for label, analysis in (('base', base_analysis), ('eval', eval_analysis)):
        if analysis['status'] != 'SUCCEEDED':
            print(f'WARNING: {label} analysis status is "{analysis["status"]}", not "SUCCEEDED" -- outputs may be incomplete.')

    output_dir = os.path.join(DIFFS_DIR, f'{base_version}_vs_{eval_version}')
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    with tempfile.TemporaryDirectory(prefix='pipseq_compare_') as tmp_dir:
        print(f'\nDownloading matching outputs for base analysis {args.base}...')
        base_files = fetch_matching_outputs(api_key, project_id, args.base, patterns, os.path.join(tmp_dir, 'base'))
        print(f'\nDownloading matching outputs for eval analysis {args.eval}...')
        eval_files = fetch_matching_outputs(api_key, project_id, args.eval, patterns, os.path.join(tmp_dir, 'eval'))

        all_relative_paths = sorted(set(base_files) | set(eval_files))

        differing = []
        only_in_base = []
        only_in_eval = []
        identical_count = 0

        for relative_path in all_relative_paths:
            if relative_path not in base_files:
                only_in_eval.append(relative_path)
                continue
            if relative_path not in eval_files:
                only_in_base.append(relative_path)
                continue

            diff_text = diff_files(base_files[relative_path], eval_files[relative_path])
            if diff_text is None:
                identical_count += 1
                continue

            differing.append(relative_path)
            diff_dest = os.path.join(output_dir, relative_path + '.diff')
            os.makedirs(os.path.dirname(diff_dest), exist_ok=True)
            with open(diff_dest, 'w') as f:
                f.write(diff_text)

        summary_lines = [
            f'base analysis: {args.base} (pipeline version {base_version}, status {base_analysis["status"]})',
            f'eval analysis: {args.eval} (pipeline version {eval_version}, status {eval_analysis["status"]})',
            '',
            f'{len(all_relative_paths)} matching file(s) found across both runs',
            f'{identical_count} identical',
            f'{len(differing)} differing',
            f'{len(only_in_base)} only in base',
            f'{len(only_in_eval)} only in eval',
        ]
        if differing:
            summary_lines += ['', 'Differing files:'] + [f'  {p}' for p in differing]
        if only_in_base:
            summary_lines += ['', 'Only in base:'] + [f'  {p}' for p in only_in_base]
        if only_in_eval:
            summary_lines += ['', 'Only in eval:'] + [f'  {p}' for p in only_in_eval]

        with open(os.path.join(output_dir, 'summary.txt'), 'w') as f:
            f.write('\n'.join(summary_lines) + '\n')

    print('')
    print('\n'.join(summary_lines))
    print(f'\nResults written to {output_dir}')

    if differing or only_in_base or only_in_eval:
        sys.exit(1)


if __name__ == '__main__':
    main()
