# MultiQC in Terra

This notebook will search through the submissions in a Terra workspace for QC files, download them to the vm, run [MultiQC](https://multiqc.info/), upload the resulting reports into the google bucket associated with this workspace, and provide a link to those reports.  Workflows will be grouped according to workflow name, and a different report will be generated for each workflow name.

This notebook must be run on a VM which has included the starup script `gs://broad-dsde-methods-ckachulis/MultiQC_in_Terra/multiqc_terra_startup.sh`.  If you don't have access to that file, you can use your own by uploading the following as a bash script into a google bucket of your choice.
```shell
#!/usr/bin/env bash

pip install git+https://github.com/kachulis/MultiQC.git@ck_gcp
```

If you can't or really don't want to use a startup script, you can uncomment the next cell, and run it instead.  Note that this cell will kill the kernel and cause it to restart.  This is not a bug, and is necessary for annoying package version reasons, and can be avoided by using the startup script instead.


```python
# !pip install git+https://github.com/kachulis/MultiQC.git@ck_gcp
# import os
# os._exit(00)
```

You can set the workspace you want to run multiqc on (by default it will run on this workspace), and any filters in the next cell.  Note that filters are cumulative.  If a workflow fails any of the filters, it will not be included.  A particular filter can be turned off by setting it to `None`.


```python
import firecloud.api as fapi
import os
namespace = os.environ.get("WORKSPACE_NAMESPACE")
workspace = os.environ.get("WORKSPACE_NAME")

submission_ids_to_include = None
worklow_ids_to_include = None
workflow_names_to_include = None
```

The next cell will find the appropriate QC files, download them, run MultiQC, and upload the resulting reports


```python
import asyncio
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from multiqc.utils import report, config
import os
from google.cloud import storage
storage_client = storage.Client()
import subprocess
import pathlib
import pytz
from datetime import datetime

config.mqc_load_userconfig()

def flat(pool):
    res = []
    for v in pool:
        if isinstance(v, list):
            res += flat(v)
        elif isinstance(v, str) or isinstance(v, tuple):
                res.append(v)
        else:
            raise RuntimeError(f'Trying to flatten and found object of unexpected type {type(v)}: {v}')
    return res  

#get submissions
print("Finding submissions...")
submissions=fapi.list_submissions(namespace, workspace).json()
submission_ids = [submission['submissionId'] for submission in submissions]
if submission_ids_to_include is not None:
    submission_ids = [sid for sid in submission_ids if sid in submission_ids_to_include]
print(f'Found {len(submission_ids)} submissions.')

#get list of tuples of workflow_ids that have succeeded and submission_ids
print("Finding workflows...")
def get_workflow_ids_for_submission(submission_id):
    submission = fapi.get_submission(namespace, workspace, submission_id).json()
    workflows = fapi.get_submission(namespace, workspace, submission_id).json()['workflows']
    workflow_ids = [(submission_id, workflow['workflowId'], submission['submissionDate']) for workflow in workflows if workflow['status'] == 'Succeeded']
    return workflow_ids

with ThreadPoolExecutor(max_workers=50) as executor:
    loop = asyncio.get_event_loop()
    tasks = [loop.run_in_executor(executor, get_workflow_ids_for_submission, submission_id) for submission_id in submission_ids]
    workflow_with_submissions = [id_tuple for id_tuples in await asyncio.gather(*tasks) for id_tuple in id_tuples]

if worklow_ids_to_include is not None:
    workflow_with_submissions = [e for e in workflow_with_sumbissions if e[1] in worklow_ids_to_include]
print(f'Found {len(workflow_with_submissions)} workflows.')


# build dictionary of outputs for specified worklows
print("Finding outputs...")
def get_outputs(submission_id, workflow_id, submission_date):
    metadata = fapi.get_workflow_metadata(namespace, workspace, submission_id, workflow_id).json()
    res = {k:(v, submission_date) for k,v in metadata['outputs'].items() if isinstance(v, list) or isinstance(v, str) and v.startswith("gs://")}
#     res = dict(filter(outputs_filter, metadata['outputs'].items()))
    return res

with ThreadPoolExecutor(max_workers=50) as executor:
    loop = asyncio.get_event_loop()
    tasks = [loop.run_in_executor(executor, get_outputs, submission_id, workflow_id, submission_date) for 
             submission_id, workflow_id, submission_date in workflow_with_submissions]
    outputs_dict = defaultdict(list)
    for d in await asyncio.gather(*tasks):
        for key, value in d.items():
            outputs_dict[key].append(value)

outputs_dict = {k: flat(v) for k,v in outputs_dict.items()}
if workflow_names_to_include is not None:
    outputs_dict = {k:v for k,v in outputs_dict.items() if k.split(".")[0] in workflow_names_to_include}
print(f'Found {sum([len(v) for v in outputs_dict.values()])} outputs')

# for each output key, check whether the output files are useful for multiqc
print("Checking outputs against MultiQC")
def check_output(key, v):
    is_found_multiqc = report.search_gcs(v)
    return key, is_found_multiqc

with ThreadPoolExecutor(max_workers=50) as executor:
    loop = asyncio.get_event_loop()
    tasks = [loop.run_in_executor(executor, check_output, output_name, output_list[0][0]) for 
             output_name, output_list in outputs_dict.items()]
    status_dict = {k:v for k,v in await asyncio.gather(*tasks)}
metrics_file_paths = flat([(k.split(".")[0], v2[0], v2[1]) for k,v in outputs_dict.items() if status_dict[k] for v2 in v])
print(f'Found {len(metrics_file_paths)} metrics files.')


workflow_names = {wn for wn, _, _ in metrics_file_paths}
now = datetime.now(pytz.timezone('US/Eastern')).strftime("%Y-%m-%d_%H:%M:%S_%Z")
for name in workflow_names:
    pathlib.Path(f'{name}_{now}_multiqc').mkdir(exist_ok=True)
    
print(f'Workflows found: {workflow_names}')

# remove duplicates
print(f'Removing duplicates/reruns...')
metrics_file_groups = defaultdict(list)
for wn, on, st in metrics_file_paths:
    metrics_file_groups[(wn, os.path.basename(on))].append((wn, on, st))

unique_metrics_files=[sorted(v,key=lambda t : datetime.strptime(t[2].split(".")[0],
                                       "%Y-%m-%dT%H:%M:%S"),
          reverse=True)[0] for v in metrics_file_groups.values()]
print(f'{len(metrics_file_paths)-len(unique_metrics_files)} duplicate/reruns removed, {len(unique_metrics_files)} unique metrics files remain.')

#download metrics files
print("Downloading metrics files...")
def get_bucket_and_blob(uri):
    return uri.replace("gs://", "").split("/", 1)

def download_file(workflow_name, uri):
    bucket_name, blob_name = get_bucket_and_blob(uri)
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    destination_name = f'{workflow_name}_{now}_multiqc/{blob_name.split("/")[-1]}'
    blob.download_to_filename(destination_name)
    
with ThreadPoolExecutor(max_workers=50) as executor:
    loop = asyncio.get_event_loop()
    tasks = [loop.run_in_executor(executor, download_file, workflow_name, uri) for 
             workflow_name, uri, _ in unique_metrics_files]
    await asyncio.gather(*tasks)
print(f'{len(unique_metrics_files)} metrics files downloaded')


workspace_bucket_name = os.environ.get("WORKSPACE_BUCKET").replace("gs://", "")
workspace_bucket = storage_client.bucket(workspace_bucket_name)
for name in workflow_names:
    directory = f'{name}_{now}_multiqc'
    report_name = f'{name}_{now}_multiqc.html'
    cmd =f'multiqc --filename {report_name} {directory}'
    !{cmd}
    print(f'Created report {report_name}')
    blob=workspace_bucket.blob(f'multiqc_reports/{report_name}')
    blob.upload_from_filename(report_name)
    print(f'Report has been uploaded to https://console.cloud.google.com/storage/browser/_details/{workspace_bucket_name}/multiqc_reports/{report_name}')
    
```
