#!/usr/bin/env python3

import argparse
from io import StringIO
import firecloud.api as fapi
from google.cloud import storage
from datetime import datetime
import pytz

def get_workspace_bucket(namespace, name):
    """Get the bucket for a workspace"""
    response = fapi.get_workspace(namespace, name)
    if not response.ok:
        raise RuntimeError(f'ERROR getting workspace {namespace}/{name}: {response.text}')
    workspace = response.json()['workspace']
    return workspace['bucketName']

def copy_files_to_delivery_workspace(source_cohort, storage_client, target_bucket_name):
    """Copy all relevant files from the source cohort to the delivery workspace"""
    target_bucket = storage_client.bucket(target_bucket_name)
    target_base_path = f'GLIMPSEImputedCohorts/{source_cohort["name"]}/'
    delivered_files = {}
    """Create md5 file for imputed_vcf"""
    imputed_vcf_path = source_cohort['attributes'].get('imputed_vcf')
    if not imputed_vcf_path:
        raise RuntimeError("Error: imputed_vcf not found in source cohort attributes")
    imputed_vcf_blob = storage_client.bucket(imputed_vcf_path.split('/')[2]).get_blob(imputed_vcf_path.split('/',3)[3])
    imputed_vcf_md5 = imputed_vcf_blob.md5_hash
    md5_file_name = f'{imputed_vcf_path.split("/")[-1]}.md5'
    with open(md5_file_name, 'w') as md5_file:
        md5_file.write(imputed_vcf_md5)
    # Upload the md5 file
    target_blob_name = f'{target_base_path}{md5_file_name}'
    target_blob = target_bucket.blob(target_blob_name)
    target_blob.upload_from_filename(md5_file_name)
    delivered_files['imputed_vcf_md5'] = f'gs://{target_bucket_name}/{target_blob_name}'
    print(f"Copied imputed_vcf_md5 to {delivered_files['imputed_vcf_md5']}")
    
    outputs_to_deliver = ['qc_report', 'imputed_vcf', 'imputed_vcf_index', 'coverage_metrics', 'qc_metrics']
    
    for output in outputs_to_deliver:
        if output not in source_cohort['attributes']:
            raise RuntimeError(f"Error: {output} not found in source cohort attributes")
            
        source_gs_path = source_cohort['attributes'][output]
        source_bucket_name = source_gs_path.split('/')[2]
        source_bucket = storage_client.bucket(source_bucket_name)
        source_blob_name = source_gs_path.split('/', 3)[3]
        source_blob = storage_client.bucket(source_bucket_name).blob(source_blob_name)
        target_blob_name = f'{target_base_path}{source_blob_name.split("/")[-1]}'
        
        # Copy the file
        source_bucket.copy_blob(source_blob, target_bucket, target_blob_name)
        delivered_files[output] = f'gs://{target_bucket_name}/{target_blob_name}'
        print(f"Copied {output} to {delivered_files[output]}")
    
    return delivered_files

def create_or_update_glimpse_table(namespace, workspace, delivered_files, source_cohort):
    """Create or update the glimpse_imputed_cohorts table in the delivery workspace"""
    table_name = "glimpse_imputed_cohorts"
    
    with StringIO() as table_io:
        table_io.write(f'entity:{table_name}_id\t' + '\t'.join(delivered_files.keys()) + '\n')
        table_io.write(f'{source_cohort["name"]}\t' + '\t'.join(delivered_files.values()) + '\n')
        response = fapi.upload_entities_tsv(namespace, workspace, table_io, model='flexible')
        if not response.ok:
            raise RuntimeError(f'ERROR adding cohort {source_cohort['name']} to {table_name} table in {namespace}/{workspace}: {response.text}')
    
    print(f"Added cohort {source_cohort['name']} to {table_name} table in {namespace}/{workspace}")

def record_delivery(source_namespace, source_workspace, cohort_name, target_namespace, target_workspace, delivered_files):
    """Record the delivery in the source workspace's delivery table"""
    table_name = "delivery"
    now = datetime.now(pytz.UTC)
    delivery_id = f"{cohort_name}_to_{target_namespace}_{target_workspace}_{now.strftime('%Y%m%d_%H%M%SZ')}"
    
    with StringIO() as table_io:
        table_io.write('entity:delivery_id\tsource_cohort\ttarget_namespace\ttarget_workspace\tdelivery_date\t' + '\t'.join(delivered_files.keys()) + '\n')
        table_io.write(f'{delivery_id}\t{cohort_name}\t{target_namespace}\t{target_workspace}\t{now.strftime('%Y_%m_%d_%H:%M:%SZ')}\t' + '\t'.join(delivered_files.values()) + '\n')
        response = fapi.upload_entities_tsv(source_namespace, source_workspace, table_io, model='flexible')
        if not response.ok:
            raise RuntimeError(f'ERROR recording delivery f{delivery_id}: {response.text}')
    print(f"Recorded delivery {delivery_id} in {table_name} table in {source_namespace}/{source_workspace}")

def get_cohort(namespace, workspace, cohort_name):
    """Get a cohort entity from a workspace"""
    all_entities_response = fapi.get_entities_with_type(namespace, workspace)
    if not all_entities_response.ok:
        raise RuntimeError(f'ERROR getting entities from {namespace}/{workspace}: {all_entities_response.text}')
    matching_entities = [e for e in all_entities_response.json() if e['name'] == cohort_name]
    if not matching_entities:
        raise RuntimeError(f'ERROR: Cohort {cohort_name} not found in {namespace}/{workspace}')
    if len(matching_entities) > 1:
        raise RuntimeError(f'ERROR: Multiple entities named {cohort_name} found in {namespace}/{workspace}')
    return matching_entities[0]

def main(source_namespace, source_workspace, cohort_name, target_namespace, target_workspace):
    # Ask for confirmation before proceeding
    print(f"\nThis will deliver cohort '{cohort_name}' from workspace '{source_namespace}/{source_workspace}' "
            f"to workspace '{target_namespace}/{target_workspace}'")
    confirmation = input("Would you like to proceed? (y/n): ").lower().strip()
    if confirmation == 'y':
        print("\nStarting delivery...")
    
        # Initialize Google Cloud Storage client
        storage_client = storage.Client()
        
        # Get the source cohort
        source_cohort = get_cohort(source_namespace, source_workspace, cohort_name)
        
        # Get the target workspace bucket
        target_bucket_name = get_workspace_bucket(target_namespace, target_workspace)
        
        # Copy files to the delivery workspace
        delivered_files = copy_files_to_delivery_workspace(source_cohort, storage_client, target_bucket_name)
        
        # Create/update the glimpse_imputed_cohorts table in the delivery workspace
        create_or_update_glimpse_table(target_namespace, target_workspace, delivered_files, source_cohort)
        
        # Record the delivery in the source workspace
        record_delivery(source_namespace, source_workspace, cohort_name, target_namespace, target_workspace, delivered_files)
        
        print("Delivery completed successfully")
    else:
        print("Delivery cancelled.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Deliver a GLIMPSE imputed cohort to another workspace')
    parser.add_argument('--source-namespace', required=True,
                      help='Namespace of the source workspace')
    parser.add_argument('--source-workspace', required=True,
                      help='Name of the source workspace')
    parser.add_argument('--cohort', required=True,
                      help='Name of the cohort to deliver')
    parser.add_argument('--target-namespace', required=True,
                      help='Namespace of the target workspace')
    parser.add_argument('--target-workspace', required=True,
                      help='Name of the target workspace')
    
    args = parser.parse_args()
    
    main(args.source_namespace,
         args.source_workspace,
         args.cohort,
         args.target_namespace,
         args.target_workspace)
