import firecloud.api as fapi
import json
import argparse
import os
from datetime import datetime, timezone
import pytz
from dataclasses import dataclass, field


@dataclass
class AggregationSet:
    lab_batch: str
    group: int = 1
    delivered: bool = False
    contains_control: bool = False
    set_id: str = field(init=False)

    def __post_init__(self):
        if self.group > 1:
            self.set_id = f'{self.lab_batch}_group_{self.group}'
        elif self.group == 1:
            self.set_id = self.lab_batch
        else:
            raise RuntimeError(
                f'Group of aggregation set for lab_batch {self.lab_batch} is {self.group}, should not be less than 1')


def pre_existing_aggregation_set(lab_batch, group, delivered):
    return AggregationSet(lab_batch, group, delivered, True)


def next_aggregation_set(agg_set):
    new_group = agg_set.group + 1
    return AggregationSet(agg_set.lab_batch, new_group)


class GroupBuilder:

    def __init__(self, workspace_namespace, workspace_name):
        self.workspace_namespace = workspace_namespace
        self.workspace_name = workspace_name

        print('Finding tables to group by lab_batch')
        entity_types_response = fapi.list_entity_types(self.workspace_namespace, self.workspace_name)
        if not entity_types_response.ok:
            raise RuntimeError(f'ERROR: {entity_types_response.text}')
        self.entity_types_dict = json.loads(entity_types_response.text)
        self.available_tables = self.entity_types_dict.keys()

    def build_groups(self):
        for table_name, description in self.entity_types_dict.items():
            if all(x in description['attributeNames'] for x in ['is_control_sample', 'lab_batch']):
                self.group_samples_into_batches(table_name)

    def group_samples_into_batches(self, table_name):
        if os.path.exists('CreateSampleSets_data'):
            os.system('rm -r CreateSampleSets_data')
        os.system('mkdir -p CreateSampleSets_data')

        existing_sample_sets = dict()  # dictionary from sample_set id to list of samples in sample_set
        samples_already_in_aggregation_sets = set()  # set of samples already in aggregation sets
        lab_batch_sample_sets_dict = dict()  # dict from lab_batch to highest group aggregation_set for that lab_batch

        if f'{table_name}_set' in self.available_tables:
            # Download current sample_set table
            print(f'Downloading {table_name}_set table...')
            sample_set_response = fapi.get_entities(self.workspace_namespace, self.workspace_name, f'{table_name}_set')
            if not sample_set_response.ok:
                raise RuntimeError(f'ERROR: {sample_set_response.text}')
            sample_sets_dict = json.loads(sample_set_response.text)
            existing_sample_sets = {s['name']: [e['entityName'] for e in s['attributes'][f'{table_name}s']['items']] for
                                    s
                                    in sample_sets_dict}

            for sample_set in sample_sets_dict:
                samples = [e['entityName'] for e in sample_set['attributes'][f'{table_name}s']['items']]
                samples_already_in_aggregation_sets.update(samples)

            for sample_set in sample_sets_dict:
                attributes = sample_set['attributes']
                lab_batch = attributes['lab_batch']
                this_aggregation_set = pre_existing_aggregation_set(lab_batch, attributes['group'],
                                                                    attributes['delivered'])
                if lab_batch in lab_batch_sample_sets_dict:
                    if lab_batch_sample_sets_dict[lab_batch].group < this_aggregation_set.group:
                        if not lab_batch_sample_sets_dict[lab_batch].delivered:
                            raise RuntimeError(
                                f'Aggregation set {lab_batch_sample_sets_dict[lab_batch].set_id} has not been delivered, '
                                f'but later set {this_aggregation_set.set_id} also exists')
                        lab_batch_sample_sets_dict[attributes['lab_batch']] = this_aggregation_set
                else:
                    lab_batch_sample_sets_dict[attributes['lab_batch']] = this_aggregation_set

        # Read samples from samples table
        print(f'Reading {table_name} table...')
        sample_response = fapi.get_entities(self.workspace_namespace, self.workspace_name, f'{table_name}')
        if not sample_response.ok:
            raise RuntimeError(f'ERROR: {sample_response.text}')

        samples = json.loads(sample_response.text)
        # Writing new sample_set_membership.tsv
        added_sample_sets_dict = dict()  # dictionary from lab_batch to aggregation sets with added samples
        control_samples_dict = dict()  # dictionary from lab_batch to sample id of control sample
        added_samples_dict = dict()  # dictionary from set_id to list of samples to be added to the set
        with open(f'CreateSampleSets_data/new_{table_name}_set_membership.tsv', 'w') as new_membership_file:
            # Write header
            new_membership_file.write(f'membership:{table_name}_set_id\t{table_name}\n')
            for sample in samples:
                if 'lab_batch' not in sample['attributes']:
                    continue
                sample_name = sample['name']
                lab_batch = sample['attributes']['lab_batch']
                is_control_sample = sample['attributes']['is_control_sample']
                rework = sample['attributes'].get('rework', False)
                if is_control_sample:
                    # do we already have a control sample for this lab batch?  that would be bad...
                    if lab_batch in control_samples_dict:
                        raise RuntimeError(
                            f'Multiple control samples for lab_bath {lab_batch}: {sample_name}, {control_samples_dict[lab_batch]}')
                    # store control sample name in dictionary
                    control_samples_dict[lab_batch] = sample_name
                    # we do not create aggregation sets if we only have the control sample.  We will add control samples to newly created aggregation sets later
                    continue
                if rework or sample_name not in samples_already_in_aggregation_sets:
                    # this (non-control) sample needs to be added to an aggregation set.
                    # Find (or create) aggregation set to add it to.
                    if lab_batch not in added_sample_sets_dict:
                        if lab_batch in lab_batch_sample_sets_dict:
                            previous_aggregation_set = lab_batch_sample_sets_dict[lab_batch]
                            if previous_aggregation_set.delivered:
                                # we need to create the next aggregation set for this lab batch
                                added_sample_sets_dict[lab_batch] = next_aggregation_set(previous_aggregation_set)
                            else:
                                # we can add to the previous aggregation set
                                added_sample_sets_dict[lab_batch] = previous_aggregation_set
                        else:
                            # we need to create the first aggregation set for this lab batch
                            added_sample_sets_dict[lab_batch] = AggregationSet(lab_batch)
                    set_id = added_sample_sets_dict[lab_batch].set_id
                    if set_id in added_samples_dict:
                        added_samples_dict[set_id].append(sample_name)
                    else:
                        added_samples_dict[set_id] = [sample_name]
            # loop through added_sample_sets_dict and write set membership if we have control sample for set
            lab_batches_without_controls = list()
            for lab_batch, agg_set in added_sample_sets_dict.items():
                set_id = agg_set.set_id
                if agg_set.contains_control:
                    # if this aggregation set already contains control, we can simply add new samples
                    for sample in added_samples_dict[set_id]:
                        new_membership_file.write(f'{set_id}\t{sample}\n')
                elif lab_batch in control_samples_dict:
                    # found control sample, so this aggregation set can be added
                    # add control sample to this aggregation set
                    added_samples_dict[set_id].append(control_samples_dict[lab_batch])
                    # write samples, including controls
                    for sample in added_samples_dict[set_id]:
                        new_membership_file.write(f'{set_id}\t{sample}\n')
                else:
                    # no control sample for this aggregation set found, so will not aggregate yet
                    del added_samples_dict[set_id]
                    lab_batches_without_controls.append(lab_batch)
            for lab_batch in lab_batches_without_controls:
                del added_sample_sets_dict[lab_batch]
        if len(added_samples_dict) == 0:
            print(f'No new {table_name}_sets to be added.')
        else:
            if f'{table_name}_set' not in self.available_tables:
                print(f'Creating new table {table_name}_set')
                # Need to upload tsv to create new table
                with open(f'CreateSampleSets_data/new_{table_name}_set.tsv', 'w') as new_set_table:
                    new_set_table.write(f'entity:{table_name}_set_id\n')
                    for set_id in added_samples_dict:
                        new_set_table.write(f'{set_id}\n')
                upload_new_table_response = fapi.upload_entities_tsv(self.workspace_namespace, self.workspace_name,
                                                                     f'CreateSampleSets_data/new_{table_name}_set.tsv',
                                                                     "flexible")
                if not upload_new_table_response.ok:
                    raise RuntimeError(f'ERROR: {upload_new_table_response.text}')
            print(f'Uploading new {table_name}_set table... ')
            upload_response = fapi.upload_entities_tsv(self.workspace_namespace, self.workspace_name,
                                                       f'CreateSampleSets_data/new_{table_name}_set_membership.tsv',
                                                       "flexible")
            if not upload_response.ok:
                raise RuntimeError(f'ERROR: {upload_response.text}')
            # Add date and time created to sample_set
            print(f'Adding date and time to newly created {table_name}_sets...')

            now = str(datetime.now(pytz.timezone('US/Eastern')))
            for i, (this_lab_batch, this_aggregation_set) in enumerate(added_sample_sets_dict.items()):
                update_response = fapi.update_entity(self.workspace_namespace, self.workspace_name,
                                                     f'{table_name}_set', this_aggregation_set.set_id,
                                                     [{"op": "AddUpdateAttribute",
                                                       "attributeName": "time_sample_set_updated",
                                                       "addUpdateAttribute": now},
                                                      {"op": "AddUpdateAttribute", "attributeName": "delivered",
                                                       "addUpdateAttribute": False},
                                                      {"op": "AddUpdateAttribute", "attributeName": "redeliver",
                                                       "addUpdateAttribute": False},
                                                      {"op": "AddUpdateAttribute", "attributeName": "group",
                                                       "addUpdateAttribute": this_aggregation_set.group},
                                                      {"op": "AddUpdateAttribute", "attributeName": "lab_batch",
                                                       "addUpdateAttribute": this_aggregation_set.lab_batch}
                                                      ])
                if not update_response.ok:
                    raise RuntimeError(f'ERROR: {update_response.text}')
                print(f'    Completed {i + 1}/{len(added_samples_dict)}')

            print(f'Updating rework field in {table_name} table')
            for samples_list in added_samples_dict.values():
                for sample in samples_list:
                    update_response = fapi.update_entity(self.workspace_namespace, self.workspace_name,
                                                         f'{table_name}', sample,
                                                         [{"op": "AddUpdateAttribute", "attributeName": "rework",
                                                           "addUpdateAttribute": False}])
                    if not update_response.ok:
                        raise RuntimeError(f'ERROR: {update_response.text}')
            # Uploading new sample_set table
            print('SUCCESS')
            print(f'Printing update {table_name}_set_membership.tsv:')
            os.system(f'cat CreateSampleSets_data/new_{table_name}_set_membership.tsv')
        os.system('rm -r CreateSampleSets_data')


def main(workspace_namespace, workspace_name):
    group_builder = GroupBuilder(workspace_namespace, workspace_name)
    group_builder.build_groups()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace_namespace", dest="workspace_namespace", required=True)
    parser.add_argument("--worspace_name", dest="workspace_name", required=True)
    parser.parse_args()
    main(parser.workspace_namespace, parser.workspace_name)
