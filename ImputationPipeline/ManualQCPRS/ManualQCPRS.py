import ipywidgets as widgets
from IPython.core.display_functions import display
from IPython.display import clear_output
import pandas as pd
import firecloud.api as fapi
import json
from google.cloud import storage
from datetime import datetime
import pytz
from abc import ABC, abstractmethod
import argparse


def get_bucket_and_blob(uri):
    return uri.replace("gs://", "").split("/", 1)


class color:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


class TableSelectionGUI:
    def __init__(self, workspace_namespace, workspace_name, workspace_bucket_name, available_tables):
        self.available_tables = available_tables
        self.workspace_namespace = workspace_namespace
        self.workspace_name = workspace_name
        self.workspace_bucket_name = workspace_bucket_name

    def run(self):
        self.build_table_selection_dropdown()

    def build_table_selection_dropdown(self):
        self.table_selection_dropdown = widgets.Dropdown(options=[None] + self.available_tables,
                                                         description='Select Table to Load Batches From',
                                                         style={"description_width": 'initial'},
                                                         layout=widgets.Layout(width='auto')
                                                         )
        self.table_selection_dropdown.observe(self.select_table, names="value")
        display(self.table_selection_dropdown)

    def select_table(self, change):
        if change['old'] != None:
            self.app.close_widgets()
        if change['new'] != None:
            self.app = ResultsModificationGUI(self.workspace_namespace, self.workspace_name, self.workspace_bucket_name,
                                              change['new'])
            self.app.run()


class WidgetGUI(ABC):
    def __init__(self):
        self.registered_widgets = []

    @abstractmethod
    def run(self):
        pass

    def close_widgets(self):
        for widget in self.registered_widgets:
            widget.close()

    def register_widget(self, widget):
        # test
        self.registered_widgets.append(widget)
        display(widget)


class ResultsModificationGUI(WidgetGUI):
    eastern_tz = pytz.timezone('US/Eastern')

    def __init__(self, workspace_namespace, workspace_name, workspace_bucket_name, table_name):
        super().__init__()
        self.workspace_namespace = workspace_namespace
        self.workspace_name = workspace_name
        self.workspace_bucket_name = workspace_bucket_name
        self.delivery_bucket = None
        self.table_name = table_name
        self.lab_batch_output_box = widgets.Output(layout={'border': '1px solid black'})
        self.status_output_box = widgets.Output(layout={'border': '1px solid black'})
        self.selected_batch_gui = None
        self.register_widget(self.status_output_box)
        self.register_widget(self.lab_batch_output_box)
        self.initialize_workspace_info_and_gcloud()

    def run(self):
        self.build_lab_batch_selection_section()

    def close_widgets(self):
        super().close_widgets()

        if self.selected_batch_gui != None:
            self.selected_batch_gui.close_widgets()

    def initialize_workspace_info_and_gcloud(self):
        with self.status_output_box:
            print(
                f"Getting {self.table_name} information for workspace " + self.workspace_namespace + "/" + self.workspace_name + "...")

        sample_set_response = fapi.get_entities(self.workspace_namespace, self.workspace_name, self.table_name)
        if not sample_set_response.ok:
            raise RuntimeError(f'ERROR: {sample_set_response.text}')

        with self.status_output_box:
            print("Converting json to dictionary...")

        sample_sets = json.loads(sample_set_response.text)

        with self.status_output_box:
            print("Finding batches with aggregated results...")

        # lab_batch selection and import/export buttons
        self.agg_batch_map = {s['name']: s for s in sample_sets if 'batch_all_results' in s['attributes'] and
                              (not s['attributes']['delivered'] or s['attributes']["redeliver"])}

        with self.status_output_box:
            print("Done")

        self.lookup_delivery_bucket()

        self.storage_client = storage.Client()

    def lookup_delivery_bucket(self):
        with self.status_output_box:
            print("Looking up delivery-bucket")
        workspace_response = fapi.get_workspace(self.workspace_namespace, self.workspace_name)
        if not workspace_response.ok:
            raise RuntimeError(f'ERROR: {workspace_response.text}')

        workspace_attributes = workspace_response.json()['workspace']['attributes']
        if 'delivery-bucket' in workspace_attributes:
            self.delivery_bucket = workspace_attributes['delivery-bucket']
            with self.status_output_box:
                print(f"Delivery bucket set to {color.BOLD} {color.BLUE} {self.delivery_bucket} {color.END}")
        else:
            with self.status_output_box:
                print("No delivery bucket found")

    def build_lab_batch_selection_section(self):
        self.agg_batch_selection_dropdown = widgets.Combobox(options=list(self.agg_batch_map.keys()),
                                                             description='Select Lab Batch',
                                                             style={"description_width": 'initial'},
                                                             layout=widgets.Layout(width='auto'),
                                                             continuous_update=False,
                                                             ensure_option=True
                                                             )
        self.agg_batch_download_button = widgets.Button(description="Load Initial Batch Results",
                                                        disabled=True,
                                                        layout=widgets.Layout(width='auto')
                                                        )
        self.agg_batch_upload_button = widgets.Button(description="Save QC'd Batch Results",
                                                      disabled=True,
                                                      layout=widgets.Layout(width='auto')
                                                      )

        self.discard_button = widgets.Button(description="Discard Loaded Results",
                                             disabled=True,
                                             layout=widgets.Layout(width='auto')
                                             )

        self.lab_batch_hbox = widgets.HBox([self.agg_batch_selection_dropdown,
                                            self.agg_batch_download_button,
                                            self.agg_batch_upload_button,
                                            self.discard_button
                                            ])
        self.register_widget(self.lab_batch_hbox)

        self.agg_batch_selection_dropdown.observe(self.lab_batch_selected, names='value')
        self.agg_batch_download_button.on_click(self.lab_batch_download_button_clicked)
        self.agg_batch_upload_button.on_click(self.save_modified_results_button_clicked)
        self.discard_button.on_click(self.discard_button_clicked)

    def lab_batch_selected(self, change):
        lab_batch_name = change['new']
        if lab_batch_name == '':
            return
        self.print_lab_batch_result_info(lab_batch_name)
        self.agg_batch_download_button.disabled = False

    def get_time_created(self, uri):
        bucket_name, path = get_bucket_and_blob(uri)
        bucket = self.storage_client.get_bucket(bucket_name)
        blob = bucket.get_blob(path)
        return blob.time_created

    def print_lab_batch_result_info(self, lab_batch_name):
        with self.lab_batch_output_box:
            clear_output()
            results_uri = self.agg_batch_map[lab_batch_name]['attributes']['batch_all_results']
            time_results_created = self.get_time_created(results_uri)
            print("Aggregated results for " + lab_batch_name + " were created at " +
                  time_results_created.astimezone(tz=self.eastern_tz).strftime("%Y-%m-%d %H:%M:%S %Z"))
            if 'qcd_batch_results' in self.agg_batch_map[lab_batch_name]['attributes']:
                prev_qcd_results_uri = self.agg_batch_map[lab_batch_name]['attributes']['qcd_batch_results']
                time_previous_qcd_results_created = self.get_time_created(prev_qcd_results_uri)
                print("This batch already has a qcd results file, which was created at " +
                      time_previous_qcd_results_created.astimezone(tz=self.eastern_tz).strftime("%Y-%m-%d %H:%M:%S %Z"))

    def lab_batch_download_button_clicked(self, button):

        self.agg_batch_download_button.disabled = True
        self.agg_batch_selection_dropdown.disabled = True
        self.discard_button.disabled = False
        self.agg_batch_upload_button.disabled = False
        self.agg_batch = self.agg_batch_selection_dropdown.value
        results_uri = self.agg_batch_map[self.agg_batch]['attributes']['batch_all_results']
        results = self.download_aggregated_results(results_uri)
        lab_batch = self.agg_batch_map[self.agg_batch]['attributes']['lab_batch']
        self.selected_batch_gui = SelectedBatchModificationGui(results, lab_batch)
        self.selected_batch_gui.run()

    def download_aggregated_results(self, uri):
        return pd.read_csv(uri, delimiter='\t').fillna("NA").set_index('sample_id')

    def discard_button_clicked(self, button):
        self.agg_batch_upload_button.disabled = True
        self.agg_batch_download_button.disabled = True
        self.discard_button.disabled = True
        self.selected_batch_gui.close_widgets()
        self.agg_batch_selection_dropdown.disabled = False
        self.agg_batch_selection_dropdown.value = ''
        self.agg_batch = None
        self.selected_batch_gui = None

    def save_modified_results_button_clicked(self, button):
        self.agg_batch_upload_button.disabled = True
        self.agg_batch_download_button.disabled = True
        self.discard_button.disabled = True
        self.selected_batch_gui.close_widgets()
        path = self.save_modified_batch_results()
        self.update_sample_sets_table_with_modified_results(path)
        self.agg_batch_selection_dropdown.disabled = False
        self.agg_batch_selection_dropdown.value = ''
        self.selected_batch_gui = None
        if self.delivery_bucket:
            self.delivery_confirmation_box = DeliveryConfirmationBox(path, self)
            self.delivery_confirmation_box.run()

        self.agg_batch = None

    def update_sample_sets_table_with_modified_results(self, path):
        with self.status_output_box:
            print(f"Updating {self.table_name} data table for " + self.agg_batch + "...")
        update_response = fapi.update_entity(self.workspace_namespace, self.workspace_name, self.table_name,
                                             f'{self.agg_batch}',
                                             [{"op": "AddUpdateAttribute",
                                               "attributeName": "qcd_batch_results", "addUpdateAttribute": path}]
                                             )
        if not update_response.ok:
            raise RuntimeError(f'ERROR: {update_response.text}')
        with self.status_output_box:
            print("Done")

    def save_modified_batch_results(self):
        modified_batch_results = self.selected_batch_gui.modified_results.append(
            self.selected_batch_gui.failed_imputation_samples).fillna("NA")
        now = datetime.now(self.eastern_tz).strftime("%Y-%m-%d_%H:%M:%S_%Z")
        path = "/".join([self.workspace_bucket_name,
                         "manually_qcd_prs_results",
                         "_".join([
                             self.agg_batch,
                             now,
                             "manually_qcd_prs_results.csv"
                         ])
                         ])
        with self.status_output_box:
            print("Saving qcd results to " + path)
        modified_batch_results.to_csv(path)
        return path


class DeliveryConfirmationBox(WidgetGUI):

    def __init__(self, source_path, requesting_results_modification_gui: ResultsModificationGUI):
        super().__init__()
        self.source_path = source_path
        self.delivery_bucket = requesting_results_modification_gui.delivery_bucket
        self.agg_batch = requesting_results_modification_gui.agg_batch
        self.workspace_namespace = requesting_results_modification_gui.workspace_namespace
        self.workspace_name = requesting_results_modification_gui.workspace_name
        self.table_name = requesting_results_modification_gui.table_name
        self.requesting_results_modification_gui = requesting_results_modification_gui
        self.requesting_results_modification_gui.agg_batch_selection_dropdown.disabled = True
        self.delivery_output_box = widgets.Output(layout={'border': '1px solid black'})
        self.register_widget(self.delivery_output_box)

    def run(self):
        self.build_confirmation_box()

    def build_confirmation_box(self):
        with self.delivery_output_box:
            print(
                f'{color.BOLD} Are you sure you want to deliver {self.source_path} to bucket {color.BLUE} {self.delivery_bucket}? {color.END}')
        self.deliver_button = widgets.Button(description="Yes, deliver",
                                             layout=widgets.Layout(width='auto'),
                                             button_style='success'
                                             )

        self.dont_deliver_button = widgets.Button(description="No, do not deliver",
                                                  layout=widgets.Layout(width='auto'),
                                                  button_style='danger'
                                                  )

        self.button_hbox = widgets.HBox([self.deliver_button,
                                         self.dont_deliver_button
                                         ])
        self.deliver_button.on_click(self.deliver_button_clicked)
        self.dont_deliver_button.on_click(self.dont_deliver_button_clicked)

        self.register_widget(self.button_hbox)

    def deliver_button_clicked(self, button):
        self.deliver_modified_batch_results()
        self.close_widgets()
        self.update_data_table()
        self.remove_set_from_lab_selection_dropdown()
        self.requesting_results_modification_gui.agg_batch_selection_dropdown.disabled = False

    def dont_deliver_button_clicked(self, button):
        self.close_widgets()
        self.requesting_results_modification_gui.agg_batch_selection_dropdown.disabled = False

    def deliver_modified_batch_results(self):
        source_bucket_name, source_blob_name = get_bucket_and_blob(self.source_path)
        storage_client = storage.Client()
        source_bucket = storage_client.bucket(source_bucket_name)
        source_blob = source_bucket.blob(source_blob_name)

        destination_bucket = storage_client.bucket(self.delivery_bucket.replace("gs://", ""))
        destination_blob_name = source_blob_name.split("/")[-1]
        source_bucket.copy_blob(source_blob, destination_bucket, destination_blob_name)
        with self.requesting_results_modification_gui.status_output_box:
            print(f'{color.BOLD} {color.BLUE} Delivering qcd results to bucket {destination_bucket.name} {color.END}')

    def update_data_table(self):
        update_response = fapi.update_entity(self.workspace_namespace, self.workspace_name,self.table_name, self.agg_batch,
                           [fapi._attr_set("delivered", True),
                            fapi._attr_set("redeliver", False)
                            ]
                           )
        fapi._check_response_code(update_response, 200)

    def remove_set_from_lab_selection_dropdown(self):
        self.requesting_results_modification_gui.agg_batch_selection_dropdown.options = \
            tuple(option for option in
                  self.requesting_results_modification_gui.agg_batch_selection_dropdown.options if
                  option != self.agg_batch)

class SelectedBatchModificationGui(WidgetGUI):
    result_to_status_dict = {"HIGH": "PASS",
                             "NOT_HIGH": "PASS",
                             "NOT_RESULTED": "FAIL",
                             "NA": "INFO"}

    status_to_style_dict = {"PASS": "success",
                            "FAIL": "danger",
                            "INFO": "info"}

    def __init__(self, results, lab_batch):
        super().__init__()
        self.manual_fails = {}
        self.drop_down_label_dict = {}
        self.drop_down_custom_dict = {}
        self.custom_failure_label_dict = {}
        self.failed_conditions = {}
        self.manual_failure_hbox_dict = {}
        self.results = results
        self.lab_batch = lab_batch
        self.modified_results = results.copy()
        self.modified_results = self.modified_results.reindex(columns=self.modified_results.columns.tolist() +
                                                                      ['poly_test_not_performed_reason', 'notes'])
        self.validate_results()

    def validate_results(self):
        batches_in_results = set(self.results['lab_batch'])
        if len(batches_in_results) > 1:
            raise RuntimeError(f'Error: more than one lab_batch in the same results table.  Batches in table are: {batches_in_results}')
        if self.lab_batch not in batches_in_results:
            raise RuntimeError(
                f'Error: lab_batch {self.lab_batch } not seen in results table, which contains batch {batches_in_results}')

    def run(self):
        self.build_failed_imputation_section()

    def build_failed_imputation_section(self):
        self.out_failed_imputation = widgets.Output(layout={'border': '1px solid black'})

        self.conditions = [col.rsplit("_", 1)[0] for col in self.results.columns if col.endswith("risk")]

        self.failed_imputation_samples = pd.DataFrame(
            columns=['sample_id', 'lab_batch', 'is_control_sample', 'poly_test_not_performed_reason',
                     'notes']).set_index('sample_id')
        self.printFailedImputationSamples()

        self.failed_imputation_sample_text_box = widgets.Text(continuous_update=False,
                                                              placeholder='Add Sample Which Failed Imputation...')
        self.failed_imputation_notes_text_box = widgets.Text(disabled=True)
        self.failed_imputation_sample_add_button = widgets.Button(description='Add/Update Failed Sample', disabled=True,
                                                                  layout=widgets.Layout(width='auto'))
        self.failed_imputation_sample_remove_button = widgets.Button(description='Remove Failed Sample', disabled=True,
                                                                     layout=widgets.Layout(width='auto')
                                                                     )
        self.finished_addition_imputation_failures_button = widgets.Button(
            description='Finished Adding Imputation Failures',
            layout=widgets.Layout(width='auto')
        )

        self.failed_imputation_sample_text_box.observe(self.enterTextForSampleFailedImputation, names='value')
        self.failed_imputation_sample_add_button.on_click(self.clickAddFailedImputationSampleButton)
        self.failed_imputation_sample_remove_button.on_click(self.clickRemoveFailedImputationSampleButton)
        self.finished_addition_imputation_failures_button.on_click(self.clickFinishedImputationFailuresButton)

        self.failed_imputation_hbox = widgets.HBox([self.failed_imputation_sample_text_box,
                                                    self.failed_imputation_notes_text_box,
                                                    self.failed_imputation_sample_add_button,
                                                    self.failed_imputation_sample_remove_button
                                                    ],
                                                   layout=widgets.Layout(display='inline-flex', flex_flow='row wrap')
                                                   )

        self.failed_imputation_vbox = widgets.VBox([self.out_failed_imputation,
                                                    self.failed_imputation_hbox,
                                                    self.finished_addition_imputation_failures_button
                                                    ])
        self.register_widget(self.failed_imputation_vbox)

    def printFailedImputationSamples(self):
        with self.out_failed_imputation:
            clear_output()
            if self.failed_imputation_samples.empty:
                print('No samples will be added as failing imputation')
            else:
                print('The following samples will be added as failing imputation')
                for row in self.failed_imputation_samples.iterrows():
                    print(row[0] + ", notes: " + row[1]['notes'])

    def addFailedImputationSample(self, sample, notes):
        self.failed_imputation_samples.loc[sample, 'poly_test_not_performed_reason'] = 'Failed Imputation'
        self.failed_imputation_samples.loc[sample, 'notes'] = notes
        self.failed_imputation_samples.loc[sample, 'lab_batch'] = self.lab_batch
        self.failed_imputation_samples.loc[sample, 'is_control_sample'] = False
        self.printFailedImputationSamples()

    def removeFailedImputationSample(self, sample):
        self.failed_imputation_samples.drop(sample, inplace=True)
        self.printFailedImputationSamples()

    def enterTextForSampleFailedImputation(self, change):
        if change['new'] == '':
            return
        if change['new'] in self.results.index:
            with self.out_failed_imputation:
                print(f'{change["new"]} already included as scored sample, '
                      f'cannot be added as sample which failed imputation')
            self.failed_imputation_sample_text_box.value = ''
            return
        self.failed_imputation_notes_text_box.disabled = False
        self.failed_imputation_notes_text_box.placeholder = 'notes...'
        if change['new'] in self.failed_imputation_samples.index:
            self.failed_imputation_sample_remove_button.disabled = False
            self.failed_imputation_notes_text_box.value = self.failed_imputation_samples.loc[change['new'], 'notes']
        self.failed_imputation_sample_add_button.disabled = False

    def resetFailedImputationWidgets(self):
        self.failed_imputation_notes_text_box.disabled = True
        self.failed_imputation_sample_add_button.disabled = True
        self.failed_imputation_sample_remove_button.disabled = True
        self.failed_imputation_sample_text_box.value = ''
        self.failed_imputation_notes_text_box.value = ''

    def clickAddFailedImputationSampleButton(self, button):
        sample = self.failed_imputation_sample_text_box.value
        notes = self.failed_imputation_notes_text_box.value
        self.addFailedImputationSample(sample, notes)
        self.resetFailedImputationWidgets()

    def clickRemoveFailedImputationSampleButton(self, button):
        sample = self.failed_imputation_sample_text_box.value
        self.removeFailedImputationSample(sample)
        self.resetFailedImputationWidgets()

    def clickFinishedImputationFailuresButton(self, button):
        self.failed_imputation_hbox.close()
        self.finished_addition_imputation_failures_button.close()
        self.build_fail_samples_all_conditions_section()
        self.printFailedSamplesAllConditions()

    def build_fail_samples_all_conditions_section(self):
        self.out_failed_sample_all_conditions = widgets.Output(layout={'border': '1px solid black'})
        self.non_control_samples = list(self.results.loc[~self.results['is_control_sample']].index)
        self.sample_failure_selection = widgets.Combobox(options=self.non_control_samples,
                                                         description='Select Sample To Fail PRS for All Conditions',
                                                         style={"description_width": 'initial'},
                                                         layout=widgets.Layout(width='auto'),
                                                         continuous_update=False,
                                                         ensure_option=True
                                                         )
        self.sample_failure_reason_selection = widgets.Combobox(options=['PCA outlier'],
                                                                continuous_update=False,
                                                                ensure_option=False,
                                                                disabled=True,
                                                                layout=widgets.Layout(width='auto')
                                                                )
        self.fail_sample_all_conditions_button = widgets.Button(description='Fail', disabled=True,
                                                                style={"description_width": 'initial'},
                                                                layout=widgets.Layout(width='auto')
                                                                )
        self.unfail_sample_all_conditions_button = widgets.Button(description='Unfail', disabled=True,
                                                                  style={"description_width": 'initial'},
                                                                  layout=widgets.Layout(width='auto')
                                                                  )

        self.finished_sample_failures_button = widgets.Button(description='Finished Failing Samples For All Conditions',
                                                              layout=widgets.Layout(width='auto')
                                                              )

        self.sample_failure_hbox = widgets.HBox([self.sample_failure_selection, self.sample_failure_reason_selection,
                                                 self.fail_sample_all_conditions_button,
                                                 self.unfail_sample_all_conditions_button
                                                 ])
        self.sample_failure_vbox = widgets.VBox([self.out_failed_sample_all_conditions,
                                                 self.sample_failure_hbox,
                                                 self.finished_sample_failures_button
                                                 ])

        self.sample_failure_selection.observe(self.selectSampleToFailAllConditions, names='value')
        self.sample_failure_reason_selection.observe(self.selectSampleFailureReason, names='value')

        self.fail_sample_all_conditions_button.on_click(self.clickFailForAllConditionsButton)
        self.unfail_sample_all_conditions_button.on_click(self.clickUnfailForAllConditionsButton)

        self.finished_sample_failures_button.on_click(self.clickFinishedSampleFailuresButton)
        self.register_widget(self.sample_failure_vbox)

    # when sample selected, enable reason combobox, and unfail button if already failed
    def selectSampleToFailAllConditions(self, change):
        sample_to_fail = change['new']
        if sample_to_fail != '':
            self.sample_failure_reason_selection.disabled = False
            if self.modified_results.loc[sample_to_fail, 'poly_test_not_performed_reason'] == 'Failed PRS':
                self.sample_failure_reason_selection.value = self.modified_results.loc[sample_to_fail, 'notes']
                self.unfail_sample_all_conditions_button.disabled = False
            else:
                self.sample_failure_reason_selection.value = ''
        else:
            self.resetSampleFailureWidgets()

    # when reason selected, enable Fail button
    def selectSampleFailureReason(self, changed):
        self.fail_sample_all_conditions_button.disabled = False

    def addManualFailureToResults(self, sample, condition, reason):
        # only need to fail if scored
        if self.results.loc[sample, f'{condition}_risk'] != "NA":
            self.modified_results.loc[sample, condition + "_raw"] = "NA"
            self.modified_results.loc[sample, condition + "_adjusted"] = "NA"
            self.modified_results.loc[sample, condition + "_percentile"] = "NA"
            self.modified_results.loc[sample, condition + "_risk"] = "NOT_RESULTED"
            self.modified_results.loc[sample, condition + "_reason_not_resulted"] = reason

    def removeManualFailureFromResults(self, sample, condition):
        self.modified_results.loc[sample, condition + "_raw"] = self.results.loc[sample, condition + "_raw"]
        self.modified_results.loc[sample, condition + "_adjusted"] = self.results.loc[sample, condition + "_adjusted"]
        self.modified_results.loc[sample, condition + "_percentile"] = self.results.loc[
            sample, condition + "_percentile"]
        self.modified_results.loc[sample, condition + "_risk"] = self.results.loc[sample, condition + "_risk"]
        self.modified_results.loc[sample, condition + "_reason_not_resulted"] = self.results.loc[
            sample, condition + "_reason_not_resulted"]

    def manuallyFailSampleForAllConditions(self, sample, reason):
        for condition in self.conditions:
            self.addManualFailureToResults(sample, condition, reason)
        self.modified_results.loc[sample, 'poly_test_not_performed_reason'] = 'Failed PRS'
        self.modified_results.loc[sample, 'notes'] = reason
        self.printFailedSamplesAllConditions()

    def manuallyUnFailSampleForAllConditions(self, sample):
        for condition in self.conditions:
            self.removeManualFailureFromResults(sample, condition)
        self.modified_results.loc[sample, 'poly_test_not_performed_reason'] = 'NA'
        self.modified_results.loc[sample, 'notes'] = 'NA'
        self.printFailedSamplesAllConditions()

    def manuallyFailConditionForSample(self, sample, condition, reason):
        self.addManualFailureToResults(sample, condition, reason)
        if sample not in self.manual_fails:
            self.manual_fails[sample] = {condition: reason}
        else:
            self.manual_fails[sample][condition] = reason

    def manuallyUnFailConditionForSample(self, sample, condition):
        self.removeManualFailureFromResults(sample, condition)
        self.manual_fails[sample].pop(condition)
        if len(self.manual_fails[sample]) == 0:
            self.manual_fails.pop(sample)

    # when buttons clicked, fail or unfail sample for all conditions
    def resetSampleFailureWidgets(self):
        self.sample_failure_selection.value = ''
        self.sample_failure_reason_selection.value = ''
        self.sample_failure_reason_selection.disabled = True
        self.fail_sample_all_conditions_button.disabled = True
        self.unfail_sample_all_conditions_button.disabled = True

    def clickFailForAllConditionsButton(self, b):
        sample = self.sample_failure_selection.value
        reason = self.sample_failure_reason_selection.value
        self.manuallyFailSampleForAllConditions(sample, reason)
        self.resetSampleFailureWidgets()

    def clickUnfailForAllConditionsButton(self, b):
        sample = self.sample_failure_selection.value
        self.manuallyUnFailSampleForAllConditions(sample)
        self.resetSampleFailureWidgets()

    def printFailedSamplesAllConditions(self):
        with self.out_failed_sample_all_conditions:
            clear_output()
            failed_samples_all_conditions_results = self.modified_results.loc[
                self.modified_results['poly_test_not_performed_reason'] == 'Failed PRS']
            if failed_samples_all_conditions_results.empty:
                print('No Samples Will Be Failed For All Conditions')
            else:
                print('The following samples will be marked as Failed PRS and failed for all conditions')
                for row in failed_samples_all_conditions_results.iterrows():
                    print(row[0] + ", notes: " + row[1]['notes'])

    # third section, fail condition for all samples
    # third section displays when second section completed
    def clickFinishedSampleFailuresButton(self, button):
        self.sample_failure_hbox.close()
        self.finished_sample_failures_button.close()
        self.build_fail_condition_section()
        self.printFailedConditionsAllSamples()

    def build_fail_condition_section(self):
        self.out_failed_condition_all_samples = widgets.Output(layout={'border': '1px solid black'})
        self.condition_failure_selection = widgets.Dropdown(options=[None] + self.conditions,
                                                            description='Select Condition to Fail for All Samples',
                                                            style={"description_width": 'initial'},
                                                            layout=widgets.Layout(width='auto')
                                                            )
        self.condition_failure_reason_selection = widgets.Combobox(
            options=['control sample score outside expected range'],
            continuous_update=False,
            ensure_option=False,
            disabled=True,
            layout=widgets.Layout(width='auto')
        )
        self.fail_condition_all_samples_button = widgets.Button(description='Fail', disabled=True,
                                                                style={"description_width": 'initial'},
                                                                layout=widgets.Layout(width='auto')
                                                                )
        self.unfail_condition_all_samples_button = widgets.Button(description='Unfail', disabled=True,
                                                                  style={"description_width": 'initial'},
                                                                  layout=widgets.Layout(width='auto')
                                                                  )

        self.finished_condition_failures_button = widgets.Button(
            description='Finished Failing Conditions For All Samples',
            layout=widgets.Layout(width='auto')
        )

        self.condition_failure_hbox = widgets.HBox(
            [self.condition_failure_selection, self.condition_failure_reason_selection,
             self.fail_condition_all_samples_button, self.unfail_condition_all_samples_button
             ])
        self.condition_failure_vbox = widgets.VBox([self.out_failed_condition_all_samples,
                                                    self.condition_failure_hbox,
                                                    self.finished_condition_failures_button
                                                    ])
        self.condition_failure_selection.observe(self.selectConditionToFailAllConditions, names='value')
        self.condition_failure_reason_selection.observe(self.selectConditionFailureReason, names='value')
        self.fail_condition_all_samples_button.on_click(self.clickFailForAllSamplesButton)
        self.unfail_condition_all_samples_button.on_click(self.clickUnfailForAllSamplesButton)
        self.finished_condition_failures_button.on_click(self.clickFinishedConditionFailuresButton)
        self.register_widget(self.condition_failure_vbox)

    # when condition selected, enable reason combobox
    def selectConditionToFailAllConditions(self, change):
        condition_to_fail = change['new']
        if condition_to_fail is not None:
            self.condition_failure_reason_selection.disabled = False
            if condition_to_fail in self.failed_conditions:
                self.unfail_condition_all_samples_button.disabled = False
                self.condition_failure_reason_selection.value = self.failed_conditions[condition_to_fail]
            else:
                self.condition_failure_reason_selection.value = ''
        else:
            self.resetConditionFailureWidgets()

    # when reason selected, enable Fail button
    def selectConditionFailureReason(self, changed):
        self.fail_condition_all_samples_button.disabled = False

    # when buttons clicked, fail or unfail condition for all samples
    def resetConditionFailureWidgets(self):
        self.condition_failure_selection.value = None
        self.condition_failure_reason_selection.value = ''
        self.condition_failure_reason_selection.disabled = True
        self.fail_condition_all_samples_button.disabled = True
        self.unfail_condition_all_samples_button.disabled = True

    def clickFailForAllSamplesButton(self, b):
        condition = self.condition_failure_selection.value
        reason = self.condition_failure_reason_selection.value
        self.manuallyFailConditionForAllSamples(condition, reason)
        self.resetConditionFailureWidgets()

    def clickUnfailForAllSamplesButton(self, b):
        condition = self.condition_failure_selection.value
        self.manuallyUnFailConditionForAllSamples(condition)
        self.resetConditionFailureWidgets()

    def manuallyFailConditionForAllSamples(self, condition, reason):
        for sample in self.non_control_samples:
            self.addManualFailureToResults(sample, condition, reason)
        self.failed_conditions[condition] = reason
        self.printFailedConditionsAllSamples()

    def manuallyUnFailConditionForAllSamples(self, condition):
        for sample in self.non_control_samples:
            self.removeManualFailureFromResults(sample, condition)
        self.failed_conditions.pop(condition)
        self.printFailedConditionsAllSamples()

    def printFailedConditionsAllSamples(self):
        with self.out_failed_condition_all_samples:
            clear_output()
            if not self.failed_conditions:
                print('No Conditions Will Be Failed For All Samples')
            else:
                print('The following conditions will be failed for all samples')
                for condition, reason in self.failed_conditions.items():
                    print(condition + ", reason: " + reason)

    # final section, fail particular sample for particular condition
    # display when finished failing conditions for all samples
    def clickFinishedConditionFailuresButton(self, button):
        self.condition_failure_hbox.close()
        self.finished_condition_failures_button.close()
        self.build_fail_sample_for_condition_section()
        self.printManualFailureStatus()

    def build_fail_sample_for_condition_section(self):
        self.out = widgets.Output(layout={'border': '1px solid black', 'text_color': 'gray'})
        self.sample_selection = widgets.Combobox(options=self.non_control_samples,
                                                 description='Select Sample To Fail Certain Conditions',
                                                 style={"description_width": 'initial'},
                                                 layout=widgets.Layout(width='auto'),
                                                 ensure_options=True,
                                                 continuous_update=False)

        observed_failure_reasons = {reason for sublist in
                                    [self.modified_results[c + "_reason_not_resulted"].unique() for c in
                                     self.conditions] for reason in sublist if reason != 'NA'}
        self.standard_failure_reasons = list({"Too Many Scoring Sites Missing"}.union(observed_failure_reasons))
        self.standard_failure_options = [('PASS', 'PASS')] + [("FAIL: " + reason, reason) for reason in
                                                              self.standard_failure_reasons] + [
                                            ('FAIL: Other Reason', 'Other')]

        self.sample_condition_failure_vbox = widgets.VBox([self.out, self.sample_selection])
        self.sample_selection.observe(self.select_sample, names='value')
        self.register_widget(self.sample_condition_failure_vbox)

    def select_failure_reason(self, change):
        condition_label = self.drop_down_label_dict[change['owner']]
        condition = condition_label.description
        new_value = change['new']
        sample = self.sample_selection.value
        custom_reason_text = self.drop_down_custom_dict[change['owner']]
        if new_value == 'PASS':
            self.manuallyUnFailConditionForSample(sample, condition)
            condition_label.button_style = 'success'
            self.resetCustomReasonText(custom_reason_text)
            self.printManualFailureStatus()
        elif new_value == 'Other':
            custom_reason_text.disabled = False
            custom_reason_text.placeholder = 'Enter Custom Failure Reason'
        else:
            self.manuallyFailConditionForSample(sample, condition, new_value)
            condition_label.button_style = 'danger'
            self.resetCustomReasonText(custom_reason_text)
            self.printManualFailureStatus()

    def resetCustomReasonText(self, custom_reason_text):
        custom_reason_text.disabled = True
        custom_reason_text.placeholder = ''
        custom_reason_text.value = ''

    def select_custom_failure_reason(self, change):
        if change['new'] == '':
            return
        condition_label = self.custom_failure_label_dict[change['owner']]
        condition = condition_label.description
        new_value = change['new']
        sample = self.sample_selection.value
        condition_label.button_style = 'danger'
        self.manuallyFailConditionForSample(sample, condition, new_value)
        self.printManualFailureStatus()

    def printManualFailureStatus(self):
        with self.out:
            clear_output()
            if not self.manual_fails:
                print("No additional manual fails")
            else:
                print("The following manual fails will be added:")
                for s in self.manual_fails:
                    print(s + " : " + ", ".join([c + "(" + r + ")" for c, r in self.manual_fails[s].items()]))

    def build_sample_condition_VBox(self, condition, sample):
        condition_status = self.result_to_status_dict[self.modified_results.loc[sample, condition + "_risk"]]
        condition_label = widgets.Button(description=condition,
                                         button_style=self.status_to_style_dict[condition_status])

        failure_reason = widgets.Dropdown(options=["NOT CALCULATED"] if condition_status == "INFO" else
            self.standard_failure_options, disabled=condition_status != "PASS")
        if condition_status == "PASS":
            failure_reason.value = 'PASS'
        elif condition_status != "INFO":
            failure_reason.value = self.modified_results.loc[sample, condition + "_reason_not_resulted"] if \
                self.modified_results.loc[
                    sample, condition + "_reason_not_resulted"] in self.standard_failure_reasons else 'Other'
        failure_reason.observe(self.select_failure_reason, names='value')

        custom_failure_reason = widgets.Text(disabled=True, continuous_update=False)
        if failure_reason.value == 'Other':
            custom_failure_reason.value = self.modified_results.loc[sample, condition + "_reason_not_resulted"]
        custom_failure_reason.observe(self.select_custom_failure_reason, names='value')
        self.drop_down_label_dict[failure_reason] = condition_label
        self.drop_down_custom_dict[failure_reason] = custom_failure_reason
        self.custom_failure_label_dict[custom_failure_reason] = condition_label
        box = widgets.VBox([condition_label, failure_reason, custom_failure_reason])
        return box

    def select_sample(self, change):
        old_sample = change['old']
        new_sample = change['new']
        if old_sample != '':
            self.manual_failure_hbox_dict[old_sample].layout.display = 'none'
        if new_sample == '':
            return
        if new_sample not in self.manual_failure_hbox_dict:
            vboxes = [self.build_sample_condition_VBox(c, change['new']) for c in self.conditions]
            hbox = widgets.HBox(vboxes, layout=widgets.Layout(display='inline-flex', flex_flow='row wrap'))
            self.manual_failure_hbox_dict[new_sample] = hbox
            self.register_widget(hbox)
        else:
            self.manual_failure_hbox_dict[new_sample].layout.display = 'inline-flex'


def main(workspace_namespace, workspace_name, workspace_bucket_name):
    entity_types_response = fapi.list_entity_types(workspace_namespace, workspace_name)
    if not entity_types_response.ok:
        raise RuntimeError(f'ERROR: {entity_types_response.text}')

    entity_types_dict = json.loads(entity_types_response.text)
    available_tables = [t for t in entity_types_dict if all(
        x in entity_types_dict[t]['attributeNames'] for x in ['batch_all_results', 'time_sample_set_updated'])]

    if len(available_tables) == 1:
        app = ResultsModificationGUI(workspace_namespace, workspace_name, workspace_bucket_name, available_tables[0])
        app.run()
    elif len(available_tables) > 1:
        app = TableSelectionGUI(workspace_namespace, workspace_name, workspace_bucket_name, available_tables)
        app.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace_namespace", dest="workspace_namespace", required=True)
    parser.add_argument("--worspace_name", dest="workspace_name", required=True)
    parser.add_argument("--worspace_bucket_name", dest="workspace_bucket_name", required=True)
    parser.parse_args()
    main(parser.workspace_namespace, parser.workspace_name, parser.workspace_bucket_name)
