import random
import ipywidgets
import pandas as pd
import ImputationPipeline.ManualQCPRS.ManualQCPRS as ManualQCPRS
import os
import itertools
import unittest
from unittest.mock import patch, MagicMock
import string
import firecloud.api as fapi
from google.cloud import storage
from io import StringIO
from datetime import datetime
import pytz
import logging
from traitlets import log

current_dir = os.path.dirname(__file__)
resources_dir = os.path.join(current_dir, "resources")


# exceptions thrown in callbacks will be logged as errors by ipywidgets
# we want them to be raised, so need to add a handler

class RaiseWarningAsErrorHandler(logging.Handler):

    def __init__(self):
        logging.Handler.__init__(self)
        self.level = logging.WARNING

    def emit(self, record):
        msg = self.format(record)
        raise RuntimeError(msg)


log = log.get_logger()
log.addHandler(RaiseWarningAsErrorHandler())


class ResultsModificationTest(unittest.TestCase):
    workspace_namespace = None
    workspace_name = None
    maxDiff = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.workspace_namespace = os.getenv("TEST_WORKSPACE_NAMESPACE")
        if cls.workspace_namespace is None:
            raise (RuntimeError("Environment variable TEST_WORKSPACE_NAMESPACE must be set to run tests"))
        user = fapi.whoami().split("@")[0]
        cls.workspace_name = f'test_workspace_ManualQCPRSTests_{user}'
        print(f"Creating test workspace {cls.workspace_namespace}/{cls.workspace_name}")
        fapi.unlock_workspace(cls.workspace_namespace, cls.workspace_name)
        fapi.delete_workspace(cls.workspace_namespace, cls.workspace_name)
        r = fapi.create_workspace(cls.workspace_namespace, cls.workspace_name)
        fapi._check_response_code(r, 201)
        cls.workspace_bucket_name = r.json()["bucketName"]
        # set delivery bucket to workspace bucket
        r = fapi.update_workspace_attributes(cls.workspace_namespace, cls.workspace_name,
                                             [fapi._attr_set("delivery-bucket", cls.workspace_bucket_name)])
        fapi._check_response_code(r, 200)

        # upload test_results into bucket
        cls.storage_client = storage.Client()
        bucket = cls.storage_client.bucket(cls.workspace_bucket_name)
        cls.test_results_sets = []
        input_prefix = "input_results"
        for test_file, lab_batch, group_n in zip(["test_results.tsv", "test_results_2.tsv", "test_results_3.tsv"],
                                                 ["BATCH_12345", "BATCH_12345", "BATCH_123456"], [1, 2, 1]):
            blob = bucket.blob(os.path.join(input_prefix, test_file))
            blob.upload_from_filename(os.path.join(resources_dir, test_file))
            cls.test_results_sets.append(("gs://" + os.path.join(cls.workspace_bucket_name, input_prefix, test_file),
                                          lab_batch, group_n
                                          )
                                         )

    def setUp(self):
        # add sample sets table to workspace
        with StringIO() as sample_sets_io:
            sample_sets_io.write('\t'.join(['entity:sample_set_id', 'lab_batch', 'group', 'batch_all_results',
                                            'time_sample_set_updated', 'delivered', 'redeliver']) + '\n')
            now = str(datetime.now(pytz.timezone('US/Eastern')))
            for results_path, lab_batch, group in self.test_results_sets:
                set_id = lab_batch if group == 1 else f'{lab_batch}_group_{group}'
                delivered = lab_batch == "BATCH_12345" and group == 1
                sample_sets_io.write(
                    '\t'.join([set_id, lab_batch, str(group), results_path, now, str(delivered), "false"]) + '\n')

            upload_sample_sets_response = fapi.upload_entities_tsv(self.workspace_namespace,
                                                                   self.workspace_name,
                                                                   sample_sets_io,
                                                                   "flexible")
            if not upload_sample_sets_response.ok:
                raise RuntimeError(f'ERROR: {upload_sample_sets_response.text}')

        self.results_modification_gui = ManualQCPRS.ResultsModificationGUI(self.workspace_namespace,
                                                                           self.workspace_name,
                                                                           f'gs://{self.workspace_bucket_name}',
                                                                           "sample_set")
        self.results_modification_gui.run()

    @classmethod
    def tearDownClass(cls) -> None:
        print(f"Deleting test workspace {cls.workspace_namespace}/{cls.workspace_name}")
        fapi.delete_workspace(cls.workspace_namespace, cls.workspace_name)

    def tearDown(self) -> None:
        fapi.delete_entities_of_type(self.workspace_namespace, self.workspace_name, "sample_set")
        # clear the top level of the bucket (for potential delivered results),
        # and the manually_qcd_prs_results prefix (for saved results)
        # we do not want to clear the whole bucket, since we want to keep the input_results prefix
        blobs = itertools.chain(self.storage_client.list_blobs(self.workspace_bucket_name, delimiter="/"),
                                self.storage_client.list_blobs(self.workspace_bucket_name,
                                                               prefix="manually_qcd_prs_results/",
                                                               delimiter="/"
                                                               )
                                )
        for blob in blobs:
            blob.delete()

    def assert_widget_status(self, batch_selection_value='', batch_selection_enabled=True, load_enabled=True,
                             save_enabled=True, discard_enabled=True):
        self.assertEqual(batch_selection_value, self.results_modification_gui.agg_batch_selection_dropdown.value)
        self.assertEqual(batch_selection_enabled,
                         not self.results_modification_gui.agg_batch_selection_dropdown.disabled)
        self.assertEqual(load_enabled, not self.results_modification_gui.agg_batch_download_button.disabled)
        self.assertEqual(save_enabled, not self.results_modification_gui.agg_batch_upload_button.disabled)
        self.assertEqual(discard_enabled, not self.results_modification_gui.discard_button.disabled)

    def test_initialization(self):
        self.assertEqual(self.results_modification_gui.delivery_bucket, self.workspace_bucket_name)

    def test_correct_agg_batches_available(self):
        expected_agg_batches = {"BATCH_12345_group_2", "BATCH_123456"}  # delivered batch not available
        self.assertEqual(expected_agg_batches, set(self.results_modification_gui.agg_batch_selection_dropdown.options))

    def test_available_for_redelivery(self):
        # set to redeliver and check that available
        fapi.update_entity(self.workspace_namespace, self.workspace_name, "sample_set", "BATCH_12345",
                           [fapi._attr_set("redeliver", True)])
        self.results_modification_gui.initialize_workspace_info_and_gcloud()  # need to redo this to get infor again
        self.results_modification_gui.run()
        expected_agg_batches = {"BATCH_12345_group_2", "BATCH_123456", "BATCH_12345"}
        self.assertEqual(expected_agg_batches, set(self.results_modification_gui.agg_batch_selection_dropdown.options))

    def test_initial_widget_status(self):
        self.assert_widget_status(load_enabled=False, save_enabled=False, discard_enabled=False)

    def test_select_batch(self):
        batch = "BATCH_12345_group_2"
        self.results_modification_gui.agg_batch_selection_dropdown.value = batch
        self.assert_widget_status(batch_selection_value=batch, save_enabled=False, discard_enabled=False)

    @patch("ImputationPipeline.ManualQCPRS.ManualQCPRS.SelectedBatchModificationGui")
    def test_load_batch(self, mock_selected_batch_modification_gui: MagicMock):

        batch = "BATCH_12345_group_2"
        self.results_modification_gui.agg_batch_selection_dropdown.value = batch
        self.results_modification_gui.agg_batch_download_button.click()

        expected_input_df = pd.read_csv(os.path.join(resources_dir, "test_results_2.tsv"), delimiter='\t')\
            .fillna("NA").set_index('sample_id')
        expected_lab_batch = "BATCH_12345"
        # cannot use assert_called_with because dataframe equality comparison doesn't work
        # so must compare arguments individually
        pd.testing.assert_frame_equal(expected_input_df, mock_selected_batch_modification_gui.call_args[0][0])
        self.assertEqual(expected_lab_batch, mock_selected_batch_modification_gui.call_args[0][1])
        # select_batch_gui has been run
        self.results_modification_gui.selected_batch_gui.run.assert_called()

        self.assert_widget_status(batch_selection_value=batch, batch_selection_enabled=False, load_enabled=False)

    @patch("ImputationPipeline.ManualQCPRS.ManualQCPRS.DeliveryConfirmationBox")
    def test_save_batch(self, mock_delivery_confirmation_box: MagicMock):
        batch = "BATCH_12345_group_2"
        self.results_modification_gui.agg_batch_selection_dropdown.value = batch
        self.results_modification_gui.agg_batch_download_button.click()
        self.results_modification_gui.agg_batch_upload_button.click()
        self.assert_widget_status(load_enabled=False, save_enabled=False, discard_enabled=False)
        self.assertIsNone(self.results_modification_gui.selected_batch_gui)

        # check table updated
        saved_path = fapi.get_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch)\
            .json()['attributes']['qcd_batch_results']
        expected_regex = "gs://" + os.path.join(self.workspace_bucket_name, "manually_qcd_prs_results",
                                                batch + "_\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}_E[D,S]T_manually_qcd_prs_results.csv")
        self.assertRegex(saved_path, expected_regex)

        # check saved result as expected
        saved_df = pd.read_csv(saved_path).fillna("NA")
        expected_df = pd.read_csv(os.path.join(resources_dir, "test_results_2.tsv"), sep="\t")
        expected_df = expected_df.reindex(columns=expected_df.columns.tolist() + ['poly_test_not_performed_reason',
                                                                                  'notes']).fillna("NA")
        pd.testing.assert_frame_equal(expected_df, saved_df, check_dtype=False)

        # check delivery confirmation box
        mock_delivery_confirmation_box.assert_called_with(saved_path, self.results_modification_gui)
        self.results_modification_gui.delivery_confirmation_box.run.assert_called()

    def test_discarded_results(self):
        batch = "BATCH_12345_group_2"
        self.results_modification_gui.agg_batch_selection_dropdown.value = batch
        self.results_modification_gui.agg_batch_download_button.click()
        self.results_modification_gui.discard_button.click()
        self.assert_widget_status(load_enabled=False, save_enabled=False, discard_enabled=False)
        self.assertIsNone(self.results_modification_gui.selected_batch_gui)

    def test_deliver_results(self):
        batch = "BATCH_12345_group_2"
        self.results_modification_gui.agg_batch_selection_dropdown.value = batch
        self.results_modification_gui.agg_batch_download_button.click()
        self.results_modification_gui.agg_batch_upload_button.click()
        self.results_modification_gui.delivery_confirmation_box.deliver_button.click()

        # check for delivery
        saved_path = fapi.get_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch)\
            .json()['attributes']['qcd_batch_results']

        delivered_path = f'gs://{self.workspace_bucket_name}/{os.path.basename(saved_path)}'
        delivered_df = pd.read_csv(delivered_path).fillna("NA")
        expected_df = pd.read_csv(os.path.join(resources_dir, "test_results_2.tsv"), sep="\t")
        expected_df = expected_df.reindex(
            columns=expected_df.columns.tolist() + ['poly_test_not_performed_reason', 'notes']).fillna("NA")
        pd.testing.assert_frame_equal(expected_df, delivered_df)

        # check that table is updated correctly
        attributes = fapi.get_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch)\
            .json()['attributes']
        self.assertTrue(attributes['delivered'])
        self.assertFalse(attributes['redeliver'])

        # check widget status
        self.assert_widget_status(load_enabled=False, save_enabled=False, discard_enabled=False)

        # check that this batch has been removed from options
        expected_agg_batches = {"BATCH_123456"}  # delivered batch not available
        self.assertEqual(expected_agg_batches, set(self.results_modification_gui.agg_batch_selection_dropdown.options))

    def test_dont_deliver(self):
        batch = "BATCH_12345_group_2"
        self.results_modification_gui.agg_batch_selection_dropdown.value = batch
        self.results_modification_gui.agg_batch_download_button.click()
        self.results_modification_gui.agg_batch_upload_button.click()
        self.results_modification_gui.delivery_confirmation_box.dont_deliver_button.click()

        # check that not delivered
        saved_path = fapi.get_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch)\
            .json()['attributes']['qcd_batch_results']

        delivered_path = f'gs://{self.workspace_bucket_name}/{os.path.basename(saved_path)}'
        bucket = self.storage_client.bucket(self.workspace_bucket_name)
        blob_should_not_exits = bucket.blob(delivered_path)
        self.assertFalse(blob_should_not_exits.exists())

        # check that delivered is false
        delivered = fapi.get_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch).json()[
            'attributes']['delivered']
        self.assertFalse(delivered)

    def test_deliver_redelivery(self):
        batch = "BATCH_12345"
        # set to redeliver and check that available
        fapi.update_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch,
                           [fapi._attr_set("redeliver", True)])
        self.results_modification_gui.initialize_workspace_info_and_gcloud()  # need to redo this to get infor again
        self.results_modification_gui.run()
        self.results_modification_gui.agg_batch_selection_dropdown.value = batch
        self.results_modification_gui.agg_batch_download_button.click()
        self.results_modification_gui.agg_batch_upload_button.click()
        self.results_modification_gui.delivery_confirmation_box.deliver_button.click()

        # check for delivery
        saved_path = \
            fapi.get_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch)\
                .json()['attributes']['qcd_batch_results']

        delivered_path = f'gs://{self.workspace_bucket_name}/{os.path.basename(saved_path)}'
        delivered_df = pd.read_csv(delivered_path).fillna("NA")
        expected_df = pd.read_csv(os.path.join(resources_dir, "test_results.tsv"), sep="\t")
        expected_df = expected_df.reindex(
            columns=expected_df.columns.tolist() + ['poly_test_not_performed_reason', 'notes']).fillna("NA")
        pd.testing.assert_frame_equal(expected_df, delivered_df)

        # check that table is updated correctly
        attributes = fapi.get_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch).json()[
            'attributes']
        self.assertTrue(attributes['delivered'])
        self.assertFalse(attributes['redeliver'])

        # check widget status
        self.assert_widget_status(load_enabled=False, save_enabled=False, discard_enabled=False)

        # check that this batch has been removed from options
        expected_agg_batches = {"BATCH_123456", "BATCH_12345_group_2"}  # delivered batch not available
        self.assertEqual(expected_agg_batches, set(self.results_modification_gui.agg_batch_selection_dropdown.options))

    def test_dont_deliver_redelivery(self):
        batch = "BATCH_12345"
        # set to redeliver and check that available
        fapi.update_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch,
                           [fapi._attr_set("redeliver", True)])
        self.results_modification_gui.initialize_workspace_info_and_gcloud()  # need to redo this to get infor again
        self.results_modification_gui.run()
        self.results_modification_gui.agg_batch_selection_dropdown.value = batch
        self.results_modification_gui.agg_batch_download_button.click()
        self.results_modification_gui.agg_batch_upload_button.click()
        self.results_modification_gui.delivery_confirmation_box.dont_deliver_button.click()

        # check that table has not changed
        attributes = fapi.get_entity(self.workspace_namespace, self.workspace_name, "sample_set", batch).json()[
            'attributes']
        self.assertTrue(attributes['delivered'])
        self.assertTrue(attributes['redeliver'])

        # check that this batch has NOT been removed from options
        expected_agg_batches = {"BATCH_123456", "BATCH_12345_group_2", "BATCH_12345"}  # delivered batch not available
        self.assertEqual(expected_agg_batches, set(self.results_modification_gui.agg_batch_selection_dropdown.options))


class FailSampleForConditionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = pd.read_csv(os.path.join(resources_dir, "test_results.tsv"), delimiter='\t'). \
            fillna("NA").set_index('sample_id')
        cls.lab_batch = "BATCH_12345"

    def setUp(self) -> None:
        self.selected_batch_gui = ManualQCPRS.SelectedBatchModificationGui(self.results, self.lab_batch)
        self.selected_batch_gui.run()
        self.selected_batch_gui.finished_addition_imputation_failures_button.click()
        self.selected_batch_gui.finished_sample_failures_button.click()
        self.selected_batch_gui.finished_condition_failures_button.click()
        self.default_reason = self.selected_batch_gui.standard_failure_reasons[0]
        self.other_reason = 'Other'

    def assert_modified_as_expected(self, failure_dict: dict = None):
        if failure_dict is None:
            failure_dict = {}
        # first assert that all samples without manual failures remain unchanged
        expected_unfailed_samples_df = self.results.query("sample_id not in @failure_dict.keys()").copy(). \
            reindex(columns=self.results.columns.tolist() + ['poly_test_not_performed_reason', 'notes']).fillna("NA")
        pd.testing.assert_frame_equal(
            self.selected_batch_gui.modified_results.query("sample_id not in @failure_dict.keys()").copy().fillna(
                "NA"),
            expected_unfailed_samples_df, check_like=True, check_dtype=False)

        # next assert that each sample with failures is correct
        for sample, condition_and_reason in failure_dict.items():
            self.assert_sample_failed_correct_conditions(sample, condition_and_reason)

    def assert_sample_failed_correct_conditions(self, sample, failed_conditions_with_reasons: list):
        this_sample_input_series = self.results.loc[sample]
        this_sample_modified_series = self.selected_batch_gui.modified_results.loc[sample]
        # first assert that all unfailed condition columns remain unchanged
        failed_condition_columns = [f'{condition}_{metric}' for condition, metric in
                                    itertools.product(
                                        [condition for condition, reason in failed_conditions_with_reasons],
                                        ["raw", "adjusted", "percentile", "risk", "reason_not_resulted"])]
        expected_unfailed_results = this_sample_input_series[~this_sample_input_series.index.
            isin(failed_condition_columns)]
        expected_unfailed_results = pd.concat(
            [expected_unfailed_results, pd.Series({'poly_test_not_performed_reason': 'NA',
                                                   'notes': 'NA'})])
        unfailed_result = this_sample_modified_series[
            ~this_sample_modified_series.index.isin(failed_condition_columns)
        ].fillna("NA")
        # need to put in same order
        expected_unfailed_results = expected_unfailed_results.reindex(unfailed_result.index.tolist())
        pd.testing.assert_series_equal(
            unfailed_result,
            expected_unfailed_results, check_names=False)

        # finally, assert that failed conditions are correctly failed
        data = {f'{c_and_r[0]}_{metric}': "NA" for c_and_r, metric in
                itertools.product(failed_conditions_with_reasons, ["raw", "adjusted", "percentile"])
                }
        data.update({f'{condition}_risk': "NOT_RESULTED" for condition, _ in failed_conditions_with_reasons})
        data.update(
            {f'{condition}_reason_not_resulted': reason for condition, reason in failed_conditions_with_reasons})
        expected_failed_results = pd.Series(data)
        failed_results = this_sample_modified_series[failed_condition_columns]
        expected_failed_results = expected_failed_results.reindex(failed_results.index.tolist())
        pd.testing.assert_series_equal(expected_failed_results, failed_results, check_names=False)

    def fail_sample_for_condition(self, sample, condition, reason=None):
        self.selected_batch_gui.sample_selection.value = sample
        sample_hbox: ipywidgets.HBox = self.selected_batch_gui.manual_failure_hbox_dict[sample]
        sample_condition_vbox: ipywidgets.VBox = next(vbox for vbox in sample_hbox.children if
                                                      vbox.children[0].description == condition)

        self.assertIsNotNone(sample_condition_vbox)
        status_selection_dropdown: ipywidgets.Dropdown = sample_condition_vbox.children[1]
        if not reason or reason == self.default_reason:
            status_selection_dropdown.value = self.default_reason
        else:
            status_selection_dropdown.value = self.other_reason
            custom_reason_text: ipywidgets.Text = sample_condition_vbox.children[2]
            custom_reason_text.value = reason

    def unfail_sample_for_condition(self, sample, condition):
        self.selected_batch_gui.sample_selection.value = sample
        sample_hbox: ipywidgets.HBox = self.selected_batch_gui.manual_failure_hbox_dict[sample]
        sample_condition_vbox: ipywidgets.VBox = next(vbox for vbox in sample_hbox.children if
                                                      vbox.children[0].description == condition)

        self.assertIsNotNone(sample_condition_vbox)
        status_selection_dropdown: ipywidgets.Dropdown = sample_condition_vbox.children[1]
        status_selection_dropdown.value = "PASS"

    def assert_condition_widget_box_status(self, condition_vbox, condition_status_dictionary: dict):
        # each condition vbox should contain 3 children
        children = condition_vbox.children
        self.assertEqual(len(children), 3)

        # check that name_button status is correct
        name_button: ipywidgets.Button = children[0]
        # check label
        self.assertTrue(name_button.description in condition_status_dictionary)
        condition_status = condition_status_dictionary[name_button.description]
        # check style
        expected_style = condition_status[0]
        self.assertEqual(name_button.button_style, expected_style)
        # name_button is always enabled (though clicking does nothing)
        self.assertFalse(name_button.disabled)

        # check that failure reason widgets based on status of condition
        status_selection: ipywidgets.Dropdown = children[1]
        custom_failure_entry: ipywidgets.Text = children[2]
        # check options
        if expected_style == "info":
            # set to not calculated, both disabled
            self.assertEqual(status_selection.options, ("NOT CALCULATED",))
            self.assertEqual(status_selection.value, "NOT CALCULATED")
            self.assertEqual(custom_failure_entry.value, "")
            self.assertTrue(status_selection.disabled)
            self.assertTrue(custom_failure_entry.disabled)
        else:
            # check options
            self.assertEqual(status_selection.options, (('PASS', 'PASS'),
                                                        ('FAIL: Too Many Scoring Sites Missing',
                                                         'Too Many Scoring Sites Missing'),
                                                        ('FAIL: Other Reason', 'Other')
                                                        )
                             )
            # status selection enabled
            self.assertFalse(status_selection.disabled)
            # check status selection value, custom reason value, whether custom reason is enabled
            if len(condition_status) == 1:
                # no failure reason given
                # status_selection is PASS, custom reason is disabled
                self.assertEqual(status_selection.value, "PASS")
                self.assertEqual(custom_failure_entry.value, "")
                self.assertTrue(custom_failure_entry.disabled)
            else:
                expected_reason = condition_status[1]
                if expected_reason == self.default_reason:
                    # status_selection enabled, custom reason disabled
                    self.assertEqual(status_selection.value, expected_reason)
                    self.assertEqual(custom_failure_entry.value, "")
                    self.assertTrue(custom_failure_entry.disabled)
                elif expected_reason == '':
                    # empty string means that "other reason" has been selected, but not yet entered
                    self.assertEqual(status_selection.value, self.other_reason)
                    self.assertEqual(custom_failure_entry.value, "")
                    self.assertEqual(custom_failure_entry.placeholder, 'Enter Custom Failure Reason')
                    self.assertFalse(custom_failure_entry.disabled)
                else:
                    self.assertEqual(status_selection.value, self.other_reason)
                    self.assertEqual(custom_failure_entry.value, expected_reason)
                    self.assertFalse(custom_failure_entry.disabled)

    def assert_widget_status(self, selected_sample='', condition_status_dictionary: dict = None):
        if condition_status_dictionary is None:
            condition_status_dictionary = {}
        self.assertEqual(self.selected_batch_gui.sample_selection.value, selected_sample)
        self.assertFalse(self.selected_batch_gui.sample_selection.disabled)
        for sample, failure_hbox in self.selected_batch_gui.manual_failure_hbox_dict.items():
            if sample == selected_sample:
                self.assertEqual(failure_hbox.layout.display, "inline-flex")
                failure_hbox_children = failure_hbox.children
                self.assertEqual(len(failure_hbox_children), len(condition_status_dictionary))
                expected_conditions = set(condition_status_dictionary.keys())
                widget_conditions = {v.children[0].description for v in failure_hbox_children}
                self.assertEqual(widget_conditions, expected_conditions)
                # loop through children, each child is a vbox for 1 condition
                for child_vbox in failure_hbox_children:
                    self.assert_condition_widget_box_status(child_vbox, condition_status_dictionary)

            else:
                self.assertEqual(failure_hbox.layout.display, "none")

    def test_correct_samples_available(self):
        self.assertEqual(set(self.results.query("not is_control_sample").index),
                         set(self.selected_batch_gui.sample_selection.options))

    def test_initial_widget_status(self):
        self.assert_widget_status()

    def test_select_sample_widget_status(self):
        sample_id = "sample_6"
        expected_condition_status_dict = {"condition_1": ("success",), "condition_2": ("info",),
                                          "condition_3": ("success",)}
        self.selected_batch_gui.sample_selection.value = sample_id
        self.assert_widget_status(selected_sample=sample_id, condition_status_dictionary=expected_condition_status_dict)

    def test_fail_default_reason(self):
        sample_id = "sample_6"
        condition_to_fail = "condition_1"
        expected_condition_status_dict = {"condition_1": ("danger", self.default_reason), "condition_2": ("info",),
                                          "condition_3": ("success",)}
        expected_failure_dict = {sample_id: [(condition_to_fail, self.default_reason)]}
        self.fail_sample_for_condition(sample_id, condition_to_fail)
        self.assert_widget_status(selected_sample=sample_id, condition_status_dictionary=expected_condition_status_dict)
        self.assert_modified_as_expected(expected_failure_dict)

    def test_select_other_reason(self):
        sample_id = "sample_6"
        condition = "condition_1"
        # when "Other reason" has been selected, but not yet entered, condition not yet failed
        expected_condition_status_dict = {"condition_1": ("success", ""), "condition_2": ("info",),
                                          "condition_3": ("success",)}
        self.selected_batch_gui.sample_selection.value = sample_id
        condition_vbox = next(vbox for vbox in self.selected_batch_gui.manual_failure_hbox_dict[sample_id].children if
                              vbox.children[0].description == condition)
        condition_status: ipywidgets.Dropdown = condition_vbox.children[1]
        condition_status.value = self.other_reason
        self.assert_widget_status(selected_sample=sample_id, condition_status_dictionary=expected_condition_status_dict)
        self.assert_modified_as_expected()  # no modifications yet

    def test_fail_sample_other_reason(self):
        sample_id = "sample_6"
        condition = "condition_1"
        reason = "a different reason"
        self.fail_sample_for_condition(sample_id, condition, reason)
        expected_condition_status_dict = {"condition_1": ("danger", reason), "condition_2": ("info",),
                                          "condition_3": ("success",)}
        self.assert_widget_status(sample_id, expected_condition_status_dict)
        expected_failure_dict = {sample_id: [(condition, reason)]}
        self.assert_modified_as_expected(expected_failure_dict)

    def test_fail_then_change_reason(self):
        sample_id = "sample_6"
        condition = "condition_1"
        reason = "a different reason"
        self.fail_sample_for_condition(sample_id, condition, self.default_reason)
        expected_condition_status_dict = {"condition_1": ("danger", self.default_reason), "condition_2": ("info",),
                                          "condition_3": ("success",)}
        self.assert_widget_status(sample_id, expected_condition_status_dict)
        expected_failure_dict = {sample_id: [(condition, self.default_reason)]}
        self.assert_modified_as_expected(expected_failure_dict)

        # change reason
        self.fail_sample_for_condition(sample_id, condition, reason)
        expected_condition_status_dict = {"condition_1": ("danger", reason), "condition_2": ("info",),
                                          "condition_3": ("success",)}
        self.assert_widget_status(sample_id, expected_condition_status_dict)
        expected_failure_dict = {sample_id: [(condition, reason)]}
        self.assert_modified_as_expected(expected_failure_dict)

    def test_fail_then_unfail(self):
        sample_id = "sample_6"
        condition = "condition_1"
        self.fail_sample_for_condition(sample_id, condition, self.default_reason)
        expected_condition_status_dict = {"condition_1": ("danger", self.default_reason), "condition_2": ("info",),
                                          "condition_3": ("success",)}
        self.assert_widget_status(sample_id, expected_condition_status_dict)
        expected_failure_dict = {sample_id: [(condition, self.default_reason)]}
        self.assert_modified_as_expected(expected_failure_dict)

        self.unfail_sample_for_condition(sample_id, condition)
        expected_condition_status_dict = {"condition_1": ("success",), "condition_2": ("info",),
                                          "condition_3": ("success",)}
        self.assert_widget_status(sample_id, expected_condition_status_dict)
        self.assert_modified_as_expected()

    def test_fail_then_select_new_sample(self):
        sample_id = "sample_6"
        condition = "condition_1"
        self.fail_sample_for_condition(sample_id, condition, self.default_reason)

        new_sample = "sample_2"
        self.selected_batch_gui.sample_selection.value = new_sample
        expected_condition_status_dict = {"condition_1": ("success",), "condition_2": ("success",),
                                          "condition_3": ("success",)}
        self.assert_widget_status(new_sample, expected_condition_status_dict)

    def test_fail_then_new_sample_then_failed_sample(self):
        sample_id = "sample_6"
        condition = "condition_1"
        self.fail_sample_for_condition(sample_id, condition, self.default_reason)

        new_sample = "sample_2"
        self.selected_batch_gui.sample_selection.value = new_sample
        self.selected_batch_gui.sample_selection.value = sample_id
        expected_condition_status_dict = {"condition_1": ("danger", self.default_reason), "condition_2": ("info",),
                                          "condition_3": ("success",)}
        self.assert_widget_status(sample_id, expected_condition_status_dict)

    def test_sequence(self):
        action_sequence = [{"sample": "sample_4", "condition": "condition_2", "reason": self.default_reason},
                           {"sample": "sample_2", "condition": "condition_1", "reason": self.default_reason},
                           {"sample": "sample_4", "condition": "condition_1", "reason": "another reason"},
                           {"sample": "sample_1", "condition": "condition_2", "reason": self.default_reason},
                           {"sample": "sample_4", "condition": "condition_1", "op": "unfail"},
                           {"sample": "sample_2", "condition": "condition_1", "reason": "a different reason"},
                           {"sample": "sample_6", "condition": "condition_1", "reason": self.default_reason},
                           {"sample": "sample_6", "condition": "condition_3", "reason": "some other reason"},
                           {"sample": "sample_6", "condition": "condition_3", "reason": self.default_reason},
                           {"sample": "sample_6", "condition": "condition_3", "op": "unfail"},
                           {"sample": "sample_4", "condition": "condition_1", "reason": self.default_reason}]

        # add a bunch more random actions to test in sequence
        possible_actions = [{"sample": f'sample_{i}', "condition": f'condition_{j}', "reason": self.default_reason}
                            for i, j in itertools.product(range(1, 7), range(1, 4))
                            if i != 3 and (i, j) not in [(4, 3), (5, 3), (6, 2)]]

        random.seed(1234)
        reason_swap = {self.default_reason: "some other reason", "some other reason": self.default_reason}
        for i in range(100):
            action = random.choice(possible_actions)
            action_sequence.append(action.copy())
            if "reason" in action:
                action["reason"] = reason_swap[action["reason"]]
                unfail_action = {"sample": action["sample"], "condition": action["condition"], "op": "unfail"}
                if unfail_action not in possible_actions:
                    possible_actions.append(unfail_action)
            elif "op" in action:
                possible_actions.remove(action)

        expected_failure_dict = dict()
        for i, action in enumerate(action_sequence):
            sample = action.get("sample")
            condition = action.get("condition")
            op = action.get("op", "fail")
            if op == "unfail":
                self.unfail_sample_for_condition(sample, condition)
                expected_failure_dict[sample] = [c_and_r for c_and_r in expected_failure_dict[sample] if
                                                 c_and_r[0] != condition]
                if len(expected_failure_dict[sample]) == 0:
                    del expected_failure_dict[sample]
            else:
                reason = action.get("reason")
                self.fail_sample_for_condition(sample, condition, reason)
                if sample in expected_failure_dict:
                    expected_failure_dict[sample] = [c_and_r for c_and_r in expected_failure_dict[sample] if
                                                     c_and_r[0] != condition]
                    expected_failure_dict[sample].append((condition, reason))
                else:
                    expected_failure_dict[sample] = [(condition, reason)]

            self.assert_modified_as_expected(expected_failure_dict)


class FailConditionForAllSamplesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = pd.read_csv(os.path.join(resources_dir, "test_results.tsv"), delimiter='\t'). \
            fillna("NA").set_index('sample_id')
        cls.lab_batch = "BATCH_12345"

    def setUp(self) -> None:
        self.selected_batch_gui = ManualQCPRS.SelectedBatchModificationGui(self.results, self.lab_batch)
        self.selected_batch_gui.run()
        self.selected_batch_gui.finished_addition_imputation_failures_button.click()
        self.selected_batch_gui.finished_sample_failures_button.click()
        self.default_reason = self.selected_batch_gui.condition_failure_reason_selection.options[0]

    def assert_modified_as_expected(self, failed_conditions: dict = None):
        if failed_conditions is None:
            failed_conditions = {}
        # first assert that all unfailed condition columns remain unchanged
        failed_condition_columns = [f'{condition}_{metric}' for condition, metric in
                                    itertools.product(failed_conditions.keys(),
                                                      ["raw", "adjusted", "percentile", "risk", "reason_not_resulted"])]
        expected_unfailed_results_df = self.results.drop(failed_condition_columns, axis=1).copy()
        expected_unfailed_results_df = expected_unfailed_results_df. \
            reindex(
            columns=expected_unfailed_results_df.columns.tolist() + ['poly_test_not_performed_reason', 'notes']). \
            fillna("NA")
        pd.testing.assert_frame_equal(
            self.selected_batch_gui.modified_results.drop(failed_condition_columns, axis=1).copy().fillna("NA"),
            expected_unfailed_results_df, check_dtype=False)

        # next assert that each failed condition is correctly failed
        if failed_conditions:
            self.assert_conditions_failed_all_samples(failed_conditions)

    def assert_conditions_failed_all_samples(self, failed_conditions_with_reasons: dict):
        failed_condition_columns = [f'{condition}_{metric}' for condition, metric in
                                    itertools.product(failed_conditions_with_reasons.keys(),
                                                      ["raw", "adjusted", "percentile", "risk", "reason_not_resulted"])]
        failed_conditions_results = self.selected_batch_gui.modified_results[failed_condition_columns]

        data = {f'{condition}_{metric}': "NA" for condition, metric in
                itertools.product(failed_conditions_with_reasons.keys(),
                                  ["raw", "adjusted", "percentile"])}

        data.update({f'{condition}_risk': "NOT_RESULTED" for condition in failed_conditions_with_reasons.keys()})
        data.update({f'{condition}_reason_not_resulted': reason for condition, reason in
                     failed_conditions_with_reasons.items()})
        # note, control sample is never failed
        data["sample_id"] = self.results.query("not is_control_sample").index
        expected_results = pd.DataFrame(data).set_index("sample_id")
        # add control sample back in unmodified
        expected_results = pd.concat(
            [expected_results, self.results.query("is_control_sample")[failed_condition_columns]])
        # unscored conditions need risk and reason reset to NA
        for sample, condition in itertools.product(self.results.query("not is_control_sample").index,
                                                   failed_conditions_with_reasons.keys()):
            if self.results.loc[sample, f'{condition}_risk'] == "NA":
                expected_results.loc[sample, f'{condition}_risk'] = "NA"
                expected_results.loc[sample, f'{condition}_reason_not_resulted'] = "NA"
        pd.testing.assert_frame_equal(expected_results, failed_conditions_results, check_like=True, check_dtype=False)

    def assert_widget_status(self, condition_value=None, reason_value='',
                             condition_enabled=True, reason_enabled=True, fail_enabled=True, unfail_enabled=True):
        self.assertEqual(self.selected_batch_gui.condition_failure_selection.value, condition_value)
        self.assertEqual(self.selected_batch_gui.condition_failure_reason_selection.value, reason_value)
        self.assertEqual(self.selected_batch_gui.condition_failure_selection.disabled, not condition_enabled)
        self.assertEqual(self.selected_batch_gui.condition_failure_reason_selection.disabled, not reason_enabled)
        self.assertEqual(self.selected_batch_gui.fail_condition_all_samples_button.disabled, not fail_enabled)
        self.assertEqual(self.selected_batch_gui.unfail_condition_all_samples_button.disabled, not unfail_enabled)
        # finish button always enabled
        self.assertFalse(self.selected_batch_gui.finished_condition_failures_button.disabled)

    def fail_condition(self, condition, reason):
        self.selected_batch_gui.condition_failure_selection.value = condition
        self.selected_batch_gui.condition_failure_reason_selection.value = reason
        self.selected_batch_gui.fail_condition_all_samples_button.click()

    def unfail_condition(self, condition):
        self.selected_batch_gui.condition_failure_selection.value = condition
        self.assertFalse(self.selected_batch_gui.unfail_condition_all_samples_button.disabled)
        self.selected_batch_gui.unfail_condition_all_samples_button.click()

    def test_correct_conditions_available(self):
        expected_conditions = {"condition_1", "condition_2", "condition_3", None}
        self.assertEqual(set(self.selected_batch_gui.condition_failure_selection.options), expected_conditions)

    def test_initial_widget_status(self):
        self.assert_widget_status(reason_enabled=False, fail_enabled=False, unfail_enabled=False)

    def test_select_condition(self):
        condition_to_fail = "condition_2"
        self.selected_batch_gui.condition_failure_selection.value = condition_to_fail
        # condition selection and reason should be enabled
        self.assert_widget_status(condition_value=condition_to_fail, fail_enabled=False, unfail_enabled=False)

    def test_select_reason(self):
        condition_to_fail = "condition_2"
        self.selected_batch_gui.condition_failure_selection.value = condition_to_fail
        self.selected_batch_gui.condition_failure_reason_selection.value = self.default_reason
        self.assert_widget_status(condition_value=condition_to_fail, reason_value=self.default_reason,
                                  unfail_enabled=False)

    def test_fail_condition(self):
        condition_to_fail = "condition_2"
        self.fail_condition(condition_to_fail, self.default_reason)
        self.assert_widget_status(reason_enabled=False, fail_enabled=False, unfail_enabled=False)
        self.assert_modified_as_expected({condition_to_fail: self.default_reason})

    def test_widget_status_failed_condition_selected(self):
        condition_to_fail = "condition_2"
        self.fail_condition(condition_to_fail, self.default_reason)
        self.selected_batch_gui.condition_failure_selection.value = condition_to_fail
        self.assert_widget_status(condition_value=condition_to_fail, reason_value=self.default_reason)

    def test_fail_then_change_reason(self):
        condition_to_fail = "condition_2"
        new_reason = "a different reason"
        self.fail_condition(condition_to_fail, self.default_reason)
        self.assert_modified_as_expected({condition_to_fail: self.default_reason})
        self.fail_condition(condition_to_fail, new_reason)
        self.assert_modified_as_expected({condition_to_fail: new_reason})

    def test_fail_then_unfail(self):
        condition_to_fail = "condition_2"
        self.fail_condition(condition_to_fail, self.default_reason)
        self.assert_modified_as_expected({condition_to_fail: self.default_reason})
        self.unfail_condition(condition_to_fail)
        self.assert_modified_as_expected()

    def test_fail_custom_reason_reselect(self):
        condition_to_fail = "condition_2"
        new_reason = "a different reason"
        self.fail_condition(condition_to_fail, new_reason)
        self.selected_batch_gui.condition_failure_selection.value = condition_to_fail
        self.assert_widget_status(condition_value=condition_to_fail, reason_value=new_reason)

    def test_sequence(self):
        action_sequence = [{"condition": "condition_3", "reason": self.default_reason},
                           {"condition": "condition_1", "reason": "a new reason"},
                           {"condition": "condition_1", "op": "unfail"},
                           {"condition": "condition_2", "reason": self.default_reason},
                           {"condition": "condition_3", "reason": "some other reason"},
                           {"condition": "condition_3", "op": "unfail"},
                           {"condition": "condition_3", "reason": self.default_reason}]

        condition_failures = dict()

        # add a bunch more random actions to test in sequence
        possible_actions = [{"condition": f'condition_{j}', "reason": self.default_reason} for j in range(1, 4)]

        random.seed(1234)
        reason_swap = {self.default_reason: "some other reason", "some other reason": self.default_reason}

        for i in range(100):
            action = random.choice(possible_actions)
            action_sequence.append(action.copy())
            if "reason" in action:
                action["reason"] = reason_swap[action["reason"]]
                unfail_action = {"condition": action["condition"], "op": "unfail"}
                if unfail_action not in possible_actions:
                    possible_actions.append(unfail_action)
            elif "op" in action:
                possible_actions.remove(action)

        for action in action_sequence:
            condition = action.get("condition")
            op = action.get("op", "fail")
            if op == "unfail":
                self.unfail_condition(condition)
                del condition_failures[condition]
            else:
                reason = action.get("reason")
                self.fail_condition(condition, reason)
                condition_failures[condition] = reason
            self.assert_modified_as_expected(condition_failures)

    def test_finish_button(self):
        with patch.object(self.selected_batch_gui.finished_condition_failures_button, 'close') as \
                finished_condition_failures_button_close_mock, \
                patch.object(self.selected_batch_gui.condition_failure_hbox, 'close') as \
                        condition_failure_hbox_close_mock, \
                patch.object(self.selected_batch_gui, 'build_fail_sample_for_condition_section') as \
                        build_fail_sample_for_condition_section_mock, \
                patch.object(self.selected_batch_gui, 'print_manual_failure_status') as \
                        print_manual_failure_status_mock:
            self.selected_batch_gui.finished_condition_failures_button.click()

        finished_condition_failures_button_close_mock.assert_called()
        condition_failure_hbox_close_mock.assert_called()
        build_fail_sample_for_condition_section_mock.assert_called()
        print_manual_failure_status_mock.assert_called()


class FailAllConditionsForSampleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = pd.read_csv(os.path.join(resources_dir, "test_results.tsv"), delimiter='\t'). \
            fillna("NA").set_index('sample_id')
        cls.lab_batch = "BATCH_12345"

    def setUp(self) -> None:
        self.selected_batch_gui = ManualQCPRS.SelectedBatchModificationGui(self.results, self.lab_batch)
        self.selected_batch_gui.run()
        self.selected_batch_gui.finished_addition_imputation_failures_button.click()
        self.default_reason = self.selected_batch_gui.sample_failure_reason_selection.options[0]

    def assert_samples_failed_all_conditions(self, failed_samples_with_reasons):
        failed_samples_results = self.selected_batch_gui.modified_results. \
            query("sample_id in @failed_samples_with_reasons.keys()").copy().fillna("NA")
        data = {f'condition_{i}_{metric}': "NA" for i, metric in
                itertools.product(range(1, 4), ["raw", "adjusted", "percentile"])
                }
        data.update({f'condition_{i}_risk': "NOT_RESULTED" for i in range(1, 4)})
        data["sample_id"] = failed_samples_with_reasons.keys()
        data["poly_test_not_performed_reason"] = "Failed PRS"
        data["lab_batch"] = self.lab_batch
        data["is_control_sample"] = False
        expected_results = pd.DataFrame(data).set_index("sample_id")
        expected_results["notes"] = expected_results.index.map(failed_samples_with_reasons)
        expected_results = expected_results.join(self.results[['PC1', 'PC2']])
        for i in range(1, 4):
            expected_results[f'condition_{i}_reason_not_resulted'] = \
                expected_results.index.map(failed_samples_with_reasons)

        # unscored conditions need to have risk and reason reset to "NA"
        for sample, i in itertools.product(failed_samples_with_reasons.keys(), range(1, 4)):
            if self.results.loc[sample, f'condition_{i}_risk'] == "NA":
                expected_results.loc[sample, f'condition_{i}_risk'] = "NA"
                expected_results.loc[sample, f'condition_{i}_reason_not_resulted'] = "NA"

        pd.testing.assert_frame_equal(expected_results, failed_samples_results, check_like=True, check_dtype=False)

    def assert_modified_as_expected(self, failed_samples=None):
        if failed_samples is None:
            failed_samples = {}
        # first assert that all entries for unfailed samples remain the same
        expected_unfailed_samples_df = self.results.query("sample_id not in @failed_samples.keys()").copy(). \
            reindex(columns=self.results.columns.tolist() + ['poly_test_not_performed_reason', 'notes']).fillna("NA")
        pd.testing.assert_frame_equal(
            self.selected_batch_gui.modified_results.query("sample_id not in @failed_samples.keys()").copy().fillna(
                "NA"),
            expected_unfailed_samples_df, check_like=True, check_dtype=False)

        # next assert that each failed sample is correctly failed
        if failed_samples:
            self.assert_samples_failed_all_conditions(failed_samples)

    def fail_sample(self, sample, reason):
        self.selected_batch_gui.sample_failure_selection.value = sample
        self.selected_batch_gui.sample_failure_reason_selection.value = reason
        self.selected_batch_gui.fail_sample_all_conditions_button.click()

    def unfail_sample(self, sample):
        self.selected_batch_gui.sample_failure_selection.value = sample
        self.assertFalse(self.selected_batch_gui.unfail_sample_all_conditions_button.disabled)
        self.selected_batch_gui.unfail_sample_all_conditions_button.click()

    def test_correct_samples_available(self):
        self.assertEqual(set(self.results.query("not is_control_sample").index),
                         set(self.selected_batch_gui.sample_failure_selection.options))

    def assert_widget_status(self, sample_value='', reason_value='',
                             sample_enabled=True, reason_enabled=True, fail_enabled=True, unfail_enabled=True):
        self.assertEqual(self.selected_batch_gui.sample_failure_selection.value, sample_value)
        self.assertEqual(self.selected_batch_gui.sample_failure_reason_selection.value, reason_value)
        self.assertEqual(self.selected_batch_gui.sample_failure_selection.disabled, not sample_enabled)
        self.assertEqual(self.selected_batch_gui.sample_failure_reason_selection.disabled, not reason_enabled)
        self.assertEqual(self.selected_batch_gui.fail_sample_all_conditions_button.disabled, not fail_enabled)
        self.assertEqual(self.selected_batch_gui.unfail_sample_all_conditions_button.disabled, not unfail_enabled)
        # finish button always enabled
        self.assertFalse(self.selected_batch_gui.finished_sample_failures_button.disabled)

    def test_initial_widget_status(self):
        self.assert_widget_status(reason_enabled=False, fail_enabled=False, unfail_enabled=False)

    def test_select_sample(self):
        sample_to_fail = "sample_2"
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail
        self.assert_widget_status(sample_value=sample_to_fail, fail_enabled=False, unfail_enabled=False)

    def test_select_reason(self):
        sample_to_fail = "sample_2"
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail
        self.selected_batch_gui.sample_failure_reason_selection.value = self.default_reason
        self.assert_widget_status(sample_value=sample_to_fail, unfail_enabled=False,
                                  reason_value=self.default_reason)

    def test_fail_sample(self):
        sample_to_fail = "sample_2"
        self.fail_sample(sample_to_fail, self.default_reason)
        self.assert_widget_status(reason_enabled=False, fail_enabled=False, unfail_enabled=False)

        self.assert_modified_as_expected({sample_to_fail: self.default_reason})

    def test_widget_status_failed_sample_selected(self):
        sample_to_fail = "sample_2"
        self.fail_sample(sample_to_fail, self.default_reason)

        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail
        self.assert_widget_status(sample_value=sample_to_fail,
                                  reason_value=self.default_reason)

    def test_fail_then_unfail(self):
        sample_to_fail = "sample_2"
        self.fail_sample(sample_to_fail, self.default_reason)
        self.unfail_sample(sample_to_fail)
        self.assert_widget_status(reason_enabled=False, fail_enabled=False, unfail_enabled=False)
        self.assert_modified_as_expected()

    def test_fail_then_change_reason(self):
        sample_to_fail = "sample_2"
        new_reason = "a new reason"
        self.fail_sample(sample_to_fail, self.default_reason)
        self.assert_modified_as_expected({sample_to_fail: self.default_reason})
        self.fail_sample(sample_to_fail, new_reason)
        self.assert_widget_status(reason_enabled=False, fail_enabled=False, unfail_enabled=False)
        self.assert_modified_as_expected({sample_to_fail: new_reason})

    def test_fail_custom_reason_reselect(self):
        sample_to_fail = "sample_2"
        new_reason = "a new reason"
        self.fail_sample(sample_to_fail, new_reason)
        self.assert_widget_status(reason_enabled=False, fail_enabled=False, unfail_enabled=False)
        self.assert_modified_as_expected({sample_to_fail: new_reason})
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail

        self.assert_widget_status(sample_value=sample_to_fail,
                                  reason_value=new_reason)

    def test_sequence(self):
        action_sequence = [{"sample": "sample_2", "reason": self.default_reason},
                           {"sample": "sample_1", "reason": "a special reason"},
                           {"sample": "sample_5", "reason": "a different_reason"},
                           {"sample": "sample_1", "op": "unfail"},
                           {"sample": "sample_6", "reason": self.default_reason},
                           {"sample": "sample_5", "reason": self.default_reason},
                           {"sample": "sample_2", "op": "unfail"}]

        sample_failures = dict()

        # add a bunch more random actions to test in sequence
        possible_actions = [{"sample": f'sample_{i}', "reason": self.default_reason}
                            for i in range(1, 7) if i != 3]

        random.seed(1234)
        reason_swap = {self.default_reason: "some other reason", "some other reason": self.default_reason}

        for i in range(100):
            action = random.choice(possible_actions)
            action_sequence.append(action.copy())
            if "reason" in action:
                action["reason"] = reason_swap[action["reason"]]
                unfail_action = {"sample": action["sample"], "op": "unfail"}
                if unfail_action not in possible_actions:
                    possible_actions.append(unfail_action)
            elif "op" in action:
                possible_actions.remove(action)

        for action in action_sequence:
            sample = action.get("sample")
            op = action.get("op", "fail")
            if op == "unfail":
                self.unfail_sample(sample)
                del sample_failures[sample]
            else:
                reason = action.get("reason")
                self.fail_sample(sample, reason)
                sample_failures[sample] = reason
            self.assert_modified_as_expected(sample_failures)

    def test_finish_button(self):
        with patch.object(self.selected_batch_gui.finished_sample_failures_button, 'close') as \
                finished_sample_failures_button_close_mock, \
                patch.object(self.selected_batch_gui.sample_failure_hbox, 'close') as sample_failure_hbox_close_mock, \
                patch.object(self.selected_batch_gui, 'build_fail_condition_section') as \
                        build_fail_condition_section_mock, \
                patch.object(self.selected_batch_gui, 'print_failed_conditions_all_samples') as \
                        print_failed_conditions_all_samples_mock:
            self.selected_batch_gui.finished_sample_failures_button.click()

        finished_sample_failures_button_close_mock.assert_called()
        sample_failure_hbox_close_mock.assert_called()
        build_fail_condition_section_mock.assert_called()
        print_failed_conditions_all_samples_mock.assert_called()


class AddSampleFailedImputationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = pd.read_csv(os.path.join(resources_dir, "test_results.tsv"), delimiter='\t'). \
            fillna("NA").set_index('sample_id')
        cls.lab_batch = "BATCH_12345"

    def setUp(self) -> None:
        self.selected_batch_gui = ManualQCPRS.SelectedBatchModificationGui(self.results, self.lab_batch)
        self.selected_batch_gui.run()

    def add_failed_imputation_sample(self, sample, note=''):
        self.selected_batch_gui.failed_imputation_sample_text_box.value = sample
        self.selected_batch_gui.failed_imputation_notes_text_box.value = note
        self.selected_batch_gui.failed_imputation_sample_add_button.click()

    def remove_failed_sample(self, sample):
        self.selected_batch_gui.failed_imputation_sample_text_box.value = sample
        self.assertFalse(self.selected_batch_gui.failed_imputation_sample_remove_button.disabled)
        self.selected_batch_gui.failed_imputation_sample_remove_button.click()

    def assert_correct_failed_samples_added(self, failed_samples_with_notes: dict = None):
        if failed_samples_with_notes is None:
            failed_samples_with_notes = {}
        expected_failed_samples_data = {"sample_id": [sample for sample in failed_samples_with_notes.keys()],
                                        "lab_batch": self.lab_batch,
                                        "is_control_sample": False,
                                        "poly_test_not_performed_reason": "Failed Imputation",
                                        "notes": [note for note in failed_samples_with_notes.values()]}
        expected_failed_samples_df = pd.DataFrame(expected_failed_samples_data).set_index("sample_id")
        if not expected_failed_samples_df.empty or not self.selected_batch_gui.failed_imputation_samples.empty:
            pd.testing.assert_frame_equal(expected_failed_samples_df, self.selected_batch_gui.failed_imputation_samples,
                                          check_like=True, check_dtype=False)

    def assert_widget_status(self, sample_value='', note_value='',
                             note_enabled=True, add_updated_enabled=True, delete_enabled=True):
        self.assertEqual(self.selected_batch_gui.failed_imputation_sample_text_box.value, sample_value)
        self.assertEqual(self.selected_batch_gui.failed_imputation_notes_text_box.value, note_value)
        self.assertFalse(self.selected_batch_gui.failed_imputation_sample_text_box.disabled)
        self.assertEqual(self.selected_batch_gui.failed_imputation_notes_text_box.disabled, not note_enabled)
        self.assertEqual(self.selected_batch_gui.failed_imputation_sample_add_button.disabled, not add_updated_enabled)
        self.assertEqual(self.selected_batch_gui.failed_imputation_sample_remove_button.disabled, not delete_enabled)
        self.assertFalse(self.selected_batch_gui.finished_addition_imputation_failures_button.disabled)
        if sample_value == '' and note_enabled:
            self.assertEqual(self.selected_batch_gui.failed_imputation_sample_text_box.placeholder, "notes...")

    def test_initial_widget_status(self):
        self.assert_widget_status(note_enabled=False, add_updated_enabled=False, delete_enabled=False)

    def test_enter_sample_already_in_input(self):
        self.selected_batch_gui.failed_imputation_sample_text_box.value = "sample_1"
        # should not be able to input scored sample, so should reset to initial widget status
        self.assert_widget_status(note_enabled=False, add_updated_enabled=False, delete_enabled=False)

    def test_enter_sample(self):
        sample = "sample_7"
        self.selected_batch_gui.failed_imputation_sample_text_box.value = sample
        self.assert_widget_status(sample_value=sample, delete_enabled=False)

    def test_enter_note(self):
        sample = "sample_7"
        note = "this is a note"
        self.selected_batch_gui.failed_imputation_sample_text_box.value = sample
        self.selected_batch_gui.failed_imputation_notes_text_box.value = note
        self.assert_widget_status(sample_value=sample, note_value=note, delete_enabled=False)

    def test_add_failed_sample(self):
        sample = "sample_7"
        note = "this is a note"
        self.add_failed_imputation_sample(sample, note)
        self.assert_widget_status(note_enabled=False, add_updated_enabled=False, delete_enabled=False)
        self.assert_correct_failed_samples_added({sample: note})

    def test_reenter_failed_sample(self):
        sample = "sample_7"
        note = "this is a note"
        self.add_failed_imputation_sample(sample, note)
        self.assert_widget_status(note_enabled=False, add_updated_enabled=False, delete_enabled=False)
        self.selected_batch_gui.failed_imputation_sample_text_box.value = sample
        self.assert_widget_status(sample_value=sample, note_value=note)

    def test_update_failed_sample(self):
        sample = "sample_7"
        note = "this is a note"
        self.add_failed_imputation_sample(sample)
        self.assert_correct_failed_samples_added({sample: ''})
        self.add_failed_imputation_sample(sample, note)
        self.assert_correct_failed_samples_added({sample: note})

    def test_remove_failed_sample(self):
        sample = "sample_7"
        self.add_failed_imputation_sample(sample)
        self.remove_failed_sample(sample)
        self.assert_widget_status(note_enabled=False, add_updated_enabled=False, delete_enabled=False)
        self.assert_correct_failed_samples_added()

    def test_sequence(self):
        random.seed(12345)
        added_failures = {}
        characters = string.ascii_letters + ' ' + string.punctuation
        for i in range(100):
            sample = f'sample_{random.randint(7, 17)}'
            if sample in added_failures and random.random() > 0.5:
                # remove
                self.remove_failed_sample(sample)
                del added_failures[sample]
            else:
                # add/update
                nchar = random.randint(0, 6)
                reason = ''.join(random.choice(characters) for _ in range(nchar))
                self.add_failed_imputation_sample(sample, reason)
                added_failures[sample] = reason
            self.assert_correct_failed_samples_added(added_failures)

    def test_finish_button(self):
        with patch.object(self.selected_batch_gui.finished_addition_imputation_failures_button, 'close') as \
                finished_addition_imputation_failures_close_mock, \
                patch.object(self.selected_batch_gui.failed_imputation_hbox, 'close') as \
                        failed_imputation_hbox_close_mock, \
                patch.object(self.selected_batch_gui, 'build_fail_samples_all_conditions_section') as \
                        build_fail_samples_all_conditions_section_mock, \
                patch.object(self.selected_batch_gui, 'print_failed_samples_all_conditions') as \
                        print_failed_samples_all_conditions_mock:
            self.selected_batch_gui.finished_addition_imputation_failures_button.click()

        finished_addition_imputation_failures_close_mock.assert_called()
        failed_imputation_hbox_close_mock.assert_called()
        build_fail_samples_all_conditions_section_mock.assert_called()
        print_failed_samples_all_conditions_mock.assert_called()


if __name__ == '__main__':
    unittest.main()
