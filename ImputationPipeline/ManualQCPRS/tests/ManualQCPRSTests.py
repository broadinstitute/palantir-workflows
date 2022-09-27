import pandas as pd
import ImputationPipeline.ManualQCPRS.ManualQCPRS as ManualQCPRS
import os
import itertools
import unittest
from unittest.mock import patch

current_dir = os.path.dirname(__file__)
resources_dir = os.path.join(current_dir, "resources")


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
            columns=expected_unfailed_results_df.columns.tolist() + ['poly_test_not_performed_reason', 'notes']).fillna(
            "NA")
        pd.testing.assert_frame_equal(
            self.selected_batch_gui.modified_results.drop(failed_condition_columns, axis=1).copy().fillna("NA"),
            expected_unfailed_results_df)

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
        pd.testing.assert_frame_equal(expected_results, failed_conditions_results, check_like=True)

    def assert_widget_status(self, condition_value=None, reason_value='',
                             condition_enabled=True, reason_enabled=True, fail_enabled=True, unfail_enabled=True):
        self.assertEqual(self.selected_batch_gui.condition_failure_selection.value, condition_value)
        self.assertEqual(self.selected_batch_gui.condition_failure_reason_selection.value, reason_value)
        self.assertEqual(self.selected_batch_gui.condition_failure_selection.disabled, not condition_enabled)
        self.assertEqual(self.selected_batch_gui.condition_failure_reason_selection.disabled, not reason_enabled)
        self.assertEqual(self.selected_batch_gui.fail_condition_all_samples_button.disabled, not fail_enabled)
        self.assertEqual(self.selected_batch_gui.unfail_condition_all_samples_button.disabled, not unfail_enabled)

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
                patch.object(self.selected_batch_gui.condition_failure_hbox,'close') as \
                condition_failure_hbox_close_mock, \
                patch.object(self.selected_batch_gui, 'build_fail_sample_for_condition_section') as \
                build_fail_sample_for_condition_section_mock, \
                patch.object(self.selected_batch_gui, 'printManualFailureStatus') as \
                printManualFailureStatus_mock:
            self.selected_batch_gui.finished_condition_failures_button.click()

        finished_condition_failures_button_close_mock.assert_called()
        condition_failure_hbox_close_mock.assert_called()
        build_fail_sample_for_condition_section_mock.assert_called()
        printManualFailureStatus_mock.assert_called()


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
        for i in range(1, 4):
            expected_results[f'condition_{i}_reason_not_resulted'] = \
                expected_results.index.map(failed_samples_with_reasons)

        pd.testing.assert_frame_equal(expected_results, failed_samples_results, check_like=True)

    def assert_modified_as_expected(self, failed_samples=None):
        if failed_samples is None:
            failed_samples = {}
        # first assert that all entries for unfailed samples remain the same
        expected_unfailed_samples_df = self.results.query("sample_id not in @failed_samples.keys()").copy(). \
            reindex(columns=self.results.columns.tolist() + ['poly_test_not_performed_reason', 'notes']).fillna("NA")
        pd.testing.assert_frame_equal(
            self.selected_batch_gui.modified_results.query("sample_id not in @failed_samples.keys()").copy().fillna(
                "NA"),
            expected_unfailed_samples_df, check_like=True)

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
        self.assertEqual(set(self.results[~self.results.is_control_sample].index),
                         set(self.selected_batch_gui.sample_failure_selection.options))

    def assert_widget_status(self, sample_value='', reason_value='',
                             sample_enabled=True, reason_enabled=True, fail_enabled=True, unfail_enabled=True):
        self.assertEqual(self.selected_batch_gui.sample_failure_selection.value, sample_value)
        self.assertEqual(self.selected_batch_gui.sample_failure_reason_selection.value, reason_value)
        self.assertEqual(self.selected_batch_gui.sample_failure_selection.disabled, not sample_enabled)
        self.assertEqual(self.selected_batch_gui.sample_failure_reason_selection.disabled, not reason_enabled)
        self.assertEqual(self.selected_batch_gui.fail_sample_all_conditions_button.disabled, not fail_enabled)
        self.assertEqual(self.selected_batch_gui.unfail_sample_all_conditions_button.disabled, not unfail_enabled)

    def assert_output_text_correct(self, failed_samples_with_reasons: dict = None):
        if not failed_samples_with_reasons:
            expected_text = "No Samples Will Be Failed For All Conditions"
            self.assertEqual(self.selected_batch_gui.out_failed_sample_all_conditions.outputs[0]["text"],
                             expected_text)
        else:
            output_text = self.selected_batch_gui.out_failed_sample_all_conditions.outputs.output[0]["text"].strip().split("\n")
            output_first_line = output_text[0]
            expected_first_line = 'The following samples will be marked as Failed PRS and failed for all conditions'
            self.assertEqual(output_first_line, expected_first_line)
            # not concerned with order of sample failure text
            sample_failure_text = set(output_text[1:])
            expected_failure_text = {f'{sample}, notes: {reason}' for sample, reason in failed_samples_with_reasons.items()}
            self.assertEqual(sample_failure_text, expected_failure_text)

    def test_initial_widget_status(self):
        self.assert_widget_status(reason_enabled=False, fail_enabled=False, unfail_enabled=False)
        self.assert_output_text_correct()

    def test_select_sample(self):
        sample_to_fail = "sample_2"
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail
        self.assert_widget_status(sample_value=sample_to_fail, fail_enabled=False, unfail_enabled=False)
        self.assert_output_text_correct()

    def test_select_reason(self):
        sample_to_fail = "sample_2"
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail
        self.selected_batch_gui.sample_failure_reason_selection.value = self.default_reason
        self.assert_widget_status(sample_value=sample_to_fail, unfail_enabled=False,
                                  reason_value=self.default_reason)
        self.assert_output_text_correct()

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
                patch.object(self.selected_batch_gui, 'printFailedConditionsAllSamples') as \
                printFailedConditionsAllSamples_mock:
            self.selected_batch_gui.finished_sample_failures_button.click()

        finished_sample_failures_button_close_mock.assert_called()
        sample_failure_hbox_close_mock.assert_called()
        build_fail_condition_section_mock.assert_called()
        printFailedConditionsAllSamples_mock.assert_called()


if __name__ == '__main__':
    unittest.main()
