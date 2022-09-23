import pandas as pd

import ImputationPipeline.ManualQCPRS.ManualQCPRS as ManualQCPRS
import unittest
import os

current_dir = os.path.dirname(__file__)
resources_dir = os.path.join(current_dir, "resources")


class FailAllConditionsForSampleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = pd.read_csv(os.path.join(resources_dir, "test_results.tsv"), delimiter='\t'). \
            fillna("NA").set_index('sample_id')
        cls.lab_batch = "ARRAY-21550"

    def setUp(self) -> None:
        self.selected_batch_gui = ManualQCPRS.SelectedBatchModificationGui(self.results, self.lab_batch)
        self.selected_batch_gui.run()
        self.selected_batch_gui.finished_addition_imputation_failures_button.click()

    def test_correct_samples_available(self):
        self.assertEqual(set(self.results[~self.results.is_control_sample].index),
                         set(self.selected_batch_gui.sample_failure_selection.options))

    def assert_initial_widget_status(self):
        # sample selection and reason should be cleared
        self.assertEqual(self.selected_batch_gui.sample_failure_selection.value, '')
        self.assertEqual(self.selected_batch_gui.sample_failure_reason_selection.value, '')
        # sample selection should be enabled
        self.assertFalse(self.selected_batch_gui.sample_failure_selection.disabled)
        # reason, fail and unfail buttons should be disabled
        self.assertTrue(self.selected_batch_gui.sample_failure_reason_selection.disabled)
        self.assertTrue(self.selected_batch_gui.fail_sample_all_conditions_button.disabled)
        self.assertTrue(self.selected_batch_gui.unfail_sample_all_conditions_button.disabled)

    def test_widget_status_initial(self):
        self.assert_initial_widget_status()

    def test_select_sample(self):
        sample_to_fail = "205247030027_R06C01"
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail
        # sample selection and reason should be enabled
        self.assertFalse(self.selected_batch_gui.sample_failure_selection.disabled)
        self.assertFalse(self.selected_batch_gui.sample_failure_reason_selection.disabled)
        # fail and unfail buttons should be disabled
        self.assertTrue(self.selected_batch_gui.fail_sample_all_conditions_button.disabled)
        self.assertTrue(self.selected_batch_gui.unfail_sample_all_conditions_button.disabled)

    def test_select_reason(self):
        sample_to_fail = "205247030027_R06C01"
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail
        self.selected_batch_gui.sample_failure_reason_selection.value = \
            self.selected_batch_gui.sample_failure_selection.options[0]
        # sample selection, reason, fail button should be enabled
        self.assertFalse(self.selected_batch_gui.sample_failure_selection.disabled)
        self.assertFalse(self.selected_batch_gui.sample_failure_reason_selection.disabled)
        self.assertFalse(self.selected_batch_gui.fail_sample_all_conditions_button.disabled)
        # unfail button should be disabled
        self.assertTrue(self.selected_batch_gui.unfail_sample_all_conditions_button.disabled)

    def test_click_fail(self):
        sample_to_fail = "205247030027_R06C01"
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail
        self.selected_batch_gui.sample_failure_reason_selection.value = \
            self.selected_batch_gui.sample_failure_reason_selection.options[0]
        self.selected_batch_gui.fail_sample_all_conditions_button.click()
        self.assert_initial_widget_status()
        # all samples other than failed should be unchanged
        remaining_samples = list(self.results.index)
        remaining_samples.remove(sample_to_fail)
        output_results = pd.read_csv(os.path.join(resources_dir, "single_sample_failed_all_conditions.tsv"),
                                     delimiter='\t').fillna("NA").set_index('sample_id')
        pd.testing.assert_frame_equal(self.selected_batch_gui.modified_results.fillna("NA"),
                                      output_results)

    def test_widget_status_failed_sample_selected(self):
        sample_to_fail = "205247030027_R06C01"
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail
        self.selected_batch_gui.sample_failure_reason_selection.value = \
            self.selected_batch_gui.sample_failure_reason_selection.options[0]
        self.selected_batch_gui.fail_sample_all_conditions_button.click()
        self.selected_batch_gui.sample_failure_selection.value = sample_to_fail

        self.assertEqual(self.selected_batch_gui.sample_failure_selection.value, sample_to_fail)
        self.assertEqual(self.selected_batch_gui.sample_failure_reason_selection.value,
                         self.selected_batch_gui.sample_failure_reason_selection.options[0])
        #sample selection, reason, and unfail should be enabled
        self.assertFalse(self.selected_batch_gui.sample_failure_selection.disabled)
        self.assertFalse(self.selected_batch_gui.sample_failure_reason_selection.disabled)
        self.assertFalse(self.selected_batch_gui.fail_sample_all_conditions_button.disabled)
        self.assertFalse(self.selected_batch_gui.unfail_sample_all_conditions_button.disabled)


if __name__ == '__main__':
    unittest.main()
