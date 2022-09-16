import ImputationPipeline.CreateSampleSets.CreateSampleSets as CreateSampleSets
import unittest
import os
import firecloud.api as fapi
import json
import random


def json_lists_to_sets(o):
    if isinstance(o, list):
        return {json_lists_to_sets(v) for v in o}
    elif isinstance(o, dict):
        return {k: json_lists_to_sets(v) for k, v in o.items()}
    return o


def build_sample_sets_lite(sample_sets: list):
    sample_sets_lite = dict()
    for sample_set in sample_sets:
        name = sample_set["name"]
        group = sample_set["attributes"]["group"]
        lab_batch = sample_set["attributes"]["lab_batch"]
        samples = {sample["entityName"] for sample in sample_set["attributes"]["samples"]["items"]}
        delivered = sample_set["attributes"]["delivered"]
        redeliver = sample_set["attributes"]["redeliver"]
        sample_sets_lite[name] = {"group": group, "lab_batch": lab_batch, "samples": samples, "delivered": delivered,
                                  "redeliver": redeliver}
    return sample_sets_lite


class IntegrationTests(unittest.TestCase):
    workspace_namespace = None
    workspace_name = None
    maxDiff = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.workspace_namespace = os.getenv("WORKSPACE_NAMESPACE")
        if cls.workspace_namespace is None:
            raise (RuntimeError("Environment variable WORKSPACE_NAMESPACE must be set to run tests"))
        workspace_hash = random.getrandbits(64)
        cls.workspace_name = f'test_workspace_{workspace_hash:016x}'
        print(f"Creating test workspace {cls.workspace_namespace}/{cls.workspace_name}")
        r = fapi.create_workspace(cls.workspace_namespace, cls.workspace_name)
        fapi._check_response_code(r, 201)

    @classmethod
    def tearDownClass(cls) -> None:
        print(f"Deleting test workspace {cls.workspace_namespace}/{cls.workspace_name}")
        fapi.delete_workspace(cls.workspace_namespace, cls.workspace_name)

    def tearDown(self) -> None:
        fapi.delete_entities_of_type(self.workspace_namespace, self.workspace_name, "sample_set")
        fapi.delete_entities_of_type(self.workspace_namespace, self.workspace_name, "sample")
        pass

    def run_test(self, samples_tsv=None, expected_sample_sets_json=None, sample_sets_to_mark_delivered=None,
                 samples_to_mark_rework=None):
        if sample_sets_to_mark_delivered is None:
            sample_sets_to_mark_delivered = []
        if samples_to_mark_rework is None:
            samples_to_mark_rework = []
        if samples_tsv is not None:
            fapi.upload_entities_tsv(self.workspace_namespace, self.workspace_name, samples_tsv, "flexible")
        for sample_set in sample_sets_to_mark_delivered:
            fapi.update_entity(self.workspace_namespace, self.workspace_name, "sample_set", sample_set,
                               [{"op": "AddUpdateAttribute",
                                 "attributeName": "delivered",
                                 "addUpdateAttribute": True}
                                ]
                               )
        for sample in samples_to_mark_rework:
            fapi.update_entity(self.workspace_namespace, self.workspace_name, "sample", sample,
                               [{"op": "AddUpdateAttribute",
                                 "attributeName": "rework",
                                 "addUpdateAttribute": True}
                                ]
                               )

        CreateSampleSets.main(self.workspace_namespace, self.workspace_name)
        if expected_sample_sets_json is not None:
            self.assert_expected_sample_sets(expected_sample_sets_json)
        self.assert_all_reworks_false()

    def test_simple(self):
        samples_seq = [("tests/resources/samples.tsv", "tests/resources/expected_sample_sets.json"),
                       ("tests/resources/samples2.tsv", "tests/resources/expected_sample_sets_2.json"),
                       ("tests/resources/samples3.tsv", "tests/resources/expected_sample_sets_3.json"),
                       ("tests/resources/samples4.tsv", "tests/resources/expected_sample_sets_4.json")
                       ]

        for samples_tsv, expected_sample_sets_json in samples_seq:
            with self.subTest(samples=samples_tsv, expected_samples=expected_sample_sets_json):
                self.run_test(samples_tsv, expected_sample_sets_json)

    def test_sample_sets_delivered(self):
        samples_seq = [("tests/resources/samples.tsv", "tests/resources/expected_sample_sets.json", []),
                       ("tests/resources/samples2.tsv", "tests/resources/expected_sample_sets_2_with_delivery.json", ["lb_2"]),
                       ("tests/resources/samples3.tsv", "tests/resources/expected_sample_sets_3_with_delivery.json", ["lb_1"]),
                       ("tests/resources/samples4.tsv", "tests/resources/expected_sample_sets_4_with_delivery.json", ["lb_2_group_2"])]

        for samples_tsv, expected_sample_sets_json, sample_sets_to_mark_delivered in samples_seq:
            with self.subTest(samples=samples_tsv, expected_samples=expected_sample_sets_json,
                              sample_set_to_mark_delivered=sample_sets_to_mark_delivered):
                self.run_test(samples_tsv, expected_sample_sets_json, sample_sets_to_mark_delivered)

    def test_lone_control(self):
        self.run_test("tests/resources/samples_lone_control.tsv", "tests/resources/expected_sample_sets_lb_2_only.json")

    def test_no_control(self):
        self.run_test("tests/resources/samples_no_control.tsv", "tests/resources/expected_sample_sets_lb_2_only.json")

    def test_multiple_controls(self):
        with self.assertRaisesRegex(RuntimeError, "Multiple control samples"):
            self.run_test("tests/resources/samples_multiple_controls.tsv")

    def test_multiple_undelivered(self):
        self.run_test("tests/resources/samples.tsv")
        # deliver lb_1, add new samples
        self.run_test("tests/resources/samples2.tsv", sample_sets_to_mark_delivered=["lb_1"])
        # unmark delivery of lb_1
        fapi.update_entity(self.workspace_namespace, self.workspace_name, "tests/resources/sample_set", "lb_1",
                           [{"op": "AddUpdateAttribute",
                             "attributeName": "delivered",
                             "addUpdateAttribute": False}
                            ]
                           )
        with self.assertRaisesRegex(RuntimeError, "has not been delivered"):
            self.run_test()

    def test_rework(self):
        samples_seq = [("tests/resources/samples.tsv", "tests/resources/expected_sample_sets.json", [], ["sample_id_1"]),
                       # sample_id_1 marked as rework, but not yet in sample set.
                       # should ignore rework and treat normally
                       (None, "tests/resources/expected_sample_sets.json", [], ["sample_id_1"]),
                       # mark sample_id_1 as reworked, but don't yet deliver, result is same
                       (None, "tests/resources/expected_sample_sets_rework.json", ["lb_1"], ["sample_id_1"]),
                       # deliver lb_1, now reworking sample_id_1 leads to additional aggregation set
                       ]
        for samples_tsv, expected_sample_sets_json, sample_sets_to_mark_delivered, samples_to_mark_rework in samples_seq:
            with self.subTest(samples=samples_tsv, expected_samples=expected_sample_sets_json,
                              sample_set_to_mark_delivered=sample_sets_to_mark_delivered,
                              samples_to_mark_rework=samples_to_mark_rework):
                self.run_test(samples_tsv, expected_sample_sets_json, sample_sets_to_mark_delivered,
                              samples_to_mark_rework)

    def test_rework_no_control(self):
        # the fact that the sample marked as rework should have no effect when there is no control sample
        self.run_test("tests/resources/samples_no_control.tsv", "tests/resources/expected_sample_sets_lb_2_only.json", samples_to_mark_rework="sample_id_1")

    def test_rework_and_new_mixed(self):
        samples_seq = [("tests/resources/samples.tsv", "tests/resources/expected_sample_sets.json", [], []),
                       ("tests/resources/samples2.tsv", "tests/resources/expected_sample_sets_2_mix_rework_and_new.json", ["lb_2"], ["sample_id_5"]) # deliver lb_2, rework sample 5
                       ]
        for samples_tsv, expected_sample_sets_json, sample_sets_to_mark_delivered, samples_to_mark_rework in samples_seq:
            with self.subTest(samples=samples_tsv, expected_samples=expected_sample_sets_json,
                              sample_set_to_mark_delivered=sample_sets_to_mark_delivered,
                              samples_to_mark_rework=samples_to_mark_rework):
                self.run_test(samples_tsv, expected_sample_sets_json, sample_sets_to_mark_delivered,
                              samples_to_mark_rework)

    def assert_expected_sample_sets(self, expected_sample_sets_json):
        with open(expected_sample_sets_json) as expected_sample_set:
            expected_sample_sets_lists = json.load(expected_sample_set)
            expected_sample_sets_lite = json_lists_to_sets(expected_sample_sets_lists)
        observed_sample_sets_lite = self.get_sample_sets_lite()
        self.assertEqual(expected_sample_sets_lite, observed_sample_sets_lite)

    def get_sample_sets_lite(self):
        sample_set_response = fapi.get_entities(self.workspace_namespace, self.workspace_name, 'sample_set')
        if not sample_set_response.ok:
            raise RuntimeError(f'ERROR: {sample_set_response.text}')
        return build_sample_sets_lite(json.loads(sample_set_response.text))

    def assert_all_reworks_false(self):
        samples = self.get_samples()
        for sample in samples:
            rework = sample["attributes"].get("rework", False)
            self.assertFalse(rework)

    def get_samples(self):
        sample_response = fapi.get_entities(self.workspace_namespace, self.workspace_name, 'sample')
        if not sample_response.ok:
            raise RuntimeError(f'ERROR: {sample_response.text}')
        return json.loads(sample_response.text)


if __name__ == '__main__':
    unittest.main()
