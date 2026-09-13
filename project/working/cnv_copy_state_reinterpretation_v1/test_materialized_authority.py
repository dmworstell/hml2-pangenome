#!/usr/bin/env python3
import csv
import hashlib
import json
import unittest
from collections import Counter
from pathlib import Path

from copy_state_authority import STATE_ARTIFACT, STATE_TRUE


LANE = Path(__file__).resolve().parent
RESULTS = LANE / "results"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MaterializedAuthorityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.authority = json.loads(
            (RESULTS / "current_unmeasured_copy_state_authority.v1.json").read_text()
        )
        cls.hg = json.loads(
            (RESULTS / "hg00423_1q22_copy_state_authority.v1.json").read_text()
        )
        cls.calibration = json.loads(
            (RESULTS / "real_control_calibration.v1.json").read_text()
        )
        cls.receipt = json.loads((RESULTS / "RECEIPT.v1.json").read_text())

    def test_exact_historical_and_current_universes(self):
        self.assertEqual(self.authority["source_universes"]["historical_unmeasured_rows"], 68)
        self.assertEqual(self.authority["source_universes"]["current_unmeasured_rows"], 44)
        self.assertEqual(len(self.authority["records"]), 44)
        with (RESULTS / "historical_unmeasured_inventory.v1.tsv").open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 68)

    def test_current_class_partition(self):
        observed = Counter(record["authority_class"] for record in self.authority["records"])
        self.assertEqual(observed, Counter({
            "AUTHENTICATED_SEGMENTAL_DUPLICATION_FIXED": 9,
            "NON_SEGDUP_DUPLICATION_LATENT_CNV_WEIGHTED": 26,
            "UNINFORMATIVE_NO_LOCAL_FLANK_MEASUREMENT": 9,
        }))

    def test_real_control_calibration(self):
        calibration = self.calibration["calibration"]
        self.assertEqual(calibration["n_positive_controls"], 61)
        self.assertTrue(self.calibration["control_authority"]["candidate_rows_excluded"])
        self.assertIsNone(calibration["prior"])
        scales = calibration["scales"]
        self.assertGreater(scales["narrow_q75_lower_residual"], 0)
        self.assertGreater(scales["central_rmse_lower_residual"], 0)
        self.assertGreater(scales["wide_q90_lower_residual"], 0)

    def test_hg00423_exact_state_contract_and_weights(self):
        self.assertEqual(self.hg["metrics"]["region_support_ratio"], 1 / 41)
        self.assertEqual(self.hg["assay_use_count"], 1)
        self.assertEqual(self.hg["assay_use_count_in_copy_state_weight"], 1)
        self.assertEqual(self.hg["central_likelihood_favors"], STATE_ARTIFACT)
        for weights in self.hg["state_likelihoods"].values():
            self.assertAlmostEqual(sum(weights.values()), 1.0)
            self.assertGreater(weights[STATE_ARTIFACT], weights[STATE_TRUE])
            self.assertGreater(weights[STATE_TRUE], 0)
        artifact = self.hg["state_contract"][STATE_ARTIFACT]
        true = self.hg["state_contract"][STATE_TRUE]
        self.assertEqual(len(artifact["biological_tips"]), 1)
        self.assertEqual(len(true["biological_tips"]), 2)
        self.assertEqual(artifact["terminal_hosts"], ["HG00423"])
        self.assertEqual(true["terminal_hosts"], ["HG00423"])

    def test_receipt_binds_every_result(self):
        for name, expected in self.receipt["result_sha256"].items():
            self.assertEqual(sha256(RESULTS / name), expected)
        self.assertFalse(self.receipt["rerun_or_submission_performed"])


if __name__ == "__main__":
    unittest.main()
