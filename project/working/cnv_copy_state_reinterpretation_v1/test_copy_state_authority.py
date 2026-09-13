#!/usr/bin/env python3
import math
import unittest

from copy_state_authority import (
    STATE_ARTIFACT,
    STATE_TRUE,
    calibrate_positive_controls,
    region_metrics,
    score_candidate,
)


class CopyStateAuthorityTests(unittest.TestCase):
    def setUp(self):
        controls = [
            ("c1", {"body_median_all": 20, "flank_median_all": 40, "sample_ref_cov": 40}),
            ("c2", {"body_median_all": 18, "flank_median_all": 38, "sample_ref_cov": 40}),
            ("c3", {"body_median_all": 22, "flank_median_all": 36, "sample_ref_cov": 40}),
            ("c4", {"body_median_all": 17, "flank_median_all": 42, "sample_ref_cov": 40}),
        ]
        self.calibration = calibrate_positive_controls(
            controls, body_crossmap_ambiguous=False
        )

    def test_hg00423_metrics_use_unique_flank_only(self):
        row = {"body_median_all": 6, "flank_median_all": 1, "sample_ref_cov": 41}
        metrics = region_metrics(row, body_crossmap_ambiguous=True)
        self.assertAlmostEqual(metrics["expected_haploid_body_coverage"], 20.5)
        self.assertAlmostEqual(metrics["body_ratio_to_one_copy"], 6 / 20.5)
        self.assertAlmostEqual(metrics["local_flank_ratio_to_represented_contig"], 1 / 41)
        self.assertAlmostEqual(metrics["region_support_ratio"], 1 / 41)
        self.assertEqual(metrics["likelihood_features"], ["local_flank_ratio"])

    def test_hg00423_artifact_is_likelihood_favored(self):
        result = score_candidate(
            candidate_key="HG00423|pat|1q22|alt2",
            row={"body_median_all": 6, "flank_median_all": 1, "sample_ref_cov": 41},
            calibration=self.calibration,
            always_biological_tips=["alt1"],
            candidate_tip="alt2",
            terminal_host="HG00423",
            body_crossmap_ambiguous=True,
        )
        self.assertEqual(result["assay_use_count"], 1)
        self.assertEqual(result["central_likelihood_favors"], STATE_ARTIFACT)
        self.assertGreater(result["central_state_likelihoods"][STATE_ARTIFACT], 0.5)
        self.assertEqual(result["state_contract"][STATE_ARTIFACT]["biological_tips"], ["alt1"])
        self.assertEqual(result["state_contract"][STATE_TRUE]["biological_tips"], ["alt1", "alt2"])
        self.assertEqual(result["state_contract"][STATE_TRUE]["terminal_hosts"], ["HG00423"])

    def test_likelihoods_normalize_at_every_scale(self):
        result = score_candidate(
            candidate_key="x", row={"body_median_all": 10, "flank_median_all": 20, "sample_ref_cov": 40},
            calibration=self.calibration, always_biological_tips=["base"], candidate_tip="candidate",
            terminal_host="person", body_crossmap_ambiguous=False,
        )
        for weights in result["state_likelihoods"].values():
            self.assertAlmostEqual(sum(weights.values()), 1.0)
            self.assertTrue(all(math.isfinite(value) for value in weights.values()))

    def test_candidate_cannot_be_reused_as_control(self):
        with self.assertRaisesRegex(ValueError, "must not also"):
            score_candidate(
                candidate_key="c1", row={"body_median_all": 20, "flank_median_all": 40, "sample_ref_cov": 40},
                calibration=self.calibration, always_biological_tips=["base"], candidate_tip="candidate",
                terminal_host="person", body_crossmap_ambiguous=False,
            )

    def test_crossmap_body_cannot_rescue_missing_flank(self):
        metrics = region_metrics(
            {"body_median_all": 500, "flank_median_all": 0, "sample_ref_cov": 40},
            body_crossmap_ambiguous=True,
        )
        self.assertEqual(metrics["region_support_ratio"], 0)

    def test_degenerate_calibration_fails_closed(self):
        controls = [(str(i), {"body_median_all": 20, "flank_median_all": 40, "sample_ref_cov": 40}) for i in range(3)]
        with self.assertRaisesRegex(ValueError, "unidentified"):
            calibrate_positive_controls(
                controls, body_crossmap_ambiguous=False
            )


if __name__ == "__main__":
    unittest.main()
