from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import unittest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "build_final_assayed_weights", HERE / "build_final_assayed_weights.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FinalAssayedWeightTest(unittest.TestCase):
    def test_depth_law_is_normalized_and_directional(self) -> None:
        low = MODULE.normalized_depth_law(1.0 / 41.0, 0.1653252220857001)
        self.assertAlmostEqual(sum(low.values()), 1.0)
        self.assertGreater(low[MODULE.ARTIFACT], 0.999999)
        high = MODULE.normalized_depth_law(1.0, 0.1653252220857001)
        self.assertGreater(high[MODULE.LATER], high[MODULE.ARTIFACT])
        self.assertEqual(high[MODULE.SEGDUP], 0.0)
        self.assertEqual(high[MODULE.ARRAY], 0.0)

    def test_hg00423_matches_accepted_two_state_depth_likelihood(self) -> None:
        law = MODULE.normalized_depth_law(1.0 / 41.0, 0.1653252220857001)
        self.assertTrue(
            math.isclose(
                law[MODULE.ARTIFACT],
                0.9999999722750369,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
        )
        self.assertTrue(
            math.isclose(
                law[MODULE.LATER],
                2.772496301734011e-08,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
        )

    def test_site_replacement_recomputes_copy_fraction(self) -> None:
        entry = {
            "site": {
                "copy_accounting": [
                    {"candidate_key": "base", "body_median_all": 63, "flank_median_all": 58},
                    {"candidate_key": "extra", "body_median_all": None, "flank_median_all": None},
                ]
            }
        }
        metrics = MODULE.site_metrics("extra", 6.0, 1.0, entry)
        self.assertEqual(metrics["site_measurement_completeness"], "COMPLETE")
        self.assertEqual(metrics["site_summed_flank_medians"], 59.0)
        self.assertAlmostEqual(metrics["candidate_flank_fraction_of_site"], 1.0 / 59.0)

    def test_site_index_accepts_receipt_declared_incremental_universe(self) -> None:
        payload = {
            "site_count": 1,
            "placement_count": 2,
            "sites": [{
                "site_key": "sample|new-locus",
                "copy_accounting": [
                    {"candidate_key": "sample|new-locus|region1"},
                    {"candidate_key": "sample|new-locus|region2"},
                ],
            }],
        }
        result = MODULE.site_index(payload)
        self.assertEqual(set(result), {
            "sample|new-locus|region1",
            "sample|new-locus|region2",
        })

    def test_boundary_measurements_uses_declared_target_set(self) -> None:
        rows = [{
            "candidate_id": "sample|new-locus|region2",
            "current_resolution": "INFORMATIVE_EXACT_SIDE_DEPTH",
            "raw_body_median": "21",
            "raw_local_flank_median": "19",
            "raw_numeric_override": "false",
        }, {
            "candidate_id": "sample|unrelated|region2",
            "current_resolution": "INFORMATIVE_EXACT_SIDE_DEPTH",
            "raw_body_median": "99",
            "raw_local_flank_median": "99",
            "raw_numeric_override": "false",
        }]
        result = MODULE.boundary_measurements(
            rows, {"sample|new-locus|region2"}
        )
        self.assertEqual(set(result), {"sample|new-locus|region2"})
        self.assertEqual(result["sample|new-locus|region2"]["flank"], 19.0)


if __name__ == "__main__":
    unittest.main()
