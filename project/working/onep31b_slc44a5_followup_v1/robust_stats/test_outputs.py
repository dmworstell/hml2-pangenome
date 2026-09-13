#!/usr/bin/env python3
import csv, json, pathlib, unittest

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "results"

def rows(name):
    with (OUT / name).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))

class RobustStatsTests(unittest.TestCase):
    def test_semantics_and_alignment(self):
        s = json.loads((OUT / "summary.json").read_text())
        self.assertEqual((s["n"], s["carriers"], s["solo_ltr_reference_people"]), (39, 13, 26))
        self.assertEqual(s["internal_dose_counts"], {"0": 26, "1": 12, "2": 1})
        self.assertEqual(s["array_carriers"], 0)
        self.assertEqual(s["corrected_identity_swaps_in_overlap"], 0)
        self.assertIn("solo-LTR", s["contrast"])
        self.assertIn("not insertion", s["contrast"])

    def test_primary_reproduces(self):
        p = next(x for x in rows("model_sensitivities.tsv") if x["model_id"] == "primary_reproduction")
        self.assertAlmostEqual(float(p["beta"]), 1.80632008094487, places=9)
        self.assertAlmostEqual(float(p["conventional_p"]), 9.82893291496621e-08, delta=1e-14)
        self.assertLess(float(p["hc3_p"]), 0.001)

    def test_robustness(self):
        m = rows("model_sensitivities.tsv")
        self.assertTrue(all(float(x["beta"]) > 0 for x in m))
        loo = rows("leave_one_out.tsv")
        self.assertEqual(len(loo), 39)
        self.assertTrue(all(float(x["beta"]) > 0 for x in loo))
        perm = rows("blocked_permutation.tsv")[0]
        self.assertEqual(int(perm["exceedances"]), 0)
        self.assertLessEqual(float(perm["corrected_p"]), 1/20000)

    def test_ancestry_and_pedigree_outputs(self):
        a = rows("ancestry_strata.tsv")
        self.assertEqual({x["superpopulation"] for x in a}, {"AFR","AMR","EAS","EUR","SAS"})
        self.assertEqual(sum(int(x["carriers"]) for x in a), 13)
        self.assertTrue((OUT / "recorded_related_pairs.tsv").exists())

if __name__ == "__main__": unittest.main(verbosity=2)
