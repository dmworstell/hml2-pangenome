"""Focused regression checks for the actual functional-refit BH consumer."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

path = Path(__file__).resolve().parents[1] / "build_artifact_filtered_functional_refit.py"
spec = importlib.util.spec_from_file_location("refit", path)
refit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(refit)


class TestBH(unittest.TestCase):
    def test_missing_is_not_zero_or_one(self):
        p = pd.Series([0.01, np.nan, 0.04, 0.03], index=[10, 11, 12, 13])
        q = refit.bh(p)
        np.testing.assert_allclose(q.loc[[10, 12, 13]], [0.03, 0.04, 0.04])
        self.assertTrue(np.isnan(q.loc[11]))

    def test_empty_and_all_missing(self):
        self.assertTrue(refit.bh(pd.Series(dtype=float)).empty)
        self.assertTrue(refit.bh(pd.Series([np.nan, np.nan])).isna().all())

    def test_ties_retained(self):
        np.testing.assert_allclose(refit.bh(pd.Series([0.01, 0.01, 0.5])), [0.015, 0.015, 0.5])

    def test_invalid_values_rejected(self):
        for p in [-0.01, 1.01, np.inf]:
            with self.assertRaises(ValueError):
                refit.bh(pd.Series([p]))


class TestMissingArrayDosage(unittest.TestCase):
    def test_unknown_haplotype_is_missing_not_zero(self):
        samples = ["A", "B", "C"]
        rows = [dict(ID=sample, Haplotype=hap, Locus=locus,
                     observation_state="PRESENT", Structure="Solo-LTR")
                for sample in samples for hap in ["mat", "pat"]
                for locus in ["HML-2_4p16.3a", "HML-2_15q25.2"]]
        oneq = pd.DataFrame({"sample": samples, "gag_compatible_dosage": [0, 0, 0]})
        seven = pd.DataFrame({"sample": ["A", "A", "B", "B", "C", "C"],
                              "haplotype": ["mat", "pat"] * 3,
                              "array_copy_number": [0, 2, 1, np.nan, 2, 3]})
        with patch.object(refit.pd, "read_csv", side_effect=[oneq, seven]):
            result = refit.derive_locus_exposures(rows, samples).set_index("sample")
        dosage = result["HML-2_7p22.1::structural::multi_vs_single"]
        self.assertEqual(dosage["A"], 1)
        self.assertTrue(np.isnan(dosage["B"]))
        self.assertEqual(dosage["C"], 2)


if __name__ == "__main__":
    unittest.main()
