"""Check the public Figure 1/S5 numeric fixtures without the external catalog."""

import csv
import unittest
from collections import Counter
from pathlib import Path

from structural_catalog_summary import STATES

ROOT = Path(__file__).resolve().parents[2]


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class RetainedSourceTests(unittest.TestCase):
    def test_figure_and_s5_are_identical_observations(self):
        figure = read(ROOT / "data/Figure_1_structural_source_data.tsv")
        table = read(ROOT / "Supplementary_Data/Table_S5_artifact_filtered_structural_spectrum.tsv")
        self.assertEqual(len(figure), 105)
        self.assertEqual([{key: value for key, value in row.items() if key != "shown_in_figure_1"}
                          for row in figure], table)
        selected = sorted((row for row in figure if row["label_scope"] == "physical_locus"),
                          key=lambda row: float(row["variability_score"]), reverse=True)[:24]
        self.assertEqual({row["locus"] for row in selected},
                         {row["locus"] for row in figure if row["shown_in_figure_1"] == "1"})

    def test_all_denominators_conserve_chromosome_units(self):
        summary = read(ROOT / "data/Figure_1_structural_source_data.tsv")
        scopes = Counter(row["label_scope"] for row in summary)
        self.assertEqual(scopes, {"physical_locus": 101, "unlocalized_record_bucket": 4})
        for row in summary:
            denominator = int(row["eligible_haplotypes"])
            if row["label_scope"] == "unlocalized_record_bucket":
                self.assertEqual(denominator, 0)
                self.assertEqual(row["Noncarrier call"], "0")
                continue
            self.assertEqual(sum(int(row[state]) for state in STATES), denominator)
            self.assertEqual(denominator + int(row["not_applicable_haplotypes"])
                             + int(row["unknown_ploidy_haplotypes"]), 584)
            self.assertEqual(denominator, 432 if row["locus"].startswith("HML-2_X") else
                             144 if row["locus"].startswith("HML-2_Y") else 584)

    def test_retained_sex_partition_counts(self):
        evidence = read(ROOT / "data/Figure_1_sex_chromosome_eligibility.tsv")
        counts = Counter((row["chromosome"], row["eligibility"]) for row in evidence)
        self.assertEqual(counts, {("X", "eligible"): 432, ("X", "not_applicable"): 144,
                                  ("X", "unknown_ploidy"): 8, ("Y", "eligible"): 144,
                                  ("Y", "not_applicable"): 432, ("Y", "unknown_ploidy"): 8})
        samples = {row["sample"]: row["sex"] for row in evidence}
        self.assertEqual(Counter(samples.values()), {"female": 144, "male": 144, "unknown": 4})


if __name__ == "__main__":
    unittest.main()
