"""Focused regression tests for the Figure 1/S1/S5 observation boundary."""

import csv
import json
import tempfile
import unittest
from pathlib import Path

from structural_catalog_summary import structural_state, structural_summary


def record(locus, sample, hap, state="PRESENT", structure="Provirus", original=None):
    return dict(Locus=locus, ID=sample, Haplotype=hap, observation_state=state,
                Structure=structure, orig_Locus=original or locus,
                haplotype_pansn_digit=hap[-1])


class StructuralObservationTests(unittest.TestCase):
    def test_missing_is_unknown_not_noncarrier(self):
        self.assertEqual(structural_state([]), "Unknown")

    def test_noncarrier_is_retained_but_not_labelled_absent(self):
        row = record("HML-2_1q22", "HG1", "h1", "NONCARRIER_SENTINEL_OBSERVATION")
        self.assertEqual(structural_state([row]), "Noncarrier call")

    def test_unknown_does_not_turn_into_absence(self):
        unknown = record("HML-2_1q22", "HG1", "h1", "UNKNOWN_TECHNICAL")
        noncarrier = record("HML-2_1q22", "HG1", "h1", "NONCARRIER_SENTINEL_OBSERVATION")
        self.assertEqual(structural_state([unknown, noncarrier]), "Unknown")

    def test_present_and_multicopy_precedence(self):
        present = record("HML-2_1q22", "HG1", "h1")
        noncarrier = record("HML-2_1q22", "HG1", "h1", "NONCARRIER_SENTINEL_OBSERVATION")
        self.assertEqual(structural_state([present, noncarrier]), "Provirus")
        self.assertEqual(structural_state([present, present.copy()]), "Multi-copy")

    def test_sex_denominators_and_unlocalized_records(self):
        roster = [(sample, hap) for sample in ("HG1", "HG2", "HG3") for hap in ("h1", "h2")]
        rows = [record("HML-2_1q22", sample, hap) for sample, hap in roster]
        rows += [record("HML-2_Xq12", "HG1", "h1"), record("HML-2_Xq12", "HG1", "h2"),
                 record("HML-2_Xq12", "HG2", "h2"), record("HML-2_Yp11.2", "HG2", "h1"),
                 record("HML-2_Yp11.2", "HG2", "h2"),  # wrong male partition is not a second Y
                 record("HML-2_acro_type1", "HG1", "h1", original="HML-2_13p13"),
                 record("HML-2_13p13", "HG2", "h1", "NONCARRIER_SENTINEL_OBSERVATION")]
        with tempfile.TemporaryDirectory() as tmp:
            frame = Path(tmp) / "frame.tsv"
            with frame.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerows([("sample_id", "sex"), ("HG1", "female"),
                                  ("HG2", "male"), ("HG3", "unknown")])
            partitions = Path(tmp) / "xy.json"
            partitions.write_text(json.dumps({"copy_assignment_contract": {"male_partition_assignments": [
                {"sample_id": "HG2", "X_partition": "h2", "Y_partition": "h1"}
            ]}}))
            summary, observations, _ = structural_summary(rows, roster, frame, partitions)
        by_locus = {row["locus"]: row for row in summary}
        self.assertEqual(by_locus["HML-2_1q22"]["eligible_haplotypes"], 6)
        self.assertEqual(by_locus["HML-2_Xq12"]["eligible_haplotypes"], 3)
        y = by_locus["HML-2_Yp11.2"]
        self.assertEqual((y["eligible_haplotypes"], y["Provirus"], y["Multi-copy"]), (1, 1, 0))
        self.assertEqual(y["unknown_ploidy_haplotypes"], 2)
        grouped = by_locus["HML-2_acro_type1"]
        self.assertEqual(grouped["label_scope"], "unlocalized_record_bucket")
        self.assertEqual(grouped["eligible_haplotypes"], 0)
        self.assertEqual(grouped["group_observed_haplotypes"], 1)
        self.assertEqual(grouped["Noncarrier call"], 0)
        self.assertEqual(by_locus["HML-2_13p13"]["missing_reassigned"], 1)
        moved = next(row for row in observations if row["locus"] == "HML-2_13p13"
                     and row["sample"] == "HG1" and row["haplotype"] == "h1")
        self.assertEqual(moved["state"], "Unknown")


if __name__ == "__main__":
    unittest.main()
