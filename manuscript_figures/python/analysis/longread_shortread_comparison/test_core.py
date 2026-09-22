#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest
import csv
import gzip
import tempfile


MODULE_PATH = Path(__file__).with_name("run_analysis.py")
SPEC = importlib.util.spec_from_file_location("lrsr", MODULE_PATH)
lrsr = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = lrsr
SPEC.loader.exec_module(lrsr)


class ComparisonCoreTests(unittest.TestCase):
    def row(self, **overrides):
        row = {
            "Haplotype": "h1",
            "Structure": "Provirus",
            "provirus_type": "type2",
            "gag": "Intact",
            "pro": "Intact",
            "pol": "Intact_FS_End",
            "env": "Intact",
        }
        row.update(overrides)
        return row

    def test_intact_fs_end_is_translatable(self):
        self.assertEqual(lrsr.product_states(self.row())[2], 1)

    def test_gag_pro_pol_definition_is_identical_by_type(self):
        type1 = self.row(provirus_type="type1", env="Nonsense")
        type2 = self.row(provirus_type="type2", env="Nonsense")
        self.assertEqual(lrsr.product_states(type1)[2], 1)
        self.assertEqual(lrsr.product_states(type2)[2], 1)

    def test_deletion_is_broken_not_unknown(self):
        self.assertEqual(lrsr.status_state("Deletion"), 0)

    def test_positive_copy_wins_over_unknown_copy(self):
        cell = lrsr.Cell()
        lrsr.add_row(cell, self.row(Haplotype="h1", gag="Undetermined"))
        lrsr.add_row(cell, self.row(Haplotype="h2", gag="Intact"))
        self.assertEqual(cell.product_state(0), 1)

    def test_negative_requires_both_haplotypes(self):
        cell = lrsr.Cell()
        lrsr.add_row(cell, self.row(Haplotype="h1", gag="Frameshift"))
        self.assertIsNone(cell.product_state(0))
        lrsr.add_row(cell, self.row(Haplotype="h2", gag="Frameshift"))
        self.assertEqual(cell.product_state(0), 0)

    def test_insertion_absent_is_unknown(self):
        cell = lrsr.Cell()
        lrsr.add_row(cell, self.row(Structure="Insertion_Absent"))
        self.assertIsNone(cell.feature_state("provirus"))

    def test_excluded_artifact_cannot_create_a_carrier(self):
        rows = [self.row(ID="HG00001", Locus="HML-2_1q22", Haplotype=h,
                         gag="Frameshift", analysis_include="1") for h in ("h1", "h2")]
        rows.append(self.row(ID="HG00001", Locus="HML-2_1q22", analysis_include="0"))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "catalog.tsv.gz"
            with gzip.open(path, "wt", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            cells, people, loci, count = lrsr.load_catalog(path)
        self.assertEqual(count, 2)
        self.assertEqual(cells[("HG00001", "HML-2_1q22")].product_state(0), 0)


if __name__ == "__main__":
    unittest.main()
