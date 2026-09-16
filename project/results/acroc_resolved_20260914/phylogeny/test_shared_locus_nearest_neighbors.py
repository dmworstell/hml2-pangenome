"""Actual-input regression for the focused tree comparison."""

import contextlib
import ast
import csv
import importlib.util
import io
import json
import shutil
import re
import tempfile
import types
import unittest
from pathlib import Path

from Bio import Phylo


SOURCE = Path(__file__).resolve().parent


class SharedLocusComparisonTest(unittest.TestCase):
    def test_restored_figure_selection_keeps_existing_tips(self):
        manuscript = SOURCE.parents[2] / "manuscript"
        owner_source = ast.parse((manuscript / "build_narrative_main_figures.py").read_text())
        restore_source = ast.parse((manuscript / "restore_figure3_20260914.py").read_text())
        representative = next(node for node in owner_source.body if isinstance(node, ast.FunctionDef) and node.name == "representative_tip_names")
        selection = next(node for node in restore_source.body if isinstance(node, ast.FunctionDef) and node.name == "select_tree_tips")
        trees = {
            region: Phylo.read(SOURCE / name, "newick")
            for region, name in (("LTR", "hml2_pan_ltr_expanded_tree.nwk"), ("Pol", "hml2_pan_orf_pol_tree.nwk"))
        }
        shared = set.intersection(*({tip.name.split("__")[0] for tip in tree.get_terminals()} for tree in trees.values()))
        namespace = {"re": re, "shared_loci": shared}
        exec(compile(ast.Module(body=[representative], type_ignores=[]), "representative_selection", "exec"), namespace)
        namespace["owner"] = types.SimpleNamespace(representative_tip_names=namespace["representative_tip_names"])
        exec(compile(ast.Module(body=[selection], type_ignores=[]), "restored_figure_selection", "exec"), namespace)
        manifest = json.loads((manuscript / "figures/comment_corrections_20260914/Figure_3_restored_manifest.json").read_text())
        for region, tree in trees.items():
            self.assertEqual(namespace["select_tree_tips"](tree), set(manifest["trees"][region]["selected_tips"]))

    def test_actual_inputs_preserve_trees_and_correct_comparison(self):
        spec = importlib.util.spec_from_file_location(
            "focused_shared_test_owner", SOURCE / "build_focused_main_panels.py"
        )
        owner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(owner)
        source_bytes = {
            name: (SOURCE / name).read_bytes() for name in owner.TREE_FILES.values()
        }
        with tempfile.TemporaryDirectory(prefix="hml2-shared-tree-test-") as tmp:
            destination = Path(tmp)
            for name in owner.TREE_FILES.values():
                shutil.copy2(SOURCE / name, destination / name)
            owner.OUT = destination
            with contextlib.redirect_stdout(io.StringIO()):
                owner.main()
            with (destination / "focused_main_nearest_neighbors.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            by_key = {(row["region"], row["focal_locus"]): row for row in rows}
            self.assertEqual(by_key["LTR", "10q24.2"]["nearest_locus"], "12q14.1")
            self.assertEqual(by_key["Pol", "10q24.2"]["nearest_locus"], "6q14.1")
            self.assertEqual(by_key["LTR", "19p12c"]["nearest_locus"], "12q14.1")
            self.assertEqual(by_key["Pol", "19p12c"]["nearest_locus"], "22q11.21")
            self.assertEqual(by_key["LTR", "19p12d"]["status"], "focal_not_in_shared_locus_set")
            self.assertEqual(by_key["LTR", "19p12d"]["nearest_locus"], "")
            self.assertEqual(by_key["Pol", "19p12d"]["status"], "not_in_callable_full_tree")
            summary = json.loads((destination / "focused_main_verification.json").read_text())
            self.assertEqual(summary["comparison_locus_count"], 64)
            shared = set(summary["comparison_loci"])
            self.assertNotIn("19p12d", shared)
            for row in rows:
                if row["nearest_locus"]:
                    self.assertIn(row["focal_locus"], shared)
                    self.assertIn(row["nearest_locus"], shared)
            for name, data in source_bytes.items():
                self.assertEqual((SOURCE / name).read_bytes(), data)
                self.assertEqual((destination / name).read_bytes(), data)
            # Comparison-set correction must not change either existing focused tree.
            for region in ("ltr", "pol"):
                name = f"focused_{region}_tree.nwk"
                self.assertEqual((destination / name).read_bytes(), (SOURCE / name).read_bytes())


if __name__ == "__main__":
    unittest.main()
