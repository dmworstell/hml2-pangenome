#!/usr/bin/env python3

import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("build_biologically_annotated_orf.py")
SPEC = importlib.util.spec_from_file_location("biological_orf", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class BiologicalAnnotationTest(unittest.TestCase):
    def row(self, locus, source, structure="Fragment", sample="HG00001"):
        return {
            "Locus": locus,
            "ID_Full": f"{sample}_{locus}",
            "ID": sample,
            "Source_Identifier": source,
            "Structure": structure,
        }

    def test_keeps_only_host_anchored_fourq_as_fourq(self):
        anchored = "HG00001#1#contigA"
        row = MODULE.annotate_row(
            self.row(MODULE.FOURQ, anchored + ":100-200"), {anchored}
        )
        self.assertEqual(row["Locus"], MODULE.FOURQ)
        self.assertEqual(row["insertion_event_group"], MODULE.TYPE2_EVENT)

        unresolved = MODULE.annotate_row(
            self.row(MODULE.FOURQ, "HG00001#1#contigB:100-200"), {anchored}
        )
        self.assertEqual(unresolved["Locus"], "HML-2_acro_type2")

    def test_folds_present_acro_rows_and_retains_arm_absence(self):
        present = MODULE.annotate_row(
            self.row("HML-2_21p13", "HG00001#1#contigA:100-200"), set()
        )
        absent = MODULE.annotate_row(
            self.row("HML-2_21p13", "Placeholder", "Absent"), set()
        )
        self.assertEqual(present["Locus"], "HML-2_acro_type2")
        self.assertEqual(absent["Locus"], "HML-2_21p13")
        self.assertEqual(absent["analysis_include"], "1")

    def test_retains_but_excludes_8q_alias(self):
        alias = MODULE.annotate_row(
            self.row("HML-2_8q24.3b", "HG00001#1#contigA:100-200"), set()
        )
        self.assertEqual(alias["Locus"], "HML-2_8q24.3c")
        self.assertEqual(alias["analysis_include"], "0")
        self.assertEqual(alias["record_relationship"], "alias_of_8q24.3c")

    def test_retains_but_excludes_depth_unsupported_artifact(self):
        id_full = "HG00001_HML-2_1q22"
        artifact = MODULE.annotate_row(
            self.row("HML-2_1q22", "HG00001#1#contigA:100-200"),
            set(),
            cnv_decisions={id_full: ("assembly_artifact", "0.999")},
        )
        self.assertEqual(artifact["ID_Full"], id_full)
        self.assertEqual(artifact["analysis_include"], "0")
        self.assertEqual(artifact["cnv_qc_state"], "assembly_artifact")

    def test_authenticated_copy_keeps_id_and_inherits_physical_assignment(self):
        authenticated = MODULE.annotate_row(
            self.row("HML-2_22p13", "HG00001#1#contigA:100-200"),
            set(),
            cnv_decisions={
                "HG00001_HML-2_22p13": ("authenticated_segdup", "1.0")
            },
        )
        physical = MODULE.annotate_row(
            self.row(MODULE.FOURQ, "HG00001#1#contigA:100-200"),
            set(),
            fourq_decisions={
                "HG00001_HML-2_4q35.2_hg38": (
                    "HML-2_22p13",
                    "primary_host_alignment_places_this_record_at_22p13_not_4q35.2",
                )
            },
        )
        annotated, excluded, unresolved = MODULE.annotate_duplicate_present_rows(
            [authenticated, physical]
        )
        self.assertEqual(excluded, 1)
        self.assertEqual(unresolved, 0)
        self.assertEqual(len(annotated), 2)
        retained = [row for row in annotated if row["analysis_include"] == "1"]
        self.assertEqual(len(retained), 1)
        self.assertEqual(retained[0]["ID_Full"], "HG00001_HML-2_22p13")
        self.assertEqual(retained[0]["Locus"], "HML-2_22p13")
        self.assertEqual(retained[0]["cnv_qc_state"], "authenticated_segdup")

    def test_chromosome4_contig_resolves_fourp_acro_duplicate_label(self):
        qname = "HG00001#1#chr4contig"
        fourp = MODULE.annotate_row(
            self.row(MODULE.FOURP, qname + ":100-200"), {qname}
        )
        acro = MODULE.annotate_row(
            self.row("HML-2_15p13a", qname + ":100-200"), {qname}
        )
        annotated, excluded, unresolved = MODULE.annotate_duplicate_present_rows(
            [fourp, acro]
        )
        self.assertEqual(excluded, 1)
        self.assertEqual(unresolved, 0)
        self.assertEqual({row["Locus"] for row in annotated}, {MODULE.FOURP})
        self.assertEqual(
            sum(row["analysis_include"] == "1" for row in annotated), 1
        )


if __name__ == "__main__":
    unittest.main()
