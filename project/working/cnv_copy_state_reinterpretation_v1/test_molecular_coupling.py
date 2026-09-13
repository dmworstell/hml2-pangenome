#!/usr/bin/env python3
import hashlib
import json
import unittest
from pathlib import Path

from copy_state_authority import STATE_ARTIFACT, STATE_TRUE


LANE = Path(__file__).resolve().parent
RESULTS = LANE / "results"
PRESERVED_WEIGHT_SHA256 = "597de702cbc995d39dd21ff51c5e549d4b8f83326d02ce708afb6dfd11d6bcda"
PRESERVED_CALIBRATION_SHA256 = "8b56d58036a3728069971bb167e5fc474c15f2e13bf8f86a32e7df7d917c50d8"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MolecularCouplingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hg = json.loads(
            (RESULTS / "hg00423_1q22_copy_state_authority.v1.json").read_text()
        )
        cls.amendment = json.loads(
            (RESULTS / "MOLECULAR_COUPLING_AMENDMENT.v1.json").read_text()
        )

    def test_artifact_candidate_is_measurement_only(self):
        law = self.hg["state_contract"][STATE_ARTIFACT]["molecular_coupling"]
        self.assertEqual(law["candidate_bytes_role"], "measurement_or_assembly_error_evidence_only")
        self.assertFalse(law["biological_lineage_created"])
        self.assertFalse(law["included_in_biological_sequence_likelihood"])
        self.assertFalse(law["second_insertion_root_created"])

    def test_true_duplication_birth_bases_are_coupled(self):
        law = self.hg["state_contract"][STATE_TRUE]["molecular_coupling"]
        self.assertEqual(
            law["duplication_birth_base_constraint"],
            "daughter_birth_sequence_equals_parent_sequence_at_duplication_time",
        )
        self.assertEqual(
            law["parent_tip_candidates"],
            ["HG00423_pat_hprc_r2_v1.0.1_HML-2_1q22_alt1"],
        )
        self.assertFalse(law["parent_assignment_latent"])
        self.assertFalse(law["second_insertion_root_created"])

    def test_required_joint_marginalization_is_explicit(self):
        required = set(
            self.hg["state_contract"][STATE_TRUE]["molecular_coupling"]["required_marginalization"]
        )
        self.assertEqual(required, {
            "shared_ancestral_bases_at_duplication_birth",
            "duplication_time",
            "parent_post_duplication_mutation_history",
            "daughter_post_duplication_mutation_history",
            "post_duplication_gene_conversion_history",
            "parent_sequence_measurement_error",
            "daughter_sequence_measurement_error",
        })

    def test_rare_duplication_event_error_is_not_default(self):
        sensitivity = self.hg["state_contract"][STATE_TRUE]["molecular_coupling"][
            "duplication_event_sequence_error_sensitivity"
        ]
        self.assertFalse(sensitivity["enabled"])
        self.assertIsNone(sensitivity["support_authority"])
        self.assertIn("independent evidence", sensitivity["activation_rule"])

    def test_numeric_weights_and_calibration_are_byte_preserved(self):
        self.assertEqual(
            sha256(RESULTS / "current_unmeasured_copy_state_weights.v1.tsv"),
            PRESERVED_WEIGHT_SHA256,
        )
        self.assertEqual(
            sha256(RESULTS / "real_control_calibration.v1.json"),
            PRESERVED_CALIBRATION_SHA256,
        )
        preserved = self.amendment["preserved_numeric_authorities"]
        self.assertFalse(preserved["cnv_weights_changed"])
        self.assertFalse(self.amendment["exact4_outputs_modified"])


if __name__ == "__main__":
    unittest.main()
