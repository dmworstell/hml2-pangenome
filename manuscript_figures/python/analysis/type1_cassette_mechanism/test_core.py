#!/usr/bin/env python3

import json
import math
import tempfile
import unittest
from pathlib import Path

from run_analysis import (
    Candidate,
    Window,
    audit_delta_records,
    build_coordinate_map,
    candidate_sort_key,
    leave_one_window_out_switching,
    matched_site_comparison,
    nearest_neighbor_switching,
    pair_distance,
    site_diversity,
    type1_locus_bootstrap_difference,
    wilson_interval,
    window_columns,
)


class CoordinateTests(unittest.TestCase):
    def test_exact_kcon_coordinate_map_ignores_alignment_insertions(self):
        aligned = "AC-GT--A"
        coordinates = build_coordinate_map(aligned, "ACGTA")
        self.assertEqual(coordinates, [0, 1, 3, 4, 7])
        window = Window("test", "T", ((1, 4),), "control")
        self.assertEqual(window_columns(window, coordinates), [1, 3, 4])

    def test_coordinate_map_rejects_changed_anchor(self):
        with self.assertRaises(ValueError):
            build_coordinate_map("AC-T", "ACG")


class DistanceTests(unittest.TestCase):
    def test_pair_distance_excludes_gaps_and_ambiguous_bases(self):
        result = pair_distance("ACGT-N", "ATGTAA", list(range(6)), minimum_overlap=4)
        self.assertEqual(result, (4, 1, 0.25))
        self.assertIsNone(pair_distance("ACGT-N", "ATGTAA", list(range(6)), minimum_overlap=5))

    def test_site_diversity_is_pairwise_mismatch_probability(self):
        sequences = ["A", "A", "C", "N"]
        callable_count, diversity = site_diversity(sequences, 0, minimum_callable=3)
        self.assertEqual(callable_count, 3)
        self.assertAlmostEqual(diversity, 2 / 3)


class ResamplingTests(unittest.TestCase):
    def test_matched_site_comparison_is_stratified_and_deterministic(self):
        cassette = {4: [0.0] * 20, 5: [0.01] * 10}
        control = {4: [0.10] * 15, 5: [0.11] * 12}
        first = matched_site_comparison(cassette, control, resamples=100, seed=7)
        second = matched_site_comparison(cassette, control, resamples=100, seed=7)
        self.assertEqual(first, second)
        self.assertEqual(first["matched_sites_per_window"], 25)
        self.assertLess(first["difference_cassette_minus_control"], 0)
        self.assertLess(first["bootstrap_ci_high"], 0)

    def test_locus_bootstrap_difference_uses_paired_type1_loci(self):
        from run_analysis import Representative

        representatives = [
            Representative(
                analysis_id=f"r{index}", locus=f"L{index}", type_call="TypeI",
                id_full=f"id{index}", person=f"HG0000{index}", haplotype="h1",
                strand="+", source_identifier="source", fasta=Path(f"{index}.fa"),
                root_label="test", sequence="", fasta_sha256="hash",
            )
            for index in range(3)
        ]
        # Cassette is invariant; every control pair differs at one of four sites.
        aligned = {
            "r0": "AAAA" + "AAAA",
            "r1": "AAAA" + "CAAA",
            "r2": "AAAA" + "GAAA",
        }
        result = type1_locus_bootstrap_difference(
            representatives, aligned, [0, 1, 2, 3], [4, 5, 6, 7],
            {"r0", "r1", "r2"}, {"r0", "r1", "r2"},
            resamples=200, seed=11,
        )
        self.assertEqual(result["paired_callable_type1_loci"], 3)
        self.assertLess(result["difference_cassette_minus_control"], 0)
        self.assertLessEqual(result["locus_bootstrap_ci_high"], 0)


class SelectionTests(unittest.TestCase):
    def test_hprc_root_precedes_graph_then_identifier(self):
        common = dict(
            locus="1q22", type_call="TypeI", person="HG00097", haplotype="h1",
            structure="Provirus", strand="+", source_identifier="source",
        )
        graph = Candidate(id_full="AAA", fasta=Path("graph.fa"), root_rank=1, root_label="graph", **common)
        hprc = Candidate(id_full="ZZZ", fasta=Path("hprc.fa"), root_rank=0, root_label="hprc", **common)
        self.assertLess(candidate_sort_key(hprc), candidate_sort_key(graph))


class DeltaAuditTests(unittest.TestCase):
    def test_only_aggregate_eligible_rows_define_both_call_status(self):
        record = {
            "direct_delta292_evidence": {
                "aggregate_eligibility_mask": {
                    "rows": [
                        {"row_index": 0, "aggregate_eligible": True},
                        {"row_index": 1, "aggregate_eligible": True},
                        {"row_index": 2, "aggregate_eligible": False},
                    ]
                },
                "rows": [
                    {"row_index": 0, "type_call": "TypeI"},
                    {"row_index": 1, "type_call": "TypeII"},
                    {"row_index": 2, "type_call": "conflict"},
                ],
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "1q22.json"
            path.write_text(json.dumps(record))
            rows, provenance = audit_delta_records({"1q22": "TypeI"}, Path(directory))
        self.assertEqual(rows[0]["eligible_rows"], 2)
        self.assertEqual(rows[0]["eligible_conflict_rows"], 0)
        self.assertEqual(rows[0]["both_typeI_and_typeII_present"], "true")
        self.assertEqual(len(provenance), 1)


class SwitchingTests(unittest.TestCase):
    def test_switching_uses_deterministic_backbone_mode(self):
        rows = []
        for window, neighbor in [
            ("backbone_1000_2000", "A"),
            ("backbone_2000_3000", "A"),
            ("backbone_3000_4000", "B"),
            ("backbone_4000_5000", "A"),
            ("backbone_5000_6000", "B"),
            ("backbone_7293_8293", "B"),
            ("type1_cassette_flanks", "C"),
        ]:
            rows.append({
                "status": "callable", "type1_locus": "T1", "window": window,
                "nearest_type2_locus": neighbor,
            })
        per_locus, summary = nearest_neighbor_switching(rows)
        self.assertEqual(per_locus[0]["backbone_modal_nearest_type2"], "A")
        self.assertEqual(per_locus[0]["cassette_differs_from_backbone_mode"], "true")
        mode = next(row for row in summary if row["comparison"] == "backbone_mode")
        self.assertEqual(mode["switch_fraction"], 1.0)

    def test_wilson_interval_is_bounded(self):
        low, high = wilson_interval(0, 10)
        self.assertEqual(low, 0.0)
        self.assertGreater(high, 0)
        self.assertTrue(math.isnan(wilson_interval(0, 0)[0]))

    def test_leave_one_out_backbone_is_negative_control(self):
        from run_analysis import CASSETTE_KEY, CONTROL_KEYS

        rows = []
        neighbors = ["A", "A", "B", "A", "B", "B"]
        for window, neighbor in zip(CONTROL_KEYS, neighbors):
            rows.append({
                "status": "callable", "type1_locus": "T1", "window": window,
                "nearest_type2_locus": neighbor,
            })
        rows.append({
            "status": "callable", "type1_locus": "T1", "window": CASSETTE_KEY,
            "nearest_type2_locus": "C",
        })
        per_locus, summary = leave_one_window_out_switching(rows)
        self.assertEqual(len(per_locus), 7)
        cassette = next(row for row in summary if row["focal_window"] == CASSETTE_KEY)
        self.assertEqual(cassette["negative_control_role"], "cassette_test")
        self.assertEqual(cassette["switch_fraction"], 1.0)
        backbone = [row for row in summary if row["focal_window"] in CONTROL_KEYS]
        self.assertTrue(all(row["negative_control_role"].startswith("leave_one_backbone") for row in backbone))


if __name__ == "__main__":
    unittest.main()
