"""Test exact native-source and competing-window boundaries of the graph consumer."""
from pathlib import Path
import copy
import importlib.util
import unittest

BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("acro_graph_consumer", BASE / "resolve_acrocentric_catalog.py")
consumer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(consumer)
audit = consumer.load_module("acro_graph_audit", consumer.BIOLOGICAL / "audit_exact_duplicate_locus_assignments.py")


def fixture():
    chrom, start, end, _ = audit.TARGETS[("t2t", "13p13")]
    length = end-start
    query_length = 100000+length+150000
    hit = dict(build="CHM13", query="g0", query_length=query_length, query_start=0,
        query_end=query_length, strand="+", target=chrom, target_start=start-100000,
        target_end=end+150000, mapq=60, tp="P", cg=f"{query_length}M")
    window = dict(query="g0", native_path="HG00001#1#ctg", window_start=1000000,
        window_end=1000000+query_length, window_length=query_length, full_reference_hits=[hit])
    match = dict(element_start0=1100000, element_end0=1100000+length,
        element_start_in_window=100000, element_end_in_window=100000+length,
        retained_sequence_orientation_vs_native="+")
    row = dict(Source_Identifier=f"HG00001#1#ctg:{1100000-500}-{1100000+length+500}",
        retained_sequence_length=length, native_windows=[dict(query="g0", exact_sequence_matches=[match])],
        unique_exact_sequence_positions=1, capture_status="EXACT_ELEMENT_MATCH")
    return row, {"g0": window}


class GraphAssignmentConsumerTests(unittest.TestCase):
    def test_exact_native_subspan_reaches_same_assignment_projector(self):
        row, windows = fixture()
        candidates, _, positions = consumer.graph_source_candidates(audit, row, windows)
        self.assertEqual(len(positions), 1)
        self.assertEqual([c["locus"] for c in candidates], ["13p13"])
        self.assertTrue(candidates[0]["interior100"])
        self.assertTrue(candidates[0]["both5"])

    def test_overlapping_window_representations_do_not_double_count(self):
        row, windows = fixture()
        second = copy.deepcopy(windows["g0"])
        second["query"] = "g1"
        second["full_reference_hits"][0]["query"] = "g1"
        windows["g1"] = second
        row["native_windows"].append({**copy.deepcopy(row["native_windows"][0]), "query": "g1"})
        candidates, hits, positions = consumer.graph_source_candidates(audit, row, windows)
        self.assertEqual((len(candidates), len(hits), len(positions)), (1, 2, 1))

    def test_competing_full_reference_hit_is_not_discarded(self):
        row, windows = fixture()
        other = copy.deepcopy(windows["g0"]["full_reference_hits"][0])
        chrom, start, _, _ = audit.TARGETS[("t2t", "15p13b")]
        other.update(target=chrom, target_start=start-100000, target_end=start-100000+other["query_length"])
        windows["g0"]["full_reference_hits"].append(other)
        candidates, _, _ = consumer.graph_source_candidates(audit, row, windows)
        self.assertEqual({c["locus"] for c in candidates}, {"13p13", "15p13b"})
        self.assertEqual(sum(c["interior100"] for c in candidates), 2)

    def test_multiple_exact_native_positions_are_not_assigned(self):
        row, windows = fixture()
        other = copy.deepcopy(row["native_windows"][0]["exact_sequence_matches"][0])
        for key in ("element_start0", "element_end0", "element_start_in_window", "element_end_in_window"):
            other[key] += 1
        row["native_windows"][0]["exact_sequence_matches"].append(other)
        row["unique_exact_sequence_positions"] = 2
        row["capture_status"] = "MULTIPLE_EXACT_ELEMENT_POSITIONS"
        candidates, hits, positions = consumer.graph_source_candidates(audit, row, windows)
        self.assertEqual((len(candidates), len(hits), len(positions)), (0, 0, 2))

    def test_source_match_must_be_inside_original_native_interval(self):
        row, windows = fixture()
        row["Source_Identifier"] = "HG00001#1#ctg:1-100"
        with self.assertRaises(AssertionError):
            consumer.graph_source_candidates(audit, row, windows)

    def test_reverse_paf_uses_original_query_coordinate(self):
        hit = dict(build="CHM13", query="g0", query_length=1000, query_start=50,
            query_end=850, strand="-", target="chr13", target_start=5000,
            target_end=5800, mapq=60, tp="P", cg="800M")
        alignment = consumer.graph_paf_alignment(audit, hit, "13p13")
        self.assertEqual(alignment.cigar, "150S800M50S")
        self.assertEqual(audit.source_to_target_bases(alignment, 120, 220, 5630, 5730), (100, 100))

    def test_truncated_cigar_is_rejected(self):
        row, windows = fixture()
        windows["g0"]["full_reference_hits"][0]["cg"] = "50M"
        with self.assertRaises(AssertionError):
            consumer.graph_source_candidates(audit, row, windows)


if __name__ == "__main__":
    unittest.main()
