#!/usr/bin/env python3
"""Offline regression tests for the manifest-driven CNV-v2 workflow.

The suite deliberately uses only tiny synthetic sequence/read placeholders.  It
does not contact object storage, a scheduler, or any biological data source.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SRC = ROOT / "src"
PLANNER = SRC / "cnv_plan.py"
CONTROLLER = SRC / "01_cnv_controller_v2.sh"
WORKER = SRC / "02_cnv_alignment_worker_v2.sh"
ANALYSIS = SRC / "03_cnv_analysis_v2.sh"
PLOTTER = SRC / "plot_cnv_depth_v2.py"
README = ROOT / "README.md"

_PLANNER_SPEC = importlib.util.spec_from_file_location("cnv_v2_plan_under_test", PLANNER)
assert _PLANNER_SPEC is not None and _PLANNER_SPEC.loader is not None
CNV_PLAN = importlib.util.module_from_spec(_PLANNER_SPEC)
_PLANNER_SPEC.loader.exec_module(CNV_PLAN)

EXPECTED_FILES = (PLANNER, CONTROLLER, WORKER, ANALYSIS, PLOTTER, README)

MANIFEST_FIELDS = (
    "sample_id",
    "locus_id",
    "part_id",
    "expected_part_count",
    "read_uri",
    "read_version",
    "read_etag",
    "read_size",
    "read_sha256",
    "read_format",
    "minimap2_preset",
    "assembly1_id",
    "assembly1_fasta",
    "assembly1_fai",
    "assembly1_fasta_sha256",
    "assembly1_fai_sha256",
    "assembly1_contig",
    "assembly1_body_start",
    "assembly1_body_end",
    "assembly2_id",
    "assembly2_fasta",
    "assembly2_fai",
    "assembly2_fasta_sha256",
    "assembly2_fai_sha256",
    "assembly2_contig",
    "assembly2_body_start",
    "assembly2_body_end",
    "full_flank_bp",
    "bait_min_aligned_bp",
    "kcon_fasta",
    "kcon_sha256",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def command_lines(text: str) -> list[str]:
    """Return logical shell lines with comments/blanks removed."""
    result: list[str] = []
    pending = ""
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        continued = line.endswith("\\")
        if continued:
            line = line[:-1].rstrip()
        pending = (pending + " " + line).strip()
        if not continued:
            result.append(pending)
            pending = ""
    if pending:
        result.append(pending)
    return result


class SyntheticRun:
    """A tiny, non-biological manifest fixture with two declared parts."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.inputs = root / "inputs"
        self.inputs.mkdir(parents=True)
        self.run_root = root / "runs"
        self.run_root.mkdir()

        self.a1 = self._fasta("assembly_one.fa", "ctgA", "ACGT" * 20)
        self.a2 = self._fasta("assembly_two.fa", "ctgB", "TGCA" * 20)
        self.kcon = self._fasta("kcon.fa", "KCON", "ACGT" * 5)[0]
        self.reads: list[Path] = []
        for number in (1, 2):
            path = self.inputs / f"reads.part{number}.bam"
            path.write_bytes(f"synthetic-read-part-{number}\n".encode("ascii"))
            self.reads.append(path)
        self.rows = [self._row(1), self._row(2)]
        self.manifest = self.root / "manifest.tsv"
        self.write_manifest()

    def _fasta(self, name: str, contig: str, sequence: str) -> tuple[Path, Path]:
        fasta = self.inputs / name
        fasta.write_text(f">{contig}\n{sequence}\n", encoding="ascii")
        fai = Path(str(fasta) + ".fai")
        offset = len(contig) + 2
        fai.write_text(
            f"{contig}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence) + 1}\n",
            encoding="ascii",
        )
        return fasta, fai

    def _row(self, part_id: int) -> dict[str, str]:
        read = self.reads[part_id - 1]
        a1_fasta, a1_fai = self.a1
        a2_fasta, a2_fai = self.a2
        return {
            "sample_id": "sampleA",
            "locus_id": "locusA",
            "part_id": str(part_id),
            "expected_part_count": "2",
            "read_uri": read.resolve().as_uri(),
            "read_version": f"immutable-v{part_id}",
            "read_etag": hashlib.sha256(f"etag-{part_id}".encode()).hexdigest(),
            "read_size": str(read.stat().st_size),
            "read_sha256": sha256(read),
            "read_format": "bam",
            "minimap2_preset": "map-ont",
            "assembly1_id": "hap1",
            "assembly1_fasta": str(a1_fasta.resolve()),
            "assembly1_fai": str(a1_fai.resolve()),
            "assembly1_fasta_sha256": sha256(a1_fasta),
            "assembly1_fai_sha256": sha256(a1_fai),
            "assembly1_contig": "ctgA",
            "assembly1_body_start": "3",
            "assembly1_body_end": "15",
            "assembly2_id": "hap2",
            "assembly2_fasta": str(a2_fasta.resolve()),
            "assembly2_fai": str(a2_fai.resolve()),
            "assembly2_fasta_sha256": sha256(a2_fasta),
            "assembly2_fai_sha256": sha256(a2_fai),
            "assembly2_contig": "ctgB",
            "assembly2_body_start": "70",
            "assembly2_body_end": "78",
            "full_flank_bp": "10",
            "bait_min_aligned_bp": "3",
            "kcon_fasta": str(self.kcon.resolve()),
            "kcon_sha256": sha256(self.kcon),
        }

    def write_manifest(
        self,
        rows: list[dict[str, str]] | None = None,
        fields: tuple[str, ...] = MANIFEST_FIELDS,
    ) -> None:
        with self.manifest.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(self.rows if rows is None else rows)

    def multi_target_rows(self) -> list[dict[str, str]]:
        rows = [dict(row) for row in self.rows]
        for part_id, template in enumerate(self.rows, 1):
            read = self.inputs / f"sampleB.reads.part{part_id}.bam"
            read.write_bytes(f"synthetic-sampleB-part-{part_id}\n".encode("ascii"))
            row = dict(template)
            row.update(
                {
                    "sample_id": "sampleB",
                    "locus_id": "locusB",
                    "read_uri": read.resolve().as_uri(),
                    "read_version": f"sampleB-immutable-v{part_id}",
                    "read_etag": hashlib.sha256(f"sampleB-etag-{part_id}".encode()).hexdigest(),
                    "read_size": str(read.stat().st_size),
                    "read_sha256": sha256(read),
                }
            )
            rows.append(row)
        return rows

    def render(self, expected_digest: str | None = None) -> subprocess.CompletedProcess[str]:
        digest = sha256(self.manifest) if expected_digest is None else expected_digest
        return subprocess.run(
            [
                sys.executable,
                str(PLANNER),
                "render",
                "--manifest",
                str(self.manifest),
                "--expected-manifest-sha256",
                digest,
                "--run-root",
                str(self.run_root),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )

    def rendered_json(self) -> list[tuple[Path, object]]:
        values: list[tuple[Path, object]] = []
        for path in self.run_root.rglob("*.json"):
            try:
                values.append((path, json.loads(path.read_text(encoding="utf-8"))))
            except json.JSONDecodeError:
                continue
        return values


class TestSourceInventory(unittest.TestCase):
    def test_expected_forward_workflow_files_exist(self) -> None:
        for path in EXPECTED_FILES:
            with self.subTest(path=path.name):
                self.assertTrue(path.is_file(), f"missing {path}")

    def test_original_scripts_are_not_used_as_runtime_sources(self) -> None:
        all_source = "\n".join(read_text(path) for path in EXPECTED_FILES)
        self.assertNotIn("/original/", all_source)
        self.assertNotRegex(all_source, r"(?m)^\s*(?:source|\.)\s+.*original/")


class TestPlannerManifestAndCoordinates(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = SyntheticRun(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def assert_rejected(self, result: subprocess.CompletedProcess[str]) -> None:
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_valid_manifest_renders_dry_run_without_external_tools(self) -> None:
        shim_dir = Path(self.temporary.name) / "shims"
        shim_dir.mkdir()
        marker = Path(self.temporary.name) / "external-called"
        for name in ("aws", "minimap2", "samtools", "sbatch", "srun"):
            shim = shim_dir / name
            shim.write_text(
                '#!/bin/sh\nprintf \'%s\\n\' "$0" >> "$CNV_TEST_MARKER"\nexit 97\n',
                encoding="utf-8",
            )
            shim.chmod(0o755)
        environment = os.environ.copy()
        environment["PATH"] = str(shim_dir) + os.pathsep + environment.get("PATH", "")
        environment["CNV_TEST_MARKER"] = str(marker)
        result = subprocess.run(
            [
                "bash",
                str(CONTROLLER),
                "--manifest",
                str(self.fixture.manifest),
                "--expected-manifest-sha256",
                sha256(self.fixture.manifest),
                "--run-root",
                str(self.fixture.run_root),
            ],
            cwd=str(ROOT),
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse(marker.exists(), "dry-run invoked an external workflow tool")
        self.assertTrue(self.fixture.rendered_json(), "dry-run did not leave a reviewable plan")

    def test_manifest_digest_is_mandatory_and_exact(self) -> None:
        self.assert_rejected(self.fixture.render("0" * 64))
        result = subprocess.run(
            [
                sys.executable,
                str(PLANNER),
                "render",
                "--manifest",
                str(self.fixture.manifest),
                "--run-root",
                str(self.fixture.run_root),
            ],
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )
        self.assert_rejected(result)

    def test_header_is_exact_and_unknown_column_is_rejected(self) -> None:
        self.fixture.write_manifest(fields=MANIFEST_FIELDS + ("ambient_extra",))
        self.assert_rejected(self.fixture.render())
        self.fixture.write_manifest(fields=MANIFEST_FIELDS[:-1])
        self.assert_rejected(self.fixture.render())

    def test_exact_part_set_is_required(self) -> None:
        self.fixture.write_manifest(rows=[self.fixture.rows[0]])
        self.assert_rejected(self.fixture.render())
        duplicate = [dict(self.fixture.rows[0]), dict(self.fixture.rows[0])]
        self.fixture.write_manifest(rows=duplicate)
        self.assert_rejected(self.fixture.render())

    def test_exactly_two_distinct_assemblies_are_required(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        for row in rows:
            row["assembly2_id"] = row["assembly1_id"]
        self.fixture.write_manifest(rows=rows)
        self.assert_rejected(self.fixture.render())

    def test_same_biological_assembly_cannot_be_relabelled_as_two(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        for row in rows:
            row.update(
                {
                    "assembly2_fasta": row["assembly1_fasta"],
                    "assembly2_fai": row["assembly1_fai"],
                    "assembly2_fasta_sha256": row["assembly1_fasta_sha256"],
                    "assembly2_fai_sha256": row["assembly1_fai_sha256"],
                    "assembly2_contig": row["assembly1_contig"],
                    "assembly2_body_start": row["assembly1_body_start"],
                    "assembly2_body_end": row["assembly1_body_end"],
                }
            )
        self.fixture.write_manifest(rows=rows)
        self.assert_rejected(self.fixture.render())

    def test_same_physical_read_object_cannot_fill_two_part_ids(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        for field in (
            "read_uri",
            "read_version",
            "read_etag",
            "read_size",
            "read_sha256",
        ):
            rows[1][field] = rows[0][field]
        self.fixture.write_manifest(rows=rows)
        self.assert_rejected(self.fixture.render())

    def test_same_local_uri_and_content_refuses_after_version_label_change(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        rows[1].update(
            {
                "read_uri": rows[0]["read_uri"],
                "read_version": "caller-changed-label-only",
                "read_etag": "caller-changed-etag-label",
                "read_size": rows[0]["read_size"],
                "read_sha256": rows[0]["read_sha256"],
            }
        )
        self.fixture.write_manifest(rows=rows)
        result = self.fixture.render()
        self.assert_rejected(result)
        self.assertRegex(result.stderr.lower(), r"content sha-256|canonical local")

    def test_two_local_paths_with_identical_content_refuse(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        alias = self.fixture.inputs / "same-content-under-another-name.bam"
        alias.write_bytes(self.fixture.reads[0].read_bytes())
        rows[1].update(
            {
                "read_uri": alias.resolve().as_uri(),
                "read_version": "different-arbitrary-label",
                "read_etag": "different-arbitrary-etag",
                "read_size": str(alias.stat().st_size),
                "read_sha256": sha256(alias),
            }
        )
        self.fixture.write_manifest(rows=rows)
        result = self.fixture.render()
        self.assert_rejected(result)
        self.assertIn("content sha-256", result.stderr.lower())

    def test_repeated_s3_bucket_key_version_refuses_despite_other_labels(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        rows[0].update({"read_uri": "s3://fixed-bucket/a%20key", "read_version": "VersionId-1"})
        rows[1].update(
            {
                "read_uri": "s3://FIXED-BUCKET/a%20key",
                "read_version": "VersionId-1",
                "read_etag": "caller-changed-etag",
            }
        )
        self.fixture.write_manifest(rows=rows)
        result = self.fixture.render()
        self.assert_rejected(result)
        self.assertIn("bucket/key/versionid", result.stderr.lower())

    def test_repeated_s3_content_refuses_across_distinct_object_labels(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        rows[0].update({"read_uri": "s3://fixed-bucket/key-one", "read_version": "VersionId-1"})
        rows[1].update(
            {
                "read_uri": "s3://fixed-bucket/key-two",
                "read_version": "VersionId-2",
                "read_sha256": rows[0]["read_sha256"],
            }
        )
        self.fixture.write_manifest(rows=rows)
        result = self.fixture.render()
        self.assert_rejected(result)
        self.assertIn("content sha-256", result.stderr.lower())

    def test_s3_object_uri_rejects_query_and_fragment(self) -> None:
        for suffix in ("?versionId=ambient", "#ambient"):
            with self.subTest(suffix=suffix):
                rows = [dict(row) for row in self.fixture.rows]
                rows[0]["read_uri"] = "s3://fixed-bucket/exact-key" + suffix
                self.fixture.write_manifest(rows=rows)
                self.assert_rejected(self.fixture.render())

    def test_target_configuration_must_match_across_parts(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        rows[1]["assembly2_body_start"] = "69"
        self.fixture.write_manifest(rows=rows)
        self.assert_rejected(self.fixture.render())

    def test_reference_digest_drift_is_rejected(self) -> None:
        self.fixture.a1[0].write_text(">ctgA\n" + "A" * 80 + "\n", encoding="ascii")
        self.assert_rejected(self.fixture.render())

    def test_body_must_be_within_declared_fai_contig(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        for row in rows:
            row["assembly2_body_end"] = "81"
        self.fixture.write_manifest(rows=rows)
        self.assert_rejected(self.fixture.render())

    def test_separate_beds_are_fai_clamped_and_contigs_are_unique(self) -> None:
        result = self.fixture.render()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        beds = {path.name: path for path in self.fixture.run_root.rglob("*.bed")}
        for required in ("core_body_union.bed", "full_window.bed", "outer_flanks.bed"):
            self.assertIn(required, beds)

        def rows(name: str) -> list[tuple[str, int, int]]:
            parsed = []
            for line in beds[name].read_text(encoding="utf-8").splitlines():
                fields = line.split("\t")
                parsed.append((fields[0], int(fields[1]), int(fields[2])))
            return parsed

        core = rows("core_body_union.bed")
        full = rows("full_window.bed")
        flanks = rows("outer_flanks.bed")
        self.assertEqual(len(core), 2)
        self.assertEqual(len({row[0] for row in core}), 2)
        self.assertTrue(all("__" in row[0] for row in core))
        self.assertEqual({(start, end) for _, start, end in core}, {(3, 15), (70, 78)})
        self.assertIn((core[0][0], 0, 25), full)
        self.assertIn((core[1][0], 60, 80), full)
        self.assertTrue(all(0 <= start < end <= 80 for _, start, end in full))
        self.assertGreaterEqual(len(flanks), 2)
        for contig, start, end in flanks:
            self.assertIn(contig, {row[0] for row in core})
            self.assertGreaterEqual(start, 0)
            self.assertLessEqual(end, 80)

        plan = next(
            value
            for _, value in self.fixture.rendered_json()
            if isinstance(value, dict)
            and value.get("schema_version") == "cnv_v2_plan_1"
        )
        by_id = {
            assembly["assembly_id"]: assembly
            for assembly in plan["target"]["assemblies"]
        }
        self.assertEqual(
            by_id["hap1"]["flank_completeness"],
            "TWO_SIDED_BOUNDARY_TRUNCATED_FLANK",
        )
        hap1_sides = {row["side"]: row for row in by_id["hap1"]["flanks"]}
        self.assertEqual(hap1_sides["left"]["requested_bp"], 10)
        self.assertEqual(hap1_sides["left"]["available_bp"], 3)
        self.assertEqual(hap1_sides["left"]["callable_bp"], 3)
        self.assertTrue(hap1_sides["left"]["boundary_truncated"])
        self.assertEqual(hap1_sides["right"]["available_bp"], 10)
        self.assertFalse(hap1_sides["right"]["boundary_truncated"])

    def test_one_sided_contig_boundary_flank_is_informative_and_side_typed(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        for row in rows:
            row["assembly1_body_start"] = "0"
        self.fixture.write_manifest(rows=rows)
        result = self.fixture.render()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        plan = next(
            value
            for _, value in self.fixture.rendered_json()
            if isinstance(value, dict)
            and value.get("schema_version") == "cnv_v2_plan_1"
        )
        hap1 = next(
            assembly
            for assembly in plan["target"]["assemblies"]
            if assembly["assembly_id"] == "hap1"
        )
        self.assertEqual(hap1["flank_completeness"], "ONE_SIDED_BOUNDARY_FLANK")
        sides = {row["side"]: row for row in hap1["flanks"]}
        self.assertEqual(
            sides["left"],
            {
                "side": "left",
                "start": 0,
                "end": 0,
                "requested_bp": 10,
                "available_bp": 0,
                "callable_bp": 0,
                "boundary_truncated": True,
            },
        )
        self.assertEqual(sides["right"]["available_bp"], 10)
        self.assertEqual(sides["right"]["callable_bp"], 10)
        self.assertFalse(sides["right"]["boundary_truncated"])
        flank_bed = next(self.fixture.run_root.rglob("outer_flanks.bed"))
        flank_names = [
            line.split("\t")[3]
            for line in flank_bed.read_text(encoding="utf-8").splitlines()
        ]
        self.assertNotIn("hap1_left_outer_flank", flank_names)
        self.assertIn("hap1_right_outer_flank", flank_names)

    def test_body_spanning_whole_contig_has_no_flank_and_is_rejected(self) -> None:
        rows = [dict(row) for row in self.fixture.rows]
        for row in rows:
            row["assembly1_body_start"] = "0"
            row["assembly1_body_end"] = "80"
        self.fixture.write_manifest(rows=rows)
        result = self.fixture.render()
        self.assert_rejected(result)
        self.assertIn("no callable outer flank", result.stderr.lower())

    def test_canonical_identity_changes_when_target_changes(self) -> None:
        first = self.fixture.render()
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        first_json = self.fixture.rendered_json()
        first_payload = json.dumps([value for _, value in first_json], sort_keys=True)
        first_paths = {str(path.relative_to(self.fixture.run_root)) for path, _ in first_json}
        first_target_ids = {
            value.get("target_id")
            for _, value in first_json
            if isinstance(value, dict) and value.get("schema_version") == "cnv_v2_plan_1"
        }

        other_root = Path(self.temporary.name) / "other"
        other = SyntheticRun(other_root)
        rows = [dict(row) for row in other.rows]
        for row in rows:
            row["locus_id"] = "locusB"
        other.write_manifest(rows=rows)
        second = other.render()
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        second_json = other.rendered_json()
        second_payload = json.dumps([value for _, value in second_json], sort_keys=True)
        second_paths = {str(path.relative_to(other.run_root)) for path, _ in second_json}
        second_target_ids = {
            value.get("target_id")
            for _, value in second_json
            if isinstance(value, dict) and value.get("schema_version") == "cnv_v2_plan_1"
        }
        self.assertNotEqual(first_payload, second_payload)
        self.assertNotEqual(first_paths, second_paths)
        self.assertNotEqual(first_target_ids, second_target_ids)

    def test_multi_target_manifest_renders_one_distinct_plan_per_target(self) -> None:
        self.fixture.write_manifest(rows=self.fixture.multi_target_rows())
        result = self.fixture.render()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        plans = [
            (path, value)
            for path, value in self.fixture.rendered_json()
            if isinstance(value, dict) and value.get("schema_version") == "cnv_v2_plan_1"
        ]
        self.assertEqual(len(plans), 2, plans)
        payloads = [json.dumps(value, sort_keys=True) for _, value in plans]
        self.assertTrue(any("sampleA" in value and "locusA" in value for value in payloads))
        self.assertTrue(any("sampleB" in value and "locusB" in value for value in payloads))
        self.assertEqual(len({str(path) for path, _ in plans}), 2)
        self.assertEqual(len({value.get("target_id") for _, value in plans}), 2)
        self.assertTrue(all(len(value.get("parts", [])) == 2 for _, value in plans))

        # The controller must consume every planner path, not treat multiline
        # output as one scalar filename.  Default mode remains scheduler-free.
        shim_dir = Path(self.temporary.name) / "multi-shims"
        shim_dir.mkdir()
        marker = Path(self.temporary.name) / "multi-external-called"
        for name in ("aws", "minimap2", "samtools", "sbatch", "srun"):
            shim = shim_dir / name
            shim.write_text(
                '#!/bin/sh\nprintf \'%s\\n\' "$0" >> "$CNV_TEST_MARKER"\nexit 97\n',
                encoding="utf-8",
            )
            shim.chmod(0o755)
        environment = os.environ.copy()
        environment["PATH"] = str(shim_dir) + os.pathsep + environment.get("PATH", "")
        environment["CNV_TEST_MARKER"] = str(marker)
        controller = subprocess.run(
            [
                "bash",
                str(CONTROLLER),
                "--manifest",
                str(self.fixture.manifest),
                "--expected-manifest-sha256",
                sha256(self.fixture.manifest),
                "--run-root",
                str(self.fixture.run_root),
            ],
            cwd=str(ROOT),
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )
        self.assertEqual(controller.returncode, 0, controller.stdout + controller.stderr)
        self.assertFalse(marker.exists(), "multi-target dry-run invoked an external tool")
        self.assertIn("sampleA", controller.stdout)
        self.assertIn("sampleB", controller.stdout)


class TestQnameClaimsAndNoClobber(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def inventory(self, name: str, payload: bytes) -> Path:
        path = self.root / name
        path.write_bytes(payload)
        return path

    def exact_plan(self, name: str = "exact") -> tuple[SyntheticRun, Path, dict[str, object]]:
        fixture = SyntheticRun(self.root / name)
        result = fixture.render()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        plans = [
            (path, value)
            for path, value in fixture.rendered_json()
            if isinstance(value, dict) and value.get("schema_version") == "cnv_v2_plan_1"
        ]
        self.assertEqual(len(plans), 1)
        plan_path, plan = plans[0]
        return fixture, plan_path, plan

    def materialize_prepared(self, plan: dict[str, object]) -> dict[str, object]:
        reference = plan["reference"]
        assert isinstance(reference, dict)
        for name in ("combined_fasta", "combined_fai", "combined_mmi", "bait_fasta"):
            path = Path(str(reference[name]))
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((name + "\n").encode("ascii"))
        artifacts = {
            name: CNV_PLAN.fingerprint(reference[name])
            for name in (
                "combined_fasta",
                "combined_fai",
                "combined_mmi",
                "bait_fasta",
                "full_window_bed",
                "core_body_union_bed",
                "outer_flanks_bed",
            )
        }
        payload = {
            "schema_version": "cnv_v2_prepared_1",
            "status": "PREPARED",
            "run_id": plan["run_id"],
            "target_id": plan["target_id"],
            "target_digest": plan["target_digest"],
            "plan_digest": plan["plan_digest"],
            "manifest_sha256": plan["manifest"]["sha256"],
            "artifacts": artifacts,
        }
        receipt_path = CNV_PLAN.write_receipt(plan["prepared_receipt"], payload)
        return json.loads(receipt_path.read_text(encoding="utf-8"))

    def tool_shims(self) -> tuple[Path, dict[str, str]]:
        shims = self.root / "shims"
        shims.mkdir(exist_ok=True)
        samtools = shims / "samtools"
        samtools.write_text(
            "#!/bin/sh\n"
            "command=$1\n"
            "case \"$command\" in\n"
            "  quickcheck|idxstats) exit 0 ;;\n"
            "  view)\n"
            "    last=\n"
            "    for item in \"$@\"; do last=$item; done\n"
            "    case \"$last\" in\n"
            "      *part_000001*) name=alpha ;;\n"
            "      *part_000002*) name=bravo ;;\n"
            "      *) name=unknown ;;\n"
            "    esac\n"
            "    printf '%s\\t0\\tref\\t1\\t60\\t1M\\t*\\t0\\t0\\tA\\t*\\n' \"$name\"\n"
            "    exit 0 ;;\n"
            "  *) exit 92 ;;\n"
            "esac\n",
            encoding="utf-8",
        )
        samtools.chmod(0o755)
        minimap2 = shims / "minimap2"
        minimap2.write_text("#!/bin/sh\nexit 93\n", encoding="utf-8")
        minimap2.chmod(0o755)
        environment = os.environ.copy()
        environment["PATH"] = str(shims) + os.pathsep + environment.get("PATH", "")
        return shims, environment

    def materialize_part_done(
        self,
        plan_path: Path,
        plan: dict[str, object],
        prepared: dict[str, object],
        part_index: int,
    ) -> Path:
        parts = plan["parts"]
        assert isinstance(parts, list)
        part = parts[part_index]
        assert isinstance(part, dict)
        for name in ("bam", "bai"):
            path = Path(str(part[name]))
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((name + "\n").encode("ascii"))
        inventory_name = "alpha" if part_index == 0 else "bravo"
        inventory = Path(str(part["primary_qnames"]))
        inventory.write_text(inventory_name + "\n", encoding="ascii")
        claim = CNV_PLAN.acquire_plan_claim(
            plan_path,
            "part",
            f"synthetic-worker-{part_index}",
            part_id=part["part_id"],
        )
        payload = {
            "schema_version": "cnv_v2_part_done_1",
            "status": "DONE",
            "run_id": plan["run_id"],
            "target_id": plan["target_id"],
            "target_digest": plan["target_digest"],
            "manifest_sha256": plan["manifest"]["sha256"],
            "plan_digest": plan["plan_digest"],
            "plan_file_sha256": sha256(plan_path),
            "part_id": part["part_id"],
            "part_index": part_index,
            "read_uri": part["read_uri"],
            "read_version": part["read_version"],
            "read_etag": part["read_etag"],
            "read_size": part["read_size"],
            "read_sha256": part["read_sha256"],
            "prepared_receipt_digest": prepared["receipt_digest"],
            "claim_digest": claim["claim_digest"],
            "prepared_artifacts": prepared["artifacts"],
            "outputs": {
                "bam": CNV_PLAN.fingerprint(part["bam"]),
                "bai": CNV_PLAN.fingerprint(part["bai"]),
                "primary_qnames": CNV_PLAN.inspect_qname_inventory(part["primary_qnames"]),
            },
        }
        return CNV_PLAN.write_receipt(part["receipt"], payload)

    def test_disjoint_inventory_union_is_streamed_and_receipt_ready(self) -> None:
        first = self.inventory("part1.qnames", b"alpha\ncharlie\n")
        second = self.inventory("part2.qnames", b"bravo\ndelta\n")
        union = CNV_PLAN.verify_qname_inventory_union([first, second])
        self.assertEqual(union["schema_version"], "cnv_v2_qname_union_1")
        self.assertEqual(union["verified_union_qname_count"], 4)
        self.assertEqual([item["qname_count"] for item in union["inventories"]], [2, 2])
        self.assertEqual(
            [item["sha256"] for item in union["inventories"]],
            [sha256(first), sha256(second)],
        )

    def test_within_and_cross_part_duplicate_qnames_refuse(self) -> None:
        duplicate_inside = self.inventory("within.qnames", b"alpha\nalpha\n")
        with self.assertRaisesRegex(CNV_PLAN.PlanError, "duplicate"):
            CNV_PLAN.inspect_qname_inventory(duplicate_inside)
        first = self.inventory("part1.qnames", b"alpha\ncharlie\n")
        second = self.inventory("part2.qnames", b"bravo\ncharlie\n")
        with self.assertRaisesRegex(CNV_PLAN.PlanError, "multiple parts"):
            CNV_PLAN.verify_qname_inventory_union([first, second])

    def test_tampered_unsorted_missing_and_malformed_inventories_refuse(self) -> None:
        cases = {
            "unsorted": b"bravo\nalpha\n",
            "no-terminal-lf": b"alpha",
            "whitespace": b"alpha beta\n",
        }
        for label, payload in cases.items():
            with self.subTest(label=label):
                path = self.inventory(label, payload)
                with self.assertRaises(CNV_PLAN.PlanError):
                    CNV_PLAN.inspect_qname_inventory(path)
        with self.assertRaisesRegex(CNV_PLAN.PlanError, "missing"):
            CNV_PLAN.inspect_qname_inventory(self.root / "missing.qnames")

    def test_empty_canonical_inventory_is_valid_and_bound(self) -> None:
        empty = self.inventory("empty.qnames", b"")
        record = CNV_PLAN.inspect_qname_inventory(empty)
        self.assertEqual(record["qname_count"], 0)
        self.assertEqual(record["size"], 0)
        self.assertEqual(record["sha256"], hashlib.sha256(b"").hexdigest())

    def test_bam_inconsistent_inventory_refuses_on_fresh_rederivation(self) -> None:
        fake_samtools = self.root / "samtools"
        fake_samtools.write_text(
            "#!/bin/sh\n"
            "test \"$1\" = view || exit 91\n"
            "printf 'alpha\\t0\\tref\\t1\\t60\\t1M\\t*\\t0\\t0\\tA\\t*\\n'\n"
            "printf 'charlie\\t0\\tref\\t2\\t60\\t1M\\t*\\t0\\t0\\tC\\t*\\n'\n",
            encoding="utf-8",
        )
        fake_samtools.chmod(0o755)
        bam = self.root / "synthetic.bam"
        bam.write_bytes(b"placeholder\n")
        matching = self.inventory("matching.qnames", b"alpha\ncharlie\n")
        CNV_PLAN.verify_inventory_matches_bam(matching, bam, samtools=str(fake_samtools))
        inconsistent = self.inventory("inconsistent.qnames", b"alpha\ndelta\n")
        with self.assertRaisesRegex(CNV_PLAN.PlanError, "differs"):
            CNV_PLAN.verify_inventory_matches_bam(inconsistent, bam, samtools=str(fake_samtools))

    def test_competing_worker_and_analysis_claims_fail_closed(self) -> None:
        for kind in ("part", "analysis"):
            with self.subTest(kind=kind):
                claim = self.root / f"{kind}.CLAIM.json"
                payload = {
                    "status": "CLAIMED",
                    "kind": kind,
                    "run_id": "run_exact",
                    "target_id": "target_exact",
                    "plan_digest": "a" * 64,
                    "claimant_token": f"first-{kind}",
                }
                first = CNV_PLAN.acquire_create_only_claim(claim, payload)
                self.assertEqual(CNV_PLAN.load_claim(claim), first)
                with self.assertRaisesRegex(CNV_PLAN.PlanError, "will not be stolen"):
                    CNV_PLAN.acquire_create_only_claim(
                        claim,
                        {**payload, "claimant_token": f"competing-{kind}"},
                    )
                self.assertTrue(claim.is_file(), "stale/failed claim was silently removed")

    def test_dangling_symlink_claim_is_occupied_without_following(self) -> None:
        missing_target = self.root / "must-not-be-created.json"
        claim = self.root / "literal.CLAIM.json"
        claim.symlink_to(missing_target)
        with self.assertRaisesRegex(CNV_PLAN.PlanError, "will not be stolen"):
            CNV_PLAN.acquire_create_only_claim(
                claim,
                {
                    "status": "CLAIMED",
                    "kind": "part",
                    "run_id": "run_exact",
                    "target_id": "target_exact",
                    "plan_digest": "a" * 64,
                    "claimant_token": "worker",
                },
            )
        self.assertTrue(claim.is_symlink())
        self.assertFalse(missing_target.exists())

    def test_parent_symlink_escape_refuses_load_and_claim_without_outside_write(self) -> None:
        _, plan_path, plan = self.exact_plan("parent-symlink-escape")
        target_dir = Path(str(plan["target_dir"]))
        parts_dir = target_dir / "parts"
        outside = self.root / "outside-target"
        outside.mkdir()
        self.assertFalse(os.path.lexists(parts_dir))
        parts_dir.symlink_to(outside, target_is_directory=True)

        with self.assertRaisesRegex(CNV_PLAN.PlanError, "outside plan target_dir"):
            CNV_PLAN.load_plan(plan_path)
        with self.assertRaisesRegex(CNV_PLAN.PlanError, "outside plan target_dir"):
            CNV_PLAN.acquire_plan_claim(
                plan,
                "part",
                "must-not-escape",
                part_id=1,
            )

        self.assertTrue(parts_dir.is_symlink())
        self.assertEqual(list(outside.iterdir()), [], "claim escaped through parent symlink")

    def test_lexical_dot_dot_escape_in_plan_refuses(self) -> None:
        _, plan_path, plan = self.exact_plan("lexical-dot-dot-escape")
        parts = plan["parts"]
        assert isinstance(parts, list)
        part = parts[0]
        assert isinstance(part, dict)
        part["claim"] = str(
            Path(str(plan["target_dir"]))
            / "parts"
            / ".."
            / ".."
            / "escaped.CLAIM.json"
        )
        plan["plan_digest"] = CNV_PLAN._plan_digest(plan)
        plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        with self.assertRaisesRegex(CNV_PLAN.PlanError, "outside plan target_dir"):
            CNV_PLAN.load_plan(plan_path, verify_files=False)
        self.assertFalse((Path(str(plan["target_dir"])).parent / "escaped.CLAIM.json").exists())

    def test_rendered_plan_declares_exact_claim_and_inventory_paths(self) -> None:
        fixture = SyntheticRun(self.root / "planned")
        result = fixture.render()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        plans = [
            (path, value)
            for path, value in fixture.rendered_json()
            if isinstance(value, dict) and value.get("schema_version") == "cnv_v2_plan_1"
        ]
        self.assertEqual(len(plans), 1)
        plan_path, plan = plans[0]
        self.assertTrue(Path(plan["analysis_claim"]).is_absolute())
        for part in plan["parts"]:
            self.assertTrue(Path(part["claim"]).is_absolute())
            self.assertTrue(Path(part["primary_qnames"]).is_absolute())
        part_claim = CNV_PLAN.acquire_plan_claim(plan_path, "part", "synthetic-worker", part_id=1)
        CNV_PLAN.validate_plan_claim(plan, "part", part_claim["claim_digest"], part_id=1)
        analysis_claim = CNV_PLAN.acquire_plan_claim(plan_path, "analysis", "synthetic-analysis")
        CNV_PLAN.validate_plan_claim(plan, "analysis", analysis_claim["claim_digest"])

    def test_create_only_publication_cannot_overwrite_canonical_output(self) -> None:
        staged = self.root / "stage.one"
        canonical = self.root / "canonical.out"
        staged.write_bytes(b"first\n")
        CNV_PLAN.publish_create_only(staged, canonical)
        self.assertFalse(staged.exists())
        self.assertEqual(canonical.read_bytes(), b"first\n")
        competitor = self.root / "stage.two"
        competitor.write_bytes(b"second\n")
        with self.assertRaisesRegex(CNV_PLAN.PlanError, "will not be overwritten"):
            CNV_PLAN.publish_create_only(competitor, canonical)
        self.assertEqual(canonical.read_bytes(), b"first\n")
        self.assertEqual(competitor.read_bytes(), b"second\n")

    def test_create_only_publication_rejects_dangling_symlink_destination(self) -> None:
        staged = self.root / "staged"
        staged.write_bytes(b"candidate\n")
        missing_target = self.root / "must-not-be-published"
        canonical = self.root / "canonical"
        canonical.symlink_to(missing_target)
        with self.assertRaisesRegex(CNV_PLAN.PlanError, "will not be overwritten"):
            CNV_PLAN.publish_create_only(staged, canonical)
        self.assertEqual(staged.read_bytes(), b"candidate\n")
        self.assertTrue(canonical.is_symlink())
        self.assertFalse(missing_target.exists())

    def test_missing_digest_done_refuses_public_validate_and_worker_resume(self) -> None:
        _, plan_path, plan = self.exact_plan("missing-digest")
        prepared = self.materialize_prepared(plan)
        receipt_path = self.materialize_part_done(plan_path, plan, prepared, 0)
        _, environment = self.tool_shims()
        valid = subprocess.run(
            [
                sys.executable,
                str(PLANNER),
                "validate-part",
                "--plan",
                str(plan_path),
                "--part-id",
                "1",
                "--samtools",
                "samtools",
            ],
            cwd=str(ROOT),
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )
        self.assertEqual(valid.returncode, 0, valid.stdout + valid.stderr)
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        receipt.pop("receipt_digest")
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        invalid = subprocess.run(
            [
                sys.executable,
                str(PLANNER),
                "validate-part",
                "--plan",
                str(plan_path),
                "--part-id",
                "1",
                "--samtools",
                "samtools",
            ],
            cwd=str(ROOT),
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )
        self.assertNotEqual(invalid.returncode, 0)
        self.assertIn("receipt digest mismatch", invalid.stderr.lower())
        worker = subprocess.run(
            ["bash", str(WORKER), "--plan", str(plan_path), "--part-index", "0", "--verify-only"],
            cwd=str(ROOT),
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )
        self.assertNotEqual(worker.returncode, 0)
        self.assertIn("stale or invalid", worker.stderr.lower())

    def test_worker_state_treats_dangling_symlink_as_partial(self) -> None:
        _, plan_path, plan = self.exact_plan("worker-dangling")
        self.materialize_prepared(plan)
        parts = plan["parts"]
        assert isinstance(parts, list) and isinstance(parts[0], dict)
        bam = Path(str(parts[0]["bam"]))
        bam.parent.mkdir(parents=True, exist_ok=True)
        bam.symlink_to(bam.parent / "missing-target.bam")
        _, environment = self.tool_shims()
        result = subprocess.run(
            ["bash", str(WORKER), "--plan", str(plan_path), "--part-index", "0", "--verify-only"],
            cwd=str(ROOT),
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("partial prior state exists", result.stderr.lower())

    def test_analysis_state_treats_dangling_symlink_as_partial(self) -> None:
        _, plan_path, plan = self.exact_plan("analysis-dangling")
        prepared = self.materialize_prepared(plan)
        self.materialize_part_done(plan_path, plan, prepared, 0)
        self.materialize_part_done(plan_path, plan, prepared, 1)
        outputs = plan["outputs"]
        assert isinstance(outputs, dict)
        partial = Path(str(outputs["raw_depth_tsv"]))
        partial.parent.mkdir(parents=True, exist_ok=True)
        partial.symlink_to(partial.parent / "missing-depth.tsv")
        _, environment = self.tool_shims()
        result = subprocess.run(
            ["bash", str(ANALYSIS), "--plan", str(plan_path), "--verify-only"],
            cwd=str(ROOT),
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("partial analysis state exists", result.stderr.lower())

    def test_worker_and_analysis_publish_done_last_without_overwriting_move(self) -> None:
        for script, receipt_variable in ((WORKER, "$PART_RECEIPT"), (ANALYSIS, "$ANALYSIS_RECEIPT")):
            with self.subTest(script=script.name):
                text = read_text(script)
                publish_lines = [line for line in command_lines(text) if " publish " in f" {line} "]
                self.assertTrue(publish_lines)
                self.assertIn(receipt_variable, publish_lines[-1])
                self.assertFalse(any(re.search(r"(?:^|\s)mv\s", line) for line in command_lines(text)))
                self.assertRegex(text.lower(), r"claim.*(?:already exists|stale|acquir)")
        self.assertIn("partial prior state exists", read_text(WORKER).lower())
        self.assertIn("partial analysis state exists", read_text(ANALYSIS).lower())
        self.assertIn('"qname_union":qname_union', read_text(ANALYSIS).replace(" ", ""))


class TestWorkerSemantics(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.text = read_text(WORKER)
        cls.lines = command_lines(cls.text)

    def test_bait_capture_has_no_cap_or_downsampling(self) -> None:
        lower = self.text.lower()
        self.assertIn("bait", lower)
        self.assertNotRegex(lower, r"\b(?:head|shuf)\b[^\n]*(?:read|name|fastq|bait)")
        self.assertNotRegex(lower, r"(?:max_reads|read_cap|downsampl|subsampl)")
        self.assertRegex(lower, r"(?:aligned|alignment|paf).*(?:bp|bases|length)|(?:bp|bases|length).*(?:aligned|alignment|paf)")

    def test_captured_reads_are_aligned_once_to_combined_reference(self) -> None:
        align_lines = [
            line
            for line in self.lines
            if "minimap2" in line and re.search(r"(?:^|\s)-a(?:x)?(?:\s|\"|'|$)", line)
        ]
        self.assertEqual(align_lines, [align_lines[0]] if align_lines else [], self.text)
        self.assertEqual(len(align_lines), 1, self.text)
        self.assertRegex(align_lines[0].lower(), r"combined|diploid")
        self.assertNotRegex(self.text, r"(?m)for\s+.*assembl.*;\s*do[\s\S]{0,400}minimap2\s+.*-a")

    def test_bam_filter_uses_full_analysis_window_not_core(self) -> None:
        view_lines = [line for line in self.lines if re.search(r"\bsamtools\s+view\b", line)]
        window_filters = [line for line in view_lines if "-L" in line and "full_window" in line.lower()]
        self.assertTrue(window_filters, self.text)
        self.assertFalse(any("core_body" in line.lower() and "-L" in line for line in view_lines), self.text)

    def test_primary_only_filter_uses_uppercase_exclusion_flag(self) -> None:
        filtered = [
            line
            for line in self.lines
            if re.search(r"\bsamtools\s+view\b", line) and "0x900" in line
        ]
        self.assertGreaterEqual(len(filtered), 4, self.text)
        self.assertTrue(all(re.search(r"(?:^|\s)-F\s+0x900(?:\s|$)", line) for line in filtered), self.text)
        self.assertFalse(any(re.search(r"(?:^|\s)-f\s+0x900(?:\s|$)", line) for line in filtered), self.text)

    def test_download_is_private_and_immutable_identity_is_checked(self) -> None:
        lower = self.text.lower()
        self.assertRegex(self.text, r"mktemp\s+(?:[^\n]*\s)?-d\b|-d\s+[^\n]*mktemp")
        self.assertIn("head-object", lower)
        for token in ("version", "etag", "contentlength", "sha256", "quickcheck"):
            with self.subTest(token=token):
                self.assertIn(token, lower.replace("_", ""))

    def test_resume_requires_exact_done_identity_and_bam_integrity(self) -> None:
        lower = self.text.lower()
        for token in ("done", "manifest_sha256", "quickcheck"):
            with self.subTest(token=token):
                self.assertIn(token, lower)
        self.assertRegex(lower, r"target_(?:identity|digest)")
        self.assertRegex(lower, r"(?:bai|index)")
        self.assertRegex(lower, r"resume|existing|already")
        self.assertIn("receipt_digest", lower)

    def test_worker_creates_only_the_canonical_part_output_directory(self) -> None:
        self.assertRegex(
            self.text.lower(),
            r"mkdir\s+-p\s+--\s+(?:\"\$part_output_dir\"|\"\$\(dirname[^\n]*\$part_bam[^\n]*\)\")",
        )
        self.assertRegex(self.text.lower(), r"part_(?:output_)?dir|dirname[^\n]*part_bam")

    def test_done_binds_prepared_receipt_and_all_seven_fingerprints(self) -> None:
        lower = self.text.lower()
        self.assertRegex(lower, r"prepared_receipt_(?:sha256|digest)|prepared[^\n]{0,80}receipt_digest")
        for artifact in (
            "combined_fasta",
            "combined_fai",
            "combined_mmi",
            "bait_fasta",
            "full_window_bed",
            "core_body_union_bed",
            "outer_flanks_bed",
        ):
            with self.subTest(artifact=artifact):
                self.assertGreaterEqual(lower.count(artifact), 2, self.text)


class TestReferencePreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.text = read_text(PLANNER)
        cls.lower = cls.text.lower()

    def test_one_run_scoped_combined_index_is_built(self) -> None:
        flags = list(re.finditer(r"[\"']-d[\"']", self.text))
        self.assertEqual(len(flags), 1, self.text)
        match = flags[0]
        context = self.text[max(0, match.start() - 300) : match.end() + 300].lower()
        self.assertIn("minimap2", context)
        self.assertRegex(context, r"combined|diploid")
        self.assertIn(".mmi", self.lower)
        self.assertRegex(self.lower, r"run_dir|run_id")

    def test_combined_contig_namespace_is_collision_checked(self) -> None:
        self.assertRegex(self.lower, r"combined_contig|combined.*contig|contig.*combined")
        self.assertRegex(self.lower, r"duplicate|distinct|unique|collision")
        self.assertRegex(self.text, r"__")

    def test_prepared_receipt_is_identity_bound_and_immutable(self) -> None:
        for token in ("prepared", "manifest_sha256", "sha256"):
            with self.subTest(token=token):
                self.assertIn(token, self.lower)
        self.assertRegex(self.lower, r"target_(?:identity|digest)")
        self.assertRegex(self.lower, r"exist|resume|already")
        self.assertRegex(self.lower, r"mismatch|refus|fail|raise")

    def test_prepared_receipt_closes_the_same_reference_and_bed_artifacts(self) -> None:
        worker = read_text(WORKER).lower()
        planner_validation = self.lower[
            self.lower.index("def validate_prepared") : self.lower.index("def _append_prefixed_fasta")
        ]
        planner_preparation = self.lower[
            self.lower.index("def prepare_reference") : self.lower.index("def get_part")
        ]
        worker_validation = worker[
            worker.index("verify_prepared_receipt") : worker.index("verify_part_receipt")
        ]
        artifacts = (
            "combined_fasta",
            "combined_fai",
            "combined_mmi",
            "bait_fasta",
            "full_window_bed",
            "core_body_union_bed",
            "outer_flanks_bed",
        )
        for artifact in artifacts:
            with self.subTest(artifact=artifact, location="planner validation"):
                self.assertIn(artifact, planner_validation)
            with self.subTest(artifact=artifact, location="planner preparation"):
                self.assertIn(artifact, planner_preparation)
            with self.subTest(artifact=artifact, location="worker validation"):
                self.assertIn(artifact, worker_validation)

    def test_part_receipt_schema_matches_worker_outputs_map(self) -> None:
        worker = read_text(WORKER).lower()
        planner_part = self.lower[
            self.lower.index("def validate_part_receipt") : self.lower.index("def _parser")
        ]
        worker_part = worker[worker.index("verify_part_receipt") :]
        self.assertIn('"outputs"', planner_part)
        self.assertIn('"outputs"', worker_part)
        self.assertNotIn('receipt.get("artifacts")', planner_part)
        for artifact in ("bam", "bai"):
            with self.subTest(artifact=artifact):
                self.assertIn(f'"{artifact}"', planner_part)
                self.assertIn(f'"{artifact}"', worker_part)
        self.assertRegex(planner_part, r"prepared_receipt_(?:sha256|digest)|prepared[^\n]{0,80}receipt_digest")


class TestAnalysisAndReporting(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.text = read_text(ANALYSIS)
        cls.lower = cls.text.lower()
        cls.lines = command_lines(cls.text)
        cls.plotter = read_text(PLOTTER)
        cls.readme = read_text(README)
        cls.controller = read_text(CONTROLLER)

    def test_merge_is_exactly_the_declared_parts(self) -> None:
        merge_lines = [line for line in self.lines if re.search(r"\bsamtools\s+merge\b", line)]
        self.assertTrue(merge_lines, self.text)
        self.assertFalse(any("*" in line or "?" in line for line in merge_lines), self.text)
        self.assertIn("expected_part", self.lower)
        self.assertIn("part_id", self.lower)
        self.assertIn("done", self.lower)

    def test_depth_filters_are_raw_and_mapq_10(self) -> None:
        depth_lines = [line for line in self.lines if re.search(r"\bsamtools\s+depth\b", line)]
        self.assertGreaterEqual(len(depth_lines), 2, self.text)
        self.assertTrue(any(re.search(r"(?:^|\s)-Q\s+10(?:\s|$)", line) for line in depth_lines), self.text)
        self.assertTrue(any(not re.search(r"(?:^|\s)-Q(?:\s|=)", line) for line in depth_lines), self.text)
        self.assertFalse(any(re.search(r"(?:^|\s)-q\s+10(?:\s|$)", line) for line in depth_lines), self.text)
        self.assertTrue(all("full_window" in line.lower() for line in depth_lines), self.text)

    def test_boundary_flank_analysis_is_side_aware_and_not_two_flank_required(self) -> None:
        self.assertIn("ONE_SIDED_BOUNDARY_FLANK", self.text)
        self.assertIn("flank_sides", self.text)
        for token in (
            "requested_bp",
            "available_bp",
            "callable_bp",
            "boundary_truncated",
            "positive_depth_bp",
        ):
            with self.subTest(token=token):
                self.assertIn(token, self.text)
        self.assertNotRegex(self.text, r"len\(fi\)\s*!=\s*2")
        self.assertRegex(
            self.text,
            r"len\(fi\)\s*!=\s*len\(planned_flanks\)",
        )
        self.assertRegex(
            self.text.lower(),
            r"missing boundary side[^\\n]*(?:neither padded nor imputed)|"
            r"(?:neither padded nor imputed)[^\\n]*missing boundary side",
        )
        self.assertRegex(
            self.text.lower(),
            r"body depth alone never rescues a no-flank observation",
        )

    def test_plotter_is_static_and_has_no_create_remove_race(self) -> None:
        self.assertRegex(self.text, r"python3\s+\"\$PLOTTER\"")
        self.assertRegex(self.text.lower(), r"plotter_sha256|plotter[^\n]{0,80}sha256")
        self.assertNotRegex(self.text, r"(?:cat|printf|echo)[^\n]*plot_cnv_depth_v2\.py")
        self.assertNotRegex(self.text, r"\brm\b[^\n]*plot_cnv_depth_v2\.py")
        self.assertNotRegex(self.text, r"\b(?:cp|mv)\b[^\n]*plot_cnv_depth_v2\.py")
        self.assertIn("sha256", self.lower)

    def test_done_follows_nonzero_target_and_valid_pdf_gates(self) -> None:
        done_position = self.lower.rfind("done")
        self.assertGreater(done_position, 0)
        nonzero_position = max(
            self.lower.find("nonzero"),
            self.lower.find("non-zero"),
            self.lower.find("target_depth"),
            self.lower.find("body_depth"),
        )
        self.assertGreater(nonzero_position, -1, self.text)
        self.assertLess(nonzero_position, done_position, self.text)
        pdf_position = max(self.lower.find("%pdf"), self.lower.find("pdf"))
        self.assertGreater(pdf_position, -1, self.text)
        self.assertLess(pdf_position, done_position, self.text)
        self.assertRegex(self.lower, r"pdf[^\n]*(?:valid|magic|eof|nonempty|non-empty|size)|(?:valid|magic|eof|nonempty|non-empty|size)[^\n]*pdf")

    def test_optional_nucfreq_is_separate_from_coverage_done(self) -> None:
        contract = (self.text + "\n" + self.controller + "\n" + self.readme).lower()
        self.assertRegex(contract, r"optional|not (?:submitted|required)|separate")
        self.assertRegex(contract, r"nucfreq|nucleotide-frequency")
        self.assertRegex(
            contract,
            r"(?:nucfreq|nucleotide-frequency)[^\n]*(?:separate|receipt|independent|not submitted|not required)|"
            r"(?:separate|receipt|independent|not submitted|not required)[^\n]*(?:nucfreq|nucleotide-frequency)",
        )
        self.assertNotRegex(self.lower, r"(?:coverage_done|done_receipt)[^\n]*&&[^\n]*nucfreq")

    def test_scientific_labels_are_descriptive_and_do_not_overclaim(self) -> None:
        reporting = (self.plotter + "\n" + self.readme).lower()
        required = (
            "raw depth",
            "mapq >= 10",
            "pooled non-overlapping outer target-flank baseline",
            "descriptive",
            "own-flank",
            "not estimated",
            "diagnostic",
        )
        for label in required:
            with self.subTest(label=label):
                self.assertIn(label, reporting)
        for forbidden in ("deduplicated", "diploid baseline", "haploid copies"):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, reporting)
        self.assertRegex(reporting, r"assembly[- ]support")
        self.assertNotRegex(reporting, r"(?:unique|uniquely mapped) (?:depth|coverage|reads?)")
        self.assertNotRegex(reporting, r"(?:automatic|inferred|estimated) copy[- ]number\s*[:=]\s*[0-9]")

    def test_own_flank_is_one_copy_scale_without_division_by_two(self) -> None:
        reporting = (self.text + "\n" + self.plotter + "\n" + self.readme).lower()
        self.assertRegex(reporting, r"one[- ]copy depth scale")
        self.assertRegex(reporting, r"body\s*/\s*own-flank|body[^\n]{0,80}own-flank[^\n]{0,80}normali")
        self.assertNotRegex(reporting, r"body[^\n]{0,120}(?:/\s*2|divide[^\n]{0,20}(?:by\s+)?2)")
        self.assertRegex(reporting, r"pooled[^\n]{0,160}(?:cannot|must not|does not)[^\n]{0,100}(?:rescue|replace)")

    def test_copy_number_is_explicitly_not_estimated(self) -> None:
        combined = (self.text + "\n" + self.plotter + "\n" + self.readme).lower()
        self.assertRegex(combined, r"read_inferred_copy_number[\"']?\s*[:=]\s*[\"']?not_estimated")
        self.assertRegex(combined, r"coverage[^\n]{0,120}(?:diagnostic|assembly-support)")

    def test_kcon_and_mapping_quality_warnings_remain_visible(self) -> None:
        combined = (self.text + "\n" + self.plotter + "\n" + self.readme).lower()
        self.assertIn("kcon", combined)
        self.assertRegex(combined, r"assembly[- ]resolved[^\n]{0,80}kcon|kcon[^\n]{0,80}assembly[- ]resolved")
        self.assertRegex(combined, r"low[- ]mapq")
        self.assertRegex(combined, r"ambigu(?:ity|ous)")

    def test_kcon_ownership_is_exact_body_overlap_and_empty_paf_is_allowed(self) -> None:
        # Prefix-only assignment can incorrectly attribute an alignment to a
        # different contig from the same assembly.  Ownership is the exact
        # planner-bound combined contig with positive half-open body overlap.
        self.assertNotRegex(
            self.text,
            r"\.startswith\(a\[\"assembly_id\"\]\s*\+\s*\"__\"\)",
        )
        marker = self.text.index("KCON")
        kcon_logic = self.text[marker:]
        for token in ("combined_contig", "body_start", "body_end", "f[7]", "f[8]"):
            with self.subTest(token=token):
                self.assertIn(token, kcon_logic)
        self.assertRegex(
            kcon_logic.lower(),
            r"(?:overlap|max\([^\n]*body_start|min\([^\n]*body_end)",
        )
        self.assertRegex(
            kcon_logic,
            r"f\[5\]\s*(?:==|!=)\s*a\[\"combined_contig\"\]|"
            r"a\[\"combined_contig\"\]\s*(?:==|!=)\s*f\[5\]",
        )

        # No KCON hit is a reportable diagnostic outcome, not a coverage-stage
        # failure.  The empty PAF is retained and produces per-assembly warnings.
        self.assertNotRegex(self.text, r"\[\[[^\n]*-s\s+\"\$kcon\"")
        self.assertRegex(self.text, r"No KCON alignment was reported")


class TestPlotterOffline(unittest.TestCase):
    @staticmethod
    def _run_plotter(
        depth: Path,
        output: Path,
        metric: str,
        summary: Path | None = None,
    ) -> subprocess.CompletedProcess[str]:
        command = [
            sys.executable,
            str(PLOTTER),
            "--depth-tsv",
            str(depth),
            "--output-pdf",
            str(output),
            "--metric",
            metric,
            "--title",
            "Synthetic depth diagnostic",
        ]
        if summary is not None:
            command.extend(("--summary-json", str(summary)))
        return subprocess.run(
            command,
            cwd=str(ROOT),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )

    def test_plotter_generates_a_structurally_valid_pdf(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            depth = root / "depth.tsv"
            depth.write_text(
                "hap1__ctgA\t1\t2\n"
                "hap1__ctgA\t2\t3\n"
                "hap2__ctgB\t71\t4\n"
                "hap2__ctgB\t72\t5\n",
                encoding="ascii",
            )
            output = root / "depth.pdf"
            result = self._run_plotter(depth, output, "raw")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            payload = output.read_bytes()
            self.assertTrue(payload.startswith(b"%PDF-"))
            self.assertTrue(payload.rstrip().endswith(b"%%EOF"))
            self.assertGreater(len(payload), 200)

    def test_summary_text_is_selected_for_the_requested_metric(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            depth = root / "depth.tsv"
            depth.write_text(
                "hap1__ctgA\t1\t2\n"
                "hap1__ctgA\t2\t3\n"
                "hap2__ctgB\t1\t4\n"
                "hap2__ctgB\t2\t5\n",
                encoding="ascii",
            )
            summary = root / "summary.json"
            summary.write_text(
                json.dumps(
                    {
                        "read_inferred_copy_number": "not_estimated",
                        "pooled_outer_flank_baselines": {
                            "raw_zero_inclusive_median_depth": 13,
                            "mapq10_zero_inclusive_median_depth": 7,
                            "descriptive_only": True,
                            "may_rescue_failed_haplotype": False,
                        },
                        "assembly_summaries": [
                            {
                                "assembly_id": "hap1",
                                "raw": {
                                    "own_outer_flank_zero_inclusive_median_depth": 11,
                                    "body_own_flank_normalized_depth": 2.5,
                                },
                                "mapq10": {
                                    "own_outer_flank_zero_inclusive_median_depth": 5,
                                    "body_own_flank_normalized_depth": 1.25,
                                },
                            },
                            {
                                "assembly_id": "hap2",
                                "raw": {
                                    "own_outer_flank_zero_inclusive_median_depth": 17,
                                    "body_own_flank_normalized_depth": 3.5,
                                },
                                "mapq10": {
                                    "own_outer_flank_zero_inclusive_median_depth": 3,
                                    "body_own_flank_normalized_depth": 0.75,
                                },
                            },
                        ],
                        "warnings": [],
                    },
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            raw_pdf = root / "raw.pdf"
            mapq_pdf = root / "mapq10.pdf"
            raw = self._run_plotter(depth, raw_pdf, "raw", summary)
            mapq = self._run_plotter(depth, mapq_pdf, "mapq10", summary)
            self.assertEqual(raw.returncode, 0, raw.stdout + raw.stderr)
            self.assertEqual(mapq.returncode, 0, mapq.stdout + mapq.stderr)

            raw_text = raw_pdf.read_bytes().decode("latin-1")
            mapq_text = mapq_pdf.read_bytes().decode("latin-1")
            for expected in (
                "Raw depth pooled outer target-flank baseline: 13 depth",
                "hap1; own-flank scale=11; body/own-flank=2.5",
                "hap2; own-flank scale=17; body/own-flank=3.5",
            ):
                with self.subTest(metric="raw", expected=expected):
                    self.assertIn(expected, raw_text)
                    self.assertNotIn(expected, mapq_text)
            for expected in (
                "MAPQ >= 10 depth pooled outer target-flank baseline: 7 depth",
                "hap1; own-flank scale=5; body/own-flank=1.25",
                "hap2; own-flank scale=3; body/own-flank=0.75",
            ):
                with self.subTest(metric="mapq10", expected=expected):
                    self.assertIn(expected, mapq_text)
                    self.assertNotIn(expected, raw_text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
