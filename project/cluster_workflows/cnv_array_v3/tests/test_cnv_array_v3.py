#!/usr/bin/env python3
"""Focused offline tests for the Rocky9 structural CNV/array closure."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import pathlib
import shutil
import subprocess
import tempfile
import unittest


PACKAGE = pathlib.Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parents[1]
MANIFEST_PATH = PACKAGE / "frozen" / "cnv_array_exact4_manifest.v1.json"


def load_module(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


plan = load_module("cnv_array_plan_v3", PACKAGE / "src" / "cnv_array_plan.py")
builder = load_module("build_cnv_array_manifest_v3", PACKAGE / "build_manifest.py")


def redigest(manifest: dict) -> dict:
    manifest = copy.deepcopy(manifest)
    manifest["plan_sha256"] = ""
    manifest["plan_sha256"] = plan.digest_without(manifest, "plan_sha256")
    return manifest


class ManifestTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_manifest_is_deterministic_and_valid(self) -> None:
        self.assertEqual(plan.canonical_bytes(builder.build()), MANIFEST_PATH.read_bytes())
        validated = plan.load_manifest(MANIFEST_PATH)
        self.assertEqual([u["sample_id"] for u in validated["units"]], list(plan.SAMPLES))
        self.assertEqual(sum(len(u["haplotypes"]) for u in validated["units"]), 8)

    def test_exact_copy_truth_and_alias_are_frozen(self) -> None:
        observed = {
            (u["sample_id"], h["assembly_haplotype"], h["truth_haplotype"]): h["physical_copy_number"]
            for u in self.manifest["units"] for h in u["haplotypes"]
        }
        self.assertEqual(observed[("HG02178", "hap1", "h1")], 1)
        self.assertEqual(observed[("HG02178", "hap2", "h2")], 4)
        self.assertEqual(sorted(observed.values()), [1, 1, 1, 1, 3, 3, 3, 4])

    def test_exclusion_sets_match_staging_receipt(self) -> None:
        staging = json.loads((ROOT / "inputs" / "cluster_execution_plans" / "cnv_array_staging_plan.json").read_text(encoding="utf-8"))
        self.assertEqual(self.manifest["selection"]["accepted_historical_samples_excluded"], staging["historical_results"]["complete_samples"])
        expected_missing_source = [x for x in staging["historical_results"]["incomplete_samples"] if x not in plan.SAMPLES]
        self.assertEqual(self.manifest["selection"]["missing_source_samples_excluded"], expected_missing_source)

    def test_forbidden_execution_actions_refuse(self) -> None:
        for key in self.manifest["execution_policy"]:
            bad = copy.deepcopy(self.manifest)
            bad["execution_policy"][key] = True
            with self.assertRaises(plan.PlanError, msg=key):
                plan.validate_manifest(redigest(bad))

    def test_afterany_refuses(self) -> None:
        bad = copy.deepcopy(self.manifest)
        bad["scheduling"]["closure"]["dependency"] = "afterany:unit_array"
        with self.assertRaises(plan.PlanError):
            plan.validate_manifest(redigest(bad))

    def test_copy_truth_mutation_refuses(self) -> None:
        bad = copy.deepcopy(self.manifest)
        bad["units"][0]["haplotypes"][1]["physical_copy_number"] = 2
        with self.assertRaises(plan.PlanError):
            plan.validate_manifest(redigest(bad))

    def test_destination_cannot_enter_protected_results(self) -> None:
        bad = copy.deepcopy(self.manifest)
        bad["units"][0]["outputs"]["evidence_tsv"] = bad["storage"]["protected_roots"][1] + "/overwrite.tsv"
        with self.assertRaises(plan.PlanError):
            plan.validate_manifest(redigest(bad))

    def test_scheduler_recipe_is_receipt_conditioned(self) -> None:
        rendered = plan.render_sbatch(MANIFEST_PATH)
        self.assertIn("--array=0-3%4", rendered)
        self.assertNotIn("--dependency", rendered)
        self.assertIn("03_cnv_array_receipt_controller.sh", rendered)
        self.assertEqual(rendered.count("HML2_CLUSTER_PYTHON="), 2)
        self.assertIn(plan.CLUSTER_PYTHON, rendered)
        self.assertNotIn("afterany", rendered.lower())
        self.assertEqual(rendered.count("/usr/bin/sbatch"), 2)
        self.assertIn("/run_20260718_003/bundle/src/01_cnv_array_unit.sh", rendered)
        self.assertNotIn(str(PACKAGE.resolve()), rendered)

    def test_shared_cpython311_is_manifest_bound(self) -> None:
        runtime = self.manifest["runtime"]
        selected = runtime["required_environment"]["HML2_CLUSTER_PYTHON"]
        self.assertEqual(selected, plan.CLUSTER_PYTHON)
        self.assertEqual(runtime["tools"]["python"]["path"], selected)
        for name in (
                "01_cnv_array_unit.sh", "02_cnv_array_closure.sh",
                "03_cnv_array_receipt_controller.sh"):
            source = (PACKAGE / "src" / name).read_text(encoding="utf-8")
            self.assertIn("HML2_CLUSTER_PYTHON", source)
            self.assertNotIn("/usr/bin/python3", source)

    def test_spooled_script_resolves_planner_only_from_absolute_manifest_binding(self) -> None:
        manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as script_tmp, tempfile.TemporaryDirectory() as cwd_tmp:
            script_dir = pathlib.Path(script_tmp) / "unrelated-slurm-spool"
            script_dir.mkdir()
            copied = script_dir / "job-script"
            shutil.copyfile(PACKAGE / "src" / "01_cnv_array_unit.sh", copied)
            copied.chmod(0o755)
            manifest_path = pathlib.Path(script_tmp) / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            proc = subprocess.run(
                [str(copied), "--dry-run", str(manifest_path.resolve())],
                cwd=cwd_tmp,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn(manifest["execution_paths"]["planner"], proc.stdout)
        self.assertIn(str(manifest_path.resolve()), proc.stdout)
        self.assertNotIn(str(script_dir), proc.stdout)
        self.assertNotIn("dirname", (PACKAGE / "src" / "01_cnv_array_unit.sh").read_text(encoding="utf-8"))

    def test_source_has_no_submission_or_data_acquisition_client(self) -> None:
        worker = (PACKAGE / "src" / "cnv_array_plan.py").read_text(encoding="utf-8")
        unit_shell = (PACKAGE / "src" / "01_cnv_array_unit.sh").read_text(encoding="utf-8")
        self.assertNotIn("subprocess.run([\"sbatch\"", worker)
        self.assertNotIn("aws s3", (worker + unit_shell).lower())
        self.assertNotIn("module load modtree/deprecated", (worker + unit_shell).lower())
        self.assertNotIn("minimap2 -d", (worker + unit_shell).lower())
        self.assertIn("APPTAINER_BIND=/cluster", unit_shell)


class LiveSeamTests(unittest.TestCase):
    def test_fake_samtools_quickcheck_and_index_query(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            bam = root / "assembly.bam"
            bai = root / "assembly.bam.bai"
            bam.write_bytes(b"assembly-bam-fixture\n")
            bai.write_bytes(b"assembly-bai-fixture\n")
            fake = root / "samtools"
            fake.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "if [[ $1 == quickcheck ]]; then exit 0; fi\n"
                "if [[ $1 == view && $2 == -c ]]; then printf '7\\n'; exit 0; fi\n"
                "if [[ $1 == --version ]]; then printf 'samtools 1.21\\n'; exit 0; fi\n"
                "exit 9\n",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
            manifest["runtime"]["tools"]["samtools"]["path"] = str(fake)
            hap = copy.deepcopy(manifest["units"][0]["haplotypes"][0])
            hap["assembly_alignment_bam"] = {
                "path": str(bam), "size_bytes": bam.stat().st_size,
                "sha256": hashlib.sha256(bam.read_bytes()).hexdigest(),
                "semantic_type": "assembly_contig_to_t2t_alignment_not_read_depth",
            }
            hap["assembly_alignment_bai"] = {
                "path": str(bai), "size_bytes": bai.stat().st_size,
                "sha256": hashlib.sha256(bai.read_bytes()).hexdigest(),
            }
            hap["expected_indexed_window_record_count"] = 7
            observed = plan._validate_haplotype_live(manifest, hap)
            self.assertTrue(observed["quickcheck_passed"])
            self.assertEqual(observed["indexed_window_record_count"], 7)
            self.assertEqual(observed["bam"]["sha256"], hap["assembly_alignment_bam"]["sha256"])

    def test_checkpoint_validation_and_bad_output_detection(self) -> None:
        manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            evidence = root / "assembly_structural_evidence.tsv"
            evidence.write_bytes(b"fixture\n")
            checkpoint_path = root / "CHECKPOINT.json"
            manifest["units"][0]["outputs"] = {"evidence_tsv": str(evidence), "checkpoint_json": str(checkpoint_path)}
            checkpoint = {
                "schema_version": plan.UNIT_SCHEMA,
                "run_id": manifest["run_id"],
                "plan_sha256": manifest["plan_sha256"],
                "unit_index": 0,
                "sample_id": "HG02027",
                "output": {"path": str(evidence), "size_bytes": evidence.stat().st_size, "sha256": hashlib.sha256(evidence.read_bytes()).hexdigest()},
                "read_depth_claim": "none",
                "copy_number_inferred_from_alignment_count": False,
                "checkpoint_sha256": "",
            }
            checkpoint["checkpoint_sha256"] = plan.digest_without(checkpoint, "checkpoint_sha256")
            checkpoint_path.write_bytes(plan.canonical_bytes(checkpoint))
            self.assertEqual(plan.validate_unit(manifest, 0)["sample_id"], "HG02027")
            evidence.write_bytes(b"changed\n")
            with self.assertRaises(plan.PlanError):
                plan.validate_unit(manifest, 0)


if __name__ == "__main__":
    unittest.main()
