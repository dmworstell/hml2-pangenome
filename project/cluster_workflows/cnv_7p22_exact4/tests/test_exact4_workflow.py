#!/usr/bin/env python3
"""Isolated, offline tests for the exact-four 7p22.1 candidate bundle."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
import exact4_plan as plan
import plot_exact4_depth as plot


def canonical(value: object) -> bytes:
    return plan.canonical_json_bytes(value)


def local_mechanics_manifest() -> dict[str, object]:
    """Build local bundle identity while truthfully preserving live pin drift.

    This helper is test-only.  The real builder continues to fail closed until
    the amended contract/acceptance hashes are frozen by the coordinator.
    """
    current_pins = tuple((relative, plan.sha256_file(PROJECT / relative)) for relative, _ in plan.PREREQUISITE_PINS)
    with mock.patch.object(plan, "PREREQUISITE_PINS", current_pins):
        return plan.build_deployment_manifest(ROOT)


def fake_tools(tool_root: Path) -> list[dict[str, object]]:
    outputs = {
        "bash": "GNU bash, version 5.2\n",
        "python": "Python 3.11.9\n",
        "aws": "aws-cli/1.16.308 Python/3.9.9\n",
        "minimap2": "2.26-r1175\n",
        "samtools": "samtools 1.21\n",
        "sbatch": "slurm 24.11\n",
    }
    rows = []
    tool_root.mkdir()
    for name in ("bash", "python", "aws", "minimap2", "samtools", "sbatch"):
        version = outputs[name]
        help_output = "-Q, --min-MQ minimum mapping quality\n-g FLAGS include flag bits\n" if name == "samtools" else f"{name} help\n"
        executable = tool_root / name
        executable.write_bytes(name.encode("ascii")); executable.chmod(0o755)
        path = str(executable)
        digest = hashlib.sha256(name.encode()).hexdigest()
        if name == "minimap2":
            digest = plan.MINIMAP2_RESIDENT_SHA256
        elif name == "samtools":
            digest = plan.SAMTOOLS_RESIDENT_SHA256
        rows.append({
            "name": name, "path": path, "size_bytes": len(name),
            "sha256": digest,
            "version": {"argv": [path, "--version"], "output": version, "output_sha256": hashlib.sha256(version.encode()).hexdigest(), "returncode": 0},
            "help": {"argv": [path, "depth", "--help"] if name == "samtools" else [path, "--help"], "output": help_output, "output_sha256": hashlib.sha256(help_output.encode()).hexdigest(), "returncode": 0},
        })
    return rows


def _artifact(path: Path, payload: bytes) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {"path": str(path), "size_bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}


def synthetic_fixture_environment(
    manifest: dict[str, object], tools: list[dict[str, object]], root: Path
) -> tuple[dict[str, object], dict[str, Path]]:
    """Create a fully replayable, externally authorized synthetic receipt."""
    by_name = {row["name"]: row for row in tools}
    minimap, samtools = by_name["minimap2"]["path"], by_name["samtools"]["path"]
    patterns = {
        plan.FIXTURE_IDS[0]: [[samtools, "faidx", "toy.fa"], [minimap, "-d", "toy.mmi", "toy.fa"]],
        plan.FIXTURE_IDS[1]: [[samtools, "view", "-b", "-o", "toy.bam", "toy.sam"], [samtools, "view", "-F", "0x900", "toy.bam"], [samtools, "fasta", "-F", "0x900", "toy.bam"]],
        plan.FIXTURE_IDS[2]: [[minimap, "-ax", "map-ont", "--secondary=no", "toy.fa", "reads.fa"]],
        plan.FIXTURE_IDS[3]: [[samtools, "view", "-u", "-F", "0x900", "-L", "toy.bed", "toy.bam"], [samtools, "sort", "-@", "4", "toy.bam"], [samtools, "index", "-@", "4", "toy.bam"], [samtools, "quickcheck", "toy.bam"]],
        plan.FIXTURE_IDS[4]: [[samtools, "merge", "-@", "4", "merged.bam", "a.bam", "b.bam"], [samtools, "index", "-@", "4", "merged.bam"], [samtools, "quickcheck", "merged.bam"]],
        plan.FIXTURE_IDS[5]: [[samtools, "depth", "-a", "-g", "0x600", "-b", "toy.bed", "toy.bam"], [samtools, "depth", "-a", "-g", "0x600", "-Q", "10", "-b", "toy.bed", "toy.bam"]],
        plan.FIXTURE_IDS[6]: [[minimap, "-c", "-p", "0.1", "-N", "10", "target.fa", "query.fa"]],
    }
    manifest_path = root / "exact4_deployment_manifest.v1.json"
    manifest_bytes = canonical(manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_bytes(manifest_bytes)
    auth_path = root / "fixture-authorization.json"
    acceptance_path = root / "fixture-authorization-acceptance.json"
    auth = {
        "schema_version": plan.TOOL_FIXTURE_AUTH_SCHEMA,
        "deployment_bundle_id": manifest["bundle_id"],
        "deployment_manifest": {"path": str(manifest_path), "size_bytes": len(manifest_bytes), "sha256": hashlib.sha256(manifest_bytes).hexdigest()},
        "correction_addendum": {"path": plan.V2_CORRECTION_ADDENDUM_PATH, "size_bytes": plan.V2_CORRECTION_ADDENDUM_SIZE, "sha256": plan.V2_CORRECTION_ADDENDUM_SHA256},
        "correction_acceptance": {"path": plan.V2_CORRECTION_ACCEPTANCE_PATH, "size_bytes": plan.V2_CORRECTION_ACCEPTANCE_SIZE, "sha256": plan.V2_CORRECTION_ACCEPTANCE_SHA256},
        "host_class": "login.pax.tufts.edu",
        "tool_bindings": {name: {key: by_name[name][key] for key in ("path", "size_bytes", "sha256", "version", "help")} for name in ("minimap2", "samtools")},
        "accepted_roots": {
            "resident_rehash_receipt": {"path": plan.RESIDENT_REHASH_RECEIPT_PATH, "size_bytes": 4571, "sha256": plan.RESIDENT_REHASH_RECEIPT_SHA256},
            "resident_rehash_acceptance": {"path": plan.RESIDENT_REHASH_ACCEPTANCE_PATH, "size_bytes": 1929, "sha256": plan.RESIDENT_REHASH_ACCEPTANCE_SHA256},
            "aws_tool_authority": {"path": plan.AWS_TOOL_AUTHORITY_PATH, "size_bytes": 1165, "sha256": plan.AWS_TOOL_AUTHORITY_SHA256},
            "aws_tool_acceptance": {"path": plan.AWS_TOOL_ACCEPTANCE_PATH, "size_bytes": 1855, "sha256": plan.AWS_TOOL_ACCEPTANCE_SHA256},
        },
        "fixture_order": list(plan.FIXTURE_IDS),
        "scratch_root": str(root / "authorized-one-shot"),
        "read_only_action_scope": plan.FIXTURE_READ_ONLY_SCOPE,
        "artifact_policy": {"create_only_local_fixture_artifacts": True, "reuse_authorized": False, "retry_loop_authorized": False},
        "action_gates": dict(plan.ACTION_GATES),
    }
    auth_bytes = canonical(auth)
    auth_path.write_bytes(auth_bytes)
    acceptance = {
        "schema_version": plan.TOOL_FIXTURE_AUTH_ACCEPTANCE_SCHEMA,
        "decision": "ACCEPT",
        "authorization": {"path": str(auth_path), "size_bytes": len(auth_bytes), "sha256": hashlib.sha256(auth_bytes).hexdigest()},
    }
    acceptance_bytes = canonical(acceptance)
    acceptance_path.write_bytes(acceptance_bytes)
    authorization_binding = {
        "authorization_path": str(auth_path), "authorization_size_bytes": len(auth_bytes),
        "authorization_sha256": hashlib.sha256(auth_bytes).hexdigest(),
        "acceptance_path": str(acceptance_path), "acceptance_size_bytes": len(acceptance_bytes),
        "acceptance_sha256": hashlib.sha256(acceptance_bytes).hexdigest(),
    }

    fixtures = []
    for fixture_id in plan.FIXTURE_IDS:
        directory = root / fixture_id
        directory.mkdir(parents=True, exist_ok=True)
        inputs: list[dict[str, object]] = []
        outputs: list[dict[str, object]] = []
        if fixture_id == plan.FIXTURE_IDS[0]:
            inputs.append(_artifact(directory / "combined.fa", b">hapA\nACGTACGTACGT\n>hapB\nTTGCTTGCTTGC\n"))
            outputs.extend((_artifact(directory / "combined.fa.fai", b"hapA\t12\t6\t12\t13\nhapB\t12\t25\t12\t13\n"), _artifact(directory / "bait.fa", b">hapA:2-9\nCGTACGTA\n"), _artifact(directory / "combined.mmi", b"toy-index\n")))
        elif fixture_id == plan.FIXTURE_IDS[1]:
            rows = ["@HD\tVN:1.6", "@SQ\tSN:hapA\tLN:20"] + [f"{name}\t{flag}\thapA\t1\t20\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF" for name, flag in (("primary", 0), ("secondary", 256), ("supplementary", 2048), ("duplicate", 1024), ("qcfail", 512))]
            inputs.append(_artifact(directory / "flags.sam", ("\n".join(rows) + "\n").encode("ascii")))
            selected = [row for row in rows[2:] if int(row.split("\t")[1]) & 0x900 == 0]
            outputs.append(_artifact(directory / "primary.sam", ("\n".join(selected) + "\n").encode("ascii")))
            outputs.append(_artifact(directory / "primary.fa", b">primary\nACGTACGT\n>duplicate\nACGTACGT\n>qcfail\nACGTACGT\n"))
        elif fixture_id == plan.FIXTURE_IDS[2]:
            inputs.extend((_artifact(directory / "combined.fa", b">hapA\nACGTACGTACGT\n>hapB\nTTGCTTGCTTGC\n"), _artifact(directory / "reads.fa", b">readA\nACGTACGT\n>readB\nTTGCTTGC\n")))
            outputs.append(_artifact(directory / "aligned.sam", b"@HD\tVN:1.6\n@SQ\tSN:hapA\tLN:12\n@SQ\tSN:hapB\tLN:12\nreadA\t0\thapA\t1\t60\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF\nreadB\t0\thapB\t1\t60\t8M\t*\t0\t0\tTTGCTTGC\tFFFFFFFF\n"))
        elif fixture_id == plan.FIXTURE_IDS[3]:
            rows = ["@HD\tVN:1.6", "@SQ\tSN:hapA\tLN:20"] + [f"{name}\t{flag}\thapA\t1\t20\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF" for name, flag in (("primary", 0), ("secondary", 256), ("supplementary", 2048), ("duplicate", 1024), ("qcfail", 512))]
            inputs.append(_artifact(directory / "window.sam", ("\n".join(rows) + "\n").encode("ascii")))
            selected = [row for row in rows[2:] if int(row.split("\t")[1]) & 0x900 == 0]
            outputs.append(_artifact(directory / "selected.sam", ("\n".join(selected) + "\n").encode("ascii")))
        elif fixture_id == plan.FIXTURE_IDS[4]:
            inputs.extend((_artifact(directory / "partA.qnames", b"a1\n"), _artifact(directory / "partB.qnames", b"b1\n"), _artifact(directory / "duplicate.qnames", b"a1\n")))
            outputs.append(_artifact(directory / "merged.sam", b"@HD\tVN:1.6\n@SQ\tSN:hapA\tLN:20\na1\t0\thapA\t1\t20\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF\nb1\t0\thapA\t1\t20\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF\n"))
        elif fixture_id == plan.FIXTURE_IDS[5]:
            raw = b"".join(f"hapA\t{i}\t{4 if i <= 10 else 0}\n".encode("ascii") for i in range(1, 21))
            q10 = b"".join(f"hapA\t{i}\t{3 if i <= 10 else 0}\n".encode("ascii") for i in range(1, 21))
            outputs.extend((_artifact(directory / "raw.depth", raw), _artifact(directory / "mapq10.depth", q10)))
        else:
            authority = {"sample_id": "FIXTURE", "haplotypes": ["hapA", "hapB"], "bed_authorities": [{"haplotype": "hapA", "contig": "toy", "full_intervals": [[0, 12000]]}], "element_identities": [{"haplotype": "hapA", "contig": "toy", "source_start": 1500, "source_end": 2500}]}
            inputs.extend((_artifact(directory / "full-window.fa", b">toy:1-12000\n" + b"A" * 12000 + b"\n"), _artifact(directory / "type2_KCON.fa", b">type2_KCON\n" + b"A" * 9472 + b"\n")))
            inputs.append(_artifact(directory / "selector_authority.json", canonical(authority)))
            boundary = {}
            for name, payload in (("boundary.bam", b"synthetic-bam-boundary\n"), ("boundary.depth", b"synthetic-depth-boundary\n"), ("boundary.callability", b"synthetic-callability-boundary\n")):
                row = _artifact(directory / name, payload); inputs.append(row); boundary[name] = {"size_bytes": row["size_bytes"], "sha256": row["sha256"]}
            inputs.append(_artifact(directory / "boundary_before.json", canonical(boundary)))
            outputs.append(_artifact(directory / "resident_minimap.paf", b"type2_KCON\t9472\t0\t9472\t+\ttoy:1-12000\t12000\t500\t9972\t9472\t9472\t60\n"))
            paf = b"".join((
                b"type2_KCON\t9472\t0\t499\t+\ttoy:1-12000\t12000\t1500\t1999\t499\t499\t60\n",
                b"type2_KCON\t9472\t0\t500\t+\ttoy:1-12000\t12000\t1500\t2000\t500\t500\t60\n",
                b"type2_KCON\t9472\t50\t750\t+\ttoy:1-12000\t12000\t1700\t2400\t680\t700\t50\n",
                b"type2_KCON\t9472\t40\t740\t+\ttoy:1-12000\t12000\t1700\t2400\t690\t700\t40\n",
                b"type2_KCON\t9472\t30\t730\t+\ttoy:1-12000\t12000\t1700\t2400\t690\t700\t50\n",
                b"type2_KCON\t9472\t20\t720\t+\ttoy:1-12000\t12000\t1600\t2300\t690\t700\t50\n",
                b"type2_KCON\t9472\t10\t710\t+\ttoy:1-12000\t12000\t1600\t2300\t690\t700\t50\tzz:Z:b\n",
                b"type2_KCON\t9472\t10\t710\t+\ttoy:1-12000\t12000\t1600\t2300\t690\t700\t50\tzz:Z:a\n",
                b"type2_KCON\t9472\t0\t700\t+\ttoy:1-12000\t12000\t2600\t3300\t700\t700\t60\n",
                b"wrong_query\t9472\t0\t700\t+\ttoy:1-12000\t12000\t1600\t2300\t700\t700\t60\n",
                b"type2_KCON\t9471\t0\t700\t+\ttoy:1-12000\t12000\t1600\t2300\t700\t700\t60\n",
                b"type2_KCON\t9472\t0\t700\t+\twrong:1-12000\t12000\t1600\t2300\t700\t700\t60\n",
            ))
            outputs.extend((_artifact(directory / "kcon.paf", paf), _artifact(directory / "boundary_after.json", canonical(boundary))))
        fixture = {
            "fixture_id": fixture_id,
            "scratch_directory": str(directory),
            "input_files": inputs,
            "commands": [{"argv": argv, "returncode": 0, "stdout_sha256": "b" * 64, "stderr_sha256": "c" * 64} for argv in patterns[fixture_id]],
            "output_files": outputs,
            "record_counts": {},
            "semantic_assertions": {},
            "passed": True,
        }
        counts, assertions = plan._recompute_fixture_oracle(fixture, manifest)
        fixture["record_counts"] = counts
        fixture["semantic_assertions"] = assertions
        fixtures.append(fixture)
    receipt = {
        "schema_version": plan.TOOL_FIXTURE_SCHEMA,
        "host_class": "login.pax.tufts.edu",
        "deployment_manifest": {"path": str(manifest_path), "size_bytes": len(manifest_bytes), "sha256": hashlib.sha256(manifest_bytes).hexdigest()},
        "authorization_binding": authorization_binding,
        "tool_bindings": {name: {key: by_name[name][key] for key in ("path", "size_bytes", "sha256", "version", "help")} for name in ("minimap2", "samtools")},
        "fixture_order": list(plan.FIXTURE_IDS),
        "fixtures": fixtures,
        "action_gates": dict(plan.ACTION_GATES),
    }
    path_map = {
        plan.DEPLOYMENT_MANIFEST_PATH: manifest_path,
        plan.FIXTURE_AUTHORIZATION_PATH: auth_path,
        plan.FIXTURE_AUTHORIZATION_ACCEPTANCE_PATH: acceptance_path,
    }
    return receipt, path_map


def synthetic_resident(
    manifest: dict[str, object], tool_root: Path,
    *, tools: list[dict[str, object]] | None = None,
    tool_fixture_receipt: dict[str, object] | None = None,
) -> dict[str, object]:
    assemblies = []
    elements = []
    beds = []
    for assembly in manifest["assembly_authorities"]:
        sample, hap = assembly["sample_id"], assembly["haplotype"]
        contig = f"{sample}_{hap}_ctg"
        assemblies.append({
            "sample_id": sample, "haplotype": hap,
            "fasta": {"path": assembly["fasta_path"], "size_bytes": 100, "sha256": assembly["fasta_sha256"]},
            "fai": {"path": assembly["fai_path"], "size_bytes": 20, "sha256": assembly["fai_sha256"]},
            "contigs": {contig: 1000},
        })
        elements.append({"sample_id": sample, "haplotype": hap, "contig": contig, "source_start": 200, "source_end": 300, "core_start": 100, "core_end": 400, "full_start": 0, "full_end": 600})
        beds.append({"sample_id": sample, "haplotype": hap, "contig": contig, "core_intervals": [[100, 400]], "full_intervals": [[0, 600]], "flank_intervals": [[0, 100], [400, 600]]})
    combined_row_sets = [{"sample_id": sample, "locus": plan.LOCUS, "header": ["Locus", "ID", "Haplotype", "Source_Identifier"], "rows": []} for sample in plan.TARGETS]
    counts = {
        "by_assembly": [{"sample_id": row["sample_id"], "haplotype": row["haplotype"], "count": 1} for row in assemblies],
        "by_sample": [{"sample_id": sample, "count": 2} for sample in plan.TARGETS],
    }
    sample_beds = []
    for sample in plan.TARGETS:
        core = plan._bed_bytes(beds, sample, "core_intervals", "_core")
        full = plan._bed_bytes(beds, sample, "full_intervals", "_full")
        sample_beds.append({"sample_id": sample, "core_bed": {"size_bytes": len(core), "sha256": hashlib.sha256(core).hexdigest(), "text": core.decode("ascii")}, "full_window_bed": {"size_bytes": len(full), "sha256": hashlib.sha256(full).hexdigest(), "text": full.decode("ascii")}})
    tools = tools or fake_tools(tool_root)
    if tool_fixture_receipt is None:
        tool_fixture_receipt, _ = synthetic_fixture_environment(manifest, tools, tool_root / "fixtures")
    return {
        "schema_version": plan.RESIDENT_SCHEMA,
        "deployment_bundle_id": manifest["bundle_id"],
        "deployment_manifest_sha256": plan.digest_value(manifest),
        "resident_files": [{"role": row["role"], "path": row["path"], "size_bytes": row["size_bytes"], "sha256": row["sha256"]} for row in manifest["resident_authorities"]],
        "assemblies": assemblies,
        "element_identities": elements,
        "bed_authorities": beds,
        "combined_row_sets": combined_row_sets,
        "assembly_encoded_counts": counts,
        "sample_bed_authorities": sample_beds,
        "tools": tools,
        "tool_fixture_receipt": tool_fixture_receipt,
        "storage": {"path": "/fixture", "block_size": 4096, "blocks_available": 1000, "inodes_available": 1000},
        "resources": manifest["resources"],
        "capture_metadata": {"observation_label": "fixture", "network_contact_performed": False, "scheduler_query_performed": False, "filesystem_write_performed": False},
    }


def synthetic_provider(manifest: dict[str, object]) -> dict[str, object]:
    objects = []
    for row in manifest["expected_provider_objects"]:
        objects.append({
            "sample_id": row["sample_id"], "part_index": row["part_index"], "uri": row["uri"],
            "version_id": row["version_id"], "etag": row["etag"],
            "content_length_bytes": row["content_length_bytes"], "last_modified": row["last_modified"],
            "observation_time": "fixture-time", "provider_command": ["fixture", "head"],
            "response_digest": hashlib.sha256((row["uri"] + row["version_id"]).encode()).hexdigest(),
        })
    aws_authority = json.loads((ROOT / "frozen" / "exact4_aws_tool_authority.v1.json").read_text())
    amended = {
        "authorization_path": plan.AMENDED_AUTHORIZATION_PATH,
        "authorization_size_bytes": 1,
        "authorization_sha256": "2" * 64,
        "acceptance_path": plan.AMENDED_AUTHORIZATION_ACCEPTANCE_PATH,
        "acceptance_size_bytes": 1,
        "acceptance_sha256": "3" * 64,
        "deployment_bundle_id": manifest["bundle_id"],
        "deployment_manifest_path": plan.DEPLOYMENT_MANIFEST_PATH,
        "deployment_manifest_size_bytes": len(canonical(manifest)),
        "deployment_manifest_sha256": plan.digest_value(manifest),
        "provider_receipt_schema": plan.PROVIDER_SCHEMA,
        "biological_target": plan.LOCUS,
        "ordered_samples": list(plan.TARGETS),
        "provider_part_counts": {sample: count for sample, count in zip(plan.TARGETS, plan.PART_COUNTS)},
        "read_only_action_scope": plan.AMENDED_READ_ONLY_SCOPE,
    }
    value = {
        "schema_version": plan.PROVIDER_SCHEMA, "deployment_bundle_id": manifest["bundle_id"],
        "objects": objects,
        "provider_tool_receipt": {"path": aws_authority["normalized_absolute_path"], "size_bytes": 1, "sha256": aws_authority["sha256"], "version_output": aws_authority["version_output"], "version_output_sha256": hashlib.sha256((aws_authority["version_output"] + "\n").encode()).hexdigest()},
        "capture_metadata": {"observation_host": "login-3.cluster.example", "capture_label": plan.PROVIDER_CAPTURE_LABEL, "object_body_downloaded": False},
        "authority_bindings": {
            "expected_provider_metadata_authority": {"path": "/fixture/expected.json", "size_bytes": 5801, "sha256": plan.EXPECTED_PROVIDER_AUTHORITY_SHA256, "acceptance_path": "/fixture/expected-acceptance.md", "acceptance_size_bytes": 2089, "acceptance_sha256": plan.EXPECTED_PROVIDER_ACCEPTANCE_SHA256, "deployment_bundle_id": "exact4-bundle-cee202a9f2a520062573fe50", "capture_authorization_sha256": plan.LEGACY_CAPTURE_AUTHORIZATION_SHA256},
            "aws_tool_authority": {"path": "/fixture/aws.json", "size_bytes": 1165, "sha256": plan.AWS_TOOL_AUTHORITY_SHA256, "acceptance_path": "/fixture/aws-acceptance.md", "acceptance_size_bytes": 1855, "acceptance_sha256": plan.AWS_TOOL_ACCEPTANCE_SHA256, "normalized_absolute_path": aws_authority["normalized_absolute_path"], "executable_sha256": aws_authority["sha256"], "version_output": aws_authority["version_output"], "module_name": aws_authority["module_name"], "cluster_host_observation_class": aws_authority["cluster_host_observation_class"], "read_only_resident_receipt_sha256": plan.RESIDENT_REHASH_RECEIPT_SHA256, "capture_authorization_sha256": plan.LEGACY_CAPTURE_AUTHORIZATION_SHA256},
            "amended_live_capture_authorization": amended,
        },
        "action_gates": dict(plan.ACTION_GATES),
    }
    aws_path = value["provider_tool_receipt"]["path"]
    for row in value["objects"]:
        uri = row["uri"][5:]
        bucket, key = uri.split("/", 1)
        row["provider_command"] = [aws_path, "s3api", "head-object", "--bucket", bucket, "--key", key, "--version-id", row["version_id"], "--no-sign-request"]
    return plan._validate_provider_receipt(value, manifest, revalidate_external_files=False)


class ExactFourFixture(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.tmp = Path(self.temporary.name).resolve()
        self.patch_root = mock.patch.object(plan, "RUN_ROOT", str(self.tmp / "runs"))
        self.patch_root.start()
        self.manifest = local_mechanics_manifest()
        self.tools = fake_tools(self.tmp / "tools")
        fixture_receipt, self.fixture_path_map = synthetic_fixture_environment(
            self.manifest, self.tools, self.tmp / "fixture-environment"
        )
        self.patch_project_paths = mock.patch.object(
            plan,
            "_project_path_from_relative",
            side_effect=lambda relative: self.fixture_path_map.get(relative, PROJECT / relative),
        )
        self.patch_project_paths.start()
        self.patch_fixture_host = mock.patch.object(
            plan.os,
            "uname",
            return_value=type("Uname", (), {"nodename": "login-03.cluster.example"})(),
        )
        self.patch_fixture_host.start()
        self.resident = synthetic_resident(
            self.manifest, self.tmp / "tools", tools=self.tools,
            tool_fixture_receipt=fixture_receipt,
        )
        derived = tuple(copy.deepcopy(self.resident[key]) for key in ("assemblies", "element_identities", "bed_authorities", "combined_row_sets", "assembly_encoded_counts", "sample_bed_authorities"))
        self.patch_derive = mock.patch.object(plan, "_derive_resident_biology", return_value=derived)
        self.patch_derive.start()
        self.patch_current_source = mock.patch.object(plan, "_validate_current_source_bytes", return_value=None)
        self.current_source_mock = self.patch_current_source.start()
        self.patch_stable_authorities = mock.patch.object(
            plan,
            "_stable_reopen_pinned_file",
            side_effect=lambda relative, *_: canonical(self.manifest) if relative == plan.DEPLOYMENT_MANIFEST_PATH else b"x",
        )
        self.patch_stable_authorities.start()
        real_fingerprint = plan._fingerprint
        minimap_path = Path(self.resident["tools"][3]["path"])
        samtools_path = Path(self.resident["tools"][4]["path"])
        def fixture_fingerprint(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
            observed = real_fingerprint(path, allow_empty=allow_empty)
            if path == minimap_path and path.read_bytes() == b"minimap2":
                observed["sha256"] = plan.MINIMAP2_RESIDENT_SHA256
            elif path == samtools_path and path.read_bytes() == b"samtools":
                observed["sha256"] = plan.SAMTOOLS_RESIDENT_SHA256
            return observed
        self.patch_fixture_fingerprint = mock.patch.object(plan, "_fingerprint", side_effect=fixture_fingerprint)
        self.patch_fixture_fingerprint.start()
        self.resident = plan._validate_resident_receipt(self.resident, self.manifest)
        self.provider = synthetic_provider(self.manifest)
        self.plan = plan.build_target_only_plan(self.manifest, self.resident, self.provider)
        self.plan_path = self.tmp / "plan.json"
        self.plan_path.write_bytes(canonical(self.plan))
        self.release = self.make_release()
        self.release_path = self.tmp / "release.json"
        self.release_path.write_bytes(canonical(self.release))

    def tearDown(self) -> None:
        self.patch_fixture_fingerprint.stop()
        self.patch_stable_authorities.stop()
        self.patch_current_source.stop()
        self.patch_derive.stop()
        self.patch_fixture_host.stop()
        self.patch_project_paths.stop()
        self.patch_root.stop()
        self.temporary.cleanup()

    def make_release(self) -> dict[str, object]:
        return {
            "schema_version": plan.RELEASE_SCHEMA,
            "plan_file_sha256": plan.sha256_file(self.plan_path), "plan_digest": self.plan["plan_digest"],
            "run_id": self.plan["run_id"], "bundle_id": self.plan["bundle_id"],
            "deployment_manifest_sha256": self.plan["deployment_manifest_sha256"],
            "resident_receipt_sha256": self.plan["resident_receipt_sha256"],
            "provider_receipt_sha256": self.plan["provider_receipt_sha256"],
            **plan.CONTRACT_BINDINGS,
            "coverage_ratio_schema_identity": plan.COVERAGE_RATIO_SCHEMA_IDENTITY,
            "coverage_ratio_schema_sha256": plan.COVERAGE_RATIO_SCHEMA_SHA256,
            "no_drift_receipt_sha256": "2" * 64, "no_drift_identity_sha256": "3" * 64,
            "entrypoint_hashes": self.plan["entrypoint_hashes"], "tool_receipts": self.resident["tools"],
            "phase_maps": plan._phase_maps(self.plan), "resources": self.plan["resources"],
            "allowed_phases": ["preparation", "alignment", "analysis"], "one_attempt": True,
            "attempt_number": 1, "allowed_claims": plan._allowed_claims(self.plan),
            "allowed_destinations": plan._allowed_destinations(self.plan), "execution_authorized": True,
            "production_authorized": False, "copy_number_inference_authorized": False,
            "model_use_authorized": False, "manuscript_use_authorized": False,
        }


class TestStrictJson(unittest.TestCase):
    def test_duplicate_and_nonfinite_rejected(self) -> None:
        with self.assertRaisesRegex(plan.ExactFourPlanError, "duplicate"):
            plan.strict_json_bytes(b'{"a":1,"a":2}\n', "fixture")
        for token in (b"NaN", b"Infinity", b"-Infinity"):
            with self.assertRaises(plan.ExactFourPlanError):
                plan.strict_json_bytes(token, "fixture")

    def test_wrong_scalar_and_float_rejected(self) -> None:
        with self.assertRaisesRegex(plan.ExactFourPlanError, "floats"):
            plan.strict_json_bytes(b'{"x":1.5}\n', "fixture")
        with self.assertRaises(plan.ExactFourPlanError):
            plan._require_int(True, "integer")

    def test_canonical_utf8_sorted_lf(self) -> None:
        self.assertEqual(plan.canonical_json_bytes({"z": "é", "a": 1}), b'{"a":1,"z":"\xc3\xa9"}\n')
        self.assertEqual(plan.strict_json_bytes(b'{"a":1}\n', "fixture"), {"a": 1})


class TestCurrentResidentIdentity(unittest.TestCase):
    def test_only_phase_relevant_source_bytes_are_freshly_rehashed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            resident_files = []
            roles = ("assembly_authority", "combined_table", "type2_kcon", "ont_index")
            for index, role in enumerate(roles):
                path = root / f"resident-{index}.dat"
                path.write_bytes(f"resident-{index}".encode("ascii"))
                resident_files.append({"role": role, **plan._fingerprint(path)})

            assemblies = []
            for index, (sample, haplotype, *_) in enumerate(plan.ASSEMBLY_AUTHORITIES):
                fasta = root / f"assembly-{index}.fa"
                fai = root / f"assembly-{index}.fa.fai"
                fasta.write_bytes(f">ctg-{index}\nACGT\n".encode("ascii"))
                fai.write_bytes(f"ctg-{index}\t4\t7\t4\t5\n".encode("ascii"))
                assemblies.append({
                    "sample_id": sample,
                    "haplotype": haplotype,
                    "fasta": plan._fingerprint(fasta),
                    "fai": plan._fingerprint(fai),
                    "contigs": {f"ctg-{index}": 4},
                })
            identity = {"resident_files": resident_files, "assemblies": assemblies}
            sample = plan.TARGETS[0]
            plan._validate_current_source_bytes(identity, "prepare", sample)
            plan._validate_current_source_bytes(identity, "analysis", sample)

            for record, kind, expected_message in (
                (assemblies[0]["fasta"], "prepare", "current assembly fasta bytes drifted"),
                (assemblies[1]["fai"], "prepare", "current assembly fai bytes drifted"),
                (resident_files[2], "analysis", "current resident file bytes drifted"),
            ):
                path = Path(record["path"])
                original = path.read_bytes()
                try:
                    path.write_bytes(original + b"drift")
                    with self.assertRaisesRegex(plan.ExactFourPlanError, expected_message):
                        plan._validate_current_source_bytes(identity, kind, sample)
                finally:
                    path.write_bytes(original)

            unrelated = Path(assemblies[-1]["fasta"]["path"])
            unrelated.write_bytes(unrelated.read_bytes() + b"drift")
            plan._validate_current_source_bytes(identity, "prepare", sample)


class TestStableAuthorityWalkerConcurrency(unittest.TestCase):
    """Real filesystem races at the stable-reader seam, without subprocesses."""

    def _tree(self) -> tuple[tempfile.TemporaryDirectory[str], Path, Path, bytes]:
        temporary = tempfile.TemporaryDirectory()
        ancestor = Path(temporary.name).resolve() / "authority-parent"
        ancestor.mkdir()
        leaf = ancestor / "accepted.bin"
        payload = b"accepted immutable authority bytes\n"
        leaf.write_bytes(payload)
        return temporary, ancestor, leaf, payload

    def _readers(self, leaf: Path, payload: bytes) -> dict[str, object]:
        digest = hashlib.sha256(payload).hexdigest()

        def absolute() -> bytes:
            return plan._stable_reopen_absolute_pinned_file(
                str(leaf), len(payload), digest, "absolute authority"
            )

        def relative() -> bytes:
            with mock.patch.object(plan, "_project_path_from_relative", return_value=leaf):
                return plan._stable_reopen_pinned_file(
                    "fixture/accepted.bin", len(payload), digest, "relative authority"
                )

        return {"absolute": absolute, "relative": relative}

    def _during_first_read(self, mutation: object, reader: object) -> bytes:
        real_read = os.read
        fired = False

        def read_once(descriptor: int, count: int) -> bytes:
            nonlocal fired
            if not fired:
                fired = True
                mutation()
            return real_read(descriptor, count)

        with mock.patch.object(plan.os, "read", side_effect=read_once):
            result = reader()
        self.assertTrue(fired)
        return result

    def test_absolute_and_relative_accept_unrelated_sibling_create_and_remove(self) -> None:
        for operation in ("create", "remove"):
            for reader_name in ("absolute", "relative"):
                with self.subTest(operation=operation, reader=reader_name):
                    temporary, ancestor, leaf, payload = self._tree()
                    try:
                        sibling = ancestor / "unrelated.tmp"
                        if operation == "remove":
                            sibling.write_bytes(b"unrelated\n")
                        before_mtime = ancestor.stat().st_mtime_ns

                        def churn() -> None:
                            if operation == "create":
                                sibling.write_bytes(b"unrelated\n")
                            else:
                                sibling.unlink()
                            current = ancestor.stat()
                            os.utime(
                                ancestor,
                                ns=(current.st_atime_ns, max(current.st_mtime_ns, before_mtime) + 1_000_000_000),
                            )

                        reader = self._readers(leaf, payload)[reader_name]
                        self.assertEqual(self._during_first_read(churn, reader), payload)
                        self.assertNotEqual(ancestor.stat().st_mtime_ns, before_mtime)
                        self.assertEqual(leaf.read_bytes(), payload)
                    finally:
                        temporary.cleanup()

    def test_absolute_and_relative_refuse_real_ancestor_and_leaf_changes(self) -> None:
        mutations = ("ancestor_directory", "ancestor_symlink", "leaf_mutation", "leaf_replacement")
        for mutation_kind in mutations:
            for reader_name in ("absolute", "relative"):
                with self.subTest(mutation=mutation_kind, reader=reader_name):
                    temporary, ancestor, leaf, payload = self._tree()
                    try:
                        parked = ancestor.with_name("parked-authority-parent")

                        def mutate() -> None:
                            if mutation_kind.startswith("ancestor_"):
                                ancestor.rename(parked)
                                if mutation_kind == "ancestor_symlink":
                                    ancestor.symlink_to(parked, target_is_directory=True)
                                else:
                                    ancestor.mkdir()
                                    (ancestor / leaf.name).write_bytes(payload)
                            elif mutation_kind == "leaf_mutation":
                                leaf.write_bytes(b"X" * len(payload))
                            else:
                                leaf.rename(leaf.with_suffix(".parked"))
                                leaf.write_bytes(payload)

                        reader = self._readers(leaf, payload)[reader_name]
                        with self.assertRaises(plan.ExactFourPlanError):
                            self._during_first_read(mutate, reader)
                    finally:
                        temporary.cleanup()


class TestManifestAndScope(unittest.TestCase):
    def test_prerequisites_and_bundle_deterministic(self) -> None:
        first = local_mechanics_manifest()
        second = local_mechanics_manifest()
        self.assertEqual(canonical(first), canonical(second))
        self.assertEqual(first["targets"], [{"sample_id": s, "locus": plan.LOCUS, "part_count": c} for s, c in zip(plan.TARGETS, plan.PART_COUNTS)])
        self.assertEqual(sum(row["content_length_bytes"] for row in first["expected_provider_objects"]), 643_422_845_702)
        self.assertEqual(len(first["assembly_authorities"]), 8)
        self.assertTrue(all(value is False for value in first["action_gates"].values()))
        self.assertTrue(all(first[key] == expected for key, expected in plan.CONTRACT_BINDINGS.items()))

    def test_manifest_binds_complete_accepted_ratio_schema(self) -> None:
        frozen = local_mechanics_manifest()
        encoded = plan.canonical_json_bytes(frozen["coverage_ratio_schema_identity"])
        self.assertEqual(len(encoded), 2549)
        self.assertEqual(hashlib.sha256(encoded).hexdigest(), plan.COVERAGE_RATIO_SCHEMA_SHA256)
        self.assertEqual(frozen["scientific_laws"]["minimum_callable_flank_mean"], {"numerator": 5, "denominator": 1})
        forbidden = canonical(frozen).lower()
        for phrase in (b"coverage-compatible", b"coverage-discordant", b"ratio_min", b"ratio_max", b"compatibility_tolerance"):
            self.assertNotIn(phrase, forbidden)

    def test_frozen_manifest_is_current_canonical_builder_output(self) -> None:
        frozen_path = ROOT / "frozen" / "exact4_deployment_manifest.v1.json"
        raw = frozen_path.read_bytes()
        frozen = plan.strict_json_bytes(raw, "frozen deployment manifest")
        self.assertEqual(raw, canonical(frozen))
        self.assertEqual(frozen, plan.build_deployment_manifest(ROOT))

    def test_extra_missing_reordered_and_stale_scope_rejected(self) -> None:
        base = local_mechanics_manifest()
        for mutate in (
            lambda value: value["targets"].append({"sample_id": "HG00512", "locus": plan.LOCUS, "part_count": 1}),
            lambda value: value["targets"].pop(),
            lambda value: value["targets"].reverse(),
            lambda value: value["expected_provider_objects"].append(copy.deepcopy(value["expected_provider_objects"][0])),
        ):
            value = copy.deepcopy(base); mutate(value); value["bundle_id"] = ""; value["bundle_id"] = "exact4-bundle-" + plan.digest_value(plan._manifest_identity(value))[:24]
            with self.assertRaises(plan.ExactFourPlanError):
                plan._validate_deployment_manifest(value)


class TestPlanReleaseAndMutation(ExactFourFixture):
    def depth_fixture(self, first_depth: int = 8, second_depth: int = 24) -> tuple[Path, Path, Path]:
        sample = self.plan["analysis_map"][0]["sample_id"]
        beds = [row for row in self.plan["immutable_run_material"]["resident_identity"]["bed_authorities"] if row["sample_id"] == sample]
        rows = []
        paf_rows = []
        for index, bed in enumerate(beds):
            depth = first_depth if index == 0 else second_depth
            for start, end in bed["full_intervals"]:
                rows.extend(f"{bed['contig']}\t{position + 1}\t{depth}\n" for position in range(start, end))
                paf_rows.append(f"type2_KCON\t9472\t0\t500\t+\t{bed['contig']}:{start + 1}-{end}\t{end - start}\t200\t300\t480\t500\t60\n")
        raw = self.tmp / "raw.depth.tsv"; raw.write_text("".join(rows), encoding="ascii")
        mapq = self.tmp / "mapq10.depth.tsv"; mapq.write_text("".join(rows), encoding="ascii")
        paf = self.tmp / "type2_KCON.paf"; paf.write_text("".join(paf_rows), encoding="ascii")
        return raw, mapq, paf

    def test_plan_determinism_identity_and_false_gates(self) -> None:
        again = plan.build_target_only_plan(self.manifest, self.resident, self.provider)
        self.assertEqual(canonical(self.plan), canonical(again))
        immutable_sha = plan.digest_value(self.plan["immutable_run_material"])
        self.assertEqual(self.plan["run_id"], "exact4-" + immutable_sha[:24])
        self.assertEqual(len(self.plan["alignment_map"]), 13)
        self.assertEqual([row["sample_id"] for row in self.plan["alignment_map"]], [s for s, c in zip(plan.TARGETS, plan.PART_COUNTS) for _ in range(c)])
        self.assertTrue(self.plan["dry_run_only"])
        self.assertTrue(all(value is False for value in self.plan["action_gates"].values()))
        self.assertEqual(set(self.plan["action_gates"]), {
            "cluster_contact_performed", "network_contact_performed", "deployment_authorized", "deployment_performed",
            "download_authorized", "download_performed", "submission_authorized", "submission_performed",
            "execution_authorized", "execution_performed", "cancellation_authorized", "cancellation_performed",
            "retry_authorized", "retry_performed", "resubmission_authorized", "resubmission_performed",
            "data_access_authorized", "production_authorized", "copy_number_inference_authorized",
            "result_use_authorized", "model_use_authorized", "manuscript_use_authorized",
        })
        self.assertTrue(all(self.plan[key] == expected == self.plan["immutable_run_material"][key] for key, expected in plan.CONTRACT_BINDINGS.items()))

    def test_resident_nested_coordinate_and_identity_drift_rejected(self) -> None:
        mutations = []
        for mutate in (
            lambda value: value["element_identities"][0].__setitem__("source_start", 201),
            lambda value: value["bed_authorities"][0].__setitem__("full_intervals", [[0, 500]]),
            lambda value: value["assemblies"][0]["contigs"].__setitem__(next(iter(value["assemblies"][0]["contigs"])), 999),
            lambda value: value["assemblies"][0].__setitem__("unknown", "forbidden"),
            lambda value: value["sample_bed_authorities"].reverse(),
        ):
            bad = copy.deepcopy(self.resident); mutate(bad); mutations.append(bad)
        for bad in mutations:
            with self.assertRaises(plan.ExactFourPlanError):
                plan._validate_resident_receipt(bad, self.manifest)

    def test_nested_phase_map_drift_rejected_even_with_refreshed_plan_digest(self) -> None:
        mutations = (
            lambda value: value["preparation_map"][0]["haplotypes"].reverse(),
            lambda value: value["alignment_map"][0].__setitem__("part_index", 1),
            lambda value: value["alignment_map"][0]["provider_object"].__setitem__("uri", "s3://foreign/object"),
            lambda value: value["analysis_map"][0]["expected_global_part_indices"].reverse(),
            lambda value: value["preparation_map"][0].__setitem__("unknown", "forbidden"),
        )
        for mutate in mutations:
            bad = copy.deepcopy(self.plan); mutate(bad)
            material = dict(bad); material.pop("plan_digest")
            bad["plan_digest"] = plan.digest_value(material)
            with self.assertRaisesRegex(plan.ExactFourPlanError, "phase maps"):
                plan.validate_target_only_plan(bad)

    def test_provider_order_version_etag_and_repeat_drift(self) -> None:
        for field in ("version_id", "etag", "content_length_bytes", "uri"):
            bad = copy.deepcopy(self.provider)
            bad["objects"][3][field] = "bad" if field != "content_length_bytes" else 1
            with self.assertRaises(plan.ExactFourPlanError):
                plan._validate_provider_receipt(bad, self.manifest)
        bad = copy.deepcopy(self.provider); bad["objects"][0], bad["objects"][1] = bad["objects"][1], bad["objects"][0]
        with self.assertRaises(plan.ExactFourPlanError):
            plan._validate_provider_receipt(bad, self.manifest)

    def test_provider_v2_full_schema_false_gates_and_external_scope_refuse(self) -> None:
        self.assertEqual(set(self.provider), {"schema_version", "deployment_bundle_id", "objects", "provider_tool_receipt", "capture_metadata", "authority_bindings", "action_gates"})
        for gate in plan.ACTION_GATES:
            for bad_value in (True, None, 0, "false"):
                bad = copy.deepcopy(self.provider); bad["action_gates"][gate] = bad_value
                with self.assertRaises(plan.ExactFourPlanError):
                    plan._validate_provider_receipt(bad, self.manifest)
        for field, bad_value in (
            ("authorization_path", "handoff/wrong.md"),
            ("acceptance_path", "handoff/wrong-acceptance.md"),
            ("deployment_bundle_id", "exact4-bundle-cee202a9f2a520062573fe50"),
            ("deployment_manifest_path", "cluster_workflows/wrong.json"),
            ("provider_receipt_schema", "hml2_7p22_exact4_provider_receipt_1"),
            ("biological_target", "other"),
            ("ordered_samples", list(reversed(plan.TARGETS))),
            ("provider_part_counts", {sample: 1 for sample in plan.TARGETS}),
            ("read_only_action_scope", plan.AMENDED_READ_ONLY_SCOPE + "_expanded"),
        ):
            bad = copy.deepcopy(self.provider)
            bad["authority_bindings"]["amended_live_capture_authorization"][field] = bad_value
            with self.assertRaises(plan.ExactFourPlanError):
                plan._validate_provider_receipt(bad, self.manifest)

    def test_full_authorization_binding_changes_provider_plan_and_run_identity(self) -> None:
        first = self.plan
        changed = copy.deepcopy(self.provider)
        changed["authority_bindings"]["amended_live_capture_authorization"]["authorization_sha256"] = "f" * 64
        with mock.patch.object(plan, "_stable_reopen_pinned_file", side_effect=lambda relative, *_: canonical(self.manifest) if relative == plan.DEPLOYMENT_MANIFEST_PATH else b"x"):
            second = plan.build_target_only_plan(self.manifest, self.resident, changed)
        self.assertNotEqual(first["provider_receipt_sha256"], second["provider_receipt_sha256"])
        self.assertNotEqual(first["immutable_run_material"]["provider_identity"], second["immutable_run_material"]["provider_identity"])
        self.assertNotEqual(first["run_id"], second["run_id"])

    def test_release_wrong_plan_phase_attempt_and_destination_rejected(self) -> None:
        plan.validate_execution_release(self.release, self.plan, plan.sha256_file(self.plan_path))
        cases = [
            ("plan_file_sha256", "f" * 64), ("run_id", "exact4-" + "f" * 24),
            ("allowed_phases", ["preparation"]), ("attempt_number", 2),
            ("implementation_contract_acceptance_sha256", "f" * 64),
            ("coverage_ratio_schema_sha256", "f" * 64),
            ("production_authorized", True),
        ]
        for key, value in cases:
            bad = copy.deepcopy(self.release); bad[key] = value
            with self.assertRaises(plan.ExactFourPlanError):
                plan.validate_execution_release(bad, self.plan, plan.sha256_file(self.plan_path))
        bad = copy.deepcopy(self.release); bad["tool_receipts"] = [None] * 6
        with self.assertRaises(plan.ExactFourPlanError):
            plan.validate_execution_release(bad, self.plan, plan.sha256_file(self.plan_path))
        bad = copy.deepcopy(self.release); bad["tool_receipts"][0]["sha256"] = "f" * 64
        with self.assertRaises(plan.ExactFourPlanError):
            plan.validate_execution_release(bad, self.plan, plan.sha256_file(self.plan_path))
        tool_path = Path(self.release["tool_receipts"][0]["path"]); original_tool = tool_path.read_bytes()
        try:
            tool_path.write_bytes(b"drift")
            with self.assertRaisesRegex(plan.ExactFourPlanError, "current tool bytes drifted"):
                plan.validate_execution_release(self.release, self.plan, plan.sha256_file(self.plan_path))
        finally:
            tool_path.write_bytes(original_tool); tool_path.chmod(0o755)

    def test_claim_create_only_and_stale_release_before_write(self) -> None:
        claim = plan._create_claim(self.plan_path, self.release_path, "part", 0, "fixture-token")
        self.assertEqual(claim["kind"], "part")
        stage_root = Path(self.plan["alignment_map"][0]["staging_root"])
        self.assertTrue(stage_root.is_dir())
        (stage_root / "work").mkdir()
        with self.assertRaises(plan.ExactFourPlanError):
            plan._create_claim(self.plan_path, self.release_path, "part", 0, "fixture-token")
        stale = copy.deepcopy(self.release); stale["run_id"] = "exact4-" + "e" * 24
        stale_path = self.tmp / "stale.json"; stale_path.write_bytes(canonical(stale))
        before = set(self.tmp.rglob("*"))
        with self.assertRaises(plan.ExactFourPlanError):
            plan._create_claim(self.plan_path, stale_path, "part", 1, "stale-token")
        self.assertEqual(before, set(self.tmp.rglob("*")))

    def test_claim_gates_only_phase_relevant_live_sources_before_write(self) -> None:
        identity = self.plan["immutable_run_material"]["resident_identity"]
        plan._create_claim(self.plan_path, self.release_path, "prepare", 0, "prepare-source-gate")
        self.current_source_mock.assert_called_once_with(identity, "prepare", plan.TARGETS[0])
        self.current_source_mock.reset_mock()
        plan._create_claim(self.plan_path, self.release_path, "part", 0, "part-no-source-gate")
        self.current_source_mock.assert_not_called()
        plan._create_claim(self.plan_path, self.release_path, "analysis", 0, "analysis-source-gate")
        self.current_source_mock.assert_called_once_with(identity, "analysis", plan.TARGETS[0])

    def test_publish_file_no_clobber_symlink_hardlink_and_outside(self) -> None:
        mapping = self.plan["alignment_map"][0]
        stage = Path(mapping["staging_root"]); stage.mkdir(parents=True)
        source = stage / "candidate"; source.write_bytes(b"candidate")
        destination = Path(mapping["unfiltered_bam"])
        plan._publish(self.plan_path, self.release_path, source, destination)
        self.assertEqual(destination.read_bytes(), b"candidate")
        source2 = stage / "candidate2"; source2.write_bytes(b"other")
        with self.assertRaises(plan.ExactFourPlanError):
            plan._publish(self.plan_path, self.release_path, source2, destination)
        source3 = stage / "candidate3"; source3.write_bytes(b"linked")
        alias = stage / "alias"; os.link(source3, alias)
        with self.assertRaisesRegex(plan.ExactFourPlanError, "hard-link"):
            plan._publish(self.plan_path, self.release_path, source3, Path(mapping["primary_bam"]))
        with self.assertRaises(plan.ExactFourPlanError):
            plan._publish(self.plan_path, self.release_path, source2, self.tmp / "outside")

    def test_parent_symlink_lexical_dotdot_and_dangling_final_block(self) -> None:
        root = Path(self.plan["run_root"]); root.mkdir(parents=True, exist_ok=True)
        trap = root / "trap"; trap.symlink_to(self.tmp)
        with self.assertRaises(plan.ExactFourPlanError):
            plan._assert_existing_components_safe(root, trap / "x")
        with self.assertRaises(plan.ExactFourPlanError):
            plan._require_plan_path(self.plan, str(root / "a" / ".." / "x"), "poison")
        mapping = self.plan["alignment_map"][0]; stage = Path(mapping["staging_root"]); stage.mkdir(parents=True, exist_ok=True)
        source = stage / "dangling-source"; source.write_bytes(b"x")
        destination = Path(mapping["primary_bai"]); destination.parent.mkdir(parents=True, exist_ok=True); destination.symlink_to(self.tmp / "missing")
        with self.assertRaises(plan.ExactFourPlanError):
            plan._publish(self.plan_path, self.release_path, source, destination)

    def test_render_plan_is_pure(self) -> None:
        before = {str(path.relative_to(self.tmp)): (path.is_file(), hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else "") for path in self.tmp.rglob("*")}
        with mock.patch("subprocess.run", side_effect=AssertionError("process forbidden")), mock.patch("pathlib.Path.write_bytes", side_effect=AssertionError("write forbidden")):
            rendered = plan.build_target_only_plan(self.manifest, self.resident, self.provider)
        after = {str(path.relative_to(self.tmp)): (path.is_file(), hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else "") for path in self.tmp.rglob("*")}
        self.assertEqual(rendered, self.plan)
        self.assertEqual(before, after)

    def test_shell_rejects_unpinned_active_python_bash_before_mutation(self) -> None:
        before = {str(path.relative_to(self.tmp)) for path in self.tmp.rglob("*")}
        environment = dict(os.environ); environment["PYTHON"] = sys.executable
        completed = subprocess.run(
            [str(SRC / "02_exact4_alignment_worker.sh"), "--plan", str(self.plan_path), "--part-index", "0", "--execution-release", str(self.release_path), "--verify-only"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=environment,
            check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertEqual(before, {str(path.relative_to(self.tmp)) for path in self.tmp.rglob("*")})

    def test_plan_bound_summary_builder_exact_ratios_and_depth_order(self) -> None:
        raw, mapq, paf = self.depth_fixture()
        summary = plot._build_summary_from_plan(self.plan, 0, raw, mapq, paf)
        self.assertEqual(len(summary["haplotype_flank_records"]), 4)
        self.assertEqual(len(summary["element_ratio_records"]), 4)
        self.assertEqual([row["F_exact"] for row in summary["haplotype_flank_records"]], [plot._rational(8, 1), plot._rational(8, 1), plot._rational(24, 1), plot._rational(24, 1)])
        self.assertTrue(all(row["R_exact"] == plot._rational(1, 1) for row in summary["element_ratio_records"]))
        payload = mapq.read_text(encoding="ascii").splitlines()
        fields = payload[0].split("\t"); fields[2] = "9"; payload[0] = "\t".join(fields)
        mapq.write_text("\n".join(payload) + "\n", encoding="ascii")
        with self.assertRaisesRegex(plot.ExactFourPlotError, "raw depth is below MAPQ10"):
            plot._build_summary_from_plan(self.plan, 0, raw, mapq, paf)


class TestQnameAndScientificFixtures(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(); self.tmp = Path(self.temporary.name).resolve()
        self.manifest = local_mechanics_manifest()
        self.tools = fake_tools(self.tmp / "tools")
        self.receipt, self.path_map = synthetic_fixture_environment(
            self.manifest, self.tools, self.tmp / "fixture-environment"
        )
        self.patch_paths = mock.patch.object(
            plan, "_project_path_from_relative",
            side_effect=lambda relative: self.path_map.get(relative, PROJECT / relative),
        )
        self.patch_paths.start()
        self.patch_host = mock.patch.object(
            plan.os, "uname",
            return_value=type("Uname", (), {"nodename": "login-03.cluster.example"})(),
        )
        self.patch_host.start()

    def tearDown(self) -> None:
        self.patch_host.stop()
        self.patch_paths.stop()
        self.temporary.cleanup()

    def test_qname_inventory_and_cross_part_uniqueness(self) -> None:
        first, second = self.tmp / "a.txt", self.tmp / "b.txt"
        first.write_bytes(b"a\nb\n"); second.write_bytes(b"c\nd\n")
        self.assertEqual(plan.verify_qname_inventory_union([first, second])["count"], 4)
        second.write_bytes(b"b\nc\n")
        with self.assertRaisesRegex(plan.ExactFourPlanError, "cross-part"):
            plan.verify_qname_inventory_union([first, second])
        first.write_bytes(b"b\na\n")
        with self.assertRaises(plan.ExactFourPlanError):
            plan.inspect_qname_inventory(first)

    def test_production_worker_qname_cli_rejects_equal_count_wrong_name(self) -> None:
        source = self.tmp / "source.txt"; source.write_bytes(b"a\nb\n")
        decoded = self.tmp / "decoded.txt"; decoded.write_bytes(b"a\nb\n")
        aligned = self.tmp / "aligned.txt"; aligned.write_bytes(b"a\nc\n")
        supplementary = self.tmp / "supplementary.txt"; supplementary.write_bytes(b"")
        command = [sys.executable, str(SRC / "exact4_plan.py"), "validate-qname-ownership", "--source", str(source), "--decoded", str(decoded), "--aligned-primary", str(aligned), "--supplementary", str(supplementary)]
        completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        self.assertNotEqual(completed.returncode, 0)
        aligned.write_bytes(b"a\nb\n")
        completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        receipt = plan.strict_json_bytes(completed.stdout.encode("utf-8"), "QNAME CLI receipt")
        self.assertEqual(receipt["qname_set_cardinality"], 2)
        self.assertIn('validate-qname-ownership', (SRC / "02_exact4_alignment_worker.sh").read_text())

    def test_qname_seam_duplicate_reorder_malformed_and_unknown_supplementary(self) -> None:
        source = self.tmp / "s.txt"; decoded = self.tmp / "d.txt"; aligned = self.tmp / "a.txt"; supplementary = self.tmp / "u.txt"
        for path in (source, decoded, aligned): path.write_bytes(b"a\nb\n")
        supplementary.write_bytes(b"b\n")
        self.assertTrue(plan.validate_complete_qname_ownership(source, decoded, aligned, supplementary)["supplementary_qnames_subset"])
        for path, payload in (
            (decoded, b"b\na\n"), (decoded, b"a\na\n"), (decoded, b"a\n"),
            (decoded, b"a\nb\nc\n"), (decoded, b"a\tbad\nb\n"), (decoded, b"a\nb"),
            (aligned, b"a\n"), (aligned, b"a\nb\nc\n"), (aligned, b"a\na\n"),
            (aligned, b"a\nc\n"), (supplementary, b"c\n"),
        ):
            with self.subTest(path=path.name, payload=payload):
                originals = {item: item.read_bytes() for item in (source, decoded, aligned, supplementary)}
                path.write_bytes(payload)
                with self.assertRaises(plan.ExactFourPlanError):
                    plan.validate_complete_qname_ownership(source, decoded, aligned, supplementary)
                for item, original in originals.items(): item.write_bytes(original)
        decoded.write_bytes(b"b\na\n")
        decoded.write_bytes(b"".join(name + b"\n" for name in sorted(decoded.read_bytes().splitlines())))
        self.assertTrue(plan.validate_complete_qname_ownership(source, decoded, aligned, supplementary)["exact_qname_set_and_multiplicity_equal"])

    def test_seven_fixture_receipt_is_all_or_nothing_and_false_gated(self) -> None:
        receipt, tools = self.receipt, self.tools
        with mock.patch.object(plan.subprocess, "run", side_effect=AssertionError("replay subprocess forbidden")):
            self.assertEqual(plan.validate_resident_tool_fixture_receipt(receipt, tools), receipt)
        for index in range(7):
            bad = copy.deepcopy(receipt); bad["fixtures"][index]["passed"] = False
            with self.assertRaises(plan.ExactFourPlanError):
                plan.validate_resident_tool_fixture_receipt(bad, tools)
        bad = copy.deepcopy(receipt); bad["fixtures"][0]["commands"][0]["argv"] = ["/usr/bin/aws", "--version"]
        with self.assertRaises(plan.ExactFourPlanError):
            plan.validate_resident_tool_fixture_receipt(bad, tools)
        bad = copy.deepcopy(receipt); bad["action_gates"]["network_authorized"] = True
        with self.assertRaises(plan.ExactFourPlanError):
            plan.validate_resident_tool_fixture_receipt(bad, tools)

    def test_real_fixture_runner_requires_exact_separate_authorization_before_subprocess(self) -> None:
        tools = self.tools
        binding = self.receipt["authorization_binding"]
        with mock.patch.object(plan.subprocess, "run", side_effect=AssertionError("subprocess forbidden")) as runner, mock.patch.object(plan.os, "uname", return_value=type("Uname", (), {"nodename": "login-03.cluster.example"})()):
            with self.assertRaises(plan.ExactFourPlanError):
                plan.run_resident_tool_fixtures(
                    tools, self.manifest,
                    str(Path(binding["authorization_path"]).parent) + "/./" + Path(binding["authorization_path"]).name,
                    binding["authorization_size_bytes"], binding["authorization_sha256"],
                    binding["acceptance_path"], binding["acceptance_size_bytes"], binding["acceptance_sha256"],
                )
        runner.assert_not_called()
        with self.assertRaises(TypeError):
            plan.run_resident_tool_fixtures(tools, self.tmp / "caller-scratch", {"fixture_execution_authorized": True})

    def test_fixture_authorization_and_acceptance_aliases_all_refuse_at_zero_subprocess(self) -> None:
        binding = self.receipt["authorization_binding"]
        auth, acceptance = binding["authorization_path"], binding["acceptance_path"]
        cases = []
        for role, original in (("authorization", auth), ("acceptance", acceptance)):
            parent, name = original.rsplit("/", 1)
            cases.extend((
                (role, "relative"),
                (role, parent + "/./" + name),
                (role, parent + "/alias/../" + name),
                (role, original.replace("/", "//", 1)),
                (role, original + "/"),
            ))
        auth_link = self.tmp / "auth-link"; auth_link.symlink_to(auth)
        auth_copy = self.tmp / "auth-copy"; auth_copy.write_bytes(Path(auth).read_bytes())
        cases.extend((("authorization", str(auth_link)), ("authorization", str(auth_copy)), ("authorization", acceptance), ("acceptance", auth)))
        for role, alias in cases:
            with self.subTest(role=role, alias=alias), mock.patch.object(plan.subprocess, "run", side_effect=AssertionError("subprocess forbidden")) as runner:
                arguments = [self.tools, self.manifest, auth, binding["authorization_size_bytes"], binding["authorization_sha256"], acceptance, binding["acceptance_size_bytes"], binding["acceptance_sha256"]]
                arguments[2 if role == "authorization" else 5] = alias
                with self.assertRaises(plan.ExactFourPlanError):
                    plan.run_resident_tool_fixtures(*arguments)
                runner.assert_not_called()

    def test_fixture_cli_requires_explicit_external_authorities(self) -> None:
        help_text = subprocess.run([sys.executable, str(SRC / "exact4_plan.py"), "run-fixtures", "--help"], check=True, stdout=subprocess.PIPE, text=True).stdout
        for option in ("--fixture-authorization", "--fixture-authorization-acceptance", "--deployment-manifest", "--tools-receipt"):
            self.assertIn(option, help_text)
        with mock.patch.object(sys, "stderr", io.StringIO()), self.assertRaises(SystemExit):
            plan._parse_cli(["run-fixtures"])

    def test_process_cli_manifest_aliases_and_identical_copy_refuse_before_subprocess_or_write(self) -> None:
        manifest_path = self.path_map[plan.DEPLOYMENT_MANIFEST_PATH]
        tools_path = self.tmp / "tools-receipt.json"; tools_path.write_bytes(canonical({"tools": self.tools}))
        fixture_path = self.tmp / "fixture-receipt.json"; fixture_path.write_bytes(canonical(self.receipt))
        copy_path = self.tmp / "manifest-copy.json"; copy_path.write_bytes(manifest_path.read_bytes())
        link_path = self.tmp / "manifest-link.json"; link_path.symlink_to(manifest_path)
        self.assertEqual(plan._load_fixed_deployment_manifest_argument(str(manifest_path), "test manifest"), self.manifest)
        parent, name = str(manifest_path).rsplit("/", 1)
        aliases = (
            "relative-manifest.json", parent + "/./" + name,
            parent + "/alias/../" + name, str(manifest_path).replace("/", "//", 1),
            str(manifest_path) + "/", str(copy_path), str(link_path),
        )
        auth = self.receipt["authorization_binding"]
        tool_paths = {row["name"]: row["path"] for row in self.tools}
        for command in ("run-fixtures", "capture-resident"):
            for alias in aliases:
                with self.subTest(command=command, alias=alias):
                    if command == "run-fixtures":
                        argv = [command, "--deployment-manifest", alias, "--tools-receipt", str(tools_path), "--fixture-authorization", auth["authorization_path"], "--fixture-authorization-size-bytes", str(auth["authorization_size_bytes"]), "--fixture-authorization-sha256", auth["authorization_sha256"], "--fixture-authorization-acceptance", auth["acceptance_path"], "--fixture-authorization-acceptance-size-bytes", str(auth["acceptance_size_bytes"]), "--fixture-authorization-acceptance-sha256", auth["acceptance_sha256"]]
                    else:
                        argv = [command, "--deployment-manifest", alias, "--tool-fixture-receipt", str(fixture_path)]
                        for tool_name in ("bash", "python", "aws", "minimap2", "samtools", "sbatch"):
                            argv.extend((f"--{tool_name}", tool_paths[tool_name]))
                    parsed = plan._parse_cli(argv)
                    self.assertIs(type(parsed.deployment_manifest), str)
                    self.assertEqual(parsed.deployment_manifest, alias)
                    with mock.patch.object(plan.subprocess, "run", side_effect=AssertionError("subprocess forbidden")) as runner, mock.patch.object(Path, "write_bytes", side_effect=AssertionError("write forbidden")) as writer, mock.patch.object(sys, "stderr", io.StringIO()):
                        self.assertEqual(plan.main(argv), 2)
                    runner.assert_not_called(); writer.assert_not_called()

    def test_retained_artifact_replay_rejects_each_fixture_even_with_updated_hash(self) -> None:
        for index in range(7):
            with self.subTest(fixture=plan.FIXTURE_IDS[index]):
                bad = copy.deepcopy(self.receipt)
                row = bad["fixtures"][index]["input_files"][0] if bad["fixtures"][index]["input_files"] else bad["fixtures"][index]["output_files"][0]
                path = Path(row["path"])
                original = path.read_bytes()
                mutated = original + b"X"
                path.write_bytes(mutated)
                row["size_bytes"] = len(mutated)
                row["sha256"] = hashlib.sha256(mutated).hexdigest()
                with mock.patch.object(plan.subprocess, "run", side_effect=AssertionError("replay subprocess forbidden")):
                    with self.assertRaises((plan.ExactFourPlanError, plot.ExactFourPlotError, UnicodeError, ValueError)):
                        plan.validate_resident_tool_fixture_receipt(bad, self.tools)
                path.write_bytes(original)

    def test_fixture7_replay_uses_pinned_core_without_plan_or_receipt_recursion(self) -> None:
        original_validator = plan.validate_resident_tool_fixture_receipt
        fixture = self.receipt["fixtures"][6]
        with mock.patch.object(plan, "validate_target_only_plan", side_effect=AssertionError("target-plan recursion forbidden")), mock.patch.object(plan, "validate_resident_tool_fixture_receipt", side_effect=AssertionError("fixture-receipt recursion forbidden")), mock.patch.object(plan.subprocess, "run", side_effect=AssertionError("replay subprocess forbidden")):
            counts, assertions = plan._recompute_fixture_oracle(fixture, self.manifest)
            self.assertEqual(counts, fixture["record_counts"])
            self.assertEqual(assertions, fixture["semantic_assertions"])
            self.assertEqual(original_validator(self.receipt, self.tools), self.receipt)

    def test_fixture7_ignores_alternate_or_monkeypatched_import(self) -> None:
        fixture = self.receipt["fixtures"][6]
        with mock.patch.object(plot, "_select_kcon_annotations_core", side_effect=AssertionError("alternate module used")):
            self.assertEqual(plan._recompute_fixture_oracle(fixture, self.manifest)[0], fixture["record_counts"])

    def test_cross_fixture_output_swap_refuses(self) -> None:
        bad = copy.deepcopy(self.receipt)
        bad["fixtures"][0]["output_files"][0], bad["fixtures"][1]["output_files"][0] = bad["fixtures"][1]["output_files"][0], bad["fixtures"][0]["output_files"][0]
        with mock.patch.object(plan.subprocess, "run", side_effect=AssertionError("replay subprocess forbidden")):
            with self.assertRaises(plan.ExactFourPlanError):
                plan.validate_resident_tool_fixture_receipt(bad, self.tools)

    def test_fixture7_resident_minimap_output_is_required_and_authority_bound(self) -> None:
        bad = copy.deepcopy(self.receipt)
        row = next(item for item in bad["fixtures"][6]["output_files"] if Path(item["path"]).name == "resident_minimap.paf")
        path = Path(row["path"]); original = path.read_bytes()
        mutated = b"type2_KCON\t9472\t0\t9472\t+\twrong:1-12000\t12000\t500\t9972\t9472\t9472\t60\n"
        path.write_bytes(mutated); row["size_bytes"] = len(mutated); row["sha256"] = hashlib.sha256(mutated).hexdigest()
        with mock.patch.object(plan.subprocess, "run", side_effect=AssertionError("replay subprocess forbidden")):
            with self.assertRaises(ValueError):
                plan.validate_resident_tool_fixture_receipt(bad, self.tools)
        path.write_bytes(original)

    def test_kcon_production_core_exact_seams_and_ranking_tiebreakers(self) -> None:
        authority = plot._validate_kcon_selector_authority({
            "sample_id": "FIXTURE", "haplotypes": ["hapA", "hapB"],
            "bed_authorities": [{"haplotype": "hapA", "contig": "toy", "full_intervals": [[0, 12000]]}],
            "element_identities": [{"haplotype": "hapA", "contig": "toy", "source_start": 1500, "source_end": 2500}],
        })
        def row(query: str = "type2_KCON", qlen: int = 9472, qstart: int = 0, tname: str = "toy:1-12000", tstart: int = 1500, block: int = 500, matches: int = 500, mapq: int = 60, tag: str | None = None) -> tuple[list[str], bytes]:
            fields = [query, str(qlen), str(qstart), str(qstart + block), "+", tname, "12000", str(tstart), str(tstart + block), str(matches), str(block), str(mapq)]
            if tag is not None: fields.append(tag)
            raw = "\t".join(fields).encode("ascii")
            return fields, raw
        disqualified = [row(block=499, matches=499), row(tstart=3000), row(query="wrong"), row(qlen=9471), row(tname="wrong")]
        candidates = [
            row(qstart=50, tstart=1700, block=700, matches=680, mapq=50),
            row(qstart=40, tstart=1700, block=700, matches=690, mapq=40),
            row(qstart=30, tstart=1700, block=700, matches=690, mapq=50),
            row(qstart=20, tstart=1600, block=700, matches=690, mapq=50),
            row(qstart=10, tstart=1600, block=700, matches=690, mapq=50, tag="zz:Z:b"),
            row(qstart=10, tstart=1600, block=700, matches=690, mapq=50, tag="zz:Z:a"),
        ]
        all_rows = disqualified + candidates
        selected = plot._select_kcon_annotations_core(authority, all_rows, "a" * 64)
        self.assertEqual(selected[0]["query_start_0based"], 10)
        self.assertEqual(selected[0]["paf_row_sha256"], hashlib.sha256(candidates[-1][1]).hexdigest())
        self.assertEqual(selected, plan._independent_kcon_oracle(authority, all_rows, "a" * 64))

    def test_production_wrapper_validates_plan_and_parses_paf_once_then_delegates(self) -> None:
        paf = self.tmp / "wrapper.paf"
        paf.write_bytes(b"type2_KCON\t9472\t0\t500\t+\ttoy:1-12000\t12000\t1500\t2000\t500\t500\t60\n")
        authority = ("FIXTURE", ["hapA", "hapB"], [{"haplotype": "hapA", "contig": "toy", "full_intervals": [[0, 12000]]}], [{"haplotype": "hapA", "contig": "toy", "source_start": 1500, "source_end": 2500}])
        validated = {"validated": True}
        with mock.patch.object(plan, "validate_target_only_plan", return_value=validated) as validate_once, mock.patch.object(plot, "_declared_sample_authority", return_value=authority) as authority_once, mock.patch.object(plot, "_parse_paf", wraps=plot._parse_paf) as parse_once, mock.patch.object(plot, "_select_kcon_annotations_core", wraps=plot._select_kcon_annotations_core) as core_once:
            observed = plot._select_kcon_annotations_from_plan({}, 0, paf)
        self.assertEqual(len(observed), 1)
        validate_once.assert_called_once_with({})
        authority_once.assert_called_once_with(validated, 0)
        parse_once.assert_called_once_with(paf)
        core_once.assert_called_once()

    @staticmethod
    def rational(numerator: int, denominator: int = 1) -> dict[str, object]:
        return plot._rational(numerator, denominator)

    def make_summary(self) -> dict[str, object]:
        run_id, sample = "exact4-" + "a" * 24, "HG02027"
        flanks, elements = [], []
        for name, depth in (("h1", 8), ("h2", 24)):
            for product in plot.PRODUCTS:
                flank_key = [sample, plot.LOCUS, name, product]
                flanks.append({"schema_id": plot.FLANK_SCHEMA, "schema_version": 1, "run_id": run_id, "sample_id": sample, "locus_id": plot.LOCUS, "haplotype_id": name, "product": product, "flank_record_key": flank_key, "flank_depth_sum": depth * 4, "flank_base_count": 4, "F_exact": self.rational(depth), "callability": "callable", "callability_reason": "callable", "read_inferred_copy_number": "not_estimated"})
                element_key = [sample, plot.LOCUS, name, name, 2, 4, product]
                elements.append({"schema_id": plot.ELEMENT_SCHEMA, "schema_version": 1, "run_id": run_id, "sample_id": sample, "locus_id": plot.LOCUS, "haplotype_id": name, "contig": name, "source_start_0based": 2, "source_end_0based": 4, "product": product, "element_record_key": element_key, "flank_record_key": flank_key, "source_depth_sum": depth * 2, "source_base_count": 2, "S_exact": self.rational(depth), "R_exact": self.rational(1), "callability": "callable", "callability_reason": "callable", "read_inferred_copy_number": "not_estimated"})
        return {"schema_id": plot.SUMMARY_SCHEMA, "schema_version": 1, "run_id": run_id, "sample_id": sample, "locus_id": plot.LOCUS, "haplotype_flank_records": flanks, "element_ratio_records": elements, "read_inferred_copy_number": "not_estimated"}

    def test_unequal_haplotypes_unit_ratios_and_pdf_closure(self) -> None:
        summary = self.make_summary(); plot.validate_summary(summary)
        raw_path = self.tmp / "raw.tsv"; raw_path.write_bytes(b"h1\t1\t8\nh1\t2\t8\nh1\t3\t8\nh1\t4\t8\nh1\t5\t8\nh1\t6\t8\nh2\t1\t24\nh2\t2\t24\nh2\t3\t24\nh2\t4\t24\nh2\t5\t24\nh2\t6\t24\n")
        summary_path = self.tmp / "summary.json"; summary_path.write_bytes(canonical(summary))
        pdf = self.tmp / "raw.pdf"; plot.render_plot(raw_path, summary_path, pdf, "raw")
        plot.validate_pdf(pdf, {"summary": summary, "metric": "raw"})
        self.assertIn(b"read_inferred_copy_number=not_estimated", pdf.read_bytes())

    def test_ratio_two_remains_callable_descriptive_evidence(self) -> None:
        summary = self.make_summary()
        for element in summary["element_ratio_records"][:2]:
            element["source_depth_sum"] = 32
            element["S_exact"] = self.rational(16)
            element["R_exact"] = self.rational(2)
        plot.validate_summary(summary)
        self.assertTrue(all(element["callability"] == "callable" for element in summary["element_ratio_records"][:2]))
        self.assertNotIn("status", summary["element_ratio_records"][0])

    def test_zero_versus_ten_denominator_no_pooled_rescue(self) -> None:
        summary = self.make_summary()
        for flank in summary["haplotype_flank_records"][:2]:
            flank.update({"flank_depth_sum": 0, "F_exact": self.rational(0), "callability": "uncallable", "callability_reason": "own_flank_mean_below_5x"})
        for element in summary["element_ratio_records"][:2]:
            element.update({"source_depth_sum": 0, "S_exact": self.rational(0), "R_exact": None, "callability": "uncallable", "callability_reason": "own_flank_mean_below_5x"})
        plot.validate_summary(summary)
        self.assertTrue(all(element["R_exact"] is None for element in summary["element_ratio_records"][:2]))

    def test_equal_numeric_coordinates_distinct_contigs_remain_distinct(self) -> None:
        summary = self.make_summary(); plot.validate_summary(summary)
        identities = {(e["haplotype_id"], e["contig"], e["source_start_0based"], e["source_end_0based"]) for e in summary["element_ratio_records"]}
        self.assertEqual(len(identities), 2)

    def test_element_product_pair_and_declared_haplotype_order_rejected(self) -> None:
        mismatched = self.make_summary()
        record = mismatched["element_ratio_records"][1]
        record["contig"] = "h1b"
        record["element_record_key"] = [record["sample_id"], record["locus_id"], record["haplotype_id"], "h1b", record["source_start_0based"], record["source_end_0based"], record["product"]]
        with self.assertRaisesRegex(plot.ExactFourPlotError, "product pair"):
            plot.validate_summary(mismatched)

        reversed_summary = self.make_summary()
        reversed_summary["haplotype_flank_records"] = reversed_summary["haplotype_flank_records"][2:] + reversed_summary["haplotype_flank_records"][:2]
        reversed_summary["element_ratio_records"] = reversed_summary["element_ratio_records"][2:] + reversed_summary["element_ratio_records"][:2]
        with self.assertRaisesRegex(plot.ExactFourPlotError, "declared assembly authority"):
            plot.validate_summary(reversed_summary, expected_haplotype_order=["h1", "h2"])

    def test_kcon_paf_numeric_bounds_rejected(self) -> None:
        paf = self.tmp / "bad.paf"
        paf.write_text("type2_KCON\t9472\t-1\t99999\t+\th1:1-10\t10\t0\t10\t10\t10\t60\n", encoding="ascii")
        with self.assertRaisesRegex(plot.ExactFourPlotError, "malformed KCON PAF"):
            plot._parse_paf(paf)

    def test_empty_flank_and_exact_five_boundary(self) -> None:
        summary = self.make_summary()
        for flank in summary["haplotype_flank_records"][:2]:
            flank.update({"flank_depth_sum": 0, "flank_base_count": 0, "F_exact": None, "callability": "uncallable", "callability_reason": "empty_own_flank_set"})
        for element in summary["element_ratio_records"][:2]:
            element.update({"R_exact": None, "callability": "uncallable", "callability_reason": "empty_own_flank_set"})
        plot.validate_summary(summary)
        boundary = self.make_summary()
        for flank in boundary["haplotype_flank_records"][:2]:
            flank.update({"flank_depth_sum": 20, "F_exact": self.rational(5)})
        for element in boundary["element_ratio_records"][:2]:
            element.update({"source_depth_sum": 10, "S_exact": self.rational(5), "R_exact": self.rational(1)})
        plot.validate_summary(boundary)

    def test_one_versus_four_element_cardinality(self) -> None:
        summary = self.make_summary()
        additions = []
        template_raw, template_mapq = copy.deepcopy(summary["element_ratio_records"][2]), copy.deepcopy(summary["element_ratio_records"][3])
        for suffix in ("b", "c", "d"):
            for template in (template_raw, template_mapq):
                row = copy.deepcopy(template)
                row["contig"] = f"h2_{suffix}"
                row["element_record_key"] = [row["sample_id"], row["locus_id"], row["haplotype_id"], row["contig"], row["source_start_0based"], row["source_end_0based"], row["product"]]
                additions.append(row)
        summary["element_ratio_records"].extend(additions)
        summary["element_ratio_records"].sort(key=lambda row: ((0 if row["haplotype_id"] == "h1" else 1), row["contig"].encode(), row["source_start_0based"], row["source_end_0based"], plot.PRODUCTS.index(row["product"])))
        self.assertEqual(len(summary["element_ratio_records"]), 10)
        plot.validate_summary(summary)

    def test_nonreduced_bool_and_invalid_pdf_rejected(self) -> None:
        summary = self.make_summary(); summary["haplotype_flank_records"][0]["F_exact"] = {"schema_id": plot.RATIONAL_SCHEMA, "schema_version": 1, "numerator": 16, "denominator": 2}
        with self.assertRaisesRegex(plot.ExactFourPlotError, "canonical reduced"):
            plot.validate_summary(summary)
        summary = self.make_summary(); summary["haplotype_flank_records"][0]["flank_depth_sum"] = True
        with self.assertRaisesRegex(plot.ExactFourPlotError, "exact integer"):
            plot.validate_summary(summary)
        bad = self.tmp / "bad.pdf"; bad.write_bytes(b"%PDF-1.4\n")
        with self.assertRaises(plot.ExactFourPlotError):
            plot.validate_pdf(bad, {"summary": self.make_summary(), "metric": "raw"})


class TestStaticAndRegressionPins(unittest.TestCase):
    def test_shell_syntax_interfaces_and_no_control_capability(self) -> None:
        scripts = [SRC / "01_exact4_prepare.sh", SRC / "02_exact4_alignment_worker.sh", SRC / "03_exact4_analysis.sh"]
        for script in scripts:
            subprocess.run(["bash", "-n", str(script)], check=True)
            text = script.read_text()
            self.assertIn("--verify-only", text)
            for forbidden in ("s" + "batch", "s" + "cancel", "s" + "queue", "s" + "acct", "s" + "sh ", "s" + "cp "):
                self.assertNotIn(forbidden, text)
            self.assertNotIn("01_cnv_controller_v2.sh", text)
            self.assertNotIn("bait_min_aligned_bp", text)
            self.assertIn("PYTHONDONTWRITEBYTECODE=1", text)
            self.assertIn("sys.executable", text)
            self.assertIn('"$BASH"', text)
            self.assertNotIn("${PYTHON:-python3}", text)
            self.assertNotIn("mkdir -p", text)

    def test_alignment_and_analysis_exact_commands(self) -> None:
        worker = (SRC / "02_exact4_alignment_worker.sh").read_text()
        analysis = (SRC / "03_exact4_analysis.sh").read_text()
        self.assertIn("-ax map-ont --secondary=no", worker)
        self.assertIn("fasta -F 0x900", worker)
        self.assertIn("view -u -F 0x900 -L", worker)
        self.assertNotIn("downsample", worker.lower())
        self.assertIn("depth -a -g 0x600 -Q 10", analysis)
        self.assertIn("depth -a -g 0x600 -b", analysis)
        self.assertNotIn("depth -a -q 10", analysis)
        self.assertIn("-c -p 0.1 -N 10", analysis)

    def test_generic_classification_and_document_pins(self) -> None:
        expected = dict(plan.PREREQUISITE_PINS)
        drift = {relative: plan.sha256_file(PROJECT / relative) for relative, digest in expected.items() if plan.sha256_file(PROJECT / relative) != digest}
        self.assertEqual(drift, {})
        self.assertEqual(canonical(plan.build_deployment_manifest(ROOT)), canonical(local_mechanics_manifest()))
        self.assertEqual(plan.sha256_file(PROJECT / "cluster_workflows/cnv_v2/src/cnv_plan.py"), "392a9a8bc63e323c06f0abdb0d013e190082a7b0fed0fbfaba7f42480d75363c")

    def test_python_compile_and_cli_has_only_contract_commands(self) -> None:
        with tempfile.TemporaryDirectory() as pycache:
            environment = dict(os.environ)
            environment["PYTHONPYCACHEPREFIX"] = pycache
            subprocess.run([sys.executable, "-m", "py_compile", str(SRC / "exact4_plan.py"), str(SRC / "plot_exact4_depth.py")], check=True, env=environment)
        completed = subprocess.run([sys.executable, str(SRC / "exact4_plan.py"), "--help"], check=True, stdout=subprocess.PIPE, text=True)
        for command in ("freeze-bundle", "capture-resident", "render-plan", "validate-plan", "validate-release", "claim", "validate-prepared", "validate-part", "validate-analysis", "validate-run", "publish"):
            self.assertIn(command, completed.stdout)


if __name__ == "__main__":
    unittest.main()
