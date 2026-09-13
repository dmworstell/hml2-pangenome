#!/usr/bin/env python3
"""Fail-closed structural CNV/array unit planner and validator.

This program has no network, provider, download, cancellation, or scheduler
submission client.  Its only writes are create-only artifacts beneath the
manifest's new versioned run root.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import pathlib
import re
import shlex
import socket
import stat
import subprocess
import sys
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


sys.dont_write_bytecode = True

SCHEMA = "hml2.cnv-array-structural-closure-manifest.v1"
PREFLIGHT_SCHEMA = "hml2.cnv-array-structural-preflight.v1"
UNIT_SCHEMA = "hml2.cnv-array-structural-unit-checkpoint.v1"
CLOSURE_SCHEMA = "hml2.cnv-array-structural-run-closure.v1"
SHA_RE = re.compile(r"^[0-9a-f]{64}$")
SAMPLES = ("HG02027", "HG02178", "HG03669", "HG03834")
CLUSTER_PYTHON = (
    "/path/to/hml2_workspace/condaenv/"
    "repeatmaskerenv/bin/python3.11"
)
EXPECTED_TRUTH = {
    "HG02027": (("mat", "mat", 1), ("pat", "pat", 3)),
    "HG02178": (("hap1", "h1", 1), ("hap2", "h2", 4)),
    "HG03669": (("mat", "mat", 3), ("pat", "pat", 1)),
    "HG03834": (("mat", "mat", 1), ("pat", "pat", 3)),
}
EVIDENCE_COLUMNS = (
    "sample_id",
    "assembly_haplotype",
    "truth_haplotype",
    "array_state",
    "physical_copy_number",
    "assembly_bam_sha256",
    "assembly_bai_sha256",
    "indexed_validation_region",
    "indexed_window_record_count",
    "quickcheck_passed",
    "evidence_channel",
    "read_depth_claim",
    "copy_number_source",
    "physical_copy_policy",
)


class PlanError(RuntimeError):
    pass


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")


def hash_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hash_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def read_json(path: pathlib.Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PlanError("cannot read canonical JSON %s: %s" % (path, exc))
    if not isinstance(value, dict):
        raise PlanError("JSON root must be an object: %s" % path)
    return value


def require_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or not SHA_RE.fullmatch(value):
        raise PlanError("invalid SHA-256 for %s" % label)
    return value


def digest_without(record: Mapping[str, Any], field: str) -> str:
    material = dict(record)
    material.pop(field, None)
    return hash_bytes(canonical_bytes(material))


def _walk_strings(value: object) -> List[str]:
    out: List[str] = []
    if isinstance(value, str):
        out.append(value)
    elif isinstance(value, dict):
        for item in value.values():
            out.extend(_walk_strings(item))
    elif isinstance(value, list):
        for item in value:
            out.extend(_walk_strings(item))
    return out


def validate_manifest(manifest: Dict[str, Any], manifest_path: Optional[pathlib.Path] = None) -> Dict[str, Any]:
    if manifest.get("schema_version") != SCHEMA:
        raise PlanError("unexpected manifest schema")
    if require_sha(manifest.get("plan_sha256"), "plan_sha256") != digest_without(manifest, "plan_sha256"):
        raise PlanError("manifest self-digest mismatch")

    scope = manifest.get("scope")
    if not isinstance(scope, dict) or scope.get("ordered_samples") != list(SAMPLES):
        raise PlanError("exact-four sample order changed")
    if scope.get("locus") != "HML-2_7p22.1" or scope.get("reference_build") != "T2T-CHM13v2.0":
        raise PlanError("locus/reference scope changed")
    if scope.get("core_region_1based") != "chr7:4699540-4717514" or scope.get("validation_region_1based") != "chr7:4649540-4767514":
        raise PlanError("frozen coordinate scope changed")

    units = manifest.get("units")
    if not isinstance(units, list) or len(units) != 4:
        raise PlanError("manifest must contain four units")
    if [u.get("index") for u in units] != list(range(4)) or [u.get("sample_id") for u in units] != list(SAMPLES):
        raise PlanError("unit order/index mismatch")
    for unit in units:
        sample = unit["sample_id"]
        haps = unit.get("haplotypes")
        if not isinstance(haps, list) or len(haps) != 2:
            raise PlanError("each sample must bind exactly two haplotypes")
        observed = []
        for hap in haps:
            observed.append((hap.get("assembly_haplotype"), hap.get("truth_haplotype"), hap.get("physical_copy_number")))
            for key in ("assembly_fasta", "assembly_fai", "assembly_alignment_bam", "assembly_alignment_bai"):
                item = hap.get(key)
                if not isinstance(item, dict) or not isinstance(item.get("path"), str):
                    raise PlanError("missing %s for %s" % (key, sample))
                require_sha(item.get("sha256"), "%s/%s" % (sample, key))
            if hap["assembly_alignment_bam"].get("semantic_type") != "assembly_contig_to_t2t_alignment_not_read_depth":
                raise PlanError("assembly BAM cannot be presented as read depth")
            if not isinstance(hap.get("expected_indexed_window_record_count"), int) or hap["expected_indexed_window_record_count"] < 1:
                raise PlanError("indexed window must have a positive frozen record count")
        if tuple(observed) != EXPECTED_TRUTH[sample]:
            raise PlanError("exact copy truth or haplotype alias changed for %s" % sample)

    runtime = manifest.get("runtime")
    if not isinstance(runtime, dict) or runtime.get("os") != "Rocky Linux 9.6" or runtime.get("partition") != "batch":
        raise PlanError("Rocky9/batch runtime binding changed")
    modules = runtime.get("modules")
    if modules != ["samtools/1.21"] or any(str(x) != str(x).lower() for x in modules):
        raise PlanError("only lowercase samtools/1.21 is allowed")
    if runtime.get("required_environment", {}).get("APPTAINER_BIND") != "/cluster":
        raise PlanError("APPTAINER_BIND=/cluster is mandatory")
    if (runtime.get("required_environment", {}).get(
            "HML2_CLUSTER_PYTHON") != CLUSTER_PYTHON):
        raise PlanError("the shared CPython 3.11 binding changed")
    for name in ("bash", "python", "samtools", "sbatch"):
        tool = runtime.get("tools", {}).get(name)
        if not isinstance(tool, dict) or not isinstance(tool.get("path"), str):
            raise PlanError("missing runtime tool %s" % name)
        require_sha(tool.get("sha256"), "tool/%s" % name)
    if runtime["tools"]["python"].get("path") != CLUSTER_PYTHON:
        raise PlanError("runtime Python differs from HML2_CLUSTER_PYTHON")

    policy = manifest.get("execution_policy")
    if not isinstance(policy, dict):
        raise PlanError("execution policy missing")
    forbidden = (
        "network_allowed", "aws_allowed", "download_allowed", "minimap2_allowed",
        "new_reference_or_index_construction_allowed", "combined_reference_construction_allowed",
        "overwrite_allowed", "submission_authorized", "submission_tool_in_package",
    )
    if any(policy.get(key) is not False for key in forbidden):
        raise PlanError("a forbidden execution action is enabled")
    strings = [x.lower() for x in _walk_strings(manifest)]
    if any("afterany" in x for x in strings if x != "afterany_allowed"):
        # The only allowed occurrence is the literal policy key, not any value.
        raise PlanError("afterany is forbidden")
    scheduling = manifest.get("scheduling")
    if not isinstance(scheduling, dict) or scheduling.get("afterany_allowed") is not False:
        raise PlanError("afterany policy missing")
    if scheduling.get("unit_array") != {"array": "0-3%4", "task_count": 4, "cpus_per_task": 1, "memory_gib": 4, "wall_minutes": 20, "dependency": None}:
        raise PlanError("unit resources changed")
    if scheduling.get("closure") != {"task_count": 1, "cpus_per_task": 1, "memory_gib": 1, "wall_minutes": 5, "dependency": "receipt_census:unit_checkpoints"}:
        raise PlanError("closure resources/dependency changed")
    if scheduling.get("total_scheduled_tasks") != 5:
        raise PlanError("scheduled task count changed")

    authorities = manifest.get("authority_bindings")
    if not isinstance(authorities, dict):
        raise PlanError("authority bindings missing")
    source_roles = (
        "exact_copy_truth", "structural_reaudit_truth", "structural_reaudit_summary",
        "physical_copy_adjudication", "physical_copy_evidence", "cluster_assembly_authority",
    )
    for role in source_roles:
        item = authorities.get(role)
        if not isinstance(item, dict):
            raise PlanError("authority missing: %s" % role)
        require_sha(item.get("sha256"), "authority/%s" % role)
    physical_policy = authorities.get("policy", {})
    expected_policy = {
        "unit": "physical_copy",
        "seven_p22_tandem_units_count_once": True,
        "merge_with_paralog_or_segdup_families": False,
        "alignment_record_count_is_copy_number": False,
        "read_depth_inference_authorized": False,
    }
    if physical_policy != expected_policy:
        raise PlanError("physical-copy policy changed")

    storage = manifest.get("storage")
    if not isinstance(storage, dict):
        raise PlanError("storage policy missing")
    run_root = pathlib.PurePosixPath(str(storage.get("run_root")))
    if str(run_root) != "/path/to/hml2_workspace/HML2_project/cnv_array_v3/runs/run_20260718_003":
        raise PlanError("run root changed")
    expected_bundle = str(run_root / "bundle")
    expected_paths = {
        "bundle_root": expected_bundle,
        "manifest": expected_bundle + "/frozen/cnv_array_exact4_manifest.v1.json",
        "planner": expected_bundle + "/src/cnv_array_plan.py",
        "unit_script": expected_bundle + "/src/01_cnv_array_unit.sh",
        "closure_script": expected_bundle + "/src/02_cnv_array_closure.sh",
        "receipt_controller": expected_bundle + "/src/03_cnv_array_receipt_controller.sh",
    }
    if manifest.get("execution_paths") != expected_paths:
        raise PlanError("absolute deployed execution paths changed")
    protected = [pathlib.PurePosixPath(str(x)) for x in storage.get("protected_roots", [])]
    if len(protected) != 4:
        raise PlanError("protected result roots changed")
    destinations = [storage.get("preflight_receipt"), storage.get("closure_tsv"), storage.get("closure_json")]
    for unit in units:
        destinations.extend(unit["outputs"].values())
    for raw in destinations:
        dest = pathlib.PurePosixPath(str(raw))
        if run_root not in dest.parents:
            raise PlanError("destination escapes run root: %s" % dest)
        for root in protected:
            if dest == root or root in dest.parents:
                raise PlanError("destination overlaps protected results: %s" % dest)

    if manifest_path is not None:
        package_root = manifest_path.resolve().parent.parent
        for item in manifest.get("bundle_files", []):
            path = package_root / item["path"]
            _validate_regular_file(path, int(item["size_bytes"]), require_sha(item["sha256"], "bundle file"), allow_symlink=False)
    return manifest


def load_manifest(path: pathlib.Path) -> Dict[str, Any]:
    return validate_manifest(read_json(path), path)


def _validate_regular_file(path: pathlib.Path, expected_size: Optional[int], expected_sha: str, allow_symlink: bool = False) -> Dict[str, Any]:
    try:
        lst = path.lstat()
    except OSError as exc:
        raise PlanError("missing input %s: %s" % (path, exc))
    if stat.S_ISLNK(lst.st_mode) and not allow_symlink:
        raise PlanError("symlink input forbidden: %s" % path)
    try:
        st = path.stat()
    except OSError as exc:
        raise PlanError("cannot stat input %s: %s" % (path, exc))
    if not stat.S_ISREG(st.st_mode):
        raise PlanError("input is not a regular file: %s" % path)
    if expected_size is not None and st.st_size != expected_size:
        raise PlanError("size mismatch: %s" % path)
    observed = hash_file(path)
    if observed != expected_sha:
        raise PlanError("SHA-256 mismatch: %s" % path)
    return {"path": str(path), "size_bytes": st.st_size, "sha256": observed, "device": st.st_dev, "inode": st.st_ino}


def _validate_tool(tool: Mapping[str, Any], name: str) -> Dict[str, Any]:
    path = pathlib.Path(str(tool["path"]))
    item = _validate_regular_file(path, None, require_sha(tool["sha256"], "tool/%s" % name), allow_symlink=(name == "python"))
    item["realpath"] = str(path.resolve())
    if name == "python" and item["realpath"] != tool.get("realpath"):
        raise PlanError("Python realpath mismatch")
    return item


def _samtools_env(manifest: Mapping[str, Any]) -> Dict[str, str]:
    env = dict(os.environ)
    for key, value in manifest["runtime"]["required_environment"].items():
        env[str(key)] = str(value)
    return env


def _run_samtools(manifest: Mapping[str, Any], args: Sequence[str], capture: bool = True) -> subprocess.CompletedProcess:
    executable = str(manifest["runtime"]["tools"]["samtools"]["path"])
    proc = subprocess.run([executable] + list(args), env=_samtools_env(manifest), stdout=subprocess.PIPE if capture else subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, check=False)
    if proc.returncode != 0:
        raise PlanError("samtools failed (%s): %s" % (" ".join(args), proc.stderr.strip()))
    return proc


def _validate_haplotype_live(manifest: Mapping[str, Any], hap: Mapping[str, Any]) -> Dict[str, Any]:
    bam = hap["assembly_alignment_bam"]
    bai = hap["assembly_alignment_bai"]
    bam_obs = _validate_regular_file(pathlib.Path(bam["path"]), int(bam["size_bytes"]), require_sha(bam["sha256"], "bam"))
    bai_obs = _validate_regular_file(pathlib.Path(bai["path"]), int(bai["size_bytes"]), require_sha(bai["sha256"], "bai"))
    _run_samtools(manifest, ["quickcheck", "-v", bam["path"]])
    region = str(manifest["scope"]["validation_region_1based"])
    proc = _run_samtools(manifest, ["view", "-c", bam["path"], region])
    try:
        count = int(proc.stdout.strip())
    except ValueError:
        raise PlanError("samtools returned a non-integer region count")
    if count != hap["expected_indexed_window_record_count"]:
        raise PlanError("indexed window count changed for %s" % bam["path"])
    return {"bam": bam_obs, "bai": bai_obs, "quickcheck_passed": True, "indexed_window_record_count": count}


def _record_digest(record: Dict[str, Any], field: str) -> Dict[str, Any]:
    record[field] = ""
    record[field] = digest_without(record, field)
    return record


def _validate_record_digest(record: Mapping[str, Any], schema: str, field: str) -> None:
    if record.get("schema_version") != schema:
        raise PlanError("unexpected receipt schema")
    if require_sha(record.get(field), field) != digest_without(record, field):
        raise PlanError("receipt self-digest mismatch")


def _write_exclusive(path: pathlib.Path, payload: bytes) -> None:
    if path.exists() or path.is_symlink():
        raise PlanError("create-only destination already exists: %s" % path)
    with path.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def _preflight_path(manifest: Mapping[str, Any]) -> pathlib.Path:
    return pathlib.Path(str(manifest["storage"]["preflight_receipt"]))


def validate_preflight(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    path = _preflight_path(manifest)
    record = read_json(path)
    _validate_record_digest(record, PREFLIGHT_SCHEMA, "receipt_sha256")
    if record.get("run_id") != manifest["run_id"] or record.get("plan_sha256") != manifest["plan_sha256"]:
        raise PlanError("preflight is for another run/plan")
    if record.get("all_eight_quickcheck_and_index_queries_passed") is not True or len(record.get("haplotypes", [])) != 8:
        raise PlanError("preflight is incomplete")
    return record


def preflight(manifest_path: pathlib.Path) -> Dict[str, Any]:
    manifest = load_manifest(manifest_path)
    receipt_path = _preflight_path(manifest)
    if receipt_path.exists():
        return validate_preflight(manifest)
    if socket.getfqdn() != manifest["runtime"]["observed_login_fqdn"]:
        raise PlanError("preflight must run on the authenticated upgraded login FQDN")
    tools = {name: _validate_tool(item, name) for name, item in manifest["runtime"]["tools"].items()}
    version = _run_samtools(manifest, ["--version"]).stdout.splitlines()[0]
    if version != "samtools 1.21":
        raise PlanError("samtools version mismatch")
    authority = manifest["authority_bindings"]["cluster_assembly_authority"]
    authority_obs = _validate_regular_file(pathlib.Path(authority["path"]), int(authority["size_bytes"]), require_sha(authority["sha256"], "cluster assembly authority"))
    observations: List[Dict[str, Any]] = []
    for unit in manifest["units"]:
        for hap in unit["haplotypes"]:
            live = _validate_haplotype_live(manifest, hap)
            observations.append({"sample_id": unit["sample_id"], "assembly_haplotype": hap["assembly_haplotype"], "truth_haplotype": hap["truth_haplotype"], "observation": live})
    run_root = pathlib.Path(str(manifest["storage"]["run_root"]))
    if not run_root.is_dir() or run_root.is_symlink():
        raise PlanError("versioned run root must already exist as a real directory")
    for name in ("units",):
        child = run_root / name
        if child.exists() or child.is_symlink():
            raise PlanError("preflight child already exists: %s" % child)
        child.mkdir(mode=0o750)
    receipt_path.parent.mkdir(mode=0o750)
    receipt = _record_digest({
        "schema_version": PREFLIGHT_SCHEMA,
        "run_id": manifest["run_id"],
        "plan_sha256": manifest["plan_sha256"],
        "login_fqdn": socket.getfqdn(),
        "os": manifest["runtime"]["os"],
        "partition": manifest["runtime"]["partition"],
        "required_environment": manifest["runtime"]["required_environment"],
        "tools": tools,
        "samtools_version": version,
        "cluster_assembly_authority": authority_obs,
        "haplotypes": observations,
        "all_eight_quickcheck_and_index_queries_passed": True,
        "network_provider_download_index_build_performed": False,
        "submission_performed": False,
    }, "receipt_sha256")
    _write_exclusive(receipt_path, canonical_bytes(receipt))
    return validate_preflight(manifest)


def _evidence_bytes(manifest: Mapping[str, Any], unit: Mapping[str, Any], observations: Sequence[Mapping[str, Any]]) -> bytes:
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(EVIDENCE_COLUMNS), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for hap, live in zip(unit["haplotypes"], observations):
        writer.writerow({
            "sample_id": unit["sample_id"],
            "assembly_haplotype": hap["assembly_haplotype"],
            "truth_haplotype": hap["truth_haplotype"],
            "array_state": hap["array_state"],
            "physical_copy_number": hap["physical_copy_number"],
            "assembly_bam_sha256": live["bam"]["sha256"],
            "assembly_bai_sha256": live["bai"]["sha256"],
            "indexed_validation_region": manifest["scope"]["validation_region_1based"],
            "indexed_window_record_count": live["indexed_window_record_count"],
            "quickcheck_passed": "true",
            "evidence_channel": "assembly_contig_to_t2t_structural_index_validation",
            "read_depth_claim": "none",
            "copy_number_source": "accepted_assembly_resolved_truth_not_alignment_count",
            "physical_copy_policy": "7p22_tandem_units_count_once_no_paralog_family_merge",
        })
    return out.getvalue().encode("utf-8")


def _unit_paths(unit: Mapping[str, Any]) -> Tuple[pathlib.Path, pathlib.Path]:
    return pathlib.Path(unit["outputs"]["evidence_tsv"]), pathlib.Path(unit["outputs"]["checkpoint_json"])


def validate_unit(manifest: Mapping[str, Any], index: int) -> Dict[str, Any]:
    unit = manifest["units"][index]
    evidence_path, checkpoint_path = _unit_paths(unit)
    checkpoint = read_json(checkpoint_path)
    _validate_record_digest(checkpoint, UNIT_SCHEMA, "checkpoint_sha256")
    if checkpoint.get("run_id") != manifest["run_id"] or checkpoint.get("plan_sha256") != manifest["plan_sha256"]:
        raise PlanError("unit checkpoint is for another run/plan")
    if checkpoint.get("unit_index") != index or checkpoint.get("sample_id") != unit["sample_id"]:
        raise PlanError("unit checkpoint identity mismatch")
    observed = _validate_regular_file(evidence_path, checkpoint["output"]["size_bytes"], require_sha(checkpoint["output"]["sha256"], "unit output"))
    if checkpoint.get("read_depth_claim") != "none" or checkpoint.get("copy_number_inferred_from_alignment_count") is not False:
        raise PlanError("unit overclaims alignment evidence")
    checkpoint["validated_output"] = observed
    return checkpoint


def run_unit(manifest_path: pathlib.Path, index: int) -> Dict[str, Any]:
    manifest = load_manifest(manifest_path)
    if index < 0 or index >= len(manifest["units"]):
        raise PlanError("unit index out of range")
    validate_preflight(manifest)
    unit = manifest["units"][index]
    evidence_path, checkpoint_path = _unit_paths(unit)
    if checkpoint_path.exists() or checkpoint_path.is_symlink():
        # Successful reruns validate current input identity before becoming a no-op.
        for hap in unit["haplotypes"]:
            _validate_haplotype_live(manifest, hap)
        return validate_unit(manifest, index)
    if evidence_path.exists() or evidence_path.is_symlink() or evidence_path.parent.exists() or evidence_path.parent.is_symlink():
        raise PlanError("partial/ambiguous unit output blocks create-only rerun")
    observations = [_validate_haplotype_live(manifest, hap) for hap in unit["haplotypes"]]
    payload = _evidence_bytes(manifest, unit, observations)
    units_root = pathlib.Path(str(manifest["storage"]["run_root"])) / "units"
    stage = units_root / (".staging-%s" % unit["sample_id"])
    if stage.exists() or stage.is_symlink():
        raise PlanError("stale staging directory blocks unit")
    stage.mkdir(mode=0o750)
    staged_evidence = stage / evidence_path.name
    staged_checkpoint = stage / checkpoint_path.name
    _write_exclusive(staged_evidence, payload)
    checkpoint = _record_digest({
        "schema_version": UNIT_SCHEMA,
        "run_id": manifest["run_id"],
        "plan_sha256": manifest["plan_sha256"],
        "unit_index": index,
        "sample_id": unit["sample_id"],
        "truth": [{"assembly_haplotype": h["assembly_haplotype"], "truth_haplotype": h["truth_haplotype"], "physical_copy_number": h["physical_copy_number"], "array_state": h["array_state"]} for h in unit["haplotypes"]],
        "input_observations": observations,
        "output": {"path": str(evidence_path), "size_bytes": len(payload), "sha256": hash_bytes(payload)},
        "read_depth_claim": "none",
        "copy_number_inferred_from_alignment_count": False,
        "physical_copy_policy": manifest["authority_bindings"]["policy"],
        "network_provider_download_index_build_performed": False,
    }, "checkpoint_sha256")
    _write_exclusive(staged_checkpoint, canonical_bytes(checkpoint))
    os.rename(str(stage), str(evidence_path.parent))
    return validate_unit(manifest, index)


def _read_evidence(path: pathlib.Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != EVIDENCE_COLUMNS:
            raise PlanError("evidence columns changed: %s" % path)
        rows = list(reader)
    if len(rows) != 2:
        raise PlanError("each unit evidence must contain two rows")
    return rows


def validate_closure(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    closure_path = pathlib.Path(str(manifest["storage"]["closure_json"]))
    tsv_path = pathlib.Path(str(manifest["storage"]["closure_tsv"]))
    closure = read_json(closure_path)
    _validate_record_digest(closure, CLOSURE_SCHEMA, "closure_sha256")
    if closure.get("run_id") != manifest["run_id"] or closure.get("plan_sha256") != manifest["plan_sha256"]:
        raise PlanError("closure is for another run/plan")
    _validate_regular_file(tsv_path, closure["combined_output"]["size_bytes"], require_sha(closure["combined_output"]["sha256"], "closure TSV"))
    if closure.get("sample_count") != 4 or closure.get("haplotype_count") != 8 or closure.get("read_depth_claim") != "none":
        raise PlanError("closure scope/evidence claim changed")
    return closure


def close_run(manifest_path: pathlib.Path) -> Dict[str, Any]:
    manifest = load_manifest(manifest_path)
    validate_preflight(manifest)
    closure_path = pathlib.Path(str(manifest["storage"]["closure_json"]))
    tsv_path = pathlib.Path(str(manifest["storage"]["closure_tsv"]))
    if closure_path.exists() or closure_path.is_symlink():
        return validate_closure(manifest)
    if tsv_path.exists() or tsv_path.is_symlink() or closure_path.parent.exists() or closure_path.parent.is_symlink():
        raise PlanError("partial/ambiguous closure blocks create-only rerun")
    checkpoints = [validate_unit(manifest, i) for i in range(4)]
    rows: List[Dict[str, str]] = []
    for unit in manifest["units"]:
        evidence_path, _ = _unit_paths(unit)
        rows.extend(_read_evidence(evidence_path))
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(EVIDENCE_COLUMNS), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    payload = out.getvalue().encode("utf-8")
    closure_path.parent.mkdir(mode=0o750)
    _write_exclusive(tsv_path, payload)
    closure = _record_digest({
        "schema_version": CLOSURE_SCHEMA,
        "run_id": manifest["run_id"],
        "plan_sha256": manifest["plan_sha256"],
        "sample_count": 4,
        "haplotype_count": 8,
        "unit_checkpoints": [{"sample_id": c["sample_id"], "path": manifest["units"][i]["outputs"]["checkpoint_json"], "sha256": hash_file(pathlib.Path(manifest["units"][i]["outputs"]["checkpoint_json"]))} for i, c in enumerate(checkpoints)],
        "combined_output": {"path": str(tsv_path), "size_bytes": len(payload), "sha256": hash_bytes(payload)},
        "copy_numbers": {row["sample_id"] + ":" + row["truth_haplotype"]: int(row["physical_copy_number"]) for row in rows},
        "read_depth_claim": "none",
        "copy_number_inferred_from_alignment_count": False,
        "physical_copy_policy": manifest["authority_bindings"]["policy"],
        "historical_results_overwritten": False,
        "network_provider_download_index_build_performed": False,
    }, "closure_sha256")
    _write_exclusive(closure_path, canonical_bytes(closure))
    return validate_closure(manifest)


def render_sbatch(manifest_path: pathlib.Path) -> str:
    manifest = load_manifest(manifest_path)
    paths = manifest["execution_paths"]
    unit_script = paths["unit_script"]
    closure_script = paths["closure_script"]
    receipt_controller = paths["receipt_controller"]
    deployed_manifest = paths["manifest"]
    unit = manifest["scheduling"]["unit_array"]
    closure = manifest["scheduling"]["closure"]
    cluster_python = manifest["runtime"]["required_environment"][
        "HML2_CLUSTER_PYTHON"]
    export_python = "HML2_CLUSTER_PYTHON=%s" % cluster_python
    first = [
        "/usr/bin/sbatch", "--parsable", "--partition=batch", "--array=%s" % unit["array"],
        "--cpus-per-task=%s" % unit["cpus_per_task"], "--mem=%sG" % unit["memory_gib"],
        "--time=00:%02d:00" % unit["wall_minutes"], "--job-name=cnv7p22_struct4",
        "--export=%s" % export_python, str(unit_script), deployed_manifest,
    ]
    second = [
        "/usr/bin/sbatch", "--parsable", "--partition=batch",
        "--cpus-per-task=%s" % closure["cpus_per_task"], "--mem=%sG" % closure["memory_gib"],
        "--time=00:%02d:00" % closure["wall_minutes"], "--job-name=cnv7p22_close",
        "--export=%s" % export_python, str(receipt_controller),
        deployed_manifest, str(closure_script),
    ]
    return "\n".join(" ".join(shlex.quote(x) for x in command) for command in (first, second)) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("validate-manifest", "preflight", "run-unit", "validate-unit", "close-run", "validate-closure", "render-sbatch"):
        child = sub.add_parser(name)
        child.add_argument("--manifest", type=pathlib.Path, required=True)
        if name in ("run-unit", "validate-unit"):
            child.add_argument("--unit-index", type=int, required=True)
    args = parser.parse_args(argv)
    if args.command == "validate-manifest":
        result = load_manifest(args.manifest)
    elif args.command == "preflight":
        result = preflight(args.manifest)
    elif args.command == "run-unit":
        result = run_unit(args.manifest, args.unit_index)
    elif args.command == "validate-unit":
        result = validate_unit(load_manifest(args.manifest), args.unit_index)
    elif args.command == "close-run":
        result = close_run(args.manifest)
    elif args.command == "validate-closure":
        result = validate_closure(load_manifest(args.manifest))
    elif args.command == "render-sbatch":
        sys.stdout.write(render_sbatch(args.manifest))
        return 0
    else:
        raise PlanError("unknown command")
    sys.stdout.write(canonical_bytes(result).decode("utf-8"))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PlanError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
