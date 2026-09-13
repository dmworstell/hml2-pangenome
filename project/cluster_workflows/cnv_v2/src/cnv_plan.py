#!/usr/bin/env python3
"""Strict, target-bound planning and receipt validation for CNV workflow v2.

The default ``render`` operation is deliberately offline: it validates a TSV
manifest and immutable local reference inputs, then writes deterministic plans
and BED files.  It never submits work or invokes an aligner.  ``prepare`` is a
separate, explicit operation that builds the run-scoped combined reference.
"""

from __future__ import annotations

import argparse
import csv
import errno
import heapq
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_VERSION = "cnv_v2_plan_1"
PREPARED_SCHEMA = "cnv_v2_prepared_1"
PART_RECEIPT_SCHEMA = "cnv_v2_part_done_1"
CLAIM_SCHEMA = "cnv_v2_claim_1"
QNAME_UNION_SCHEMA = "cnv_v2_qname_union_1"
QNAME_ENCODING = "SAM-QNAME ASCII bytes, LF terminated"
QNAME_ORDER = "LC_ALL=C bytewise ascending, unique"
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")

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

TARGET_FIELDS = tuple(
    field
    for field in MANIFEST_FIELDS
    if field
    not in {
        "part_id",
        "read_uri",
        "read_version",
        "read_etag",
        "read_size",
        "read_sha256",
        "read_format",
        "minimap2_preset",
    }
)

DEFAULT_CODE_NAMES = (
    "cnv_plan.py",
    "01_cnv_controller_v2.sh",
    "02_cnv_alignment_worker_v2.sh",
    "03_cnv_analysis_v2.sh",
    "plot_cnv_depth_v2.py",
)


class PlanError(RuntimeError):
    """A fail-closed manifest, plan, or receipt validation error."""


def _literal_absolute_path(path: os.PathLike[str] | str, label: str) -> Path:
    """Normalize lexical components without following the final path entry."""
    expanded = os.path.expanduser(os.fspath(path))
    if not os.path.isabs(expanded):
        raise PlanError(f"{label} must be absolute")
    return Path(os.path.normpath(expanded))


def _resolved_target_root(plan: Mapping[str, Any]) -> Path:
    literal = _literal_absolute_path(str(plan.get("target_dir", "")), "plan target_dir")
    try:
        return literal.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise PlanError(f"plan target_dir cannot be resolved: {literal}: {error}") from error


def _resolved_parent_with_literal_basename(
    path: os.PathLike[str] | str,
    label: str,
) -> Path:
    literal = _literal_absolute_path(path, label)
    try:
        resolved_parent = literal.parent.resolve(strict=False)
    except (OSError, RuntimeError) as error:
        raise PlanError(f"{label} parent cannot be resolved: {literal.parent}: {error}") from error
    return resolved_parent / literal.name


def _validate_plan_bound_destination(
    plan: Mapping[str, Any],
    path: os.PathLike[str] | str,
    label: str,
) -> Path:
    root = _resolved_target_root(plan)
    candidate = _resolved_parent_with_literal_basename(path, label)
    try:
        within = os.path.commonpath((str(root), str(candidate))) == str(root)
    except ValueError as error:
        raise PlanError(f"{label} is not on the plan target path domain: {candidate}") from error
    if not within:
        raise PlanError(
            f"{label} resolves outside plan target_dir through an existing parent component: {candidate}"
        )
    return candidate


def _canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode()


def digest_value(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def sha256_file(path: os.PathLike[str] | str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path: os.PathLike[str] | str) -> dict[str, Any]:
    item = Path(path)
    if not item.is_file():
        raise PlanError(f"required file is missing: {item}")
    size = item.stat().st_size
    if size <= 0:
        raise PlanError(f"required file is empty: {item}")
    return {"path": str(item.resolve()), "size": size, "sha256": sha256_file(item)}


def fingerprint_allow_empty(path: os.PathLike[str] | str) -> dict[str, Any]:
    """Fingerprint a required regular file whose empty byte string is meaningful."""
    item = Path(path)
    if not item.is_file():
        raise PlanError(f"required file is missing: {item}")
    return {
        "path": str(item.resolve()),
        "size": item.stat().st_size,
        "sha256": sha256_file(item),
    }


def _require_sha(value: str, label: str) -> str:
    lowered = value.lower()
    if not SHA256.fullmatch(lowered):
        raise PlanError(f"{label} must be a 64-character SHA-256 hex digest")
    return lowered


def _require_id(value: str, label: str) -> str:
    if not SAFE_ID.fullmatch(value):
        raise PlanError(f"{label} is not a safe identifier: {value!r}")
    return value


def _positive_int(value: str, label: str, *, allow_zero: bool = False) -> int:
    if not re.fullmatch(r"[0-9]+", value):
        raise PlanError(f"{label} must be an unsigned decimal integer")
    parsed = int(value)
    if parsed < (0 if allow_zero else 1):
        raise PlanError(f"{label} is out of range")
    return parsed


def _absolute_file(value: str, label: str) -> Path:
    item = Path(value)
    if not item.is_absolute():
        raise PlanError(f"{label} must be an absolute path: {value!r}")
    if not item.is_file():
        raise PlanError(f"{label} does not exist: {value!r}")
    if item.stat().st_size <= 0:
        raise PlanError(f"{label} is empty: {value!r}")
    return item.resolve()


def _verify_file(path: str, expected_sha256: str, label: str) -> Path:
    item = _absolute_file(path, label)
    expected = _require_sha(expected_sha256, f"{label} SHA-256")
    observed = sha256_file(item)
    if observed != expected:
        raise PlanError(f"{label} SHA-256 mismatch: expected {expected}, observed {observed}")
    return item


def _read_fai(path: Path) -> dict[str, int]:
    records: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            fields = line.split("\t")
            if len(fields) < 5:
                raise PlanError(f"malformed FAI line {line_number}: {path}")
            name = fields[0]
            if not name or name in records:
                raise PlanError(f"empty or duplicate FAI contig {name!r}: {path}")
            records[name] = _positive_int(fields[1], f"FAI length at {path}:{line_number}")
    if not records:
        raise PlanError(f"FAI contains no records: {path}")
    return records


def _fasta_lengths(path: Path) -> dict[str, int]:
    """Stream FASTA identifiers and lengths without retaining sequences."""
    lengths: dict[str, int] = {}
    current: str | None = None
    length = 0
    with path.open("r", encoding="ascii") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = length
                header = line[1:].split(None, 1)
                if not header or not header[0] or header[0] in lengths:
                    raise PlanError(f"empty or duplicate FASTA identifier at {path}:{line_number}")
                current = header[0]
                length = 0
            else:
                if current is None:
                    raise PlanError(f"sequence before first FASTA header at {path}:{line_number}")
                if not line or not re.fullmatch(r"[A-Za-z*.-]+", line):
                    raise PlanError(f"invalid FASTA sequence at {path}:{line_number}")
                length += len(line)
    if current is not None:
        if current in lengths:
            raise PlanError(f"duplicate FASTA identifier: {current!r} in {path}")
        lengths[current] = length
    if not lengths or any(value <= 0 for value in lengths.values()):
        raise PlanError(f"FASTA contains no records or an empty sequence: {path}")
    return lengths


def _canonical_read_uri(value: str, label: str) -> str:
    parsed = urllib.parse.urlsplit(value)
    if parsed.scheme == "file":
        if parsed.netloc not in {"", "localhost"} or parsed.query or parsed.fragment:
            raise PlanError(f"{label} file URI must have no remote authority, query, or fragment")
        decoded = urllib.parse.unquote(parsed.path, errors="strict")
        if not Path(decoded).is_absolute():
            raise PlanError(f"{label} file URI must contain an absolute path")
        # Resolve lexical aliases and symlinks when present.  The worker decodes
        # the resulting canonical URI; caller labels never participate in local
        # physical-object identity.
        return Path(decoded).resolve(strict=False).as_uri()
    if parsed.scheme == "s3" and parsed.netloc and parsed.path.lstrip("/"):
        try:
            port = parsed.port
        except ValueError as error:
            raise PlanError(f"{label} has an invalid S3 authority") from error
        if (
            parsed.query
            or parsed.fragment
            or parsed.username is not None
            or parsed.password is not None
            or port is not None
        ):
            raise PlanError(f"{label} S3 URI must not contain userinfo, a port, query, or fragment")
        bucket = parsed.netloc.lower()
        try:
            key = urllib.parse.unquote(parsed.path.lstrip("/"), errors="strict")
        except UnicodeDecodeError as error:
            raise PlanError(f"{label} contains an invalid percent-encoded key") from error
        if not key:
            raise PlanError(f"{label} S3 URI has an empty decoded key")
        return f"s3://{bucket}/{urllib.parse.quote(key, safe='/-_.~')}"
    raise PlanError(f"{label} must use a valid s3:// or file:// URI")


def load_manifest(
    manifest: os.PathLike[str] | str,
    expected_manifest_sha256: str,
) -> tuple[Path, str, list[dict[str, str]]]:
    """Read a strict TSV manifest after checking its caller-pinned digest."""
    path = _absolute_file(str(Path(manifest).resolve()), "manifest")
    expected = _require_sha(expected_manifest_sha256, "expected manifest SHA-256")
    observed = sha256_file(path)
    if observed != expected:
        raise PlanError(f"manifest SHA-256 mismatch: expected {expected}, observed {observed}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", restkey="__extra__", restval=None)
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            raise PlanError(
                "manifest header must exactly equal: " + "\t".join(MANIFEST_FIELDS)
            )
        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, 2):
            if row.get("__extra__") is not None:
                raise PlanError(f"manifest line {line_number} has extra columns")
            if any(row[field] is None or row[field] == "" for field in MANIFEST_FIELDS):
                raise PlanError(f"manifest line {line_number} has an empty field")
            if any("\r" in row[field] or "\n" in row[field] for field in MANIFEST_FIELDS):
                raise PlanError(f"manifest line {line_number} contains a line break in a field")
            normalized = {field: str(row[field]) for field in MANIFEST_FIELDS}
            normalized["__line__"] = str(line_number)
            rows.append(normalized)
    if not rows:
        raise PlanError("manifest has no data rows")
    return path, observed, rows


def _validated_groups(rows: Sequence[Mapping[str, str]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[Mapping[str, str]]] = {}
    seen_read_identity: set[tuple[str, str, str, str]] = set()
    seen_local_uris: set[tuple[str, str, str]] = set()
    seen_s3_objects: set[tuple[str, str, str, str, str]] = set()
    seen_content_sha256: set[tuple[str, str, str]] = set()
    verified_files: dict[tuple[str, str], Path] = {}
    fai_cache: dict[Path, dict[str, int]] = {}
    fasta_length_cache: dict[Path, dict[str, int]] = {}

    def verify_cached(path: str, expected_sha256: str, label: str) -> Path:
        key = (str(Path(path).expanduser()), expected_sha256.lower())
        if key not in verified_files:
            verified_files[key] = _verify_file(path, expected_sha256, label)
        return verified_files[key]

    for row in rows:
        line = row.get("__line__", "?")
        sample = _require_id(row["sample_id"], f"sample_id at line {line}")
        locus = _require_id(row["locus_id"], f"locus_id at line {line}")
        part = _positive_int(row["part_id"], f"part_id at line {line}")
        _positive_int(row["expected_part_count"], f"expected_part_count at line {line}")
        _positive_int(row["read_size"], f"read_size at line {line}")
        _require_sha(row["read_sha256"], f"read_sha256 at line {line}")
        if row["read_format"] != "bam":
            raise PlanError(f"read_format at line {line} must be exactly 'bam'")
        if row["minimap2_preset"] not in {"map-ont", "map-hifi"}:
            raise PlanError(f"unsupported minimap2_preset at line {line}")
        canonical_read_uri = _canonical_read_uri(row["read_uri"], f"read_uri at line {line}")
        identity = (sample, locus, row["part_id"], row["read_uri"])
        if identity in seen_read_identity:
            raise PlanError(f"duplicate part row at line {line}")
        seen_read_identity.add(identity)
        content_identity = (sample, locus, row["read_sha256"].lower())
        if content_identity in seen_content_sha256:
            raise PlanError(
                f"target {sample}/{locus} repeats read content SHA-256 at line {line}; "
                "a label or alternate path cannot make identical content a new part"
            )
        seen_content_sha256.add(content_identity)
        parsed_read = urllib.parse.urlsplit(canonical_read_uri)
        if parsed_read.scheme == "file":
            local_identity = (sample, locus, canonical_read_uri)
            if local_identity in seen_local_uris:
                raise PlanError(
                    f"target {sample}/{locus} repeats one canonical local read URI at line {line}; "
                    "read_version is not local object identity"
                )
            seen_local_uris.add(local_identity)
        else:
            s3_identity = (
                sample,
                locus,
                parsed_read.netloc,
                urllib.parse.unquote(parsed_read.path.lstrip("/"), errors="strict"),
                row["read_version"],
            )
            if s3_identity in seen_s3_objects:
                raise PlanError(
                    f"target {sample}/{locus} repeats one immutable S3 bucket/key/VersionId "
                    f"at line {line}; labels and ETags cannot make it a new part"
                )
            seen_s3_objects.add(s3_identity)
        grouped.setdefault((sample, locus), []).append(row)

    output: list[dict[str, Any]] = []
    for (sample, locus), members in sorted(grouped.items()):
        first = members[0]
        for field in TARGET_FIELDS:
            values = {row[field] for row in members}
            if len(values) != 1:
                raise PlanError(f"target {sample}/{locus} has inconsistent {field}")
        expected_count = _positive_int(first["expected_part_count"], "expected_part_count")
        part_ids = [_positive_int(row["part_id"], "part_id") for row in members]
        expected_ids = list(range(1, expected_count + 1))
        if sorted(part_ids) != expected_ids or len(part_ids) != len(set(part_ids)):
            raise PlanError(
                f"target {sample}/{locus} parts must be exactly 1..{expected_count}; "
                f"observed {sorted(part_ids)}"
            )

        assemblies: list[dict[str, Any]] = []
        for index in (1, 2):
            prefix = f"assembly{index}"
            assembly_id = _require_id(first[f"{prefix}_id"], f"{prefix}_id")
            fasta = verify_cached(
                first[f"{prefix}_fasta"], first[f"{prefix}_fasta_sha256"], f"{prefix}_fasta"
            )
            fai = verify_cached(
                first[f"{prefix}_fai"], first[f"{prefix}_fai_sha256"], f"{prefix}_fai"
            )
            contig = first[f"{prefix}_contig"]
            if not contig or any(char.isspace() for char in contig):
                raise PlanError(f"{prefix}_contig must be a non-whitespace FASTA identifier")
            if fai not in fai_cache:
                fai_cache[fai] = _read_fai(fai)
            if fasta not in fasta_length_cache:
                fasta_length_cache[fasta] = _fasta_lengths(fasta)
            fai_records = fai_cache[fai]
            fasta_lengths = fasta_length_cache[fasta]
            if contig not in fai_records:
                raise PlanError(f"{prefix}_contig {contig!r} is absent from {fai}")
            if contig not in fasta_lengths:
                raise PlanError(f"{prefix}_contig {contig!r} is absent from {fasta}")
            if fasta_lengths[contig] != fai_records[contig]:
                raise PlanError(
                    f"{prefix}_contig length differs between FASTA ({fasta_lengths[contig]}) "
                    f"and FAI ({fai_records[contig]})"
                )
            start = _positive_int(first[f"{prefix}_body_start"], f"{prefix}_body_start", allow_zero=True)
            end = _positive_int(first[f"{prefix}_body_end"], f"{prefix}_body_end")
            if not start < end <= fai_records[contig]:
                raise PlanError(
                    f"{prefix} body interval [{start}, {end}) exceeds contig length {fai_records[contig]}"
                )
            assemblies.append(
                {
                    "assembly_id": assembly_id,
                    "fasta": str(fasta),
                    "fai": str(fai),
                    "fasta_sha256": first[f"{prefix}_fasta_sha256"].lower(),
                    "fai_sha256": first[f"{prefix}_fai_sha256"].lower(),
                    "contig": contig,
                    "combined_contig": f"{assembly_id}__{contig}",
                    "contig_length": fai_records[contig],
                    "body_start": start,
                    "body_end": end,
                }
            )
        if assemblies[0]["assembly_id"] == assemblies[1]["assembly_id"]:
            raise PlanError(f"target {sample}/{locus} must declare two distinct assembly IDs")
        if assemblies[0]["combined_contig"] == assemblies[1]["combined_contig"]:
            raise PlanError(f"target {sample}/{locus} produces duplicate combined contig IDs")
        physical_assembly_targets = {
            (assembly["fasta_sha256"], assembly["contig"]) for assembly in assemblies
        }
        if len(physical_assembly_targets) != 2:
            raise PlanError(
                f"target {sample}/{locus} repeats one biological assembly target under "
                "different assembly IDs; two physical assembly targets are required"
            )

        flank = _positive_int(first["full_flank_bp"], "full_flank_bp")
        for assembly in assemblies:
            body_start = int(assembly["body_start"])
            body_end = int(assembly["body_end"])
            contig_length = int(assembly["contig_length"])
            flank_records = [
                {
                    "side": "left",
                    "start": max(0, body_start - flank),
                    "end": body_start,
                    "requested_bp": flank,
                    "available_bp": min(flank, body_start),
                    "callable_bp": min(flank, body_start),
                    "boundary_truncated": body_start < flank,
                },
                {
                    "side": "right",
                    "start": body_end,
                    "end": min(contig_length, body_end + flank),
                    "requested_bp": flank,
                    "available_bp": min(flank, contig_length - body_end),
                    "callable_bp": min(flank, contig_length - body_end),
                    "boundary_truncated": contig_length - body_end < flank,
                },
            ]
            nonempty = [record for record in flank_records if record["available_bp"] > 0]
            if not nonempty:
                raise PlanError(
                    f"assembly {assembly['assembly_id']} has no callable outer flank "
                    "on either side after FAI boundary clamping"
                )
            if len(nonempty) == 1:
                completeness = "ONE_SIDED_BOUNDARY_FLANK"
            elif any(record["boundary_truncated"] for record in flank_records):
                completeness = "TWO_SIDED_BOUNDARY_TRUNCATED_FLANK"
            else:
                completeness = "TWO_SIDED_FULL_FLANK"
            assembly["flanks"] = flank_records
            assembly["flank_completeness"] = completeness
            assembly["boundary_truncated_sides"] = [
                record["side"]
                for record in flank_records
                if record["boundary_truncated"]
            ]
        bait_min = _positive_int(first["bait_min_aligned_bp"], "bait_min_aligned_bp")
        kcon = verify_cached(first["kcon_fasta"], first["kcon_sha256"], "kcon_fasta")
        parts = []
        for row in sorted(members, key=lambda member: int(member["part_id"])):
            parts.append(
                {
                    "part_id": int(row["part_id"]),
                    "read_uri": _canonical_read_uri(row["read_uri"], "read_uri"),
                    "read_version": row["read_version"],
                    "read_etag": row["read_etag"],
                    "read_size": int(row["read_size"]),
                    "read_sha256": row["read_sha256"].lower(),
                    "read_format": row["read_format"],
                    "minimap2_preset": row["minimap2_preset"],
                }
            )
        output.append(
            {
                "sample_id": sample,
                "locus_id": locus,
                "expected_part_count": expected_count,
                "assemblies": assemblies,
                "full_flank_bp": flank,
                "bait_min_aligned_bp": bait_min,
                "kcon_fasta": str(kcon),
                "kcon_sha256": first["kcon_sha256"].lower(),
                "parts": parts,
            }
        )
    return output


def _bed_texts(target: Mapping[str, Any]) -> dict[str, str]:
    core: list[str] = []
    full: list[str] = []
    outer: list[str] = []
    for assembly in target["assemblies"]:
        contig = assembly["combined_contig"]
        body_start = int(assembly["body_start"])
        body_end = int(assembly["body_end"])
        flanks = assembly["flanks"]
        nonempty = [record for record in flanks if int(record["available_bp"]) > 0]
        if not nonempty:
            raise PlanError(f"assembly {assembly['assembly_id']} has no callable outer flank")
        window_start = min([body_start, *[int(record["start"]) for record in nonempty]])
        window_end = max([body_end, *[int(record["end"]) for record in nonempty]])
        core.append(f"{contig}\t{body_start}\t{body_end}\t{assembly['assembly_id']}_body\n")
        full.append(f"{contig}\t{window_start}\t{window_end}\t{assembly['assembly_id']}_full_window\n")
        for record in nonempty:
            outer.append(
                f"{contig}\t{record['start']}\t{record['end']}\t"
                f"{assembly['assembly_id']}_{record['side']}_outer_flank\n"
            )
    return {"core_body_union": "".join(core), "full_window": "".join(full), "outer_flanks": "".join(outer)}


def _default_code_paths() -> list[Path]:
    source = Path(__file__).resolve().parent
    return [source / name for name in DEFAULT_CODE_NAMES]


def _code_fingerprints(code_paths: Sequence[os.PathLike[str] | str] | None) -> dict[str, dict[str, Any]]:
    paths = [Path(path).resolve() for path in (code_paths or _default_code_paths())]
    names = [path.name for path in paths]
    if len(names) != len(set(names)):
        raise PlanError("code paths must have unique basenames")
    return {path.name: fingerprint(path) for path in sorted(paths, key=lambda item: item.name)}


def _write_exact(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if not path.is_file() or path.read_bytes() != content:
            raise PlanError(f"refusing to replace non-identical immutable artifact: {path}")
        return
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("xb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _plan_digest(plan: Mapping[str, Any]) -> str:
    material = dict(plan)
    material.pop("plan_digest", None)
    return digest_value(material)


def render_plans(
    manifest: os.PathLike[str] | str,
    expected_manifest_sha256: str,
    run_root: os.PathLike[str] | str,
    *,
    code_paths: Sequence[os.PathLike[str] | str] | None = None,
) -> list[Path]:
    """Validate and render one immutable target plan per sample/locus."""
    manifest_path, manifest_sha, rows = load_manifest(manifest, expected_manifest_sha256)
    groups = _validated_groups(rows)
    codes = _code_fingerprints(code_paths)
    if "cnv_plan.py" not in codes:
        raise PlanError("code paths must include cnv_plan.py for claims, publication, and fresh resume validation")
    semantic_rows = [{key: value for key, value in row.items() if key != "__line__"} for row in rows]
    semantic_sha = digest_value(semantic_rows)
    run_material = {
        "schema_version": SCHEMA_VERSION,
        "manifest_sha256": manifest_sha,
        "manifest_semantic_sha256": semantic_sha,
        "code_sha256": {name: value["sha256"] for name, value in codes.items()},
    }
    run_id = "run_" + digest_value(run_material)[:24]
    root = Path(run_root).expanduser().resolve()
    if root.exists() and not root.is_dir():
        raise PlanError(f"run_root is not a directory: {root}")
    run_dir = root / run_id
    plan_paths: list[Path] = []

    for target in groups:
        bed_texts = _bed_texts(target)
        target_material = {
            "run_id": run_id,
            "sample_id": target["sample_id"],
            "locus_id": target["locus_id"],
            "expected_part_count": target["expected_part_count"],
            "assemblies": target["assemblies"],
            "full_flank_bp": target["full_flank_bp"],
            "bait_min_aligned_bp": target["bait_min_aligned_bp"],
            "kcon_sha256": target["kcon_sha256"],
            "parts": target["parts"],
            "beds_sha256": {name: hashlib.sha256(text.encode()).hexdigest() for name, text in bed_texts.items()},
        }
        target_digest = digest_value(target_material)
        target_id = f"{target['sample_id']}__{target['locus_id']}__{target_digest[:16]}"
        target_dir = run_dir / "targets" / target_id
        bed_dir = target_dir / "beds"
        reference_dir = target_dir / "reference"
        part_dir = target_dir / "parts"
        analysis_dir = target_dir / "analysis"
        prefix = target_id
        bed_paths = {
            # The enclosing immutable target directory supplies the canonical
            # identity; stable basenames make plan review and tooling simpler.
            "core_body_union_bed": bed_dir / "core_body_union.bed",
            "full_window_bed": bed_dir / "full_window.bed",
            "outer_flanks_bed": bed_dir / "outer_flanks.bed",
        }
        _write_exact(bed_paths["core_body_union_bed"], bed_texts["core_body_union"].encode())
        _write_exact(bed_paths["full_window_bed"], bed_texts["full_window"].encode())
        _write_exact(bed_paths["outer_flanks_bed"], bed_texts["outer_flanks"].encode())

        parts = []
        for part in target["parts"]:
            part_prefix = f"{prefix}.part_{part['part_id']:06d}"
            part_output = part_dir / f"part_{part['part_id']:06d}"
            parts.append(
                {
                    **part,
                    "bam": str(part_output / f"{part_prefix}.full_window.bam"),
                    "bai": str(part_output / f"{part_prefix}.full_window.bam.bai"),
                    "primary_qnames": str(part_output / f"{part_prefix}.primary_qnames.txt"),
                    "claim": str(part_output / f"{part_prefix}.CLAIM.json"),
                    "receipt": str(part_output / f"{part_prefix}.DONE.json"),
                }
            )
        plotter = codes.get("plot_cnv_depth_v2.py")
        if plotter is None:
            raise PlanError("code paths must include plot_cnv_depth_v2.py")
        plan: dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "workflow_status": "provisional",
            "production_authorized": False,
            "run_id": run_id,
            "run_dir": str(run_dir),
            "target_id": target_id,
            "target_dir": str(target_dir),
            "target_digest": target_digest,
            "manifest": {"path": str(manifest_path), "sha256": manifest_sha, "semantic_sha256": semantic_sha},
            "code": codes,
            "plan_digest": "",
            "prepared_receipt": str(reference_dir / f"{prefix}.PREPARED.json"),
            "analysis_claim": str(analysis_dir / f"{prefix}.ANALYSIS_CLAIM.json"),
            "bait_min_aligned_bp": target["bait_min_aligned_bp"],
            "target": {
                "sample_id": target["sample_id"],
                "locus_id": target["locus_id"],
                "expected_part_count": target["expected_part_count"],
                "full_flank_bp": target["full_flank_bp"],
                "assemblies": target["assemblies"],
            },
            "reference": {
                "combined_fasta": str(reference_dir / f"{prefix}.combined.fa"),
                "combined_fai": str(reference_dir / f"{prefix}.combined.fa.fai"),
                "combined_mmi": str(reference_dir / f"{prefix}.combined.mmi"),
                "bait_fasta": str(reference_dir / f"{prefix}.full_window_bait.fa"),
                **{key: str(value) for key, value in bed_paths.items()},
                "core_body_union_bed_sha256": hashlib.sha256(
                    bed_texts["core_body_union"].encode()
                ).hexdigest(),
                "full_window_bed_sha256": hashlib.sha256(
                    bed_texts["full_window"].encode()
                ).hexdigest(),
                "outer_flanks_bed_sha256": hashlib.sha256(
                    bed_texts["outer_flanks"].encode()
                ).hexdigest(),
                "kcon_fasta": target["kcon_fasta"],
                "kcon_sha256": target["kcon_sha256"],
            },
            "plotter": {"path": plotter["path"], "sha256": plotter["sha256"]},
            "parts": parts,
            "outputs": {
                "merged_bam": str(analysis_dir / f"{prefix}.merged.full_window.bam"),
                "merged_bai": str(analysis_dir / f"{prefix}.merged.full_window.bam.bai"),
                "raw_depth_tsv": str(analysis_dir / f"{prefix}.raw_depth.tsv"),
                "mapq10_depth_tsv": str(analysis_dir / f"{prefix}.mapq10_depth.tsv"),
                "summary_json": str(analysis_dir / f"{prefix}.coverage_summary.json"),
                "raw_pdf": str(analysis_dir / f"{prefix}.raw_depth.pdf"),
                "mapq10_pdf": str(analysis_dir / f"{prefix}.mapq10_depth.pdf"),
                "kcon_paf": str(analysis_dir / f"{prefix}.assembly_resolved_kcon.paf"),
                "analysis_receipt": str(analysis_dir / f"{prefix}.ANALYSIS_DONE.json"),
            },
        }
        plan["plan_digest"] = _plan_digest(plan)
        plan_path = target_dir / f"{prefix}.plan.json"
        _write_exact(plan_path, json.dumps(plan, indent=2, sort_keys=True).encode() + b"\n")
        plan_paths.append(plan_path)

    run_manifest = {
        "schema_version": "cnv_v2_run_1",
        "run_id": run_id,
        "manifest": {"path": str(manifest_path), "sha256": manifest_sha, "semantic_sha256": semantic_sha},
        "code": codes,
        "plans": [str(path) for path in plan_paths],
    }
    _write_exact(run_dir / f"{run_id}.json", json.dumps(run_manifest, indent=2, sort_keys=True).encode() + b"\n")
    return plan_paths


def load_plan(path: os.PathLike[str] | str, *, verify_files: bool = True) -> dict[str, Any]:
    plan_path = Path(path).expanduser().resolve()
    if not plan_path.is_file():
        raise PlanError(f"plan does not exist: {plan_path}")
    try:
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PlanError(f"cannot read plan {plan_path}: {error}") from error
    if not isinstance(plan, dict) or plan.get("schema_version") != SCHEMA_VERSION:
        raise PlanError(f"unsupported plan schema: {plan_path}")
    expected = plan.get("plan_digest")
    if not isinstance(expected, str) or not SHA256.fullmatch(expected) or _plan_digest(plan) != expected:
        raise PlanError(f"plan digest mismatch: {plan_path}")
    if not SAFE_ID.fullmatch(str(plan.get("target_id", ""))):
        raise PlanError(f"unsafe target_id in plan: {plan_path}")
    parts = plan.get("parts")
    expected_count = plan.get("target", {}).get("expected_part_count")
    if not isinstance(parts, list) or [part.get("part_id") for part in parts] != list(range(1, expected_count + 1)):
        raise PlanError(f"plan does not declare exact ordered parts: {plan_path}")
    destination_paths: list[str] = [
        str(plan.get("prepared_receipt", "")),
        str(plan.get("analysis_claim", "")),
    ]
    reference = plan.get("reference")
    if not isinstance(reference, dict):
        raise PlanError(f"plan reference output map is missing: {plan_path}")
    destination_paths.extend(
        str(reference.get(name, ""))
        for name in (
            "combined_fasta",
            "combined_fai",
            "combined_mmi",
            "bait_fasta",
            "core_body_union_bed",
            "full_window_bed",
            "outer_flanks_bed",
        )
    )
    for part in parts:
        destination_paths.extend(
            str(part.get(name, ""))
            for name in ("bam", "bai", "primary_qnames", "claim", "receipt")
        )
    outputs = plan.get("outputs")
    if not isinstance(outputs, dict):
        raise PlanError(f"plan analysis output map is missing: {plan_path}")
    destination_paths.extend(str(value) for value in outputs.values())
    if any(not value or not Path(value).is_absolute() for value in destination_paths):
        raise PlanError(f"plan contains a missing or non-absolute destination path: {plan_path}")
    resolved_parent_paths = [
        _validate_plan_bound_destination(plan, value, "plan-bound destination path")
        for value in destination_paths
    ]
    if len({str(value) for value in resolved_parent_paths}) != len(resolved_parent_paths):
        raise PlanError(f"plan contains duplicate destination paths: {plan_path}")
    if verify_files:
        manifest = plan["manifest"]
        _verify_file(manifest["path"], manifest["sha256"], "plan manifest")
        for name, item in plan["code"].items():
            _verify_file(item["path"], item["sha256"], f"plan code {name}")
        _verify_file(plan["plotter"]["path"], plan["plotter"]["sha256"], "plan plotter")
        _verify_file(plan["reference"]["kcon_fasta"], plan["reference"]["kcon_sha256"], "plan KCON")
        for assembly in plan["target"]["assemblies"]:
            _verify_file(assembly["fasta"], assembly["fasta_sha256"], "plan assembly FASTA")
            _verify_file(assembly["fai"], assembly["fai_sha256"], "plan assembly FAI")
        for bed_name in ("core_body_union_bed", "full_window_bed", "outer_flanks_bed"):
            bed = fingerprint(plan["reference"][bed_name])
            expected_bed_sha = plan["reference"].get(f"{bed_name}_sha256")
            if bed["sha256"] != expected_bed_sha:
                raise PlanError(f"plan BED SHA-256 mismatch: {bed_name}")
    return plan


def write_receipt(path: os.PathLike[str] | str, payload: Mapping[str, Any]) -> Path:
    """Write an immutable JSON receipt after its payload is fully assembled."""
    receipt_path = Path(path).expanduser().resolve()
    document = dict(payload)
    document["receipt_digest"] = digest_value(document)
    _write_exact(receipt_path, json.dumps(document, indent=2, sort_keys=True).encode() + b"\n")
    return receipt_path


def acquire_create_only_claim(
    path: os.PathLike[str] | str,
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    """Atomically create one persistent, self-digested claim; never steal it."""
    claim_path = _literal_absolute_path(path, "claim path")
    claim_path.parent.mkdir(parents=True, exist_ok=True)
    if os.path.lexists(claim_path):
        raise PlanError(f"claim already exists and will not be stolen: {claim_path}")
    document = dict(payload)
    document["schema_version"] = CLAIM_SCHEMA
    document["claim_digest"] = digest_value(document)
    encoded = json.dumps(document, indent=2, sort_keys=True).encode() + b"\n"
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    try:
        descriptor = os.open(claim_path, flags, 0o444)
    except FileExistsError as error:
        raise PlanError(f"claim already exists and will not be stolen: {claim_path}") from error
    try:
        with os.fdopen(descriptor, "wb", closefd=True) as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        # A short or otherwise damaged claim is deliberately left visible.  A
        # later invocation must fail closed rather than silently deleting it.
        raise
    return document


def load_claim(path: os.PathLike[str] | str) -> dict[str, Any]:
    claim_path = _literal_absolute_path(path, "claim path")
    if claim_path.is_symlink():
        raise PlanError(f"claim path must be a literal regular file, not a symlink: {claim_path}")
    try:
        claim = json.loads(claim_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PlanError(f"cannot read claim {claim_path}: {error}") from error
    if not isinstance(claim, dict) or claim.get("schema_version") != CLAIM_SCHEMA:
        raise PlanError(f"claim schema mismatch: {claim_path}")
    claimed = claim.get("claim_digest")
    material = dict(claim)
    material.pop("claim_digest", None)
    if not isinstance(claimed, str) or digest_value(material) != claimed:
        raise PlanError(f"claim digest mismatch: {claim_path}")
    return claim


def _claim_payload(
    plan: Mapping[str, Any],
    kind: str,
    claimant_token: str,
    *,
    part_id: int | None = None,
) -> tuple[Path, dict[str, Any]]:
    if not claimant_token or any(char in claimant_token for char in "\r\n\x00"):
        raise PlanError("claimant token must be nonempty and single-line")
    common: dict[str, Any] = {
        "status": "CLAIMED",
        "kind": kind,
        "run_id": plan["run_id"],
        "target_id": plan["target_id"],
        "target_digest": plan["target_digest"],
        "plan_digest": plan["plan_digest"],
        "manifest_sha256": plan["manifest"]["sha256"],
        "claimant_token": claimant_token,
    }
    if kind == "part":
        if part_id is None:
            raise PlanError("part claim requires part_id")
        part = get_part(plan, part_id)
        common.update(
            {
                "part_id": part_id,
                "read_uri": part["read_uri"],
                "read_version": part["read_version"],
                "read_etag": part["read_etag"],
                "read_size": part["read_size"],
                "read_sha256": part["read_sha256"],
            }
        )
        path = Path(part["claim"])
        _validate_plan_bound_destination(plan, path, "part claim path")
        return path, common
    if kind == "analysis":
        if part_id is not None:
            raise PlanError("analysis claim must not declare part_id")
        path = Path(plan["analysis_claim"])
        _validate_plan_bound_destination(plan, path, "analysis claim path")
        return path, common
    raise PlanError(f"unsupported claim kind: {kind}")


def acquire_plan_claim(
    plan_or_path: Mapping[str, Any] | os.PathLike[str] | str,
    kind: str,
    claimant_token: str,
    *,
    part_id: int | None = None,
) -> dict[str, Any]:
    plan = dict(plan_or_path) if isinstance(plan_or_path, Mapping) else load_plan(plan_or_path)
    path, payload = _claim_payload(plan, kind, claimant_token, part_id=part_id)
    # Recheck immediately before any mkdir/O_EXCL action.  A parent directory
    # replaced by a symlink after plan rendering must not redirect the claim.
    _validate_plan_bound_destination(plan, path, f"{kind} claim path")
    return acquire_create_only_claim(path, payload)


def validate_plan_claim(
    plan: Mapping[str, Any],
    kind: str,
    expected_digest: str,
    *,
    part_id: int | None = None,
) -> dict[str, Any]:
    path, expected = _claim_payload(plan, kind, "placeholder", part_id=part_id)
    claim = load_claim(path)
    expected.pop("claimant_token")
    for key, value in expected.items():
        if claim.get(key) != value:
            raise PlanError(f"{kind} claim binding mismatch for {key}")
    if claim.get("claim_digest") != expected_digest:
        raise PlanError(f"{kind} claim digest mismatch")
    token = claim.get("claimant_token")
    if not isinstance(token, str) or not token or any(char in token for char in "\r\n\x00"):
        raise PlanError(f"{kind} claim has an invalid claimant token")
    return claim


def publish_create_only(
    staged: os.PathLike[str] | str,
    destination: os.PathLike[str] | str,
) -> Path:
    """Atomically publish a same-filesystem file without replacement."""
    source = Path(staged).expanduser().resolve()
    target = _literal_absolute_path(destination, "publication destination")
    if not source.is_file():
        raise PlanError(f"staged publication file is missing: {source}")
    if not target.parent.is_dir():
        raise PlanError(f"publication directory is missing: {target.parent}")
    if os.path.lexists(target):
        raise PlanError(f"canonical output already exists and will not be overwritten: {target}")
    if source.stat().st_dev != target.parent.stat().st_dev:
        raise PlanError("staged and canonical publication paths are on different filesystems")
    try:
        os.link(source, target)
    except FileExistsError as error:
        raise PlanError(f"canonical output already exists and will not be overwritten: {target}") from error
    except OSError as error:
        if error.errno == errno.EXDEV:
            raise PlanError("staged and canonical publication paths are on different filesystems") from error
        raise PlanError(f"create-only publication failed for {target}: {error}") from error
    source.unlink()
    return target


def _iter_canonical_qnames(path: os.PathLike[str] | str) -> Iterable[bytes]:
    inventory = Path(path)
    if not inventory.is_file():
        raise PlanError(f"primary-QNAME inventory is missing: {inventory}")
    previous: bytes | None = None
    with inventory.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.endswith(b"\n"):
                raise PlanError(f"primary-QNAME inventory lacks terminal LF at line {line_number}: {inventory}")
            qname = raw[:-1]
            if not qname or any(byte <= 32 or byte >= 127 for byte in qname):
                raise PlanError(f"primary-QNAME inventory has an invalid QNAME at line {line_number}: {inventory}")
            if previous is not None and qname <= previous:
                reason = "duplicate" if qname == previous else "unsorted"
                raise PlanError(f"primary-QNAME inventory is {reason} at line {line_number}: {inventory}")
            previous = qname
            yield qname


def inspect_qname_inventory(path: os.PathLike[str] | str) -> dict[str, Any]:
    count = sum(1 for _ in _iter_canonical_qnames(path))
    return {
        **fingerprint_allow_empty(path),
        "qname_count": count,
        "encoding": QNAME_ENCODING,
        "order": QNAME_ORDER,
    }


def verify_qname_inventory_union(paths: Sequence[os.PathLike[str] | str]) -> dict[str, Any]:
    """Stream sorted inventories and prove a pairwise-disjoint exact union."""
    if not paths:
        raise PlanError("at least one primary-QNAME inventory is required")
    inventories = [inspect_qname_inventory(path) for path in paths]
    iterators = [iter(_iter_canonical_qnames(record["path"])) for record in inventories]
    heap: list[tuple[bytes, int]] = []
    for index, iterator in enumerate(iterators):
        try:
            heapq.heappush(heap, (next(iterator), index))
        except StopIteration:
            pass
    union_count = 0
    previous: bytes | None = None
    previous_source: int | None = None
    try:
        while heap:
            qname, source = heapq.heappop(heap)
            if previous == qname:
                raise PlanError(
                    "primary QNAME occurs in multiple parts: "
                    f"{qname.decode('ascii')} (inventory indexes {previous_source} and {source})"
                )
            union_count += 1
            previous = qname
            previous_source = source
            try:
                heapq.heappush(heap, (next(iterators[source]), source))
            except StopIteration:
                pass
    finally:
        for iterator in iterators:
            close = getattr(iterator, "close", None)
            if close is not None:
                close()
    for initial in inventories:
        if inspect_qname_inventory(initial["path"]) != initial:
            raise PlanError("primary-QNAME inventory changed during union verification")
    return {
        "schema_version": QNAME_UNION_SCHEMA,
        "encoding": QNAME_ENCODING,
        "order": QNAME_ORDER,
        "inventories": inventories,
        "verified_union_qname_count": union_count,
    }


def _derive_primary_qname_inventory(
    bam: os.PathLike[str] | str,
    output: os.PathLike[str] | str,
    *,
    samtools: str = "samtools",
) -> Path:
    """Stream primary BAM QNAMEs through a bytewise external unique sort."""
    bam_path = Path(bam).resolve()
    output_path = Path(output).resolve()
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    try:
        viewer = subprocess.Popen(
            [samtools, "view", "-F", "0x900", str(bam_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        with output_path.open("xb") as destination:
            sorter = subprocess.Popen(
                ["sort", "-u"],
                stdin=subprocess.PIPE,
                stdout=destination,
                stderr=subprocess.PIPE,
                env=environment,
            )
            assert viewer.stdout is not None
            assert sorter.stdin is not None
            try:
                for raw in viewer.stdout:
                    qname = raw.split(b"\t", 1)[0]
                    if not qname or any(byte <= 32 or byte >= 127 for byte in qname):
                        raise PlanError("final BAM emitted an invalid primary QNAME")
                    sorter.stdin.write(qname + b"\n")
            finally:
                sorter.stdin.close()
                viewer.stdout.close()
            sorter_stderr = sorter.stderr.read() if sorter.stderr is not None else b""
            if sorter.stderr is not None:
                sorter.stderr.close()
            sorter_status = sorter.wait()
        viewer_stderr = viewer.stderr.read() if viewer.stderr is not None else b""
        if viewer.stderr is not None:
            viewer.stderr.close()
        viewer_status = viewer.wait()
    except FileNotFoundError as error:
        raise PlanError(f"required executable not found while deriving QNAME inventory: {error.filename}") from error
    if viewer_status != 0:
        raise PlanError(f"samtools primary-QNAME derivation failed: {viewer_stderr.decode(errors='replace').strip()}")
    if sorter_status != 0:
        raise PlanError(f"bytewise primary-QNAME sort failed: {sorter_stderr.decode(errors='replace').strip()}")
    inspect_qname_inventory(output_path)
    return output_path


def verify_inventory_matches_bam(
    inventory: os.PathLike[str] | str,
    bam: os.PathLike[str] | str,
    *,
    samtools: str = "samtools",
) -> dict[str, Any]:
    expected = Path(inventory).resolve()
    record = inspect_qname_inventory(expected)
    with tempfile.TemporaryDirectory(prefix="cnv_v2.qnames.") as temporary:
        observed = Path(temporary) / "fresh.primary_qnames.txt"
        _derive_primary_qname_inventory(bam, observed, samtools=samtools)
        with expected.open("rb") as left, observed.open("rb") as right:
            while True:
                left_block = left.read(1024 * 1024)
                right_block = right.read(1024 * 1024)
                if left_block != right_block:
                    raise PlanError("primary-QNAME inventory differs from the final part BAM")
                if not left_block:
                    break
    if inspect_qname_inventory(expected) != record:
        raise PlanError("primary-QNAME inventory changed during BAM comparison")
    return record


def load_receipt(path: os.PathLike[str] | str, expected_schema: str) -> dict[str, Any]:
    receipt_path = _literal_absolute_path(path, "receipt path")
    if receipt_path.is_symlink():
        raise PlanError(f"receipt path must be a literal regular file, not a symlink: {receipt_path}")
    if not receipt_path.is_file():
        raise PlanError(f"receipt is missing: {receipt_path}")
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PlanError(f"cannot read receipt {receipt_path}: {error}") from error
    if not isinstance(receipt, dict) or receipt.get("schema_version") != expected_schema:
        raise PlanError(f"receipt schema mismatch: {receipt_path}")
    claimed = receipt.get("receipt_digest")
    material = dict(receipt)
    material.pop("receipt_digest", None)
    if not isinstance(claimed, str) or not SHA256.fullmatch(claimed) or digest_value(material) != claimed:
        raise PlanError(f"receipt digest mismatch: {receipt_path}")
    return receipt


def _assert_fingerprint(record: Mapping[str, Any], expected_path: str, label: str) -> None:
    observed = fingerprint(expected_path)
    if record != observed:
        raise PlanError(f"{label} fingerprint mismatch")


def validate_prepared(plan_or_path: Mapping[str, Any] | os.PathLike[str] | str) -> dict[str, Any]:
    plan = dict(plan_or_path) if isinstance(plan_or_path, Mapping) else load_plan(plan_or_path)
    receipt = load_receipt(plan["prepared_receipt"], PREPARED_SCHEMA)
    if not SHA256.fullmatch(str(receipt.get("receipt_digest", ""))):
        raise PlanError("prepared receipt is missing its canonical receipt digest")
    if receipt.get("run_id") != plan["run_id"] or receipt.get("target_id") != plan["target_id"]:
        raise PlanError("prepared receipt belongs to a different run or target")
    if receipt.get("status") != "PREPARED" or receipt.get("target_digest") != plan["target_digest"]:
        raise PlanError("prepared receipt status/target binding mismatch")
    if receipt.get("plan_digest") != plan["plan_digest"] or receipt.get("manifest_sha256") != plan["manifest"]["sha256"]:
        raise PlanError("prepared receipt plan/manifest binding mismatch")
    expected = {
        "combined_fasta": plan["reference"]["combined_fasta"],
        "combined_fai": plan["reference"]["combined_fai"],
        "combined_mmi": plan["reference"]["combined_mmi"],
        "bait_fasta": plan["reference"]["bait_fasta"],
        "full_window_bed": plan["reference"]["full_window_bed"],
        "core_body_union_bed": plan["reference"]["core_body_union_bed"],
        "outer_flanks_bed": plan["reference"]["outer_flanks_bed"],
    }
    artifacts = receipt.get("artifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != set(expected):
        raise PlanError("prepared receipt artifact set mismatch")
    for name, path in expected.items():
        _assert_fingerprint(artifacts[name], path, f"prepared {name}")
    return receipt


def _append_prefixed_fasta(
    source: Path,
    output: Any,
    *,
    assembly_id: str,
    target_contig: str,
    target_length: int,
    window_start: int,
    window_end: int,
    combined_seen: set[str],
) -> str:
    """Stream one FASTA into the combined namespace and return its bait slice."""
    source_seen: set[str] = set()
    current: str | None = None
    current_length = 0
    bait_chunks: list[str] = []
    target_seen = False

    def finish_record() -> None:
        nonlocal target_seen
        if current is None:
            return
        if current_length <= 0:
            raise PlanError(f"FASTA contains an empty sequence for {current!r}: {source}")
        if current == target_contig:
            target_seen = True
            if current_length != target_length:
                raise PlanError(
                    f"FASTA/FAI length mismatch for {assembly_id} target contig: "
                    f"{current_length} != {target_length}"
                )

    with source.open("r", encoding="ascii") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                finish_record()
                header = line[1:]
                split = header.split(None, 1)
                if not split or not split[0] or split[0] in source_seen:
                    raise PlanError(f"empty or duplicate FASTA identifier at {source}:{line_number}")
                current = split[0]
                source_seen.add(current)
                current_length = 0
                combined = f"{assembly_id}__{current}"
                if combined in combined_seen:
                    raise PlanError(f"duplicate combined FASTA identifier: {combined}")
                combined_seen.add(combined)
                description = split[1] if len(split) == 2 else ""
                output.write(f">{combined}{(' ' + description) if description else ''}\n")
                continue
            if current is None:
                raise PlanError(f"sequence before first FASTA header at {source}:{line_number}")
            if not line or not re.fullmatch(r"[A-Za-z*.-]+", line):
                raise PlanError(f"invalid FASTA sequence at {source}:{line_number}")
            if current == target_contig:
                overlap_start = max(window_start, current_length)
                overlap_end = min(window_end, current_length + len(line))
                if overlap_start < overlap_end:
                    bait_chunks.append(line[overlap_start - current_length : overlap_end - current_length])
            current_length += len(line)
            output.write(line + "\n")
    finish_record()
    if not source_seen or not target_seen:
        raise PlanError(f"target contig {target_contig!r} is absent from {source}")
    bait = "".join(bait_chunks)
    if len(bait) != window_end - window_start:
        raise PlanError(f"could not extract the exact bait interval from {source}")
    return bait


def _write_fasta_record(handle: Any, name: str, sequence: str, description: str = "") -> None:
    handle.write(f">{name}{(' ' + description) if description else ''}\n")
    for offset in range(0, len(sequence), 60):
        handle.write(sequence[offset : offset + 60] + "\n")


def _run_checked(command: Sequence[str], *, stdout: Any = None) -> None:
    try:
        subprocess.run(command, check=True, stdout=stdout)
    except FileNotFoundError as error:
        raise PlanError(f"required executable not found: {command[0]}") from error
    except subprocess.CalledProcessError as error:
        raise PlanError(f"command failed with status {error.returncode}: {' '.join(command)}") from error


def prepare_reference(
    plan_path: os.PathLike[str] | str,
    *,
    samtools: str = "samtools",
    minimap2: str = "minimap2",
) -> Path:
    """Build the combined reference only after an explicit prepare request."""
    plan = load_plan(plan_path)
    receipt_path = Path(plan["prepared_receipt"])
    if receipt_path.exists():
        # Resume accepts only an exact identity-bound PREPARED receipt.  Any
        # mismatch fails closed; existing partial artifacts without that
        # receipt are refused below rather than overwritten or reused.
        validate_prepared(plan)
        return receipt_path
    artifact_paths = [
        Path(plan["reference"][name])
        for name in ("combined_fasta", "combined_fai", "combined_mmi", "bait_fasta")
    ]
    if any(path.exists() for path in artifact_paths):
        raise PlanError("partial prepared artifacts exist without an exact PREPARED receipt")
    reference_dir = artifact_paths[0].parent
    reference_dir.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{plan['target_id']}.prepare.", dir=reference_dir))
    try:
        combined_temp = temporary / artifact_paths[0].name
        combined_seen: set[str] = set()
        bait_records: list[tuple[str, int, int, str]] = []
        with combined_temp.open("x", encoding="ascii", newline="\n") as output:
            for assembly in plan["target"]["assemblies"]:
                nonempty_flanks = [
                    record
                    for record in assembly["flanks"]
                    if int(record["available_bp"]) > 0
                ]
                start = min(
                    int(assembly["body_start"]),
                    *[int(record["start"]) for record in nonempty_flanks],
                )
                end = max(
                    int(assembly["body_end"]),
                    *[int(record["end"]) for record in nonempty_flanks],
                )
                bait = _append_prefixed_fasta(
                    Path(assembly["fasta"]),
                    output,
                    assembly_id=assembly["assembly_id"],
                    target_contig=assembly["contig"],
                    target_length=int(assembly["contig_length"]),
                    window_start=start,
                    window_end=end,
                    combined_seen=combined_seen,
                )
                bait_records.append((assembly["combined_contig"], start, end, bait))
        bait_temp = temporary / artifact_paths[3].name
        with bait_temp.open("x", encoding="ascii", newline="\n") as output:
            for contig, start, end, bait in bait_records:
                _write_fasta_record(output, f"{contig}:{start}-{end}", bait)
        fai_temp = temporary / artifact_paths[1].name
        mmi_temp = temporary / artifact_paths[2].name
        _run_checked([samtools, "faidx", str(combined_temp)])
        generated_fai = Path(str(combined_temp) + ".fai")
        if not generated_fai.is_file():
            raise PlanError("samtools faidx did not create the expected FAI")
        if generated_fai != fai_temp:
            generated_fai.rename(fai_temp)
        _run_checked([minimap2, "-d", str(mmi_temp), str(combined_temp)])
        for staged, final in zip((combined_temp, fai_temp, mmi_temp, bait_temp), artifact_paths):
            if not staged.is_file() or staged.stat().st_size <= 0:
                raise PlanError(f"prepare did not produce a non-empty artifact: {staged}")
            os.replace(staged, final)
        payload = {
            "schema_version": PREPARED_SCHEMA,
            "status": "PREPARED",
            "run_id": plan["run_id"],
            "target_id": plan["target_id"],
            "target_digest": plan["target_digest"],
            "plan_digest": plan["plan_digest"],
            "manifest_sha256": plan["manifest"]["sha256"],
            "artifacts": {
                name: fingerprint(plan["reference"][name])
                for name in (
                    "combined_fasta",
                    "combined_fai",
                    "combined_mmi",
                    "bait_fasta",
                    "full_window_bed",
                    "core_body_union_bed",
                    "outer_flanks_bed",
                )
            },
        }
        write_receipt(receipt_path, payload)
        validate_prepared(plan)
        return receipt_path
    finally:
        shutil.rmtree(temporary, ignore_errors=True)


def get_part(plan: Mapping[str, Any], part_id: int) -> dict[str, Any]:
    if part_id < 1:
        raise PlanError("part_id must be positive")
    matches = [part for part in plan["parts"] if part["part_id"] == part_id]
    if len(matches) != 1:
        raise PlanError(f"part_id {part_id} is not declared exactly once")
    return dict(matches[0])


def validate_part_receipt(
    plan_or_path: Mapping[str, Any] | os.PathLike[str] | str,
    part_id: int,
    *,
    samtools: str = "samtools",
) -> dict[str, Any]:
    plan = dict(plan_or_path) if isinstance(plan_or_path, Mapping) else load_plan(plan_or_path)
    prepared = validate_prepared(plan)
    part = get_part(plan, part_id)
    part_index = next(index for index, value in enumerate(plan["parts"]) if value["part_id"] == part_id)
    receipt = load_receipt(part["receipt"], PART_RECEIPT_SCHEMA)
    bindings = {
        "run_id": plan["run_id"],
        "target_id": plan["target_id"],
        "plan_digest": plan["plan_digest"],
        "plan_file_sha256": sha256_file(Path(plan_or_path).expanduser().resolve())
        if not isinstance(plan_or_path, Mapping)
        else receipt.get("plan_file_sha256"),
        "manifest_sha256": plan["manifest"]["sha256"],
        "target_digest": plan["target_digest"],
        "part_id": part_id,
        "part_index": part_index,
        "read_uri": part["read_uri"],
        "read_version": part["read_version"],
        "read_etag": part["read_etag"],
        "read_size": part["read_size"],
        "read_sha256": part["read_sha256"],
        "prepared_receipt_digest": prepared["receipt_digest"],
    }
    for key, expected in bindings.items():
        if receipt.get(key) != expected:
            raise PlanError(f"part {part_id} receipt binding mismatch for {key}")
    if receipt.get("status") != "DONE":
        raise PlanError(f"part {part_id} receipt is not DONE")
    claim_digest = receipt.get("claim_digest")
    if not isinstance(claim_digest, str) or not SHA256.fullmatch(claim_digest):
        raise PlanError(f"part {part_id} receipt has no exact claim digest")
    validate_plan_claim(plan, "part", claim_digest, part_id=part_id)
    if receipt.get("prepared_artifacts") != prepared["artifacts"]:
        raise PlanError(f"part {part_id} receipt does not bind the exact prepared artifacts")
    outputs = receipt.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != {"bam", "bai", "primary_qnames"}:
        raise PlanError(f"part {part_id} receipt artifact set mismatch")
    _assert_fingerprint(outputs["bam"], part["bam"], f"part {part_id} BAM")
    _assert_fingerprint(outputs["bai"], part["bai"], f"part {part_id} BAI")
    inventory = inspect_qname_inventory(part["primary_qnames"])
    if outputs["primary_qnames"] != inventory:
        raise PlanError(f"part {part_id} primary-QNAME inventory fingerprint mismatch")
    _run_checked([samtools, "quickcheck", "-v", part["bam"]])
    _run_checked([samtools, "idxstats", part["bam"]], stdout=subprocess.DEVNULL)
    verify_inventory_matches_bam(part["primary_qnames"], part["bam"], samtools=samtools)
    return receipt


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    render = subparsers.add_parser("render", help="validate and render deterministic target plans (offline)")
    render.add_argument("--manifest", required=True)
    render.add_argument("--expected-manifest-sha256", required=True)
    render.add_argument("--run-root", required=True)
    render.add_argument("--code-path", action="append", default=None)
    render.add_argument("--print-plan-path", action="store_true")

    prepare = subparsers.add_parser("prepare", help="explicitly build one plan's combined reference")
    prepare.add_argument("--plan", required=True)
    prepare.add_argument("--samtools", default="samtools")
    prepare.add_argument("--minimap2", default="minimap2")

    validate_plan_parser = subparsers.add_parser("validate-plan", help="validate an immutable plan")
    validate_plan_parser.add_argument("--plan", required=True)
    validate_plan_parser.add_argument("--print-json", action="store_true")

    prepared = subparsers.add_parser("validate-prepared", help="validate a PREPARED receipt")
    prepared.add_argument("--plan", required=True)

    part = subparsers.add_parser("validate-part", help="validate an exact part DONE receipt")
    part.add_argument("--plan", required=True)
    part.add_argument("--part-id", required=True, type=int)
    part.add_argument("--samtools", default="samtools")

    claim = subparsers.add_parser("claim", help="atomically acquire one persistent exact-identity claim")
    claim.add_argument("--plan", required=True)
    claim.add_argument("--kind", choices=("part", "analysis"), required=True)
    claim.add_argument("--part-id", type=int)
    claim.add_argument("--claimant-token", required=True)

    validate_claim_parser = subparsers.add_parser("validate-claim", help="validate a persistent exact-identity claim")
    validate_claim_parser.add_argument("--plan", required=True)
    validate_claim_parser.add_argument("--kind", choices=("part", "analysis"), required=True)
    validate_claim_parser.add_argument("--part-id", type=int)
    validate_claim_parser.add_argument("--expected-digest", required=True)

    publish = subparsers.add_parser("publish", help="same-filesystem create-only atomic publication")
    publish.add_argument("--staged", required=True)
    publish.add_argument("--destination", required=True)

    union = subparsers.add_parser(
        "verify-qname-union",
        help="stream canonical primary-QNAME inventories and prove pairwise disjointness",
    )
    union.add_argument("--inventory", action="append", required=True)
    union.add_argument("--output-json")

    match = subparsers.add_parser(
        "verify-inventory-bam",
        help="freshly rederive a primary-QNAME set and compare it with its inventory",
    )
    match.add_argument("--inventory", required=True)
    match.add_argument("--bam", required=True)
    match.add_argument("--samtools", default="samtools")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "render":
            paths = render_plans(
                args.manifest,
                args.expected_manifest_sha256,
                args.run_root,
                code_paths=args.code_path,
            )
            if args.print_plan_path:
                for path in paths:
                    print(path)
            else:
                print(f"validated {len(paths)} target(s); dry-run plans rendered beneath {Path(args.run_root).resolve()}")
        elif args.command == "prepare":
            print(prepare_reference(args.plan, samtools=args.samtools, minimap2=args.minimap2))
        elif args.command == "validate-plan":
            plan = load_plan(args.plan)
            print(json.dumps(plan, sort_keys=True) if args.print_json else plan["plan_digest"])
        elif args.command == "validate-prepared":
            print(validate_prepared(args.plan)["receipt_digest"])
        elif args.command == "validate-part":
            print(validate_part_receipt(args.plan, args.part_id, samtools=args.samtools)["receipt_digest"])
        elif args.command == "claim":
            if args.kind == "part" and args.part_id is None:
                raise PlanError("part claim requires --part-id")
            if args.kind == "analysis" and args.part_id is not None:
                raise PlanError("analysis claim forbids --part-id")
            claim = acquire_plan_claim(
                args.plan,
                args.kind,
                args.claimant_token,
                part_id=args.part_id,
            )
            print(claim["claim_digest"])
        elif args.command == "publish":
            print(publish_create_only(args.staged, args.destination))
        elif args.command == "validate-claim":
            plan = load_plan(args.plan)
            validate_plan_claim(
                plan,
                args.kind,
                args.expected_digest,
                part_id=args.part_id,
            )
            print(args.expected_digest)
        elif args.command == "verify-qname-union":
            union = verify_qname_inventory_union(args.inventory)
            encoded = json.dumps(union, indent=2, sort_keys=True) + "\n"
            if args.output_json:
                output = Path(args.output_json).expanduser().resolve()
                with output.open("x", encoding="utf-8") as handle:
                    handle.write(encoded)
            else:
                print(encoded, end="")
        elif args.command == "verify-inventory-bam":
            record = verify_inventory_matches_bam(
                args.inventory,
                args.bam,
                samtools=args.samtools,
            )
            print(json.dumps(record, sort_keys=True))
        else:  # pragma: no cover - argparse prevents this
            raise PlanError(f"unsupported command: {args.command}")
    except PlanError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
