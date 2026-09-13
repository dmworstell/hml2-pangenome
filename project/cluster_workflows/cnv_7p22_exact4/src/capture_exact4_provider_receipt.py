#!/usr/bin/env python3
"""Emit the exact-four provider receipt using metadata-only S3 HEAD requests.

This helper has one read-only purpose.  It validates the frozen exact-four
deployment manifest and a separately accepted expected-metadata authority,
then asks the explicitly named AWS executable for metadata for exactly the 13
versioned objects.  It never reads an object body and has no filesystem-output,
scheduler, SSH, deployment, submission, retry, or cleanup interface.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
from typing import Mapping, Sequence
from urllib.parse import urlsplit

import exact4_plan as plan


sys.dont_write_bytecode = True


EXPECTED_METADATA_SCHEMA = "hml2_7p22_exact4_expected_provider_metadata_1"
AWS_TOOL_AUTHORITY_SCHEMA = "hml2_7p22_exact4_aws_tool_authority_1"
PROVIDER_RECEIPT_SCHEMA = "hml2_7p22_exact4_provider_receipt_2"
CAPTURE_AUTHORIZATION_SHA256 = (
    "b7d2a267e44291590868e2f679beb4c7f62790f1b89778e804f105266ab04d58"
)
CAPTURE_LABEL = f"read_only_exact13_metadata_head:{CAPTURE_AUTHORIZATION_SHA256}"
EXPECTED_METADATA_AUTHORITY_SIZE_BYTES = 5801
EXPECTED_METADATA_AUTHORITY_SHA256 = (
    "5406e35fed4d44233591b978914b8befbb90d11f6c647584fa6464edca9697ff"
)
EXPECTED_METADATA_ACCEPTANCE_SIZE_BYTES = 2089
EXPECTED_METADATA_ACCEPTANCE_SHA256 = (
    "cd18055697f2433933d2ea154718237d41e70265300c024d5d27d35a24e4117b"
)
AWS_TOOL_AUTHORITY_SIZE_BYTES = 1165
AWS_TOOL_AUTHORITY_SHA256 = (
    "23728802f3249480adb3d45041a75fb38e0c6c52cadb24e4e03580e89bfeba63"
)
AWS_TOOL_ACCEPTANCE_SIZE_BYTES = 1855
AWS_TOOL_ACCEPTANCE_SHA256 = (
    "59710dc0a77ca92f52b4a531a6e974b3706a1d45649da48fba193896ab994177"
)
UNAUTHORIZED_ACTION_GATES = {
    "cancellation_authorized": False,
    "cluster_contact_authorized": False,
    "copy_number_inference_authorized": False,
    "data_access_authorized": False,
    "deployment_authorized": False,
    "download_authorized": False,
    "execution_authorized": False,
    "manuscript_use_authorized": False,
    "model_use_authorized": False,
    "network_authorized": False,
    "production_authorized": False,
    "resubmission_authorized": False,
    "result_use_authorized": False,
    "retry_authorized": False,
    "submission_authorized": False,
}
AWS_CHILD_ENVIRONMENT = {
    "AWS_CONFIG_FILE": "/dev/null",
    "AWS_EC2_METADATA_DISABLED": "true",
    "AWS_MAX_ATTEMPTS": "1",
    "AWS_PAGER": "",
    "AWS_SHARED_CREDENTIALS_FILE": "/dev/null",
    "LANG": "C",
    "LC_ALL": "C",
}
EXPECTED_SAMPLE_PARTS = (
    ("HG02027", 0),
    ("HG02027", 1),
    ("HG02027", 2),
    ("HG02178", 0),
    ("HG02178", 1),
    ("HG02178", 2),
    ("HG03669", 0),
    ("HG03669", 1),
    ("HG03669", 2),
    ("HG03669", 3),
    ("HG03834", 0),
    ("HG03834", 1),
    ("HG03834", 2),
)
EXPECTED_METADATA_AUTHORITY_PATH = "cluster_workflows/cnv_7p22_exact4/frozen/exact4_expected_provider_metadata.v1.json"
EXPECTED_METADATA_ACCEPTANCE_PATH = "handoff/7P22_EXACT_FOUR_EXPECTED_PROVIDER_METADATA_AUTHORITY_ACCEPTANCE.md"
AWS_TOOL_AUTHORITY_PATH = "cluster_workflows/cnv_7p22_exact4/frozen/exact4_aws_tool_authority.v1.json"
AWS_TOOL_ACCEPTANCE_PATH = "handoff/7P22_EXACT_FOUR_AWS_TOOL_AUTHORITY_ACCEPTANCE.md"


class ProviderReceiptCaptureError(RuntimeError):
    """A provider authority, metadata response, or command gate failed."""


def _require_lexical_normalized_absolute(
    path_text: object, label: str, *, expected: Path | None = None
) -> str:
    """Reject every alias spelling before any process-capable operation."""
    text = _require_exact_str(path_text, label)
    if "\x00" in text or not os.path.isabs(text):
        raise ProviderReceiptCaptureError(f"{label} path must be normalized absolute")
    normalized = os.path.normpath(text)
    if text != normalized or text.endswith(os.sep):
        raise ProviderReceiptCaptureError(f"{label} path must use its unique normalized absolute spelling")
    parts = Path(text).parts[1:]
    if any(part in (".", "..", "") for part in parts) or os.sep * 2 in text:
        raise ProviderReceiptCaptureError(f"{label} path contains a lexical alias")
    if expected is not None and text != str(expected):
        raise ProviderReceiptCaptureError(f"{label} path differs from the fixed authority role")
    return text


def _require_exact_str(value: object, label: str) -> str:
    if type(value) is not str or not value:
        raise ProviderReceiptCaptureError(f"{label} must be an exact nonempty string")
    return value


def _require_exact_int(value: object, label: str) -> int:
    if type(value) is not int or value <= 0:
        raise ProviderReceiptCaptureError(f"{label} must be an exact positive integer")
    return value


def _require_sha256(value: object, label: str) -> str:
    text = _require_exact_str(value, label)
    if plan.SHA256_RE.fullmatch(text) is None:
        raise ProviderReceiptCaptureError(f"{label} must be lowercase 64-hex SHA-256")
    return text


def _exact_keys(value: Mapping[str, object], expected: Sequence[str], label: str) -> None:
    observed_keys = set(value)
    expected_keys = set(expected)
    if observed_keys != expected_keys:
        raise ProviderReceiptCaptureError(
            f"{label} keys differ; missing={sorted(expected_keys-observed_keys)}, "
            f"unknown={sorted(observed_keys-expected_keys)}"
        )


def _load_canonical_object(path: str, label: str) -> dict[str, object]:
    item = Path(path)
    if not item.is_absolute():
        raise ProviderReceiptCaptureError(f"{label} path must be absolute")
    if not item.is_file() or item.is_symlink():
        raise ProviderReceiptCaptureError(f"{label} must be a regular non-symlink file")
    payload = item.read_bytes()
    value = plan.strict_json_bytes(payload, label)
    if type(value) is not dict:
        raise ProviderReceiptCaptureError(f"{label} must be an exact JSON object")
    if plan.canonical_json_bytes(value) != payload:
        raise ProviderReceiptCaptureError(f"{label} must use canonical JSON bytes")
    return value


def _stable_read_accepted_bytes(
    path_text: str,
    label: str,
    expected_size_bytes: int,
    expected_sha256: str,
) -> tuple[bytes, dict[str, object]]:
    """Read one accepted file through one descriptor before any decoding."""

    path = Path(_require_lexical_normalized_absolute(path_text, label))
    leaf_fields = (
        "st_dev", "st_ino", "st_mode", "st_nlink", "st_uid", "st_gid",
        "st_size", "st_mtime_ns", "st_ctime_ns",
    )

    def walk() -> tuple[tuple[object, ...], ...]:
        identities = []
        current = Path(path.anchor)
        for component in ((), *[(part,) for part in path.parts[1:]]):
            if component:
                current /= component[0]
            info = os.lstat(current)
            leaf = current == path
            if stat.S_ISLNK(info.st_mode) or (leaf and not stat.S_ISREG(info.st_mode)) or (not leaf and not stat.S_ISDIR(info.st_mode)):
                raise ProviderReceiptCaptureError(f"{label} path walk found unsafe component")
            if leaf:
                identities.append((str(current), *(getattr(info, field) for field in leaf_fields)))
            else:
                # Do not bind child-sensitive directory timestamps: unrelated
                # sibling churn does not alter this accepted leaf's authority.
                identities.append((str(current), info.st_dev, info.st_ino, info.st_mode))
        return tuple(identities)

    before_walk = walk()
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise ProviderReceiptCaptureError(f"{label} cannot be opened safely") from error
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ProviderReceiptCaptureError(f"{label} must be a regular file")
        if before.st_size != expected_size_bytes:
            raise ProviderReceiptCaptureError(
                f"{label} size drift: {before.st_size} != {expected_size_bytes}"
            )
        chunks: list[bytes] = []
        remaining = expected_size_bytes + 1
        while remaining:
            chunk = os.read(descriptor, min(1024 * 1024, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if any(getattr(before, field) != getattr(after, field) for field in leaf_fields):
        raise ProviderReceiptCaptureError(f"{label} changed during stable read")
    if len(payload) != expected_size_bytes:
        raise ProviderReceiptCaptureError(f"{label} byte count drift")
    digest = hashlib.sha256(payload).hexdigest()
    if digest != expected_sha256:
        raise ProviderReceiptCaptureError(f"{label} SHA-256 drift")
    after_walk = walk()
    if before_walk != after_walk:
        raise ProviderReceiptCaptureError(f"{label} root/ancestor/leaf identity drift")
    try:
        reopened = os.open(path, flags)
    except OSError as error:
        raise ProviderReceiptCaptureError(f"{label} cannot be reopened safely") from error
    try:
        reopened_info = os.fstat(reopened)
    finally:
        os.close(reopened)
    if any(getattr(reopened_info, field) != getattr(before, field) for field in leaf_fields):
        raise ProviderReceiptCaptureError(f"{label} changed before post-read reopen")
    return payload, {
        "path": str(path),
        "size_bytes": expected_size_bytes,
        "sha256": expected_sha256,
    }


def _decode_canonical_authority(payload: bytes, label: str) -> dict[str, object]:
    value = plan.strict_json_bytes(payload, label)
    if type(value) is not dict:
        raise ProviderReceiptCaptureError(f"{label} must be an exact JSON object")
    if plan.canonical_json_bytes(value) != payload:
        raise ProviderReceiptCaptureError(f"{label} must use canonical JSON bytes")
    return value


def _manifest_provider_order(manifest: Mapping[str, object]) -> list[tuple[str, int, str]]:
    rows = manifest["expected_provider_objects"]
    if type(rows) is not list or len(rows) != 13:
        raise ProviderReceiptCaptureError("deployment manifest must contain exactly 13 provider rows")
    result: list[tuple[str, int, str]] = []
    for index, row_value in enumerate(rows):
        if type(row_value) is not dict:
            raise ProviderReceiptCaptureError(f"deployment provider row {index} is not an object")
        sample = _require_exact_str(row_value.get("sample_id"), f"manifest sample {index}")
        part = row_value.get("part_index")
        if type(part) is not int or part < 0:
            raise ProviderReceiptCaptureError(f"manifest part index {index} is invalid")
        uri = _require_exact_str(row_value.get("uri"), f"manifest URI {index}")
        result.append((sample, part, uri))
    if tuple((sample, part) for sample, part, _ in result) != EXPECTED_SAMPLE_PARTS:
        raise ProviderReceiptCaptureError("deployment manifest provider order is not exact 3/3/4/3")
    if len({uri for _, _, uri in result}) != 13:
        raise ProviderReceiptCaptureError("deployment manifest provider URIs are not unique")
    return result


def validate_expected_provider_metadata(
    value: object,
    manifest: Mapping[str, object],
) -> dict[str, object]:
    if type(value) is not dict:
        raise ProviderReceiptCaptureError("expected provider metadata must be an exact object")
    _exact_keys(
        value,
        (
            "schema_version",
            "deployment_bundle_id",
            "capture_authorization_sha256",
            "contract_bindings",
            "objects",
        ),
        "expected provider metadata",
    )
    if value["schema_version"] != EXPECTED_METADATA_SCHEMA:
        raise ProviderReceiptCaptureError("expected provider metadata schema drift")
    if value["deployment_bundle_id"] != "exact4-bundle-cee202a9f2a520062573fe50":
        raise ProviderReceiptCaptureError("expected provider metadata source-bundle drift")
    if value["capture_authorization_sha256"] != CAPTURE_AUTHORIZATION_SHA256:
        raise ProviderReceiptCaptureError("provider capture authorization drift")

    bindings = value["contract_bindings"]
    if type(bindings) is not dict:
        raise ProviderReceiptCaptureError("expected contract bindings must be an exact object")
    _exact_keys(bindings, tuple(plan.CONTRACT_BINDINGS), "expected contract bindings")
    manifest_bindings = {key: manifest[key] for key in plan.CONTRACT_BINDINGS}
    for key, expected in plan.CONTRACT_BINDINGS.items():
        _require_sha256(bindings[key], f"expected binding {key}")
        if bindings[key] != expected or bindings[key] != manifest_bindings[key]:
            raise ProviderReceiptCaptureError(f"expected provider metadata binding drift: {key}")

    objects = value["objects"]
    if type(objects) is not list or len(objects) != 13:
        raise ProviderReceiptCaptureError("expected provider metadata must contain exactly 13 rows")
    manifest_order = _manifest_provider_order(manifest)
    normalized: list[dict[str, object]] = []
    for index, (row_value, manifest_key) in enumerate(zip(objects, manifest_order)):
        if type(row_value) is not dict:
            raise ProviderReceiptCaptureError(f"expected provider row {index} is not an object")
        _exact_keys(
            row_value,
            (
                "sample_id",
                "part_index",
                "uri",
                "version_id",
                "etag",
                "content_length_bytes",
                "last_modified",
            ),
            f"expected provider row {index}",
        )
        sample = _require_exact_str(row_value["sample_id"], f"expected sample {index}")
        part = row_value["part_index"]
        if type(part) is not int or part < 0:
            raise ProviderReceiptCaptureError(f"expected part index {index} is invalid")
        uri = _require_exact_str(row_value["uri"], f"expected URI {index}")
        if (sample, part, uri) != manifest_key:
            raise ProviderReceiptCaptureError(f"expected provider URI/order drift at row {index}")
        version_id = _require_exact_str(row_value["version_id"], f"expected version ID {index}")
        etag = _require_exact_str(row_value["etag"], f"expected ETag {index}")
        if not (len(etag) >= 2 and etag.startswith('"') and etag.endswith('"')):
            raise ProviderReceiptCaptureError(f"expected ETag {index} must retain HTTP quotes")
        content_length = _require_exact_int(
            row_value["content_length_bytes"], f"expected content length {index}"
        )
        last_modified = _require_exact_str(
            row_value["last_modified"], f"expected last-modified {index}"
        )
        normalized.append(
            {
                "sample_id": sample,
                "part_index": part,
                "uri": uri,
                "version_id": version_id,
                "etag": etag,
                "content_length_bytes": content_length,
                "last_modified": last_modified,
            }
        )
    return {
        "schema_version": EXPECTED_METADATA_SCHEMA,
        "deployment_bundle_id": value["deployment_bundle_id"],
        "capture_authorization_sha256": CAPTURE_AUTHORIZATION_SHA256,
        "contract_bindings": dict(plan.CONTRACT_BINDINGS),
        "objects": normalized,
    }


def validate_aws_tool_authority(value: object) -> dict[str, object]:
    if type(value) is not dict:
        raise ProviderReceiptCaptureError("AWS tool authority must be an exact object")
    _exact_keys(
        value,
        (
            "schema_version",
            "normalized_absolute_path",
            "sha256",
            "version_output",
            "module_name",
            "read_only_resident_receipt_sha256",
            "capture_authorization_sha256",
            "cluster_host_observation_class",
            "action_gates",
        ),
        "AWS tool authority",
    )
    if value["schema_version"] != AWS_TOOL_AUTHORITY_SCHEMA:
        raise ProviderReceiptCaptureError("AWS tool authority schema drift")
    path = _require_exact_str(value["normalized_absolute_path"], "AWS authority path")
    if not Path(path).is_absolute() or str(Path(path)) != path:
        raise ProviderReceiptCaptureError("AWS authority path is not normalized absolute")
    _require_sha256(value["sha256"], "AWS authority executable SHA-256")
    _require_exact_str(value["version_output"], "AWS authority version output")
    _require_exact_str(value["module_name"], "AWS authority module")
    _require_sha256(
        value["read_only_resident_receipt_sha256"],
        "AWS authority resident receipt SHA-256",
    )
    if value["capture_authorization_sha256"] != CAPTURE_AUTHORIZATION_SHA256:
        raise ProviderReceiptCaptureError("AWS authority capture authorization drift")
    _require_exact_str(
        value["cluster_host_observation_class"], "AWS authority host class"
    )
    gates = value["action_gates"]
    if type(gates) is not dict:
        raise ProviderReceiptCaptureError("AWS authority action gates must be an exact object")
    _exact_keys(gates, tuple(UNAUTHORIZED_ACTION_GATES), "AWS authority action gates")
    if gates != UNAUTHORIZED_ACTION_GATES or any(type(gates[key]) is not bool or gates[key] is not False for key in UNAUTHORIZED_ACTION_GATES):
        raise ProviderReceiptCaptureError("AWS authority action gates are not all exactly false")
    return copy.deepcopy(value)


def parse_s3_uri(uri: str) -> tuple[str, str]:
    parsed = urlsplit(uri)
    if (
        parsed.scheme != "s3"
        or not parsed.netloc
        or parsed.username is not None
        or parsed.password is not None
        or parsed.port is not None
        or parsed.query
        or parsed.fragment
        or not parsed.path.startswith("/")
        or parsed.path.startswith("//")
    ):
        raise ProviderReceiptCaptureError(f"invalid exact S3 URI: {uri!r}")
    bucket = parsed.netloc
    key = parsed.path[1:]
    if not key or any(character in uri for character in ("\x00", "\n", "\r")):
        raise ProviderReceiptCaptureError(f"invalid exact S3 URI: {uri!r}")
    return bucket, key


def _validate_aws_executable(path_text: str) -> tuple[Path, dict[str, object]]:
    path = Path(path_text)
    if not path.is_absolute():
        raise ProviderReceiptCaptureError("AWS executable path must be absolute")
    if not path.is_file() or path.is_symlink() or not os.access(path, os.X_OK):
        raise ProviderReceiptCaptureError("AWS executable must be executable, regular, and non-symlinked")
    if path.resolve(strict=True) != path:
        raise ProviderReceiptCaptureError("AWS executable path must be absolute and normalized")
    fingerprint = plan._fingerprint(path)
    return path, fingerprint


def _command_environment() -> dict[str, str]:
    # Deliberately do not inherit profiles, retry modes, endpoints, proxy
    # settings, credential locations, locale, or other ambient AWS settings.
    return dict(AWS_CHILD_ENVIRONMENT)


def _run_aws(command: list[str], label: str) -> subprocess.CompletedProcess[bytes]:
    completed = subprocess.run(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=_command_environment(),
    )
    if completed.returncode != 0:
        raise ProviderReceiptCaptureError(
            f"{label} failed with return code {completed.returncode}"
        )
    return completed


def _decode_utf8(payload: bytes, label: str) -> str:
    try:
        return payload.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise ProviderReceiptCaptureError(f"{label} is not valid UTF-8") from error


def _validated_file_binding(
    value: object,
    label: str,
    expected_size_bytes: int,
    expected_sha256: str,
) -> dict[str, object]:
    if type(value) is not dict:
        raise ProviderReceiptCaptureError(f"{label} binding must be an exact object")
    _exact_keys(value, ("path", "size_bytes", "sha256"), f"{label} binding")
    path = _require_exact_str(value["path"], f"{label} path")
    if not Path(path).is_absolute():
        raise ProviderReceiptCaptureError(f"{label} path must be absolute")
    if value["size_bytes"] != expected_size_bytes or type(value["size_bytes"]) is not int:
        raise ProviderReceiptCaptureError(f"{label} accepted size binding drift")
    if value["sha256"] != expected_sha256:
        raise ProviderReceiptCaptureError(f"{label} accepted SHA-256 binding drift")
    return dict(value)


def capture_provider_receipt(
    manifest: dict[str, object],
    expected_metadata: dict[str, object],
    expected_metadata_file_binding: dict[str, object],
    expected_metadata_acceptance_binding: dict[str, object],
    aws_tool_authority: dict[str, object],
    aws_tool_authority_file_binding: dict[str, object],
    aws_tool_acceptance_binding: dict[str, object],
    aws_executable: str,
    aws_module: str,
    observation_time: str,
    amended_live_capture_authorization: dict[str, object],
) -> dict[str, object]:
    manifest = plan._validate_deployment_manifest(manifest)
    amended_binding = plan._validate_amended_live_capture_authorization(
        amended_live_capture_authorization, manifest, revalidate_files=True
    )
    expected = validate_expected_provider_metadata(expected_metadata, manifest)
    expected_file_binding = _validated_file_binding(
        expected_metadata_file_binding,
        "expected provider metadata authority",
        EXPECTED_METADATA_AUTHORITY_SIZE_BYTES,
        EXPECTED_METADATA_AUTHORITY_SHA256,
    )
    expected_acceptance_binding = _validated_file_binding(
        expected_metadata_acceptance_binding,
        "expected provider metadata acceptance",
        EXPECTED_METADATA_ACCEPTANCE_SIZE_BYTES,
        EXPECTED_METADATA_ACCEPTANCE_SHA256,
    )
    tool_authority = validate_aws_tool_authority(aws_tool_authority)
    tool_file_binding = _validated_file_binding(
        aws_tool_authority_file_binding,
        "AWS tool authority",
        AWS_TOOL_AUTHORITY_SIZE_BYTES,
        AWS_TOOL_AUTHORITY_SHA256,
    )
    tool_acceptance_binding = _validated_file_binding(
        aws_tool_acceptance_binding,
        "AWS tool authority acceptance",
        AWS_TOOL_ACCEPTANCE_SIZE_BYTES,
        AWS_TOOL_ACCEPTANCE_SHA256,
    )
    observation_time = _require_exact_str(observation_time, "observation time")
    observation_host = os.uname().nodename
    if type(observation_host) is not str or re.fullmatch(r"login-[0-9]+\.cluster\.example", observation_host) is None:
        raise ProviderReceiptCaptureError("kernel nodename is not an accepted login FQDN")
    observation_host_class = re.sub(r"^login-[0-9]+", "login", observation_host)
    aws_module = _require_exact_str(aws_module, "AWS module")
    if aws_executable != tool_authority["normalized_absolute_path"]:
        raise ProviderReceiptCaptureError("supplied AWS path differs from accepted authority")
    if aws_module != tool_authority["module_name"]:
        raise ProviderReceiptCaptureError("supplied AWS module differs from accepted authority")
    if observation_host_class != tool_authority["cluster_host_observation_class"]:
        raise ProviderReceiptCaptureError("supplied host class differs from accepted authority")
    aws_path, aws_fingerprint = _validate_aws_executable(aws_executable)
    aws_text = str(aws_path)
    if aws_text != tool_authority["normalized_absolute_path"]:
        raise ProviderReceiptCaptureError("validated AWS path differs from accepted authority")
    if aws_fingerprint.get("path") != aws_text:
        raise ProviderReceiptCaptureError("AWS fingerprint path drift")
    if aws_fingerprint.get("sha256") != tool_authority["sha256"]:
        raise ProviderReceiptCaptureError("AWS executable SHA-256 differs from accepted authority")

    version_command = [aws_text, "--version"]
    version_result = _run_aws(version_command, "AWS version receipt")
    version_bytes = version_result.stdout + version_result.stderr
    if not version_bytes:
        raise ProviderReceiptCaptureError("AWS version receipt is empty")
    _decode_utf8(version_bytes, "AWS version receipt")
    accepted_version_bytes = (tool_authority["version_output"] + "\n").encode("utf-8")
    if version_bytes != accepted_version_bytes:
        raise ProviderReceiptCaptureError("AWS version output differs from accepted authority")

    objects: list[dict[str, object]] = []
    for index, row in enumerate(expected["objects"]):
        bucket, key = parse_s3_uri(row["uri"])
        command = [
            aws_text,
            "s3api",
            "head-object",
            "--bucket",
            bucket,
            "--key",
            key,
            "--version-id",
            row["version_id"],
            "--no-sign-request",
        ]
        result = _run_aws(command, f"provider HEAD row {index}")
        if result.stderr:
            raise ProviderReceiptCaptureError(f"provider HEAD row {index} emitted stderr")
        response = plan.strict_json_bytes(result.stdout, f"provider HEAD response {index}")
        if type(response) is not dict:
            raise ProviderReceiptCaptureError(f"provider HEAD response {index} is not an object")
        required = {
            "VersionId": row["version_id"],
            "ETag": row["etag"],
            "ContentLength": row["content_length_bytes"],
            "LastModified": row["last_modified"],
        }
        for key_name, expected_value in required.items():
            if key_name not in response:
                raise ProviderReceiptCaptureError(
                    f"provider HEAD response {index} is missing {key_name}"
                )
            if type(response[key_name]) is not type(expected_value) or response[key_name] != expected_value:
                raise ProviderReceiptCaptureError(
                    f"provider HEAD response {index} drift: {key_name}"
                )
        objects.append(
            {
                **row,
                "observation_time": observation_time,
                "provider_command": command,
                "response_digest": hashlib.sha256(result.stdout).hexdigest(),
            }
        )

    authority_bindings = {
        "expected_provider_metadata_authority": {
            **expected_file_binding,
            "acceptance_path": expected_acceptance_binding["path"],
            "acceptance_size_bytes": expected_acceptance_binding["size_bytes"],
            "acceptance_sha256": expected_acceptance_binding["sha256"],
            "deployment_bundle_id": expected["deployment_bundle_id"],
            "capture_authorization_sha256": expected["capture_authorization_sha256"],
        },
        "aws_tool_authority": {
            **tool_file_binding,
            "acceptance_path": tool_acceptance_binding["path"],
            "acceptance_size_bytes": tool_acceptance_binding["size_bytes"],
            "acceptance_sha256": tool_acceptance_binding["sha256"],
            "normalized_absolute_path": tool_authority["normalized_absolute_path"],
            "executable_sha256": tool_authority["sha256"],
            "version_output": tool_authority["version_output"],
            "module_name": tool_authority["module_name"],
            "cluster_host_observation_class": tool_authority["cluster_host_observation_class"],
            "read_only_resident_receipt_sha256": tool_authority["read_only_resident_receipt_sha256"],
            "capture_authorization_sha256": tool_authority["capture_authorization_sha256"],
        },
        "amended_live_capture_authorization": amended_binding,
    }
    receipt = {
        "schema_version": PROVIDER_RECEIPT_SCHEMA,
        "deployment_bundle_id": manifest["bundle_id"],
        "objects": objects,
        "provider_tool_receipt": {
            **aws_fingerprint,
            "version_output": tool_authority["version_output"],
            "version_output_sha256": hashlib.sha256(version_bytes).hexdigest(),
        },
        "capture_metadata": {
            "observation_host": observation_host,
            "capture_label": CAPTURE_LABEL,
            "object_body_downloaded": False,
        },
        "authority_bindings": authority_bindings,
        "action_gates": dict(UNAUTHORIZED_ACTION_GATES),
    }
    return plan._validate_provider_receipt(receipt, manifest, revalidate_external_files=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Emit one canonical exact-13 metadata-HEAD-only provider receipt."
    )
    parser.add_argument("--deployment-manifest", required=True)
    parser.add_argument("--expected-provider-metadata", required=True)
    parser.add_argument("--expected-provider-metadata-acceptance", required=True)
    parser.add_argument("--aws-tool-authority", required=True)
    parser.add_argument("--aws-tool-authority-acceptance", required=True)
    parser.add_argument("--aws", required=True)
    parser.add_argument("--aws-module", required=True)
    parser.add_argument("--observation-time", required=True)
    parser.add_argument("--amended-authorization-size-bytes", required=True, type=int)
    parser.add_argument("--amended-authorization-sha256", required=True)
    parser.add_argument("--amended-authorization-acceptance-size-bytes", required=True, type=int)
    parser.add_argument("--amended-authorization-acceptance-sha256", required=True)
    parser.add_argument("--deployment-manifest-size-bytes", required=True, type=int)
    parser.add_argument("--deployment-manifest-sha256", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        expected_manifest_path = plan._project_path_from_relative(plan.DEPLOYMENT_MANIFEST_PATH)
        fixed_roles = {
            "deployment_manifest": expected_manifest_path,
            "expected_provider_metadata": plan._project_path_from_relative(EXPECTED_METADATA_AUTHORITY_PATH),
            "expected_provider_metadata_acceptance": plan._project_path_from_relative(EXPECTED_METADATA_ACCEPTANCE_PATH),
            "aws_tool_authority": plan._project_path_from_relative(AWS_TOOL_AUTHORITY_PATH),
            "aws_tool_authority_acceptance": plan._project_path_from_relative(AWS_TOOL_ACCEPTANCE_PATH),
        }
        for name, fixed_path in fixed_roles.items():
            _require_lexical_normalized_absolute(getattr(args, name), name.replace("_", " "), expected=fixed_path)
        manifest_payload, _ = _stable_read_accepted_bytes(
            args.deployment_manifest,
            "amended deployment manifest",
            args.deployment_manifest_size_bytes,
            args.deployment_manifest_sha256,
        )
        manifest = _decode_canonical_authority(manifest_payload, "amended deployment manifest")
        expected_payload, expected_binding = _stable_read_accepted_bytes(
            args.expected_provider_metadata,
            "expected provider metadata authority",
            EXPECTED_METADATA_AUTHORITY_SIZE_BYTES,
            EXPECTED_METADATA_AUTHORITY_SHA256,
        )
        _, expected_acceptance_binding = _stable_read_accepted_bytes(
            args.expected_provider_metadata_acceptance,
            "expected provider metadata acceptance",
            EXPECTED_METADATA_ACCEPTANCE_SIZE_BYTES,
            EXPECTED_METADATA_ACCEPTANCE_SHA256,
        )
        tool_payload, tool_binding = _stable_read_accepted_bytes(
            args.aws_tool_authority,
            "AWS tool authority",
            AWS_TOOL_AUTHORITY_SIZE_BYTES,
            AWS_TOOL_AUTHORITY_SHA256,
        )
        _, tool_acceptance_binding = _stable_read_accepted_bytes(
            args.aws_tool_authority_acceptance,
            "AWS tool authority acceptance",
            AWS_TOOL_ACCEPTANCE_SIZE_BYTES,
            AWS_TOOL_ACCEPTANCE_SHA256,
        )
        expected = _decode_canonical_authority(
            expected_payload, "expected provider metadata authority"
        )
        tool_authority = _decode_canonical_authority(tool_payload, "AWS tool authority")
        amended_binding = {
            "authorization_path": plan.AMENDED_AUTHORIZATION_PATH,
            "authorization_size_bytes": args.amended_authorization_size_bytes,
            "authorization_sha256": args.amended_authorization_sha256,
            "acceptance_path": plan.AMENDED_AUTHORIZATION_ACCEPTANCE_PATH,
            "acceptance_size_bytes": args.amended_authorization_acceptance_size_bytes,
            "acceptance_sha256": args.amended_authorization_acceptance_sha256,
            "deployment_bundle_id": manifest["bundle_id"],
            "deployment_manifest_path": plan.DEPLOYMENT_MANIFEST_PATH,
            "deployment_manifest_size_bytes": args.deployment_manifest_size_bytes,
            "deployment_manifest_sha256": args.deployment_manifest_sha256,
            "provider_receipt_schema": PROVIDER_RECEIPT_SCHEMA,
            "biological_target": plan.LOCUS,
            "ordered_samples": list(plan.TARGETS),
            "provider_part_counts": {sample: count for sample, count in zip(plan.TARGETS, plan.PART_COUNTS)},
            "read_only_action_scope": plan.AMENDED_READ_ONLY_SCOPE,
        }
        receipt = capture_provider_receipt(
            manifest,
            expected,
            expected_binding,
            expected_acceptance_binding,
            tool_authority,
            tool_binding,
            tool_acceptance_binding,
            args.aws,
            args.aws_module,
            args.observation_time,
            amended_binding,
        )
        sys.stdout.buffer.write(plan.canonical_json_bytes(receipt))
        return 0
    except (ProviderReceiptCaptureError, plan.ExactFourPlanError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
