#!/usr/bin/env python3
"""Offline tests for accepted-authority exact-four provider receipt capture."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from urllib.parse import urlsplit


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
MANIFEST_PATH = ROOT / "frozen" / "exact4_deployment_manifest.v1.json"
EXPECTED_PATH = ROOT / "frozen" / "exact4_expected_provider_metadata.v1.json"
AWS_AUTHORITY_PATH = ROOT / "frozen" / "exact4_aws_tool_authority.v1.json"
HANDOFF = ROOT.parents[1] / "handoff"
EXPECTED_ACCEPTANCE_PATH = (
    HANDOFF / "7P22_EXACT_FOUR_EXPECTED_PROVIDER_METADATA_AUTHORITY_ACCEPTANCE.md"
)
AWS_ACCEPTANCE_PATH = HANDOFF / "7P22_EXACT_FOUR_AWS_TOOL_AUTHORITY_ACCEPTANCE.md"
sys.path.insert(0, str(SRC))
import capture_exact4_provider_receipt as capture
import exact4_plan as plan


def binding(path: Path) -> dict[str, object]:
    payload = path.read_bytes()
    return {
        "path": str(path),
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


class TestStableAcceptedAuthorityConcurrency(unittest.TestCase):
    def _tree(self) -> tuple[tempfile.TemporaryDirectory[str], Path, Path, bytes]:
        temporary = tempfile.TemporaryDirectory()
        ancestor = Path(temporary.name).resolve() / "authority-parent"
        ancestor.mkdir()
        leaf = ancestor / "accepted.bin"
        payload = b"accepted immutable authority bytes\n"
        leaf.write_bytes(payload)
        return temporary, ancestor, leaf, payload

    def _read(self, leaf: Path, payload: bytes) -> bytes:
        result, binding_value = capture._stable_read_accepted_bytes(
            str(leaf), "accepted authority", len(payload), hashlib.sha256(payload).hexdigest()
        )
        self.assertEqual(binding_value["path"], str(leaf))
        return result

    def _during_first_read(self, mutation: object, leaf: Path, payload: bytes) -> bytes:
        real_read = os.read
        fired = False

        def read_once(descriptor: int, count: int) -> bytes:
            nonlocal fired
            if not fired:
                fired = True
                mutation()
            return real_read(descriptor, count)

        with mock.patch.object(capture.os, "read", side_effect=read_once):
            result = self._read(leaf, payload)
        self.assertTrue(fired)
        return result

    def test_accepts_unrelated_sibling_create_and_remove(self) -> None:
        for operation in ("create", "remove"):
            with self.subTest(operation=operation):
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

                    self.assertEqual(self._during_first_read(churn, leaf, payload), payload)
                    self.assertNotEqual(ancestor.stat().st_mtime_ns, before_mtime)
                    self.assertEqual(leaf.read_bytes(), payload)
                finally:
                    temporary.cleanup()

    def test_refuses_real_ancestor_and_leaf_changes(self) -> None:
        mutations = ("ancestor_directory", "ancestor_symlink", "leaf_mutation", "leaf_replacement")
        for mutation_kind in mutations:
            with self.subTest(mutation=mutation_kind):
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

                    with self.assertRaises(capture.ProviderReceiptCaptureError):
                        self._during_first_read(mutate, leaf, payload)
                finally:
                    temporary.cleanup()


class ProviderCaptureFixture(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.tmp = Path(self.temporary.name).resolve()
        current_pins = tuple((relative, plan.sha256_file(ROOT.parents[1] / relative)) for relative, _ in plan.PREREQUISITE_PINS)
        with mock.patch.object(plan, "PREREQUISITE_PINS", current_pins):
            self.manifest = plan.build_deployment_manifest(ROOT)
        self.manifest_path = self.tmp / "exact4_deployment_manifest.v1.json"
        self.manifest_path.write_bytes(plan.canonical_json_bytes(self.manifest))
        self.authorization_path = self.tmp / "amended-authorization.md"
        self.authorization_path.write_bytes(b"fixture amended authorization\n")
        self.authorization_acceptance_path = self.tmp / "amended-authorization-acceptance.md"
        self.authorization_acceptance_path.write_bytes(b"fixture amended authorization acceptance\n")
        self.expected = plan._load_json_file(EXPECTED_PATH, "expected provider metadata")
        self.tool = plan._load_json_file(AWS_AUTHORITY_PATH, "AWS tool authority")
        self.expected_path = self._copy(EXPECTED_PATH, "expected.json")
        self.expected_acceptance_path = self._copy(
            EXPECTED_ACCEPTANCE_PATH, "expected-acceptance.md"
        )
        self.tool_path = self._copy(AWS_AUTHORITY_PATH, "aws-authority.json")
        self.tool_acceptance_path = self._copy(
            AWS_ACCEPTANCE_PATH, "aws-acceptance.md"
        )
        self.fingerprint_sha = self.tool["sha256"]
        self.version_bytes = (self.tool["version_output"] + "\n").encode("utf-8")
        self.version_failure = False
        self.head_failure_version: str | None = None
        self.response_drift: str | None = None
        self.nodename = "login-03.cluster.example"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _copy(self, source: Path, name: str) -> Path:
        destination = self.tmp / name
        destination.write_bytes(source.read_bytes())
        return destination

    def command(self) -> list[str]:
        return [
            "--deployment-manifest",
            str(self.manifest_path),
            "--expected-provider-metadata",
            str(self.expected_path),
            "--expected-provider-metadata-acceptance",
            str(self.expected_acceptance_path),
            "--aws-tool-authority",
            str(self.tool_path),
            "--aws-tool-authority-acceptance",
            str(self.tool_acceptance_path),
            "--aws",
            self.tool["normalized_absolute_path"],
            "--aws-module",
            self.tool["module_name"],
            "--observation-time",
            "2026-07-15T21:00:00Z",
            "--amended-authorization-size-bytes",
            str(self.authorization_path.stat().st_size),
            "--amended-authorization-sha256",
            hashlib.sha256(self.authorization_path.read_bytes()).hexdigest(),
            "--amended-authorization-acceptance-size-bytes",
            str(self.authorization_acceptance_path.stat().st_size),
            "--amended-authorization-acceptance-sha256",
            hashlib.sha256(self.authorization_acceptance_path.read_bytes()).hexdigest(),
            "--deployment-manifest-size-bytes",
            str(self.manifest_path.stat().st_size),
            "--deployment-manifest-sha256",
            hashlib.sha256(self.manifest_path.read_bytes()).hexdigest(),
        ]

    def fake_validate_aws(self, path_text: str) -> tuple[Path, dict[str, object]]:
        path = Path(path_text)
        return path, {
            "path": str(path),
            "size_bytes": 987654,
            "sha256": self.fingerprint_sha,
        }

    def fake_run(self, command: list[str], **_: object) -> subprocess.CompletedProcess[bytes]:
        if command[1:] == ["--version"]:
            return subprocess.CompletedProcess(
                command,
                7 if self.version_failure else 0,
                stdout=b"",
                stderr=self.version_bytes,
            )
        self.assertEqual(command[1:3], ["s3api", "head-object"])
        self.assertEqual(len(command), 10)
        self.assertEqual(command[3], "--bucket")
        self.assertEqual(command[5], "--key")
        self.assertEqual(command[7], "--version-id")
        self.assertEqual(command[9], "--no-sign-request")
        version = command[8]
        row = next(item for item in self.expected["objects"] if item["version_id"] == version)
        parsed = urlsplit(row["uri"])
        self.assertEqual(command[4], parsed.netloc)
        self.assertEqual(command[6], parsed.path[1:])
        if version == self.head_failure_version:
            return subprocess.CompletedProcess(command, 13, stdout=b"", stderr=b"failed\n")
        response = {
            "AcceptRanges": "bytes",
            "ContentLength": row["content_length_bytes"],
            "ETag": row["etag"],
            "LastModified": row["last_modified"],
            "Metadata": {},
            "VersionId": row["version_id"],
        }
        if self.response_drift == "ETag":
            response["ETag"] = '"drift"'
        elif self.response_drift == "ContentLength":
            response["ContentLength"] += 1
        elif self.response_drift == "LastModified":
            response["LastModified"] = "drift"
        elif self.response_drift == "VersionId":
            response["VersionId"] = "drift"
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=plan.canonical_json_bytes(response),
            stderr=b"",
        )

    def run_main(
        self, arguments: list[str] | None = None
    ) -> tuple[int, bytes, bytes, mock.Mock]:
        stdout_bytes = io.BytesIO()
        stderr_bytes = io.BytesIO()
        stdout = io.TextIOWrapper(stdout_bytes, encoding="utf-8", write_through=True)
        stderr = io.TextIOWrapper(stderr_bytes, encoding="utf-8", write_through=True)
        with (
            mock.patch.object(capture, "_validate_aws_executable", self.fake_validate_aws),
            mock.patch.object(capture.subprocess, "run", side_effect=self.fake_run) as runner,
            mock.patch.object(capture.os, "uname", return_value=type("Uname", (), {"nodename": self.nodename})()),
            mock.patch.object(plan, "_project_path_from_relative", side_effect=lambda relative: {
                plan.DEPLOYMENT_MANIFEST_PATH: self.manifest_path,
                plan.AMENDED_AUTHORIZATION_PATH: self.authorization_path,
                plan.AMENDED_AUTHORIZATION_ACCEPTANCE_PATH: self.authorization_acceptance_path,
                capture.EXPECTED_METADATA_AUTHORITY_PATH: self.expected_path,
                capture.EXPECTED_METADATA_ACCEPTANCE_PATH: self.expected_acceptance_path,
                capture.AWS_TOOL_AUTHORITY_PATH: self.tool_path,
                capture.AWS_TOOL_ACCEPTANCE_PATH: self.tool_acceptance_path,
            }[relative]),
            mock.patch.object(sys, "stdout", stdout),
            mock.patch.object(sys, "stderr", stderr),
        ):
            return_code = capture.main(arguments or self.command())
            stdout.flush()
            stderr.flush()
            output = stdout_bytes.getvalue()
            errors = stderr_bytes.getvalue()
        stdout.detach()
        stderr.detach()
        return return_code, output, errors, runner


class TestAcceptedExactThirteenCapture(ProviderCaptureFixture):
    def test_canonical_receipt_exact_order_commands_and_child_environment(self) -> None:
        return_code, output, errors, runner = self.run_main()
        self.assertEqual(return_code, 0, errors.decode())
        receipt = plan.strict_json_bytes(output, "provider receipt")
        self.assertEqual(output, plan.canonical_json_bytes(receipt))
        self.assertEqual(len(receipt["objects"]), 13)
        self.assertEqual(
            [(row["sample_id"], row["part_index"]) for row in receipt["objects"]],
            list(capture.EXPECTED_SAMPLE_PARTS),
        )
        self.assertEqual(runner.call_count, 14)
        commands = [call.args[0] for call in runner.call_args_list]
        self.assertEqual(commands[0], [self.tool["normalized_absolute_path"], "--version"])
        for index, (command, row) in enumerate(zip(commands[1:], receipt["objects"])):
            bucket, key = capture.parse_s3_uri(row["uri"])
            self.assertEqual(
                command,
                [
                    self.tool["normalized_absolute_path"],
                    "s3api",
                    "head-object",
                    "--bucket",
                    bucket,
                    "--key",
                    key,
                    "--version-id",
                    row["version_id"],
                    "--no-sign-request",
                ],
                index,
            )
            self.assertEqual(row["provider_command"], command)
        for call in runner.call_args_list:
            self.assertEqual(call.kwargs["env"], capture.AWS_CHILD_ENVIRONMENT)
            self.assertIs(call.kwargs["stdin"], subprocess.DEVNULL)
            self.assertIs(call.kwargs["stdout"], subprocess.PIPE)
            self.assertIs(call.kwargs["stderr"], subprocess.PIPE)
            self.assertIs(call.kwargs["check"], False)
        forbidden_environment = {
            "AWS_PROFILE",
            "AWS_DEFAULT_PROFILE",
            "AWS_RETRY_MODE",
            "AWS_ENDPOINT_URL",
            "AWS_CA_BUNDLE",
            "HTTP_PROXY",
            "HTTPS_PROXY",
            "HOME",
        }
        self.assertTrue(forbidden_environment.isdisjoint(capture.AWS_CHILD_ENVIRONMENT))
        self.assertEqual(capture.AWS_CHILD_ENVIRONMENT["AWS_MAX_ATTEMPTS"], "1")

    def test_receipt_hard_binds_all_four_accepted_files_tool_and_false_gates(self) -> None:
        return_code, output, errors, _ = self.run_main()
        self.assertEqual(return_code, 0, errors.decode())
        receipt = plan.strict_json_bytes(output, "provider receipt")
        bindings = receipt["authority_bindings"]
        expected_binding = bindings["expected_provider_metadata_authority"]
        self.assertEqual(expected_binding["size_bytes"], 5801)
        self.assertEqual(expected_binding["sha256"], capture.EXPECTED_METADATA_AUTHORITY_SHA256)
        self.assertEqual(expected_binding["acceptance_size_bytes"], 2089)
        self.assertEqual(
            expected_binding["acceptance_sha256"],
            capture.EXPECTED_METADATA_ACCEPTANCE_SHA256,
        )
        tool_binding = bindings["aws_tool_authority"]
        self.assertEqual(tool_binding["size_bytes"], 1165)
        self.assertEqual(tool_binding["sha256"], capture.AWS_TOOL_AUTHORITY_SHA256)
        self.assertEqual(tool_binding["acceptance_size_bytes"], 1855)
        self.assertEqual(
            tool_binding["acceptance_sha256"], capture.AWS_TOOL_ACCEPTANCE_SHA256
        )
        self.assertEqual(
            tool_binding["normalized_absolute_path"], self.tool["normalized_absolute_path"]
        )
        self.assertEqual(tool_binding["executable_sha256"], self.tool["sha256"])
        self.assertEqual(tool_binding["version_output"], self.tool["version_output"])
        self.assertEqual(tool_binding["module_name"], self.tool["module_name"])
        self.assertEqual(
            tool_binding["cluster_host_observation_class"],
            self.tool["cluster_host_observation_class"],
        )
        self.assertEqual(receipt["action_gates"], capture.UNAUTHORIZED_ACTION_GATES)
        self.assertEqual(len(receipt["action_gates"]), 15)
        self.assertTrue(all(value is False for value in receipt["action_gates"].values()))
        self.assertIs(receipt["capture_metadata"]["object_body_downloaded"], False)
        self.assertEqual(
            receipt["provider_tool_receipt"]["version_output"], self.tool["version_output"]
        )

    def test_full_receipt_passes_directly_to_planner_and_all_identity_layers(self) -> None:
        return_code, output, errors, _ = self.run_main()
        self.assertEqual(return_code, 0, errors.decode())
        full_receipt = plan.strict_json_bytes(output, "provider receipt")
        module_spec = importlib.util.spec_from_file_location("exact4_workflow_fixture_helpers", ROOT / "tests" / "test_exact4_workflow.py")
        assert module_spec is not None and module_spec.loader is not None
        helpers = importlib.util.module_from_spec(module_spec)
        module_spec.loader.exec_module(helpers)
        tool_root = self.tmp / "planner-tools"
        tools = helpers.fake_tools(tool_root)
        fixture_receipt, fixture_paths = helpers.synthetic_fixture_environment(
            self.manifest, tools, self.tmp / "planner-fixture-environment"
        )
        resident = helpers.synthetic_resident(
            self.manifest, tool_root, tools=tools,
            tool_fixture_receipt=fixture_receipt,
        )
        derived = tuple(copy.deepcopy(resident[key]) for key in ("assemblies", "element_identities", "bed_authorities", "combined_row_sets", "assembly_encoded_counts", "sample_bed_authorities"))
        path_map = {
            plan.DEPLOYMENT_MANIFEST_PATH: self.manifest_path,
            plan.AMENDED_AUTHORIZATION_PATH: self.authorization_path,
            plan.AMENDED_AUTHORIZATION_ACCEPTANCE_PATH: self.authorization_acceptance_path,
            **fixture_paths,
        }
        with (
            mock.patch.object(plan, "_derive_resident_biology", return_value=derived),
            mock.patch.object(plan, "_project_path_from_relative", side_effect=lambda relative: path_map.get(relative, ROOT.parents[1] / relative)),
            mock.patch.object(plan.os, "uname", return_value=type("Uname", (), {"nodename": "login-03.cluster.example"})()),
        ):
            resident = plan._validate_resident_receipt(resident, self.manifest)
            target_plan = plan.build_target_only_plan(self.manifest, resident, full_receipt)
        self.assertEqual(target_plan["provider_receipt_sha256"], plan.digest_value(full_receipt))
        self.assertEqual(target_plan["immutable_run_material"]["provider_identity"]["authority_bindings"], full_receipt["authority_bindings"])
        self.assertEqual(target_plan["immutable_run_material"]["provider_identity"]["action_gates"], full_receipt["action_gates"])

    def test_output_is_deterministic_and_no_output_option_exists(self) -> None:
        first_code, first, first_errors, _ = self.run_main()
        second_code, second, second_errors, _ = self.run_main()
        self.assertEqual((first_code, first_errors), (0, b""))
        self.assertEqual((second_code, second_errors), (0, b""))
        self.assertEqual(first, second)
        with mock.patch.object(sys, "stderr", io.StringIO()), self.assertRaises(SystemExit):
            capture.build_parser().parse_args([*self.command(), "--output", "receipt.json"])


class TestAcceptedAuthorityRefusals(ProviderCaptureFixture):
    def test_every_authority_role_alias_refuses_before_first_subprocess(self) -> None:
        options = (
            "--deployment-manifest",
            "--expected-provider-metadata",
            "--expected-provider-metadata-acceptance",
            "--aws-tool-authority",
            "--aws-tool-authority-acceptance",
        )
        for option in options:
            original = self.command()[self.command().index(option) + 1]
            aliases = (
                "relative-authority",
                original.rsplit("/", 1)[0] + "/./" + original.rsplit("/", 1)[1],
                original.rsplit("/", 1)[0] + "/alias/../" + original.rsplit("/", 1)[1],
                original.replace("/", "//", 1),
                original + "/",
            )
            for alias in aliases:
                with self.subTest(option=option, alias=alias):
                    arguments = self.command()
                    arguments[arguments.index(option) + 1] = alias
                    return_code, output, _, runner = self.run_main(arguments)
                    self.assertNotEqual(return_code, 0)
                    self.assertEqual(output, b"")
                    runner.assert_not_called()

    def test_role_swaps_symlink_and_identical_bytes_other_path_refuse_before_subprocess(self) -> None:
        for left, right in (
            ("--expected-provider-metadata", "--expected-provider-metadata-acceptance"),
            ("--aws-tool-authority", "--aws-tool-authority-acceptance"),
        ):
            arguments = self.command()
            arguments[arguments.index(left) + 1], arguments[arguments.index(right) + 1] = arguments[arguments.index(right) + 1], arguments[arguments.index(left) + 1]
            code, output, _, runner = self.run_main(arguments)
            self.assertNotEqual(code, 0); self.assertEqual(output, b""); runner.assert_not_called()
        link = self.tmp / "expected-link.json"; link.symlink_to(self.expected_path)
        copied = self.tmp / "identical-elsewhere.json"; copied.write_bytes(self.expected_path.read_bytes())
        for bad_path in (link, copied):
            arguments = self.command()
            arguments[arguments.index("--expected-provider-metadata") + 1] = str(bad_path)
            code, output, _, runner = self.run_main(arguments)
            self.assertNotEqual(code, 0); self.assertEqual(output, b""); runner.assert_not_called()

    def test_rewritten_metadata_rejected_before_fake_matching_head(self) -> None:
        rewritten = copy.deepcopy(self.expected)
        old = rewritten["objects"][0]["last_modified"]
        rewritten["objects"][0]["last_modified"] = "X" * len(old)
        payload = plan.canonical_json_bytes(rewritten)
        self.assertEqual(len(payload), capture.EXPECTED_METADATA_AUTHORITY_SIZE_BYTES)
        self.expected_path.write_bytes(payload)
        self.expected = rewritten  # The fake HEAD would match the rewrite if it were called.
        return_code, output, _, runner = self.run_main()
        self.assertNotEqual(return_code, 0)
        self.assertEqual(output, b"")
        self.assertEqual(runner.call_count, 0)

    def test_each_acceptance_and_tool_authority_rewrite_rejected_before_aws(self) -> None:
        paths = (
            self.expected_acceptance_path,
            self.tool_path,
            self.tool_acceptance_path,
        )
        for path in paths:
            with self.subTest(path=path.name):
                original = path.read_bytes()
                changed = bytes([original[0] ^ 1]) + original[1:]
                self.assertEqual(len(changed), len(original))
                path.write_bytes(changed)
                return_code, output, _, runner = self.run_main()
                self.assertNotEqual(return_code, 0)
                self.assertEqual(output, b"")
                self.assertEqual(runner.call_count, 0)
                path.write_bytes(original)

    def test_relative_and_symlinked_accepted_authorities_rejected(self) -> None:
        with self.assertRaisesRegex(capture.ProviderReceiptCaptureError, "absolute"):
            capture._stable_read_accepted_bytes(
                "relative.json",
                "test authority",
                capture.EXPECTED_METADATA_AUTHORITY_SIZE_BYTES,
                capture.EXPECTED_METADATA_AUTHORITY_SHA256,
            )
        link = self.tmp / "metadata-link.json"
        link.symlink_to(self.expected_path)
        with self.assertRaises(capture.ProviderReceiptCaptureError):
            capture._stable_read_accepted_bytes(
                str(link),
                "test authority",
                capture.EXPECTED_METADATA_AUTHORITY_SIZE_BYTES,
                capture.EXPECTED_METADATA_AUTHORITY_SHA256,
            )

    def test_supplied_path_and_module_must_match_authority(self) -> None:
        mutations = (
            ("--aws", "/absolute/wrong/aws"),
            ("--aws-module", "awscli/wrong"),
        )
        for option, replacement in mutations:
            with self.subTest(option=option):
                arguments = self.command()
                arguments[arguments.index(option) + 1] = replacement
                return_code, output, _, runner = self.run_main(arguments)
                self.assertNotEqual(return_code, 0)
                self.assertEqual(output, b"")
                self.assertEqual(runner.call_count, 0)

    def test_kernel_nodename_is_sole_host_authority_and_refuses_before_aws(self) -> None:
        for hostname in ("login-3", "LOGIN-3.cluster.example", "login-3.cluster.example.", "compute-3.cluster.example", "login-x.cluster.example"):
            with self.subTest(hostname=hostname):
                self.nodename = hostname
                return_code, output, _, runner = self.run_main()
                self.assertNotEqual(return_code, 0)
                self.assertEqual(output, b"")
                self.assertEqual(runner.call_count, 0)
        self.assertNotIn("--observation-host", capture.build_parser().format_help())
        self.assertNotIn("--observation-host-class", capture.build_parser().format_help())

    def test_each_amended_authorization_binding_drift_refuses_before_aws(self) -> None:
        for option, replacement in (
            ("--amended-authorization-sha256", "0" * 64),
            ("--amended-authorization-acceptance-sha256", "1" * 64),
            ("--deployment-manifest-sha256", "2" * 64),
            ("--amended-authorization-size-bytes", "1"),
        ):
            with self.subTest(option=option):
                arguments = self.command()
                arguments[arguments.index(option) + 1] = replacement
                return_code, output, _, runner = self.run_main(arguments)
                self.assertNotEqual(return_code, 0)
                self.assertEqual(output, b"")
                self.assertEqual(runner.call_count, 0)

    def test_executable_hash_must_match_before_version_or_head(self) -> None:
        self.fingerprint_sha = "0" * 64
        return_code, output, _, runner = self.run_main()
        self.assertNotEqual(return_code, 0)
        self.assertEqual(output, b"")
        self.assertEqual(runner.call_count, 0)


class TestNoRetryAndResponseRefusals(ProviderCaptureFixture):
    def test_version_mismatch_uses_exactly_one_subprocess_and_refuses(self) -> None:
        self.version_bytes = b"aws-cli/drift\n"
        return_code, output, _, runner = self.run_main()
        self.assertNotEqual(return_code, 0)
        self.assertEqual(output, b"")
        self.assertEqual(runner.call_count, 1)
        self.assertEqual(runner.call_args_list[0].args[0][1:], ["--version"])

    def test_version_failure_uses_exactly_one_subprocess_and_refuses(self) -> None:
        self.version_failure = True
        return_code, output, _, runner = self.run_main()
        self.assertNotEqual(return_code, 0)
        self.assertEqual(output, b"")
        self.assertEqual(runner.call_count, 1)

    def test_head_failure_has_one_source_attempt_and_no_retry(self) -> None:
        failed_version = self.expected["objects"][4]["version_id"]
        self.head_failure_version = failed_version
        return_code, output, _, runner = self.run_main()
        self.assertNotEqual(return_code, 0)
        self.assertEqual(output, b"")
        commands = [call.args[0] for call in runner.call_args_list]
        self.assertEqual(len(commands), 6)  # one version plus rows 0..4
        self.assertEqual(
            sum(failed_version in command for command in commands),
            1,
        )

    def test_each_provider_metadata_drift_refuses_without_second_attempt(self) -> None:
        for field in ("ETag", "ContentLength", "LastModified", "VersionId"):
            with self.subTest(field=field):
                self.response_drift = field
                return_code, output, _, runner = self.run_main()
                self.assertNotEqual(return_code, 0)
                self.assertEqual(output, b"")
                self.assertEqual(runner.call_count, 2)
                self.response_drift = None


if __name__ == "__main__":
    unittest.main()
