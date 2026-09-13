from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ROOT.parents[1]


class ReceiptControllerTests(unittest.TestCase):
    def test_active_shell_uses_explicit_cpython311_runtime(self) -> None:
        for path in (
            ROOT / "src/01_cnv_controller_v2.sh",
            ROOT / "src/02_cnv_alignment_worker_v2.sh",
            ROOT / "src/03_cnv_analysis_v2.sh",
        ):
            text = path.read_text()
            self.assertIn("PYTHON=${HML2_PYTHON:-python3.11}", text, path)
            self.assertIn('command -v "$PYTHON"', text, path)
            executable_lines = [
                line for line in text.splitlines()
                if line.strip().startswith("python3 ")
            ]
            self.assertEqual(executable_lines, [], path)

    def test_controller_has_no_scheduler_dependency(self) -> None:
        submitter = (ROOT / "src/01_cnv_controller_v2.sh").read_text()
        controller = (ROOT / "src/04_cnv_receipt_controller_v2.sh").read_text()
        self.assertNotIn("--dependency", submitter)
        self.assertNotIn("--dependency", controller)
        self.assertIn("part.get(\"receipt\"", controller)
        self.assertIn("--begin=now+5minutes", controller)
        self.assertIn("03_cnv_analysis_v2.sh", submitter)

    def test_legacy_dependency_controllers_remain_unreachable(self) -> None:
        legacy_paths = (
            ROOT / "original/01_cnv_controller.sh",
            PROJECT_ROOT
            / "generative/cluster_cnv_retry/live_snapshot/01_cnv_controller.sh",
            PROJECT_ROOT
            / (
                "generative/cluster_cnv_retry/live_snapshot/"
                "01_cnv_controller_retry_v1.reviewed.sh"
            ),
        )
        for legacy in legacy_paths:
            self.assertTrue(legacy.is_file(), legacy)
            legacy_text = legacy.read_text()
            self.assertTrue(
                "--dependency" in legacy_text
                or "afterok" in legacy_text
                or "afterany" in legacy_text,
                legacy,
            )

        active_paths = tuple(sorted((ROOT / "src").glob("*.sh"))) + (
            PROJECT_ROOT
            / (
                "generative/cluster_cnv_retry/live_snapshot/"
                "01_cnv_controller_retry_v2.reviewed.sh"
            ),
        )
        for active in active_paths:
            self.assertTrue(active.is_file(), active)
            active_text = active.read_text()
            self.assertNotIn("--dependency", active_text, active)
            self.assertNotIn("afterok", active_text.lower(), active)
            self.assertNotIn("afterany", active_text.lower(), active)
            for legacy in legacy_paths:
                self.assertNotIn(
                    legacy.relative_to(PROJECT_ROOT).as_posix(),
                    active_text,
                    active,
                )

        retry_v2_text = active_paths[-1].read_text()
        self.assertIn("default-closed", retry_v2_text)
        self.assertNotIn("sbatch", retry_v2_text)


if __name__ == "__main__":
    unittest.main()
