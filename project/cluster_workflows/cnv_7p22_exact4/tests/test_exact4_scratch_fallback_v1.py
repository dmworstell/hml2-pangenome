from pathlib import Path
import unittest


WORKER = Path(__file__).parents[1] / "src" / "02_exact4_alignment_worker.sh"


class ExactFourScratchFallbackTests(unittest.TestCase):
    def test_worker_uses_private_fallback_and_bounded_cleanup(self) -> None:
        text = WORKER.read_text(encoding="utf-8")
        self.assertIn("NODE_TMP=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}", text)
        self.assertIn("install -d -m 0700", text)
        self.assertIn("SCRATCH_ROOT=${NODE_TMP%/}/$WORKER_USER", text)
        self.assertIn("trap cleanup_scratch EXIT INT TERM", text)
        self.assertIn('rm -f -- "$SCRATCH"', text)
        self.assertNotIn("SLURM_TMPDIR:?", text)


if __name__ == "__main__":
    unittest.main()
