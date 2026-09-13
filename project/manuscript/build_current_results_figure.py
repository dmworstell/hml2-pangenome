#!/usr/bin/env python3
"""Build the manuscript functional-results figure from frozen result tables."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from figure_style import BLUE, GREEN, INK, LIGHT, PURPLE, apply_style, finish_axis, panel_title, save_figure


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
REGISTER = ROOT / "working/claude_safe_functional_register_v1/results/full_functional_test_register.tsv"
ONEQ22 = ROOT / "manuscript/supplement/Table_S9_1q22_Gag_artifact_corrected_truth.tsv"
CNV_SUPPLEMENT = ROOT / "manuscript/supplement/Table_S4_CNV_assembly_artifact_qc.tsv"
OUTDIR = ROOT / "manuscript/figures/current_results"
PNG = OUTDIR / "functional_cnv_summary.png"
PDF = OUTDIR / "functional_cnv_summary.pdf"
PROVENANCE = OUTDIR / "functional_cnv_summary.provenance.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def extract(pattern: str, text: str) -> float:
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"missing pattern {pattern!r}")
    return float(match.group(1).rstrip(".;"))


def forest(ax, rows, *, xlabel: str, zero: bool = True) -> None:
    y = np.arange(len(rows))[::-1]
    for yi, row in zip(y, rows):
        label, beta, se, color, note = row
        ax.errorbar(beta, yi, xerr=1.96 * se, fmt="o", color=color, ecolor=color, capsize=2.2, ms=4.5, lw=1.1)
        ax.text(beta, yi - 0.18, note, color=color, ha="center", fontsize=7.2)
    if zero:
        ax.axvline(0, color=INK, lw=0.8)
    ax.set_yticks(y, [row[0] for row in rows])
    ax.set_ylim(-0.5, len(rows) - 0.5)
    ax.set_xlabel(xlabel)
    finish_axis(ax, grid="x")


def main() -> None:
    apply_style()
    register = {row["row_id"]: row for row in read_tsv(REGISTER)}

    gag_counts = {0: 0, 1: 0, 2: 0}
    for row in read_tsv(ONEQ22):
        gag_counts[int(row["gag_compatible_dosage"])] += 1

    reg21 = register["REG-021"]["effect_uncertainty_or_power"]
    reg22 = register["REG-022"]["effect_uncertainty_or_power"]
    reg52 = register["REG-052"]["effect_uncertainty_or_power"]
    nulls = [
        ("1q22 Gag → growth\n(n=34)", extract(r"beta=([^;]+)", reg21), extract(r"se_hc3=([^;]+)", reg21), BLUE, "P=0.61"),
        ("1q22 Gag → EBV load\n(n=116)", extract(r"beta=([^;]+)", reg22), extract(r"se_hc3=([^;]+)", reg22), BLUE, "P=0.41"),
        ("7p22.1 copy no. → EBV load\n(n=116)", extract(r"E3_total_cn: beta=([^,]+)", reg52), extract(r"se_hc3=([^,]+)", reg52), PURPLE, "P=0.92"),
    ]

    # The association panels are already present in main Figure 6.  Keep only
    # the supplementary dosage distribution and additional screens.
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.1, 3.05), constrained_layout=True)

    x = np.array(sorted(gag_counts))
    values = np.array([gag_counts[v] for v in x])
    ax_a.bar(x, values, color=[LIGHT, BLUE, GREEN], edgecolor="white")
    for xx, value in zip(x, values):
        ax_a.text(xx, value + 3, str(value), ha="center", fontweight="bold")
    ax_a.set_xticks(x)
    ax_a.set_xlabel("Sequence-compatible 1q22 Gag copies per person")
    ax_a.set_ylabel("People")
    panel_title(ax_a, "A", "1q22 Gag dosage")
    finish_axis(ax_a, grid="y")

    forest(ax_b, nulls, xlabel="Effect estimate")
    panel_title(ax_b, "B", "Additional functional screens")

    save_figure(fig, PNG, PDF)
    plt.close(fig)

    artifact_rows = [row for row in read_tsv(CNV_SUPPLEMENT) if row["manuscript_analysis_decision"] == "EXCLUDE_EXACT_ASSEMBLY_RECORD"]
    if len(artifact_rows) != 35:
        raise ValueError("terminal artifact list changed")
    provenance = {
        "schema": "hml2.functional-figure.v2",
        "inputs": {str(path.relative_to(WORKSPACE)): sha256(path) for path in (REGISTER, ONEQ22, CNV_SUPPLEMENT)},
        "outputs": {str(PNG.relative_to(WORKSPACE)): sha256(PNG), str(PDF.relative_to(WORKSPACE)): sha256(PDF)},
        "scope_guard": "Thirty-five exact assembly records are excluded; no sample or locus is excluded wholesale.",
    }
    PROVENANCE.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(PNG)
    print(PDF)


if __name__ == "__main__":
    main()
