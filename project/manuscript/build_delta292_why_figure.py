#!/usr/bin/env python3
"""Build the manuscript synthesis figure for Delta292 history and 'why' tests."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from figure_style import GREEN, GRAY, INK, LIGHT, ORANGE, apply_style, finish_axis, panel_title, save_figure


WORKSPACE = Path(__file__).resolve().parents[2]
PROJECT = WORKSPACE / "project"
MECHANISM = (
    WORKSPACE
    / "manuscript_figures/python/analysis/type1_cassette_mechanism/results"
)
COUNTERFACTUAL = PROJECT / "manuscript/delta292_ape_counterfactuals_v1"
OUT = PROJECT / "manuscript/figures/delta292_why"
PNG = OUT / "delta292_why_synthesis.png"
PDF = OUT / "delta292_why_synthesis.pdf"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))




def main() -> None:
    apply_style()
    switches = read_tsv(MECHANISM / "nearest_neighbor_leave_one_out_summary.tsv")
    ape = json.loads((COUNTERFACTUAL / "summary.json").read_text())

    order = ["B1", "B2", "B3", "B4", "B5", "B6", "Cass."]

    # The regional-divergence panel is already shown in the main manuscript.
    # The supplementary figure retains only analyses that add new information.
    fig = plt.figure(figsize=(7.1, 5.25), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    ax_a = fig.add_subplot(grid[0, :])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    switch_by = {row["focal_window_label"]: row for row in switches}
    vals = np.array([float(switch_by[label]["switch_fraction"]) for label in order])
    lo = np.array([float(switch_by[label]["wilson_95ci_low"]) for label in order])
    hi = np.array([float(switch_by[label]["wilson_95ci_high"]) for label in order])
    x = np.arange(len(order))
    colors = [GRAY] * 6 + [GREEN]
    ax_a.errorbar(x, vals, yerr=[vals - lo, hi - vals], fmt="none", ecolor=INK, lw=0.8, capsize=2)
    ax_a.scatter(x, vals, c=colors, s=25, zorder=3)
    ax_a.set_xticks(x, order, rotation=45, ha="right")
    ax_a.set_ylim(0, 1.05)
    ax_a.set_ylabel("Nearest-neighbour switch fraction")
    panel_title(ax_a, "A", "Nearest-neighbour switching")
    finish_axis(ax_a, grid="y")


    cf = ape["counterfactuals"]
    matrix = np.array(
        [
            [26 - cf["delta292_without_6331T_6492T"], cf["delta292_without_6331T_6492T"]],
            [cf["non_delta_with_6331T_6492T"], 57 - cf["non_delta_with_6331T_6492T"]],
        ]
    )
    from matplotlib.colors import LinearSegmentedColormap

    muted_teal = LinearSegmentedColormap.from_list(
        "muted_teal", ["#F2F3F2", "#AFC7C1", GREEN]
    )
    ax_c.imshow(matrix, cmap=muted_teal, vmin=0, vmax=matrix.max())
    ax_c.set_xticks([0, 1], ["6331T+6492T", "Other pair"], rotation=18, ha="right")
    ax_c.set_yticks([0, 1], ["Δ292", "No Δ292"])
    for row in range(2):
        for col in range(2):
            ax_c.text(col, row, str(matrix[row, col]), ha="center", va="center", fontweight="bold", color="white" if matrix[row, col] > 12 else INK)
    panel_title(ax_c, "B", "Linked substitutions")
    ax_c.tick_params(length=0)

    lesions = [
        (112, 1, "112 nt", GRAY),
        (292, 10, "Δ292", GREEN),
        (556, 1, "556 nt", ORANGE),
    ]
    for length, units, label, color in lesions:
        ax_d.scatter(length, units, s=45, color=color, zorder=3)
        ax_d.annotate(label, (length, units), xytext=(0, 7), textcoords="offset points", ha="center", fontweight="bold")
    ax_d.set_yscale("log")
    ax_d.set_ylim(0.7, 18)
    ax_d.set_xlim(50, 620)
    ax_d.set_xlabel("Deletion length (nt)")
    ax_d.set_ylabel("Independent lesion groups")
    ax_d.set_yticks([1, 2, 5, 10], ["1", "2", "5", "10"])
    panel_title(ax_d, "C", "Natural deletion recurrence")
    finish_axis(ax_d, grid="both")

    save_figure(fig, PNG, PDF)
    plt.close(fig)
    print(PNG)
    print(PDF)


if __name__ == "__main__":
    main()
