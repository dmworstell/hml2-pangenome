#!/usr/bin/env python3
"""Build the supplementary helper-dependence sensitivity figure."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

from figure_style import BLUE, GREEN, INK, LIGHT, PURPLE, RED, apply_style, finish_axis, panel_title, save_figure


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "manuscript/delta292_helper_abc_v1"
OUT = ROOT / "manuscript/figures/delta292_helper_abc"
PNG = OUT / "delta292_helper_abc_summary.png"
PDF = OUT / "delta292_helper_abc_summary.pdf"


def read_tsv(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))




def main() -> None:
    apply_style()
    equilibrium = read_tsv("helper_equilibrium_requirements.tsv")
    neutral = read_tsv("neutral_founder_hitting_bounds.tsv")

    fig, (ax_b, ax_c) = plt.subplots(2, 1, figsize=(7.1, 5.8), constrained_layout=True)


    for target, label, color in [(0.625, "10/16 floor", BLUE), (0.7391304347826086, "17/23 ceiling", RED)]:
        rows = [r for r in equilibrium if abs(float(r["target_focal_fraction"]) - target) < 1e-8]
        ax_b.plot([float(r["moi"]) for r in rows], [float(r["minimum_replication_gain_for_equilibrium"]) for r in rows], marker="o", color=color, label=label)
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlabel("Average HML-2 genomes per infected cell")
    ax_b.set_ylabel("Minimum gain with helper virus")
    ax_b.legend()
    panel_title(ax_b, "A", "Helper availability")
    finish_axis(ax_b, grid="both")

    for copies, label, color in [(1, "1 starting copy", PURPLE), (5, "5 starting copies", GREEN)]:
        rows = [r for r in neutral if int(r["starting_copies"]) == copies and abs(float(r["target_fraction"]) - 0.625) < 1e-8]
        ax_c.plot([int(r["effective_population_size"]) for r in rows], [float(r["martingale_upper_bound_probability_ever_reach_target"]) for r in rows], marker="o", color=color, label=label)
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_xlabel("Effective viral population size")
    ax_c.set_ylabel("Maximum probability of reaching 10/16")
    ax_c.legend()
    panel_title(ax_c, "B", "Neutral founder probability")
    finish_axis(ax_c, grid="both")


    save_figure(fig, PNG, PDF)
    plt.close(fig)
    print(PNG)
    print(PDF)


if __name__ == "__main__":
    main()
