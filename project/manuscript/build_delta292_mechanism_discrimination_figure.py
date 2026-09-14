#!/usr/bin/env python3
"""Build the supplementary source-versus-effect and gain-budget figure."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from figure_style import BLUE, GOLD, GREEN, LIGHT, ORANGE, PURPLE, RED, apply_style, finish_axis, panel_title, save_figure


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "manuscript/delta292_source_vs_effect_v1"
BUDGET = ROOT / "manuscript/delta292_budget_abc_v1"
OUT = ROOT / "manuscript/figures/delta292_mechanism_discrimination"
PNG = OUT / "delta292_mechanism_discrimination.png"
PDF = OUT / "delta292_mechanism_discrimination.pdf"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    apply_style()
    source = read_tsv(SOURCE / "model_comparison.tsv")
    requirement = read_tsv(SOURCE / "source_only_opportunity_requirement.tsv")
    budget = read_tsv(BUDGET / "budget_model_comparison.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.45), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.flat

    keys = [
        ("orthology_aware_floor", "narrow", "10/16 narrow", BLUE, "o"),
        ("orthology_aware_floor", "broad", "10/16 broad", BLUE, "s"),
        ("inclusive_ceiling", "narrow", "17/23 narrow", RED, "o"),
        ("inclusive_ceiling", "broad", "17/23 broad", RED, "s"),
    ]
    for obs, asc, label, color, marker in keys:
        if obs == "orthology_aware_floor" and asc == "narrow":
            continue  # This single primary curve is shown in main Figure 5D.
        rows = [r for r in source if r["observation"] == obs and r["ascertainment_scenario"] == asc and r["model"] == "source_plus_focal_effect"]
        rows.sort(key=lambda r: float(r["source_opportunity_cv"]), reverse=True)
        ax_a.plot([float(r["source_opportunity_cv"]) for r in rows], [float(r["likelihood_ratio_effect_vs_source_only"]) for r in rows], marker=marker, color=color, label=label)
    ax_a.axhline(1, color="#172A3A", lw=0.8)
    ax_a.set_yscale("log")
    ax_a.set_xlabel("Source-weight coefficient of variation")
    ax_a.set_ylabel("Bayes factor for an\nadded Δ292 advantage")
    ax_a.legend(ncol=1, fontsize=8)
    panel_title(ax_a, "A", "Δ292 effect under alternative\ncounts and ascertainment")
    finish_axis(ax_a, grid="both")

    for obs, asc, label, color, marker in keys:
        rows = [r for r in requirement if r["observation"] == obs and r["ascertainment_scenario"] == asc]
        rows.sort(key=lambda r: float(r["source_opportunity_cv"]))
        ax_b.plot([float(r["source_opportunity_cv"]) for r in rows], [float(r["posterior_source_ratio_median"]) for r in rows], marker=marker, color=color, label=label)
    ax_b.axhline(1, color="#172A3A", lw=0.8)
    ax_b.set_xlabel("Source-weight coefficient of variation")
    ax_b.set_ylabel("Δ292 source contribution /\nother sources")
    panel_title(ax_b, "B", "Δ292 source contribution\nwithout an added effect")
    finish_axis(ax_b, grid="both")

    model_order = ["length_only", "rna_redirection", "producer_cost_avoidance", "packaging_per_genomic_rna", "rna_plus_producer", "rna_plus_packaging", "producer_plus_packaging", "all_three"]
    labels = ["Length", "RNA", "Producer", "Packaging", "RNA + producer", "RNA + packaging", "Producer + packaging", "All three"]
    colors = [LIGHT, BLUE, GREEN, GOLD, PURPLE, ORANGE, "#8B6D4F", RED]
    rows = [r for r in budget if r["prior"] == "conservative" and r["observation"] == "orthology_aware_floor"]
    by = {r["model"]: r for r in rows}
    y = np.arange(len(model_order))[::-1]
    vals = [float(by[m]["likelihood_ratio_vs_length_only"]) for m in model_order]
    ax_c.barh(y, vals, color=colors)
    ax_c.set_xscale("log")
    ax_c.set_yticks(y, labels)
    ax_c.set_xlabel("Support relative to length alone")
    panel_title(ax_c, "C", "RNA supply, producer cost\nand packaging")
    finish_axis(ax_c, grid="x")

    categories = {
        "Length only": ["length_only"],
        "One component": ["rna_redirection", "producer_cost_avoidance", "packaging_per_genomic_rna"],
        "Two components": ["rna_plus_producer", "rna_plus_packaging", "producer_plus_packaging"],
        "All three": ["all_three"],
    }
    scenarios = [("conservative", "orthology_aware_floor", "Conserv.\n10/16"), ("conservative", "inclusive_ceiling", "Conserv.\n17/23"), ("broad", "orthology_aware_floor", "Broad\n10/16"), ("broad", "inclusive_ceiling", "Broad\n17/23")]
    x = np.arange(len(scenarios))
    bottom = np.zeros(len(scenarios))
    cat_colors = [LIGHT, BLUE, PURPLE, RED]
    for (cat, models), color in zip(categories.items(), cat_colors):
        vals = []
        for prior, obs, _ in scenarios:
            subset = [r for r in budget if r["prior"] == prior and r["observation"] == obs and r["model"] in models]
            vals.append(sum(float(r["equal_model_prior_share"]) for r in subset))
        ax_d.bar(x, vals, bottom=bottom, color=color, width=0.72, label=cat)
        bottom += vals
    ax_d.set_xticks(x, [s[2] for s in scenarios])
    ax_d.set_ylim(0, 1)
    ax_d.set_ylabel("Model weight under\nequal model priors")
    # Keep the category key outside the stacked bars.  The former in-panel
    # legend obscured both the bar heights and its own labels.
    ax_d.legend(
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        borderaxespad=0,
    )
    panel_title(ax_d, "D", "Model weights by\ncomponent combination")
    finish_axis(ax_d, grid="y")

    save_figure(fig, PNG, PDF)
    plt.close(fig)
    print(PNG)
    print(PDF)


if __name__ == "__main__":
    main()
