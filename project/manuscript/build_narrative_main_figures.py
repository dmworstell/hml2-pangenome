#!/usr/bin/env python3
"""Build the seven discovery-led manuscript figures from accepted current inputs."""

from __future__ import annotations

import csv
import itertools
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageOps
from structural_catalog_summary import structural_summary, write_structural_tables

from retained_panels import draw_eightq_network, draw_type_state_counts, draw_duplicated_groups

from figure_style import (
    BLUE,
    GOLD,
    GREEN,
    GRAY,
    GRID,
    INK,
    LIGHT,
    MUTED,
    ORANGE,
    PURPLE,
    RED,
    SKY,
    apply_style,
    finish_axis,
    panel_title,
    save_figure,
)


PROJECT = Path(__file__).resolve().parents[1]
WORKSPACE = PROJECT.parent
OUT = PROJECT / "manuscript/figures/narrative"
SUPPLEMENT = PROJECT / "manuscript/supplement"
CATALOG = (
    PROJECT
    / "results/resolved_manuscript_catalog_20260914/"
    "combined_hml2_orf_analysis.RESOLVED.tsv"
)
SEVENP22 = PROJECT / "working/sevenp22_proxy_resolution_agent/haplotype_copy_number_truth.tsv"
ONEP31 = PROJECT / "working/onep31b_array_recovery_agent/haplotype_array_reconciliation.tsv"
IGSR_SAMPLES = WORKSPACE / "HML2_ProjectResources/data/ref/igsr_samples.tsv"
PHY = OUT / "phylogeny"
PHY_SOURCE = (
    PROJECT
    / "manuscript/source_snapshot/figures/Main/Fig6_Phylogeny_recombination"
)
PROCESSED_LOCI = WORKSPACE / "HML2_project_data/processed_loci"
SUBFAMILY_AUTHORITY = (
    PROJECT
    / "working/type1_ltr_authority_resolution_agent/subfamily_authority.tsv"
)
TREE_SUBFAMILY_FALLBACKS = {
    # Same duplicated locus groups.
    "11p15.4b": "LTR5A",
    "1p36.21a": "LTR5B",
    # Direct LTR-reference matches in the resident ape-state audit.
    "20q11.22": "LTR5B",
    "21p13": "LTR5A",
    "22p13": "LTR5A",
    "22q11.23": "LTR5B",
    "9q34.11": "LTR5A",
    "9q34.3": "LTR5B",
    "Xq28b": "LTR5B",
    # Both recovered resident ape LTR tracts match the LTR5A reference.
    "15p13a": "LTR5A",
}
MECHANISM = (
    WORKSPACE
    / "manuscript_figures/python/analysis/"
    "type1_cassette_mechanism/results"
)
TYPE1_HYPOTHESES = (
    PROJECT / "working/type1_selection_hypotheses_v1/results/locus_inputs.tsv"
)
SOURCE_EFFECT = PROJECT / "manuscript/delta292_source_vs_effect_v1"
BIORENDER = (
    PROJECT
    / "manuscript/source_snapshot/figures/Main/Fig1_Structural_landscape/"
    "A_Motivation_Longreads.png"
)
LEADS = (
    PROJECT
    / "working/claude_safe_positive_lead_followup_v1/results/"
    "lead_reanalysis_summary.tsv"
)
SLC_INPUTS = (
    PROJECT
    / "working/hml2_functional_evidence_synthesis_claude_v2/results/"
    "synthesis_slc44a5_inputs.tsv"
)
GROWTH_GROUPS = (
    PROJECT
    / "working/direct_ebv_fitness_screen_v1/results/Im2012_intrinsic_growth/"
    "growth_6q14_group_summary.tsv"
)
TWELVE_MODELS = (
    PROJECT
    / "working/twelveq13_2_anti_cd20_followup_v1/state_robustness/results/"
    "adjusted_state_models.tsv"
)
TWELVE = (
    PROJECT
    / "working/claude_safe_positive_lead_followup_v1/results/"
    "twelveq13_complement_causal_discrimination.tsv"
)
SIXQ_LOO = (
    PROJECT
    / "working/direct_ebv_fitness_screen_v1/results/Im2012_intrinsic_growth/"
    "growth_6q14_leave_one_person_out.tsv"
)
TYPE1_LOCUS_AUDIT = (
    PROJECT / "working/type1_causal_explanation_v1/typeI_locus_type_audit.tsv"
)

TYPE_I = "#356F6A"
TYPE_II = "#765A78"
TREE_LTR5HS = "#0072B2"
TREE_LTR5A = "#7A5195"
TREE_LTR5B = "#D55E00"
PUBLIC_ID = re.compile(r"^(?:HG|NA)\d+$")
COMPATIBLE = {"Intact", "Frameshift_at_end", "Intact_FS_End"}
PROVIRUS = {"Provirus", "Provirus_from_Multi"}
STATE_COLORS = {
    "Noncarrier call": "#E8ECEF",
    "Solo-LTR": ORANGE,
    "Fragment": SKY,
    "Provirus": GREEN,
    "Multi-copy": PURPLE,
    "Unknown": "#8A949E",
}
FEATURES = ["gag", "gag_pro_route", "gag_pro_pol_route", "env", "accessory"]
FEATURE_LABELS = ["Gag", "Gag–Pro", "Gag–Pro–Pol", "Env", "Np9 / Rec"]
MAIN_TREE_LABELS = {
    "1p31.1b",
    "1p36.21a",
    "1p36.21b",
    "1p36.21c",
    "1q22",
    "3q12.3",
    "3q13.2",
    "3q27.2",
    "4p16.1a",
    "4p16.1b",
    "5q33.3",
    "6q14.1",
    "7p22.1",
    "8p23.1a",
    "8p23.1b",
    "8p23.1c",
    "8p23.1d",
    "8p23.1e",
    "8q11.23_new",
    "10p12.1",
    "12q13.2",
    "12q24.33",
    "19p12a",
    "19p12c",
    "22q11.21",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def display_locus_name(value: str) -> str:
    """Normalize a locus name for figures without changing analysis keys."""
    label = value.removeprefix("HML-2_")
    for suffix in ("_new", "_hg38"):
        label = label.removesuffix(suffix)
    return {
        "acro_type1": "Acrocentric Type I",
        "acro_type2": "Acrocentric Type II",
        "8p23.1_duplicate_group_unresolved": "8p23.1 duplicate group",
    }.get(label, label)


def finding_title(
    ax,
    letter: str,
    title: str,
    *,
    letter_x: float = -0.105,
) -> None:
    """Short, centered panel title with a separate panel letter."""
    ax.set_title(title, loc="center", fontsize=11.4, fontweight="normal", pad=9)
    ax.text(
        letter_x,
        1.025,
        letter.upper(),
        transform=ax.transAxes,
        fontsize=11.5,
        fontweight="bold",
        color=INK,
        va="bottom",
    )


def load_catalog() -> tuple[list[dict[str, str]], list[tuple[str, str]]]:
    with CATALOG.open(newline="", encoding="utf-8") as handle:
        public_rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if PUBLIC_ID.fullmatch(row["ID"])
        ]
    excluded = Counter(
        row["analysis_exclusion_reason"]
        for row in public_rows
        if row["analysis_include"] == "0"
    )
    expected_exclusions = Counter(
        {
            "alias_duplicate_of_8q24.3c": 584,
            "assembly_artifact_not_supported_by_CNV_depth": 35,
            "duplicate_catalog_label_for_same_assembled_interval": 78,
            "non_HML2_HML11_sequence_identity": 1168,
        }
    )
    if excluded != expected_exclusions:
        raise ValueError(f"unexpected biological exclusions: {excluded}")
    rows = [row for row in public_rows if row["analysis_include"] == "1"]
    roster = sorted({(row["ID"], row["Haplotype"]) for row in rows})
    if len(rows) != 59_656 or len(roster) != 584:
        raise ValueError(f"unexpected current catalog dimensions: {len(rows)}, {len(roster)}")
    return rows, roster


def gene_call(row: dict[str, str], feature: str) -> tuple[bool, bool]:
    if feature == "gag_pro_route":
        values = [row["gag"], row["pro"]]
        callable_ = all(value not in {"", "NA"} for value in values)
        return callable_, callable_ and all(value in COMPATIBLE for value in values)
    if feature == "gag_pro_pol_route":
        values = [row["gag"], row["pro"], row["pol"]]
        callable_ = all(value not in {"", "NA"} for value in values)
        return callable_, callable_ and all(value in COMPATIBLE for value in values)
    field = "np9" if feature == "accessory" and row["provirus_type"] == "type1" else (
        "rec" if feature == "accessory" else feature
    )
    value = row[field]
    callable_ = value not in {"", "NA"}
    return callable_, callable_ and value in COMPATIBLE


def build_figure_1(rows, roster) -> Path:
    summary, observations, sex_evidence = structural_summary(rows, roster)
    write_structural_tables(summary, observations, sex_evidence, SUPPLEMENT)
    selected = sorted(
        (row for row in summary if row["label_scope"] == "physical_locus"),
        key=lambda row: row["variability_score"], reverse=True,
    )[:24][::-1]
    with (OUT / "Figure_1_structural_source_data.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]) + ["shown_in_figure_1"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        selected_loci = {row["locus"] for row in selected}
        writer.writerows({**row, "shown_in_figure_1": int(row["locus"] in selected_loci)} for row in summary)

    apply_style()
    fig = plt.figure(figsize=(7.1, 8.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.45, 1.35])
    ax_a = fig.add_subplot(gs[0])
    ax_a.imshow(crop_white(BIORENDER, pad=8))
    ax_a.axis("off")
    finding_title(ax_a, "A", "Long reads resolve tandem HML-2 arrays")
    ax = fig.add_subplot(gs[1])
    left = np.zeros(len(selected))
    y = np.arange(len(selected))
    for state, color in STATE_COLORS.items():
        if state == "Unknown":
            continue
        values = np.array([row[state] / row["known_haplotypes"] for row in selected])
        ax.barh(y, values, left=left, color=color, label=state, height=0.78)
        left += values
    if not np.allclose(left, 1):
        raise ValueError("Structural-state fractions do not sum to one over recovered calls")
    for index, row in enumerate(selected):
        ax.text(1.015, index, f'{row["known_haplotypes"]}/{row["eligible_haplotypes"]}',
                transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=8.5, clip_on=False)
    ax.text(1.015, 1.015, "Calls / total", transform=ax.transAxes, ha="left", fontsize=8.5)
    ax.set_yticks(y, [display_locus_name(row["locus"]) for row in selected])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of recovered locus-level calls")
    ax.set_title(
        "Structural states among recovered haplotypes",
        loc="center",
        fontsize=10.5,
        fontweight="normal",
        y=1.14,
    )
    ax.text(
        -0.105, 1.14, "B", transform=ax.transAxes, fontsize=11.5,
        fontweight="bold", color=INK, va="bottom",
    )
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.015))
    finish_axis(ax, grid="x")
    path = OUT / "Figure_1_long_read_structural_pangenome.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def load_array_counts(rows, roster):
    seven = Counter(int(row["array_copy_number"]) for row in read_tsv(SEVENP22))
    one_authority = {
        (row["sample"], row["haplotype"]): int(row["copy_number"])
        for row in read_tsv(ONEP31)
    }
    one_cells: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["Locus"] == "HML-2_1p31.1b":
            one_cells[(row["ID"], row["Haplotype"])].append(row)
    one = Counter()
    for key in roster:
        if key in one_authority:
            cn = one_authority[key]
        else:
            present = [
                row for row in one_cells.get(key, [])
                if row["observation_state"] == "PRESENT"
            ]
            cn = 1 if any(row["Structure"] == "Fragment" for row in present) else 0
        one[cn] += 1
    if sum(seven.values()) != 584 or sum(one.values()) != 584:
        raise ValueError("array authorities do not span 584 haplotypes")
    return seven, one


def build_figure_2(rows, roster) -> Path:
    seven, one = load_array_counts(rows, roster)
    cells: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        cells[(row["Locus"], row["ID"], row["Haplotype"])].append(row)
    array_loci = sorted({
        row["Locus"]
        for row in rows
        if row["cnv_support_class"] == "AUTHENTICATED_ARRAY_MULTICOPY"
    })
    per_locus: dict[str, list[int]] = {}
    for locus in array_loci:
        counts = []
        for sample, haplotype in roster:
            observed = [
                row for row in cells.get((locus, sample, haplotype), [])
                if row["observation_state"] == "PRESENT"
            ]
            if not observed:
                counts.append(0)
                continue
            counts.append(max(int(float(row["expected_biological_copy_count"])) for row in observed))
        per_locus[locus] = counts
    # Use the two locus-specific authorities for the headline arrays.
    per_locus["HML-2_7p22.1"] = [
        cn for cn, n in sorted(seven.items()) for _ in range(n)
    ]
    per_locus["HML-2_1p31.1b"] = [
        cn for cn, n in sorted(one.items()) for _ in range(n)
    ]
    order = sorted(
        array_loci,
        key=lambda locus: (
            -max(per_locus[locus]),
            -sum(cn > 1 for cn in per_locus[locus]),
            locus,
        ),
    )
    array_units: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if (
            row["observation_state"] != "PRESENT"
            or row["cnv_support_class"] != "AUTHENTICATED_ARRAY_MULTICOPY"
            or "_part" not in row["ID_Full"]
        ):
            continue
        array_base = re.sub(r"_part\d+.*$", "", row["ID_Full"])
        array_units[
            (row["Locus"], row["ID"], row["Haplotype"], array_base)
        ].append(row)
    resolved_arrays = [
        units for units in array_units.values() if len(units) >= 2
    ]
    if len(resolved_arrays) != 290:
        raise ValueError(
            f"within-array comparison denominator changed: {len(resolved_arrays)}"
        )
    features = ("gag", "pro", "pol", "env")
    arrays_by_locus: dict[str, list[list[dict[str, str]]]] = defaultdict(list)
    for units in resolved_arrays:
        arrays_by_locus[units[0]["Locus"]].append(units)
    variation_classes = {
        locus: {
            feature: Counter(
                (
                    "orf_state"
                    if len({row[feature] for row in units}) > 1
                    else "sequence_only"
                )
                for units in units_by_locus
                if len(
                    {
                        (row[feature], row[f"missense_{feature}"])
                        for row in units
                    }
                )
                > 1
            )
            for feature in features
        }
        for locus, units_by_locus in arrays_by_locus.items()
    }
    variation_counts = {
        locus: {
            feature: sum(variation_classes[locus][feature].values())
            for feature in features
        }
        for locus in arrays_by_locus
    }
    variation_table = OUT / "Figure_2_within_array_variation_source_data.tsv"
    with variation_table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "locus",
                "region",
                "multi_copy_haplotypes",
                "haplotypes_with_multiple_profiles",
                "arrays_with_orf_state_difference",
                "arrays_with_sequence_difference_only",
                "fraction_with_multiple_profiles",
            ]
        )
        for locus in sorted(arrays_by_locus):
            denominator = len(arrays_by_locus[locus])
            for feature in features:
                numerator = variation_counts[locus][feature]
                writer.writerow(
                    [
                        locus,
                        feature,
                        denominator,
                        numerator,
                        variation_classes[locus][feature]["orf_state"],
                        variation_classes[locus][feature]["sequence_only"],
                        numerator / denominator,
                    ]
                )

    apply_style()
    fig = plt.figure(figsize=(7.1, 5.2))
    ax_a = fig.add_axes([0.075, 0.635, 0.515, 0.285])
    rng = np.random.default_rng(19)
    for index, locus in enumerate(order):
        carrier_cn = np.array([cn for cn in per_locus[locus] if cn > 0])
        x = index + rng.uniform(-0.18, 0.18, size=len(carrier_cn))
        color = TYPE_I if locus.endswith("7p22.1") else (
            GREEN if locus.endswith("1p31.1b") else BLUE
        )
        ax_a.scatter(
            x, carrier_cn, s=8, color=color, alpha=0.28, linewidths=0
        )
        if len(carrier_cn):
            ax_a.plot(
                [index - 0.25, index + 0.25],
                [np.median(carrier_cn)] * 2,
                color=INK,
                lw=1.2,
            )
    ax_a.set_ylim(0.70, 6.30)
    ax_a.set_xticks(
        np.arange(len(order)),
        [display_locus_name(locus) for locus in order],
        rotation=42,
        ha="right",
    )
    ax_a.set_ylabel("Copies")
    finding_title(ax_a, "A", "Array copy number per haplotype")
    finish_axis(ax_a, grid="y")

    ax_b_high = fig.add_axes([0.075, 0.405, 0.515, 0.075])
    ax_b = fig.add_axes([0.075, 0.145, 0.515, 0.225], sharex=ax_b_high)
    x = np.arange(len(order))
    bottom = np.zeros(len(order))
    palette = {2: SKY, 3: GREEN, 4: GOLD, 5: ORANGE, 6: TYPE_II}
    for cn in range(2, 7):
        values = np.array([
            sum(value == cn for value in per_locus[locus]) / len(roster)
            for locus in order
        ])
        for axis in (ax_b_high, ax_b):
            axis.bar(
                x,
                values,
                bottom=bottom,
                color=palette[cn],
                label=str(cn),
            )
        bottom += values
    ax_b_high.set_ylim(0.382, 0.426)
    ax_b.set_ylim(0, 0.026)
    ax_b_high.spines["bottom"].set_visible(False)
    ax_b.spines["top"].set_visible(False)
    ax_b_high.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_b_high.set_yticks([0.40])
    diagonal = 0.008
    break_style = dict(color=INK, clip_on=False, lw=0.85)
    for axis, y_values in (
        (ax_b_high, (-diagonal, +diagonal)),
        (ax_b, (1 - diagonal, 1 + diagonal)),
    ):
        axis.plot(
            (-diagonal, +diagonal),
            y_values,
            transform=axis.transAxes,
            **break_style,
        )
    ax_b.set_xticks(
        x,
        [display_locus_name(locus) for locus in order],
        rotation=42,
        ha="right",
    )
    ax_b.set_ylabel("Frequency")
    ax_b_high.legend(title="Copies", ncol=5, loc="upper right")
    finding_title(ax_b_high, "B", "Expanded haplotype frequency")
    finish_axis(ax_b_high, grid="y")
    finish_axis(ax_b, grid="y")

    variation_positions = [
        [0.700, 0.700, 0.275, 0.130],
        [0.700, 0.500, 0.275, 0.130],
        [0.700, 0.300, 0.275, 0.130],
        [0.700, 0.100, 0.275, 0.130],
    ]
    feature_titles = {
        "gag": "Gag",
        "pro": "Pro",
        "pol": "Pol",
        "env": "Env",
    }
    variation_colors = {
        "orf_state": "#5D7684",
        "sequence_only": "#A48A63",
    }
    for feature_index, feature in enumerate(features):
        axis = fig.add_axes(variation_positions[feature_index])
        loci_with_variation = sorted(
            (
                locus
                for locus in arrays_by_locus
                if variation_counts[locus][feature] > 0
            ),
            key=lambda locus: (
                -len(arrays_by_locus[locus]),
                -variation_counts[locus][feature]
                / len(arrays_by_locus[locus]),
                locus,
            ),
        )
        y = np.arange(len(loci_with_variation))
        denominators = np.array(
            [len(arrays_by_locus[locus]) for locus in loci_with_variation]
        )
        state_fractions = np.array(
            [
                variation_classes[locus][feature]["orf_state"]
                / len(arrays_by_locus[locus])
                for locus in loci_with_variation
            ]
        )
        sequence_fractions = np.array(
            [
                variation_classes[locus][feature]["sequence_only"]
                / len(arrays_by_locus[locus])
                for locus in loci_with_variation
            ]
        )
        fractions = state_fractions + sequence_fractions
        axis.barh(
            y,
            state_fractions,
            color=variation_colors["orf_state"],
            height=0.62,
            zorder=2,
        )
        axis.barh(
            y,
            sequence_fractions,
            left=state_fractions,
            color=variation_colors["sequence_only"],
            height=0.62,
            zorder=2,
        )
        for row_index, (locus, fraction, denominator) in enumerate(
            zip(loci_with_variation, fractions, denominators)
        ):
            numerator = variation_counts[locus][feature]
            label_to_left = fraction >= 0.72
            axis.text(
                fraction - 0.025 if label_to_left else fraction + 0.025,
                row_index,
                f"{numerator}/{denominator}",
                ha="right" if label_to_left else "left",
                va="center",
                fontsize=7.2,
                color=INK,
                bbox={
                    "boxstyle": "square,pad=0.03",
                    "fc": "white",
                    "ec": "none",
                    "alpha": 0.90,
                },
            )
        axis.set_xlim(0, 1.08)
        axis.set_ylim(-0.55, len(loci_with_variation) - 0.45)
        axis.invert_yaxis()
        axis.set_yticks(
            y,
            [
                display_locus_name(locus)
                for locus in loci_with_variation
            ],
            fontsize=7.0,
        )
        axis.set_xticks([0, 0.5, 1.0], ["0", "50", "100"])
        if feature_index < 3:
            axis.tick_params(axis="x", labelbottom=False)
        finding_title(
            axis,
            chr(ord("C") + feature_index),
            feature_titles[feature],
            letter_x=-0.16,
        )
        finish_axis(axis, grid="x")
    legend_handles = [
        patches.Patch(
            facecolor=variation_colors["orf_state"],
            label="ΔORF",
        ),
        patches.Patch(
            facecolor=variation_colors["sequence_only"],
            label="Mutation only",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.838, 0.980),
        frameon=False,
        fontsize=7.2,
        ncol=2,
        handlelength=1.2,
        handletextpad=0.5,
        labelspacing=0.25,
    )
    fig.text(
        0.838,
        0.025,
        "Arrays (%)",
        ha="center",
        va="center",
        fontsize=9.5,
        color=INK,
    )
    path = OUT / "Figure_2_population_variable_arrays.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def crop_white(path: Path, pad: int = 18) -> Image.Image:
    image = Image.open(path).convert("RGB")
    diff = ImageChops.difference(image, Image.new("RGB", image.size, "white")).convert("L")
    box = diff.point(lambda value: 255 if value > 12 else 0).getbbox()
    if box is None:
        return image
    left, top, right, bottom = box
    return image.crop((
        max(0, left - pad), max(0, top - pad),
        min(image.width, right + pad), min(image.height, bottom + pad),
    ))


def fit_panel(source: Path, size: tuple[int, int]) -> Image.Image:
    image = crop_white(source)
    image.thumbnail(size, Image.Resampling.LANCZOS)
    panel = Image.new("RGB", size, "white")
    panel.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return panel


def load_subfamily_authority() -> dict[str, str]:
    # The original tree build persisted a consensus-anchored assignment for
    # every locus in that topology. Published or directly sequence-audited
    # assignments override those tree-local calls.
    subfamilies = {}
    local_map = PHY / "locus_subfamily.tsv"
    if local_map.exists():
        for line in local_map.read_text().splitlines():
            if not line.strip():
                continue
            locus, subfamily = line.split("\t")
            subfamilies[locus] = subfamily
    for row in read_tsv(SUBFAMILY_AUTHORITY):
        subfamily = row["final_subfamily"]
        if subfamily == "LTR5Hs":
            subfamily = "LTR5_Hs"
        if subfamily and subfamily != "not_applicable":
            subfamilies[row["locus"]] = subfamily
    subfamilies.update(TREE_SUBFAMILY_FALLBACKS)
    expanded_audit = PHY / "expanded_ltr_subfamily_assignments.tsv"
    if expanded_audit.exists():
        for row in read_tsv(expanded_audit):
            subfamilies.setdefault(row["locus"], row["subfamily"])
    return subfamilies


def kcon_ltr_sequence(sequence: str) -> str | None:
    """Use the 3' LTR when callable, otherwise the 5' LTR."""
    sequence = sequence.upper()
    five_prime = sequence[:968]
    three_prime = sequence[8504:9472] if len(sequence) >= 9472 else ""

    def canonical_count(value: str) -> int:
        return sum(base in "ACGT" for base in value)

    if canonical_count(three_prime) >= 600:
        return three_prime
    if canonical_count(five_prime) >= 600:
        return five_prime
    return None


def aligned_p_distance(first: str, second: str) -> float:
    valid = [
        (left, right)
        for left, right in zip(first, second)
        if left in "ACGT" and right in "ACGT"
    ]
    if not valid:
        return 1.0
    return sum(left != right for left, right in valid) / len(valid)


def build_expanded_ltr_tree(
    rows: list[dict[str, str]],
) -> tuple[str, Counter]:
    """Build a representative LTR topology across every callable locus.

    Up to the two most frequent exact 3'/solo-LTR sequence clusters are
    retained per locus. The full artifact-filtered catalog supplies the
    admissible sample-haplotype set and the label denominators.
    """
    allowed = {
        (
            row["Locus"].removeprefix("HML-2_"),
            row["ID"],
            row["Haplotype"],
        )
        for row in rows
        if row["observation_state"] == "PRESENT"
    }

    def normalized_haplotype(value: str) -> str:
        return {"hap1": "h1", "hap2": "h2"}.get(value.lower(), value.lower())

    per_locus: dict[str, Counter] = {}
    for path in sorted((PROJECT / "data").glob("*_carrier_kcon.aln.fa")):
        locus = path.name.removesuffix("_carrier_kcon.aln.fa")
        sequences = Counter()
        for record in SeqIO.parse(path, "fasta"):
            match = re.fullmatch(
                r"(.+)_(hap1|hap2|h1|h2|pat|mat)",
                record.id,
            )
            if not match:
                continue
            sample = match.group(1)
            haplotype = normalized_haplotype(match.group(2))
            if (locus, sample, haplotype) not in allowed:
                continue
            ltr = kcon_ltr_sequence(str(record.seq))
            if ltr is not None:
                sequences[ltr] += 1
        if sequences:
            per_locus[locus] = sequences

    if len(per_locus) < 70:
        raise ValueError(
            f"expanded LTR panel unexpectedly contains only {len(per_locus)} loci"
        )

    known_subfamilies = load_subfamily_authority()
    modal_sequences = {
        locus: sequences.most_common(1)[0][0]
        for locus, sequences in per_locus.items()
    }
    assignment_rows = []
    for locus in sorted(per_locus):
        if locus in known_subfamilies:
            assignment_rows.append(
                (locus, known_subfamilies[locus], "existing_authority", "", "")
            )
            continue
        candidates = [
            (
                aligned_p_distance(
                    modal_sequences[locus],
                    modal_sequences[other],
                ),
                other,
                known_subfamilies[other],
            )
            for other in modal_sequences
            if other in known_subfamilies
        ]
        distance, nearest_locus, subfamily = min(candidates)
        known_subfamilies[locus] = subfamily
        assignment_rows.append(
            (
                locus,
                subfamily,
                "nearest_labeled_resident_ltr",
                nearest_locus,
                f"{distance:.6f}",
            )
        )

    audit_path = PHY / "expanded_ltr_subfamily_assignments.tsv"
    with audit_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "locus",
                "subfamily",
                "assignment_source",
                "nearest_labeled_locus",
                "p_distance",
            ]
        )
        writer.writerows(assignment_rows)

    names = []
    sequences = []
    locus_counts = Counter()
    for locus in sorted(per_locus):
        locus_counts[locus] = sum(per_locus[locus].values())
        for cluster_index, (sequence, count) in enumerate(
            per_locus[locus].most_common(2),
            start=1,
        ):
            names.append(f"{locus}__hap{cluster_index}__n{count}")
            sequences.append(sequence)

    array = np.array([list(sequence) for sequence in sequences])
    canonical = np.isin(array, list("ACGT"))
    distance_rows = []
    for row_index in range(len(names)):
        row = []
        for column_index in range(row_index + 1):
            valid = canonical[row_index] & canonical[column_index]
            if not np.any(valid):
                distance = 1.0
            else:
                distance = float(
                    np.mean(array[row_index, valid] != array[column_index, valid])
                )
            row.append(distance)
        distance_rows.append(row)
    tree = DistanceTreeConstructor().nj(
        DistanceMatrix(names=names, matrix=distance_rows)
    )
    for clade in tree.find_clades():
        if clade.branch_length is not None and clade.branch_length < 0:
            clade.branch_length = 0.0
    tree_name = "hml2_pan_ltr_expanded_tree.nwk"
    Phylo.write(tree, PHY / tree_name, "newick")
    return tree_name, locus_counts


def representative_tip_names(tree) -> set[str]:
    """Return the most frequent haplotype cluster for each named locus."""
    representatives: dict[str, tuple[int, str]] = {}
    for terminal in tree.get_terminals():
        locus = terminal.name.split("__", 1)[0]
        count_match = re.search(r"__n(\d+)$", terminal.name)
        count = int(count_match.group(1)) if count_match else 1
        if locus not in representatives or count > representatives[locus][0]:
            representatives[locus] = (count, terminal.name)
    return {name for _, name in representatives.values()}


def observation_count_from_tip(name: str) -> int:
    count_match = re.search(r"__n(\d+)$", name)
    return int(count_match.group(1)) if count_match else 1


def locus_observation_counts(tree) -> Counter:
    """Sum every collapsed tip count at each locus."""
    counts = Counter()
    for terminal in tree.get_terminals():
        locus = terminal.name.split("__", 1)[0]
        counts[locus] += observation_count_from_tip(terminal.name)
    return counts


def build_representative_phylogeny(
    tree_name: str,
    title: str,
    output_name: str,
    *,
    mark_ltr_clades: bool = False,
    count_authority: Counter | None = None,
    label_loci: set[str] | None = None,
    figsize: tuple[float, float] = (3.35, 7.8),
    label_fontsize: float = 6.0,
    title_fontsize: float = 11.5,
    axis_fontsize: float = 9.5,
    clade_label_fontsize: float | None = None,
    star_size: float = 105,
) -> Path:
    """Render the full topology while labeling one cluster per locus."""
    subfamilies = load_subfamily_authority()
    colors = {
        "LTR5_Hs": TREE_LTR5HS,
        "LTR5A": TREE_LTR5A,
        "LTR5B": TREE_LTR5B,
        "LTR5": PURPLE,
    }
    clade_fontsize = (
        max(6.8, label_fontsize - 0.7)
        if clade_label_fontsize is None
        else clade_label_fontsize
    )
    tree = Phylo.read(PHY / tree_name, "newick")
    # Consensus sequences are annotation references, not observations. Remove
    # them from the topology and annotate subfamily clades separately.
    for terminal in list(tree.get_terminals()):
        if terminal.name.endswith("(ref)"):
            tree.prune(terminal)
    representatives = representative_tip_names(tree)
    locus_counts = locus_observation_counts(tree)
    missing_subfamilies = sorted(
        {
            terminal.name.split("__", 1)[0]
            for terminal in tree.get_terminals()
            if terminal.name.split("__", 1)[0] not in subfamilies
        }
    )
    if missing_subfamilies:
        raise ValueError(
            "missing LTR-subfamily assignments for displayed loci: "
            + ", ".join(missing_subfamilies)
        )
    terminals = tree.get_terminals()
    terminal_by_name = {terminal.name: terminal for terminal in terminals}

    maxheight = tree.count_terminals()
    y_positions = {
        tip: maxheight - index
        for index, tip in enumerate(reversed(terminals))
    }

    def assign_internal_y(clade):
        for child in clade:
            if child not in y_positions:
                assign_internal_y(child)
        y_positions[clade] = (
            y_positions[clade.clades[0]] + y_positions[clade.clades[-1]]
        ) / 2

    if tree.root.clades:
        assign_internal_y(tree.root)
    x_positions = tree.depths()
    if not max(x_positions.values()):
        x_positions = tree.depths(unit_branch_lengths=True)

    for clade in tree.find_clades():
        clade.color = "#A8ADB0"
        clade.width = 0.55

    apply_style()
    fig, ax = plt.subplots(figsize=figsize)
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda _clade: None,
    )

    labels = []
    for name in representatives:
        terminal = terminal_by_name[name]
        locus = name.split("__", 1)[0]
        if label_loci is not None and locus not in label_loci:
            continue
        display_count = locus_counts[locus]
        if count_authority and count_authority.get(locus, 0) > 0:
            display_count = count_authority[locus]
        labels.append(
            (
                y_positions[terminal],
                x_positions[terminal],
                f"{display_locus_name(locus)} (n={display_count})",
                colors.get(subfamilies.get(locus, ""), GRAY),
            )
        )
    labels.sort(key=lambda item: item[0])
    x_tree = max(x_positions.values())
    x_label = x_tree * 1.07
    label_y = np.linspace(maxheight * 0.025, maxheight * 0.975, len(labels))
    for (tip_y, tip_x, label, color), target_y in zip(labels, label_y):
        ax.plot(
            [tip_x, x_label * 0.985],
            [tip_y, target_y],
            color="#D4D6D8",
            lw=0.35,
            zorder=0,
        )
        ax.text(
            x_label,
            target_y,
            label,
            fontsize=label_fontsize,
            color=color,
            va="center",
        )

    if mark_ltr_clades:
        for subfamily in ("LTR5_Hs", "LTR5A", "LTR5B", "LTR5"):
            candidates = []
            for clade in tree.find_clades(order="postorder"):
                descendants = clade.get_terminals()
                assigned = [
                    subfamilies.get(tip.name.split("__", 1)[0], "")
                    for tip in descendants
                ]
                target_count = sum(value == subfamily for value in assigned)
                if target_count < 2:
                    continue
                resolved = [
                    value
                    for value in assigned
                    if value in {"LTR5_Hs", "LTR5A", "LTR5B", "LTR5"}
                ]
                if not resolved:
                    continue
                purity = target_count / len(resolved)
                threshold = 0.98 if subfamily == "LTR5_Hs" else 0.90
                if purity >= threshold:
                    candidates.append(
                        (
                            target_count,
                            purity,
                            x_positions[clade],
                            frozenset(descendants),
                            clade,
                        )
                    )
            if not candidates:
                continue
            candidates.sort(key=lambda item: item[:3], reverse=True)
            selected = []
            for candidate in candidates:
                descendant_set = candidate[3]
                if any(descendant_set & prior[3] for prior in selected):
                    continue
                selected.append(candidate)
                if subfamily != "LTR5_Hs":
                    break
            if subfamily == "LTR5_Hs":
                # The Hs subfamily can occupy several distinct portions of the
                # topology. Merge only immediately adjacent pure spans so the
                # brackets identify each visible group without redundant,
                # stacked labels.
                spans = sorted(
                    (
                        min(y_positions[tip] for tip in clade.get_terminals()),
                        max(y_positions[tip] for tip in clade.get_terminals()),
                        clade,
                    )
                    for _, _, _, _, clade in selected
                )
                merged = []
                merge_gap = max(1.5, maxheight * 0.015)
                for span_low, span_high, clade in spans:
                    if merged and span_low <= merged[-1][1] + merge_gap:
                        prior_low, prior_high, prior_clades = merged[-1]
                        merged[-1] = (
                            prior_low,
                            max(prior_high, span_high),
                            prior_clades + [clade],
                        )
                    else:
                        merged.append((span_low, span_high, [clade]))
                selected_for_drawing = sorted(
                    [
                        (
                            max(
                                clades,
                                key=lambda item: len(item.get_terminals()),
                            ),
                            span_low,
                            span_high,
                        )
                        for span_low, span_high, clades in merged
                    ],
                    key=lambda item: item[2] - item[1],
                    reverse=True,
                )[:2]
                label_spans = selected_for_drawing
            else:
                selected_for_drawing = [
                    (candidate[4], None, None) for candidate in selected
                ]
                label_spans = []
            for clade, merged_low, merged_high in selected_for_drawing:
                star_x = x_positions[clade]
                star_y = y_positions[clade]
                descendant_y = [y_positions[tip] for tip in clade.get_terminals()]
                if subfamily != "LTR5_Hs":
                    ax.scatter(
                        star_x,
                        star_y,
                        marker="*",
                        s=star_size,
                        facecolor=colors[subfamily],
                        edgecolor=INK,
                        linewidth=0.65,
                        zorder=5,
                    )
                    ax.annotate(
                        subfamily,
                        (star_x, star_y),
                        xytext=(4, 5),
                        textcoords="offset points",
                        fontsize=clade_fontsize,
                        fontfamily="Arial",
                        fontweight="bold",
                        color=INK,
                    )
                    continue
                span_low = (
                    merged_low if merged_low is not None else min(descendant_y)
                )
                span_high = (
                    merged_high if merged_high is not None else max(descendant_y)
                )
                bracket_x = -0.030 * x_tree
                tick_x = -0.010 * x_tree
                span_midpoint = (span_low + span_high) / 2
                ax.plot(
                    [bracket_x, bracket_x],
                    [span_low, span_high],
                    color=INK,
                    lw=1.5,
                    clip_on=False,
                    zorder=6,
                )
                ax.plot(
                    [bracket_x, tick_x],
                    [span_low, span_low],
                    color=INK,
                    lw=1.5,
                    clip_on=False,
                    zorder=6,
                )
                ax.plot(
                    [bracket_x, tick_x],
                    [span_high, span_high],
                    color=INK,
                    lw=1.5,
                    clip_on=False,
                    zorder=6,
                )
                ax.scatter(
                    0.005 * x_tree,
                    span_midpoint,
                    marker="*",
                    s=star_size,
                    facecolor=colors[subfamily],
                    edgecolor=INK,
                    linewidth=0.65,
                    clip_on=False,
                    zorder=7,
                )
                if any(
                    clade is label_clade
                    and merged_low == label_low
                    and merged_high == label_high
                    for label_clade, label_low, label_high in label_spans
                ):
                    ax.text(
                        -0.050 * x_tree,
                        span_midpoint,
                        "LTR5Hs",
                        rotation=90,
                        ha="center",
                        va="center",
                        fontsize=clade_fontsize,
                        fontfamily="Arial",
                        fontweight="bold",
                        color=INK,
                        clip_on=False,
                    )

    ax.set_title(title, fontsize=title_fontsize, fontweight="normal", pad=8)
    ax.set_xlabel("Substitutions per site", fontsize=axis_fontsize)
    ax.tick_params(axis="x", labelsize=max(8.0, axis_fontsize - 1.2))
    ax.set_ylabel("")
    ax.set_xlim(-0.065 * x_tree, x_tree * 1.45)
    ax.set_ylim(0, maxheight + 1)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    fig.tight_layout(pad=0.35)
    path = PHY / output_name
    fig.savefig(path, dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def build_figure_3(rows, recombination_path: Path) -> Path:
    """Show fully labeled LTR/pol topologies and acrocentric exchange."""
    expanded_ltr_tree, ltr_counts = build_expanded_ltr_tree(rows)
    gag_counts = Counter()
    pro_counts = Counter()
    pol_counts = Counter()
    env_counts = Counter()
    missing_pol = {
        "", "NA", "N", "Protein_Missing", "ComparisonError_ProteinMissing",
    }
    for row in proviral_rows(rows):
        if row["gag"] not in missing_pol:
            gag_counts[row["Locus"].removeprefix("HML-2_")] += 1
        if row["pro"] not in missing_pol:
            pro_counts[row["Locus"].removeprefix("HML-2_")] += 1
        if row["pol"] not in missing_pol:
            pol_counts[row["Locus"].removeprefix("HML-2_")] += 1
        if row["env"] not in missing_pol:
            env_counts[row["Locus"].removeprefix("HML-2_")] += 1
    build_representative_phylogeny(
        "hml2_pan_ltr_subfamily_tree.nwk",
        "LTR",
        "Supplementary_full_ltr_phylogeny.png",
        mark_ltr_clades=True,
        figsize=(7.1, 12.5),
        label_fontsize=7.5,
    )
    build_representative_phylogeny(
        "hml2_pan_orf_gag_tree.nwk",
        "Gag",
        "Supplementary_full_gag_phylogeny.png",
        mark_ltr_clades=True,
        count_authority=gag_counts,
        figsize=(7.1, 12.5),
        label_fontsize=7.5,
    )
    build_representative_phylogeny(
        "hml2_pan_orf_pro_tree.nwk",
        "Pro",
        "Supplementary_full_pro_phylogeny.png",
        mark_ltr_clades=True,
        count_authority=pro_counts,
        figsize=(7.1, 12.5),
        label_fontsize=7.5,
    )
    build_representative_phylogeny(
        "hml2_pan_orf_pol_tree.nwk",
        "Pol",
        "Supplementary_full_pol_phylogeny.png",
        mark_ltr_clades=True,
        count_authority=pol_counts,
        figsize=(7.1, 12.5),
        label_fontsize=7.5,
    )
    build_representative_phylogeny(
        "hml2_pan_orf_env_tree.nwk",
        "Env",
        "Supplementary_full_env_phylogeny.png",
        mark_ltr_clades=True,
        count_authority=env_counts,
        figsize=(7.1, 16.0),
        label_fontsize=7.5,
    )
    build_representative_phylogeny(
        expanded_ltr_tree,
        "LTR",
        "Supplementary_expanded_ltr_phylogeny.png",
        mark_ltr_clades=True,
        count_authority=ltr_counts,
        figsize=(7.1, 16.0),
        label_fontsize=7.2,
    )
    ltr_path = build_representative_phylogeny(
        expanded_ltr_tree,
        "LTR",
        "hml2_pan_ltr_representative_phylogeny.png",
        mark_ltr_clades=True,
        count_authority=ltr_counts,
        figsize=(5.0, 10.2),
        # Eighty-two locus labels share this panel. At 11 pt the rendered
        # labels overlap after the tree is fitted into the composite figure.
        label_fontsize=9.2,
        title_fontsize=16.0,
        axis_fontsize=13.5,
        clade_label_fontsize=12.0,
        star_size=310,
    )
    pol_path = build_representative_phylogeny(
        "hml2_pan_orf_pol_tree.nwk",
        "Pol",
        "hml2_pan_orf_pol_representative_phylogeny.png",
        mark_ltr_clades=True,
        count_authority=pol_counts,
        figsize=(5.0, 10.2),
        label_fontsize=11.0,
        title_fontsize=16.0,
        axis_fontsize=13.5,
        clade_label_fontsize=12.0,
        star_size=310,
    )
    acrocentric_path = draw_duplicated_groups(OUT / "Duplicated_HML2_groups.png")
    chromosome_path = build_chromosome_location_schematic()
    # Keep the full-width composite short enough to fit a portrait manuscript
    # page without Word clipping the lower panels.
    canvas = Image.new("RGB", (4200, 5500), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 96)
        legend_font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 72)
    except OSError:
        font = ImageFont.load_default()
        legend_font = font
    panels = [
        ("a", ltr_path, (110, 100, 1900, 3650)),
        ("b", pol_path, (2300, 100, 4090, 3650)),
        ("c", acrocentric_path, (80, 3870, 2290, 5440)),
        ("d", chromosome_path, (2350, 3870, 4140, 5440)),
    ]
    for letter, source, (left, top, right, bottom) in panels:
        letter_y = {"a": 0, "b": 0, "c": 3870, "d": 3870}[letter]
        letter_x = {"a": 25, "b": 2160, "c": 15, "d": 2260}[letter]
        draw.text((letter_x, letter_y), letter, font=font, fill=INK)
        panel = fit_panel(source, (right - left, bottom - top))
        canvas.paste(panel, (left, top))
    legend_items = [
        ("LTR5Hs", TREE_LTR5HS),
        ("LTR5A", TREE_LTR5A),
        ("LTR5B", TREE_LTR5B),
    ]
    try:
        legend_font = ImageFont.truetype(
            "/System/Library/Fonts/Helvetica.ttc",
            96,
        )
    except OSError:
        pass
    for y, (label, color) in zip((1380, 1600, 1820), legend_items):
        draw.rectangle((1920, y, 1984, y + 64), fill=color)
        draw.text((2010, y - 22), label, font=legend_font, fill=INK)
    path = OUT / "Figure_3_hml2_phylogenetic_mosaic.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(path, dpi=(450, 450), optimize=True)
    return path




CYTOBANDS_HG38 = {
    "1": [
        (2300000, "gneg"), (5300000, "gpos25"), (7100000, "gneg"),
        (9100000, "gpos25"), (12500000, "gneg"), (15900000, "gpos50"),
        (20100000, "gneg"), (23600000, "gpos25"), (27600000, "gneg"),
        (29900000, "gpos25"), (32300000, "gneg"), (34300000, "gpos25"),
        (39600000, "gneg"), (43700000, "gpos25"), (46300000, "gneg"),
        (50200000, "gpos75"), (55600000, "gneg"), (58500000, "gpos50"),
        (60800000, "gneg"), (68500000, "gpos50"), (69300000, "gneg"),
        (84400000, "gpos100"), (87900000, "gneg"), (91500000, "gpos75"),
        (94300000, "gneg"), (99300000, "gpos75"), (101800000, "gneg"),
        (106700000, "gpos100"), (111200000, "gneg"), (115500000, "gpos50"),
        (117200000, "gneg"), (120400000, "gpos50"), (121700000, "gneg"),
        (123400000, "acen"), (125100000, "acen"), (143200000, "gvar"),
        (147500000, "gneg"), (150600000, "gpos50"), (155100000, "gneg"),
        (156600000, "gpos50"), (159100000, "gneg"), (160500000, "gpos50"),
        (165500000, "gneg"), (167200000, "gpos50"), (170900000, "gneg"),
        (173000000, "gpos75"), (176100000, "gneg"), (180300000, "gpos50"),
        (185800000, "gneg"), (190800000, "gpos100"), (193800000, "gneg"),
        (198700000, "gpos100"), (207100000, "gneg"), (211300000, "gpos25"),
        (214400000, "gneg"), (223900000, "gpos100"), (224400000, "gneg"),
        (226800000, "gpos25"), (230500000, "gneg"), (234600000, "gpos50"),
        (236400000, "gneg"), (243500000, "gpos75"), (248956422, "gneg"),
    ],
    "4": [
        (4500000, "gneg"), (6000000, "gpos25"), (11300000, "gneg"),
        (15000000, "gpos50"), (17700000, "gneg"), (21300000, "gpos75"),
        (27700000, "gneg"), (35800000, "gpos100"), (41200000, "gneg"),
        (44600000, "gpos50"), (48200000, "gneg"), (50000000, "acen"),
        (51800000, "acen"), (58500000, "gneg"), (65500000, "gpos100"),
        (69400000, "gneg"), (75300000, "gpos75"), (78000000, "gneg"),
        (81500000, "gpos50"), (83200000, "gneg"), (86000000, "gpos25"),
        (87100000, "gneg"), (92800000, "gpos75"), (94200000, "gneg"),
        (97900000, "gpos75"), (100100000, "gneg"),
        (106700000, "gpos50"), (113200000, "gneg"),
        (119900000, "gpos75"), (122800000, "gneg"),
        (127900000, "gpos50"), (130100000, "gneg"),
        (138500000, "gpos100"), (140600000, "gneg"),
        (145900000, "gpos25"), (147500000, "gneg"),
        (150200000, "gpos25"), (154600000, "gneg"),
        (160800000, "gpos100"), (163600000, "gneg"),
        (169200000, "gpos100"), (171000000, "gneg"),
        (175400000, "gpos75"), (176600000, "gneg"),
        (182300000, "gpos100"), (186200000, "gneg"),
        (190214555, "gpos25"),
    ],
    "8": [
        (2300000, "gneg"), (6300000, "gpos75"), (12800000, "gneg"),
        (19200000, "gpos100"), (23500000, "gneg"), (27500000, "gpos50"),
        (29000000, "gneg"), (36700000, "gpos75"), (38500000, "gneg"),
        (39900000, "gpos25"), (43200000, "gneg"), (45200000, "acen"),
        (47200000, "acen"), (51300000, "gneg"), (51700000, "gpos75"),
        (54600000, "gneg"), (60600000, "gpos50"), (61300000, "gneg"),
        (65100000, "gpos50"), (67100000, "gneg"), (69600000, "gpos50"),
        (72000000, "gneg"), (74600000, "gpos100"), (74700000, "gneg"),
        (83500000, "gpos75"), (85900000, "gneg"), (92300000, "gpos100"),
        (97900000, "gneg"), (100500000, "gpos25"), (105100000, "gneg"),
        (109500000, "gpos75"), (111100000, "gneg"), (116700000, "gpos100"),
        (118300000, "gneg"), (121500000, "gpos50"), (126300000, "gneg"),
        (130400000, "gpos50"), (135400000, "gneg"), (138900000, "gpos75"),
        (145138636, "gneg"),
    ],
    "X": [
        (4400000, "gneg"), (6100000, "gpos50"), (9600000, "gneg"),
        (17400000, "gpos50"), (19200000, "gneg"), (21900000, "gpos50"),
        (24900000, "gneg"), (29300000, "gpos100"), (31500000, "gneg"),
        (37800000, "gpos100"), (42500000, "gneg"), (47600000, "gpos75"),
        (50100000, "gneg"), (54800000, "gpos25"), (58100000, "gneg"),
        (61000000, "acen"), (63800000, "acen"), (65400000, "gneg"),
        (68500000, "gpos50"), (73000000, "gneg"), (74700000, "gpos50"),
        (76800000, "gneg"), (85400000, "gpos100"), (87000000, "gneg"),
        (92700000, "gpos100"), (94300000, "gneg"), (99100000, "gpos75"),
        (103300000, "gneg"), (104500000, "gpos50"), (109400000, "gneg"),
        (117400000, "gpos75"), (121800000, "gneg"), (129500000, "gpos100"),
        (131300000, "gneg"), (134500000, "gpos25"), (138900000, "gneg"),
        (141200000, "gpos75"), (143000000, "gneg"), (148000000, "gpos100"),
        (156040895, "gneg"),
    ],
    "13": [
        (4600000, "gvar"), (10100000, "stalk"), (16500000, "gvar"),
        (17700000, "acen"), (18900000, "acen"), (22600000, "gneg"),
        (24900000, "gpos25"), (27200000, "gneg"), (28300000, "gpos25"),
        (31600000, "gneg"), (33400000, "gpos50"), (34900000, "gneg"),
        (39500000, "gpos75"), (44600000, "gneg"), (45200000, "gpos25"),
        (46700000, "gneg"), (50300000, "gpos50"), (54700000, "gneg"),
        (59000000, "gpos100"), (61800000, "gneg"), (65200000, "gpos75"),
        (68100000, "gneg"), (72800000, "gpos100"), (74900000, "gneg"),
        (76700000, "gpos50"), (78500000, "gneg"), (87100000, "gpos100"),
        (89400000, "gneg"), (94400000, "gpos100"), (97500000, "gneg"),
        (98700000, "gpos25"), (101100000, "gneg"), (104200000, "gpos100"),
        (106400000, "gneg"), (109600000, "gpos100"), (114364328, "gneg"),
    ],
    "15": [
        (4200000, "gvar"), (9700000, "stalk"), (17500000, "gvar"),
        (19000000, "acen"), (20500000, "acen"), (25500000, "gneg"),
        (27800000, "gpos50"), (30000000, "gneg"), (30900000, "gpos50"),
        (33400000, "gneg"), (39800000, "gpos75"), (42500000, "gneg"),
        (43300000, "gpos25"), (44500000, "gneg"), (49200000, "gpos75"),
        (52600000, "gneg"), (58800000, "gpos75"), (59000000, "gneg"),
        (63400000, "gpos25"), (66900000, "gneg"), (67000000, "gpos25"),
        (67200000, "gneg"), (72400000, "gpos25"), (74900000, "gneg"),
        (76300000, "gpos25"), (78000000, "gneg"), (81400000, "gpos50"),
        (84700000, "gneg"), (88500000, "gpos50"), (93800000, "gneg"),
        (98000000, "gpos50"), (101991189, "gneg"),
    ],
    "21": [
        (3100000, "gvar"), (7000000, "stalk"), (10900000, "gvar"),
        (12000000, "acen"), (13000000, "acen"), (15000000, "gneg"),
        (22600000, "gpos100"), (25500000, "gneg"), (30200000, "gpos75"),
        (34400000, "gneg"), (36400000, "gpos50"), (38300000, "gneg"),
        (41200000, "gpos50"), (46709983, "gneg"),
    ],
    "22": [
        (4300000, "gvar"), (9400000, "stalk"), (13700000, "gvar"),
        (15000000, "acen"), (17400000, "acen"), (21700000, "gneg"),
        (23100000, "gpos25"), (25500000, "gneg"), (29200000, "gpos50"),
        (31800000, "gneg"), (37200000, "gpos50"), (40600000, "gneg"),
        (43800000, "gpos50"), (48100000, "gneg"), (49100000, "gpos50"),
        (50818468, "gneg"),
    ],
}


def build_chromosome_location_schematic() -> Path:
    """Place HML-2 groups on GRCh38 cytoband ideograms."""
    apply_style()
    fig, ax = plt.subplots(figsize=(4.05, 3.55))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    band_colors = {
        "gneg": "#FFFFFF",
        "gpos25": "#D9DDE0",
        "gpos50": "#AEB5BA",
        "gpos75": "#737C82",
        "gpos100": "#30363A",
        "gvar": "#C8CDD0",
        "stalk": "#E3E6E8",
        "acen": "#B96868",
    }
    chromosome_edge = "#687279"
    ideogram_top = 0.82
    ideogram_bottom = 0.12
    ideogram_height = ideogram_top - ideogram_bottom
    chromosome_width = 0.050
    chromosome_x = {
        "1": 0.055,
        "4": 0.145,
        "8": 0.235,
        "X": 0.335,
        "13": 0.52,
        "15": 0.65,
        "21": 0.80,
        "22": 0.93,
    }

    def genomic_y(chromosome_name: str, coordinate: float) -> float:
        total = CYTOBANDS_HG38[chromosome_name][-1][0]
        return ideogram_top - coordinate / total * ideogram_height

    def chromosome(chromosome_name: str) -> None:
        x = chromosome_x[chromosome_name]
        bands = CYTOBANDS_HG38[chromosome_name]
        total = bands[-1][0]
        clip = patches.FancyBboxPatch(
            (x - chromosome_width / 2, ideogram_bottom),
            chromosome_width,
            ideogram_height,
            boxstyle="round,pad=0.001,rounding_size=0.008",
            facecolor="white",
            edgecolor="none",
            zorder=1,
        )
        ax.add_patch(clip)
        start = 0
        acen_limits = []
        for end, stain in bands:
            band_top = ideogram_top - start / total * ideogram_height
            band_bottom = ideogram_top - end / total * ideogram_height
            rectangle = patches.Rectangle(
                (x - chromosome_width / 2, band_bottom),
                chromosome_width,
                band_top - band_bottom,
                facecolor=band_colors[stain],
                edgecolor="#8D969B",
                linewidth=0.55,
                zorder=2,
            )
            rectangle.set_clip_path(clip)
            ax.add_patch(rectangle)
            if stain == "acen":
                acen_limits.append((band_bottom, band_top))
            start = end
        if acen_limits:
            centromere = np.mean(
                [coordinate for limits in acen_limits for coordinate in limits]
            )
            notch_height = 0.016
            for side in (-1, 1):
                outer_x = x + side * chromosome_width / 2
                inner_x = x + side * chromosome_width * 0.12
                ax.add_patch(
                    patches.Polygon(
                        [
                            (outer_x, centromere - notch_height),
                            (inner_x, centromere),
                            (outer_x, centromere + notch_height),
                        ],
                        closed=True,
                        facecolor="white",
                        edgecolor="none",
                        zorder=4,
                    )
                )
        ax.add_patch(
            patches.FancyBboxPatch(
                (x - chromosome_width / 2, ideogram_bottom),
                chromosome_width,
                ideogram_height,
                boxstyle="round,pad=0.001,rounding_size=0.008",
                facecolor="none",
                edgecolor=chromosome_edge,
            linewidth=1.1,
                zorder=5,
            )
        )
        ax.text(
            x,
            0.065,
            chromosome_name,
            ha="center",
            va="center",
            fontsize=9.6,
            color=INK,
        )

    def mark_group(
        chromosome_name: str,
        coordinate: float,
        count: int,
        color: str,
    ) -> list[float]:
        x = chromosome_x[chromosome_name]
        center = genomic_y(chromosome_name, coordinate)
        offsets = (np.arange(count) - (count - 1) / 2) * 0.009
        positions = [center + offset for offset in offsets]
        for y in positions:
            ax.plot(
                [x - 0.034, x + 0.034],
                [y, y],
                color=color,
                linewidth=1.9,
                solid_capstyle="butt",
                zorder=7,
            )
        if count > 1:
            bracket_x = x + 0.041
            ax.plot(
                [bracket_x, bracket_x],
                [min(positions), max(positions)],
                color=color,
                linewidth=1.0,
                zorder=7,
            )
        return positions

    def connect(
        first,
        second,
        color,
        *,
        curve=-0.15,
        linewidth=1.6,
        alpha=0.82,
        linestyle="solid",
        zorder=6,
        arrowstyle="-",
    ):
        ax.add_patch(
            patches.FancyArrowPatch(
                first,
                second,
                arrowstyle=arrowstyle,
                mutation_scale=8,
                shrinkA=2,
                shrinkB=2,
                connectionstyle=f"arc3,rad={curve}",
                color=color,
                linewidth=linewidth,
                alpha=alpha,
                linestyle=linestyle,
                zorder=zorder,
            )
        )

    ax.text(
        0.20,
        0.94,
        "Peritelomeric groups",
        ha="center",
        va="center",
        fontsize=10.3,
        color=INK,
    )
    ax.text(
        0.73,
        0.94,
        "Duplicated p-arm groups",
        ha="center",
        va="center",
        fontsize=10.3,
        color=INK,
    )
    ax.plot([0.41, 0.41], [0.10, 0.90], color=GRID, linewidth=0.9)

    for chromosome_name in ("1", "4", "8", "X", "13", "15", "21", "22"):
        chromosome(chromosome_name)

    mark_group("1", 14200000, 3, BLUE)
    four_positions = mark_group("4", 190109902, 1, TYPE_II)
    mark_group("8", 9550000, 4, ORANGE)
    mark_group("X", 152000000, 2, RED)
    ax.text(0.055, 0.865, "1p36.21", ha="center", fontsize=9.2, color=INK)
    ax.text(0.235, 0.865, "8p23.1", ha="center", fontsize=9.2, color=INK)
    ax.text(
        chromosome_x["4"] - 0.018,
        0.000,
        "4q35.2",
        ha="center",
        va="bottom",
        fontsize=8.7,
        color=INK,
    )
    ax.text(
        chromosome_x["X"] + 0.030,
        0.865,
        "Xq28",
        ha="center",
        fontsize=8.7,
        color=INK,
    )

    thirteen_positions = mark_group("13", 2300000, 1, TYPE_I)
    fifteen_positions = mark_group("15", 2100000, 2, TYPE_II)
    twentyone_positions = mark_group("21", 1550000, 1, TYPE_II)
    twentytwo_positions = mark_group("22", 2150000, 1, TYPE_II)
    ax.text(
        chromosome_x["13"],
        0.885,
        "13p13",
        ha="center",
        va="center",
        fontsize=8.7,
        color=INK,
    )
    ax.text(
        chromosome_x["15"],
        0.88,
        "15p13\na,b",
        ha="center",
        va="center",
        fontsize=8.4,
        linespacing=0.9,
        color=INK,
    )
    ax.text(
        chromosome_x["21"],
        0.885,
        "21p13",
        ha="center",
        va="center",
        fontsize=8.7,
        color=INK,
    )
    ax.text(
        chromosome_x["22"],
        0.885,
        "22p13",
        ha="center",
        va="center",
        fontsize=8.7,
        color=INK,
    )
    connect(
        (chromosome_x["13"], thirteen_positions[0] + 0.005),
        (chromosome_x["15"], fifteen_positions[-1] + 0.005),
        TYPE_I,
        curve=-0.24,
        linewidth=2.0,
    )
    connect(
        (chromosome_x["15"], fifteen_positions[0]),
        (chromosome_x["21"], twentyone_positions[0]),
        TYPE_II,
        curve=-0.20,
        linewidth=1.7,
    )
    connect(
        (chromosome_x["15"], fifteen_positions[0] - 0.003),
        (chromosome_x["22"], twentytwo_positions[0]),
        TYPE_II,
        curve=-0.26,
        linewidth=1.5,
    )
    connect(
        (chromosome_x["21"], twentyone_positions[0]),
        (chromosome_x["22"], twentytwo_positions[0]),
        TYPE_II,
        curve=-0.24,
        linewidth=1.8,
    )
    # Route the inferred 21p13-to-4q35.2 duplication around the chromosome
    # field. The arrowhead terminates at the 4q35.2 marker.
    source = (chromosome_x["21"] + 0.034, twentyone_positions[0])
    top_lane = 0.85
    right_lane = 0.985
    bottom_lane = 0.022
    approach_lane = 0.205
    target_y = four_positions[0]
    route = [
        source,
        (chromosome_x["21"] + 0.055, top_lane),
        (right_lane, top_lane),
        (right_lane, bottom_lane),
        (approach_lane, bottom_lane),
        (approach_lane, target_y),
    ]
    ax.plot(
        [point[0] for point in route],
        [point[1] for point in route],
        color=TYPE_II,
        linewidth=1.55,
        solid_capstyle="round",
        solid_joinstyle="round",
        zorder=6,
    )
    ax.add_patch(
        patches.FancyArrowPatch(
            route[-1],
            (chromosome_x["4"] + 0.034, target_y),
            arrowstyle="-|>",
            mutation_scale=9,
            shrinkA=0,
            shrinkB=0,
            color=TYPE_II,
            linewidth=1.55,
            zorder=7,
        )
    )

    ax.set_title(
        "Chromosomal locations",
        fontsize=12.0,
        fontweight="normal",
        pad=5,
    )
    path = OUT / "Chromosome_location_schematic.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def proviral_rows(rows):
    return [
        row for row in rows
        if row["observation_state"] == "PRESENT" and row["Structure"] in PROVIRUS
    ]


def pattern_tuple(row) -> tuple[bool, ...]:
    return tuple(gene_call(row, feature)[1] for feature in FEATURES)


def draw_upset(ax_bar, ax_matrix, rows, type_name: str, color: str, letter: str) -> None:
    if type_name == "type1":
        labels = ["Gag", "Gag–Pro", "Gag–Pro–Pol", "Np9"]
        feature_indices = (0, 1, 2, 4)
    else:
        labels = ["Gag", "Gag–Pro", "Gag–Pro–Pol", "Env", "Rec"]
        feature_indices = (0, 1, 2, 3, 4)
    counts = Counter(
        tuple(pattern_tuple(row)[index] for index in feature_indices)
        for row in rows
    )
    patterns = [item for item in counts.most_common(13) if any(item[0])]
    if not patterns:
        raise ValueError(f"no ORF patterns for {type_name}")
    x = np.arange(len(patterns))
    heights = [count for _, count in patterns]
    ax_bar.bar(x, heights, color=color, width=0.72)
    ax_bar.set_xlim(-0.5, len(patterns) - 0.5)
    ax_bar.set_ylim(0, max(heights) * 1.16)
    ax_bar.set_ylabel("Proviral observations")
    ax_bar.set_xticks([])
    for xi, height in zip(x, heights):
        ax_bar.text(
            xi,
            height + max(heights) * 0.025,
            f"{height:,}",
            ha="center",
            fontsize=7.0,
        )
    finding_title(
        ax_bar,
        letter,
        "Type I" if type_name == "type1" else "Type II",
    )
    finish_axis(ax_bar, grid="y")

    ax_matrix.set_xlim(-0.5, len(patterns) - 0.5)
    ax_matrix.set_ylim(-0.5, len(labels) - 0.5)
    ax_matrix.invert_yaxis()
    ax_matrix.set_yticks(np.arange(len(labels)), labels)
    ax_matrix.set_xticks([])
    for column, (pattern, _) in enumerate(patterns):
        on = [index for index, present in enumerate(pattern) if present]
        if len(on) > 1:
            ax_matrix.plot([column, column], [min(on), max(on)], color=INK, lw=1.15, zorder=1)
        for row_index, present in enumerate(pattern):
            if present:
                ax_matrix.scatter(
                    column,
                    row_index,
                    marker="s",
                    s=24,
                    color=INK,
                    linewidths=0,
                    zorder=2,
                )
    ax_matrix.spines[:].set_visible(False)
    ax_matrix.tick_params(axis="y", length=0)
    ax_matrix.set_xlabel("")


def build_figure_4(rows) -> Path:
    prov = proviral_rows(rows)
    burden_features = (
        "gag",
        "gag_pro_route",
        "gag_pro_pol_route",
        "env_type2",
        "np9_type1",
        "rec_type2",
    )
    burden_labels = (
        "Gag",
        "Gag–Pro",
        "Gag–Pro–Pol",
        "Env",
        "Np9",
        "Rec",
    )
    burden_c_labels = (
        "Gag",
        "Gag–\nPro",
        "Gag–\nPro–\nPol",
        "Env",
        "Np9",
        "Rec",
    )
    per_person: dict[str, dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    for row in rows:
        if not PUBLIC_ID.fullmatch(row["ID"]):
            continue
        weight = float(row["observation_weight"])
        biological_copies = row["expected_biological_copy_count"]
        route_copy_count = (
            0.0
            if biological_copies in {"", "NA"}
            else float(biological_copies)
        )
        values = {
            "gag": row["expected_intact_gag_copies"],
            "gag_pro_route": (
                route_copy_count
                if gene_call(row, "gag_pro_route")[1]
                else 0.0
            ),
            "gag_pro_pol_route": (
                route_copy_count
                if gene_call(row, "gag_pro_pol_route")[1]
                else 0.0
            ),
            "env_type2": (
                row["expected_intact_env_copies"]
                if row["provirus_type"] == "type2"
                else 0.0
            ),
            "np9_type1": (
                row["expected_intact_accessory_copies"]
                if row["provirus_type"] == "type1"
                else 0.0
            ),
            "rec_type2": (
                row["expected_intact_accessory_copies"]
                if row["provirus_type"] == "type2"
                else 0.0
            ),
        }
        for feature, value in values.items():
            if value not in {"", "NA"}:
                per_person[row["ID"]][feature] += weight * float(value)
    if len(per_person) != 292:
        raise ValueError(f"unexpected ORF-burden donor count: {len(per_person)}")

    sample_superpopulation = {}
    with IGSR_SAMPLES.open(newline="") as handle:
        for sample_row in csv.DictReader(handle, delimiter="\t"):
            sample_superpopulation[sample_row["Sample name"]] = (
                sample_row["Superpopulation code"]
            )
    superpopulation_order = ("AFR", "AMR", "EAS", "EUR", "SAS")
    superpopulation_colors = {
        "AFR": "#0072B2",
        "AMR": "#E69F00",
        "EAS": "#7A5195",
        "EUR": "#009E73",
        "SAS": "#CC79A7",
    }
    mapped_people = [
        sample
        for sample in per_person
        if sample_superpopulation.get(sample) in superpopulation_order
    ]
    if len(mapped_people) != 286:
        raise ValueError(
            f"unexpected mapped ORF-burden donor count: {len(mapped_people)}"
        )

    burden_path = OUT / "Figure_4_person_orf_burden.tsv"
    with burden_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["sample", "superpopulation", *burden_features]
        )
        for sample in sorted(per_person):
            writer.writerow(
                [
                    sample,
                    sample_superpopulation.get(sample, "Unknown"),
                    *[
                        f"{per_person[sample][feature]:.6g}"
                        for feature in burden_features
                    ],
                ]
            )

    apply_style()
    fig = plt.figure(figsize=(7.1, 5.35))
    type1_bar = fig.add_axes([0.10, 0.70, 0.37, 0.225])
    type1_matrix = fig.add_axes([0.10, 0.49, 0.37, 0.175])
    type2_bar = fig.add_axes([0.61, 0.70, 0.37, 0.225])
    type2_matrix = fig.add_axes([0.61, 0.49, 0.37, 0.175])
    draw_upset(
        type1_bar,
        type1_matrix,
        [row for row in prov if row["provirus_type"] == "type1"],
        "type1",
        TYPE_I,
        "A",
    )
    draw_upset(
        type2_bar,
        type2_matrix,
        [row for row in prov if row["provirus_type"] == "type2"],
        "type2",
        TYPE_II,
        "B",
    )

    ax_c = fig.add_axes([0.10, 0.09, 0.37, 0.275])
    feature_values = [
        np.array(
            [per_person[sample][feature] for sample in sorted(per_person)]
        )
        for feature in burden_features
    ]
    rng = np.random.default_rng(41)
    population_offsets = dict(
        zip(
            superpopulation_order,
            np.linspace(-0.22, 0.22, len(superpopulation_order)),
        )
    )
    for index, values in enumerate(feature_values, start=1):
        for superpopulation in superpopulation_order:
            selected_people = [
                sample
                for sample in mapped_people
                if sample_superpopulation[sample] == superpopulation
            ]
            selected_values = np.array([
                per_person[sample][burden_features[index - 1]]
                for sample in selected_people
            ])
            jitter = rng.uniform(-0.025, 0.025, len(selected_values))
            ax_c.scatter(
                (
                    np.full(len(selected_values), index)
                    + population_offsets[superpopulation]
                    + jitter
                ),
                selected_values,
                s=8,
                color=superpopulation_colors[superpopulation],
                alpha=0.55,
                linewidths=0,
                zorder=1,
            )
    means = [float(np.mean(values)) for values in feature_values]
    ax_c.scatter(
        np.arange(1, len(means) + 1),
        means,
        marker="D",
        s=30,
        color=INK,
        edgecolor="white",
        linewidth=0.7,
        zorder=4,
    )
    for index, (mean, values) in enumerate(
        zip(means, feature_values),
        start=1,
    ):
        ax_c.text(
            index,
            float(np.max(values)) + 1.4,
            f"{mean:.1f}",
            ha="center",
            va="bottom",
            fontsize=8.8,
            color=INK,
            bbox={
                "boxstyle": "square,pad=0.08",
                "fc": "white",
                "ec": "none",
                "alpha": 0.92,
            },
            zorder=5,
        )
    ax_c.set_xticks(
        np.arange(1, len(burden_c_labels) + 1),
        burden_c_labels,
    )
    ax_c.tick_params(axis="x", labelsize=9.5)
    ax_c.tick_params(axis="y", labelsize=9.5)
    ax_c.set_ylabel("Copies", fontsize=10.5)
    maximum_copy_count = max(float(np.max(values)) for values in feature_values)
    ax_c.set_ylim(0, math.ceil((maximum_copy_count + 5) / 10) * 10)
    ax_c.set_xlim(0.55, len(burden_features) + 0.45)
    finding_title(ax_c, "C", "Coding capacity per person")
    ax_c.title.set_fontsize(13.5)
    population_handles = [
        matplotlib.lines.Line2D(
            [],
            [],
            linestyle="none",
            marker="o",
            markersize=4.8,
            markerfacecolor=superpopulation_colors[superpopulation],
            markeredgewidth=0,
            label=("Unmapped" if superpopulation == "Unknown"
                   else superpopulation),
        )
        for superpopulation in superpopulation_order
    ]
    finish_axis(ax_c, grid="y")

    # Panel D asks a different question from panel C: which biologically
    # interpretable coding features, at named loci, show the largest observed
    # frequency differences among superpopulations? Pro and Pol are counted
    # only within continuous Gag–Pro and Gag–Pro–Pol reading paths, matching
    # panels A–C. Use an absolute 0--100% scale so the plot does not visually
    # magnify small differences in aggregate copy burden.
    locus_cells: dict[
        tuple[str, str, str], list[dict[str, str]]
    ] = defaultdict(list)
    for row in rows:
        if (
            PUBLIC_ID.fullmatch(row["ID"])
            and sample_superpopulation.get(row["ID"])
            in superpopulation_order
        ):
            locus_cells[(row["Locus"], row["ID"], row["Haplotype"])].append(row)

    haplotypes_by_superpopulation = {
        superpopulation: sorted(
            {
                (row["ID"], row["Haplotype"])
                for row in rows
                if sample_superpopulation.get(row["ID"]) == superpopulation
            }
        )
        for superpopulation in superpopulation_order
    }
    locus_feature_specs = (
        ("Gag", "gag", None),
        ("Gag–Pro", "gag_pro_route", None),
        ("Gag–Pro–Pol", "gag_pro_pol_route", None),
        ("Env", "env", "type2"),
        ("Np9", "accessory", "type1"),
        ("Rec", "accessory", "type2"),
    )
    locus_feature_summaries = []
    for feature_label, feature, required_type in locus_feature_specs:
        for locus in sorted({row["Locus"] for row in rows}):
            values_by_superpopulation = {}
            counts_by_superpopulation = {}
            for superpopulation in superpopulation_order:
                values = []
                for sample, haplotype in haplotypes_by_superpopulation[
                    superpopulation
                ]:
                    cell_rows = locus_cells.get(
                        (locus, sample, haplotype),
                        [],
                    )
                    if not cell_rows or all(
                        row["observation_state"] == "UNKNOWN_TECHNICAL"
                        for row in cell_rows
                    ):
                        continue
                    expected_feature_copies = 0.0
                    for row in cell_rows:
                        if (
                            required_type is not None
                            and row["provirus_type"] != required_type
                        ):
                            continue
                        if feature == "gag":
                            value = row["expected_intact_gag_copies"]
                        elif feature == "env":
                            value = row["expected_intact_env_copies"]
                        elif feature == "accessory":
                            value = row["expected_intact_accessory_copies"]
                        else:
                            biological_copies = row[
                                "expected_biological_copy_count"
                            ]
                            value = (
                                biological_copies
                                if gene_call(row, feature)[1]
                                else "0"
                            )
                        if value not in {"", "NA"}:
                            expected_feature_copies += (
                                float(row["observation_weight"])
                                * float(value)
                            )
                    # Report the expected fraction of haplotypes carrying at
                    # least one intact copy, rather than letting tandem arrays
                    # contribute multiple times.
                    values.append(min(1.0, expected_feature_copies))
                values_by_superpopulation[superpopulation] = (
                    100.0 * float(np.mean(values)) if values else math.nan
                )
                counts_by_superpopulation[superpopulation] = len(values)
            if min(counts_by_superpopulation.values()) < 50:
                continue
            frequencies = list(values_by_superpopulation.values())
            overall_frequency = float(
                np.average(
                    frequencies,
                    weights=list(counts_by_superpopulation.values()),
                )
            )
            if not 1.0 <= overall_frequency <= 99.0:
                continue
            locus_feature_summaries.append(
                {
                    "locus": display_locus_name(locus),
                    "feature": feature_label,
                    "frequencies": values_by_superpopulation,
                    "counts": counts_by_superpopulation,
                    "range": max(frequencies) - min(frequencies),
                }
            )

    # Show the locus with the largest observed range for each coding feature.
    displayed_locus_features = []
    for feature_label, _, _ in locus_feature_specs:
        candidates = [
            summary
            for summary in locus_feature_summaries
            if summary["feature"] == feature_label
        ]
        displayed_locus_features.append(
            max(candidates, key=lambda summary: summary["range"])
        )
    displayed_locus_features.sort(
        key=lambda summary: summary["range"],
        reverse=True,
    )

    locus_orf_path = OUT / "Figure_4_locus_orf_superpopulation_frequencies.tsv"
    with locus_orf_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "locus",
                "coding_feature",
                "superpopulation",
                "feature_haplotype_frequency_percent",
                "known_haplotypes",
                "max_min_difference_percentage_points",
            ]
        )
        for summary in displayed_locus_features:
            for superpopulation in superpopulation_order:
                writer.writerow(
                    [
                        summary["locus"],
                        summary["feature"],
                        superpopulation,
                        f'{summary["frequencies"][superpopulation]:.3f}',
                        summary["counts"][superpopulation],
                        f'{summary["range"]:.3f}',
                    ]
                )

    ax_d = fig.add_axes([0.62, 0.09, 0.355, 0.275])
    y_positions = np.arange(len(displayed_locus_features))[::-1]
    for y_position, summary in zip(y_positions, displayed_locus_features):
        frequencies = summary["frequencies"]
        ax_d.plot(
            [min(frequencies.values()), max(frequencies.values())],
            [y_position, y_position],
            color="#D7D5D0",
            linewidth=2.0,
            solid_capstyle="round",
            zorder=1,
        )
        for superpopulation in superpopulation_order:
            ax_d.scatter(
                frequencies[superpopulation],
                y_position,
                s=23,
                color=superpopulation_colors[superpopulation],
                edgecolor="white",
                linewidth=0.45,
                zorder=3,
            )
    ax_d.set_yticks(
        y_positions,
        [
            f'{summary["locus"]} {summary["feature"]}'
            for summary in displayed_locus_features
        ],
    )
    ax_d.set_xlim(0, 100)
    ax_d.set_xticks([0, 25, 50, 75, 100])
    ax_d.set_xlabel("Haplotypes carrying feature (%)", fontsize=8.8)
    ax_d.tick_params(axis="x", labelsize=8.0)
    ax_d.tick_params(axis="y", labelsize=8.2, length=0, pad=3)
    finding_title(
        ax_d,
        "D",
        "Locus-specific coding capacity",
        letter_x=-0.22,
    )
    ax_d.title.set_fontsize(12.3)
    finish_axis(ax_d, grid="x")

    # One figure-level legend, centered between the upper and lower rows,
    # explicitly applies the same superpopulation colors to panels C and D.
    fig.legend(
        handles=population_handles,
        loc="center",
        bbox_to_anchor=(0.54, 0.445),
        frameon=False,
        fontsize=8.0,
        ncol=5,
        columnspacing=1.15,
        handletextpad=0.35,
        borderaxespad=0,
    )

    path = OUT / "Figure_4_same_molecule_coding_capacity.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def normalized_locus_for_signature(locus: str) -> str:
    locus = locus.removeprefix("HML-2_").removesuffix("_hg38")
    if locus in {"13p13", "15p13b"}:
        return "acro_type1"
    if locus in {"15p13a", "21p13", "22p13"}:
        return "acro_type2"
    return locus


def clean_long_missense_signature(value: str) -> str | None:
    missing = {
        "", "NA", "N", "nan", "Protein_Missing",
        "ComparisonError_ProteinMissing",
    }
    if value in missing:
        return None
    value = re.sub(r"Undetermined:[^,]*", "", value)
    mutations = [
        token.strip()
        for token in value.split(",")
        if token.strip() and token.strip() not in missing
    ]
    if len(mutations) <= 15:
        return None
    return ",".join(mutations)


def build_peritelomeric_signature_figure(rows) -> Path:
    """Recompute exact long-signature sharing from the current corrected catalog."""
    observations: set[tuple[str, str, str]] = set()
    for row in proviral_rows(rows):
        locus = normalized_locus_for_signature(row["Locus"])
        for orf in ("gag", "pro", "pol", "env"):
            signature = clean_long_missense_signature(row[f"missense_{orf}"])
            if signature:
                observations.add((locus, orf, signature))

    signature_loci: dict[tuple[str, str], set[str]] = defaultdict(set)
    for locus, orf, signature in observations:
        signature_loci[(orf, signature)].add(locus)

    pair_counts = Counter()
    pair_orf_counts = Counter()
    for (orf, _signature), loci in signature_loci.items():
        if len(loci) < 2:
            continue
        for first, second in itertools.combinations(sorted(loci), 2):
            pair_counts[(first, second)] += 1
            pair_orf_counts[(first, second, orf)] += 1
    if not pair_counts:
        raise ValueError("no inter-locus long missense signatures in current catalog")

    # Preserve chromosome-arm labels for the focused Type-II comparison. The
    # general summary above intentionally collapses unresolved acrocentric
    # labels, but that would hide the observed sharing with 4q35.2.
    type2_loci = ("15p13a", "21p13", "22p13", "4q35.2_hg38")
    type2_observations: set[tuple[str, str, str]] = set()
    for row in rows:
        locus = row["Locus"].removeprefix("HML-2_")
        if locus not in type2_loci or row["observation_state"] != "PRESENT":
            continue
        for orf in ("gag", "pro", "pol", "env"):
            signature = clean_long_missense_signature(row[f"missense_{orf}"])
            if signature:
                type2_observations.add((locus, orf, signature))
    type2_signature_loci: dict[tuple[str, str], set[str]] = defaultdict(set)
    for locus, orf, signature in type2_observations:
        type2_signature_loci[(orf, signature)].add(locus)
    type2_pair_counts = Counter()
    for signature_members in type2_signature_loci.values():
        for first, second in itertools.combinations(sorted(signature_members), 2):
            type2_pair_counts[(first, second)] += 1
    if not any("4q35.2_hg38" in pair for pair in type2_pair_counts):
        raise ValueError("4q35.2 is absent from the Type-II signature comparison")

    table_path = OUT / "peritelomeric_shared_signature_pairs.tsv"
    with table_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["locus_1", "locus_2", "shared_signatures", "gag", "pro", "pol", "env"]
        )
        for (first, second), count in pair_counts.most_common():
            writer.writerow(
                [
                    first,
                    second,
                    count,
                    *[
                        pair_orf_counts[(first, second, orf)]
                        for orf in ("gag", "pro", "pol", "env")
                    ],
                ]
            )

    def draw_signature_matrix(ax, counts, loci, title):
        matrix = np.zeros((len(loci), len(loci)), dtype=float)
        for (first, second), count in counts.items():
            if first not in loci or second not in loci:
                continue
            i = loci.index(first)
            j = loci.index(second)
            matrix[i, j] = count
            matrix[j, i] = count
        display = matrix.copy()
        display[np.tril_indices_from(display)] = np.nan
        display[display == 0] = np.nan
        color_map = matplotlib.colormaps["viridis"].copy()
        color_map.set_bad("white")
        maximum = max(counts.values())
        image = ax.imshow(
            display,
            cmap=color_map,
            vmin=1,
            vmax=maximum,
            aspect="equal",
        )
        for i in range(len(loci)):
            for j in range(i + 1, len(loci)):
                value = int(matrix[i, j])
                if not value:
                    continue
                color = "white" if value >= maximum * 0.48 else INK
                ax.text(
                    j,
                    i,
                    str(value),
                    ha="center",
                    va="center",
                    color=color,
                    fontsize=9.2,
                )
        labels = [
            display_locus_name(locus)
            .replace("8p23.1 duplicate group", "8p23.1 unresolved")
            .replace("_hg38", "")
            for locus in loci
        ]
        ax.set_xticks(np.arange(len(loci)), labels, rotation=35, ha="right")
        ax.set_yticks(np.arange(len(loci)), labels)
        ax.tick_params(length=0, labelsize=8.5)
        ax.set_title(title, fontweight="normal", fontsize=11.5)
        for spine in ax.spines.values():
            spine.set_visible(False)
        colorbar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.025)
        colorbar.set_label("Shared signatures", fontsize=8.5)
        colorbar.ax.tick_params(labelsize=7.7)

    general_loci = sorted({locus for pair in pair_counts for locus in pair})
    apply_style()
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9.4, 4.8),
        gridspec_kw={"width_ratios": [2.15, 1]},
        constrained_layout=True,
    )
    draw_signature_matrix(
        axes[0],
        pair_counts,
        general_loci,
        "Peritelomeric duplication groups",
    )
    draw_signature_matrix(
        axes[1],
        type2_pair_counts,
        list(type2_loci),
        "Type-II group and 4q35.2",
    )
    path = OUT / "Supplementary_peritelomeric_shared_signatures.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def build_supplementary_orf_figures(rows) -> list[Path]:
    prov = proviral_rows(rows)
    by_type_locus: dict[str, dict[str, list[dict[str, str]]]] = {
        "type1": defaultdict(list),
        "type2": defaultdict(list),
    }
    for row in prov:
        if row["provirus_type"] in by_type_locus:
            by_type_locus[row["provirus_type"]][row["Locus"]].append(row)

    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 8.0), constrained_layout=True)
    for ax, type_name, letter, color in (
        (axes[0], "type1", "A", TYPE_I),
        (axes[1], "type2", "B", TYPE_II),
    ):
        summary = []
        for locus, locus_rows in by_type_locus[type_name].items():
            fractions = []
            for feature in FEATURES:
                calls = [gene_call(row, feature) for row in locus_rows]
                n = sum(callable_ for callable_, _ in calls)
                fractions.append(sum(hit for _, hit in calls) / n if n else math.nan)
            summary.append((locus, len(locus_rows), fractions))
        summary.sort(
            key=lambda item: (sum(v for v in item[2] if not math.isnan(v)), item[1]),
            reverse=True,
        )
        selected = summary[:35]
        matrix = np.array([item[2] for item in selected])
        image = ax.imshow(matrix, vmin=0, vmax=1, cmap="YlGnBu", aspect="auto")
        ax.set_xticks(np.arange(5), FEATURE_LABELS, rotation=30, ha="right")
        ax.set_yticks(
            np.arange(len(selected)),
            [display_locus_name(item[0]) for item in selected],
        )
        panel_title(ax, letter, f"{type_name.replace('type', 'Type ')} coding states")
        ax.spines[:].set_color(color)
    cbar = fig.colorbar(image, ax=axes, pad=0.02, shrink=0.55)
    cbar.set_label("Compatible fraction among evaluable proviral copies")
    path_heat = OUT / "Supplementary_ORF_by_type_current_v3r1.png"
    save_figure(fig, path_heat)
    plt.close(fig)

    pattern_counts = Counter((row["provirus_type"], pattern_tuple(row)) for row in prov)
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 7.3), constrained_layout=True)
    for ax, type_name, letter, color in (
        (axes[0], "type1", "A", TYPE_I),
        (axes[1], "type2", "B", TYPE_II),
    ):
        top = sorted(
            ((pattern, count) for (kind, pattern), count in pattern_counts.items() if kind == type_name),
            key=lambda item: item[1],
            reverse=True,
        )[:18][::-1]
        y = np.arange(len(top))
        ax.barh(y, [count for _, count in top], color=color)
        labels = []
        for pattern, _ in top:
            feature_labels = [
                label for label, present in zip(
                    [
                        "Gag",
                        "Gag–Pro",
                        "Gag–Pro–Pol",
                        "Env",
                        "Np9" if type_name == "type1" else "Rec",
                    ],
                    pattern,
                )
                if present
            ]
            labels.append(" + ".join(feature_labels) if feature_labels else "none")
        ax.set_yticks(y, labels)
        ax.set_xlabel("Artifact-filtered proviral observations")
        panel_title(ax, letter, f"{type_name.replace('type', 'Type ')} same-molecule patterns")
        finish_axis(ax, grid="x")
    path_patterns = OUT / "Supplementary_ORF_cooccurrence_current_v3r1.png"
    save_figure(fig, path_patterns)
    plt.close(fig)
    return [path_heat, path_patterns]


def build_supplementary_phylogeny_detail() -> Path:
    canvas = Image.new("RGB", (3300, 4200), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 72)
        small = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 56)
    except OSError:
        font = ImageFont.load_default()
        small = font
    for letter, label, source, left in (
        ("A", "gag", PHY / "hml2_pan_orf_gag_phylogeny.png", 40),
        ("B", "pro", PHY / "hml2_pan_orf_pro_phylogeny.png", 1670),
    ):
        draw.text((left, 20), letter, font=font, fill=INK)
        draw.text((left + 110, 30), label, font=small, fill=INK)
        panel = fit_panel(source, (1570, 4070))
        canvas.paste(panel, (left, 120))
    path = OUT / "Supplementary_gag_pro_phylogenies.png"
    canvas.save(path, dpi=(450, 450), optimize=True)
    return path


def draw_type_schematic(ax) -> None:
    """Draw one HML-2 map with the operational cassette and Δ292 nested within it."""
    ax.set_xlim(0, 9500)
    ax.set_ylim(0, 10)
    ax.axis("off")
    finding_title(ax, "A", "Type-I cassette")
    genes = [
        ("5′ LTR", 0, 968, 3.8, 2.0, "#8A949E"),
        ("gag", 1111, 2001, 5.0, 1.65, GREEN),
        ("pro", 2913, 1005, 5.0, 1.65, GOLD),
        ("pol", 3878, 2871, 5.0, 1.65, PURPLE),
        ("env", 6450, 2100, 2.75, 1.65, BLUE),
        ("3′ LTR", 8504, 968, 3.8, 2.0, "#8A949E"),
    ]
    ax.plot([970, 8530], [4.8, 4.8], color=INK, lw=1.0, zorder=0)
    for gene, x, width, y, height, color in genes:
        ax.add_patch(patches.Rectangle((x, y), width, height, color=color, alpha=0.9))
        ax.text(
            x + width / 2, y + height / 2, gene, ha="center", va="center",
            color="white", fontsize=9.0, fontweight="bold",
        )

    cassette_start, deletion_start, deletion_end, cassette_end = 6000, 6501, 6793, 7293
    cassette_polygon = [
        (cassette_start, 2.42),
        (cassette_start + 90, 3.05),
        (cassette_start, 3.68),
        (cassette_start + 90, 4.31),
        (cassette_start, 4.94),
        (cassette_start + 90, 5.57),
        (cassette_start, 6.20),
        (cassette_start + 90, 6.97),
        (cassette_end - 90, 6.97),
        (cassette_end, 6.20),
        (cassette_end - 90, 5.57),
        (cassette_end, 4.94),
        (cassette_end - 90, 4.31),
        (cassette_end, 3.68),
        (cassette_end - 90, 3.05),
        (cassette_end, 2.42),
    ]
    ax.add_patch(
        patches.Polygon(
            cassette_polygon,
            closed=True,
            facecolor="#D9E7E5",
            edgecolor=TYPE_I,
            linewidth=1.4,
            alpha=0.55,
        )
    )
    ax.add_patch(
        patches.Rectangle(
            (deletion_start, 2.22),
            deletion_end - deletion_start,
            4.95,
            facecolor="white",
            edgecolor=RED,
            linewidth=1.7,
            hatch="////",
        )
    )
    ax.annotate(
        "",
        xy=(cassette_start, 7.4),
        xytext=(cassette_end, 7.4),
        arrowprops=dict(arrowstyle="|-|", color=TYPE_I, lw=1.5),
    )
    ax.text(
        (cassette_start + cassette_end) / 2,
        8.0,
        "Cassette",
        ha="center",
        va="bottom",
        color=TYPE_I,
        fontsize=9.0,
        fontweight="normal",
    )
    ax.annotate(
        "Δ292",
        xy=((deletion_start + deletion_end) / 2, 2.18),
        xytext=((deletion_start + deletion_end) / 2, 1.15),
        ha="center",
        va="top",
        arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.1),
        color=RED,
        fontsize=9.0,
        fontweight="normal",
    )


def build_figure_5() -> Path:
    divergence = read_tsv(MECHANISM / "pairwise_divergence_summary.tsv")
    lineage_rows = [
        row for row in read_tsv(TYPE1_HYPOTHESES)
        if row["hypothesis"] == "H1_lower_host_cost"
        and row["ltr_subfamily"] in {"LTR5Hs", "LTR5A", "LTR5B"}
    ]
    lineage_counts = {
        type_name: Counter(
            row["ltr_subfamily"]
            for row in lineage_rows
            if row["direct_type"] == type_name
        )
        for type_name in ("TypeI", "TypeII")
    }
    effect_rows = [
        row for row in read_tsv(SOURCE_EFFECT / "model_comparison.tsv")
        if row["observation"] == "orthology_aware_floor"
        and row["ascertainment_scenario"] == "narrow"
        and row["model"] == "source_plus_focal_effect"
    ]
    effect_rows.sort(key=lambda row: float(row["source_opportunity_cv"]))

    apply_style()
    fig = plt.figure(figsize=(7.1, 4.75), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[0.72, 0.80])
    draw_type_schematic(fig.add_subplot(gs[0, :]))

    ax_b = fig.add_subplot(gs[1, 0])
    type_order = ["TypeI", "TypeII"]
    display = ["Type I", "Type II"]
    bottom = np.zeros(2)
    palette = {"LTR5Hs": TYPE_I, "LTR5A": GOLD, "LTR5B": ORANGE}
    for subfamily in ("LTR5Hs", "LTR5A", "LTR5B"):
        values = [lineage_counts[type_name][subfamily] for type_name in type_order]
        ax_b.bar(display, values, bottom=bottom, color=palette[subfamily], label=subfamily)
        bottom += values
    ax_b.set_ylabel("Loci")
    ax_b.legend(loc="upper left")
    finding_title(ax_b, "B", "LTR subfamilies")
    finish_axis(ax_b, grid="y")

    ax_c = fig.add_subplot(gs[1, 1])
    order = ["B1", "B2", "B3", "B4", "B5", "B6", "Cass."]
    for contrast, label, color, marker in (
        ("within_TypeI", "Type I", TYPE_I, "o"),
        ("within_TypeII", "Type II", TYPE_II, "s"),
        ("between_TypeI_TypeII", "Between", GOLD, "D"),
    ):
        selected = {r["window_label"]: r for r in divergence if r["contrast"] == contrast}
        ax_c.plot(
            order,
            [float(selected[key]["mean_p_distance"]) for key in order],
            marker=marker,
            lw=1.5,
            color=color,
            label=label,
        )
    ax_c.axvspan(5.55, 6.45, color="#ECE8DF", zorder=-2)
    ax_c.set_ylabel("Sequence divergence")
    ax_c.set_xlabel("Genome window")
    ax_c.set_ylim(0.02, 0.16)
    ax_c.legend(
        ncol=1,
        loc="upper right",
        bbox_to_anchor=(1.02, 1.02),
        columnspacing=0.8,
        handlelength=1.0,
        fontsize=8.5,
    )
    finding_title(ax_c, "C", "Regional divergence")
    finish_axis(ax_c, grid="y")

    ax_d = fig.add_subplot(gs[1, 2])
    x = np.array([float(row["source_opportunity_cv"]) for row in effect_rows])
    y = np.array([
        float(row["likelihood_ratio_effect_vs_source_only"])
        for row in effect_rows
    ])
    ax_d.plot(x, y, color=TYPE_I, marker="o", lw=1.8, ms=4)
    ax_d.axhline(1, color=INK, lw=0.9, ls="--")
    ax_d.set_yscale("log")
    ax_d.set_xscale("log")
    ax_d.set_xticks([x[0], 1.0, x[-1]], ["0.224", "1", "3.16"])
    ax_d.minorticks_off()
    ax_d.set_xlabel("Source-weight\ncoefficient of variation")
    ax_d.set_ylabel("Bayes factor")
    ax_d.text(
        x[0] * 1.06,
        1.12,
        "Equal support",
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=INK,
    )
    finding_title(ax_d, "D", "Added Δ292 effect")
    finish_axis(ax_d, grid="y")
    path = OUT / "Figure_5_type1_persistent_recombining_cassette.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def extract(pattern: str, text: str) -> float:
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"missing pattern {pattern!r}")
    return float(match.group(1).rstrip(".;"))


def forest(ax, rows, xlabel: str) -> None:
    y = np.arange(len(rows))[::-1]
    for yi, (label, beta, se, color) in zip(y, rows):
        ax.errorbar(beta, yi, xerr=1.96 * se, fmt="o", color=color, capsize=2.5, lw=1.2)
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_yticks(y, [row[0] for row in rows])
    ax.set_ylim(-0.5, len(rows) - 0.5)
    ax.set_xlabel(xlabel)
    finish_axis(ax, grid="x")


def build_figure_6() -> Path:
    slc_inputs = {
        row["cohort"]: row for row in read_tsv(SLC_INPUTS)
        if row["role"] == "PRIMARY_disjoint_pair"
    }
    growth = read_tsv(GROWTH_GROUPS)
    twelve_rows = {
        row["outcome_id"]: row for row in read_tsv(TWELVE_MODELS)
        if row["model_id"] == "structural_internal_vs_retained_noninternal"
    }
    apply_style()
    fig = plt.figure(figsize=(7.8, 2.72), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.86, 1.04, 1.42])

    ax_a = fig.add_subplot(gs[0, 0])
    slc_rows = []
    for label, cohort in (
        ("MAGE", "MAGE_v1_direct39"),
        ("GEUVADIS", "GEUVADIS_28_MAGE_disjoint"),
    ):
        source = slc_inputs[cohort]
        slc_rows.append((
            label,
            float(source["beta_std"]),
            float(source["se_std"]),
        ))
    x = np.arange(len(slc_rows))
    effects = np.array([row[1] for row in slc_rows])
    errors = 1.96 * np.array([row[2] for row in slc_rows])
    ax_a.errorbar(
        x,
        effects,
        yerr=errors,
        fmt="o",
        color=TYPE_I,
        markersize=6,
        capsize=3,
        lw=1.4,
    )
    ax_a.axhline(0, color=INK, lw=0.8)
    ax_a.set_xticks(x, [row[0] for row in slc_rows])
    ax_a.set_ylabel("Expression difference (SD)")
    ax_a.set_ylim(min(0, np.min(effects - errors)) - 0.08, np.max(effects + errors) + 0.08)
    finish_axis(ax_a, grid="y")
    finding_title(ax_a, "A", "SLC44A5 expression", letter_x=-0.20)

    ax_b = fig.add_subplot(gs[0, 1])
    populations = sorted({row["population"] for row in growth})
    growth_lookup = {
        (row["population"], row["exposure"]): float(row["mean_growth_per_10000"])
        for row in growth
    }
    for index, population in enumerate(populations):
        values = [growth_lookup[(population, exposure)] for exposure in ("0", "1")]
        ax_b.plot([0, 1], values, marker="o", color=GREEN, alpha=0.85)
        ax_b.text(1.04, values[1], population, va="center", fontsize=7.2)
    ax_b.set_xticks([0, 1], ["Provirus", "Solo LTR"])
    ax_b.set_ylabel("Growth rate (×10,000)")
    finding_title(ax_b, "B", "LCL growth", letter_x=-0.20)
    finish_axis(ax_b, grid="y")

    ax_c = fig.add_subplot(gs[0, 2])
    contexts = [
        ("Ofatumumab + serum", "Ofat_serum", RED),
        ("Obinutuzumab + serum", "Obin_serum", ORANGE),
        ("Rituximab + serum", "Ritux_serum", GOLD),
        ("Ofatumumab + media", "Ofat_media", BLUE),
        ("Obinutuzumab + media", "Obin_media", PURPLE),
    ]
    context_rows = []
    for label, outcome, color in contexts:
        row = twelve_rows[outcome]
        context_rows.append((
            label,
            float(row["beta"]),
            float(row["hc3_se"]),
            color,
        ))
    forest(ax_c, context_rows, "Difference in live-cell fraction")
    finding_title(ax_c, "C", "Anti-CD20 survival", letter_x=-0.20)
    path = OUT / "Figure_6_functional_association_leads.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def read_fasta_sequence(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text().splitlines()
        if line and not line.startswith(">")
    ).upper()


def solo_ltr_diversity(rows) -> list[dict[str, object]]:
    """Compute direct within-locus solo-LTR diversity from retained sequences."""
    sequences: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        if (
            row["observation_state"] != "PRESENT"
            or row["Structure"] != "Solo-LTR"
            or not PUBLIC_ID.fullmatch(row["ID"])
        ):
            continue
        path = PROCESSED_LOCI / row["Locus"] / f"{row['ID_Full']}.fa"
        if path.is_file():
            sequences[row["Locus"]].append(read_fasta_sequence(path))

    summaries: list[dict[str, object]] = []
    for locus, locus_sequences in sequences.items():
        if len(locus_sequences) < 20:
            continue
        lengths = {len(sequence) for sequence in locus_sequences}
        if len(lengths) != 1:
            continue
        n = len(locus_sequences)
        length = lengths.pop()
        total_pair_bases = (n * (n - 1) / 2) * length
        pair_differences = 0.0
        for column in range(length):
            counts = Counter(sequence[column] for sequence in locus_sequences)
            pair_differences += (
                n * n - sum(value * value for value in counts.values())
            ) / 2
        dominant = Counter(locus_sequences).most_common(1)[0][1]
        summaries.append(
            {
                "locus": locus.removeprefix("HML-2_"),
                "n": n,
                "sequence_length": length,
                "unique_sequences": len(set(locus_sequences)),
                "dominant_fraction": dominant / n,
                "nucleotide_diversity": pair_differences / total_pair_bases,
            }
        )
    return sorted(summaries, key=lambda row: row["nucleotide_diversity"])


def build_figure_7(rows) -> Path:
    return draw_eightq_network(WORKSPACE / "Supplementary_Data/Figure_7_solo_LTR_haplotype_counts.tsv", OUT / "Figure_7_8q11_structure_variation.png")


def build_type_state_counts() -> Path:
    return draw_type_state_counts(TYPE1_LOCUS_AUDIT, OUT / "Supplementary_type_state_counts.png")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, roster = load_catalog()
    recombination_path = build_peritelomeric_signature_figure(rows)
    outputs = [
        build_figure_1(rows, roster),
        build_figure_2(rows, roster),
        build_figure_3(rows, recombination_path),
        build_figure_4(rows),
        build_figure_5(),
        build_figure_6(),
        build_figure_7(rows),
    ]
    outputs.extend(build_supplementary_orf_figures(rows))
    outputs.append(build_supplementary_phylogeny_detail())
    outputs.append(recombination_path)
    outputs.append(build_type_state_counts())
    print("\n".join(str(path) for path in outputs))


if __name__ == "__main__":
    main()
