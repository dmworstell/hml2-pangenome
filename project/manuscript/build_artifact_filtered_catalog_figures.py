#!/usr/bin/env python3
"""Rebuild catalog, VCF-QC, CNV/array, and ORF figures after artifact removal."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from structural_catalog_summary import structural_summary, write_structural_tables

from figure_style import (
    BLUE,
    GOLD,
    GRAY,
    GREEN,
    GRID,
    INK,
    LIGHT,
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
CATALOG = (
    PROJECT
    / "results/resolved_manuscript_catalog_20260914/"
    "combined_hml2_orf_analysis.RESOLVED.tsv"
)
CATALOG_SUMMARY = CATALOG.parent / "verification.json"
SHORT_READ = (
    PROJECT
    / "working/short_read_direction_corrected_v4/results/"
    "bio5_direction_corrected_v4/short_read_numeric_authority_v1.json"
)
SEVENP22 = (
    PROJECT
    / "working/sevenp22_proxy_resolution_agent/haplotype_copy_number_truth.tsv"
)
ONEP31_ARRAYS = (
    PROJECT
    / "working/onep31b_array_recovery_agent/haplotype_array_reconciliation.tsv"
)
CNV_SUPPLEMENT = (
    PROJECT / "manuscript/supplement/Table_S4_CNV_assembly_artifact_qc.tsv"
)
ONEQ22_AUTHORITY = (
    PROJECT
    / "working/cnv_copy_state_reinterpretation_v1/results/"
    "hg00423_1q22_copy_state_authority.v1.json"
)
DIRECT_MATRIX = (
    PROJECT
    / "working/direct_ebv_fitness_screen_v1/results/person_level_direct_matrix.tsv"
)

OUTDIR = PROJECT / "manuscript/figures/artifact_filtered"
SUPPLEMENT = PROJECT / "manuscript/supplement"
PUBLIC_ID = re.compile(r"^(?:HG|NA)\d+$")
COMPATIBLE = {"Intact", "Frameshift_at_end", "Intact_FS_End"}
PROVIRUS_STRUCTURES = {"Provirus", "Provirus_from_Multi"}
COLORS = {
    "Noncarrier call": LIGHT,
    "Solo-LTR": ORANGE,
    "Fragment": SKY,
    "Provirus": GREEN,
    "Multi-copy": PURPLE,
    "Unknown": GRAY,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


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
        "8p23.1_duplicate_group_unresolved": "8p23.1 copy group",
    }.get(label, label)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


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
        raise ValueError(f"unexpected filtered catalog dimensions: {len(rows)}, {len(roster)}")
    return rows, roster


def build_structural_figure(
    rows: list[dict[str, str]], roster: list[tuple[str, str]]
) -> Path:
    summary, observations, sex_evidence = structural_summary(rows, roster)
    states = list(COLORS)
    selected = sorted(
        (row for row in summary if row["label_scope"] == "physical_locus"),
        key=lambda row: row["variability_score"], reverse=True,
    )[:30]
    selected.reverse()
    write_structural_tables(summary, observations, sex_evidence, SUPPLEMENT)

    apply_style()
    fig, ax = plt.subplots(figsize=(7.1, 7.0), constrained_layout=True)
    left = np.zeros(len(selected))
    y = np.arange(len(selected))
    for state in states:
        values = np.array([row[state] / row["eligible_haplotypes"] for row in selected])
        ax.barh(y, values, left=left, color=COLORS[state], label=state, height=0.78)
        left += values
    ax.set_yticks(y, [display_locus_name(row["locus"]) for row in selected])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of eligible chromosome copies")
    panel_title(ax, "A", "Structural variation is concentrated at a subset of loci")
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.01))
    finish_axis(ax, grid="x")
    path = OUTDIR / "Figure_1_artifact_filtered_structural_spectrum.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def build_short_read_figure() -> Path:
    authority = json.loads(SHORT_READ.read_text())
    rows = []
    endpoint_labels = {
        "E2": "false alternate",
        "E5": "structural-state/type error",
    }
    for endpoint in ("E2", "E5"):
        by_locus = authority["cluster_weighted_empirical_diagnostics"][endpoint][
            "by_locus"
        ]
        for locus, payload in by_locus.items():
            rows.append(
                {
                    "endpoint": endpoint,
                    "endpoint_label": endpoint_labels[endpoint],
                    "locus": locus,
                    "sample_locus_clusters": payload["sample_locus_clusters"],
                    "cluster_weighted_error_fraction": payload[
                        "cluster_weighted_error_fraction"
                    ],
                    "record_rows": payload["record_rows"],
                    "analysis_role": "DESCRIPTIVE_DIAGNOSTIC_ONLY",
                }
            )
    write_tsv(
        SUPPLEMENT / "Table_S6_short_read_vcf_locus_diagnostics.tsv",
        rows,
        [
            "endpoint",
            "endpoint_label",
            "locus",
            "sample_locus_clusters",
            "cluster_weighted_error_fraction",
            "record_rows",
            "analysis_role",
        ],
    )

    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 4.35), constrained_layout=True)
    for ax, endpoint in zip(axes, ("E2", "E5")):
        subset = [
            row
            for row in rows
            if row["endpoint"] == endpoint and row["sample_locus_clusters"] >= 10
        ]
        subset = sorted(
            subset, key=lambda row: row["cluster_weighted_error_fraction"], reverse=True
        )[:15]
        subset.reverse()
        y = np.arange(len(subset))
        values = [row["cluster_weighted_error_fraction"] for row in subset]
        ax.barh(y, values, color=BLUE if endpoint == "E2" else RED, height=0.72)
        ax.set_yticks(y, [display_locus_name(row["locus"]) for row in subset])
        ax.set_xlim(0, 1)
        ax.set_xlabel("Fraction of matched events that disagree")
        panel_title(
            ax,
            "A" if endpoint == "E2" else "B",
            "Alternate calls contradicted by assemblies"
            if endpoint == "E2"
            else "Structural state disagrees with assemblies",
        )
        finish_axis(ax, grid="x")
    path = OUTDIR / "Figure_2_short_read_vcf_diagnostics.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def build_cnv_figure(
    rows: list[dict[str, str]], roster: list[tuple[str, str]]
) -> Path:
    seven_rows = read_tsv(SEVENP22)
    seven_counts = Counter(int(row["array_copy_number"]) for row in seven_rows)
    if sum(seven_counts.values()) != 584:
        raise ValueError("7p22.1 authority is not a complete 584-haplotype panel")

    arrays = {
        (row["sample"], row["haplotype"]): int(row["copy_number"])
        for row in read_tsv(ONEP31_ARRAYS)
    }
    one_cells: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["Locus"] == "HML-2_1p31.1b":
            one_cells[(row["ID"], row["Haplotype"])].append(row)
    one_cn = {}
    for key in roster:
        if key in arrays:
            one_cn[key] = arrays[key]
        else:
            present = [
                row
                for row in one_cells.get(key, [])
                if row["observation_state"] == "PRESENT"
            ]
            one_cn[key] = 1 if any(row["Structure"] == "Fragment" for row in present) else 0
    one_counts = Counter(one_cn.values())

    cnv_rows = read_tsv(CNV_SUPPLEMENT)
    artifacts = [
        row
        for row in cnv_rows
        if row["manuscript_analysis_decision"] == "EXCLUDE_EXACT_ASSEMBLY_RECORD"
    ]
    artifact_loci = Counter(display_locus_name(row["locus"]) for row in artifacts)
    hg00423 = next(
        row
        for row in artifacts
        if row["candidate_key"] == "HG00423|HML-2_1q22|region2"
    )

    array_rows = []
    for locus, counts in (("7p22.1", seven_counts), ("1p31.1b", one_counts)):
        for cn in sorted(counts):
            array_rows.append(
                {"locus": locus, "copy_number": cn, "haplotypes": counts[cn]}
            )
    write_tsv(
        SUPPLEMENT / "Table_S7_array_copy_number_distributions.tsv",
        array_rows,
        ["locus", "copy_number", "haplotypes"],
    )

    apply_style()
    fig = plt.figure(figsize=(7.1, 5.65), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    for ax, locus, counts, color in (
        (ax_a, "7p22.1", seven_counts, BLUE),
        (ax_b, "1p31.1b", one_counts, GREEN),
    ):
        x = sorted(counts)
        ax.bar(x, [counts[value] for value in x], color=color)
        ax.set_xticks(x)
        ax.set_xlabel("HML-2 copies per haplotype")
        ax.set_ylabel("Haplotypes")
        panel_title(ax, "A" if locus == "7p22.1" else "B", f"{locus} arrays")
        for xv, count in zip(x, [counts[value] for value in x]):
            ax.text(xv, count, str(count), ha="center", va="bottom", fontsize=6.5)
        finish_axis(ax, grid="y")

    top = artifact_loci.most_common(10)
    top.reverse()
    ax_c.barh(
        np.arange(len(top)),
        [count for _, count in top],
        color=PURPLE,
    )
    ax_c.set_yticks(np.arange(len(top)), [locus for locus, _ in top])
    ax_c.set_xlabel("Assembly records excluded")
    panel_title(ax_c, "C", "Excluded records by locus")
    for yi, count in enumerate([count for _, count in top]):
        ax_c.text(count + 0.12, yi, str(count), va="center", fontsize=6.5)
    finish_axis(ax_c, grid="x")

    labels = ["retained paternal alt1", "excluded paternal alt2"]
    body = [63, float(hg00423["raw_body_median"])]
    flank = [58, float(hg00423["raw_local_flank_median"])]
    x = np.arange(2)
    width = 0.34
    ax_d.bar(x - width / 2, body, width, label="HML-2 body", color=GREEN)
    ax_d.bar(x + width / 2, flank, width, label="local flank", color=GOLD)
    ax_d.axhline(
        float(hg00423["sample_global_diploid_target_flank_baseline"]),
        color="#222222",
        ls="--",
        label="sample baseline",
    )
    ax_d.set_xticks(x, labels, rotation=12, ha="right")
    ax_d.set_ylabel("Median read depth")
    panel_title(ax_d, "D", "One erroneous 1q22 assembly copy")
    ax_d.legend(loc="upper right")
    ax_d.text(1, max(body[1], flank[1]) + 2.5, "excluded", ha="center", color=RED, fontweight="bold")
    finish_axis(ax_d, grid="y")

    path = OUTDIR / "Figure_3_cnv_arrays_and_artifacts.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def gene_call(row: dict[str, str], gene: str) -> tuple[bool, bool]:
    if gene == "gag_pro_route":
        values = [row["gag"], row["pro"]]
        callable_ = all(value not in {"", "NA"} for value in values)
        return callable_, callable_ and all(value in COMPATIBLE for value in values)
    if gene == "gag_pro_pol_route":
        values = [row["gag"], row["pro"], row["pol"]]
        callable_ = all(value not in {"", "NA"} for value in values)
        return callable_, callable_ and all(value in COMPATIBLE for value in values)
    if gene == "accessory":
        field = "np9" if row["provirus_type"] == "type1" else "rec"
    else:
        field = gene
    value = row[field]
    callable_ = value not in {"", "NA"}
    return callable_, callable_ and value in COMPATIBLE


def build_orf_figure(rows: list[dict[str, str]]) -> Path:
    genes = ["gag", "gag_pro_route", "gag_pro_pol_route", "env", "accessory"]
    labels = ["Gag", "Gag–Pro", "Gag–Pro–Pol", "Env frame", "Np9/Rec"]
    by_locus: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if (
            row["observation_state"] == "PRESENT"
            and row["Structure"] in PROVIRUS_STRUCTURES
        ):
            by_locus[row["Locus"]].append(row)
    summary = []
    for locus, locus_rows in by_locus.items():
        result: dict[str, object] = {
            "locus": locus,
            "provirus_type": Counter(
                row["provirus_type"] for row in locus_rows
            ).most_common(1)[0][0],
            "proviral_copy_rows": len(locus_rows),
        }
        for gene in genes:
            calls = [gene_call(row, gene) for row in locus_rows]
            callable_n = sum(callable_ for callable_, _ in calls)
            compatible_n = sum(compatible for _, compatible in calls)
            result[f"{gene}_callable"] = callable_n
            result[f"{gene}_compatible"] = compatible_n
            result[f"{gene}_fraction"] = (
                compatible_n / callable_n if callable_n else math.nan
            )
        result["max_compatible_fraction"] = max(
            value
            for value in [result[f"{gene}_fraction"] for gene in genes]
            if not math.isnan(value)
        )
        summary.append(result)
    summary.sort(
        key=lambda row: (row["max_compatible_fraction"], row["proviral_copy_rows"]),
        reverse=True,
    )
    selected = summary[:30]
    fields = ["locus", "provirus_type", "proviral_copy_rows"]
    for gene in genes:
        fields.extend(
            [f"{gene}_callable", f"{gene}_compatible", f"{gene}_fraction"]
        )
    fields.append("max_compatible_fraction")
    write_tsv(
        SUPPLEMENT / "Table_S8_artifact_filtered_orf_coding_potential.tsv",
        summary,
        fields,
    )

    matrix = np.array(
        [[row[f"{gene}_fraction"] for gene in genes] for row in selected],
        dtype=float,
    )
    apply_style()
    fig, ax = plt.subplots(figsize=(7.1, 6.65), constrained_layout=True)
    image = ax.imshow(matrix, vmin=0, vmax=1, cmap="YlGnBu", aspect="auto")
    ax.set_xticks(np.arange(len(labels)), labels, rotation=22, ha="right")
    ylabels = [
        f"{display_locus_name(row['locus'])} ({row['provirus_type']})"
        for row in selected
    ]
    ax.set_yticks(np.arange(len(selected)), ylabels)
    for y, row in enumerate(selected):
        for x, gene in enumerate(genes):
            value = row[f"{gene}_fraction"]
            if not math.isnan(value):
                ax.text(
                    x,
                    y,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    color="white" if value > 0.66 else INK,
                )
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("Fraction sequence-compatible among evaluable proviral copies")
    ax.set_title("ORF annotations at 30 selected HML-2 loci and copy groups", loc="left", pad=6)
    path = OUTDIR / "Figure_4_orf_coding_potential.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def build_coding_path_disruption_figure(rows: list[dict[str, str]]) -> Path:
    """Summarize disruption along biologically meaningful reading paths."""
    by_locus: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if (
            row["observation_state"] == "PRESENT"
            and row["Structure"] in PROVIRUS_STRUCTURES
        ):
            by_locus[row["Locus"]].append(row)

    summaries = []
    for locus, locus_rows in by_locus.items():
        unspliced = Counter()
        env = Counter()
        for row in locus_rows:
            values = [row["gag"], row["pro"], row["pol"]]
            if all(value not in {"", "NA"} for value in values):
                if values[0] not in COMPATIBLE:
                    unspliced["Gag disrupted"] += 1
                elif values[1] not in COMPATIBLE:
                    unspliced["Pro disrupted"] += 1
                elif values[2] not in COMPATIBLE:
                    unspliced["Pol disrupted"] += 1
                else:
                    unspliced["Gag–Pro–Pol compatible"] += 1
            if row["env"] not in {"", "NA"}:
                env[
                    "Env compatible"
                    if row["env"] in COMPATIBLE
                    else "Env disrupted"
                ] += 1
        if sum(unspliced.values()):
            summaries.append(
                {
                    "locus": locus,
                    "unspliced": unspliced,
                    "env": env,
                    "unspliced_n": sum(unspliced.values()),
                    "env_n": sum(env.values()),
                }
            )

    selected = sorted(
        summaries,
        key=lambda row: (row["unspliced_n"], row["env_n"]),
        reverse=True,
    )[:35][::-1]
    unspliced_order = [
        "Gag disrupted",
        "Pro disrupted",
        "Pol disrupted",
        "Gag–Pro–Pol compatible",
    ]
    env_order = ["Env disrupted", "Env compatible"]
    palette = {
        "Gag disrupted": GOLD,
        "Pro disrupted": SKY,
        "Pol disrupted": GREEN,
        "Gag–Pro–Pol compatible": LIGHT,
        "Env disrupted": ORANGE,
        "Env compatible": GRAY,
    }

    apply_style()
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(7.1, 8.2), sharey=True, constrained_layout=True
    )
    y = np.arange(len(selected))
    for ax, field, categories, columns in (
        (ax_a, "unspliced", unspliced_order, 2),
        (ax_b, "env", env_order, 1),
    ):
        left = np.zeros(len(selected))
        for category in categories:
            values = np.array(
                [
                    row[field][category] / sum(row[field].values())
                    if sum(row[field].values())
                    else 0
                    for row in selected
                ]
            )
            ax.barh(
                y,
                values,
                left=left,
                color=palette[category],
                label=category,
                height=0.78,
            )
            left += values
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, len(selected) + 3.5)
        ax.set_xlabel("Fraction of evaluable proviral copies")
        ax.legend(loc="upper right", ncol=columns)
        finish_axis(ax, grid="x")
    ax_a.set_yticks(
        y, [display_locus_name(row["locus"]) for row in selected]
    )
    ax_b.tick_params(axis="y", left=False, labelleft=False)
    panel_title(ax_a, "A", "Gag–Pro–Pol path")
    panel_title(ax_b, "B", "Env frame")
    path = OUTDIR / "Figure_S8_coding_path_disruption.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def build_apparent_copy_review_figure() -> Path:
    """Show the current review outcome for apparent additional assembly copies."""
    rows = read_tsv(CNV_SUPPLEMENT)
    state_order = [
        "assembly_artifact",
        "authenticated_segdup",
        "later_duplication",
    ]
    labels = {
        "assembly_artifact": "Assembly artifact",
        "authenticated_segdup": "Segmental duplication",
        "later_duplication": "Later duplication",
    }
    colors = {
        "assembly_artifact": RED,
        "authenticated_segdup": PURPLE,
        "later_duplication": GREEN,
    }
    totals = Counter(row["dominant_state"] for row in rows)
    by_locus: dict[str, Counter] = defaultdict(Counter)
    for row in rows:
        by_locus[row["locus"]][row["dominant_state"]] += 1
    loci = sorted(
        by_locus,
        key=lambda locus: (sum(by_locus[locus].values()), locus),
    )

    apply_style()
    fig = plt.figure(figsize=(7.1, 6.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.1, 1.3])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[:, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    x = np.arange(len(state_order))
    values = [totals[state] for state in state_order]
    ax_a.bar(x, values, color=[colors[state] for state in state_order])
    ax_a.set_xticks(x, [labels[state] for state in state_order], rotation=28, ha="right")
    ax_a.set_ylabel("Reviewed assembly records")
    for index, value in enumerate(values):
        ax_a.text(index, value + 0.7, str(value), ha="center", va="bottom")
    panel_title(ax_a, "A", "Review outcome")
    finish_axis(ax_a, grid="y")

    y = np.arange(len(loci))
    left = np.zeros(len(loci))
    for state in state_order:
        values = np.array([by_locus[locus][state] for locus in loci])
        ax_b.barh(
            y,
            values,
            left=left,
            color=colors[state],
            label=labels[state],
            height=0.76,
        )
        left += values
    ax_b.set_yticks(y, [display_locus_name(locus) for locus in loci])
    ax_b.set_xlabel("Reviewed assembly records")
    ax_b.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1)
    panel_title(ax_b, "B", "Outcome by locus")
    finish_axis(ax_b, grid="x")

    depth_source = PROJECT / "working/cnv_copy_state_reinterpretation_v1/source_snapshot/cnv_depth_summary.tsv"
    if sha256(depth_source) != "566cec9fc714c5b3f37f6cd892d800c87bcd46cd9b5b66fbe5a17497d0d12923":
        raise ValueError("HG00423 source depth table changed")
    depth = {row["region"]: row for row in read_tsv(depth_source)
             if row["sample"] == "HG00423" and row["locus"] == "HML-2_1q22"}
    if set(depth) != {"1", "2"}:
        raise ValueError("HG00423 1q22 must have exactly two measured candidate regions")
    body = [float(depth[region]["body_median_all"]) for region in ("1", "2")]
    flank = [float(depth[region]["flank_median_all"]) for region in ("1", "2")]
    baseline = {float(row["sample_ref_cov"]) for row in depth.values()}
    if len(baseline) != 1:
        raise ValueError("Inconsistent sample-wide depth baseline")
    x = np.arange(2)
    ax_c.bar(x - .17, body, .34, label="HML-2 body", color=GREEN)
    ax_c.bar(x + .17, flank, .34, label="Host flank", color=GOLD)
    ax_c.axhline(baseline.pop(), color=INK, ls="--", label="Sample baseline")
    ax_c.set_xticks(x, ["Retained copy", "Excluded copy"], rotation=12, ha="right")
    ax_c.set_ylabel("Median read depth")
    ax_c.set_ylim(0, 90)
    panel_title(ax_c, "C", "HG00423 at 1q22")
    ax_c.legend(loc="upper right", fontsize=7)
    finish_axis(ax_c, grid="y")

    path = OUTDIR / "Figure_S9_apparent_extra_copy_review.png"
    save_figure(fig, path)
    plt.close(fig)
    return path


def build_oneq22_audit(rows: list[dict[str, str]]) -> None:
    locus_rows = [
        row
        for row in rows
        if row["Locus"] == "HML-2_1q22"
        and row["observation_state"] == "PRESENT"
    ]
    by_person: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in locus_rows:
        by_person[row["ID"]].append(row)
    audit = []
    for sample, sample_rows in sorted(by_person.items()):
        provirus = [
            row for row in sample_rows if row["Structure"] in PROVIRUS_STRUCTURES
        ]
        compatible_gag = sum(row["gag"] in COMPATIBLE for row in provirus)
        audit.append(
            {
                "sample": sample,
                "observed_1q22_haplotypes": len(sample_rows),
                "provirus_haplotypes": len(provirus),
                "solo_ltr_haplotypes": sum(
                    row["Structure"] == "Solo-LTR" for row in sample_rows
                ),
                "gag_compatible_dosage": compatible_gag,
                "gag_marker_state": (
                    "gag_compatible" if compatible_gag else "gag_not_compatible"
                ),
                "artifact_policy": (
                    "HG00423 paternal alt2 exact record removed; retained "
                    "maternal and paternal alt1 remain"
                    if sample == "HG00423"
                    else "no 1q22 assembly-artifact exclusion"
                ),
            }
        )
    write_tsv(
        SUPPLEMENT / "Table_S9_1q22_Gag_artifact_corrected_truth.tsv",
        audit,
        list(audit[0]),
    )
    matrix = {row["sample"]: row for row in read_tsv(DIRECT_MATRIX)}
    hg = matrix.get("HG00423", {})
    outcome_fields = [
        key
        for key in ("ebv_load", "log2_ebv_load")
        if hg.get(key, "") not in {"", "NA"}
    ]
    summary = {
        "schema": "hml2.1q22-artifact-corrected-audit.v1",
        "people": len(audit),
        "proviral_haplotypes": sum(row["provirus_haplotypes"] for row in audit),
        "solo_ltr_haplotypes": sum(row["solo_ltr_haplotypes"] for row in audit),
        "gag_compatible_haplotypes": sum(
            row["gag_compatible_dosage"] for row in audit
        ),
        "gag_not_compatible_haplotypes": sum(
            row["provirus_haplotypes"] - row["gag_compatible_dosage"]
            for row in audit
        ),
        "person_dosage_counts": dict(
            sorted(Counter(row["gag_compatible_dosage"] for row in audit).items())
        ),
        "HG00423_current_gag_compatible_dosage": next(
            row["gag_compatible_dosage"]
            for row in audit
            if row["sample"] == "HG00423"
        ),
        "HG00423_nonmissing_direct_matrix_fields": outcome_fields,
        "coefficient_impact_statement": (
            "HG00423 has no nonmissing direct-screen phenotype outcome; restoring "
            "its two valid 1q22 haplotypes changes truth-vector denominators but "
            "does not change fitted rows in those outcome models."
        ),
    }
    (
        OUTDIR / "oneq22_Gag_artifact_corrected_audit.json"
    ).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows, roster = load_catalog()
    figures = [
        build_structural_figure(rows, roster),
        build_short_read_figure(),
        build_cnv_figure(rows, roster),
        build_orf_figure(rows),
        build_coding_path_disruption_figure(rows),
        build_apparent_copy_review_figure(),
    ]
    build_oneq22_audit(rows)
    provenance = {
        "schema": "hml2.biologically-filtered-manuscript-figures.v2",
        "catalog": str(CATALOG.relative_to(WORKSPACE)),
        "catalog_sha256": sha256(CATALOG),
        "catalog_summary_sha256": sha256(CATALOG_SUMMARY),
        "excluded_exact_assembly_records": 35,
        "excluded_8q24.3b_alias_rows_in_public_panel": 584,
        "excluded_duplicate_catalog_labels_in_public_panel": 78,
        "exclusion_scope": (
            "Rows remain in the annotated catalog. Analyses use only analysis_include=1; "
            "valid copies, samples, and loci remain."
        ),
        "short_read_policy": (
            "phased-VCF locus diagnostics only; no consensus-assembly recovery "
            "claim, no pooled rate, and no biological ORF inference"
        ),
        "figures": {
            str(path.relative_to(WORKSPACE)): sha256(path) for path in figures
        },
    }
    (OUTDIR / "artifact_filtered_figures.provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    for figure in figures:
        print(figure)


if __name__ == "__main__":
    main()
