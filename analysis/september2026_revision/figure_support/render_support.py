from __future__ import annotations

import csv

import itertools

import json

import math

import re

from collections import Counter, defaultdict

from pathlib import Path

import matplotlib

import matplotlib.pyplot as plt

import matplotlib.patches as patches

import numpy as np

from Bio import Phylo, SeqIO

from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageOps

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

PHY = OUT / "phylogeny"

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
    / "hervk-hml2-age/manuscript_figures/python/analysis/"
    "type1_cassette_mechanism/results"
)

TYPE1_HYPOTHESES = (
    PROJECT / "working/type1_selection_hypotheses_v1/results/locus_inputs.tsv"
)

SOURCE_EFFECT = PROJECT / "manuscript/delta292_source_vs_effect_v1"

TYPE_I = "#163B75"
CASSETTE = "#B44787"

TYPE_II = "#735324"

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def display_locus_name(value: str) -> str:
    """Normalize a locus name for figures without changing analysis keys."""
    label = value.removeprefix("HML-2_")
    for suffix in ("_new", "_hg38"):
        label = label.removesuffix(suffix)
    return {
        "acro_type1": "Telomeric Type I",
        "acro_type2": "Telomeric Type II",
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

def draw_type_schematic(ax) -> None:
    """Show the exact windows used in panel C, excluding Δ292 from Cass."""
    ax.set_xlim(0, 9500)
    ax.set_ylim(-0.2, 10.5)
    ax.axis("off")
    finding_title(ax, "A", "Sequence comparison windows")
    genes = [
        ("5′ LTR", 0, 968, 8.8, "#8A949E"),
        ("gag", 1111, 3112, 9.1, "#527F4C"),
        ("pro", 2913, 3918, 9.1, "#444C56"),
        ("pol", 3878, 6749, 9.1, "#735394"),
        ("env", 6450, 8550, 8.1, "#519AC4"),
        ("3′ LTR", 8504, 9472, 8.8, "#8A949E"),
    ]
    ax.plot([968, 8504], [8.9, 8.9], color=INK, lw=0.7, zorder=0)
    for label, start, end, y, color in genes:
        ax.add_patch(patches.Rectangle((start, y), end-start, .78, color=color))
        ax.text((start+end)/2, y+.39, label, ha="center", va="center",
                color="white", fontsize=8, fontweight="bold")
    # Every backbone interval is exactly 1,000 KCON positions.
    for index, start in enumerate((1000, 2000, 3000, 4000, 5000, 7293), 1):
        ax.add_patch(patches.Rectangle((start, 6.6), 1000, .65,
                                      facecolor="#E9EAEC", edgecolor="#5B6369", lw=.7))
        ax.text(start+500, 6.925, f"B{index}", ha="center", va="center", fontsize=8)
    for start, end in ((6000, 6501), (6793, 7293)):
        ax.add_patch(patches.Rectangle((start, 6.6), end-start, .65,
                                      facecolor=CASSETTE, edgecolor=CASSETTE, lw=.7))
    ax.add_patch(patches.Rectangle((6501, 6.6), 292, .65,
                                  facecolor="white", edgecolor=INK, hatch="////", lw=.8))
    ax.text(6646.5, 7.55, "Cass.", ha="center", va="center", fontsize=8, color=CASSETTE)
    ax.text(1000, 5.9, "1000", ha="center", fontsize=6.5)
    ax.text(6000, 5.9, "6000", ha="right", fontsize=6.5)
    ax.text(7293, 5.9, "7293", ha="left", fontsize=6.5)
    ax.text(8293, 5.9, "8293", ha="center", fontsize=6.5)
    # Enlarge the same interval to make the two included flanks explicit.
    zoom_left, zoom_right = 1500, 8000
    zoom = lambda p: zoom_left + (p-6000)/(7293-6000)*(zoom_right-zoom_left)
    for source, target in ((6000, zoom_left), (7293, zoom_right)):
        ax.plot([source, target], [5.55, 4.25], color="#A9AFB3", lw=.65, ls="--")
    for start, end, label in ((6000, 6501, "501 bp"), (6793, 7293, "500 bp")):
        ax.add_patch(patches.Rectangle((zoom(start), 2.8), zoom(end)-zoom(start), 1.1,
                                      facecolor=CASSETTE, edgecolor=CASSETTE, lw=.8))
        ax.text((zoom(start)+zoom(end))/2, 3.35, label, ha="center", va="center",
                fontsize=8.5, color="white", fontweight="bold")
    ax.add_patch(patches.Rectangle((zoom(6501), 2.8), zoom(6793)-zoom(6501), 1.1,
                                  facecolor="white", edgecolor=INK, hatch="////", lw=1))
    ax.text(zoom(6647), 4.45, "Δ292", ha="center", fontsize=8.5, color=INK)
    for coordinate in (6000, 6501, 6793, 7293):
        ax.text(zoom(coordinate), 2.15, str(coordinate), ha="center", fontsize=7)
    ax.text(4750, 1.05, "Cass. comparison = 501 bp + 500 bp (Δ292 excluded)",
            ha="center", fontsize=8, color=CASSETTE)
    ax.text(4750, .15, "KCON coordinates are zero-based, half-open",
            ha="center", fontsize=7, color=INK)


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
    fig = plt.figure(figsize=(7.1, 5.55), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.02, 0.80])
    draw_type_schematic(fig.add_subplot(gs[0, :]))

    ax_b = fig.add_subplot(gs[1, 0])
    type_order = ["TypeI", "TypeII"]
    display = ["Type I", "Type II"]
    bottom = np.zeros(2)
    palette = {"LTR5Hs": "#007F73", "LTR5A": "#CEAA36", "LTR5B": "#D7733F"}
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
        ("within_TypeI", "Type I", "#163B75", "o"),
        ("within_TypeII", "Type II", "#735324", "s"),
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
    ax_d.plot(x, y, color=INK, marker="o", lw=1.8, ms=4)
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

# The retained Figure 3/5 entrypoints consume only these portable inputs.
_REVISION = Path(__file__).resolve().parents[1]
MECHANISM = _REVISION/'inputs/type1'
TYPE1_HYPOTHESES = _REVISION/'inputs/figures/locus_inputs.tsv'
SOURCE_EFFECT = _REVISION/'inputs/figures'
SUBFAMILY_AUTHORITY = _REVISION/'inputs/type1/subfamily_authority.tsv'
PHY = _REVISION/'inputs/phylogeny'
