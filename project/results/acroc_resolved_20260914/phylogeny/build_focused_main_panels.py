#!/usr/bin/env python3
"""Prune existing full trees into the focused main-figure views; never refit."""

from __future__ import annotations

import copy
import csv
import itertools
import json
import os
import re
import sys
from hashlib import sha256
from pathlib import Path

OUT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
os.environ["MPLCONFIGDIR"] = str(OUT / "mplcache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

FOCAL = ("19p12c", "19p12d", "12q14.1", "10q24.2")
TREE_FILES = {"LTR": "hml2_pan_ltr_expanded_tree.nwk", "Pol": "hml2_pan_orf_pol_tree.nwk"}


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def observation_count(tip):
    return int(re.search(r"__n(\d+)$", tip.name).group(1))


def locus(tip):
    return tip.name.split("__", 1)[0]


def label(tip):
    loc, cluster, _count = tip.name.split("__")
    loc = loc.removesuffix("_new").removesuffix("_hg38")
    return f"{loc} c{cluster.removeprefix('hap')} (n={observation_count(tip)})"


def draw_tree(ax, tree):
    terminals = tree.get_terminals()
    x = tree.depths()
    y = {tip: index for index, tip in enumerate(reversed(terminals))}

    def position(node):
        for child in node.clades:
            if child not in y:
                position(child)
        y[node] = (y[node.clades[0]] + y[node.clades[-1]]) / 2

    position(tree.root)
    xmax = max(x.values())
    label_x = xmax * 1.06
    for node in tree.find_clades():
        if not node.clades:
            continue
        child_y = [y[child] for child in node.clades]
        ax.plot([x[node], x[node]], [min(child_y), max(child_y)], color="#727D85", lw=0.75)
        for child in node.clades:
            ax.plot([x[node], x[child]], [y[child], y[child]], color="#727D85", lw=0.75)
    for tip in terminals:
        focal = locus(tip) in FOCAL
        ax.plot([x[tip], label_x * 0.98], [y[tip], y[tip]], color="#C3CBD1", lw=0.45)
        ax.text(label_x, y[tip], label(tip), va="center", color="#006CAB" if focal else "#333A40",
                fontsize=9.0, weight="bold" if focal else "normal")
    ax.set_xlim(-0.015 * xmax, xmax * 2.58)
    ax.set_ylim(-0.65, len(terminals) - 0.35)
    ax.set_yticks([])
    ax.set_xlabel("Substitutions per site", fontsize=8, labelpad=3)
    ax.tick_params(axis="x", labelsize=7.5, length=3)
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.7)
    return len(terminals)


def verify_label_geometry(fig, axes):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [(text.get_text(), text.get_window_extent(renderer)) for ax in axes for text in ax.texts]
    for index, (first, box) in enumerate(boxes):
        assert box.x0 >= 0 and box.x1 <= fig.bbox.width, (first, "horizontal clip")
        assert box.y0 >= 0 and box.y1 <= fig.bbox.height, (first, "vertical clip")
        for second, other in boxes[index + 1:]:
            assert not box.overlaps(other), (first, second)


def render(tree, region, basename):
    plt.rcParams.update({"font.family": "Arial", "font.size": 9, "svg.fonttype": "none", "savefig.bbox": None})
    fig = plt.figure(figsize=(3.5, 3.25))
    ax = fig.add_axes([0.035, 0.19, 0.94, 0.69])
    tip_count = draw_tree(ax, tree)
    fig.text(0.02, 0.965, "A" if region == "LTR" else "B", va="top", fontsize=12, weight="bold")
    fig.text(0.52, 0.965, region, va="top", ha="center", fontsize=11)
    fig.text(0.50, 0.065, "c1/c2: modal sequence clusters; n: cluster observations", ha="center", fontsize=6.5)
    if region == "Pol":
        fig.text(0.50, 0.018, "19p12d: no callable Pol sequence in the retained panel", ha="center", fontsize=6.5)
    verify_label_geometry(fig, [ax])
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"{basename}.{suffix}", dpi=450, facecolor="white", bbox_inches=None)
    plt.close(fig)
    return {"panel": basename, "width_inches": 3.5, "height_inches": 3.25,
            "tip_labels": tip_count, "font_points": 9, "overlapping_labels": 0, "clipped_labels": 0}


def render_pair(trees):
    fig = plt.figure(figsize=(7.1, 2.5))
    axes = []
    for offset, region in ((0.0, "LTR"), (3.6, "Pol")):
        ax = fig.add_axes([(offset + 0.035 * 3.5) / 7.1, 0.15, 0.94 * 3.5 / 7.1, 0.785])
        axes.append(ax)
        draw_tree(ax, trees[region])
        fig.text((offset + 1.75) / 7.1, 0.995, region, va="top", ha="center", fontsize=11)
    fig.text(0.012, 0.995, "A", va="top", fontsize=12, weight="bold")
    verify_label_geometry(fig, axes)
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"Main_LTR_Pol_focus.{suffix}", dpi=450, facecolor="white", bbox_inches=None)
    plt.close(fig)
    return {"panel": "Main_LTR_Pol_focus", "width_inches": 7.1, "height_inches": 2.5,
            "tip_labels": sum(tree.count_terminals() for tree in trees.values()),
            "font_points": 9, "overlapping_labels": 0, "clipped_labels": 0}


def main():
    trees = {region: Phylo.read(OUT / name, "newick") for region, name in TREE_FILES.items()}
    original_bytes = {region: (OUT / name).read_bytes() for region, name in TREE_FILES.items()}
    shared_loci = set.intersection(*(
        {locus(tip) for tip in tree.get_terminals()} for tree in trees.values()
    ))
    nearest_rows, all_distances = [], []
    selected_loci = set(FOCAL)
    for region, tree in trees.items():
        terminals = tree.get_terminals()
        for focal in FOCAL:
            candidates = [tip for tip in terminals if locus(tip) == focal]
            if not candidates:
                nearest_rows.append({"region": region, "focal_locus": focal, "focal_tip": "",
                    "nearest_locus": "", "nearest_tip": "", "patristic_distance": "",
                    "status": "not_in_callable_full_tree", "comparison_locus_count": len(shared_loci)})
                continue
            focal_tip = sorted(candidates, key=lambda tip: (-observation_count(tip), tip.name))[0]
            distances = sorted((tree.distance(focal_tip, tip), tip.name, locus(tip)) for tip in terminals if locus(tip) != focal)
            shared_distances = [item for item in distances if item[2] in shared_loci]
            nearest = shared_distances[0][0] if focal in shared_loci else None
            if nearest is None:
                nearest_rows.append({"region": region, "focal_locus": focal, "focal_tip": focal_tip.name,
                    "nearest_locus": "", "nearest_tip": "", "patristic_distance": "",
                    "status": "focal_not_in_shared_locus_set", "comparison_locus_count": len(shared_loci)})
            for distance, tip, other_locus in distances:
                all_distances.append({"region": region, "focal_locus": focal, "focal_tip": focal_tip.name,
                    "other_locus": other_locus, "other_tip": tip, "patristic_distance": f"{distance:.15g}",
                    "included_in_shared_locus_comparison": focal in shared_loci and other_locus in shared_loci})
                if nearest is not None and other_locus in shared_loci and abs(distance - nearest) <= 1e-12:
                    nearest_rows.append({"region": region, "focal_locus": focal, "focal_tip": focal_tip.name,
                        "nearest_locus": other_locus, "nearest_tip": tip, "patristic_distance": f"{distance:.15g}",
                        "status": "nearest_among_shared_loci_outside_focal_locus", "comparison_locus_count": len(shared_loci)})
                    selected_loci.add(other_locus)
    tip_rows, pair_rows, panels, pruned_trees = [], [], [], {}
    for region, original in trees.items():
        keep = {tip.name for tip in original.get_terminals() if locus(tip) in selected_loci}
        focused = copy.deepcopy(original)
        for tip in list(focused.get_terminals()):
            if tip.name not in keep:
                focused.prune(tip)
        assert {tip.name for tip in focused.get_terminals()} == keep
        for row in nearest_rows:
            if row["region"] == region and row["nearest_tip"]:
                assert row["focal_tip"] in keep and row["nearest_tip"] in keep
        for first, second in itertools.combinations(sorted(keep), 2):
            full_distance = original.distance(first, second)
            focal_distance = focused.distance(first, second)
            difference = abs(full_distance - focal_distance)
            assert difference < 1e-10, (region, first, second, difference)
            pair_rows.append({"region": region, "tip_1": first, "tip_2": second,
                "full_tree_distance": f"{full_distance:.15g}", "focused_tree_distance": f"{focal_distance:.15g}",
                "absolute_difference": f"{difference:.15g}"})
        for tip in focused.get_terminals():
            tip_rows.append({"region": region, "tip": tip.name, "locus": locus(tip),
                "cluster_observations": observation_count(tip), "focal_locus": int(locus(tip) in FOCAL)})
        filename = f"focused_{region.lower()}_tree.nwk"
        Phylo.write(focused, OUT / filename, "newick", format_branch_length="%1.15g")
        reread = Phylo.read(OUT / filename, "newick")
        pruned_trees[region] = reread
        for first, second in itertools.combinations(sorted(keep), 2):
            assert abs(original.distance(first, second) - reread.distance(first, second)) < 1e-10
        panels.append(render(reread, region, "Main_LTR" if region == "LTR" else "Main_Pol"))
        assert original_bytes[region] == (OUT / TREE_FILES[region]).read_bytes()
    panels.append(render_pair(pruned_trees))
    write_tsv(OUT / "focused_main_selected_tips.tsv", tip_rows, list(tip_rows[0]))
    write_tsv(OUT / "focused_main_nearest_neighbors.tsv", nearest_rows, list(nearest_rows[0]))
    write_tsv(OUT / "focused_main_full_tree_distances.tsv", all_distances, list(all_distances[0]))
    write_tsv(OUT / "focused_main_branch_length_verification.tsv", pair_rows, list(pair_rows[0]))
    summary = {"status": "PASS_PRUNED_WITHOUT_REFITTING", "focal_loci": list(FOCAL),
        "comparison_loci": sorted(shared_loci), "comparison_locus_count": len(shared_loci),
        "context_loci": sorted(selected_loci - set(FOCAL)), "selected_loci": sorted(selected_loci),
        "full_tree_sha256": {region: sha256(value).hexdigest() for region, value in original_bytes.items()},
        "selected_tip_counts": {region: sum(row["region"] == region for row in tip_rows) for region in trees},
        "verified_pairwise_distances": len(pair_rows),
        "maximum_distance_difference": max(float(row["absolute_difference"]) for row in pair_rows),
        "selection_rule": "Use the most frequent retained cluster at each preselected focal locus (lexical tip-ID tie break). Restrict focal and candidate loci to those represented in both LTR and Pol trees. Search every retained other-locus cluster within this common set for its patristic nearest neighbor, including all exact ties. Distances use the existing full trees without refitting. Include both retained clusters of every preselected focal or shared-set nearest-neighbor locus in both panels when callable. c1/c2 identify frequency-ranked sequence clusters within each region, not phased cross-region haplotypes.",
        "nearest_neighbors": nearest_rows, "panels": panels,
        "interpretation": "19p12c, 12q14.1, and 10q24.2 have different nearest-locus neighborhoods between LTR and Pol after restricting candidates to the loci represented in both retained trees. 19p12d lacks callable Pol and is excluded from this paired comparison. Nearest-neighbor differences are descriptive and do not establish a recombination mechanism or statistical support."}
    (OUT / "focused_main_verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    full_summary_path = OUT / "verification.json"
    if full_summary_path.exists():
        full_summary = json.loads(full_summary_path.read_text())
        full_summary["main_panels"] = summary
        full_summary["render_geometry"] = [row for row in full_summary.get("render_geometry", []) if not row["panel"].startswith("Main_")] + panels
        full_summary_path.write_text(json.dumps(full_summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
