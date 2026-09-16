"""Recompute regional nearest neighbors using the same locus candidate set.

This checks the existing published NJ trees. It does not infer new trees or
equate frequency-ranked clusters across different gene regions.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio import Phylo


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phylogeny", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--figure-manifest", type=Path)
    args = parser.parse_args()
    filenames = {
        "LTR": "hml2_pan_ltr_expanded_tree.nwk",
        "Pol": "hml2_pan_orf_pol_tree.nwk",
    }
    trees = {
        region: Phylo.read(args.phylogeny / filename, "newick")
        for region, filename in filenames.items()
    }
    tips = {region: tree.get_terminals() for region, tree in trees.items()}
    loci = {
        region: {tip.name.split("__")[0] for tip in leaves}
        for region, leaves in tips.items()
    }
    common = loci["LTR"] & loci["Pol"]
    rows, summaries = [], []
    for focal in ("19p12c", "10q24.2"):
        if focal not in common:
            raise ValueError(f"Focal locus absent from one region: {focal}")
        for region, tree in trees.items():
            queries = [
                tip for tip in tips[region]
                if tip.name.startswith(focal + "__hap1__")
            ]
            if len(queries) != 1:
                raise ValueError(f"Expected one modal cluster: {region} {focal}")
            query = queries[0]
            candidates = sorted(
                (tree.distance(query, tip), tip.name, tip.name.split("__")[0])
                for tip in tips[region]
                if tip.name.split("__")[0] != focal
            )
            for distance, name, locus in candidates:
                rows.append({
                    "region": region,
                    "focal_locus": focal,
                    "focal_tip": query.name,
                    "candidate_locus": locus,
                    "candidate_tip": name,
                    "patristic_distance": distance,
                    "candidate_in_both_regions": locus in common,
                })
            for scope in ("all_region_loci", "shared_loci_only"):
                eligible = [
                    item for item in candidates
                    if scope == "all_region_loci" or item[2] in common
                ]
                minimum = eligible[0][0]
                nearest = [item for item in eligible if abs(item[0] - minimum) < 1e-10]
                summaries.append({
                    "region": region,
                    "focal_locus": focal,
                    "focal_tip": query.name,
                    "comparison": scope,
                    "candidate_loci": len(loci[region] - {focal}) if scope == "all_region_loci" else len(common - {focal}),
                    "nearest_tips": [item[1] for item in nearest],
                    "nearest_loci": sorted({item[2] for item in nearest}),
                    "minimum_patristic_distance": minimum,
                })
    if args.figure_manifest:
        manifest = json.loads(args.figure_manifest.read_text())
        for summary in summaries:
            if summary["comparison"] != "shared_loci_only":
                continue
            selected = set(manifest["trees"][summary["region"]]["selected_tips"])
            required = {summary["focal_tip"], *summary["nearest_tips"]}
            summary["required_tips_already_in_main_figure"] = required <= selected
            if not required <= selected:
                raise ValueError(f"Main figure does not display {required - selected}")
    # Verify the adjacent three-gene example against one common locus set too.
    all_trees = dict(trees)
    for region in ("Gag", "Pro", "Env"):
        filename = f"hml2_pan_orf_{region.lower()}_tree.nwk"
        filenames[region] = filename
        all_trees[region] = Phylo.read(args.phylogeny / filename, "newick")
    all_loci = {
        region: {tip.name.split("__")[0] for tip in tree.get_terminals()}
        for region, tree in all_trees.items()
    }
    common_five = set.intersection(*all_loci.values())
    five_region_comparison = []
    for region, tree in all_trees.items():
        query = next(tip for tip in tree.get_terminals() if tip.name.startswith("19p12c__hap1__"))
        candidates = sorted(
            (tree.distance(query, tip), tip.name)
            for tip in tree.get_terminals()
            if tip.name.split("__")[0] in common_five - {"19p12c"}
        )
        minimum = candidates[0][0]
        five_region_comparison.append({
            "region": region,
            "focal_tip": query.name,
            "candidate_loci": len(common_five - {"19p12c"}),
            "nearest_tips": [name for distance, name in candidates if abs(distance - minimum) < 1e-10],
            "minimum_patristic_distance": minimum,
        })
    result = {
        "method": "Patristic distance from each region's modal focal sequence cluster, restricting candidate loci to the intersection represented in both existing regional trees. Existing topology and branch lengths are unchanged. Sequence clusters are region-specific, not matched donor haplotypes.",
        "locus_counts": {region: len(value) for region, value in loci.items()},
        "shared_locus_count": len(common),
        "shared_loci": sorted(common),
        "LTR_only": sorted(loci["LTR"] - common),
        "Pol_only": sorted(loci["Pol"] - common),
        "five_region_shared_locus_count": len(common_five),
        "five_region_shared_loci": sorted(common_five),
        "five_region_19p12c_comparison": five_region_comparison,
        "tree_input_sha256": {
            filename: hashlib.sha256((args.phylogeny / filename).read_bytes()).hexdigest()
            for filename in filenames.values()
        },
        "comparisons": summaries,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "regional_tree_candidate_distances.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.output / "shared_locus_comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"shared_loci": len(common), "comparisons": summaries}, indent=2))


if __name__ == "__main__":
    main()
