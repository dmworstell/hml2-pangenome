#!/usr/bin/env python3
"""Apply the repaired assignment reader and rebuild the retained DNA network."""

import argparse
import csv
import importlib.util
import json
import os
import sys
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
BIOLOGICAL = PROJECT / "working/biological_orf_annotation_20260802"
EXACT_AUDIT = PROJECT / "results/biological_orf_annotation_20260802/exact_duplicate_locus_assignment_audit.tsv"
ATTACHMENTS = PROJECT / "results/duplicated_locus_attachment_audit_20260802/duplicated_locus_attachment_scan.exact_cigar.tsv"
POOLS = {"acro_type1", "acro_type2", "8p23.1_duplicate_group_unresolved", "Yq11.23_duplicate_group_unresolved"}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def graph_paf_alignment(audit, hit, locus):
    """Feed native-window PAF geometry to the existing exact CIGAR projector."""
    build = {"CHM13": "t2t", "GRCh38": "hg38"}[hit["build"]]
    length, start, end = hit["query_length"], hit["query_start"], hit["query_end"]
    assert 0 <= start < end <= length and hit["strand"] in {"+", "-"}
    left, right = (start, length-end) if hit["strand"] == "+" else (length-end, start)
    cigar = (f"{left}S" if left else "") + hit["cg"] + (f"{right}S" if right else "")
    operations = audit.parse_cigar(hit["cg"])
    assert sum(n for n, op in operations if op in "MI=X") == end-start
    assert sum(n for n, op in operations if op in "MDN=X") == hit["target_end"]-hit["target_start"]
    return audit.Alignment(build, locus, hit["query"], hit["target"], hit["target_start"],
        hit["target_end"], start, end, length, hit["strand"], hit["mapq"], hit.get("tp") == "P", cigar)


def graph_source_candidates(audit, row, windows):
    """Require one byte-bound native element position, retaining every competing hit."""
    positions = set()
    for native_window in row["native_windows"]:
        window = windows[native_window["query"]]
        source, requested_start, requested_end = audit.source_interval(row["Source_Identifier"])
        assert window["native_path"] == source
        assert window["window_end"]-window["window_start"] == window["window_length"]
        for match in native_window["exact_sequence_matches"]:
            start, end = match["element_start0"], match["element_end0"]
            assert requested_start <= start < end <= requested_end
            assert end-start == row["retained_sequence_length"]
            assert start-window["window_start"] == match["element_start_in_window"]
            assert end-window["window_start"] == match["element_end_in_window"]
            assert 0 <= match["element_start_in_window"] < match["element_end_in_window"] <= window["window_length"]
            assert match["retained_sequence_orientation_vs_native"] in {"+", "-"}
            positions.add((start, end, match["retained_sequence_orientation_vs_native"]))
    assert len(positions) == row["unique_exact_sequence_positions"]
    if len(positions) != 1:
        return [], [], positions
    assert row["capture_status"] == "EXACT_ELEMENT_MATCH"
    best = {}
    matched_hits = []
    for native_window in row["native_windows"]:
        window = windows[native_window["query"]]
        for match in native_window["exact_sequence_matches"]:
            alignments = []
            for hit in window["full_reference_hits"]:
                assert hit["query"] == window["query"] and hit["query_length"] == window["window_length"]
                build = {"CHM13": "t2t", "GRCh38": "hg38"}[hit["build"]]
                for (target_build, locus), target in audit.TARGETS.items():
                    if build != target_build or hit["target"] != target[0]:
                        continue
                    if hit["target_start"] >= target[2] or hit["target_end"] <= target[1]:
                        continue
                    alignments.append(graph_paf_alignment(audit, hit, locus))
                matched_hits.append((window, match, hit))
            candidates = audit.best_candidates(window["query"], match["element_start_in_window"],
                match["element_end_in_window"], {window["query"]: alignments})
            for candidate in candidates:
                key = candidate["build"], candidate["locus"]
                if key not in best or candidate["_score"] > best[key]["_score"]:
                    best[key] = candidate
    return sorted(best.values(), key=lambda c: c["_score"], reverse=True), matched_hits, positions


def incorporate_graph_capture(audit, exact_rows, capture_path, catalog_rows, output, reference_authority):
    """Consume the allocated graph's exact-source receipt, never its capture labels."""
    payload = capture_path.read_bytes()
    capture = json.loads(payload)
    assert not capture.get("error")
    assert capture["host"] == capture["allocation"] and capture["host"].startswith("pax")
    assert str(capture["numeric_job"]).isdigit()
    requests = json.loads((output/"omitted_graph_binding_requests.json").read_text())
    requested = {row["ID_Full"]: row for row in requests}
    assert len(requested) == len(requests) == len(capture["acro_rows"])
    assert {row["ID_Full"] for row in capture["acro_rows"]} == set(requested)
    for build, prior_build in (("CHM13", "t2t"), ("GRCh38", "hg38")):
        for key in ("path", "bytes", "sha256"):
            assert capture["references"][build][key] == reference_authority[prior_build][key]
    catalog = {row["ID_Full"]: row for row in catalog_rows}
    old_by_id = {row["ID_Full"]: row for row in exact_rows}
    mapped, reasons, evidence = {}, {}, []
    for row in capture["acro_rows"]:
        request = requested[row["ID_Full"]]
        assert all(row[key] == value for key, value in request.items())
        original = catalog[row["ID_Full"]]
        assert original["Source_Identifier"] == row["Source_Identifier"]
        assert original["orig_Locus"] == row["orig_Locus"]
        retained = Path(row["retained_fasta"]).read_bytes()
        assert sha256(retained).hexdigest() == row["retained_fasta_sha256"]
        lines = retained.decode().splitlines()
        assert sum(line.startswith(">") for line in lines) == 1
        sequence = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
        assert len(sequence) == row["retained_sequence_length"]
        assert sha256(sequence.encode()).hexdigest() == row["retained_sequence_sha256"]
        candidates, hits, positions = graph_source_candidates(audit, row, capture["source_windows"])
        if not positions:
            reason = "exact_native_source_path_not_retained_in_regional_graph" if not row["native_windows"] else "retained_graph_paths_do_not_contain_exact_source_sequence"
        elif len(positions) > 1:
            reason = "multiple_exact_source_positions_within_native_window"
        elif not candidates:
            reason = "full_reference_alignments_do_not_support_named_element_placement"
        else:
            reason = ""
        source = row["Source_Identifier"]
        reasons[source] = reason
        if len(positions) == 1:
            mapped[source] = candidates
            new = dict.fromkeys(exact_rows[0], "")
            new.update(ID_Full=row["ID_Full"], raw_locus=row["orig_Locus"], Source_Identifier=source,
                same_source_loci=row["orig_Locus"], candidate_count=str(len(candidates)),
                all_candidates=";".join(audit.format_candidate(c) for c in candidates))
            if candidates:
                best = candidates[0]
                new.update(best_build=best["build"], best_locus="HML-2_"+best["locus"], best_primary=int(best["primary"]),
                    best_mapq=best["mapq"], best_both5=int(best["both5"]), best_interior100=int(best["interior100"]),
                    best_interior1m=int(best["interior1m"]), best_left_flank=best["left_flank"],
                    best_right_flank=best["right_flank"], best_candidate=audit.format_candidate(best))
            old_by_id[row["ID_Full"]] = new
        evidence.append({"ID_Full": row["ID_Full"], "Source_Identifier": source,
            "retained_fasta": row["retained_fasta"], "retained_sequence_sha256": row["retained_sequence_sha256"],
            "capture_status": row["capture_status"], "unique_exact_sequence_positions": len(positions),
            "exact_native_element_positions": json.dumps(sorted(positions)),
            "native_windows_examined": len(row["native_windows"]), "full_reference_hits_examined": len(hits),
            "all_candidates": ";".join(audit.format_candidate(c) for c in candidates), "source_evidence_limit": reason,
            "allocated_capture": str(capture_path), "allocated_capture_sha256": sha256(payload).hexdigest()})
    write_tsv(output/"graph_source_placement_evidence.tsv", evidence)
    effective = output/"exact_assignment_audit.with_source_remap.tsv"
    write_tsv(effective, list(old_by_id.values()))
    return list(old_by_id.values()), mapped, reasons, effective


def incorporate_source_remap(audit, exact_rows, remap_path, binding_path, catalog_rows, output):
    """Bind byte-verified retained element spans to allocated full-reference CIGARs."""
    remap = json.loads(remap_path.read_text())
    bindings = json.loads(binding_path.read_text())
    assert not remap.get("error") and not bindings.get("error")
    assert remap["requested"] == remap["extracted"] == bindings["validated"] == 72
    assert remap["host"].startswith("pax") and bindings["host"].startswith("pax")
    binding_by_id = {row["ID_Full"]: row for row in bindings["rows"]}
    catalog = {row["ID_Full"]: row for row in catalog_rows}
    old_by_id = {row["ID_Full"]: row for row in exact_rows}
    assert len(binding_by_id) == len(remap["rows"]) == 72
    mapped_sources = {}
    evidence = []
    off_catalog = []
    for row in remap["rows"]:
        bound = binding_by_id[row["ID_Full"]]
        original = catalog[row["ID_Full"]]
        assert bound["mapped_window_unchanged"] and bound["matching_orientation"] in {"+", "-"}
        assert bound["remap_receipt_sha256"] == sha256(remap_path.read_bytes()).hexdigest()
        assert sha256(Path(bound["retained_fasta"]).read_bytes()).hexdigest() == bound["retained_fasta_sha256"]
        assert row["Source_Identifier"] == bound["Source_Identifier"] == original["Source_Identifier"]
        assert row["orig_Locus"] == original["orig_Locus"]
        alignments = [audit.Alignment(**item) for item in row["candidate_alignments"]]
        old_candidates = audit.best_candidates(row["query"], row["source_start_in_window"], row["source_end_in_window"], {row["query"]: alignments})
        assert ";".join(audit.format_candidate(item) for item in old_candidates) == row["all_candidates"]
        start = bound["source_start0"]-row["window_start"]
        end = bound["source_end0"]-row["window_start"]
        candidates = audit.best_candidates(row["query"], start, end, {row["query"]: alignments})
        new = dict.fromkeys(exact_rows[0], "")
        new.update(ID_Full=row["ID_Full"], raw_locus=row["orig_Locus"], Source_Identifier=row["Source_Identifier"],
            same_source_loci=row["orig_Locus"], candidate_count=str(len(candidates)),
            all_candidates=";".join(audit.format_candidate(item) for item in candidates))
        if candidates:
            best = candidates[0]
            new.update(best_build=best["build"], best_locus="HML-2_"+best["locus"], best_primary=int(best["primary"]),
                best_mapq=best["mapq"], best_both5=int(best["both5"]), best_interior100=int(best["interior100"]),
                best_interior1m=int(best["interior1m"]), best_left_flank=best["left_flank"],
                best_right_flank=best["right_flank"], best_candidate=audit.format_candidate(best))
        old_by_id[row["ID_Full"]] = new
        mapped_sources[row["Source_Identifier"]] = candidates
        if not candidates:
            for hit in row["full_reference_hits"]:
                fields = hit["paf"].split("\t")
                tags = {item[:2]: item[5:] for item in fields[12:]}
                if tags.get("tp") != "P" or int(fields[11]) < 20:
                    continue
                strand, reference_position = fields[4], int(fields[7])
                query_position = int(fields[2]) if strand == "+" else int(fields[3])
                spans = []
                for length, operation in audit.parse_cigar(tags["cg"]):
                    if operation in audit.ALIGNED:
                        low, high = (query_position, query_position+length) if strand == "+" else (query_position-length, query_position)
                        a, b = max(start, low), min(end, high)
                        if a < b:
                            spans.append((reference_position+a-low, reference_position+b-low) if strand == "+" else
                                (reference_position+high-b, reference_position+high-a))
                    if operation in audit.REFERENCE_CONSUMING:
                        reference_position += length
                    if operation in audit.QUERY_COORDINATE_CONSUMING:
                        query_position += length if strand == "+" else -length
                covered = sum(b-a for a, b in spans)
                if covered >= 0.8*(end-start):
                    off_catalog.append({"ID_Full": row["ID_Full"], "Source_Identifier": row["Source_Identifier"],
                        "reference": hit["build"], "chromosome_or_reference_contig": fields[5],
                        "reference_start0": min(a for a, b in spans), "reference_end0": max(b for a, b in spans),
                        "MAPQ": int(fields[11]), "aligned_source_bases": covered, "retained_source_bases": end-start,
                        "evidence": "byte_bound_retained_source_and_long_host_window_map_outside_every_named_catalog_element"})
        evidence.append({"ID_Full": row["ID_Full"], "Source_Identifier": row["Source_Identifier"],
            "mapped_source_start0": bound["source_start0"], "mapped_source_end0": bound["source_end0"],
            "retained_fasta": bound["retained_fasta"], "retained_fasta_sha256": bound["retained_fasta_sha256"],
            "source_sequence_sha256": bound["retained_sequence_sha256"], "matching_orientation": bound["matching_orientation"],
            "original_padded_window_candidates": row["all_candidates"], "actual_source_candidates": new["all_candidates"]})
    write_tsv(output/"remapped_source_placement_evidence.tsv", evidence)
    write_tsv(output/"off_catalog_physical_placements.tsv", off_catalog,
        ["ID_Full", "Source_Identifier", "reference", "chromosome_or_reference_contig", "reference_start0", "reference_end0", "MAPQ", "aligned_source_bases", "retained_source_bases", "evidence"])
    effective = output/"exact_assignment_audit.with_source_remap.tsv"
    write_tsv(effective, list(old_by_id.values()))
    return list(old_by_id.values()), mapped_sources, effective


def render_assignment_resolution(changed, unresolved, public_ids, output):
    """Show positive assignments and the remaining evidence limits in one cohort."""
    import matplotlib.pyplot as plt
    resolved = Counter(row["Locus"].removeprefix("HML-2_") for row in changed if row["ID_Full"] in public_ids)
    pending = Counter(row["reason"] for row in unresolved if row["ID_Full"] in public_ids)
    names = {
        "source_contig_absent_from_retained_attachment_capture": "No retained exact-source host linkage",
        "retained_element_alignments_do_not_support_this_source_interval": "Captured alignments do not support\nthe exact source interval",
        "insufficient_primary_mapping_or_bidirectional_host_flanks": "Insufficient unique host-flank support",
        "duplicated_local_block_without_unique_long_host_anchor": "Competing local duplicate-block\nalignments",
        "competing_long_host_anchors": "Competing long host anchors",
        "full_reference_host_remap_places_source_outside_named_catalog_elements": "Mapped outside named catalog loci",
        "exact_native_source_path_not_retained_in_regional_graph": "Exact native source path not captured",
        "retained_graph_paths_do_not_contain_exact_source_sequence": "Captured native paths lack the\nexact source sequence",
        "multiple_exact_source_positions_within_native_window": "Multiple exact native source positions",
        "full_reference_alignments_do_not_support_named_element_placement": "Full-reference alignment does not\nsupport a named element placement",
    }
    data = [{"scope": "newly_named_physical_locus", "category": locus, "public_records": resolved[locus]} for locus in ("13p13", "15p13b", "15p13a", "21p13", "22p13")]
    data += [{"scope": "remaining_placement", "category": reason, "public_records": count} for reason, count in pending.items()]
    write_tsv(output/"Figure_S5_assignment_resolution.tsv", data)
    fig = plt.figure(figsize=(7.1, 4.5))
    groups = [(data[:5], [.44, .56, .52, .35], "A", "Records assigned to named loci"),
        (data[5:], [.44, .105, .52, .35], "B", "Why records remain in copy groups")]
    maximum = max(row["public_records"] for row in data)*1.16
    for records, bounds, panel, title in groups:
        axis = fig.add_axes(bounds)
        y = list(range(len(records)))[::-1]
        labels = [row["category"] if panel == "A" else names[row["category"]] for row in records]
        colors = [("#356F6A" if row["category"] in {"13p13", "15p13b"} else "#765A78") if panel == "A" else "#87939B" for row in records]
        axis.barh(y, [row["public_records"] for row in records], color=colors, height=.62)
        axis.set_yticks(y, labels, fontsize=9.5)
        axis.set_xlim(0, maximum)
        axis.tick_params(axis="x", labelsize=9)
        axis.tick_params(axis="y", length=0)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.spines["bottom"].set_color("#AAB1B6")
        axis.grid(axis="x", color="#E9ECEF", linewidth=.6)
        axis.set_axisbelow(True)
        for position, row in zip(y, records):
            axis.text(row["public_records"]+maximum*.012, position, f"{row['public_records']:,}", va="center", fontsize=9.5)
        axis.set_title(title, fontsize=10, loc="left", pad=7)
        fig.text(.015, bounds[1]+bounds[3]+.025, panel, fontsize=11, weight="bold")
    axis.set_xlabel("Public-cohort catalog records", fontsize=9.5)
    for suffix in ("png", "svg"):
        fig.savefig(output/f"Figure_S5_assignment_resolution.{suffix}", dpi=450, facecolor="white")
    plt.close(fig)


def render_network(edge_rows, output, location_path):
    os.environ["MPLCONFIGDIR"] = str(output / "mplcache")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont

    plt.rcParams.update({"font.family": "Arial", "font.size": 9.5, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(7.1, 3.0))
    fig.text(.014, .988, "B", fontsize=11, weight="bold", va="top")
    groups = [
        ("Acrocentric Type I", [.035, .57, .32, .29], "#356F6A", {"13p13": (-1, 0), "15p13b": (1, 0)}),
        ("Acrocentric Type II", [.44, .49, .52, .40], "#765A78", {"15p13a": (-1, .6), "21p13": (0, -.65), "22p13": (1, .6)}),
        ("1p36.21", [.025, .035, .285, .38], "#336D9E", {"1p36.21a": (-1, .55), "1p36.21b": (0, -.55), "1p36.21c": (1, .55)}),
        ("8p23.1", [.33, .02, .40, .405], "#B96A32", {"8p23.1b": (-1, .65), "8p23.1c": (-1, -.65), "8p23.1d": (1, -.65), "8p23.1e": (1, .65)}),
        ("Xq28", [.76, .06, .22, .25], "#A83F50", {"Xq28a": (-1, 0), "Xq28b": (1, 0)}),
    ]
    plotted = set()
    for title, bounds, color, positions in groups:
        ax = fig.add_axes(bounds)
        ax.set_xlim(-1.65, 1.65)
        ax.set_ylim(-1.2, 1.2)
        ax.axis("off")
        ax.set_title(title, fontsize=9.5, pad=8)
        nodes = {}
        for name, xy in positions.items():
            label = name if title.startswith("Acrocentric") else name[-1]
            nodes[name] = ax.text(*xy, label, ha="center", va="center", fontsize=9.5,
                bbox={"boxstyle": "round,pad=.22", "fc": "white", "ec": color, "lw": 1.5}, zorder=5).get_bbox_patch()
        for row in edge_rows:
            a, b = row["locus_1"], row["locus_2"]
            if a not in positions or b not in positions:
                continue
            plotted.add((a, b))
            p, q = np.array(positions[a]), np.array(positions[b])
            rad = .24 if p[0]*q[1] - p[1]*q[0] >= 0 else -.24
            diagonal = title == "8p23.1" and abs(p[0]-q[0]) > 1 and abs(p[1]-q[1]) > 1
            if diagonal:
                rad = .12 if a.endswith("b") else -.12
            patch = FancyArrowPatch(p, q, patchA=nodes[a], patchB=nodes[b], arrowstyle="-",
                connectionstyle=f"arc3,rad={rad}", color=color, lw=2, zorder=1)
            ax.add_patch(patch)
            vertices = patch.get_path().vertices
            assert len(vertices) == 3
            t = (.28 if a.endswith("b") else .68) if diagonal else .5
            middle = (1-t)**2*vertices[0] + 2*(1-t)*t*vertices[1] + t*t*vertices[2]
            ax.text(*middle, str(row["distinct_gene_sequences"]), ha="center", va="center", fontsize=9.5,
                bbox={"fc": "white", "ec": "none", "pad": .6}, zorder=4)
    assert plotted == {(r["locus_1"], r["locus_2"]) for r in edge_rows}
    for label, name in (("B", "Figure_3B_nucleotide_sharing"), ("C", "Figure_3C_nucleotide_sharing"), ("D", "Figure_S12D_nucleotide_sharing")):
        fig.texts[0].set_text(label)
        fig.savefig(output / f"{name}.png", dpi=450, facecolor="white")
        fig.savefig(output / f"{name}.svg", facecolor="white")
        fig.savefig(output / f"{name}.pdf", facecolor="white")
    plt.close(fig)

    network = Image.open(output / "Figure_3C_nucleotide_sharing.png").convert("RGB")
    locations = Image.open(location_path).convert("RGB")
    locations = locations.crop(Image.eval(locations.convert("L"), lambda x: 255-x).getbbox())
    width = 4200
    network = network.resize((width, round(width*network.height/network.width)), Image.Resampling.LANCZOS)
    scale = min(3650/locations.width, 2180/locations.height)
    locations = locations.resize((round(scale*locations.width), round(scale*locations.height)), Image.Resampling.LANCZOS)
    combined = Image.new("RGB", (width, network.height+110+locations.height+50), "white")
    combined.paste(network, (0, 0))
    combined.paste(locations, ((width-locations.width)//2, network.height+110))
    draw = ImageDraw.Draw(combined)
    font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 105)
    draw.text((45, network.height+110), "D", fill="black", font=font)
    combined.save(output / "Figure_3CD.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-catalog", type=Path, required=True)
    parser.add_argument("--prior-nucleotide-dir", type=Path, required=True)
    parser.add_argument("--nucleotide-builder", type=Path, required=True)
    parser.add_argument("--reuse-nucleotide-dir", type=Path)
    parser.add_argument("--source-remap", type=Path)
    parser.add_argument("--source-bindings", type=Path)
    parser.add_argument("--graph-capture", type=Path)
    parser.add_argument("--public-catalog", type=Path,
        help="Final merged catalog whose public inclusion/exclusion flags govern the figures; never edited here")
    parser.add_argument("--location-figure", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    bio = load_module("biological_orf", BIOLOGICAL / "build_biologically_annotated_orf.py")
    audit = load_module("exact_assignment_audit", BIOLOGICAL / "audit_exact_duplicate_locus_assignments.py")
    rows = read_tsv(args.source_catalog)
    original_rows = [dict(row) for row in rows]
    exact_rows = read_tsv(EXACT_AUDIT)
    alignments = audit.load_alignments(ATTACHMENTS)
    replayed = 0
    for row in exact_rows:
        if row["raw_locus"] not in bio.ACRO_FAMILY:
            continue
        qname, start, end = audit.source_interval(row["Source_Identifier"])
        candidates = audit.best_candidates(qname, start, end, alignments)
        assert ";".join(audit.format_candidate(c) for c in candidates) == row["all_candidates"], row["ID_Full"]
        replayed += 1
    effective_audit = EXACT_AUDIT
    remapped_sources = {}
    assert bool(args.source_remap) == bool(args.source_bindings)
    if args.source_remap:
        exact_rows, remapped_sources, effective_audit = incorporate_source_remap(
            audit, exact_rows, args.source_remap, args.source_bindings, rows, output)
    graph_reasons = {}
    graph_mapped_sources = {}
    if args.graph_capture:
        assert args.source_remap, "Full-reference identity authority is required"
        exact_rows, graph_mapped_sources, graph_reasons, effective_audit = incorporate_graph_capture(
            audit, exact_rows, args.graph_capture, rows, output, json.loads(args.source_remap.read_text())["references"])
        assert not (remapped_sources.keys() & graph_mapped_sources.keys())
        remapped_sources.update(graph_mapped_sources)
    by_id = {r["ID_Full"]: r for r in exact_rows}
    by_source = {r["Source_Identifier"]: r for r in exact_rows}
    decisions = bio.load_exact_source_decisions(effective_audit)

    changed = []
    mutable_fields = ("Locus", "physical_locus_assignment", "physical_assignment_evidence")
    for row in rows:
        if row["Locus"] not in {"HML-2_acro_type1", "HML-2_acro_type2"}:
            continue
        decision = decisions.get(row["Source_Identifier"])
        if decision is None or decision[0] == row["Locus"]:
            continue
        support = by_id[row["ID_Full"]]
        assert support["Source_Identifier"] == row["Source_Identifier"]
        assert support["raw_locus"] == row["orig_Locus"]
        raw = dict(row, Locus=row["orig_Locus"])
        annotated = bio.annotate_row(raw, set(), exact_source_decisions=decisions)
        assert annotated["insertion_event_group"] == row["insertion_event_group"]
        delta = {key: row[key] for key in ("ID_Full", "Source_Identifier", "orig_Locus", "analysis_include")}
        delta["prior_Locus"] = row["Locus"]
        for field in mutable_fields:
            row[field] = annotated[field]
            delta[field] = row[field]
        delta["all_candidates"] = support["all_candidates"]
        changed.append(delta)
    destination = output / "combined_hml2_orf_analysis.ACROCENTRIC_RESOLVED.tsv"
    write_tsv(destination, rows)
    write_tsv(output / "acrocentric_assignment_changes.tsv", changed)
    changed_columns = {key for old, new in zip(original_rows, rows) for key in old if old[key] != new[key]}
    assert changed_columns <= set(mutable_fields), changed_columns
    assert len(rows) == len(original_rows)
    sys.path.insert(0, str(PROJECT / "manuscript"))
    manuscript = load_module("current_manuscript_catalog_consumer", PROJECT / "manuscript/build_narrative_main_figures.py")
    manuscript.CATALOG = args.public_catalog or destination
    public_rows, public_roster = manuscript.load_catalog()
    public_ids = {row["ID_Full"] for row in public_rows}
    cohort_catalog = read_tsv(args.public_catalog) if args.public_catalog else rows
    assert len(cohort_catalog) == len(rows)
    original_public_ids = set()
    for original, current in zip(rows, cohort_catalog):
        assert tuple(original[key] for key in ("ID", "Haplotype", "orig_Locus")) == tuple(
            current[key] for key in ("ID", "Haplotype", "orig_Locus"))
        if current["ID_Full"] in public_ids:
            original_public_ids.add(original["ID_Full"])
    public_exclusions = Counter(row["analysis_exclusion_reason"] for row in cohort_catalog
        if manuscript.PUBLIC_ID.fullmatch(row["ID"]) and row["analysis_include"] == "0")

    unresolved = []
    eligible_unconsumed = []
    for row in rows:
        if row["Locus"] not in {"HML-2_acro_type1", "HML-2_acro_type2"} or row["analysis_include"] != "1":
            continue
        qname, start, end = audit.source_interval(row["Source_Identifier"])
        candidates = remapped_sources.get(row["Source_Identifier"])
        if candidates is None:
            candidates = audit.best_candidates(qname, start, end, alignments)
        if row["Source_Identifier"] not in by_source and candidates:
            eligible_unconsumed.append(row["ID_Full"])
        anchors = [c for c in candidates if c["primary"] and c["mapq"] >= 20]
        if graph_reasons.get(row["Source_Identifier"]):
            reason = graph_reasons[row["Source_Identifier"]]
        elif row["Source_Identifier"] in remapped_sources and not candidates:
            reason = "full_reference_host_remap_places_source_outside_named_catalog_elements"
        elif not alignments.get(qname) and row["Source_Identifier"] not in remapped_sources:
            reason = "source_contig_absent_from_retained_attachment_capture"
        elif not candidates:
            reason = "retained_element_alignments_do_not_support_this_source_interval"
        elif len([c for c in anchors if c["interior100"]]) > 1:
            reason = "competing_long_host_anchors"
        elif len(candidates) > 1 and any(c["both5"] for c in anchors):
            reason = "duplicated_local_block_without_unique_long_host_anchor"
        else:
            reason = "insufficient_primary_mapping_or_bidirectional_host_flanks"
        unresolved.append({"ID_Full": row["ID_Full"], "Source_Identifier": row["Source_Identifier"],
            "orig_Locus": row["orig_Locus"], "Locus": row["Locus"], "Structure": row["Structure"],
            "reason": reason, "all_candidates": ";".join(audit.format_candidate(c) for c in candidates)})
    assert not eligible_unconsumed, eligible_unconsumed
    write_tsv(output / "remaining_acrocentric_placement_evidence.tsv", unresolved)

    nucleotide = load_module("retained_nucleotide_analysis", args.nucleotide_builder)
    nucleotide.CATALOG = args.public_catalog or destination
    nucleotide.OUT = output / "nucleotide_analysis"
    if args.reuse_nucleotide_dir:
        nucleotide.OUT.mkdir(parents=True, exist_ok=True)
        by_catalog_id = {r["ID_Full"]: r for r in rows}
        reused_sources = read_tsv(args.reuse_nucleotide_dir / "source_records.tsv")
        reused_observations = read_tsv(args.reuse_nucleotide_dir / "gene_sequence_observations.tsv")
        for record in reused_sources:
            catalog_row = by_catalog_id[record["ID_Full"]]
            assert record["Source_Identifier"] == catalog_row["Source_Identifier"]
            assert sha256(Path(record["path"]).read_bytes()).hexdigest() == record["source_sha256"]
            for key in ("Locus", "physical_locus_assignment"):
                record[key] = catalog_row[key]
        for observation in reused_observations:
            catalog_row = by_catalog_id[observation["ID_Full"]]
            assert observation["Source_Identifier"] == catalog_row["Source_Identifier"]
            observation["locus"] = catalog_row["Locus"].removeprefix("HML-2_")
        write_tsv(nucleotide.OUT / "source_records.tsv", reused_sources)
        write_tsv(nucleotide.OUT / "gene_sequence_observations.tsv", reused_observations)
        reused_summary = json.loads((args.reuse_nucleotide_dir / "summary.json").read_text())
        reused_summary["catalog_sha256"] = sha256(destination.read_bytes()).hexdigest()
        reused_summary["reused_byte_verified_nucleotide_evidence"] = str(args.reuse_nucleotide_dir)
        reused_summary["genes_by_locus"] = dict(Counter(r["locus"] for r in reused_observations))
        raw_loci = defaultdict(set)
        raw_observations = defaultdict(list)
        for observation in reused_observations:
            raw_loci[(observation["gene"], observation["sha256"])].add(observation["locus"])
            raw_observations[(observation["gene"], observation["sha256"])].append(observation)
        raw_edges = Counter(pair for loci in raw_loci.values() for pair in combinations(sorted(loci), 2))
        reused_summary["edges"] = [{"locus_1": a, "locus_2": b, "distinct_gene_sequences": n} for (a, b), n in sorted(raw_edges.items())]
        write_tsv(nucleotide.OUT / "exact_nucleotide_edges.tsv", reused_summary["edges"])
        raw_witnesses = []
        for key, members in raw_observations.items():
            for a, b in combinations(sorted(raw_loci[key]), 2):
                raw_witnesses.append({"locus_1": a, "locus_2": b, "gene": key[0], "length": members[0]["length"],
                    "sha256": key[1], "witnesses": json.dumps([(r["locus"], r["ID_Full"], r["Source_Identifier"]) for r in members])})
        write_tsv(nucleotide.OUT / "edge_witnesses.tsv", raw_witnesses)
        (nucleotide.OUT / "summary.json").write_text(json.dumps(reused_summary, indent=2)+"\n")
    else:
        nucleotide.main()
    observations = read_tsv(nucleotide.OUT / "gene_sequence_observations.tsv")
    previous = read_tsv(args.prior_nucleotide_dir / "gene_sequence_observations.tsv")
    identity = lambda r: tuple(r[k] for k in ("ID_Full", "Source_Identifier", "gene", "length", "sha256"))
    current_identities = Counter(map(identity, observations))
    previous_identities = Counter(map(identity, previous))
    assert not (previous_identities-current_identities), 'Previously verified DNA sequences changed'
    current_by_id = {r["ID_Full"]: r for r in cohort_catalog}
    added_identities = current_identities-previous_identities
    assert all(current_by_id[key[0]].get("v3_row_origin") == "BROAD_GRAPH_SOURCE_RECOVERY_20260914"
               for key in added_identities), 'New gene observations lack recovered-source evidence'
    sources = read_tsv(nucleotide.OUT / "source_records.tsv")
    for source in sources:
        assert sha256(Path(source["path"]).read_bytes()).hexdigest() == source["source_sha256"]

    by_sequence = defaultdict(lambda: defaultdict(list))
    for row in observations:
        if row["locus"] not in POOLS and row["ID_Full"] in public_ids:
            by_sequence[(row["gene"], row["sha256"])][row["locus"]].append(row)
    edges = defaultdict(set)
    witnesses = []
    for sequence, loci in by_sequence.items():
        for a, b in combinations(sorted(loci), 2):
            left = {r["Source_Identifier"] for r in loci[a]}
            right = {r["Source_Identifier"] for r in loci[b]}
            assert left.isdisjoint(right), (a, b, left & right)
            for left_source in left:
                left_contig, left_start, left_end = audit.source_interval(left_source)
                for right_source in right:
                    right_contig, right_start, right_end = audit.source_interval(right_source)
                    assert left_contig != right_contig or max(left_start, right_start) >= min(left_end, right_end), (a, b, left_source, right_source)
            edges[(a, b)].add(sequence)
            witnesses.append({"locus_1": a, "locus_2": b, "gene": sequence[0], "sha256": sequence[1],
                "length": loci[a][0]["length"], "witnesses": json.dumps([(r["locus"], r["ID_Full"], r["Source_Identifier"]) for r in loci[a]+loci[b]])})
    edge_rows = [{"locus_1": a, "locus_2": b, "distinct_gene_sequences": len(sequences)} for (a, b), sequences in sorted(edges.items())]
    write_tsv(output / "Figure_3C_exact_nucleotide_edges.tsv", edge_rows)
    write_tsv(output / "Figure_3C_sequence_witnesses.tsv", witnesses)
    render_network(edge_rows, output, args.location_figure)
    render_assignment_resolution(changed, unresolved, original_public_ids, output)
    summary = {
        "input_catalog": str(args.source_catalog), "input_catalog_sha256": sha256(args.source_catalog.read_bytes()).hexdigest(),
        "output_catalog": str(destination), "output_catalog_sha256": sha256(destination.read_bytes()).hexdigest(),
        "exact_assignment_audit": str(EXACT_AUDIT), "exact_assignment_audit_sha256": sha256(EXACT_AUDIT.read_bytes()).hexdigest(),
        "effective_exact_assignment_audit": str(effective_audit), "effective_exact_assignment_audit_sha256": sha256(effective_audit.read_bytes()).hexdigest(),
        "allocated_source_remap": str(args.source_remap) if args.source_remap else None,
        "allocated_source_bindings": str(args.source_bindings) if args.source_bindings else None,
        "current_assembly_source_sequences_verified": len(remapped_sources)-len(graph_mapped_sources),
        "allocated_graph_capture": str(args.graph_capture) if args.graph_capture else None,
        "graph_source_sequences_uniquely_bound": len(graph_mapped_sources),
        "graph_capture_source_status": dict(Counter(reason or "exact_source_bound_candidate_set_available" for reason in graph_reasons.values())),
        "attachment_capture": str(ATTACHMENTS), "attachment_capture_sha256": sha256(ATTACHMENTS.read_bytes()).hexdigest(),
        "acrocentric_audit_rows_reproduced_from_CIGAR": replayed,
        "assignment_changes": len(changed), "changes_by_locus": dict(Counter(r["Locus"] for r in changed)),
        "changed_catalog_columns": list(mutable_fields), "included_rows_unchanged": sum(r["analysis_include"] == "1" for r in rows),
        "public_cohort_catalog": str(manuscript.CATALOG), "public_cohort_catalog_sha256": sha256(manuscript.CATALOG.read_bytes()).hexdigest(),
        "public_included_rows_unchanged": len(public_rows), "public_haplotypes_unchanged": len(public_roster),
        "public_included_assignment_repairs": sum(row["ID_Full"] in original_public_ids for row in changed),
        "public_exclusions_unchanged": dict(public_exclusions),
        "retained_pool_rows": len(unresolved), "retained_pool_reasons": dict(Counter(r["reason"] for r in unresolved)),
        "public_retained_pool_rows": sum(row["ID_Full"] in original_public_ids for row in unresolved),
        "public_retained_pool_reasons": dict(Counter(row["reason"] for row in unresolved if row["ID_Full"] in original_public_ids)),
        "eligible_retained_source_evidence_omitted_from_old_audit": eligible_unconsumed,
        "verified_source_FASTAs": len(sources), "DNA_observations_reproduced_unchanged": len(previous),
        "additional_recovered_source_gene_observations": sum(added_identities.values()),
        "excluded_pool_gene_observations": sum(r["locus"] in POOLS for r in observations),
        "excluded_nonprimary_cohort_gene_observations": sum(r["ID_Full"] not in public_ids for r in observations),
        "physical_locus_edges": len(edge_rows), "edges": edge_rows,
        "edge_endpoint_same_source_or_overlapping_interval_aliases": 0,
        "scope": "Only pooled acrocentric placement fields changed; biological calls, ORFs, CNV weights, and duplicate exclusions were preserved. Unresolved records are evidence-limited placements, not biological absences or terminal non-identifiability."
    }
    (output / "verification.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
