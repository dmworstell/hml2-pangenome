#!/usr/bin/env python3
"""Plot retained ORF annotations and sequence-level longest open frames.

The annotation tier pools Intact and Intact_FS_End, as in the main ORF figure.
It is not an observation of a translated polyprotein. The independent sequence
scan finds the longest ATG-starting open frame in each retained provirus.
An ATG inside the gag interval is not necessarily its canonical ATG.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import re

PUBLIC_ID = re.compile(r"^(?:HG|NA)\d+$")
PASSING = {"Intact", "Intact_FS_End"}
GENES = ("gag", "pro", "pol", "env")
FEATURES = {
    "type1": {"gag": (1111, 3112), "pro": (2913, 3918),
              "pol": (3878, 6501), "env": (6512, 8258)},
    "type2": {"gag": (1111, 3112), "pro": (2913, 3918),
              "pol": (3878, 6749), "env": (6450, 8550)},
}
TIER_ORDER = ("Gag", "Gag + Pro", "Gag + Pro + Pol",
              "Gag + Pro + Pol + Env (Type I)")
REACH_ORDER = ("gag only", "gag + pro", "gag + pro + pol")
SCAN_FIELDS = ("ID_Full", "Locus", "orig_Locus", "ID", "Haplotype", "type",
               "longest_ATG_ORF_aa", "start_gene", "genes_spanned",
               "start_nt0", "end_nt0", "start_KCON_nt0", "frame0",
               "end_reason", "sequence_bases", "noncanonical_bases", "alignment_path")


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_catalog(path):
    """Apply the main figure's public-ID and analysis_include rules."""
    rows, rejected = [], collections.Counter()
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not PUBLIC_ID.fullmatch(row["ID"]):
                rejected["reference_or_nonpublic_ID"] += 1
            elif row["analysis_include"] == "0":
                rejected[row["analysis_exclusion_reason"]] += 1
            elif row["analysis_include"] == "1":
                rows.append(row)
            else:
                raise ValueError("Unrecognized inclusion state: " + row["ID_Full"])
    if len({row["ID_Full"] for row in rows}) != len(rows):
        raise ValueError("Repeated retained catalog record identity")
    return rows, rejected


def annotation_tiers(rows):
    counts, witnesses = collections.Counter(), []
    for row in rows:
        if "Provirus" not in row["Structure"] or row["gag"] not in PASSING:
            continue
        tier = 0
        if row["pro"] in PASSING:
            tier = 1
            if row["pol"] in PASSING:
                tier = 2
                if row["provirus_type"] == "type1" and row["env"] in PASSING:
                    tier = 3
        counts[TIER_ORDER[tier]] += 1
        witnesses.append({key: row[key] for key in
                          ("ID_Full", "Locus", "ID", "Haplotype", "gag", "pro", "pol", "env", "provirus_type")}
                         | {"highest_tier": TIER_ORDER[tier]})
    return counts, witnesses


def alignment_request(rows):
    """Bind current retained proviruses to their original per-locus producer."""
    requests = []
    for row in rows:
        if "Provirus" not in row["Structure"]:
            continue
        original_path = row["row_source_path"]
        alignment = re.sub(r"_orf_integrity_results_(type[12]_KCON)(?:_completed)?\.csv$",
                           r"_all_haplotypes_aligned_\1.fasta", original_path)
        if alignment == original_path:
            raise ValueError("No original KCON alignment binding: " + row["ID_Full"])
        requests.append({key: row[key] for key in
                         ("ID_Full", "Locus", "orig_Locus", "ID", "Haplotype", "provirus_type")}
                        | {"alignment_path": alignment})
    return requests


def longest_open_frame(sequence):
    """Keep the longest ATG-to-stop/boundary segment, preserving real indels.

    Ambiguous codons terminate a segment instead of being deleted. Ties retain
    the first frame and first ATG encountered, as in the historical scan.
    """
    best = (0, 0, 0, -1, "no_ATG")
    stop_codons = {"TAA", "TAG", "TGA"}
    for frame in range(3):
        start = None
        end = frame
        for offset in range(frame, len(sequence) - 2, 3):
            codon = sequence[offset:offset + 3]
            end = offset + 3
            ambiguous = any(base not in "ACGT" for base in codon)
            if codon in stop_codons or ambiguous:
                if start is not None and (offset - start) // 3 > best[0]:
                    best = ((offset - start) // 3, start, offset, frame,
                            "ambiguous_codon" if ambiguous else "stop_codon")
                start = None
            elif start is None and codon == "ATG":
                start = offset
        if start is not None and (end - start) // 3 > best[0]:
            best = ((end - start) // 3, start, end, frame, "sequence_end")
    return best


def gene_at(position, features):
    # Shared overlap bases belong to the first (upstream) interval for the
    # regional-start summary, matching the historical scan's convention.
    for gene, (start, end) in features.items():
        if start <= position < end:
            return gene
    return "-"


def scan_alignment_records(requests, alignment_root=None, min_aa=200):
    """Scan true gapped MSAs. Do not supply fixed-width KCON projections."""
    from Bio import SeqIO

    grouped = collections.defaultdict(list)
    for index, request in enumerate(requests):
        grouped[request["alignment_path"]].append((index, request))
    result = {"input_alignments": [], "rows": [], "missing": [],
              "below_minimum": [], "minimum_ORF_aa": min_aa}
    for original_path, wanted in sorted(grouped.items()):
        path = Path(original_path)
        if alignment_root is not None:
            path = Path(alignment_root) / path.parent.name / path.name
        if not path.is_file():
            result["missing"].extend([index, "alignment_file_absent"] for index, _ in wanted)
            continue
        with path.open() as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        types = {request["provirus_type"] for _, request in wanted}
        if len(types) != 1:
            raise ValueError("Conflicting type routing: " + original_path)
        typ = next(iter(types))
        references = [record for record in records if record.id == typ + "_KCON"]
        if len(references) != 1:
            raise ValueError("No unique type-correct reference in " + str(path))
        kcon = str(references[0].seq).upper()
        if len({len(record.seq) for record in records}) != 1:
            raise ValueError("Unequal aligned sequence lengths: " + str(path))
        kcon_positions, position = [], 0
        for base in kcon:
            kcon_positions.append(position)
            if base != "-":
                position += 1
        lookup = {}
        for record in records:
            identity = record.id.split("|", 1)[0]
            if identity in lookup:
                raise ValueError("Duplicated alignment record: " + identity)
            lookup[identity] = record
        result["input_alignments"].append({
            "original_path": original_path, "read_path": str(path),
            "sha256": file_sha256(path), "bytes": path.stat().st_size,
            "reference_id": references[0].id, "aligned_columns": len(kcon),
            "reference_bases": position,
            "reference_sequence_sha256": hashlib.sha256(kcon.replace("-", "").encode()).hexdigest(),
            "requested_records": len(wanted), "alignment_records": len(records)})
        for index, request in wanted:
            record = lookup.get(request["ID_Full"])
            if record is None:
                result["missing"].append([index, "exact_ID_absent_from_alignment"])
                continue
            aligned = str(record.seq).upper()
            columns = [column for column, base in enumerate(aligned) if base != "-"]
            sequence = aligned.replace("-", "")
            length, start, end, frame, reason = longest_open_frame(sequence)
            if length < min_aa:
                result["below_minimum"].append([index, length])
                continue
            start_kcon = kcon_positions[columns[start]]
            start_gene = gene_at(start_kcon, FEATURES[typ])
            genes = {gene_at(kcon_positions[columns[offset]], FEATURES[typ])
                     for offset in range(start, end)} - {"-"}
            mask = sum(1 << GENES.index(gene) for gene in genes)
            # Compact, lossless result rows keep allocated transfer bounded.
            result["rows"].append([index, length, GENES.index(start_gene) if start_gene in GENES else -1,
                mask, start, end, start_kcon, frame,
                ("stop_codon", "sequence_end", "ambiguous_codon").index(reason),
                len(sequence), sum(base not in "ACGT" for base in sequence)])
    result["rows"].sort()
    indexes = [row[0] for row in result["rows"] + result["missing"] + result["below_minimum"]]
    if sorted(indexes) != list(range(len(requests))):
        raise ValueError("Incomplete or duplicated scan accounting")
    return result


def expand_scan_rows(requests, result):
    rows = []
    for (index, length, start_code, mask, start, end, start_kcon, frame,
         end_code, bases, ambiguous) in result["rows"]:
        request = requests[index]
        start_gene = GENES[start_code] if start_code >= 0 else "-"
        rows.append({key: request[key] for key in
                     ("ID_Full", "Locus", "orig_Locus", "ID", "Haplotype", "alignment_path")}
                    | {"type": request["provirus_type"], "longest_ATG_ORF_aa": length,
                       "start_gene": start_gene,
                       "genes_spanned": "-".join(gene for i, gene in enumerate(GENES) if mask & (1 << i)) or "-",
                       "start_nt0": start, "end_nt0": end, "start_KCON_nt0": start_kcon,
                       "frame0": frame,
                       "end_reason": ("stop_codon", "sequence_end", "ambiguous_codon")[end_code],
                       "sequence_bases": bases, "noncanonical_bases": ambiguous})
    return rows


def matched_scan_rows(scan_path, retained):
    lookup = {row["ID_Full"]: row for row in retained if "Provirus" in row["Structure"]}
    used, rejected, accepted = set(), collections.Counter(), []
    with Path(scan_path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            identity = row.get("ID_Full", row.get("provirus", "")).split("|", 1)[0]
            current = lookup.get(identity)
            if current is None:
                rejected["outside_retained_public_proviruses"] += 1
                continue
            if identity in used:
                raise ValueError("Duplicate scan record: " + identity)
            if row.get("type") != current["provirus_type"]:
                raise ValueError("Scan/catalog type mismatch: " + identity)
            used.add(identity)
            row.update({key: current[key] for key in ("ID_Full", "Locus", "ID", "Haplotype")})
            accepted.append(row)
    return accepted, rejected, sorted(set(lookup) - used)


def regional_start_counts(rows):
    counts, witnesses = collections.Counter(), []
    for row in rows:
        if row["start_gene"] != "gag":
            continue
        genes = set(row["genes_spanned"].split("-")) - {"", "-"}
        if genes <= {"gag"}:
            category = REACH_ORDER[0]
        elif genes <= {"gag", "pro"}:
            category = REACH_ORDER[1]
        elif genes == {"gag", "pro", "pol"}:
            category = REACH_ORDER[2]
        else:
            raise ValueError("Additional gene-region combination requires a plotted category: " + row["ID_Full"])
        counts[category] += 1
        witnesses.append(row | {"regional_reach": category})
    return counts, witnesses


def write_tsv(path, rows, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def render_figure(tiers, reach, out, stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({"font.family": "Arial", "font.size": 10,
                               "axes.titlesize": 10.5, "axes.labelsize": 10,
                               "xtick.labelsize": 10, "ytick.labelsize": 10,
                               "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.3, 3.65), layout="constrained")
    palette = ("#858783", "#829EA3", "#506573", "#98664F")
    values = [tiers[tier] for tier in TIER_ORDER]
    left.barh(range(4), values, color=palette, edgecolor="black", linewidth=0.5)
    left.set_yticks(range(4), ["Gag", "Gag + Pro", "Gag + Pro + Pol", "Gag + Pro + Pol + Env\n(Type I)"])
    left.invert_yaxis()
    left.set_xlim(0, max(values, default=1) * 1.19)
    left.set_xlabel("Retained proviral copies")
    left.set_title("A   Highest passing ORF tier", loc="left", pad=10)
    for index, value in enumerate(values):
        left.text(value + max(values, default=1) * .018, index, f"{value:,}", va="center")
    values = [reach[key] for key in REACH_ORDER]
    right.bar(range(3), values, color=(palette[0], palette[1], palette[2]), edgecolor="black", linewidth=.5)
    right.set_xticks(range(3), ["gag only", "gag + pro", "gag + pro + pol"])
    right.set_ylim(0, max(values, default=1) * 1.15)
    right.set_ylabel("Retained proviral copies")
    right.set_title("B   Longest ORF starting within gag", loc="left", pad=10)
    for index, value in enumerate(values):
        right.text(index, value + max(values, default=1) * .02, f"{value:,}", ha="center", va="bottom")
    for axis in (left, right):
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(length=3, width=.6)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(out / (stem + "." + extension), dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tsv", required=True, type=Path)
    parser.add_argument("--scan", required=True, type=Path)
    parser.add_argument("-o", "--out_dir", type=Path, default=Path("."))
    parser.add_argument("--stem", default="hml2_orf_tiers_and_longest_frames")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    retained, excluded = read_catalog(args.tsv)
    tiers, tier_rows = annotation_tiers(retained)
    matched, scan_excluded, unscanned = matched_scan_rows(args.scan, retained)
    reach, reach_rows = regional_start_counts(matched)
    write_tsv(args.out_dir / "highest_annotation_tiers.tsv", tier_rows, list(tier_rows[0]))
    write_tsv(args.out_dir / "gag_region_longest_frame_witnesses.tsv", reach_rows, (*SCAN_FIELDS, "regional_reach"))
    write_tsv(args.out_dir / "catalog_proviruses_without_reported_scan.tsv",
              [{"ID_Full": identity} for identity in unscanned], ("ID_Full",))
    summary = {"catalog": str(args.tsv), "catalog_sha256": file_sha256(args.tsv),
        "scan": str(args.scan), "scan_sha256": file_sha256(args.scan),
        "retained_public_rows": len(retained),
        "public_haplotypes": len({(row["ID"], row["Haplotype"]) for row in retained}),
        "retained_proviruses": sum("Provirus" in row["Structure"] for row in retained),
        "catalog_exclusions": dict(excluded), "annotation_tiers": dict(tiers),
        "reported_scan_rows": len(matched), "scan_exclusions": dict(scan_excluded),
        "retained_proviruses_without_reported_scan": len(unscanned),
        "longest_frame_start_regions": dict(collections.Counter(row["start_gene"] for row in matched)),
        "gag_region_reach": dict(reach),
        "gag_region_rows_with_start_at_KCON_1111": sum(int(row.get("start_KCON_nt0", -1)) == 1111 for row in reach_rows),
        "gag_region_end_reasons": dict(collections.Counter(row.get("end_reason", "not_recorded") for row in reach_rows)),
        "annotation_rule": "Mutually exclusive highest tier, pooling Intact and Intact_FS_End. The fourth tier additionally requires Type I and a passing Env annotation.",
        "scan_rule": "Longest single ATG-starting open frame, not necessarily beginning at the canonical gene start. Coordinates map its regional overlap, not the identity of a translated protein."}
    (args.out_dir / "figure_counts_and_provenance.json").write_text(json.dumps(summary, indent=2) + "\n")
    render_figure(tiers, reach, args.out_dir, args.stem)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
