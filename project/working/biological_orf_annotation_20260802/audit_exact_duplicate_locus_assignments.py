#!/usr/bin/env python3
"""Map duplicated-locus ORF intervals through exact alignment CIGAR strings."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


PRESENT = {"Provirus", "Solo-LTR", "Fragment", "Provirus_from_Multi"}
csv.field_size_limit(sys.maxsize)
DUPLICATED_LOCI = {
    "HML-2_1p36.21a",
    "HML-2_1p36.21b",
    "HML-2_1p36.21c",
    "HML-2_8p23.1a_hg38",
    "HML-2_8p23.1b",
    "HML-2_8p23.1c",
    "HML-2_8p23.1d",
    "HML-2_8p23.1e",
    "HML-2_Xq28a",
    "HML-2_Xq28b",
    "HML-2_13p13",
    "HML-2_15p13a",
    "HML-2_15p13b",
    "HML-2_21p13",
    "HML-2_22p13",
    "HML-2_4q35.2_hg38",
    "HML-2_Yq11.23a",
    "HML-2_Yq11.23b",
}
TARGETS = {
    ("t2t", "13p13"): ("chr13", 5_381_952, 5_388_076, "right"),
    ("t2t", "15p13a"): ("chr15", 12_216, 19_499, "right"),
    ("t2t", "15p13b"): ("chr15", 2_092_086, 2_101_275, "right"),
    ("t2t", "21p13"): ("chr21", 9_160, 16_435, "right"),
    ("t2t", "22p13"): ("chr22", 10_726, 18_005, "right"),
    ("t2t", "1p36.21a"): ("chr1", 12_324_389, 12_330_209, "right"),
    ("t2t", "1p36.21b"): ("chr1", 12_529_193, 12_538_734, "right"),
    ("t2t", "1p36.21c"): ("chr1", 12_793_963, 12_803_479, "right"),
    ("t2t", "8p23.1b"): ("chr8", 7_508_401, 7_517_950, "right"),
    ("t2t", "8p23.1c"): ("chr8", 11_536_938, 11_546_458, "right"),
    ("t2t", "8p23.1d"): ("chr8", 12_201_792, 12_211_329, "right"),
    ("t2t", "8p23.1e"): ("chr8", 12_738_560, 12_748_099, "right"),
    ("t2t", "Xq28a"): ("chrX", 152_825_089, 152_832_431, "left"),
    ("t2t", "Xq28b"): ("chrX", 152_849_562, 152_852_185, "left"),
    ("hg38", "8p23.1a_hg38"): ("chr8", 7_497_875, 7_507_337, "right"),
    ("hg38", "4q35.2_hg38"): ("chr4", 190_106_259, 190_113_546, "left"),
}
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
DETAIL_RE = re.compile(
    r"([^;|]+)\|([^:|]+):(\d+)-(\d+)\|q:(\d+)-(\d+)/(\d+)"
    r"\|([+-])\|mapq:(\d+)\|(nonsecondary|secondary)\|cigar:([^;]+)"
)
REFERENCE_CONSUMING = frozenset("MDN=X")
QUERY_COORDINATE_CONSUMING = frozenset("MI=XS H".replace(" ", ""))
ALIGNED = frozenset("M=X")


@dataclass(frozen=True)
class Alignment:
    build: str
    locus: str
    qname: str
    chrom: str
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    query_length: int
    strand: str
    mapq: int
    primary: bool
    cigar: str


def source_interval(value: str) -> tuple[str, int, int] | None:
    match = re.match(r"^(.*):(\d+)-(\d+)$", value or "")
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3))


@lru_cache(maxsize=None)
def parse_cigar(cigar: str) -> tuple[tuple[int, str], ...]:
    operations = tuple(
        (int(length), operation) for length, operation in CIGAR_RE.findall(cigar)
    )
    if "".join(f"{length}{operation}" for length, operation in operations) != cigar:
        raise ValueError(f"unsupported CIGAR: {cigar}")
    return operations


@lru_cache(maxsize=None)
def aligned_reference_bases(
    alignment: Alignment, window_start: int, window_end: int
) -> int:
    ref = alignment.ref_start
    total = 0
    for length, operation in parse_cigar(alignment.cigar):
        if operation in REFERENCE_CONSUMING:
            operation_end = ref + length
            if operation in ALIGNED:
                total += max(
                    0,
                    min(operation_end, window_end) - max(ref, window_start),
                )
            ref = operation_end
    return total


def source_to_target_bases(
    alignment: Alignment,
    source_start: int,
    source_end: int,
    target_start: int,
    target_end: int,
) -> tuple[int, int]:
    """Return source bases aligned anywhere and source bases aligned in target."""
    ref = alignment.ref_start
    query = alignment.query_length if alignment.strand == "-" else 0
    source_aligned = 0
    target_aligned = 0
    for length, operation in parse_cigar(alignment.cigar):
        consumes_query_coordinate = operation in QUERY_COORDINATE_CONSUMING
        if operation in ALIGNED:
            if alignment.strand == "+":
                query_low, query_high = query, query + length
            else:
                query_low, query_high = query - length, query
            overlap_low = max(source_start, query_low)
            overlap_high = min(source_end, query_high)
            if overlap_high > overlap_low:
                overlap = overlap_high - overlap_low
                source_aligned += overlap
                if alignment.strand == "+":
                    mapped_low = ref + overlap_low - query_low
                    mapped_high = ref + overlap_high - query_low
                else:
                    mapped_low = ref + query_high - overlap_high
                    mapped_high = ref + query_high - overlap_low
                target_aligned += max(
                    0,
                    min(mapped_high, target_end) - max(mapped_low, target_start),
                )
        if operation in REFERENCE_CONSUMING:
            ref += length
        if consumes_query_coordinate:
            query += length if alignment.strand == "+" else -length
    return source_aligned, target_aligned


def windows(target: tuple[str, int, int, str]) -> dict[str, tuple[int, int]]:
    _, start, end, interior = target
    result = {
        "element": (start, end),
        "left5": (max(0, start - 5_000), start),
        "right5": (end, end + 5_000),
    }
    if interior == "right":
        result["interior100"] = (end + 80_000, end + 100_000)
        result["interior1m"] = (end + 500_000, end + 1_000_000)
    else:
        result["interior100"] = (max(0, start - 100_000), max(0, start - 80_000))
        result["interior1m"] = (max(0, start - 1_000_000), max(0, start - 500_000))
    return result


def load_alignments(path: Path) -> dict[str, list[Alignment]]:
    by_qname: dict[str, list[Alignment]] = defaultdict(list)
    seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["build"], row["locus"])
            if key not in TARGETS:
                continue
            for match in DETAIL_RE.finditer(row["details"]):
                alignment = Alignment(
                    build=row["build"],
                    locus=row["locus"],
                    qname=match.group(1),
                    chrom=match.group(2),
                    ref_start=int(match.group(3)),
                    ref_end=int(match.group(4)),
                    query_start=int(match.group(5)),
                    query_end=int(match.group(6)),
                    query_length=int(match.group(7)),
                    strand=match.group(8),
                    mapq=int(match.group(9)),
                    primary=match.group(10) == "nonsecondary",
                    cigar=match.group(11),
                )
                identity = (
                    alignment.build,
                    alignment.locus,
                    alignment.qname,
                    alignment.chrom,
                    alignment.ref_start,
                    alignment.ref_end,
                    alignment.query_start,
                    alignment.query_end,
                    alignment.strand,
                    alignment.cigar,
                )
                if identity not in seen:
                    seen.add(identity)
                    by_qname[alignment.qname].append(alignment)
    return by_qname


def candidate_metrics(
    alignment: Alignment, source_start: int, source_end: int
) -> dict[str, object] | None:
    target = TARGETS[(alignment.build, alignment.locus)]
    target_windows = windows(target)
    source_length = source_end - source_start
    source_aligned, target_aligned = source_to_target_bases(
        alignment,
        source_start,
        source_end,
        target_windows["element"][0],
        target_windows["element"][1],
    )
    if source_length <= 0 or source_aligned / source_length < 0.8:
        return None
    if target_aligned / source_length < 0.5:
        return None
    fractions = {
        name: aligned_reference_bases(alignment, start, end) / (end - start)
        for name, (start, end) in target_windows.items()
    }
    left_flank = max(0, source_start - alignment.query_start)
    right_flank = max(0, alignment.query_end - source_end)
    return {
        "build": alignment.build,
        "locus": alignment.locus,
        "mapq": alignment.mapq,
        "primary": alignment.primary,
        "source_fraction": source_aligned / source_length,
        "target_source_fraction": target_aligned / source_length,
        "element_fraction": fractions["element"],
        "both5": fractions["left5"] >= 0.5 and fractions["right5"] >= 0.5,
        "interior100": fractions["interior100"] >= 0.5,
        "interior1m": fractions["interior1m"] >= 0.5,
        "left_flank": left_flank,
        "right_flank": right_flank,
        "min_flank": min(left_flank, right_flank),
        "max_flank": max(left_flank, right_flank),
    }


def best_candidates(
    qname: str,
    source_start: int,
    source_end: int,
    alignments: dict[str, list[Alignment]],
) -> list[dict[str, object]]:
    best: dict[tuple[str, str], dict[str, object]] = {}
    for alignment in alignments.get(qname, []):
        metrics = candidate_metrics(alignment, source_start, source_end)
        if metrics is None:
            continue
        key = (alignment.build, alignment.locus)
        score = (
            int(metrics["primary"]),
            int(metrics["both5"]),
            int(metrics["interior1m"]),
            int(metrics["interior100"]),
            int(metrics["mapq"]),
            int(metrics["max_flank"]),
        )
        previous = best.get(key)
        if previous is None or score > previous["_score"]:
            metrics["_score"] = score
            best[key] = metrics
    return sorted(best.values(), key=lambda item: item["_score"], reverse=True)


def format_candidate(candidate: dict[str, object]) -> str:
    return (
        f"{candidate['build']}:{candidate['locus']}"
        f":primary={int(candidate['primary'])}"
        f":mapq={candidate['mapq']}"
        f":source={candidate['source_fraction']:.3f}"
        f":target={candidate['target_source_fraction']:.3f}"
        f":element={candidate['element_fraction']:.3f}"
        f":both5={int(candidate['both5'])}"
        f":int100={int(candidate['interior100'])}"
        f":int1m={int(candidate['interior1m'])}"
        f":flanks={candidate['left_flank']},{candidate['right_flank']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orf", required=True)
    parser.add_argument("--attachments", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    alignments = load_alignments(Path(args.attachments))
    with Path(args.orf).open(newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["Structure"] in PRESENT and row["Locus"] in DUPLICATED_LOCI
        ]
    loci_by_source: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        loci_by_source[row["Source_Identifier"]].add(row["Locus"])

    output_rows = []
    candidates_by_source = {}
    for row in rows:
        source = source_interval(row["Source_Identifier"])
        if source is None:
            continue
        qname, start, end = source
        candidates = candidates_by_source.get(row["Source_Identifier"])
        if candidates is None:
            candidates = best_candidates(qname, start, end, alignments)
            candidates_by_source[row["Source_Identifier"]] = candidates
        if len(loci_by_source[row["Source_Identifier"]]) == 1 and not candidates:
            continue
        output_rows.append(
            {
                "ID_Full": row["ID_Full"],
                "raw_locus": row["Locus"],
                "Source_Identifier": row["Source_Identifier"],
                "same_source_loci": ",".join(
                    sorted(loci_by_source[row["Source_Identifier"]])
                ),
                "candidate_count": len(candidates),
                "best_build": candidates[0]["build"] if candidates else "",
                "best_locus": (
                    f"HML-2_{candidates[0]['locus']}" if candidates else ""
                ),
                "best_primary": int(candidates[0]["primary"]) if candidates else "",
                "best_mapq": candidates[0]["mapq"] if candidates else "",
                "best_both5": int(candidates[0]["both5"]) if candidates else "",
                "best_interior100": (
                    int(candidates[0]["interior100"]) if candidates else ""
                ),
                "best_interior1m": (
                    int(candidates[0]["interior1m"]) if candidates else ""
                ),
                "best_left_flank": candidates[0]["left_flank"] if candidates else "",
                "best_right_flank": candidates[0]["right_flank"] if candidates else "",
                "best_candidate": format_candidate(candidates[0]) if candidates else "",
                "all_candidates": ";".join(format_candidate(item) for item in candidates),
            }
        )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "ID_Full",
        "raw_locus",
        "Source_Identifier",
        "same_source_loci",
        "candidate_count",
        "best_build",
        "best_locus",
        "best_primary",
        "best_mapq",
        "best_both5",
        "best_interior100",
        "best_interior1m",
        "best_left_flank",
        "best_right_flank",
        "best_candidate",
        "all_candidates",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"wrote {len(output_rows)} rows to {output_path}")


if __name__ == "__main__":
    main()
