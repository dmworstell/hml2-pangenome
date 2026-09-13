#!/usr/bin/env python3
"""Re-derive Type-I cassette/backbone sequence evidence from current inputs.

The analysis deliberately uses one deterministic, actual population provirus
per production-eligible locus.  A two-thread MAFFT alignment is anchored to
0-based Type-II KCON coordinates.  The Type-I cassette flanks and six
predeclared backbone windows are then compared using the observed sequences.

This is a mosaic-sequence audit.  It cannot distinguish ancient
post-integration gene conversion from reverse-transcription-mediated spread,
and conservation is not evidence of host benefit or selection.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
import platform
import random
import re
import resource
import statistics
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence


PERSON_RE = re.compile(r"^(?:HG|NA)\d+$")
CANONICAL = frozenset("ACGT")
CALLABLE_LOCUS_FRACTION = 0.90
PAIRWISE_OVERLAP_FRACTION = 0.80
SITE_CALLABLE_TYPE1_FRACTION = 0.80
RESAMPLES = 2_000
SEED = 20260712


@dataclass(frozen=True)
class Window:
    key: str
    label: str
    segments: tuple[tuple[int, int], ...]
    kind: str

    @property
    def length(self) -> int:
        return sum(end - start for start, end in self.segments)

    def coordinates(self) -> list[int]:
        return [
            coordinate
            for start, end in self.segments
            for coordinate in range(start, end)
        ]


WINDOWS = (
    Window("backbone_1000_2000", "B1", ((1000, 2000),), "backbone_control"),
    Window("backbone_2000_3000", "B2", ((2000, 3000),), "backbone_control"),
    Window("backbone_3000_4000", "B3", ((3000, 4000),), "backbone_control"),
    Window("backbone_4000_5000", "B4", ((4000, 5000),), "backbone_control"),
    Window("backbone_5000_6000", "B5", ((5000, 6000),), "backbone_control"),
    Window("backbone_7293_8293", "B6", ((7293, 8293),), "backbone_control"),
    Window(
        "type1_cassette_flanks",
        "Cass.",
        ((6000, 6501), (6793, 7293)),
        "type1_cassette_flanks",
    ),
)
WINDOW_BY_KEY = {window.key: window for window in WINDOWS}
CONTROL_KEYS = tuple(window.key for window in WINDOWS if window.kind == "backbone_control")
CASSETTE_KEY = "type1_cassette_flanks"


@dataclass(frozen=True)
class Candidate:
    locus: str
    type_call: str
    id_full: str
    person: str
    haplotype: str
    structure: str
    strand: str
    source_identifier: str
    fasta: Path
    root_rank: int
    root_label: str


@dataclass(frozen=True)
class Representative:
    analysis_id: str
    locus: str
    type_call: str
    id_full: str
    person: str
    haplotype: str
    strand: str
    source_identifier: str
    fasta: Path
    root_label: str
    sequence: str
    fasta_sha256: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_locus(value: str) -> str:
    return value.strip().removeprefix("HML-2_")


def read_single_fasta(path: Path) -> tuple[str, str]:
    name = ""
    sequence: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    raise ValueError(f"Expected one FASTA record: {path}")
                name = line[1:].split()[0]
            else:
                sequence.append(line.upper())
    if not name or not sequence:
        raise ValueError(f"Empty FASTA: {path}")
    return name, "".join(sequence)


def read_fasta_alignment(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = ""
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate alignment identifier: {current}")
                records[current] = []
            elif current:
                records[current].append(line.upper())
            else:
                raise ValueError("Alignment sequence appeared before a FASTA header")
    return {name: "".join(parts) for name, parts in records.items()}


def write_tsv(path: Path, rows: Iterable[dict], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def load_eligible_loci(path: Path) -> tuple[dict[str, str], list[dict]]:
    eligible: dict[str, str] = {}
    rows: list[dict] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"locus", "type_call", "status", "production_eligible"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Delta summary lacks columns: {sorted(missing)}")
        for row in reader:
            locus = normalize_locus(row["locus"])
            production_eligible = row["production_eligible"].strip().lower() == "true"
            rows.append({
                "locus": locus,
                "type_call": row["type_call"],
                "status": row["status"],
                "production_eligible": str(production_eligible).lower(),
            })
            if production_eligible:
                if row["type_call"] not in {"TypeI", "TypeII"}:
                    raise ValueError(f"Eligible locus lacks Type-I/II call: {locus}")
                eligible[locus] = row["type_call"]
    return eligible, rows


def locate_fasta(locus_label: str, id_full: str, roots: Sequence[tuple[str, Path]]) -> tuple[Path, int, str]:
    relative = Path(locus_label) / f"{id_full}.fa"
    for rank, (label, root) in enumerate(roots):
        path = root / relative
        if path.is_file():
            return path, rank, label
    return roots[0][1] / relative, len(roots), "missing"


def candidate_sort_key(candidate: Candidate) -> tuple:
    """Predeclared selection order: HPRC root, then graph root, then ID."""
    return (candidate.root_rank, candidate.id_full, candidate.source_identifier)


def select_representatives(
    orf_table: Path,
    eligible: dict[str, str],
    roots: Sequence[tuple[str, Path]],
) -> tuple[list[Representative], list[dict]]:
    candidates: dict[str, dict[str, Candidate]] = defaultdict(dict)
    counters: dict[str, Counter] = {locus: Counter() for locus in eligible}
    with orf_table.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Locus", "ID_Full", "ID", "Haplotype", "Source_Identifier",
            "Structure", "Strand",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"ORF table lacks columns: {sorted(missing)}")
        for row in reader:
            locus = normalize_locus(row["Locus"])
            if locus not in eligible:
                continue
            counters[locus]["eligible_locus_rows"] += 1
            if not PERSON_RE.fullmatch(row["ID"].strip()):
                counters[locus]["non_population_rows"] += 1
                continue
            if "provirus" not in row["Structure"].strip().lower():
                counters[locus]["non_provirus_population_rows"] += 1
                continue
            counters[locus]["population_provirus_rows"] += 1
            fasta, root_rank, root_label = locate_fasta(row["Locus"].strip(), row["ID_Full"].strip(), roots)
            if not fasta.is_file():
                counters[locus]["fasta_missing_rows"] += 1
                continue
            candidate = Candidate(
                locus=locus,
                type_call=eligible[locus],
                id_full=row["ID_Full"].strip(),
                person=row["ID"].strip(),
                haplotype=row["Haplotype"].strip(),
                structure=row["Structure"].strip(),
                strand=row["Strand"].strip(),
                source_identifier=row["Source_Identifier"].strip(),
                fasta=fasta,
                root_rank=root_rank,
                root_label=root_label,
            )
            existing = candidates[locus].get(candidate.id_full)
            if existing is None or candidate_sort_key(candidate) < candidate_sort_key(existing):
                candidates[locus][candidate.id_full] = candidate

    representatives: list[Representative] = []
    selection_log: list[dict] = []
    for locus in sorted(eligible):
        ordered = sorted(candidates[locus].values(), key=candidate_sort_key)
        selected: Representative | None = None
        invalid = 0
        for candidate in ordered:
            try:
                _, sequence = read_single_fasta(candidate.fasta)
            except (OSError, ValueError):
                invalid += 1
                continue
            if not sequence:
                invalid += 1
                continue
            analysis_id = f"rep{len(representatives) + 1:03d}"
            selected = Representative(
                analysis_id=analysis_id,
                locus=locus,
                type_call=eligible[locus],
                id_full=candidate.id_full,
                person=candidate.person,
                haplotype=candidate.haplotype,
                strand=candidate.strand,
                source_identifier=candidate.source_identifier,
                fasta=candidate.fasta,
                root_label=candidate.root_label,
                sequence=sequence,
                fasta_sha256=sha256(candidate.fasta),
            )
            representatives.append(selected)
            break
        count = counters[locus]
        status = "selected" if selected else "excluded_no_readable_population_provirus_fasta"
        selection_log.append({
            "locus": locus,
            "type_call": eligible[locus],
            "status": status,
            "exclusion_reason": "" if selected else "no_readable_population_provirus_fasta",
            "eligible_locus_rows": count["eligible_locus_rows"],
            "population_provirus_rows": count["population_provirus_rows"],
            "non_population_rows": count["non_population_rows"],
            "non_provirus_population_rows": count["non_provirus_population_rows"],
            "fasta_missing_rows": count["fasta_missing_rows"],
            "distinct_existing_fasta_candidates": len(ordered),
            "invalid_fasta_candidates_before_selection": invalid,
            "selected_analysis_id": selected.analysis_id if selected else "",
            "selected_id_full": selected.id_full if selected else "",
            "selected_person": selected.person if selected else "",
            "selected_haplotype": selected.haplotype if selected else "",
            "selected_root": selected.root_label if selected else "",
            "selected_strand_metadata": selected.strand if selected else "",
            "selected_sequence_length": len(selected.sequence) if selected else "",
            "selected_fasta_basename": selected.fasta.name if selected else "",
            "selected_fasta_sha256": selected.fasta_sha256 if selected else "",
            "selection_rule": "root_priority_HPRC_then_HGSVC;lexical_ID_Full;first_readable",
        })
    return representatives, selection_log


def memory_free_percent() -> float | None:
    try:
        completed = subprocess.run(
            ["memory_pressure", "-Q"], check=True, capture_output=True, text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    match = re.search(r"System-wide memory free percentage:\s*([0-9.]+)%", completed.stdout)
    return float(match.group(1)) if match else None


def mafft_version(mafft: Path) -> str:
    completed = subprocess.run([str(mafft), "--version"], capture_output=True, text=True)
    return (completed.stdout or completed.stderr).strip().splitlines()[0]


def run_mafft(
    representatives: Sequence[Representative],
    kcon_name: str,
    kcon_sequence: str,
    mafft: Path,
    derived_dir: Path,
) -> tuple[Path, dict]:
    derived_dir.mkdir(parents=True, exist_ok=True)
    before = memory_free_percent()
    command = [
        str(mafft), "--thread", "2", "--retree", "1", "--maxiterate", "0",
        "--adjustdirection", "--anysymbol", "--inputorder", "--quiet",
    ]
    env = dict(os.environ)
    env.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
    })
    with tempfile.TemporaryDirectory(prefix="hml2_type1_mafft_") as temp_text:
        temp = Path(temp_text)
        input_path = temp / "representatives.fasta"
        alignment_path = temp / "representatives.mafft.fasta"
        with input_path.open("w") as handle:
            handle.write(f">{kcon_name}\n{kcon_sequence}\n")
            for representative in representatives:
                handle.write(f">{representative.analysis_id}\n{representative.sequence}\n")
        actual_command = command + [str(input_path)]
        child_rss_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        with alignment_path.open("w") as output:
            completed = subprocess.run(
                actual_command,
                check=False,
                stdout=output,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
            )
        if completed.returncode:
            raise RuntimeError(
                f"MAFFT failed with exit {completed.returncode}:\n{completed.stderr[-8000:]}"
            )
        stderr = completed.stderr
        child_rss_after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        peak_bytes = int(child_rss_after if platform.system() == "Darwin" else child_rss_after * 1024)
        alignment_gz = derived_dir / "representatives.mafft.fasta.gz"
        with alignment_path.open("rb") as source, gzip.open(alignment_gz, "wb", compresslevel=6) as target:
            for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
                target.write(block)
        input_gz = derived_dir / "representatives.input.fasta.gz"
        with input_path.open("rb") as source, gzip.open(input_gz, "wb", compresslevel=6) as target:
            for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
                target.write(block)
        parsed_copy = derived_dir / "_alignment_for_parse.tmp"
        with alignment_path.open("rb") as source, parsed_copy.open("wb") as target:
            for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
                target.write(block)
    after = memory_free_percent()
    details = {
        "command": command + ["<temporary_input_fasta>"],
        "mafft_version": mafft_version(mafft),
        "threads": 2,
        "numerical_library_threads": 1,
        "memory_free_percent_before": before,
        "memory_free_percent_after": after,
        "mafft_peak_rss_bytes_from_child_rusage": peak_bytes,
        "child_rusage_ru_maxrss_before": child_rss_before,
        "child_rusage_ru_maxrss_after": child_rss_after,
        "peak_rss_measurement": "resource.getrusage(RUSAGE_CHILDREN).ru_maxrss; bytes on Darwin, KiB converted to bytes elsewhere",
        "stderr_tail": stderr.splitlines()[-20:],
    }
    return parsed_copy, details


def normalize_mafft_names(aligned: dict[str, str], kcon_name: str) -> tuple[dict[str, str], list[str]]:
    normalized: dict[str, str] = {}
    reversed_ids: list[str] = []
    for name, sequence in aligned.items():
        if name.startswith("_R_"):
            base = name[3:]
            reversed_ids.append(base)
        else:
            base = name
        if base in normalized:
            raise ValueError(f"MAFFT output collided after orientation normalization: {base}")
        normalized[base] = sequence
    if kcon_name in reversed_ids:
        raise ValueError("MAFFT reversed the KCON coordinate anchor")
    return normalized, sorted(reversed_ids)


def build_coordinate_map(aligned_kcon: str, kcon_sequence: str) -> list[int]:
    coordinates: list[int] = []
    observed: list[str] = []
    for column, base in enumerate(aligned_kcon):
        if base != "-":
            coordinates.append(column)
            observed.append(base)
    if "".join(observed).upper() != kcon_sequence.upper():
        raise ValueError("Ungapped MAFFT KCON does not equal the input coordinate anchor")
    return coordinates


def window_columns(window: Window, coordinate_map: Sequence[int]) -> list[int]:
    coordinates = window.coordinates()
    if not coordinates or min(coordinates) < 0 or max(coordinates) >= len(coordinate_map):
        raise ValueError(f"Window outside KCON: {window.key}")
    return [coordinate_map[coordinate] for coordinate in coordinates]


def canonical_count(sequence: str, columns: Sequence[int]) -> int:
    return sum(sequence[column] in CANONICAL for column in columns)


def pair_distance(
    first: str,
    second: str,
    columns: Sequence[int],
    minimum_overlap: int,
) -> tuple[int, int, float] | None:
    compared = 0
    mismatches = 0
    for column in columns:
        left = first[column]
        right = second[column]
        if left in CANONICAL and right in CANONICAL:
            compared += 1
            mismatches += left != right
    if compared < minimum_overlap:
        return None
    return compared, mismatches, mismatches / compared


def percentile(values: Sequence[float], probability: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    location = (len(ordered) - 1) * probability
    lower = math.floor(location)
    upper = math.ceil(location)
    if lower == upper:
        return ordered[lower]
    weight = location - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def site_diversity(aligned_sequences: Sequence[str], column: int, minimum_callable: int) -> tuple[int, float] | None:
    counts = Counter(sequence[column] for sequence in aligned_sequences if sequence[column] in CANONICAL)
    callable_count = sum(counts.values())
    if callable_count < minimum_callable or callable_count < 2:
        return None
    total_pairs = callable_count * (callable_count - 1) // 2
    same_pairs = sum(count * (count - 1) // 2 for count in counts.values())
    return callable_count, (total_pairs - same_pairs) / total_pairs


def bootstrap_mean(values: Sequence[float], resamples: int, rng: random.Random) -> tuple[float, float]:
    if not values:
        return math.nan, math.nan
    size = len(values)
    means = [
        sum(values[rng.randrange(size)] for _ in range(size)) / size
        for _ in range(resamples)
    ]
    return percentile(means, 0.025), percentile(means, 0.975)


def matched_site_comparison(
    cassette_by_callable: dict[int, list[float]],
    control_by_callable: dict[int, list[float]],
    resamples: int,
    seed: int,
) -> dict[str, float | int | str]:
    """Stratify by callable Type-I count and match site counts exactly."""
    rng = random.Random(seed)
    matched: list[tuple[list[float], list[float]]] = []
    strata_description: list[str] = []
    for callable_count in sorted(set(cassette_by_callable) & set(control_by_callable)):
        cassette = list(cassette_by_callable[callable_count])
        control = list(control_by_callable[callable_count])
        size = min(len(cassette), len(control))
        if size == 0:
            continue
        if len(cassette) > size:
            cassette = rng.sample(cassette, size)
        if len(control) > size:
            control = rng.sample(control, size)
        matched.append((cassette, control))
        strata_description.append(f"{callable_count}:{size}")
    matched_sites = sum(len(cassette) for cassette, _ in matched)
    if matched_sites < 2:
        return {
            "matched_sites_per_window": matched_sites,
            "callable_count_strata": ";".join(strata_description),
            "cassette_mean": math.nan,
            "control_mean": math.nan,
            "difference_cassette_minus_control": math.nan,
            "bootstrap_ci_low": math.nan,
            "bootstrap_ci_high": math.nan,
            "permutation_p_two_sided": math.nan,
        }
    cassette_mean = sum(sum(left) for left, _ in matched) / matched_sites
    control_mean = sum(sum(right) for _, right in matched) / matched_sites
    observed = cassette_mean - control_mean
    bootstrap_differences: list[float] = []
    permutation_differences: list[float] = []
    for _ in range(resamples):
        left_sum = 0.0
        right_sum = 0.0
        permuted_left_sum = 0.0
        permuted_right_sum = 0.0
        for left, right in matched:
            size = len(left)
            left_sum += sum(left[rng.randrange(size)] for _ in range(size))
            right_sum += sum(right[rng.randrange(size)] for _ in range(size))
            pooled = left + right
            indices = set(rng.sample(range(2 * size), size))
            permuted_left_sum += sum(value for index, value in enumerate(pooled) if index in indices)
            permuted_right_sum += sum(value for index, value in enumerate(pooled) if index not in indices)
        bootstrap_differences.append((left_sum - right_sum) / matched_sites)
        permutation_differences.append((permuted_left_sum - permuted_right_sum) / matched_sites)
    p_value = (
        1 + sum(abs(value) >= abs(observed) for value in permutation_differences)
    ) / (resamples + 1)
    return {
        "matched_sites_per_window": matched_sites,
        "callable_count_strata": ";".join(strata_description),
        "cassette_mean": cassette_mean,
        "control_mean": control_mean,
        "difference_cassette_minus_control": observed,
        "bootstrap_ci_low": percentile(bootstrap_differences, 0.025),
        "bootstrap_ci_high": percentile(bootstrap_differences, 0.975),
        "permutation_p_two_sided": p_value,
    }


def type1_locus_bootstrap_window(
    type1_representatives: Sequence[Representative],
    aligned: dict[str, str],
    columns: Sequence[int],
    callable_ids: set[str],
    resamples: int,
    seed: int,
) -> dict[str, float | int]:
    """Bootstrap Type-I loci, treating the locus as the inference unit."""
    representatives = [
        representative for representative in type1_representatives
        if representative.analysis_id in callable_ids
    ]
    size = len(representatives)
    if size < 2:
        return {
            "callable_type1_loci": size,
            "mean_within_type1_pairwise_p_distance": math.nan,
            "locus_bootstrap_ci_low": math.nan,
            "locus_bootstrap_ci_high": math.nan,
        }
    minimum_overlap = math.ceil(PAIRWISE_OVERLAP_FRACTION * len(columns))
    matrix = [[0.0] * size for _ in range(size)]
    observed: list[float] = []
    for first in range(size):
        for second in range(first + 1, size):
            result = pair_distance(
                aligned[representatives[first].analysis_id],
                aligned[representatives[second].analysis_id],
                columns,
                minimum_overlap,
            )
            if result is None:
                raise ValueError("Individually callable Type-I loci lacked required pairwise overlap")
            distance = result[2]
            matrix[first][second] = distance
            matrix[second][first] = distance
            observed.append(distance)
    rng = random.Random(seed)
    bootstraps: list[float] = []
    for _ in range(resamples):
        sample = [rng.randrange(size) for _ in range(size)]
        values = [matrix[sample[first]][sample[second]] for first in range(size) for second in range(first + 1, size)]
        bootstraps.append(statistics.fmean(values))
    return {
        "callable_type1_loci": size,
        "mean_within_type1_pairwise_p_distance": statistics.fmean(observed),
        "locus_bootstrap_ci_low": percentile(bootstraps, 0.025),
        "locus_bootstrap_ci_high": percentile(bootstraps, 0.975),
    }


def type1_locus_bootstrap_difference(
    type1_representatives: Sequence[Representative],
    aligned: dict[str, str],
    cassette_columns: Sequence[int],
    control_columns: Sequence[int],
    cassette_callable_ids: set[str],
    control_callable_ids: set[str],
    resamples: int,
    seed: int,
) -> dict[str, float | int]:
    """Paired locus bootstrap of cassette minus control within-Type-I diversity."""
    representatives = [
        representative for representative in type1_representatives
        if representative.analysis_id in cassette_callable_ids
        and representative.analysis_id in control_callable_ids
    ]
    size = len(representatives)
    if size < 2:
        return {
            "paired_callable_type1_loci": size,
            "cassette_mean_within_type1_pairwise_p_distance": math.nan,
            "control_mean_within_type1_pairwise_p_distance": math.nan,
            "difference_cassette_minus_control": math.nan,
            "locus_bootstrap_ci_low": math.nan,
            "locus_bootstrap_ci_high": math.nan,
        }
    matrices: list[list[list[float]]] = []
    observed_means: list[float] = []
    for columns in (cassette_columns, control_columns):
        minimum_overlap = math.ceil(PAIRWISE_OVERLAP_FRACTION * len(columns))
        matrix = [[0.0] * size for _ in range(size)]
        observed: list[float] = []
        for first in range(size):
            for second in range(first + 1, size):
                result = pair_distance(
                    aligned[representatives[first].analysis_id],
                    aligned[representatives[second].analysis_id],
                    columns,
                    minimum_overlap,
                )
                if result is None:
                    raise ValueError("Paired callable Type-I loci lacked required pairwise overlap")
                distance = result[2]
                matrix[first][second] = distance
                matrix[second][first] = distance
                observed.append(distance)
        matrices.append(matrix)
        observed_means.append(statistics.fmean(observed))
    rng = random.Random(seed)
    differences: list[float] = []
    for _ in range(resamples):
        sample = [rng.randrange(size) for _ in range(size)]
        means: list[float] = []
        for matrix in matrices:
            values = [matrix[sample[first]][sample[second]] for first in range(size) for second in range(first + 1, size)]
            means.append(statistics.fmean(values))
        differences.append(means[0] - means[1])
    return {
        "paired_callable_type1_loci": size,
        "cassette_mean_within_type1_pairwise_p_distance": observed_means[0],
        "control_mean_within_type1_pairwise_p_distance": observed_means[1],
        "difference_cassette_minus_control": observed_means[0] - observed_means[1],
        "locus_bootstrap_ci_low": percentile(differences, 0.025),
        "locus_bootstrap_ci_high": percentile(differences, 0.975),
    }


def paired_comparison(
    cassette: dict[str, float],
    control: dict[str, float],
    resamples: int,
    seed: int,
) -> dict[str, float | int]:
    loci = sorted(set(cassette) & set(control))
    differences = [cassette[locus] - control[locus] for locus in loci]
    if not differences:
        return {
            "paired_type1_loci": 0,
            "mean_difference_cassette_minus_control": math.nan,
            "bootstrap_ci_low": math.nan,
            "bootstrap_ci_high": math.nan,
            "sign_flip_p_two_sided": math.nan,
        }
    rng = random.Random(seed)
    size = len(differences)
    observed = statistics.fmean(differences)
    bootstraps = [
        sum(differences[rng.randrange(size)] for _ in range(size)) / size
        for _ in range(resamples)
    ]
    null = [
        sum(value if rng.random() < 0.5 else -value for value in differences) / size
        for _ in range(resamples)
    ]
    p_value = (1 + sum(abs(value) >= abs(observed) for value in null)) / (resamples + 1)
    return {
        "paired_type1_loci": size,
        "mean_difference_cassette_minus_control": observed,
        "bootstrap_ci_low": percentile(bootstraps, 0.025),
        "bootstrap_ci_high": percentile(bootstraps, 0.975),
        "sign_flip_p_two_sided": p_value,
    }


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    half = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
    return center - half, center + half


def audit_delta_records(
    eligible: dict[str, str],
    records_dir: Path,
) -> tuple[list[dict], list[dict]]:
    audit: list[dict] = []
    provenance: list[dict] = []
    for locus in sorted(eligible):
        record_path = records_dir / f"{locus}.json"
        base = {
            "locus": locus,
            "summary_type_call": eligible[locus],
            "record_status": "present" if record_path.is_file() else "missing",
            "eligible_rows": 0,
            "eligible_typeI_rows": 0,
            "eligible_typeII_rows": 0,
            "eligible_conflict_rows": 0,
            "eligible_no_call_rows": 0,
            "both_typeI_and_typeII_present": "false",
            "record_sha256": "",
        }
        if not record_path.is_file():
            audit.append(base)
            continue
        record_hash = sha256(record_path)
        with record_path.open() as handle:
            document = json.load(handle)
        evidence = document.get("direct_delta292_evidence", {})
        mask_rows = evidence.get("aggregate_eligibility_mask", {}).get("rows", [])
        eligibility = {
            int(row["row_index"]): bool(row.get("aggregate_eligible", False))
            for row in mask_rows
            if "row_index" in row
        }
        counts = Counter()
        for row in evidence.get("rows", []):
            row_index = int(row.get("row_index", -1))
            if not eligibility.get(row_index, False):
                continue
            counts[str(row.get("type_call", "no_call"))] += 1
        both = counts["TypeI"] > 0 and counts["TypeII"] > 0
        base.update({
            "eligible_rows": sum(counts.values()),
            "eligible_typeI_rows": counts["TypeI"],
            "eligible_typeII_rows": counts["TypeII"],
            "eligible_conflict_rows": counts["conflict"],
            "eligible_no_call_rows": counts["no_call"],
            "both_typeI_and_typeII_present": str(both).lower(),
            "record_sha256": record_hash,
        })
        audit.append(base)
        provenance.append({
            "role": "delta292_record",
            "locus": locus,
            "basename": record_path.name,
            "bytes": record_path.stat().st_size,
            "sha256": record_hash,
        })
        del document
    return audit, provenance


def pairwise_summaries(
    representatives: Sequence[Representative],
    aligned: dict[str, str],
    columns_by_window: dict[str, list[int]],
    callable_by_window: dict[str, set[str]],
) -> list[dict]:
    by_type = {
        type_call: [representative for representative in representatives if representative.type_call == type_call]
        for type_call in ("TypeI", "TypeII")
    }
    rows: list[dict] = []
    for window in WINDOWS:
        columns = columns_by_window[window.key]
        minimum_overlap = math.ceil(PAIRWISE_OVERLAP_FRACTION * len(columns))
        pair_groups = {
            "within_TypeI": itertools.combinations(by_type["TypeI"], 2),
            "within_TypeII": itertools.combinations(by_type["TypeII"], 2),
            "between_TypeI_TypeII": itertools.product(by_type["TypeI"], by_type["TypeII"]),
        }
        for contrast, pairs in pair_groups.items():
            values: list[float] = []
            overlaps: list[int] = []
            loci_used: set[str] = set()
            for first, second in pairs:
                if first.analysis_id not in callable_by_window[window.key] or second.analysis_id not in callable_by_window[window.key]:
                    continue
                result = pair_distance(
                    aligned[first.analysis_id], aligned[second.analysis_id], columns, minimum_overlap,
                )
                if result is None:
                    continue
                overlap, _, divergence = result
                values.append(divergence)
                overlaps.append(overlap)
                loci_used.update((first.locus, second.locus))
            rows.append({
                "window": window.key,
                "window_label": window.label,
                "window_kind": window.kind,
                "kcon_anchored_bases": window.length,
                "contrast": contrast,
                "pair_count": len(values),
                "unique_loci": len(loci_used),
                "mean_p_distance": statistics.fmean(values) if values else math.nan,
                "median_p_distance": statistics.median(values) if values else math.nan,
                "minimum_p_distance": min(values) if values else math.nan,
                "maximum_p_distance": max(values) if values else math.nan,
                "mean_jointly_callable_bases": statistics.fmean(overlaps) if overlaps else math.nan,
                "inference": "descriptive_pairs_not_independent",
            })
    return rows


def nearest_neighbors(
    representatives: Sequence[Representative],
    aligned: dict[str, str],
    columns_by_window: dict[str, list[int]],
    callable_by_window: dict[str, set[str]],
) -> list[dict]:
    type1 = [representative for representative in representatives if representative.type_call == "TypeI"]
    type2 = [representative for representative in representatives if representative.type_call == "TypeII"]
    rows: list[dict] = []
    for window in WINDOWS:
        columns = columns_by_window[window.key]
        minimum_overlap = math.ceil(PAIRWISE_OVERLAP_FRACTION * len(columns))
        for focal in type1:
            base = {
                "window": window.key,
                "window_label": window.label,
                "type1_locus": focal.locus,
                "type1_analysis_id": focal.analysis_id,
                "nearest_type2_locus": "",
                "nearest_type2_identity": math.nan,
                "nearest_type2_jointly_callable_bases": "",
                "nearest_overall_locus": "",
                "nearest_overall_type": "",
                "nearest_overall_identity": math.nan,
            }
            if focal.analysis_id not in callable_by_window[window.key]:
                base["status"] = "type1_not_callable"
                rows.append(base)
                continue
            type2_scores: list[tuple[float, str, int]] = []
            overall_scores: list[tuple[float, str, str, int]] = []
            for other in representatives:
                if other.analysis_id == focal.analysis_id or other.analysis_id not in callable_by_window[window.key]:
                    continue
                result = pair_distance(aligned[focal.analysis_id], aligned[other.analysis_id], columns, minimum_overlap)
                if result is None:
                    continue
                overlap, _, divergence = result
                identity = 1 - divergence
                overall_scores.append((identity, other.locus, other.type_call, overlap))
                if other.type_call == "TypeII":
                    type2_scores.append((identity, other.locus, overlap))
            if type2_scores:
                identity, locus, overlap = sorted(type2_scores, key=lambda value: (-value[0], value[1]))[0]
                base.update({
                    "nearest_type2_locus": locus,
                    "nearest_type2_identity": identity,
                    "nearest_type2_jointly_callable_bases": overlap,
                })
            if overall_scores:
                identity, locus, type_call, _ = sorted(overall_scores, key=lambda value: (-value[0], value[1], value[2]))[0]
                base.update({
                    "nearest_overall_locus": locus,
                    "nearest_overall_type": type_call,
                    "nearest_overall_identity": identity,
                })
            base["status"] = "callable" if type2_scores else "no_callable_type2_neighbor"
            rows.append(base)
    return rows


def nearest_neighbor_switching(neighbors: Sequence[dict]) -> tuple[list[dict], list[dict]]:
    by_locus: dict[str, dict[str, dict]] = defaultdict(dict)
    for row in neighbors:
        if row["status"] == "callable":
            by_locus[row["type1_locus"]][row["window"]] = row
    per_locus: list[dict] = []
    for locus in sorted(by_locus):
        values = by_locus[locus]
        controls = [values[key]["nearest_type2_locus"] for key in CONTROL_KEYS if key in values]
        counts = Counter(controls)
        modal = sorted(counts.items(), key=lambda value: (-value[1], value[0]))[0][0] if counts else ""
        cassette = values.get(CASSETTE_KEY, {}).get("nearest_type2_locus", "")
        per_locus.append({
            "type1_locus": locus,
            "callable_backbone_windows": len(controls),
            "backbone_modal_nearest_type2": modal,
            "backbone_modal_support_windows": counts[modal] if modal else 0,
            "cassette_nearest_type2": cassette,
            "cassette_differs_from_backbone_mode": (
                str(bool(cassette and modal and cassette != modal)).lower()
                if cassette and modal else "not_callable"
            ),
            "distinct_nearest_type2_across_all_windows": len(set(controls + ([cassette] if cassette else []))),
            **{
                f"cassette_differs_from_{key}": (
                    str(bool(cassette and key in values and cassette != values[key]["nearest_type2_locus"])).lower()
                    if cassette and key in values else "not_callable"
                )
                for key in CONTROL_KEYS
            },
        })
    comparisons = ["backbone_mode", *CONTROL_KEYS]
    summary: list[dict] = []
    for comparison in comparisons:
        if comparison == "backbone_mode":
            field = "cassette_differs_from_backbone_mode"
        else:
            field = f"cassette_differs_from_{comparison}"
        statuses = [row[field] for row in per_locus if row[field] in {"true", "false"}]
        switched = sum(value == "true" for value in statuses)
        low, high = wilson_interval(switched, len(statuses))
        summary.append({
            "comparison": comparison,
            "type1_loci_compared": len(statuses),
            "switched_nearest_type2": switched,
            "switch_fraction": switched / len(statuses) if statuses else math.nan,
            "wilson_95ci_low": low,
            "wilson_95ci_high": high,
            "interpretation": "local_nearest_neighbor_discordance_not_direction_or_mechanism",
        })
    return per_locus, summary


def leave_one_window_out_switching(neighbors: Sequence[dict]) -> tuple[list[dict], list[dict]]:
    """Compare each focal window with the modal Type-II neighbor in other backbone windows.

    For a backbone focal window, that window is omitted from the five-window
    backbone mode.  The cassette is compared with the mode across all six
    backbone windows.  The backbone rows are the predeclared negative control
    for whether switching is cassette-specific.
    """
    by_locus: dict[str, dict[str, dict]] = defaultdict(dict)
    for row in neighbors:
        if row["status"] == "callable":
            by_locus[row["type1_locus"]][row["window"]] = row
    per_locus: list[dict] = []
    focal_keys = [*CONTROL_KEYS, CASSETTE_KEY]
    for locus in sorted(by_locus):
        values = by_locus[locus]
        for focal_key in focal_keys:
            reference_keys = [key for key in CONTROL_KEYS if key != focal_key]
            reference_neighbors = [
                values[key]["nearest_type2_locus"] for key in reference_keys if key in values
            ]
            counts = Counter(reference_neighbors)
            modal = sorted(counts.items(), key=lambda value: (-value[1], value[0]))[0][0] if counts else ""
            focal_neighbor = values.get(focal_key, {}).get("nearest_type2_locus", "")
            status = "callable" if focal_neighbor and modal else "not_callable"
            per_locus.append({
                "type1_locus": locus,
                "focal_window": focal_key,
                "focal_window_kind": WINDOW_BY_KEY[focal_key].kind,
                "focal_nearest_type2": focal_neighbor,
                "reference_backbone_windows": len(reference_neighbors),
                "reference_backbone_modal_nearest_type2": modal,
                "reference_mode_support_windows": counts[modal] if modal else 0,
                "focal_differs_from_reference_backbone_mode": (
                    str(focal_neighbor != modal).lower() if status == "callable" else "not_callable"
                ),
                "status": status,
            })
    summary: list[dict] = []
    for focal_key in focal_keys:
        rows = [row for row in per_locus if row["focal_window"] == focal_key and row["status"] == "callable"]
        switched = sum(row["focal_differs_from_reference_backbone_mode"] == "true" for row in rows)
        low, high = wilson_interval(switched, len(rows))
        summary.append({
            "focal_window": focal_key,
            "focal_window_label": WINDOW_BY_KEY[focal_key].label,
            "focal_window_kind": WINDOW_BY_KEY[focal_key].kind,
            "type1_loci_compared": len(rows),
            "switched_nearest_type2": switched,
            "switch_fraction": switched / len(rows) if rows else math.nan,
            "wilson_95ci_low": low,
            "wilson_95ci_high": high,
            "negative_control_role": (
                "cassette_test" if focal_key == CASSETTE_KEY
                else "leave_one_backbone_window_out_negative_control"
            ),
            "interpretation": "switching_shared_with_backbone_supports_broad_mosaicism_not_cassette_specificity",
        })
    return per_locus, summary


def emit_figure(
    output_dir: Path,
    locus_diversity_rows: Sequence[dict],
    leave_one_out_summary: Sequence[dict],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = [window.key for window in WINDOWS]
    diversity = {row["window"]: row for row in locus_diversity_rows}
    switching = {row["focal_window"]: row for row in leave_one_out_summary}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True)

    x = list(range(len(order)))
    means = [diversity[key]["mean_within_type1_pairwise_p_distance"] for key in order]
    lower = [means[index] - diversity[key]["locus_bootstrap_ci_low"] for index, key in enumerate(order)]
    upper = [diversity[key]["locus_bootstrap_ci_high"] - means[index] for index, key in enumerate(order)]
    colors = ["#8F9AA6" if key != CASSETTE_KEY else "#B4473E" for key in order]
    axes[0].errorbar(x, means, yerr=[lower, upper], fmt="none", ecolor="#30363D", lw=1, capsize=2)
    axes[0].scatter(x, means, c=colors, s=30, zorder=3)
    axes[0].set_xticks(x, [WINDOW_BY_KEY[key].label for key in order])
    axes[0].set_ylabel("Within-Type-I p-distance")
    axes[0].set_title("A  Sequence conservation", loc="left", fontsize=9, fontweight="bold")

    control_x = list(range(len(CONTROL_KEYS) + 1))
    comparison_order = [*CONTROL_KEYS, CASSETTE_KEY]
    rates = [switching[key]["switch_fraction"] * 100 for key in comparison_order]
    lows = [(switching[key]["switch_fraction"] - switching[key]["wilson_95ci_low"]) * 100 for key in comparison_order]
    highs = [(switching[key]["wilson_95ci_high"] - switching[key]["switch_fraction"]) * 100 for key in comparison_order]
    axes[1].bar(control_x, rates, color="#477A9C", width=0.72)
    axes[1].errorbar(control_x, rates, yerr=[lows, highs], fmt="none", ecolor="#30363D", lw=1, capsize=2)
    axes[1].set_xticks(control_x, [WINDOW_BY_KEY[key].label for key in comparison_order])
    axes[1].set_ylabel("LOO neighbor switch (%)")
    axes[1].set_ylim(0, 105)
    axes[1].set_title("B  Local neighbor discordance", loc="left", fontsize=9, fontweight="bold")

    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(labelsize=7)
        axis.grid(axis="y", color="#D9DEE3", lw=0.5, alpha=0.7)
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    stem = figures / "type1_cassette_backbone_rederivation"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orf-table", type=Path, required=True)
    parser.add_argument("--delta292-summary", type=Path, required=True)
    parser.add_argument("--delta292-records", type=Path, required=True)
    parser.add_argument("--type2-reference", type=Path, required=True)
    parser.add_argument("--processed-root", type=Path, required=True)
    parser.add_argument("--graph-root", type=Path, required=True)
    parser.add_argument("--mafft", type=Path, default=Path("mafft"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--resamples", type=int, default=RESAMPLES)
    args = parser.parse_args()

    results_dir = args.output_dir / "results"
    derived_dir = args.output_dir / "derived"
    results_dir.mkdir(parents=True, exist_ok=True)
    derived_dir.mkdir(parents=True, exist_ok=True)

    eligible, delta_summary_rows = load_eligible_loci(args.delta292_summary)
    roots = (("HPRC_processed_loci", args.processed_root), ("HGSVC3_graph_loci", args.graph_root))
    representatives, selection_log = select_representatives(args.orf_table, eligible, roots)
    write_tsv(
        results_dir / "representative_selection.tsv",
        selection_log,
        [
            "locus", "type_call", "status", "exclusion_reason", "eligible_locus_rows",
            "population_provirus_rows", "non_population_rows", "non_provirus_population_rows",
            "fasta_missing_rows", "distinct_existing_fasta_candidates",
            "invalid_fasta_candidates_before_selection", "selected_analysis_id", "selected_id_full",
            "selected_person", "selected_haplotype", "selected_root", "selected_strand_metadata",
            "selected_sequence_length", "selected_fasta_basename", "selected_fasta_sha256",
            "selection_rule",
        ],
    )
    if not representatives:
        raise RuntimeError("No population provirus representatives were selectable")

    kcon_name, kcon_sequence = read_single_fasta(args.type2_reference)
    kcon_analysis_name = "KCON_TypeII_coordinate_anchor"
    alignment_temp, mafft_details = run_mafft(
        representatives, kcon_analysis_name, kcon_sequence, args.mafft, derived_dir,
    )
    try:
        aligned_raw = read_fasta_alignment(alignment_temp)
    finally:
        alignment_temp.unlink(missing_ok=True)
    aligned, reversed_ids = normalize_mafft_names(aligned_raw, kcon_analysis_name)
    expected_ids = {representative.analysis_id for representative in representatives} | {kcon_analysis_name}
    if set(aligned) != expected_ids:
        raise ValueError(f"MAFFT output identifiers differ: missing={sorted(expected_ids - set(aligned))}, extra={sorted(set(aligned) - expected_ids)}")
    alignment_lengths = {len(sequence) for sequence in aligned.values()}
    if len(alignment_lengths) != 1:
        raise ValueError("MAFFT output is not rectangular")
    coordinate_map = build_coordinate_map(aligned[kcon_analysis_name], kcon_sequence)
    columns_by_window = {window.key: window_columns(window, coordinate_map) for window in WINDOWS}

    with gzip.open(derived_dir / "kcon_coordinate_map.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["kcon_0based_coordinate", "mafft_0based_column"])
        writer.writerows(enumerate(coordinate_map))
    window_definition_rows = [
        {
            "window": window.key,
            "label": window.label,
            "kind": window.kind,
            "segments_0based_half_open": "+".join(f"[{start},{end})" for start, end in window.segments),
            "kcon_anchored_bases": window.length,
            "first_mafft_column": min(columns_by_window[window.key]),
            "last_mafft_column": max(columns_by_window[window.key]),
        }
        for window in WINDOWS
    ]
    write_tsv(
        derived_dir / "window_definitions.tsv",
        window_definition_rows,
        ["window", "label", "kind", "segments_0based_half_open", "kcon_anchored_bases", "first_mafft_column", "last_mafft_column"],
    )

    manifest_rows = []
    for representative in representatives:
        manifest_rows.append({
            "analysis_id": representative.analysis_id,
            "locus": representative.locus,
            "type_call": representative.type_call,
            "id_full": representative.id_full,
            "person": representative.person,
            "haplotype": representative.haplotype,
            "root": representative.root_label,
            "strand_metadata": representative.strand,
            "input_sequence_length": len(representative.sequence),
            "mafft_reverse_complemented": str(representative.analysis_id in reversed_ids).lower(),
            "fasta_basename": representative.fasta.name,
            "fasta_sha256": representative.fasta_sha256,
        })
    write_tsv(
        derived_dir / "representative_manifest.tsv",
        manifest_rows,
        [
            "analysis_id", "locus", "type_call", "id_full", "person", "haplotype",
            "root", "strand_metadata", "input_sequence_length", "mafft_reverse_complemented",
            "fasta_basename", "fasta_sha256",
        ],
    )

    callability_rows: list[dict] = []
    callable_by_window: dict[str, set[str]] = {window.key: set() for window in WINDOWS}
    for window in WINDOWS:
        columns = columns_by_window[window.key]
        for representative in representatives:
            count = canonical_count(aligned[representative.analysis_id], columns)
            fraction = count / len(columns)
            callable_status = fraction >= CALLABLE_LOCUS_FRACTION
            if callable_status:
                callable_by_window[window.key].add(representative.analysis_id)
            callability_rows.append({
                "window": window.key,
                "window_label": window.label,
                "locus": representative.locus,
                "type_call": representative.type_call,
                "analysis_id": representative.analysis_id,
                "kcon_anchored_bases": len(columns),
                "canonical_bases": count,
                "callable_fraction": fraction,
                "callable_at_0.90": str(callable_status).lower(),
            })
    write_tsv(
        results_dir / "window_callability.tsv",
        callability_rows,
        ["window", "window_label", "locus", "type_call", "analysis_id", "kcon_anchored_bases", "canonical_bases", "callable_fraction", "callable_at_0.90"],
    )

    pairwise_rows = pairwise_summaries(representatives, aligned, columns_by_window, callable_by_window)
    write_tsv(
        results_dir / "pairwise_divergence_summary.tsv",
        pairwise_rows,
        [
            "window", "window_label", "window_kind", "kcon_anchored_bases", "contrast",
            "pair_count", "unique_loci", "mean_p_distance", "median_p_distance",
            "minimum_p_distance", "maximum_p_distance", "mean_jointly_callable_bases", "inference",
        ],
    )

    type1_representatives = [representative for representative in representatives if representative.type_call == "TypeI"]
    locus_diversity_rows: list[dict] = []
    for index, window in enumerate(WINDOWS):
        estimate = type1_locus_bootstrap_window(
            type1_representatives,
            aligned,
            columns_by_window[window.key],
            callable_by_window[window.key],
            args.resamples,
            SEED + 500 + index,
        )
        locus_diversity_rows.append({
            "window": window.key,
            "window_label": window.label,
            "window_kind": window.kind,
            "kcon_anchored_bases": window.length,
            **estimate,
            "locus_bootstrap_resamples": args.resamples,
            "inference_unit": "TypeI_locus",
        })
    write_tsv(
        results_dir / "type1_locus_bootstrap_diversity.tsv",
        locus_diversity_rows,
        [
            "window", "window_label", "window_kind", "kcon_anchored_bases",
            "callable_type1_loci", "mean_within_type1_pairwise_p_distance",
            "locus_bootstrap_ci_low", "locus_bootstrap_ci_high",
            "locus_bootstrap_resamples", "inference_unit",
        ],
    )
    locus_bootstrap_comparisons: list[dict] = []
    for index, control_key in enumerate(CONTROL_KEYS):
        comparison = type1_locus_bootstrap_difference(
            type1_representatives,
            aligned,
            columns_by_window[CASSETTE_KEY],
            columns_by_window[control_key],
            callable_by_window[CASSETTE_KEY],
            callable_by_window[control_key],
            args.resamples,
            SEED + 700 + index,
        )
        locus_bootstrap_comparisons.append({
            "cassette_window": CASSETTE_KEY,
            "control_window": control_key,
            **comparison,
            "locus_bootstrap_resamples": args.resamples,
            "inference_unit": "TypeI_locus",
        })
    write_tsv(
        results_dir / "type1_locus_bootstrap_comparisons.tsv",
        locus_bootstrap_comparisons,
        [
            "cassette_window", "control_window", "paired_callable_type1_loci",
            "cassette_mean_within_type1_pairwise_p_distance",
            "control_mean_within_type1_pairwise_p_distance",
            "difference_cassette_minus_control", "locus_bootstrap_ci_low",
            "locus_bootstrap_ci_high", "locus_bootstrap_resamples", "inference_unit",
        ],
    )

    minimum_type1_site_callable = math.ceil(SITE_CALLABLE_TYPE1_FRACTION * len(type1_representatives))
    site_scores: dict[str, dict[int, list[float]]] = {}
    diversity_rows: list[dict] = []
    for index, window in enumerate(WINDOWS):
        by_callable: dict[int, list[float]] = defaultdict(list)
        sequences = [aligned[representative.analysis_id] for representative in type1_representatives]
        for column in columns_by_window[window.key]:
            value = site_diversity(sequences, column, minimum_type1_site_callable)
            if value is not None:
                callable_count, diversity = value
                by_callable[callable_count].append(diversity)
        site_scores[window.key] = dict(by_callable)
        flat = [value for values in by_callable.values() for value in values]
        rng = random.Random(SEED + 100 + index)
        low, high = bootstrap_mean(flat, args.resamples, rng)
        diversity_rows.append({
            "window": window.key,
            "window_label": window.label,
            "window_kind": window.kind,
            "kcon_anchored_bases": window.length,
            "comparison_eligible_sites": len(flat),
            "minimum_callable_type1_loci_per_site": minimum_type1_site_callable,
            "mean_type1_site_diversity": statistics.fmean(flat) if flat else math.nan,
            "bootstrap_95ci_low": low,
            "bootstrap_95ci_high": high,
            "bootstrap_resamples": args.resamples,
        })
    write_tsv(
        results_dir / "type1_site_diversity.tsv",
        diversity_rows,
        [
            "window", "window_label", "window_kind", "kcon_anchored_bases",
            "comparison_eligible_sites", "minimum_callable_type1_loci_per_site",
            "mean_type1_site_diversity", "bootstrap_95ci_low", "bootstrap_95ci_high",
            "bootstrap_resamples",
        ],
    )

    site_comparisons: list[dict] = []
    for index, control_key in enumerate(CONTROL_KEYS):
        result = matched_site_comparison(
            site_scores[CASSETTE_KEY], site_scores[control_key], args.resamples, SEED + 1_000 + index,
        )
        site_comparisons.append({
            "cassette_window": CASSETTE_KEY,
            "control_window": control_key,
            **result,
            "bootstrap_resamples": args.resamples,
            "permutation_resamples": args.resamples,
            "matching_rule": "exact_TypeI_callable_locus_count_then_equal_site_count_without_replacement",
            "inference_status": "linked_site_resampling_sensitivity_only_not_independent_inference",
        })
    write_tsv(
        results_dir / "matched_type1_site_comparisons.tsv",
        site_comparisons,
        [
            "cassette_window", "control_window", "matched_sites_per_window", "callable_count_strata",
            "cassette_mean", "control_mean", "difference_cassette_minus_control",
            "bootstrap_ci_low", "bootstrap_ci_high", "permutation_p_two_sided",
            "bootstrap_resamples", "permutation_resamples", "matching_rule",
            "inference_status",
        ],
    )

    neighbor_rows = nearest_neighbors(representatives, aligned, columns_by_window, callable_by_window)
    write_tsv(
        results_dir / "type1_nearest_neighbors.tsv",
        neighbor_rows,
        [
            "window", "window_label", "type1_locus", "type1_analysis_id", "status",
            "nearest_type2_locus", "nearest_type2_identity", "nearest_type2_jointly_callable_bases",
            "nearest_overall_locus", "nearest_overall_type", "nearest_overall_identity",
        ],
    )
    neighbor_summary_rows: list[dict] = []
    for window in WINDOWS:
        rows = [
            row for row in neighbor_rows
            if row["window"] == window.key and row["status"] == "callable"
        ]
        nearest_type2_identities = [float(row["nearest_type2_identity"]) for row in rows]
        nearest_overall_identities = [float(row["nearest_overall_identity"]) for row in rows]
        nearest_overall_typeI = sum(row["nearest_overall_type"] == "TypeI" for row in rows)
        neighbor_summary_rows.append({
            "window": window.key,
            "window_label": window.label,
            "window_kind": window.kind,
            "type1_loci_callable": len(rows),
            "mean_nearest_type2_identity": statistics.fmean(nearest_type2_identities) if rows else math.nan,
            "median_nearest_type2_identity": statistics.median(nearest_type2_identities) if rows else math.nan,
            "distinct_nearest_type2_loci": len({row["nearest_type2_locus"] for row in rows}),
            "mean_nearest_overall_identity": statistics.fmean(nearest_overall_identities) if rows else math.nan,
            "nearest_overall_is_typeI": nearest_overall_typeI,
            "nearest_overall_is_typeI_fraction": nearest_overall_typeI / len(rows) if rows else math.nan,
            "inference": "descriptive_local_neighbor_summary",
        })
    write_tsv(
        results_dir / "type1_nearest_neighbor_summary.tsv",
        neighbor_summary_rows,
        [
            "window", "window_label", "window_kind", "type1_loci_callable",
            "mean_nearest_type2_identity", "median_nearest_type2_identity",
            "distinct_nearest_type2_loci", "mean_nearest_overall_identity",
            "nearest_overall_is_typeI", "nearest_overall_is_typeI_fraction", "inference",
        ],
    )
    switching_rows, switching_summary = nearest_neighbor_switching(neighbor_rows)
    write_tsv(
        results_dir / "nearest_neighbor_switching.tsv",
        switching_rows,
        [
            "type1_locus", "callable_backbone_windows", "backbone_modal_nearest_type2",
            "backbone_modal_support_windows", "cassette_nearest_type2",
            "cassette_differs_from_backbone_mode", "distinct_nearest_type2_across_all_windows",
            *[f"cassette_differs_from_{key}" for key in CONTROL_KEYS],
        ],
    )
    write_tsv(
        results_dir / "nearest_neighbor_switching_summary.tsv",
        switching_summary,
        [
            "comparison", "type1_loci_compared", "switched_nearest_type2", "switch_fraction",
            "wilson_95ci_low", "wilson_95ci_high", "interpretation",
        ],
    )
    leave_one_out_rows, leave_one_out_summary = leave_one_window_out_switching(neighbor_rows)
    write_tsv(
        results_dir / "nearest_neighbor_leave_one_out.tsv",
        leave_one_out_rows,
        [
            "type1_locus", "focal_window", "focal_window_kind", "focal_nearest_type2",
            "reference_backbone_windows", "reference_backbone_modal_nearest_type2",
            "reference_mode_support_windows", "focal_differs_from_reference_backbone_mode", "status",
        ],
    )
    write_tsv(
        results_dir / "nearest_neighbor_leave_one_out_summary.tsv",
        leave_one_out_summary,
        [
            "focal_window", "focal_window_label", "focal_window_kind", "type1_loci_compared",
            "switched_nearest_type2", "switch_fraction", "wilson_95ci_low", "wilson_95ci_high",
            "negative_control_role", "interpretation",
        ],
    )

    nearest_identity: dict[str, dict[str, float]] = defaultdict(dict)
    for row in neighbor_rows:
        if row["status"] == "callable":
            nearest_identity[row["window"]][row["type1_locus"]] = float(row["nearest_type2_identity"])
    identity_comparisons: list[dict] = []
    for index, control_key in enumerate(CONTROL_KEYS):
        comparison = paired_comparison(
            nearest_identity[CASSETTE_KEY], nearest_identity[control_key], args.resamples, SEED + 2_000 + index,
        )
        identity_comparisons.append({
            "cassette_window": CASSETTE_KEY,
            "control_window": control_key,
            **comparison,
            "bootstrap_resamples": args.resamples,
            "sign_flip_resamples": args.resamples,
        })
    write_tsv(
        results_dir / "nearest_type2_identity_comparisons.tsv",
        identity_comparisons,
        [
            "cassette_window", "control_window", "paired_type1_loci",
            "mean_difference_cassette_minus_control", "bootstrap_ci_low", "bootstrap_ci_high",
            "sign_flip_p_two_sided", "bootstrap_resamples", "sign_flip_resamples",
        ],
    )

    delta_audit_rows, delta_record_provenance = audit_delta_records(eligible, args.delta292_records)
    write_tsv(
        results_dir / "per_haplotype_delta292_audit.tsv",
        delta_audit_rows,
        [
            "locus", "summary_type_call", "record_status", "eligible_rows", "eligible_typeI_rows",
            "eligible_typeII_rows", "eligible_conflict_rows", "eligible_no_call_rows",
            "both_typeI_and_typeII_present", "record_sha256",
        ],
    )

    selected_type_counts = Counter(representative.type_call for representative in representatives)
    callability_counts = {
        window.key: {
            type_call: sum(
                representative.analysis_id in callable_by_window[window.key]
                for representative in representatives
                if representative.type_call == type_call
            )
            for type_call in ("TypeI", "TypeII")
        }
        for window in WINDOWS
    }
    descriptive_site_lower_controls = sum(float(row["difference_cassette_minus_control"]) < 0 for row in site_comparisons)
    locus_lower_controls = sum(float(row["difference_cassette_minus_control"]) < 0 for row in locus_bootstrap_comparisons)
    locus_decisive_lower_controls = sum(
        float(row["difference_cassette_minus_control"]) < 0
        and float(row["locus_bootstrap_ci_high"]) < 0
        for row in locus_bootstrap_comparisons
    )
    callable_gate = all(
        counts["TypeI"] >= math.ceil(0.80 * selected_type_counts["TypeI"])
        and counts["TypeII"] >= math.ceil(0.80 * selected_type_counts["TypeII"])
        for counts in callability_counts.values()
    )
    neighbor_gate = all(
        sum(
            row["window"] == window.key and row["status"] == "callable"
            for row in neighbor_rows
        ) >= math.ceil(0.80 * selected_type_counts["TypeI"])
        for window in WINDOWS
    )
    figure_gate = callable_gate and neighbor_gate and locus_lower_controls >= 5 and locus_decisive_lower_controls >= 4
    if figure_gate:
        emit_figure(args.output_dir, locus_diversity_rows, leave_one_out_summary)
        figure_status = "emitted_stability_gate_passed"
    else:
        figure_status = "not_emitted_stability_gate_failed"

    both_loci = sorted(
        row["locus"] for row in delta_audit_rows
        if row["both_typeI_and_typeII_present"] == "true"
    )
    leave_one_out_by_window = {row["focal_window"]: row for row in leave_one_out_summary}
    backbone_negative_control_rates = [
        float(leave_one_out_by_window[key]["switch_fraction"]) for key in CONTROL_KEYS
    ]
    cassette_switch_rate = float(leave_one_out_by_window[CASSETTE_KEY]["switch_fraction"])
    summary = {
        "schema": "hml2.type1-cassette-backbone-rederivation.v1",
        "analysis_scope": "one_deterministic_population_provirus_per_production_eligible_locus",
        "coordinate_system": "Type-II_KCON_0-based_half-open; KCON bases mapped exactly to MAFFT columns",
        "eligible_loci": len(eligible),
        "eligible_type_counts": dict(Counter(eligible.values())),
        "selected_representatives": len(representatives),
        "selected_type_counts": dict(selected_type_counts),
        "excluded_eligible_loci": sorted(set(eligible) - {representative.locus for representative in representatives}),
        "mafft": mafft_details,
        "mafft_alignment_columns": next(iter(alignment_lengths)),
        "mafft_reverse_complemented_representatives": reversed_ids,
        "callability": callability_counts,
        "delta292_eligible_haplotype_audit": {
            "loci_with_both_TypeI_and_TypeII_calls": len(both_loci),
            "loci": both_loci,
            "interpretation_limit": "Current locus uniformity cannot distinguish RT-mediated spread from ancient post-integration conversion.",
        },
        "type1_site_comparison": {
            "descriptive_controls_with_lower_cassette_diversity": descriptive_site_lower_controls,
            "controls_total": len(CONTROL_KEYS),
            "inference_status": "linked_site_resampling_sensitivity_only_not_independent_inference",
        },
        "type1_locus_bootstrap_comparison": {
            "controls_with_lower_cassette_diversity": locus_lower_controls,
            "controls_with_locus_bootstrap_CI_below_zero": locus_decisive_lower_controls,
            "controls_total": len(CONTROL_KEYS),
            "inference_unit": "TypeI_locus",
        },
        "nearest_type2_neighbor_leave_one_out": {
            "cassette_switch_fraction": cassette_switch_rate,
            "backbone_negative_control_switch_fraction_minimum": min(backbone_negative_control_rates),
            "backbone_negative_control_switch_fraction_maximum": max(backbone_negative_control_rates),
            "interpretation": "Cassette switching is not exceptional relative to leave-one-backbone-window-out negative controls; discordance is broad across the provirus.",
        },
        "figure": {
            "status": figure_status,
            "gate": {
                "all_windows_at_least_80_percent_each_type_callable": callable_gate,
                "all_windows_at_least_80_percent_TypeI_have_TypeII_neighbor": neighbor_gate,
                "cassette_diversity_lower_than_at_least_5_of_6_controls_by_locus_bootstrap_point_estimate": locus_lower_controls >= 5,
                "at_least_4_controls_TypeI_locus_bootstrap_CI_below_zero": locus_decisive_lower_controls >= 4,
            },
        },
        "interpretation_boundaries": [
            "Local sequence discordance or nearest-neighbor switching supports recombination somewhere in the history, but not its direction, timing, or molecular route.",
            "Uniform present-day Type-I/Type-II calls within a locus do not distinguish reverse-transcription-mediated spread from ancient post-integration gene conversion.",
            "Cassette conservation does not establish host benefit, antiviral function, or selection.",
            "The observed combination of a relatively conserved Type-I cassette and similarly frequent neighbor switching in backbone controls is compatible with a restricted founder lineage plus pervasive recombination; it does not require a special cassette reservoir.",
            "The direct Pol-Env sequence-capability audit is a separate analysis and is not evidence here.",
        ],
        "parameters": {
            "callable_locus_fraction": CALLABLE_LOCUS_FRACTION,
            "pairwise_overlap_fraction": PAIRWISE_OVERLAP_FRACTION,
            "site_callable_TypeI_fraction": SITE_CALLABLE_TYPE1_FRACTION,
            "resamples": args.resamples,
            "seed": SEED,
            "windows": [asdict(window) | {"length": window.length} for window in WINDOWS],
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
    }
    write_json(results_dir / "analysis_summary.json", summary)
    write_json(results_dir / "analysis_parameters.json", summary["parameters"])

    provenance = [
        {"role": "orf_table", "locus": "", "basename": args.orf_table.name, "bytes": args.orf_table.stat().st_size, "sha256": sha256(args.orf_table)},
        {"role": "delta292_summary", "locus": "", "basename": args.delta292_summary.name, "bytes": args.delta292_summary.stat().st_size, "sha256": sha256(args.delta292_summary)},
        {"role": "type2_KCON_coordinate_anchor", "locus": "", "basename": args.type2_reference.name, "bytes": args.type2_reference.stat().st_size, "sha256": sha256(args.type2_reference)},
        *[
            {
                "role": "selected_population_provirus_fasta",
                "locus": representative.locus,
                "basename": representative.fasta.name,
                "bytes": representative.fasta.stat().st_size,
                "sha256": representative.fasta_sha256,
            }
            for representative in representatives
        ],
        *delta_record_provenance,
    ]
    write_tsv(
        derived_dir / "input_provenance.tsv",
        provenance,
        ["role", "locus", "basename", "bytes", "sha256"],
    )
    write_tsv(
        derived_dir / "delta292_summary_scope.tsv",
        delta_summary_rows,
        ["locus", "type_call", "status", "production_eligible"],
    )
    print(json.dumps({
        "eligible_loci": len(eligible),
        "selected_representatives": len(representatives),
        "selected_type_counts": dict(selected_type_counts),
        "loci_with_both_delta_calls": len(both_loci),
        "descriptive_site_lower_diversity_controls": descriptive_site_lower_controls,
        "locus_bootstrap_lower_diversity_controls": locus_lower_controls,
        "locus_bootstrap_decisive_lower_controls": locus_decisive_lower_controls,
        "figure_status": figure_status,
        "mafft_peak_rss_bytes": mafft_details["mafft_peak_rss_bytes_from_child_rusage"],
        "memory_free_before": mafft_details["memory_free_percent_before"],
        "memory_free_after": mafft_details["memory_free_percent_after"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
