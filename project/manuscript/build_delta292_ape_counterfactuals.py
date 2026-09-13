#!/usr/bin/env python3
"""Test whether ape HML-2 sequences separate Delta292 from its human-linked haplotype.

Coordinates are zero-based KCON coordinates.  The frozen candidate ledger owns
the Type-I/Type-II classification; this script does not reclassify candidates.
Raw species observations are reported separately from ancestral integrations.
Orthologous observations are never interpreted as independent integrations.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[2]
PROJECT = WORKSPACE / "project"
OUT = PROJECT / "manuscript" / "delta292_ape_counterfactuals_v1"
FASTA = (
    PROJECT
    / "working"
    / "ape_recombination_origin_v1"
    / "inputs"
    / "ape_candidates_kcon_projection.fa"
)
LEDGER = (
    PROJECT
    / "working"
    / "ape_genomewide_hml2_scan_v1"
    / "results"
    / "candidate_ledger.tsv"
)
UNITS = (
    PROJECT
    / "working"
    / "claude_safe_delta292_why_execution_v1"
    / "results"
    / "exact_orthology_collapsed_lesion_units.tsv"
)

SITES = (6328, 6331, 6485, 6492, 6495, 6496, 6498)
FIVE_SITES = (6328, 6331, 6492, 6496, 6498)
HUMAN_TYPEI_FIVE_SITE = "TTTAC"
HUMAN_MINIMAL_PAIR = "TT"
CANONICAL = frozenset("ACGT")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = ""
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        else:
            if not name:
                raise ValueError(f"sequence before header in {path}")
            records[name].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    sequences = read_fasta(FASTA)
    ledger = read_tsv(LEDGER)
    if len(sequences) != 83 or len(ledger) != 83:
        raise ValueError("expected the frozen 83-candidate ape panel")
    if set(sequences) != {row["candidate_id"] for row in ledger}:
        raise ValueError("FASTA headers and candidate ledger do not match")
    if {len(sequence) for sequence in sequences.values()} != {9472}:
        raise ValueError("ape projections must all be 9472 KCON columns")

    record_to_unit: dict[str, str] = {}
    unit_status: dict[str, str] = {}
    for unit in read_tsv(UNITS):
        unit_status[unit["ancestral_integration_unit_id"]] = unit["resolution_status"]
        for record in unit["member_records"].split(";"):
            if record:
                record_to_unit[record] = unit["ancestral_integration_unit_id"]

    calls: list[dict[str, object]] = []
    for source in ledger:
        candidate = source["candidate_id"]
        sequence = sequences[candidate]
        site_calls = {position: sequence[position] for position in SITES}
        five_site = "".join(site_calls[position] for position in FIVE_SITES)
        pair = site_calls[6331] + site_calls[6492]
        local_positions = tuple(range(6300, 6501)) + tuple(range(6793, 6851))
        local_callable = sum(sequence[position] in CANONICAL for position in local_positions)
        delta = source["provirus_type"] == "TypeI_canonical_Delta292"
        full_haplotype = five_site == HUMAN_TYPEI_FIVE_SITE
        minimal_pair = pair == HUMAN_MINIMAL_PAIR
        unit = record_to_unit.get(candidate, "")
        calls.append(
            {
                "candidate_id": candidate,
                "species": source["species"],
                "provirus_type": source["provirus_type"],
                "best_identity": source["best_identity"],
                "query_union_coverage": source["query_union_coverage"],
                "deletion_events_kcon_start_length": source[
                    "deletion_events_kcon_start_length"
                ],
                "ancestral_integration_unit_id": unit or "not_resolved_in_lesion_unit_table",
                "unit_resolution_status": unit_status.get(
                    unit, "not_resolved_in_lesion_unit_table"
                ),
                "alleles_6328_6331_6485_6492_6495_6496_6498": "".join(
                    site_calls[position] for position in SITES
                ),
                "five_site_haplotype_6328_6331_6492_6496_6498": five_site,
                "pair_6331_6492": pair,
                "carries_human_typeI_five_site_haplotype": str(full_haplotype).lower(),
                "carries_human_minimal_pair": str(minimal_pair).lower(),
                "exact_delta292_frozen_call": str(delta).lower(),
                "retained_local_callable_sites": local_callable,
                "retained_local_total_sites": len(local_positions),
                "retained_local_callable_fraction": f"{local_callable / len(local_positions):.6f}",
                "natural_counterfactual_class": (
                    "delta292_without_human_minimal_pair"
                    if delta and not minimal_pair
                    else (
                        "delta292_without_human_five_site_haplotype"
                        if delta and not full_haplotype
                        else (
                            "human_minimal_pair_without_delta292"
                            if not delta and minimal_pair
                            else (
                                "human_five_site_haplotype_without_delta292"
                                if not delta and full_haplotype
                                else "none"
                            )
                        )
                    )
                ),
            }
        )

    state_counts: list[dict[str, object]] = []
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in calls:
        grouped[(str(row["provirus_type"]), str(row["species"]))].append(row)
    for (provirus_type, species), rows in sorted(grouped.items()):
        haplotypes = Counter(
            str(row["five_site_haplotype_6328_6331_6492_6496_6498"])
            for row in rows
        )
        state_counts.append(
            {
                "provirus_type": provirus_type,
                "species": species,
                "n_observations": len(rows),
                "n_human_five_site_haplotype": sum(
                    row["carries_human_typeI_five_site_haplotype"] == "true"
                    for row in rows
                ),
                "n_human_minimal_pair": sum(
                    row["carries_human_minimal_pair"] == "true" for row in rows
                ),
                "five_site_haplotype_counts": ";".join(
                    f"{haplotype}:{count}"
                    for haplotype, count in sorted(haplotypes.items())
                ),
                "count_semantics": (
                    "species observations; orthologous copies are not independent integrations"
                ),
            }
        )

    counterfactuals = [
        row for row in calls if row["natural_counterfactual_class"] != "none"
    ]
    delta_rows = [
        row for row in calls if row["provirus_type"] == "TypeI_canonical_Delta292"
    ]
    type2_rows = [row for row in calls if row["provirus_type"] == "TypeII_retained"]
    alternative_rows = [
        row for row in calls if row["provirus_type"] == "alternative_pol_env_deletion"
    ]

    # These assertions freeze the empirical linkage breaks that motivated this analysis.
    pair_without_delta = [
        row
        for row in calls
        if row["exact_delta292_frozen_call"] == "false"
        and row["carries_human_minimal_pair"] == "true"
    ]
    delta_without_pair = [
        row
        for row in calls
        if row["exact_delta292_frozen_call"] == "true"
        and row["carries_human_minimal_pair"] == "false"
    ]
    if len(pair_without_delta) != 4 or len(delta_without_pair) != 1:
        raise ValueError("ape natural-counterfactual counts changed")

    conclusions = [
        {
            "question": "Is the human 6331T+6492T pair sufficient for Delta292?",
            "result": "no",
            "evidence": (
                "Three orangutan Type-II candidates and one orangutan 112-nt "
                "alternative-deletion candidate carry 6331T+6492T without Delta292."
            ),
            "boundary": (
                "These are distinct same-species genomic coordinates, not repeated "
                "observations of one orthologous integration."
            ),
        },
        {
            "question": "Is the human 6331T+6492T pair universally required for a resident Delta292 element?",
            "result": "no",
            "evidence": (
                "The siamang exact-Delta292 candidate carries 6331C+6492G and the "
                "five-site haplotype CCGAC."
            ),
            "boundary": (
                "This candidate has lower identity/coverage and unresolved cross-species "
                "orthology, but all tested retained sites are directly callable. Its state "
                "can reflect ancestral sequence or post-integration substitution and therefore "
                "does not by itself identify the founding haplotype."
            ),
        },
        {
            "question": "Does the full human five-site haplotype occur naturally without Delta292 in this ape panel?",
            "result": "no_observation",
            "evidence": (
                "No Type-II or alternative-deletion candidate carries the full TTTAC "
                "five-site state."
            ),
            "boundary": (
                "Absence in 57 non-Delta292 ape candidates is an ascertainment-limited "
                "observation, not proof of functional dependence."
            ),
        },
        {
            "question": "Does the ape panel remove the human-only deletion-versus-linked-pair identifiability barrier?",
            "result": "yes_for_the_6331_6492_pair",
            "evidence": (
                "Both directions of natural separation are observed: pair without Delta292 "
                "and Delta292 without the pair."
            ),
            "boundary": (
                "The panel still does not isolate the founding viral genotype or measure a "
                "replicative effect of the deletion."
            ),
        },
    ]

    OUT.mkdir(parents=True, exist_ok=True)
    call_fields = list(calls[0])
    write_tsv(OUT / "candidate_genotypes.tsv", calls, call_fields)
    write_tsv(
        OUT / "natural_counterfactuals.tsv",
        counterfactuals,
        call_fields,
    )
    write_tsv(
        OUT / "species_stratified_state_counts.tsv",
        state_counts,
        list(state_counts[0]),
    )
    write_tsv(
        OUT / "validated_conclusions.tsv",
        conclusions,
        list(conclusions[0]),
    )

    summary = {
        "analysis": "ape-wide Delta292/local-haplotype natural counterfactual search",
        "coordinate_system": "zero-based KCON projection",
        "input_hashes": {
            str(FASTA.relative_to(WORKSPACE)): sha256(FASTA),
            str(LEDGER.relative_to(WORKSPACE)): sha256(LEDGER),
            str(UNITS.relative_to(WORKSPACE)): sha256(UNITS),
        },
        "panel": {
            "total": len(calls),
            "delta292": len(delta_rows),
            "typeII": len(type2_rows),
            "alternative_pol_env_deletion": len(alternative_rows),
        },
        "counterfactuals": {
            "non_delta_with_6331T_6492T": len(pair_without_delta),
            "typeII_non_delta_with_6331T_6492T": sum(
                row["provirus_type"] == "TypeII_retained" for row in pair_without_delta
            ),
            "delta292_without_6331T_6492T": len(delta_without_pair),
            "non_delta_with_full_TTTAC": sum(
                row["exact_delta292_frozen_call"] == "false"
                and row["carries_human_typeI_five_site_haplotype"] == "true"
                for row in calls
            ),
        },
        "main_result": (
            "The ape panel breaks the human-only alias between Delta292 and the "
            "6331T+6492T pair. The pair is neither sufficient for Delta292 nor "
            "universally required for a resident Delta292 element."
        ),
        "interpretive_boundary": (
            "This rejects the pair as a universal deterministic explanation. It does "
            "not establish whether the exact deletion, a broader founding haplotype, "
            "or source-lineage opportunity caused the historical expansion."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
