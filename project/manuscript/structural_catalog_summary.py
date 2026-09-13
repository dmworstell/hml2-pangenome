"""Structural observations for Figure 1 and Tables S1/S5, with explicit units.

The catalog's noncarrier sentinel is not a confirmed empty-site allele. Missing
records are not noncarrier calls. Unlocalized copy buckets are not single loci.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

STATES = ("Noncarrier call", "Solo-LTR", "Fragment", "Provirus", "Multi-copy", "Unknown")
GROUP_LABELS = {
    "HML-2_acro_type1",
    "HML-2_acro_type2",
    "HML-2_8p23.1_duplicate_group_unresolved",
    "HML-2_Yq11.23_duplicate_group_unresolved",
}
PROJECT = Path(__file__).resolve().parents[1]
SAMPLE_FRAME = PROJECT / "inputs/sampling_frame/figure1_sample_sex.tsv"
XY_PARTITIONS = PROJECT / "inputs/population_authority/figure1_xy_partitions.json"


def structural_state(rows: list[dict[str, str]]) -> str:
    present = [row for row in rows if row["observation_state"] == "PRESENT"]
    if len(present) > 1:
        return "Multi-copy"
    if present:
        structure = present[0]["Structure"]
        if structure in {"Solo-LTR", "Fragment"}:
            return structure
        if structure in {"Provirus", "Provirus_from_Multi"}:
            return "Provirus"
        raise ValueError(f"unrecognized present structural state: {structure}")
    if not rows or any(row["observation_state"] == "UNKNOWN_TECHNICAL" for row in rows):
        return "Unknown"
    if all(row["observation_state"] == "NONCARRIER_SENTINEL_OBSERVATION" for row in rows):
        return "Noncarrier call"
    raise ValueError(f"unrecognized catalog observation states: {set(row['observation_state'] for row in rows)}")


def chromosome_eligibility(rows, roster, sample_frame=SAMPLE_FRAME, xy_partitions=XY_PARTITIONS):
    """Use documented sex and sequence-based partitions, not pat/mat semantics.

    The full descriptive panel retains relatives. Only the male X/Y partition
    evidence is used from the older authority, not its founder-only counts or
    its different structural-state calls.
    """
    with sample_frame.open(newline="") as handle:
        metadata = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    source = json.loads(xy_partitions.read_text())
    assignments = source["copy_assignment_contract"]["male_partition_assignments"]
    source_by_sample = defaultdict(list)
    for assignment in assignments:
        source_by_sample[assignment["sample_id"]].append(assignment)
    by_sample = defaultdict(set)
    digits = defaultdict(set)
    observed_x = defaultdict(set)
    for sample, hap in roster:
        by_sample[sample].add(hap)
    for row in rows:
        value = row.get("haplotype_pansn_digit", "")
        if value in {"1", "2"}:
            digits[(row["ID"], row["Haplotype"])].add(value)
        if row["Locus"].startswith("HML-2_X") and row["observation_state"] == "PRESENT":
            observed_x[row["ID"]].add(row["Haplotype"])
    eligibility = {}
    evidence = []
    for sample, haps in sorted(by_sample.items()):
        if sample not in metadata:
            raise ValueError(f"sample missing from retained metadata: {sample}")
        sex = metadata[sample]["sex"]
        if len(haps) != 2:
            raise ValueError(f"descriptive panel does not have two partitions: {sample}")
        if sex == "female":
            selected = {"X": set(haps), "Y": set()}
            method = "retained_female_metadata_two_X_no_Y"
        elif sex == "male":
            selected = {}
            for chrom in ("X", "Y"):
                candidates = set()
                for assignment in source_by_sample[sample]:
                    partition = assignment[f"{chrom}_partition"]
                    if partition in haps:
                        candidates.add(partition)
                    elif partition in {"h1", "hap1", "h2", "hap2"}:
                        candidates.update(hap for hap in haps if digits[(sample, hap)] == {partition[-1]})
                if len(candidates) != 1:
                    raise ValueError(f"unresolved retained {chrom} partition: {sample}, {candidates}")
                selected[chrom] = candidates
            if selected["X"] & selected["Y"]:
                raise ValueError(f"male X and Y assigned to same partition: {sample}")
            if observed_x[sample] != selected["X"]:
                raise ValueError(f"retained X partition disagrees with current catalog: {sample}")
            method = "retained_sequence_bearing_X_partition_and_complementary_Y_partition"
        elif sex == "unknown":
            selected = {"X": set(), "Y": set()}
            method = "sex_unknown_no_chromosome_count_inferred"
        else:
            raise ValueError(f"unrecognized metadata sex: {sample}, {sex}")
        for hap in sorted(haps):
            for chrom in ("X", "Y"):
                status = ("eligible" if hap in selected[chrom] else
                          "unknown_ploidy" if sex == "unknown" else "not_applicable")
                eligibility[(sample, hap, chrom)] = status
                evidence.append(dict(sample=sample, haplotype=hap, chromosome=chrom,
                                     sex=sex, eligibility=status, evidence=method))
    return eligibility, evidence


def structural_summary(rows, roster, sample_frame=SAMPLE_FRAME, xy_partitions=XY_PARTITIONS):
    eligibility, sex_evidence = chromosome_eligibility(rows, roster, sample_frame, xy_partitions)
    cells = defaultdict(list)
    original_cells = defaultdict(list)
    for row in rows:
        key = (row["Locus"], row["ID"], row["Haplotype"])
        cells[key].append(row)
        original_cells[(row.get("orig_Locus", row["Locus"]), row["ID"], row["Haplotype"])].append(row)
    summaries = []
    observations = []
    for locus in sorted({row["Locus"] for row in rows}):
        grouped = locus in GROUP_LABELS
        chromosome = locus.removeprefix("HML-2_")[0]
        counts = Counter()
        outside = Counter()
        reasons = Counter()
        max_copies = 0
        record_count = 0
        group_haplotypes = 0
        for sample, hap in roster:
            cell = cells.get((locus, sample, hap), [])
            status = eligibility[(sample, hap, chromosome)] if chromosome in {"X", "Y"} else "eligible"
            present = sum(row["observation_state"] == "PRESENT" for row in cell)
            if grouped:
                # This bucket collects copies that cannot be placed at one locus.
                # A sample without a bucket record has no defined "absent" call.
                if cell:
                    state = structural_state(cell)
                    counts[state] += 1
                    max_copies = max(max_copies, present)
                    record_count += len(cell)
                    group_haplotypes += 1
                observations.append(dict(locus=locus, sample=sample, haplotype=hap,
                                         eligibility="unlocalized_record_bucket", state=state if cell else "No bucket record",
                                         reason="not_a_locus_frequency_unit", present_records=present))
                continue
            if status != "eligible":
                outside[status] += 1
                observations.append(dict(locus=locus, sample=sample, haplotype=hap,
                                         eligibility=status, state="Outside denominator",
                                         reason="retained_sex_and_sequence_based_partition", present_records=present))
                continue
            state = structural_state(cell)
            counts[state] += 1
            max_copies = max(max_copies, present)
            if not cell:
                moved = original_cells.get((locus, sample, hap), [])
                reason = "records_reassigned_to_other_or_unlocalized_label" if moved else "no_retained_source_record"
            elif state == "Unknown":
                reason = "explicit_technical_unknown"
            elif state == "Noncarrier call":
                reason = "pipeline_noncarrier_sentinel_not_confirmed_empty_site"
            else:
                reason = "retained_structural_observation"
            reasons[reason] += 1
            observations.append(dict(locus=locus, sample=sample, haplotype=hap,
                                     eligibility=status, state=state, reason=reason, present_records=present))
        denominator = sum(counts.values()) if not grouped else 0
        called = denominator - counts["Unknown"] if not grouped else 0
        variability = 1 - max(counts[state] for state in STATES if state != "Unknown") / called if called else 0
        if grouped:
            entries = [f"{state.lower().replace('multi-copy', 'multiple records')} {counts[state]}" for state in STATES if counts[state]]
            text = f"Unlocalized records in {group_haplotypes} haplotypes, " + ", ".join(entries) + ". Not a single-locus frequency."
        else:
            entries = [f"{state.lower()} {counts[state]}/{denominator}" for state in STATES if counts[state]]
            text = ", ".join(entries)
            if outside:
                text += f". Outside denominator: {outside['not_applicable']} inapplicable partitions, {outside['unknown_ploidy']} unknown-sex partitions."
        summaries.append(dict(locus=locus, label_scope="unlocalized_record_bucket" if grouped else "physical_locus",
                              **{state: counts[state] for state in STATES}, eligible_haplotypes=denominator,
                              known_haplotypes=called, not_applicable_haplotypes=outside["not_applicable"],
                              unknown_ploidy_haplotypes=outside["unknown_ploidy"],
                              missing_reassigned=reasons["records_reassigned_to_other_or_unlocalized_label"],
                              missing_no_retained_record=reasons["no_retained_source_record"],
                              explicit_technical_unknown=reasons["explicit_technical_unknown"],
                              group_record_count=record_count, group_observed_haplotypes=group_haplotypes,
                              maximum_retained_copies_per_eligible_haplotype=max_copies,
                              variability_score=variability, structural_observations=text))
    return summaries, observations, sex_evidence


def write_structural_tables(summary, observations, sex_evidence, output):
    output.mkdir(parents=True, exist_ok=True)
    for filename, records in (
        ("Table_S5_artifact_filtered_structural_spectrum.tsv", summary),
        ("Table_S1_structural_observations_corrected.tsv", summary),
        ("Figure_1_structural_observation_cells.tsv", observations),
        ("Figure_1_sex_chromosome_eligibility.tsv", sex_evidence),
    ):
        with (output / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
