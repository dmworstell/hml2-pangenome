#!/usr/bin/env python3
"""Build the ORF table with duplicated loci assigned only when host flanks support them."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


PRESENT = {"Provirus", "Solo-LTR", "Fragment", "Provirus_from_Multi"}
ACRO_FAMILY = {
    "HML-2_13p13": "HML-2_acro_type1",
    "HML-2_15p13b": "HML-2_acro_type1",
    "HML-2_15p13a": "HML-2_acro_type2",
    "HML-2_21p13": "HML-2_acro_type2",
    "HML-2_22p13": "HML-2_acro_type2",
}
TYPE1_EVENT = "type1_acrocentric_ancestral_integration"
TYPE2_EVENT = "type2_acrocentric_4q35.2_ancestral_integration"
FOURQ = "HML-2_4q35.2_hg38"
FOURP = "HML-2_4p16.3a"
FOURQ_ANCHOR_EVIDENCE = (
    "same_alignment_reaches_80kb_or_more_into_chromosome4_interior"
)
FOURQ_GROUP = {"15p13a", "21p13", "22p13", "4q35.2_hg38"}
DETAIL_RE = re.compile(
    r"([^;|]+)\|chr([^:]+):(\d+)-(\d+)\|q:(\d+)-(\d+)/(\d+)"
    r"\|([+-])\|mapq:(\d+)\|(nonsecondary|secondary)"
)
EIGHTQ_ALIASES = {"HML-2_8q24.3b", "HML-2_8q24.3b_hg38"}
ASM_DUP_ID = "HG00658_pat_hprc_r2_v1.0.1_HML-2_7p22.1_asmdup"
DROP_ROW = "__DROP_ROW__"
EXACT_CANDIDATE_RE = re.compile(
    r"(?:t2t|hg38):([^:;]+):primary=(\d):mapq=(\d+)"
    r":source=[^:;]+:target=[^:;]+:element=[^:;]+"
    r":both5=(\d):int100=(\d):int1m=(\d):flanks=(\d+),(\d+)"
)
UNRESOLVED_DUPLICATE_FAMILIES = {
    frozenset({"HML-2_8p23.1b", "HML-2_8p23.1c", "HML-2_8p23.1d", "HML-2_8p23.1e"}): (
        "HML-2_8p23.1_duplicate_group_unresolved",
        "8p23.1_segmental_duplication_group",
    ),
    frozenset({"HML-2_Yq11.23a", "HML-2_Yq11.23b"}): (
        "HML-2_Yq11.23_duplicate_group_unresolved",
        "Yq11.23_duplication_group",
    ),
}


def source_qname(source_identifier: str) -> str:
    if not source_identifier or source_identifier in {
        "Placeholder",
        "Validator_Script",
        "Successful_Empty_Site",
    }:
        return ""
    return source_identifier.rsplit(":", 1)[0]


def locus_name(locus: str) -> str:
    return locus[6:] if locus.startswith("HML-2_") else locus


def source_interval(source_identifier: str) -> tuple[str, int, int] | None:
    try:
        qname, span = source_identifier.rsplit(":", 1)
        start_text, end_text = span.split("-", 1)
        return qname, int(start_text), int(end_text)
    except (AttributeError, ValueError):
        return None


def load_attachment_evidence(
    path: Path,
) -> tuple[set[str], dict[str, list[dict[str, object]]]]:
    qnames = set()
    alignments: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["locus"] == "4q35.2_hg38":
                for column in ("interior100_qnames", "interior1m_qnames"):
                    qnames.update(qname for qname in row[column].split(",") if qname)
            if row["locus"] not in FOURQ_GROUP:
                continue
            for match in DETAIL_RE.finditer(row["details"]):
                qname = match.group(1)
                alignments[qname].append(
                    {
                        "locus": row["locus"],
                        "ref_start": int(match.group(3)),
                        "ref_end": int(match.group(4)),
                        "query_start": int(match.group(5)),
                        "query_end": int(match.group(6)),
                        "mapq": int(match.group(9)),
                        "primary": match.group(10) == "nonsecondary",
                    }
                )
    if not qnames:
        raise ValueError("attachment audit contains no chromosome 4 anchor QNAMEs")
    return qnames, alignments


def load_cnv_decisions(
    weights_path: Path, receipt_path: Path
) -> tuple[set[str], dict[str, tuple[str, str]]]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    id_by_candidate = receipt["authority_exact_controller_replay_map"]
    artifact_ids = {ASM_DUP_ID}
    decisions = {ASM_DUP_ID: ("assembly_artifact", "")}
    with weights_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            id_full = id_by_candidate[row["candidate_key"]]
            state = row["dominant_state"]
            decisions[id_full] = (state, row.get(state, ""))
            if row["dominant_state"] == "assembly_artifact":
                artifact_ids.add(id_full)
    return artifact_ids, decisions


def load_exact_source_decisions(path: Path) -> dict[str, tuple[str, str]]:
    """Resolve exact duplicate labels only when host alignments distinguish them."""
    rows_by_source: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if "," in row["same_source_loci"]:
                rows_by_source[row["Source_Identifier"]].append(row)

    decisions = {}
    for source, rows in rows_by_source.items():
        loci = set(rows[0]["same_source_loci"].split(","))
        if loci <= {"HML-2_Yq11.23a", "HML-2_Yq11.23b"}:
            decisions[source] = (
                "HML-2_Yq11.23_duplicate_group_unresolved",
                "same_interval_has_no_host_flank_alignment_that_distinguishes_Yq11.23a_from_Yq11.23b",
            )
            continue

        candidates = []
        for candidate in rows[0]["all_candidates"].split(";"):
            match = EXACT_CANDIDATE_RE.fullmatch(candidate)
            if not match:
                continue
            candidates.append(
                {
                    "locus": f"HML-2_{match.group(1)}",
                    "primary": match.group(2) == "1",
                    "mapq": int(match.group(3)),
                    "both5": match.group(4) == "1",
                    "interior100": match.group(5) == "1",
                    "interior1m": match.group(6) == "1",
                    "left_flank": int(match.group(7)),
                    "right_flank": int(match.group(8)),
                }
            )

        if all("8p23.1" in locus for locus in loci):
            long_anchors = [
                candidate
                for candidate in candidates
                if candidate["primary"] and candidate["interior1m"]
            ]
            local_anchors = [
                candidate
                for candidate in candidates
                if candidate["primary"]
                and candidate["both5"]
                and min(candidate["left_flank"], candidate["right_flank"]) >= 5_000
            ]
            resolved = None
            evidence = ""
            if len(long_anchors) == 1:
                resolved = long_anchors[0]
                evidence = "one_alignment_links_the_element_to_at_least_1Mb_of_locus_specific_host_sequence"
            elif len(candidates) == 1 and len(local_anchors) == 1:
                resolved = local_anchors[0]
                evidence = "one_alignment_links_the_element_and_both_5kb_flanks_to_a_single_8p23.1_copy"
            if resolved is not None:
                decisions[source] = (str(resolved["locus"]), evidence)
            else:
                decisions[source] = (
                    "HML-2_8p23.1_duplicate_group_unresolved",
                    "same_interval_maps_to_more_than_one_8p23.1_duplicate_block_without_one_decisive_long_host_anchor",
                )
            continue

        t2t_anchors = [
            candidate
            for candidate in candidates
            if candidate["locus"] != FOURQ
            and candidate["primary"]
            and candidate["both5"]
        ]
        t2t_loci = {candidate["locus"] for candidate in t2t_anchors}
        if len(t2t_loci) == 1:
            target = next(iter(t2t_loci))
            decisions[source] = (
                target,
                "exact_CIGAR_alignment_of_the_element_and_both_host_flanks_to_"
                + locus_name(target),
            )
            continue

        family = unresolved_duplicate_family(loci)
        if family is not None:
            generic_locus, _ = family
            decisions[source] = (
                generic_locus,
                "same_interval_cannot_be_assigned_to_one_member_of_the_duplicate_group",
            )
    return decisions


def build_source_loci(rows: list[dict[str, str]]) -> dict[str, set[str]]:
    loci_by_source: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["Structure"] in PRESENT and source_interval(row["Source_Identifier"]):
            loci_by_source[row["Source_Identifier"]].add(locus_name(row["Locus"]))
    return loci_by_source


def fourq_alignment_candidates(
    row: dict[str, str], alignments: dict[str, list[dict[str, object]]]
) -> list[tuple[str, tuple[int, int, int]]]:
    parsed = source_interval(row.get("Source_Identifier", ""))
    if not parsed:
        return []
    qname, source_start, source_end = parsed
    source_length = source_end - source_start
    best_by_locus: dict[str, tuple[int, int, int]] = {}
    for alignment in alignments.get(qname, []):
        if not alignment["primary"]:
            continue
        query_start = int(alignment["query_start"])
        query_end = int(alignment["query_end"])
        overlap = max(
            0,
            min(source_end, query_end) - max(source_start, query_start),
        )
        if source_length <= 0 or overlap < 0.8 * source_length:
            continue
        ref_span = abs(int(alignment["ref_end"]) - int(alignment["ref_start"]))
        flank_span = max(0, source_start - query_start) + max(
            0, query_end - source_end
        )
        host_span = max(flank_span, ref_span - source_length)
        score = (int(alignment["mapq"]), host_span, ref_span)
        locus = str(alignment["locus"])
        if score > best_by_locus.get(locus, (-1, -1, -1)):
            best_by_locus[locus] = score
    return sorted(best_by_locus.items(), key=lambda item: item[1], reverse=True)


def resolve_fourq_rows(
    rows: list[dict[str, str]],
    alignments: dict[str, list[dict[str, object]]],
) -> dict[str, tuple[str, str]]:
    loci_by_source = build_source_loci(rows)
    decisions: dict[str, tuple[str, str]] = {}
    for row in rows:
        if row["Locus"] != FOURQ or row["Structure"] not in PRESENT:
            continue
        source = row.get("Source_Identifier", "")
        competing_loci = loci_by_source.get(source, set()) - {locus_name(FOURQ)}
        candidates = fourq_alignment_candidates(row, alignments)
        score_by_locus = dict(candidates)

        competing_acro = sorted(competing_loci & (FOURQ_GROUP - {"4q35.2_hg38"}))
        supported_competitors = [
            locus for locus in competing_acro if locus in score_by_locus
        ]
        if len(supported_competitors) == 1:
            target = supported_competitors[0]
            decisions[row["ID_Full"]] = (
                f"HML-2_{target}",
                f"same_element_interval_and_primary_host_sequence_place_it_at_{target}",
            )
            continue

        fourq_score = score_by_locus.get("4q35.2_hg38")
        if competing_loci - FOURQ_GROUP and (
            fourq_score is None or fourq_score[0] < 20
        ):
            decisions[row["ID_Full"]] = (
                DROP_ROW,
                "same_element_interval_belongs_to_another_catalog_locus",
            )
            continue

        if candidates and candidates[0][1][0] >= 20:
            best_locus, best_score = candidates[0]
            second_score = candidates[1][1] if len(candidates) > 1 else None
            clear = second_score is None or (
                best_score[0] >= second_score[0] + 10
                or best_score[1] >= second_score[1] + 20_000
            )
            if clear:
                if best_locus == "4q35.2_hg38":
                    evidence = "primary_alignment_of_element_and_host_sequence_to_chromosome4"
                else:
                    evidence = (
                        f"primary_host_alignment_places_this_record_at_{best_locus}_not_4q35.2"
                    )
                decisions[row["ID_Full"]] = (f"HML-2_{best_locus}", evidence)
                continue

        decisions[row["ID_Full"]] = (
            "HML-2_acro_type2",
            "type2_element_block_without_a_unique_chromosome_assignment",
        )
    return decisions


def annotate_row(
    row: dict[str, str],
    fourq_anchors: set[str],
    fourq_decisions: dict[str, tuple[str, str]] | None = None,
    cnv_decisions: dict[str, tuple[str, str]] | None = None,
    exact_source_decisions: dict[str, tuple[str, str]] | None = None,
) -> dict[str, str]:
    original_locus = row["Locus"]
    present = row["Structure"] in PRESENT
    row = dict(row)
    row["orig_Locus"] = original_locus
    row["physical_locus_assignment"] = locus_name(original_locus)
    row["physical_assignment_evidence"] = "original_catalog_assignment"
    row["insertion_event_group"] = locus_name(original_locus)
    state, probability = (cnv_decisions or {}).get(
        row["ID_Full"], ("not_assayed", "")
    )
    row["cnv_qc_state"] = state
    row["cnv_qc_probability"] = probability
    row["analysis_include"] = "1"
    row["analysis_exclusion_reason"] = ""
    row["record_relationship"] = "primary_record"
    row["representative_ID_Full"] = row["ID_Full"]

    if state == "assembly_artifact":
        row["analysis_include"] = "0"
        row["analysis_exclusion_reason"] = "assembly_artifact_not_supported_by_CNV_depth"
        row["record_relationship"] = "assembly_artifact"

    if original_locus in EIGHTQ_ALIASES:
        row["Locus"] = "HML-2_8q24.3c"
        row["physical_locus_assignment"] = "8q24.3c"
        row["physical_assignment_evidence"] = (
            "8q24.3b_is_an_alias_of_8q24.3c_not_a_separate_insertion"
        )
        row["insertion_event_group"] = "8q24.3c"
        row["analysis_include"] = "0"
        row["analysis_exclusion_reason"] = "alias_duplicate_of_8q24.3c"
        row["record_relationship"] = "alias_of_8q24.3c"
        return row

    if present:
        exact_decision = (exact_source_decisions or {}).get(
            row.get("Source_Identifier", "")
        )
        if exact_decision:
            target, evidence = exact_decision
            row["Locus"] = target
            row["physical_locus_assignment"] = locus_name(target)
            row["physical_assignment_evidence"] = evidence
            if target in {"HML-2_13p13", "HML-2_15p13b", "HML-2_acro_type1"}:
                row["insertion_event_group"] = TYPE1_EVENT
            elif target in {
                "HML-2_15p13a",
                "HML-2_21p13",
                "HML-2_22p13",
                FOURQ,
                "HML-2_acro_type2",
            }:
                row["insertion_event_group"] = TYPE2_EVENT
            elif "8p23.1" in target:
                row["insertion_event_group"] = "8p23.1_segmental_duplication_group"
            elif "Yq11.23" in target:
                row["insertion_event_group"] = "Yq11.23_duplication_group"
            return row

    if original_locus == FOURQ and present:
        decision = (fourq_decisions or {}).get(row["ID_Full"])
        if decision:
            target, evidence = decision
            if target == DROP_ROW:
                target = "HML-2_acro_type2"
                evidence = "same_element_interval_belongs_to_the_type2_duplicate_family_not_4q35.2"
            row["Locus"] = target
            row["physical_locus_assignment"] = locus_name(target)
            row["physical_assignment_evidence"] = evidence
            row["insertion_event_group"] = TYPE2_EVENT
            return row
        qname = source_qname(row.get("Source_Identifier", ""))
        if row.get("ID") == "GCA":
            row["physical_assignment_evidence"] = "GRCh38_reference_coordinate"
        elif qname in fourq_anchors:
            row["physical_assignment_evidence"] = FOURQ_ANCHOR_EVIDENCE
        else:
            row["Locus"] = "HML-2_acro_type2"
            row["physical_locus_assignment"] = "acrocentric_type2_unresolved"
            row["physical_assignment_evidence"] = (
                "element_and_local_duplicate_block_only_without_chromosome4_interior_anchor"
            )
        row["insertion_event_group"] = TYPE2_EVENT
        return row

    if original_locus == FOURP and present:
        qname = source_qname(row.get("Source_Identifier", ""))
        if qname in fourq_anchors:
            row["physical_assignment_evidence"] = (
                "chromosome4_contig_reaches_from_this_p_arm_record_into_q_arm_interior"
            )
        return row

    if original_locus in ACRO_FAMILY:
        if not present:
            row["insertion_event_group"] = (
                TYPE1_EVENT
                if ACRO_FAMILY[original_locus] == "HML-2_acro_type1"
                else TYPE2_EVENT
            )
            return row
        row["Locus"] = ACRO_FAMILY[original_locus]
        row["physical_locus_assignment"] = (
            "acrocentric_type1_unresolved"
            if row["Locus"] == "HML-2_acro_type1"
            else "acrocentric_type2_unresolved"
        )
        row["physical_assignment_evidence"] = (
            "distal_acrocentric_duplicate_family_without_reliable_short_arm_assignment"
        )
        row["insertion_event_group"] = (
            TYPE1_EVENT if row["Locus"] == "HML-2_acro_type1" else TYPE2_EVENT
        )
    return row


def unresolved_duplicate_family(
    loci: set[str],
) -> tuple[str, str] | None:
    for members, family in UNRESOLVED_DUPLICATE_FAMILIES.items():
        if loci.issubset(members):
            return family
    return None


def annotate_duplicate_present_rows(
    rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], int, int]:
    present_by_source = defaultdict(list)
    for index, row in enumerate(rows):
        if row["Structure"] in PRESENT and row.get("Source_Identifier"):
            present_by_source[row["Source_Identifier"]].append(index)

    duplicate_labels_excluded = 0
    unresolved_cross_family = 0
    for indices in present_by_source.values():
        if len(indices) < 2:
            continue
        eligible = [
            index
            for index in indices
            if rows[index]["analysis_include"] == "1"
        ]
        if not eligible:
            continue
        loci = {rows[index]["Locus"] for index in eligible}
        if len(loci) > 1:
            family = unresolved_duplicate_family(loci)
            if family is not None:
                generic_locus, event_group = family
                for index in eligible:
                    rows[index]["Locus"] = generic_locus
                    rows[index]["physical_locus_assignment"] = locus_name(generic_locus)
                    rows[index]["physical_assignment_evidence"] = (
                        "same_interval_cannot_be_assigned_to_one_member_of_the_duplicate_group"
                    )
                    rows[index]["insertion_event_group"] = event_group
            else:
                anchored = [
                    index
                    for index in eligible
                    if rows[index]["physical_assignment_evidence"].startswith(
                        (
                            "exact_CIGAR",
                            "one_alignment",
                            "same_element_interval",
                            "primary_host_alignment",
                            "primary_alignment",
                            "chromosome4_contig",
                        )
                    )
                ]
                if len({rows[index]["Locus"] for index in anchored}) != 1:
                    unresolved_cross_family += 1

        anchored = [
            index
            for index in eligible
            if rows[index]["physical_assignment_evidence"].startswith(
                (
                    "exact_CIGAR",
                    "one_alignment",
                    "same_element_interval",
                    "primary_host_alignment",
                    "primary_alignment",
                    "chromosome4_contig",
                )
            )
        ]
        if len({rows[index]["Locus"] for index in anchored}) == 1 and anchored:
            resolved = rows[anchored[0]]
            for index in eligible:
                rows[index]["Locus"] = resolved["Locus"]
                rows[index]["physical_locus_assignment"] = resolved[
                    "physical_locus_assignment"
                ]
                rows[index]["physical_assignment_evidence"] = resolved[
                    "physical_assignment_evidence"
                ]
                rows[index]["insertion_event_group"] = resolved[
                    "insertion_event_group"
                ]

        cnv_supported = [
            index
            for index in eligible
            if rows[index]["cnv_qc_state"] in {"authenticated_segdup", "later_duplication"}
        ]
        pool = cnv_supported if cnv_supported else eligible
        keep = max(
            pool,
            key=lambda index: (
                rows[index]["orig_Locus"] == rows[index]["Locus"],
                rows[index]["physical_assignment_evidence"].startswith(
                    (
                        "exact_CIGAR",
                        "one_alignment",
                        "same_element_interval",
                        "primary_host_alignment",
                        "primary_alignment",
                        "chromosome4_contig",
                    )
                ),
                rows[index]["physical_assignment_evidence"] == FOURQ_ANCHOR_EVIDENCE,
                rows[index]["ID_Full"],
            ),
        )
        representative = rows[keep]["ID_Full"]
        resolved = rows[keep]
        for index in indices:
            rows[index]["Locus"] = resolved["Locus"]
            rows[index]["physical_locus_assignment"] = resolved[
                "physical_locus_assignment"
            ]
            rows[index]["physical_assignment_evidence"] = resolved[
                "physical_assignment_evidence"
            ]
            rows[index]["insertion_event_group"] = resolved["insertion_event_group"]
            rows[index]["representative_ID_Full"] = representative
            if index == keep or rows[index]["analysis_include"] == "0":
                continue
            rows[index]["analysis_include"] = "0"
            rows[index]["analysis_exclusion_reason"] = (
                "duplicate_catalog_label_for_same_assembled_interval"
            )
            rows[index]["record_relationship"] = "duplicate_label"
            duplicate_labels_excluded += 1
    return rows, duplicate_labels_excluded, unresolved_cross_family


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--attachment-audit", required=True)
    parser.add_argument("--exact-assignment-audit", required=True)
    parser.add_argument("--weights", required=True)
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    summary_path = Path(args.summary)
    fourq_anchors, alignments = load_attachment_evidence(Path(args.attachment_audit))
    artifact_ids, cnv_decisions = load_cnv_decisions(
        Path(args.weights), Path(args.receipt)
    )
    exact_source_decisions = load_exact_source_decisions(
        Path(args.exact_assignment_audit)
    )

    with input_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_fields = list(reader.fieldnames or ())
        input_rows = list(reader)
    required = {"Locus", "ID_Full", "ID", "Source_Identifier", "Structure"}
    if not required.issubset(input_fields):
        raise ValueError(f"ORF table missing columns: {sorted(required - set(input_fields))}")

    artifact_rows = [row for row in input_rows if row["ID_Full"] in artifact_ids]
    clean_for_assignment = [row for row in input_rows if row["Locus"] not in EIGHTQ_ALIASES]
    fourq_decisions = resolve_fourq_rows(clean_for_assignment, alignments)
    rows = []
    alias_rows = 0
    for row in input_rows:
        annotated = annotate_row(
            row,
            fourq_anchors,
            fourq_decisions=fourq_decisions,
            cnv_decisions=cnv_decisions,
            exact_source_decisions=exact_source_decisions,
        )
        if row["Locus"] in EIGHTQ_ALIASES:
            alias_rows += 1
        rows.append(annotated)
    rows, duplicate_rows, unresolved_source_conflicts = annotate_duplicate_present_rows(rows)

    output_fields = input_fields + [
        field
        for field in (
            "orig_Locus",
            "physical_locus_assignment",
            "physical_assignment_evidence",
            "insertion_event_group",
            "cnv_qc_state",
            "cnv_qc_probability",
            "analysis_include",
            "analysis_exclusion_reason",
            "record_relationship",
            "representative_ID_Full",
        )
        if field not in input_fields
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    present_fourq = [
        row for row in rows if row["Locus"] == FOURQ and row["Structure"] in PRESENT
    ]
    unresolved_type2_from_fourq = [
        row
        for row in rows
        if row["orig_Locus"] == FOURQ
        and row["Locus"] == "HML-2_acro_type2"
        and row["Structure"] in PRESENT
    ]
    summary = {
        "schema": "hml2.biologically-annotated-orf.v3",
        "input_rows": len(input_rows),
        "output_rows": len(rows),
        "exact_cnv_artifact_rows_retained_and_excluded": len(artifact_rows),
        "artifact_ids_not_present_in_input": len(
            artifact_ids - {row["ID_Full"] for row in artifact_rows}
        ),
        "8q24.3b_alias_rows_retained_and_excluded": alias_rows,
        "duplicate_present_labels_retained_and_excluded": duplicate_rows,
        "unresolved_exact_source_conflicts_preserved": unresolved_source_conflicts,
        "analysis_included_rows": sum(
            row["analysis_include"] == "1" for row in rows
        ),
        "analysis_excluded_rows_by_reason": dict(
            sorted(
                Counter(
                    row["analysis_exclusion_reason"]
                    for row in rows
                    if row["analysis_include"] == "0"
                ).items()
            )
        ),
        "resolved_4q35.2_present_rows": len(present_fourq),
        "unresolved_type2_rows_previously_labeled_4q35.2": len(
            unresolved_type2_from_fourq
        ),
        "fourq_rows_reassigned_by_host_sequence": Counter(
            row["Locus"]
            for row in rows
            if row["orig_Locus"] == FOURQ and row["Locus"] != FOURQ
        ),
        "present_rows_by_biological_locus": dict(
            sorted(
                Counter(
                    row["Locus"]
                    for row in rows
                    if row["Structure"] in PRESENT
                    and row["analysis_include"] == "1"
                ).items()
            )
        ),
        "insertion_event_policy": {
            TYPE1_EVENT: ["13p13", "15p13b"],
            TYPE2_EVENT: ["15p13a", "21p13", "22p13", "4q35.2_hg38"],
        },
        "fourq_directional_interpretation": (
            "4q35.2 most likely arose from the 21p13 lineage. The direction is supported "
            "by the near-fixed acrocentric Type II family, polymorphic anchored 4q35.2, "
            "and long gag signatures shared with 21p13 but not 15p13a or 22p13."
        ),
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
