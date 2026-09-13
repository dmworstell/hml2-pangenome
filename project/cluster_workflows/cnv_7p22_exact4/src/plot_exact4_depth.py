#!/usr/bin/env python3
"""Exact, threshold-free coverage-ratio validation and diagnostic PDFs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import re
import sys
from pathlib import Path
from typing import Mapping, Sequence

sys.dont_write_bytecode = True


class ExactFourPlotError(ValueError):
    """An exact ratio summary, depth table, KCON PAF, or PDF failed closure."""


SUMMARY_SCHEMA = "hml2_7p22_exact4_coverage_ratio_summary_1"
FLANK_SCHEMA = "hml2_7p22_exact4_haplotype_flank_record_1"
ELEMENT_SCHEMA = "hml2_7p22_exact4_element_ratio_record_1"
RATIONAL_SCHEMA = "hml2_7p22_exact4_reduced_rational_1"
LOCUS = "HML-2_7p22.1"
PRODUCTS = ("raw", "mapq10")
CALLABILITY = {("callable", "callable"), ("uncallable", "empty_own_flank_set"), ("uncallable", "own_flank_mean_below_5x")}
SUMMARY_FIELDS = ("schema_id", "schema_version", "run_id", "sample_id", "locus_id", "haplotype_flank_records", "element_ratio_records", "read_inferred_copy_number")
FLANK_FIELDS = ("schema_id", "schema_version", "run_id", "sample_id", "locus_id", "haplotype_id", "product", "flank_record_key", "flank_depth_sum", "flank_base_count", "F_exact", "callability", "callability_reason", "read_inferred_copy_number")
ELEMENT_FIELDS = ("schema_id", "schema_version", "run_id", "sample_id", "locus_id", "haplotype_id", "contig", "source_start_0based", "source_end_0based", "product", "element_record_key", "flank_record_key", "source_depth_sum", "source_base_count", "S_exact", "R_exact", "callability", "callability_reason", "read_inferred_copy_number")
RATIONAL_FIELDS = ("schema_id", "schema_version", "numerator", "denominator")
FORBIDDEN_RATIO_CATEGORIES = ("coverage-compatible", "coverage-discordant", "compatibility", "discordance", "ratio_tolerance")


def _strict(value: object, label: str) -> None:
    if type(value) in {str, int, bool} or value is None:
        return
    if type(value) is float:
        raise ExactFourPlotError(f"{label} floats/nonfinite values are forbidden")
    if type(value) is list:
        for index, item in enumerate(value):
            _strict(item, f"{label}[{index}]")
        return
    if type(value) is dict:
        for key, item in value.items():
            if type(key) is not str:
                raise ExactFourPlotError(f"{label} has a non-string key")
            _strict(item, f"{label}.{key}")
        return
    raise ExactFourPlotError(f"{label} contains forbidden type {type(value).__name__}")


def _keys(value: Mapping[str, object], expected: Sequence[str], label: str) -> None:
    if set(value) != set(expected):
        raise ExactFourPlotError(f"{label} keys drifted")


def _integer(value: object, label: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ExactFourPlotError(f"{label} must be an exact integer >= {minimum}")
    return value


def _string(value: object, label: str) -> str:
    if type(value) is not str or not value:
        raise ExactFourPlotError(f"{label} must be an exact nonempty string")
    return value


def _rational(numerator: int, denominator: int) -> dict[str, object]:
    if type(numerator) is not int or type(denominator) is not int or numerator < 0 or denominator <= 0:
        raise ExactFourPlotError("invalid exact rational inputs")
    divisor = math.gcd(numerator, denominator)
    numerator, denominator = numerator // divisor, denominator // divisor
    return {"schema_id": RATIONAL_SCHEMA, "schema_version": 1, "numerator": numerator, "denominator": denominator}


def _validate_rational(value: object, label: str, *, nullable: bool) -> dict[str, object] | None:
    if value is None:
        if nullable:
            return None
        raise ExactFourPlotError(f"{label} cannot be null")
    if type(value) is not dict:
        raise ExactFourPlotError(f"{label} must be a reduced rational object")
    _keys(value, RATIONAL_FIELDS, label)
    if value["schema_id"] != RATIONAL_SCHEMA or value["schema_version"] != 1:
        raise ExactFourPlotError(f"{label} rational schema drift")
    numerator = _integer(value["numerator"], f"{label} numerator")
    denominator = _integer(value["denominator"], f"{label} denominator", 1)
    if math.gcd(numerator, denominator) != 1 or (numerator == 0 and denominator != 1):
        raise ExactFourPlotError(f"{label} is not a canonical reduced rational")
    return value


def _callability_pair(value: Mapping[str, object], label: str) -> tuple[str, str]:
    pair = (_string(value["callability"], f"{label} callability"), _string(value["callability_reason"], f"{label} reason"))
    if pair not in CALLABILITY:
        raise ExactFourPlotError(f"{label} has a forbidden callability pair")
    return pair


def validate_summary(
    summary: dict[str, object],
    *,
    expected_haplotype_order: Sequence[str] | None = None,
) -> dict[str, object]:
    if type(summary) is not dict:
        raise ExactFourPlotError("summary must be an exact object")
    _strict(summary, "summary")
    _keys(summary, SUMMARY_FIELDS, "summary")
    serialized = json.dumps(summary, ensure_ascii=False, sort_keys=True, separators=(",", ":")).lower()
    if any(phrase in serialized for phrase in FORBIDDEN_RATIO_CATEGORIES):
        raise ExactFourPlotError("ratio-derived compatibility/discordance category is forbidden")
    if summary["schema_id"] != SUMMARY_SCHEMA or summary["schema_version"] != 1 or summary["locus_id"] != LOCUS:
        raise ExactFourPlotError("summary schema/version/locus drift")
    run_id, sample = _string(summary["run_id"], "run_id"), _string(summary["sample_id"], "sample_id")
    if summary["read_inferred_copy_number"] != "not_estimated":
        raise ExactFourPlotError("read_inferred_copy_number must remain not_estimated")
    flanks, elements = summary["haplotype_flank_records"], summary["element_ratio_records"]
    if type(flanks) is not list or len(flanks) != 4 or type(elements) is not list or not elements or len(elements) % 2:
        raise ExactFourPlotError("per-sample summary requires four flank and 2*k element records")
    flank_by_key: dict[tuple[object, ...], dict[str, object]] = {}
    hap_order: list[str] = []
    expected_flank_order = []
    for record in flanks:
        if type(record) is not dict:
            raise ExactFourPlotError("flank record must be an exact object")
        _keys(record, FLANK_FIELDS, "flank record")
        if record["schema_id"] != FLANK_SCHEMA or record["schema_version"] != 1 or record["run_id"] != run_id or record["sample_id"] != sample or record["locus_id"] != LOCUS or record["read_inferred_copy_number"] != "not_estimated":
            raise ExactFourPlotError("flank record identity drift")
        hap = _string(record["haplotype_id"], "haplotype_id")
        product = _string(record["product"], "product")
        if product not in PRODUCTS:
            raise ExactFourPlotError("product must be raw or mapq10")
        if hap not in hap_order:
            hap_order.append(hap)
        key = [sample, LOCUS, hap, product]
        if record["flank_record_key"] != key:
            raise ExactFourPlotError("flank_record_key drift")
        tuple_key = tuple(key)
        if tuple_key in flank_by_key:
            raise ExactFourPlotError("duplicate flank record key")
        depth_sum = _integer(record["flank_depth_sum"], "flank depth sum")
        bases = _integer(record["flank_base_count"], "flank base count")
        rational = _validate_rational(record["F_exact"], "F_exact", nullable=True)
        pair = _callability_pair(record, "flank record")
        if bases == 0:
            if depth_sum != 0 or rational is not None or pair != ("uncallable", "empty_own_flank_set"):
                raise ExactFourPlotError("empty own-flank law drift")
        else:
            expected = _rational(depth_sum, bases)
            if rational != expected:
                raise ExactFourPlotError("F_exact does not reduce flank sum/count")
            callable_value = expected["numerator"] >= 5 * expected["denominator"]
            expected_pair = ("callable", "callable") if callable_value else ("uncallable", "own_flank_mean_below_5x")
            if pair != expected_pair:
                raise ExactFourPlotError("flank callability does not use exact 5/1 boundary")
        flank_by_key[tuple_key] = record
        expected_flank_order.append((hap, product))
    if len(hap_order) != 2 or expected_flank_order != [(hap, product) for hap in hap_order for product in PRODUCTS]:
        raise ExactFourPlotError("flank records must follow declared haplotype and raw/mapq10 order")
    if expected_haplotype_order is not None:
        declared = list(expected_haplotype_order)
        if len(declared) != 2 or len(set(declared)) != 2 or hap_order != declared:
            raise ExactFourPlotError("summary haplotype order differs from declared assembly authority")
    for hap in hap_order:
        paired = [flank_by_key[(sample, LOCUS, hap, product)] for product in PRODUCTS]
        if [record["product"] for record in paired] != list(PRODUCTS):
            raise ExactFourPlotError("flank raw/mapq10 product pairing drift")
        if paired[0]["flank_base_count"] != paired[1]["flank_base_count"]:
            raise ExactFourPlotError("flank raw/mapq10 base-count pairing drift")
    seen_elements = set()
    observed_order = []
    element_products: dict[tuple[object, ...], list[str]] = {}
    for record in elements:
        if type(record) is not dict:
            raise ExactFourPlotError("element record must be an exact object")
        _keys(record, ELEMENT_FIELDS, "element record")
        if record["schema_id"] != ELEMENT_SCHEMA or record["schema_version"] != 1 or record["run_id"] != run_id or record["sample_id"] != sample or record["locus_id"] != LOCUS or record["read_inferred_copy_number"] != "not_estimated":
            raise ExactFourPlotError("element record identity drift")
        hap, contig, product = _string(record["haplotype_id"], "element haplotype"), _string(record["contig"], "element contig"), _string(record["product"], "element product")
        start, end = _integer(record["source_start_0based"], "source start"), _integer(record["source_end_0based"], "source end", 1)
        if hap not in hap_order or product not in PRODUCTS or start >= end:
            raise ExactFourPlotError("element identity/product drift")
        key = [sample, LOCUS, hap, contig, start, end, product]
        flank_key = [sample, LOCUS, hap, product]
        if record["element_record_key"] != key or record["flank_record_key"] != flank_key:
            raise ExactFourPlotError("element/flank record key drift")
        if tuple(key) in seen_elements or tuple(flank_key) not in flank_by_key:
            raise ExactFourPlotError("duplicate element or missing unique flank reference")
        seen_elements.add(tuple(key))
        identity = (sample, LOCUS, hap, contig, start, end)
        element_products.setdefault(identity, []).append(product)
        source_sum = _integer(record["source_depth_sum"], "source depth sum")
        source_bases = _integer(record["source_base_count"], "source base count", 1)
        if source_bases != end - start:
            raise ExactFourPlotError("source base count differs from exact interval")
        if _validate_rational(record["S_exact"], "S_exact", nullable=False) != _rational(source_sum, source_bases):
            raise ExactFourPlotError("S_exact does not reduce source sum/count")
        flank = flank_by_key[tuple(flank_key)]
        pair = _callability_pair(record, "element record")
        if pair != (flank["callability"], flank["callability_reason"]):
            raise ExactFourPlotError("element does not inherit unique flank callability")
        ratio = _validate_rational(record["R_exact"], "R_exact", nullable=True)
        if pair == ("callable", "callable"):
            expected_ratio = _rational(source_sum * flank["flank_base_count"], source_bases * flank["flank_depth_sum"])
            if ratio != expected_ratio:
                raise ExactFourPlotError("R_exact is not exact same-haplotype S/F")
        elif ratio is not None:
            raise ExactFourPlotError("uncallable element must have literal-null R_exact")
        observed_order.append((hap_order.index(hap), contig.encode("utf-8"), start, end, PRODUCTS.index(product)))
    if observed_order != sorted(observed_order):
        raise ExactFourPlotError("element records violate canonical haplotype/contig/interval/product order")
    if any(products != list(PRODUCTS) for products in element_products.values()):
        raise ExactFourPlotError("every element must have one exact raw/mapq10 product pair")
    if len(elements) != 2 * len(element_products):
        raise ExactFourPlotError("element product-pair cardinality drift")
    return summary


def _depth_rows(path: Path) -> tuple[list[tuple[str, int, int]], str]:
    if not path.is_file() or path.is_symlink():
        raise ExactFourPlotError(f"depth table is not a regular file: {path}")
    payload = path.read_bytes()
    if not payload or not payload.endswith(b"\n"):
        raise ExactFourPlotError("depth table must be nonempty and LF terminated")
    rows = []
    try:
        for number, raw in enumerate(payload.decode("ascii").splitlines(), 1):
            fields = raw.split("\t")
            if len(fields) != 3 or not fields[0] or not fields[1].isdigit() or not fields[2].isdigit():
                raise ExactFourPlotError(f"malformed depth row {number}")
            position, depth = int(fields[1]), int(fields[2])
            if position < 1:
                raise ExactFourPlotError("depth position must be one-based positive")
            rows.append((fields[0], position, depth))
    except UnicodeError as error:
        raise ExactFourPlotError("depth table is not strict ASCII") from error
    keys = [(contig, position) for contig, position, _ in rows]
    if keys != sorted(keys) or len(keys) != len(set(keys)):
        raise ExactFourPlotError("depth coordinates must be sorted and unique")
    return rows, hashlib.sha256(payload).hexdigest()


def _parse_paf(path: Path) -> tuple[list[tuple[list[str], bytes]], str]:
    if not path.is_file() or path.is_symlink():
        raise ExactFourPlotError("KCON PAF is missing or symlinked")
    payload = path.read_bytes()
    if not payload or not payload.endswith(b"\n"):
        raise ExactFourPlotError("KCON PAF must be nonempty and LF terminated")
    rows = []
    for number, raw in enumerate(payload.splitlines(), 1):
        try:
            fields = raw.decode("ascii").split("\t")
            numeric_fields = (1, 2, 3, 6, 7, 8, 9, 10, 11)
            if len(fields) < 12 or not fields[0] or fields[4] not in {"+", "-"} or not fields[5] or not all(fields[index].isdigit() for index in numeric_fields):
                raise ValueError
            query_length, query_start, query_end = (int(fields[index]) for index in (1, 2, 3))
            target_length, target_start, target_end = (int(fields[index]) for index in (6, 7, 8))
            matches, block_length, mapq = (int(fields[index]) for index in (9, 10, 11))
            if not (
                query_length > 0
                and 0 <= query_start < query_end <= query_length
                and target_length > 0
                and 0 <= target_start < target_end <= target_length
                and 0 <= matches <= block_length
                and block_length > 0
                and 0 <= mapq <= 255
            ):
                raise ValueError
        except (UnicodeError, ValueError) as error:
            raise ExactFourPlotError(f"malformed KCON PAF row {number}") from error
        rows.append((fields, raw))
    return rows, hashlib.sha256(payload).hexdigest()


def _declared_sample_authority(
    plan: Mapping[str, object], sample_index: int
) -> tuple[str, list[str], list[dict[str, object]], list[dict[str, object]]]:
    """Return the exact sample authority in its declared assembly order."""
    if type(sample_index) is not int or not 0 <= sample_index < 4:
        raise ExactFourPlotError("sample index must be 0..3")
    analysis = plan["analysis_map"][sample_index]
    preparation = plan["preparation_map"][sample_index]
    sample = analysis["sample_id"]
    if preparation["sample_id"] != sample or preparation["index"] != sample_index or analysis["index"] != sample_index:
        raise ExactFourPlotError("sample preparation/analysis authority drift")
    resident = plan["immutable_run_material"]["resident_identity"]
    assemblies = [row for row in resident["assemblies"] if row["sample_id"] == sample]
    haplotypes = [row["haplotype"] for row in assemblies]
    if len(haplotypes) != 2 or len(set(haplotypes)) != 2 or preparation["haplotypes"] != haplotypes:
        raise ExactFourPlotError("declared assembly-authority haplotype order drift")
    beds = [row for row in resident["bed_authorities"] if row["sample_id"] == sample]
    elements = [row for row in resident["element_identities"] if row["sample_id"] == sample]
    if not beds or not elements:
        raise ExactFourPlotError("sample lacks authenticated bed or element authority")
    assembly_contigs = {
        (row["haplotype"], contig)
        for row in assemblies
        for contig in row["contigs"]
    }
    if any((row["haplotype"], row["contig"]) not in assembly_contigs for row in beds + elements):
        raise ExactFourPlotError("bed/element authority is absent from the declared assembly")
    expected_element_order = sorted(
        elements,
        key=lambda row: (
            haplotypes.index(row["haplotype"]),
            row["contig"].encode("utf-8"),
            row["source_start"],
            row["source_end"],
        ),
    )
    return sample, haplotypes, beds, expected_element_order


def _select_kcon_annotations_from_plan(
    plan: dict[str, object], sample_index: int, kcon_paf_path: Path
) -> list[dict[str, object]]:
    """Select one exact, annotation-only KCON PAF row per authority element."""
    try:
        import exact4_plan
    except ImportError as error:
        raise ExactFourPlotError("exact4_plan must be beside the frozen plotter") from error
    plan = exact4_plan.validate_target_only_plan(plan)
    sample, haplotypes, beds, elements = _declared_sample_authority(plan, sample_index)
    paf_rows, paf_sha256 = _parse_paf(kcon_paf_path)
    authority = {
        "sample_id": sample,
        "haplotypes": list(haplotypes),
        "bed_authorities": [dict(row) for row in beds],
        "element_identities": [dict(row) for row in elements],
    }
    return _select_kcon_annotations_core(authority, paf_rows, paf_sha256)


def _validate_kcon_selector_authority(value: object) -> dict[str, object]:
    """Validate the plan-free authority accepted by the pure KCON selector."""
    if type(value) is not dict or set(value) != {
        "sample_id", "haplotypes", "bed_authorities", "element_identities"
    }:
        raise ExactFourPlotError("KCON selector authority has an inexact schema")
    sample = value["sample_id"]
    haplotypes = value["haplotypes"]
    beds = value["bed_authorities"]
    elements = value["element_identities"]
    if type(sample) is not str or not sample or type(haplotypes) is not list or len(haplotypes) != 2 or any(type(item) is not str or not item for item in haplotypes) or len(set(haplotypes)) != 2:
        raise ExactFourPlotError("KCON selector sample/haplotype authority is invalid")
    if type(beds) is not list or not beds or type(elements) is not list or not elements:
        raise ExactFourPlotError("KCON selector authority lacks beds or elements")
    normalized_beds: list[dict[str, object]] = []
    for row in beds:
        if type(row) is not dict or set(row) != {"haplotype", "contig", "full_intervals"}:
            raise ExactFourPlotError("KCON selector bed authority has an inexact schema")
        hap, contig, intervals = row["haplotype"], row["contig"], row["full_intervals"]
        if hap not in haplotypes or type(contig) is not str or not contig or type(intervals) is not list or not intervals:
            raise ExactFourPlotError("KCON selector bed authority is invalid")
        checked: list[list[int]] = []
        for interval in intervals:
            if type(interval) not in (list, tuple) or len(interval) != 2 or any(type(x) is not int for x in interval) or not 0 <= interval[0] < interval[1]:
                raise ExactFourPlotError("KCON selector full interval is invalid")
            checked.append([interval[0], interval[1]])
        if checked != sorted(checked) or any(left[1] > right[0] for left, right in zip(checked, checked[1:])):
            raise ExactFourPlotError("KCON selector full intervals are not ordered/disjoint")
        normalized_beds.append({"haplotype": hap, "contig": contig, "full_intervals": checked})
    normalized_elements: list[dict[str, object]] = []
    for row in elements:
        if type(row) is not dict or set(row) != {"haplotype", "contig", "source_start", "source_end"}:
            raise ExactFourPlotError("KCON selector element authority has an inexact schema")
        hap, contig, start, end = row["haplotype"], row["contig"], row["source_start"], row["source_end"]
        if hap not in haplotypes or type(contig) is not str or not contig or type(start) is not int or type(end) is not int or not 0 <= start < end:
            raise ExactFourPlotError("KCON selector element authority is invalid")
        containing = [
            (bed["contig"], interval[0], interval[1])
            for bed in normalized_beds
            if bed["haplotype"] == hap and bed["contig"] == contig
            for interval in bed["full_intervals"]
            if interval[0] <= start and end <= interval[1]
        ]
        if len(containing) != 1:
            raise ExactFourPlotError("each authority element must have one exact containing full window")
        normalized_elements.append({"haplotype": hap, "contig": contig, "source_start": start, "source_end": end})
    expected_order = sorted(normalized_elements, key=lambda row: (haplotypes.index(row["haplotype"]), row["contig"].encode("utf-8"), row["source_start"], row["source_end"]))
    if normalized_elements != expected_order or len({(row["haplotype"], row["contig"], row["source_start"], row["source_end"]) for row in normalized_elements}) != len(normalized_elements):
        raise ExactFourPlotError("KCON selector elements are not canonical and unique")
    return {"sample_id": sample, "haplotypes": list(haplotypes), "bed_authorities": normalized_beds, "element_identities": normalized_elements}


def _select_kcon_annotations_core(
    authority: dict[str, object],
    paf_rows: list[tuple[list[str], bytes]],
    paf_sha256: str,
) -> list[dict[str, object]]:
    """Pure selector over already-validated authority and already-parsed PAF."""
    sample = authority["sample_id"]
    haplotypes = authority["haplotypes"]
    beds = authority["bed_authorities"]
    elements = authority["element_identities"]
    annotations: list[dict[str, object]] = []
    for element in elements:
        containing = [
            (row["contig"], start, end)
            for row in beds
            if row["haplotype"] == element["haplotype"] and row["contig"] == element["contig"]
            for start, end in row["full_intervals"]
            if start <= element["source_start"] and element["source_end"] <= end
        ]
        if len(containing) != 1:
            raise ExactFourPlotError("each authority element must have one exact containing full window")
        contig, window_start, window_end = containing[0]
        expected_target = f"{contig}:{window_start + 1}-{window_end}"
        expected_target_length = window_end - window_start
        qualifying: list[tuple[list[str], bytes]] = []
        for fields, raw in paf_rows:
            if fields[0] != "type2_KCON" or int(fields[1]) != 9472:
                continue
            if fields[5] != expected_target or int(fields[6]) != expected_target_length:
                continue
            absolute_target_start = window_start + int(fields[7])
            absolute_target_end = window_start + int(fields[8])
            if int(fields[10]) >= 500 and absolute_target_start < element["source_end"] and absolute_target_end > element["source_start"]:
                qualifying.append((fields, raw))
        if not qualifying:
            raise ExactFourPlotError("authenticated source window lacks a qualifying Type-II KCON row")
        fields, raw = min(
            qualifying,
            key=lambda item: (
                -int(item[0][10]),
                -int(item[0][9]),
                -int(item[0][11]),
                int(item[0][7]),
                int(item[0][2]),
                item[1],
            ),
        )
        annotations.append({
            "sample_id": sample,
            "haplotype_id": element["haplotype"],
            "contig": contig,
            "source_start_0based": element["source_start"],
            "source_end_0based": element["source_end"],
            "query_name": fields[0],
            "query_length": int(fields[1]),
            "query_start_0based": int(fields[2]),
            "query_end_0based": int(fields[3]),
            "strand": fields[4],
            "target_name": fields[5],
            "target_length": int(fields[6]),
            "target_start_0based": int(fields[7]),
            "target_end_0based": int(fields[8]),
            "matches": int(fields[9]),
            "aligned_block_length": int(fields[10]),
            "mapq": int(fields[11]),
            "paf_row_sha256": hashlib.sha256(raw).hexdigest(),
            "paf_file_sha256": paf_sha256,
            "annotation_only": True,
            "read_inferred_copy_number": "not_estimated",
        })
    observed = [
        (row["haplotype_id"], row["contig"], row["source_start_0based"], row["source_end_0based"])
        for row in annotations
    ]
    expected = [
        (row["haplotype"], row["contig"], row["source_start"], row["source_end"])
        for row in elements
    ]
    if observed != expected or any(row["haplotype_id"] not in haplotypes for row in annotations):
        raise ExactFourPlotError("KCON annotation authority order drift")
    return annotations


def _build_summary_and_kcon_annotations_from_plan(
    plan: dict[str, object],
    sample_index: int,
    raw_depth_path: Path,
    mapq10_depth_path: Path,
    kcon_paf_path: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Build the exact summary and its plan-bound, annotation-only KCON rows."""
    try:
        import exact4_plan
    except ImportError as error:
        raise ExactFourPlotError("exact4_plan must be beside the frozen plotter") from error
    plan = exact4_plan.validate_target_only_plan(plan)
    sample, hap_order, beds, elements = _declared_sample_authority(plan, sample_index)
    raw_rows, _ = _depth_rows(raw_depth_path)
    mapq_rows, _ = _depth_rows(mapq10_depth_path)
    if [(c, p) for c, p, _ in raw_rows] != [(c, p) for c, p, _ in mapq_rows]:
        raise ExactFourPlotError("raw and MAPQ10 coordinate inventories differ")
    raw_map = {(c, p - 1): d for c, p, d in raw_rows}
    mapq_map = {(c, p - 1): d for c, p, d in mapq_rows}
    if any(raw_map[key] < mapq_map[key] for key in raw_map):
        raise ExactFourPlotError("raw depth is below MAPQ10 depth")
    expected = [(row["contig"], pos) for row in beds for start, end in row["full_intervals"] for pos in range(start, end)]
    if sorted(raw_map) != sorted(expected) or len(raw_map) != len(expected):
        raise ExactFourPlotError("depth inventory does not equal full-window union")
    kcon_annotations = _select_kcon_annotations_from_plan(plan, sample_index, kcon_paf_path)
    flank_records, element_records = [], []
    for hap in hap_order:
        hap_beds = [row for row in beds if row["haplotype"] == hap]
        flank_coords = [(row["contig"], pos) for row in hap_beds for start, end in row["flank_intervals"] for pos in range(start, end)]
        for product, depths in (("raw", raw_map), ("mapq10", mapq_map)):
            depth_sum, base_count = sum(depths[key] for key in flank_coords), len(flank_coords)
            if base_count == 0:
                f_exact, callability, reason = None, "uncallable", "empty_own_flank_set"
            else:
                f_exact = _rational(depth_sum, base_count)
                if f_exact["numerator"] >= 5 * f_exact["denominator"]:
                    callability, reason = "callable", "callable"
                else:
                    callability, reason = "uncallable", "own_flank_mean_below_5x"
            flank_key = [sample, LOCUS, hap, product]
            flank_records.append({"schema_id": FLANK_SCHEMA, "schema_version": 1, "run_id": plan["run_id"], "sample_id": sample, "locus_id": LOCUS, "haplotype_id": hap, "product": product, "flank_record_key": flank_key, "flank_depth_sum": depth_sum, "flank_base_count": base_count, "F_exact": f_exact, "callability": callability, "callability_reason": reason, "read_inferred_copy_number": "not_estimated"})
    flank_lookup = {tuple(record["flank_record_key"]): record for record in flank_records}
    for element in elements:
        source_coords = [(element["contig"], pos) for pos in range(element["source_start"], element["source_end"])]
        for product, depths in (("raw", raw_map), ("mapq10", mapq_map)):
            source_sum, source_count = sum(depths[key] for key in source_coords), len(source_coords)
            s_exact = _rational(source_sum, source_count)
            flank_key = [sample, LOCUS, element["haplotype"], product]
            flank = flank_lookup[tuple(flank_key)]
            r_exact = _rational(source_sum * flank["flank_base_count"], source_count * flank["flank_depth_sum"]) if flank["callability"] == "callable" else None
            element_key = [sample, LOCUS, element["haplotype"], element["contig"], element["source_start"], element["source_end"], product]
            element_records.append({"schema_id": ELEMENT_SCHEMA, "schema_version": 1, "run_id": plan["run_id"], "sample_id": sample, "locus_id": LOCUS, "haplotype_id": element["haplotype"], "contig": element["contig"], "source_start_0based": element["source_start"], "source_end_0based": element["source_end"], "product": product, "element_record_key": element_key, "flank_record_key": flank_key, "source_depth_sum": source_sum, "source_base_count": source_count, "S_exact": s_exact, "R_exact": r_exact, "callability": flank["callability"], "callability_reason": flank["callability_reason"], "read_inferred_copy_number": "not_estimated"})
    summary = {"schema_id": SUMMARY_SCHEMA, "schema_version": 1, "run_id": plan["run_id"], "sample_id": sample, "locus_id": LOCUS, "haplotype_flank_records": flank_records, "element_ratio_records": element_records, "read_inferred_copy_number": "not_estimated"}
    return validate_summary(summary, expected_haplotype_order=hap_order), kcon_annotations


def _build_summary_from_plan(plan: dict[str, object], sample_index: int, raw_depth_path: Path, mapq10_depth_path: Path, kcon_paf_path: Path) -> dict[str, object]:
    """Build the exact per-sample ratio summary in a private analysis stage."""
    summary, _ = _build_summary_and_kcon_annotations_from_plan(
        plan, sample_index, raw_depth_path, mapq10_depth_path, kcon_paf_path
    )
    return summary


def _pdf_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _minimal_pdf(lines: Sequence[str]) -> bytes:
    commands = ["BT", "/F1 8 Tf", "36 756 Td"]
    for index, line in enumerate(lines):
        if index:
            commands.append("0 -11 Td")
        commands.append(f"({_pdf_escape(line)}) Tj")
    commands.append("ET")
    stream = ("\n".join(commands) + "\n").encode("latin-1", "replace")
    objects = [b"<< /Type /Catalog /Pages 2 0 R >>", b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>", b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 5 0 R >> >> /Contents 4 0 R >>", b"<< /Length " + str(len(stream)).encode() + b" >>\nstream\n" + stream + b"endstream", b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>"]
    output = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n"); offsets = [0]
    for number, obj in enumerate(objects, 1):
        offsets.append(len(output)); output.extend(f"{number} 0 obj\n".encode()); output.extend(obj + b"\nendobj\n")
    xref = len(output); output.extend(f"xref\n0 {len(objects)+1}\n".encode()); output.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]: output.extend(f"{offset:010d} 00000 n \n".encode())
    output.extend(f"trailer\n<< /Size {len(objects)+1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n".encode()); return bytes(output)


def _rational_text(value: object) -> str:
    return "null" if value is None else f"{value['numerator']}/{value['denominator']}"


def _validate_kcon_annotations(
    summary: Mapping[str, object], annotations: Sequence[Mapping[str, object]]
) -> list[Mapping[str, object]]:
    if not annotations:
        return []
    expected_identities = []
    for record in summary["element_ratio_records"]:
        if record["product"] == PRODUCTS[0]:
            expected_identities.append((
                record["sample_id"],
                record["haplotype_id"],
                record["contig"],
                record["source_start_0based"],
                record["source_end_0based"],
            ))
    observed_identities = []
    for annotation in annotations:
        required = {
            "sample_id", "haplotype_id", "contig", "source_start_0based", "source_end_0based",
            "query_name", "query_length", "query_start_0based", "query_end_0based", "strand",
            "target_name", "target_length", "target_start_0based", "target_end_0based",
            "matches", "aligned_block_length", "mapq", "paf_row_sha256", "paf_file_sha256",
            "annotation_only", "read_inferred_copy_number",
        }
        if type(annotation) is not dict or set(annotation) != required:
            raise ExactFourPlotError("KCON annotation fields drifted")
        if annotation["query_name"] != "type2_KCON" or annotation["query_length"] != 9472:
            raise ExactFourPlotError("KCON annotation query identity drifted")
        if annotation["aligned_block_length"] < 500 or annotation["annotation_only"] is not True:
            raise ExactFourPlotError("KCON annotation scientific law drifted")
        if annotation["read_inferred_copy_number"] != "not_estimated":
            raise ExactFourPlotError("KCON annotation cannot estimate read copy number")
        for key in ("paf_row_sha256", "paf_file_sha256"):
            if type(annotation[key]) is not str or re.fullmatch(r"[0-9a-f]{64}", annotation[key]) is None:
                raise ExactFourPlotError(f"KCON annotation {key} is not an exact SHA-256")
        observed_identities.append((
            annotation["sample_id"],
            annotation["haplotype_id"],
            annotation["contig"],
            annotation["source_start_0based"],
            annotation["source_end_0based"],
        ))
    if observed_identities != expected_identities:
        raise ExactFourPlotError("KCON annotations do not pair exactly with summary elements")
    if annotations and len({row["paf_file_sha256"] for row in annotations}) != 1:
        raise ExactFourPlotError("KCON annotations do not share one exact PAF identity")
    return list(annotations)


def _labels(
    summary: Mapping[str, object],
    metric: str,
    kcon_annotations: Sequence[Mapping[str, object]] = (),
) -> list[str]:
    labels = ["HML-2 exact-four threshold-free diagnostic ratio", f"run={summary['run_id']}", f"sample={summary['sample_id']}", f"locus={summary['locus_id']}", f"product={metric}", "read_inferred_copy_number=not_estimated", "competitive chosen-primary; not phased truth", "ratio categories/thresholds absent"]
    for record in summary["haplotype_flank_records"]:
        if record["product"] == metric:
            labels.append(f"flank={record['haplotype_id']} sum={record['flank_depth_sum']} bases={record['flank_base_count']} F_exact={_rational_text(record['F_exact'])} callability={record['callability']} reason={record['callability_reason']}")
    for record in summary["element_ratio_records"]:
        if record["product"] == metric:
            labels.append(f"element={record['sample_id']}|{record['haplotype_id']}|{record['contig']}:{record['source_start_0based']}-{record['source_end_0based']} S_exact={_rational_text(record['S_exact'])} R_exact={_rational_text(record['R_exact'])} callability={record['callability']} reason={record['callability_reason']}")
    for annotation in _validate_kcon_annotations(summary, kcon_annotations):
        identity = f"{annotation['sample_id']}|{annotation['haplotype_id']}|{annotation['contig']}:{annotation['source_start_0based']}-{annotation['source_end_0based']}"
        labels.extend([
            f"KCON selected annotation={identity} annotation_only=true",
            f"KCON query={annotation['query_name']}:{annotation['query_start_0based']}-{annotation['query_end_0based']}/{annotation['query_length']} strand={annotation['strand']}",
            f"KCON target={annotation['target_name']}:{annotation['target_start_0based']}-{annotation['target_end_0based']}/{annotation['target_length']}",
            f"KCON evidence matches={annotation['matches']} block={annotation['aligned_block_length']} MAPQ={annotation['mapq']}",
            f"KCON selected_row_sha256={annotation['paf_row_sha256']}",
            f"KCON paf_file_sha256={annotation['paf_file_sha256']}",
        ])
    return labels


def _load_target_only_plan(path: Path) -> dict[str, object]:
    try:
        import exact4_plan
    except ImportError as error:
        raise ExactFourPlotError("exact4_plan must be beside the frozen plotter") from error
    if not path.is_file() or path.is_symlink():
        raise ExactFourPlotError("target-only plan must be one regular file")
    raw = path.read_bytes()
    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=lambda pairs: _reject_duplicate_pairs(pairs, "target-only plan"),
            parse_constant=lambda token: (_ for _ in ()).throw(ExactFourPlotError(f"nonfinite token {token}")),
        )
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ExactFourPlotError(f"cannot load strict target-only plan: {error}") from error
    if exact4_plan.canonical_json_bytes(value) != raw:
        raise ExactFourPlotError("target-only plan must be canonical JSON")
    return exact4_plan.validate_target_only_plan(value)


def _reject_duplicate_pairs(pairs: list[tuple[str, object]], label: str) -> dict[str, object]:
    output: dict[str, object] = {}
    for key, value in pairs:
        if key in output:
            raise ExactFourPlotError(f"duplicate {label} key: {key}")
        output[key] = value
    return output


def render_plot(
    depth_path: pathlib.Path,
    summary_path: pathlib.Path,
    output_path: pathlib.Path,
    metric: str,
    *,
    target_only_plan_path: pathlib.Path | None = None,
    sample_index: int | None = None,
    kcon_paf_path: pathlib.Path | None = None,
) -> None:
    if metric not in PRODUCTS:
        raise ExactFourPlotError("metric must be exactly raw or mapq10")
    try:
        raw = summary_path.read_bytes(); summary = json.loads(raw.decode("utf-8"), object_pairs_hook=lambda pairs: _reject_duplicate_pairs(pairs, "summary"), parse_constant=lambda token: (_ for _ in ()).throw(ExactFourPlotError(f"nonfinite token {token}")))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ExactFourPlotError(f"cannot load strict summary: {error}") from error
    if (json.dumps(summary, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode() != raw:
        raise ExactFourPlotError("summary must be canonical JSON")
    authority_arguments = (target_only_plan_path, sample_index, kcon_paf_path)
    if any(value is not None for value in authority_arguments) and not all(value is not None for value in authority_arguments):
        raise ExactFourPlotError("plan, sample index, and KCON PAF must be supplied together")
    kcon_annotations: list[dict[str, object]] = []
    if all(value is not None for value in authority_arguments):
        plan = _load_target_only_plan(target_only_plan_path)
        sample, haplotypes, _, _ = _declared_sample_authority(plan, sample_index)
        if summary["run_id"] != plan["run_id"] or summary["sample_id"] != sample:
            raise ExactFourPlotError("summary differs from target-only plan run/sample authority")
        validate_summary(summary, expected_haplotype_order=haplotypes)
        kcon_annotations = _select_kcon_annotations_from_plan(plan, sample_index, kcon_paf_path)
    else:
        validate_summary(summary)
    rows, _ = _depth_rows(depth_path)
    if not rows:
        raise ExactFourPlotError("depth table has no coordinate inventory")
    payload = _minimal_pdf(_labels(summary, metric, kcon_annotations)); output_path.parent.mkdir(parents=False, exist_ok=True)
    try:
        with output_path.open("xb") as handle: handle.write(payload); handle.flush()
    except FileExistsError as error:
        raise ExactFourPlotError("refusing to replace existing PDF") from error
    expected = {"summary": summary, "metric": metric}
    if kcon_annotations:
        expected["kcon_annotations"] = kcon_annotations
    validate_pdf(output_path, expected)


def validate_pdf(path: pathlib.Path, expected: dict[str, object]) -> None:
    if type(expected) is not dict or set(expected) not in ({"summary", "metric"}, {"summary", "metric", "kcon_annotations"}):
        raise ExactFourPlotError("PDF authority must contain summary/metric and optional bound KCON annotations")
    summary = validate_summary(expected["summary"]); metric = expected["metric"]
    kcon_annotations = expected.get("kcon_annotations", [])
    if metric not in PRODUCTS or not path.is_file() or path.is_symlink():
        raise ExactFourPlotError("PDF metric/path closure failed")
    payload = path.read_bytes()
    if not payload.startswith(b"%PDF-") or not payload.rstrip().endswith(b"%%EOF") or len(re.findall(rb"/Type\s*/Page\b", payload)) != 1:
        raise ExactFourPlotError("PDF structural closure failed")
    for label in _labels(summary, metric, kcon_annotations):
        if _pdf_escape(label).encode("latin-1", "replace") not in payload:
            raise ExactFourPlotError(f"PDF lacks exact diagnostic label: {label}")


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--depth", type=Path, required=True); parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True); parser.add_argument("--metric", choices=PRODUCTS, required=True)
    parser.add_argument("--target-only-plan", type=Path, required=True)
    parser.add_argument("--sample-index", type=int, choices=range(4), required=True)
    parser.add_argument("--kcon-paf", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        render_plot(
            args.depth,
            args.summary,
            args.output,
            args.metric,
            target_only_plan_path=args.target_only_plan,
            sample_index=args.sample_index,
            kcon_paf_path=args.kcon_paf,
        )
        return 0
    except ExactFourPlotError as error: print(f"plot_exact4_depth: {error}", file=sys.stderr); return 2


if __name__ == "__main__":
    raise SystemExit(main())
