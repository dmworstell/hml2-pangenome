#!/usr/bin/env python3
"""Finalize normalized CNV measurement weights for materialized assays."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics
from typing import Any

ARTIFACT = "assembly_artifact"
LATER = "later_duplication"
SEGDUP = "authenticated_segdup"
ARRAY = "array_multipart"
STATES = (ARTIFACT, LATER, SEGDUP, ARRAY)
LATENT = "NON_SEGDUP_DUPLICATION_LATENT_CNV_WEIGHTED"
FIXED_SEGDUP = "AUTHENTICATED_SEGMENTAL_DUPLICATION_FIXED"
UNMEASURED = "UNINFORMATIVE_NO_LOCAL_FLANK_MEASUREMENT"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        rows = list(reader)
    if not reader.fieldnames or any(None in row for row in rows):
        raise SystemExit(f"malformed TSV: {path}")
    return rows


def number(value: Any, label: str) -> float:
    if value in (None, "", "NA", "None", "nan"):
        raise ValueError(f"missing {label}")
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"invalid {label}: {value}")
    return result


def normalized_depth_law(ratio: float, sigma: float) -> dict[str, float]:
    if not 0 <= ratio <= 1 or sigma <= 0:
        raise ValueError("invalid depth-likelihood input")
    logs = {
        ARTIFACT: -0.5 * (ratio / sigma) ** 2,
        LATER: -0.5 * ((ratio - 1.0) / sigma) ** 2,
    }
    anchor = max(logs.values())
    raw = {state: math.exp(value - anchor) for state, value in logs.items()}
    total = sum(raw.values())
    return {
        ARTIFACT: raw[ARTIFACT] / total,
        LATER: raw[LATER] / total,
        SEGDUP: 0.0,
        ARRAY: 0.0,
    }


def boundary_measurements(
    rows: list[dict[str, str]], targets: set[str]
) -> dict[str, dict[str, Any]]:
    result = {}
    accepted = {
        "INFORMATIVE_RAW_ONE_OR_TWO_SIDED_DEPTH",
        "INFORMATIVE_EXACT_SIDE_DEPTH",
    }
    for row in rows:
        key = row["candidate_id"]
        if key not in targets or row["current_resolution"] not in accepted:
            continue
        result[key] = {
            "body": number(row["raw_body_median"], "raw body median"),
            "flank": number(row["raw_local_flank_median"], "raw local flank"),
            "source": (
                "AUTHENTICATED_NUMERIC_OVERRIDE"
                if row["raw_numeric_override"] == "true"
                else "AUTHORITY_BOUND_EXACT_SIDE_DEPTH"
            ),
        }
    return result


def exact_measurements(
    path: Path | None, targets: set[str]
) -> dict[str, dict[str, Any]]:
    if path is None:
        return {}
    result = {}
    for row in read_tsv(path):
        key = row["candidate_id"]
        if key not in targets or key in result:
            raise SystemExit(f"unexpected exact candidate: {key}")
        if row["flank_measurement_status"] != "INFORMATIVE_LOCAL_FLANK_DEPTH":
            continue
        result[key] = {
            "body": number(row["raw_body_median"], "raw body median"),
            "flank": number(row["raw_local_flank_median"], "raw local flank"),
            "source": row["measurement_source"],
        }
    return result


def key_value_receipt(path: Path) -> dict[str, str]:
    result = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t", strict=True):
            if len(row) != 2 or row[0] in result:
                raise SystemExit(f"malformed key-value receipt: {path}")
            result[row[0]] = row[1]
    return result


def exact_depth_measurement(path: Path, authority: dict[str, str]) -> tuple[float, float]:
    selected = []
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if reader.fieldnames != [
            "region_id", "contig", "position_1based", "raw_depth", "mapq10_depth"
        ]:
            raise SystemExit(f"unexpected exact-depth schema: {path}")
        for row in reader:
            if row["region_id"] != authority["region_id"]:
                continue
            if row["contig"] != authority["raw_contig"]:
                raise SystemExit(f"exact-depth contig mismatch: {authority['candidate_id']}")
            selected.append((int(row["position_1based"]), int(row["raw_depth"])))
    if not selected or len(selected) != len({position for position, _ in selected}):
        raise SystemExit(f"missing or duplicate exact depth: {authority['candidate_id']}")
    start = int(authority["raw_body_start0"])
    end = min(int(authority["raw_body_end0"]), int(authority["raw_contig_length"]))
    flank_bp = int(authority["flank_requested_bp"])
    left_min = max(1, start - flank_bp + 1)
    right_max = min(int(authority["raw_contig_length"]), int(authority["raw_body_end0"]) + flank_bp)
    left = [depth for position, depth in selected if left_min <= position <= start]
    body = [depth for position, depth in selected if start < position <= end]
    right = [depth for position, depth in selected if int(authority["raw_body_end0"]) < position <= right_max]
    if len(body) != end - start or not left + right:
        raise SystemExit(f"incomplete exact body/flank depth: {authority['candidate_id']}")
    return float(statistics.median(body)), float(statistics.median(left + right))


def ready_remote_measurements(
    manifest_path: Path | None,
    authority_rows: list[dict[str, str]],
    targets: set[str],
) -> dict[str, dict[str, Any]]:
    if manifest_path is None:
        return {}
    manifest = read_tsv(manifest_path)
    for row in manifest:
        path = manifest_path.parent / row["path"]
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise SystemExit(f"ready-remote manifest mismatch: {path}")
    authority_by_sample_locus: dict[
        tuple[str, str], list[dict[str, str]]
    ] = {}
    for row in authority_rows:
        if row["candidate_id"] in targets:
            authority_by_sample_locus.setdefault(
                (row["sample_id"], row["locus"]), []
            ).append(row)
    result = {}
    summaries = [
        manifest_path.parent / row["path"]
        for row in manifest if row["path"].endswith("/coverage_summary.json")
    ]
    for summary_path in summaries:
        summary = json.loads(summary_path.read_text())
        identity = (summary["sample_id"], summary["locus"])
        authorities = authority_by_sample_locus.get(identity)
        if not authorities:
            raise SystemExit(f"ready-remote identity not authoritative: {summary_path}")
        regions = {
            (region["region_id"], region["contig"]): region
            for region in summary.get("regions", [])
        }
        if len(regions) != len(summary.get("regions", [])):
            raise SystemExit(f"ready-remote regions repeat: {summary_path}")
        locus_root = summary_path.parent
        depth_path = locus_root / "depth.raw_and_mapq10.tsv.gz"
        receipt_path = locus_root / "receipt.tsv"
        receipt = key_value_receipt(receipt_path)
        if (
            receipt.get("status") not in {"COMPLETE", "NO_CALL_ZERO_BODY_DEPTH"}
            or receipt.get("sample_id") != identity[0]
            or receipt.get("locus") != identity[1]
            or receipt.get("depth_sha256") != sha256(depth_path)
            or receipt.get("summary_sha256") != sha256(summary_path)
        ):
            raise SystemExit(f"ready-remote receipt mismatch: {summary_path}")
        for authority in authorities:
            identity_key = (authority["region_id"], authority["raw_contig"])
            if identity_key not in regions:
                raise SystemExit(
                    f"ready-remote region mismatch: {authority['candidate_id']}"
                )
            candidate_id = authority["candidate_id"]
            if candidate_id in result:
                raise SystemExit(f"duplicate ready-remote candidate: {candidate_id}")
            body, flank = exact_depth_measurement(depth_path, authority)
            result[candidate_id] = {
                "body": body,
                "flank": flank,
                "source": "AUTHORITY_BOUND_TRANSFERRED_PER_BASE_DEPTH",
                "evidence_sha256s": [
                    sha256(summary_path),
                    sha256(depth_path),
                    sha256(receipt_path),
                ],
            }
    return result


def site_index(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    sites = payload.get("sites")
    if not isinstance(sites, list) or payload.get("site_count") != len(sites):
        raise SystemExit("site accounting count does not match its sites")
    result = {}
    placement_count = 0
    for site in sites:
        copies = site.get("copy_accounting")
        if not isinstance(copies, list) or not copies:
            raise SystemExit("site accounting contains an empty placement set")
        placement_count += len(copies)
        for copy in copies:
            if copy["candidate_key"] in result:
                raise SystemExit("duplicate site candidate")
            result[copy["candidate_key"]] = {"site": site, "copy": copy}
    if (
        payload.get("placement_count") != placement_count
        or len(result) != placement_count
    ):
        raise SystemExit("site placement count does not match its copy accounting")
    return result


def site_metrics(key: str, body: float, flank: float,
                 entry: dict[str, Any] | None,
                 measurement_overrides: dict[str, dict[str, Any]] | None = None,
                 ) -> dict[str, Any]:
    empty = {
        "site_measurement_completeness": "PARTIAL_CONTEXT_ONLY",
        "site_summed_body_medians": None,
        "site_summed_flank_medians": None,
        "candidate_body_fraction_of_site": None,
        "candidate_flank_fraction_of_site": None,
    }
    if entry is None:
        return {**empty, "site_measurement_completeness": "NOT_IN_BASE_SITE_ACCOUNTING"}
    bodies, flanks = [], []
    overrides = measurement_overrides or {}
    for copy in entry["site"]["copy_accounting"]:
        copy_key = copy["candidate_key"]
        replacement = overrides.get(copy_key)
        b = (body if copy_key == key else
             replacement["body"] if replacement is not None else
             copy["body_median_all"])
        f = (flank if copy_key == key else
             replacement["flank"] if replacement is not None else
             copy["flank_median_all"])
        if b is None or f is None:
            return empty
        bodies.append(float(b))
        flanks.append(float(f))
    total_body, total_flank = sum(bodies), sum(flanks)
    return {
        "site_measurement_completeness": "COMPLETE",
        "site_summed_body_medians": total_body,
        "site_summed_flank_medians": total_flank,
        "candidate_body_fraction_of_site": body / total_body if total_body else None,
        "candidate_flank_fraction_of_site": flank / total_flank if total_flank else None,
    }


def final_row(authority: dict[str, Any], scales: dict[str, float],
              sites: dict[str, dict[str, Any]],
              measurement: dict[str, Any] | None = None,
              all_measurements: dict[str, dict[str, Any]] | None = None,
              ) -> dict[str, Any]:
    key = authority["candidate_key"]
    cls = authority["authority_class"]
    if cls == FIXED_SEGDUP:
        central = {ARTIFACT: 0.0, LATER: 0.0, SEGDUP: 1.0, ARRAY: 0.0}
        laws = None
        envelope = {state: {"minimum": value, "maximum": value}
                    for state, value in central.items()}
        body = float(authority.get("retained_body_median_all") or 0)
        flank = float(authority.get("retained_local_flank_median_all") or 0)
        baseline = authority.get("retained_sample_ref_cov")
        ratio = None if baseline in (None, 0) else flank / float(baseline)
        source = "AUTHENTICATED_STRUCTURAL_SEGDUP_AUTHORITY"
    elif cls == LATENT or (cls == UNMEASURED and measurement is not None):
        if measurement is None:
            body = number(authority["retained_body_median_all"], "body")
            flank = number(authority["retained_local_flank_median_all"], "flank")
            source = "RETAINED_AUTHORITY_LOCAL_FLANK_DEPTH"
        else:
            body, flank = measurement["body"], measurement["flank"]
            source = measurement["source"]
        baseline = number(authority["retained_sample_ref_cov"], "sample baseline")
        ratio = min(1.0, flank / baseline)
        laws = {name: normalized_depth_law(ratio, sigma)
                for name, sigma in scales.items()}
        central = laws["central_rmse_lower_residual"]
        envelope = {state: {
            "minimum": min(law[state] for law in laws.values()),
            "maximum": max(law[state] for law in laws.values()),
        } for state in STATES}
    else:
        raise ValueError(f"unmaterialized candidate: {key}")
    row = {
        "schema": "hml2.cnv-final-assayed-candidate-weight.v1",
        "candidate_key": key,
        "sample": authority["sample"],
        "locus": authority["locus"],
        "region": authority["region"],
        "authority_class": (
            LATENT if cls == UNMEASURED and measurement is not None else cls
        ),
        "measurement_source": source,
        "measurement_evidence_sha256s": (
            [] if measurement is None else measurement.get("evidence_sha256s", [])
        ),
        "raw_body_median": body,
        "raw_local_flank_median": flank,
        "sample_global_diploid_target_flank_baseline": baseline,
        "local_flank_ratio_to_sample_global_baseline": ratio,
        **site_metrics(
            key, body, flank, sites.get(key), all_measurements),
        "normalized_measurement_likelihoods": central,
        "control_scale_state_likelihoods": laws,
        "control_scale_uncertainty_envelope": envelope,
        "dominant_state": max(central, key=central.get),
        "hard_call": False,
        "weight_status": "FINAL_FOR_MATERIALIZED_ASSAY",
        "weight_semantics": "normalized_measurement_likelihood_not_population_posterior",
        "probability_owner": "cnv_depth_weight_finalization_v1_exactly_once",
        "array_requires_direct_shared_ltr_or_junction_authority": True,
        "segdup_requires_authenticated_structural_and_paralog_authority": True,
        "assembly_error_is_not_locus_absence": True,
    }
    row["row_receipt_sha256"] = canonical_sha256(row)
    return row


def write_outputs(outdir: Path, rows: list[dict[str, Any]], pending: list[dict[str, str]],
                  inputs: dict[str, str], ledger: dict[str, Any],
                  marker_receipt: dict[str, Any],
                  expected_candidate_count: int) -> None:
    outdir.mkdir(parents=True, exist_ok=False)
    json_path = outdir / "final_assayed_candidate_weights.v1.json"
    json_path.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    tsv_path = outdir / "final_assayed_candidate_weights.v1.tsv"
    fields = [
        "candidate_key", "sample", "locus", "region", "authority_class",
        "measurement_source", "raw_body_median", "raw_local_flank_median",
        "sample_global_diploid_target_flank_baseline",
        "local_flank_ratio_to_sample_global_baseline",
        "site_measurement_completeness", "site_summed_body_medians",
        "site_summed_flank_medians", "candidate_body_fraction_of_site",
        "candidate_flank_fraction_of_site", *STATES,
        "artifact_min", "artifact_max", "later_min", "later_max",
        "segdup_min", "segdup_max", "array_min", "array_max",
        "dominant_state", "hard_call", "weight_status", "row_receipt_sha256",
    ]
    with tsv_path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            law, env = row["normalized_measurement_likelihoods"], row["control_scale_uncertainty_envelope"]
            writer.writerow({
                **{field: row.get(field) for field in fields}, **law,
                "artifact_min": env[ARTIFACT]["minimum"], "artifact_max": env[ARTIFACT]["maximum"],
                "later_min": env[LATER]["minimum"], "later_max": env[LATER]["maximum"],
                "segdup_min": env[SEGDUP]["minimum"], "segdup_max": env[SEGDUP]["maximum"],
                "array_min": env[ARRAY]["minimum"], "array_max": env[ARRAY]["maximum"],
            })
    prototype = next(row for row in rows if row["candidate_key"] == "HG00423|HML-2_1q22|region2")
    panel_complete = not pending and len(rows) == expected_candidate_count
    receipt = {
        "schema": "hml2.cnv-final-assayed-weight-receipt.v1",
        "status": (
            f"COMPLETE_ALL_{expected_candidate_count}_CANDIDATE_AUTHORITIES_FINAL"
            if panel_complete else
            "FINAL_FOR_ALL_MATERIALIZED_ASSAYS_PANEL_STILL_OPEN"
        ),
        "panel_complete": panel_complete,
        "assayed_final_row_count": len(rows),
        "latent_depth_final_row_count": sum(row["authority_class"] != FIXED_SEGDUP for row in rows),
        "authenticated_segdup_final_row_count": sum(row["authority_class"] == FIXED_SEGDUP for row in rows),
        "pending_unmaterialized_candidate_count": len(pending),
        "pending_candidates": pending,
        "state_universe": list(STATES),
        "uncertainty": "min/max across 61-control narrow, central, and wide empirical residual scales",
        "weight_semantics": "normalized_measurement_likelihood_not_population_posterior",
        "hard_calls_emitted": 0,
        "inputs_sha256": inputs,
        "outputs_sha256": {"weights_tsv": sha256(tsv_path), "weights_json": sha256(json_path)},
        "ordered_row_receipt_sha256s": [row["row_receipt_sha256"] for row in rows],
        "hg00423_1q22_prototype": {
            "final_weight_row": prototype,
            "physical_ledger_adjudication": ledger["adjudication"],
            "accepted_latent_marker_receipt": marker_receipt["hg00423_1q22"],
            "conflict_resolution": "Duplicate assembly contigs are retained as evidence, but no unique junction is materialized; local-flank depth supplies the artifact-versus-later-duplication likelihood and cannot create a segdup or array route.",
        },
        "gates": {"production": False, "result": False, "manuscript": False},
    }
    (outdir / "FINAL_ASSAYED_WEIGHT_RECEIPT.v1.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    for name in ("base_authority", "site_accounting", "calibration",
                 "boundary_authority", "ready_index", "hg00423_ledger",
                 "latent_receipt"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--exact-evidence", type=Path)
    parser.add_argument("--ready-remote-manifest", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    base = json.loads(args.base_authority.read_text())
    records = base["records"]
    by_key = {row["candidate_key"]: row for row in records}
    if not records or len(by_key) != len(records):
        raise SystemExit("base authority has no records or repeats candidate keys")
    counts = {}
    for row in records:
        counts[row["authority_class"]] = counts.get(row["authority_class"], 0) + 1
    if not set(counts) <= {FIXED_SEGDUP, LATENT, UNMEASURED}:
        raise SystemExit("base authority contains an unsupported candidate class")
    targets = {
        row["candidate_key"]
        for row in records
        if row["authority_class"] == UNMEASURED
    }
    calibration = json.loads(args.calibration.read_text())["calibration"]
    if calibration["n_positive_controls"] != 61 or calibration["prior"] is not None:
        raise SystemExit("calibration is not the accepted prior-free 61-control model")
    scales = {key: float(value) for key, value in calibration["scales"].items()}
    boundary_rows = read_tsv(args.boundary_authority)
    measurements = boundary_measurements(boundary_rows, targets)
    fresh = exact_measurements(args.exact_evidence, targets)
    transferred = ready_remote_measurements(
        args.ready_remote_manifest, boundary_rows, targets
    )
    for key in set(fresh) & set(transferred):
        if fresh[key] != transferred[key]:
            raise SystemExit(f"conflicting transferred evidence for {key}")
    fresh.update(transferred)
    for key in set(measurements) & set(fresh):
        if measurements[key] != fresh[key]:
            raise SystemExit(f"conflicting evidence for {key}")
    measurements.update(fresh)
    ready = {row["candidate_id"]: row for row in read_tsv(args.ready_index)}
    if set(ready) != targets:
        raise SystemExit("ready index does not match the unmeasured candidate universe")
    sites = site_index(json.loads(args.site_accounting.read_text()))
    rows = [final_row(row, scales, sites, all_measurements=measurements)
            for row in records
            if row["authority_class"] in {LATENT, FIXED_SEGDUP}]
    rows.extend(final_row(
                    by_key[key], scales, sites, measurement,
                    all_measurements=measurements)
                for key, measurement in measurements.items())
    rows.sort(key=lambda row: row["candidate_key"])
    if len(rows) != len({row["candidate_key"] for row in rows}):
        raise SystemExit("duplicate output candidate")
    for row in rows:
        if not math.isclose(sum(row["normalized_measurement_likelihoods"].values()), 1.0,
                            rel_tol=0.0, abs_tol=1e-12):
            raise SystemExit(f"non-normalized row: {row['candidate_key']}")
    pending = [{
        "candidate_key": key,
        "scheduler_or_materialization_status": ready[key]["status"],
        "remote_assay_result_root": ready[key]["assay_result_root"],
    } for key in sorted(targets - set(measurements))]
    paths = {
        "base_authority": args.base_authority,
        "site_accounting": args.site_accounting,
        "calibration": args.calibration,
        "boundary_authority": args.boundary_authority,
        "ready_index": args.ready_index,
        "hg00423_ledger": args.hg00423_ledger,
        "latent_receipt": args.latent_receipt,
    }
    if args.exact_evidence:
        paths["exact_evidence"] = args.exact_evidence
    if args.ready_remote_manifest:
        paths["ready_remote_manifest"] = args.ready_remote_manifest
    write_outputs(args.outdir, rows, pending,
                  {name: sha256(path) for name, path in paths.items()},
                  json.loads(args.hg00423_ledger.read_text()),
                  json.loads(args.latent_receipt.read_text()),
                  len(records))


if __name__ == "__main__":
    main()
