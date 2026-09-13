#!/usr/bin/env python3
"""Audit whether extant Type-I loci retain greater Gag-Pro-Pol autonomy.

Pol is not treated as a standalone product.  The primary endpoint is the
canonical Gag -> Pro -> Pol translation route on the same proviral copy.
This script extracts the already-computed artifact-filtered, locus-level
comparisons and freezes the exact input identities used for interpretation.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMPARISONS = (
    ROOT
    / "working/type1_adjusted_matched_analysis_agent/results"
    / "adjusted_type_comparisons.tsv"
)
ORF_INPUT = (
    ROOT
    / "inputs/orf_analysis/current_v3r1"
    / "combined_hml2_orf_analysis.CNV_WEIGHTED.v3r1.tsv"
)
ARTIFACTS = ROOT / "manuscript/supplement/Table_S_CNV_excluded_ORF_records.tsv"
OUT = ROOT / "manuscript/delta292_autonomy_hypothesis_audit_v1"

ENDPOINTS = (
    "poly_e2_canonical_translation",
    "poly_combined_canonical_pol",
    "integrase_catalytic_motif",
    "integrase_core_sequence_compatibility",
)

SCOPES = (
    "all_typed_unadjusted",
    "ltr5hs_frozen_authority",
    "ltr5hs_carriage_callable_adjusted",
    "resolved_subfamily_carriage_callable_adjusted",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
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
    OUT.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(COMPARISONS)
    selected: list[dict[str, object]] = []
    keyed = {(r["analysis"], r["endpoint"]): r for r in rows}

    for scope in SCOPES:
        for endpoint in ENDPOINTS:
            row = keyed.get((scope, endpoint))
            if row is None:
                raise ValueError(f"missing required comparison: {scope} / {endpoint}")
            difference = float(row["typeI_minus_typeII"])
            selected.append(
                {
                    "analysis_scope": scope,
                    "endpoint": endpoint,
                    "pol_translation_semantics": (
                        "same-copy Gag-Pro-Pol route"
                        if endpoint == "poly_e2_canonical_translation"
                        else (
                            "combined Pol compatibility conditional on upstream route"
                            if endpoint == "poly_combined_canonical_pol"
                            else "sequence compatibility; not translation evidence"
                        )
                    ),
                    "n_loci_typeI": int(row["n_loci_typeI"]),
                    "n_loci_typeII": int(row["n_loci_typeII"]),
                    "mean_typeI": float(row["mean_typeI"]),
                    "mean_typeII": float(row["mean_typeII"]),
                    "typeI_minus_typeII": difference,
                    "direction_supports_greater_typeI_autonomy": difference > 0,
                    "permutation_p": float(row["permutation_p"]),
                    "permutation_bh_q": (
                        row["permutation_bh_q_within_family_and_analysis"]
                    ),
                    "model": row["model"],
                    "covariates": row["covariates"],
                }
            )

    primary = next(
        row
        for row in selected
        if row["analysis_scope"] == "ltr5hs_carriage_callable_adjusted"
        and row["endpoint"] == "poly_e2_canonical_translation"
    )
    combined_pol = next(
        row
        for row in selected
        if row["analysis_scope"] == "ltr5hs_carriage_callable_adjusted"
        and row["endpoint"] == "poly_combined_canonical_pol"
    )

    # The v3r1 ORF table contains all records; the downstream comparison input
    # was built after the exact-record exclusion table was applied.  Freeze
    # both identities so later reruns cannot silently switch inputs.
    n_orf_rows = sum(1 for _ in ORF_INPUT.open()) - 1
    n_artifact_rows = sum(1 for _ in ARTIFACTS.open()) - 1
    if n_orf_rows != 61_936 or n_artifact_rows != 35:
        raise ValueError(
            f"unexpected current input sizes: {n_orf_rows=} {n_artifact_rows=}"
        )

    summary = {
        "analysis": "delta292_autonomy_hypothesis_audit_v1",
        "question": (
            "Do extant Type-I loci show greater retention of the same-copy "
            "Gag-Pro-Pol translation route or Pol compatibility?"
        ),
        "primary_adjusted_result": primary,
        "combined_pol_adjusted_result": combined_pol,
        "verdict": (
            "No current fossil-sequence support for greater retained Type-I "
            "Gag-Pro-Pol autonomy. Both primary adjusted contrasts point in "
            "the opposite direction and are nonsignificant."
        ),
        "claim_boundary": (
            "This falsifies present-day ORF preservation as supporting evidence "
            "for the autonomy branch. It does not reconstruct ancestral founder "
            "ORFs and therefore does not prove that the propagating ancestor "
            "lacked a competent Gag-Pro-Pol route."
        ),
        "artifact_filter": {
            "terminal_orf_rows": n_orf_rows,
            "exact_excluded_records": n_artifact_rows,
            "biological_analysis_rows": n_orf_rows - n_artifact_rows,
        },
        "input_sha256": {
            str(COMPARISONS.relative_to(ROOT)): sha256(COMPARISONS),
            str(ORF_INPUT.relative_to(ROOT)): sha256(ORF_INPUT),
            str(ARTIFACTS.relative_to(ROOT)): sha256(ARTIFACTS),
        },
    }

    write_tsv(OUT / "autonomy_scope_summary.tsv", selected)
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    (OUT / "MODEL_SPECIFICATION.md").write_text(
        "# Type-I retained-autonomy audit\n\n"
        "The primary endpoint is `poly_e2_canonical_translation`, which requires "
        "a compatible Gag-to-Pro-to-Pol translation route on the same proviral "
        "copy. Pol is never treated as independently translated. Loci, not "
        "haplotype records, are the inferential units. The principal matched "
        "scope is frozen-authority LTR5Hs with carrier probability and callable "
        "denominator adjustment. Exact CNV/assembly-artifact records listed in "
        "`Table_S_CNV_excluded_ORF_records.tsv` are excluded upstream; samples "
        "and loci are retained.\n\n"
        "A negative extant-locus result is evidence against superior retained "
        "fossil ORF preservation. It is not a reconstruction of the ancestral "
        "propagating virus and is not used to force its ancestral state.\n"
    )
    print(OUT)


if __name__ == "__main__":
    main()
