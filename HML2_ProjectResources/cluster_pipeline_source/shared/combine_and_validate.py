import pandas as pd
import argparse
import os
import sys
import logging
import re
import tempfile
import json
from pathlib import Path


# --- Canonical non-carriage states --------------------------------------------------------------
# The master table encodes exactly two ways a provirus can fail to appear, and they are NOT
# interchangeable:
#
#   Insertion_Absent -- the sample/haplotype/locus WAS successfully interrogated and the provirus is
#       genuinely not there. A real biological observation (a non-carrier), retained in the output,
#       and the ONLY state that may enter a biological-absence denominator.
#   Technical_Error  -- the sample/haplotype/locus was NOT successfully interrogated (missing BAM or
#       index, graph-extraction failure, ORF failure, missing/empty upstream result). Absence here is
#       UNKNOWN. It always carries an explicit cause, is never relabelled as absence, and never
#       enters the absence denominator.
BIOLOGICAL_ABSENCE = "Insertion_Absent"
TECHNICAL_ERROR = "Technical_Error"

# Upstream Structure spellings meaning "interrogated, nothing here" -> canonicalized to
# Insertion_Absent.
BIOLOGICAL_ABSENCE_FLAGS = ["Absent", "Absent_or_Unaligned", "Insertion_Absent"]
# Upstream Structure spellings meaning "not interrogated" -> canonicalized to Technical_Error, each
# keeping its cause.
TECHNICAL_ERROR_FLAGS = ["BAM_or_Index_Error", "Graph_Extraction_Error", "ORF_Analysis_Error",
                         "Python_Mapping_Error", "Processing_Error", "Empty_File",
                         "No_FASTA_Records", "Technical_Error", "Technical_Missing"]

# The locus type selects the KCON reference and the accessory ORF (np9 vs rec). It is a closed set:
# anything else is a defect and must fail closed rather than be absorbed by a default.
VALID_LOCUS_TYPES = {"typei": "type1", "typeii": "type2"}
COPY_STATE_WEIGHT_COLS = ["copy_state_single_weight", "copy_state_artifact_weight", "copy_state_true_dup_weight"]
COPY_EXPECTATION_COLS = [
    "expected_intact_gag_copies", "expected_intact_pro_copies", "expected_intact_pol_copies",
    "expected_intact_env_copies", "expected_intact_accessory_copies",
    "expected_all_orfs_intact_copies",
]
COPY_STATE_LEDGER_COLUMNS = [
    "locus", "ID_Full", "observation_weight", *COPY_STATE_WEIGHT_COLS,
    "expected_biological_copy_count", *COPY_EXPECTATION_COLS,
    "copy_state_weight_provenance", "artifact_evidence_sha256", "artifact_evidence_bytes",
    "duplication_birth_constraint", "divergence_marginalization_model",
    "consequence_marginalization_provenance",
]


def reverse_complement(value):
    return value.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def normalize_present_day_tsd_semantics(rows):
    """Keep sample-specific TSD observations in insertion orientation."""
    result = rows.copy()
    added = (
        "oriented_tsd_pair", "tsd_length",
        "present_day_tsd_mismatch_count", "tsd_observable_state",
        "tsd_call_status", "tsd_pair_admission",
        "tsd_call_provenance", "tsd_source_interval",
        "tsd_evidence_role", "latent_insertion_tsd_relationship",
        "genomic_forward_left_tsd", "genomic_forward_right_tsd",
        "TSD_Orientation",
        "integrated_reference_build", "integrated_reference_interval",
        "integrated_sequence_authority", "orthologous_empty_site_authority",
    )
    for field in added:
        if field not in result.columns:
            result[field] = ""
    # ERROR/MAP_ERROR are old producer spellings for an unobservable junction,
    # not nucleotide strings.  Treat them exactly like a deleted/unspanned TSD
    # side: retain the provirus row, contribute no paired-TSD likelihood.
    unavailable = {
        "", "NA", "N/A", "NONE", "UNKNOWN", "UNOBSERVED", "-",
        "ERROR", "MAP_ERROR",
    }

    def censor_unproven_pair(index):
        result.at[index, "5'_TSD"] = "NONE"
        result.at[index, "3'_TSD"] = "NONE"
        result.at[index, "oriented_tsd_pair"] = "NONE/NONE"
        result.at[index, "tsd_length"] = ""
        result.at[index, "present_day_tsd_mismatch_count"] = ""
        result.at[index, "tsd_observable_state"] = "UNOBSERVABLE"
        result.at[index, "tsd_call_status"] = "NO_OBSERVABLE_TSD"
        result.at[index, "tsd_pair_admission"] = "UNAVAILABLE"
        result.at[index, "tsd_evidence_role"] = "NO_PAIRED_TSD_OBSERVATION"
        result.at[index, "latent_insertion_tsd_relationship"] = (
            "IDENTICAL_COPIES_AT_INSERTION"
        )
        result.at[index, "genomic_forward_left_tsd"] = "NONE"
        result.at[index, "genomic_forward_right_tsd"] = "NONE"

    for index, row in result.iterrows():
        # ``presence_call`` is assigned directly from the ORF result immediately
        # before this reducer runs.  Older ORF result CSVs do not carry the later
        # ``observation_state`` column, so requiring that extra field here erased
        # every real TSD in the combined table.  When a state is present it must
        # agree; when it is absent, the explicit presence call is authoritative.
        observation_state = str(row.get("observation_state", "")).upper()
        present = (
            str(row.get("presence_call", "")).lower() == "present"
            and observation_state in {"", "PRESENT"}
        )
        left = str(row.get("5'_TSD", "")).upper()
        right = str(row.get("3'_TSD", "")).upper()
        left = "" if left in unavailable else left
        right = "" if right in unavailable else right
        if not present:
            result.at[index, "5'_TSD"] = ""
            result.at[index, "3'_TSD"] = ""
            result.at[index, "oriented_tsd_pair"] = ""
            result.at[index, "tsd_length"] = ""
            result.at[index, "present_day_tsd_mismatch_count"] = ""
            result.at[index, "tsd_observable_state"] = "NOT_APPLICABLE"
            result.at[index, "tsd_call_status"] = "NO_COPY_OR_UNRESOLVED"
            result.at[index, "tsd_pair_admission"] = "UNAVAILABLE"
            result.at[index, "TSD_Orientation"] = "NOT_APPLICABLE"
            result.at[index, "tsd_call_provenance"] = (
                "ORF_ANALYSIS_ADJACENT_LTR_HOST_JUNCTIONS"
            )
            result.at[index, "tsd_source_interval"] = str(
                row.get("Source_Identifier", "")
            )
            result.at[index, "tsd_evidence_role"] = "NO_PAIRED_TSD_OBSERVATION"
            result.at[index, "latent_insertion_tsd_relationship"] = (
                "NOT_APPLICABLE_NO_COPY_OR_UNRESOLVED"
            )
            continue
        strand = str(row.get("Strand", "")).upper()
        strand = (
            "REVERSE" if strand in {"-", "REVERSE"} else
            "FORWARD" if strand in {"+", "FORWARD"} else "UNKNOWN"
        )
        result.at[index, "Strand"] = strand
        tsd_orientation = str(row.get("TSD_Orientation", "")).upper()
        if tsd_orientation == "INSERTION":
            if strand == "REVERSE":
                genomic_left, genomic_right = (
                    reverse_complement(right) if right else "",
                    reverse_complement(left) if left else "",
                )
            else:
                genomic_left, genomic_right = left, right
        else:
            genomic_left, genomic_right = left, right
            if strand == "REVERSE":
                left, right = (
                    reverse_complement(right) if right else "",
                    reverse_complement(left) if left else "",
                )
        result.at[index, "genomic_forward_left_tsd"] = (
            genomic_left if genomic_left else "NONE"
        )
        result.at[index, "genomic_forward_right_tsd"] = (
            genomic_right if genomic_right else "NONE"
        )
        result.at[index, "TSD_Orientation"] = "INSERTION"
        # Missing/deleted/unspanned sides are biological no-observations, not
        # empty serialization.  Keep the explicit sentinel in the published
        # combined table while the local booleans below drive admission.
        result.at[index, "5'_TSD"] = left if left else "NONE"
        result.at[index, "3'_TSD"] = right if right else "NONE"
        result.at[index, "oriented_tsd_pair"] = (
            f"{left or 'NONE'}/{right or 'NONE'}"
        )
        result.at[index, "tsd_call_provenance"] = (
            "ORF_ANALYSIS_ADJACENT_LTR_HOST_JUNCTIONS;"
            "LTR_TERMINAL_REFINEMENT_MAX_SHIFT_2BP;"
            "NO_HOST_KMER_BOUNDARY_SELECTION"
        )
        source_interval = str(row.get("tsd_source_interval", ""))
        if source_interval.upper() in unavailable:
            source_interval = str(row.get("Source_Identifier", ""))
        result.at[index, "tsd_source_interval"] = source_interval
        if left and right:
            if (
                len(left) != len(right)
                or len(left) not in (4, 5, 6)
                or set(left + right) - set("ACGT")
            ):
                censor_unproven_pair(index)
                continue
            mismatch = sum(a != b for a, b in zip(left, right))
            if mismatch > 2:
                censor_unproven_pair(index)
                continue
            result.at[index, "tsd_length"] = str(len(left))
            result.at[index, "present_day_tsd_mismatch_count"] = str(mismatch)
            result.at[index, "tsd_observable_state"] = "TWO_SIDED"
            result.at[index, "tsd_call_status"] = "PAIRED_ACCEPTED"
            result.at[index, "tsd_pair_admission"] = "PAIRED_ACCEPTED"
            result.at[index, "tsd_evidence_role"] = (
                "POST_INSERTION_EVOLUTION_OBSERVATION"
            )
        elif bool(left) != bool(right):
            result.at[index, "tsd_length"] = ""
            result.at[index, "present_day_tsd_mismatch_count"] = ""
            result.at[index, "tsd_observable_state"] = "ONE_SIDED"
            result.at[index, "tsd_call_status"] = "ONE_SIDED_NO_PAIR"
            result.at[index, "tsd_pair_admission"] = "ONE_SIDED"
            result.at[index, "tsd_evidence_role"] = "NO_PAIRED_TSD_OBSERVATION"
        else:
            result.at[index, "tsd_length"] = ""
            result.at[index, "present_day_tsd_mismatch_count"] = ""
            result.at[index, "tsd_observable_state"] = "UNOBSERVABLE"
            result.at[index, "tsd_call_status"] = "NO_OBSERVABLE_TSD"
            result.at[index, "tsd_pair_admission"] = "UNAVAILABLE"
            result.at[index, "tsd_evidence_role"] = "NO_PAIRED_TSD_OBSERVATION"
        result.at[index, "latent_insertion_tsd_relationship"] = (
            "IDENTICAL_COPIES_AT_INSERTION"
        )
    return result


def load_copy_state_ledger(path):
    """Load exact observation weights; validate probabilities without inventing states."""
    ledger = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    # The current weighted combined ORF table is itself the direct biological
    # copy-state authority.  Accept it without materializing a second ledger;
    # extra ORF/TSD columns are ignored and loci outside the current 101-locus
    # source partition are never consumed.
    combined_source = "locus" not in ledger.columns and "Locus" in ledger.columns
    if combined_source:
        ledger = ledger.rename(columns={"Locus": "locus"})
        if "presence_call" in ledger.columns:
            ledger = ledger[ledger["presence_call"].eq("present")].copy()
    missing = [c for c in COPY_STATE_LEDGER_COLUMNS if c not in ledger.columns]
    if missing:
        raise ValueError(f"copy-state ledger is missing biological fields: {missing}")
    ledger = ledger[COPY_STATE_LEDGER_COLUMNS].copy()
    if ledger.duplicated(["locus", "ID_Full"]).any():
        dup = ledger.loc[ledger.duplicated(["locus", "ID_Full"], keep=False), ["locus", "ID_Full"]]
        raise ValueError(f"duplicate exact copy-state rows: {dup.to_dict('records')[:10]}")
    if ledger[["locus", "ID_Full"]].apply(lambda s: s.str.strip().eq("")).any().any():
        raise ValueError("copy-state locus and ID_Full must be nonempty exact identifiers")
    numeric = ["observation_weight", *COPY_STATE_WEIGHT_COLS, "expected_biological_copy_count",
               *COPY_EXPECTATION_COLS, "artifact_evidence_bytes"]
    for col in numeric:
        ledger[col] = pd.to_numeric(ledger[col], errors="raise")
    state_sum = ledger[COPY_STATE_WEIGHT_COLS].sum(axis=1)
    if ((ledger[["observation_weight", *COPY_STATE_WEIGHT_COLS]] < 0).any().any()
            or (ledger[["observation_weight", *COPY_STATE_WEIGHT_COLS]] > 1).any().any()
            or (state_sum - 1.0).abs().gt(1e-9).any()):
        raise ValueError("copy-state and observation weights must be within [0,1]; state weights must sum to 1 per observation")
    if (ledger["expected_biological_copy_count"] < 0).any():
        raise ValueError("expected_biological_copy_count must be nonnegative")
    for col in COPY_EXPECTATION_COLS:
        if (ledger[col] < 0).any() or (ledger[col] > ledger["expected_biological_copy_count"] + 1e-9).any():
            raise ValueError(f"{col} must be within [0, expected_biological_copy_count]")
    for col in ("copy_state_weight_provenance", "consequence_marginalization_provenance"):
        if ledger[col].str.strip().isin(["", "NA", "UNKNOWN", "UNMEASURED"]).any():
            raise ValueError(f"{col} must name measured model provenance")
    true_dup = ledger["copy_state_true_dup_weight"] > 0
    artifact = ledger["copy_state_artifact_weight"] > 0
    if (ledger.loc[artifact, "artifact_evidence_bytes"] <= 0).any():
        raise ValueError("artifact weight requires positive retained artifact_evidence_bytes")
    if (ledger.loc[true_dup, "duplication_birth_constraint"] != "IDENTICAL_AT_DUPLICATION_BIRTH").any():
        raise ValueError("true-dup weight requires IDENTICAL_AT_DUPLICATION_BIRTH coupling")
    if ledger.loc[true_dup, "divergence_marginalization_model"].str.strip().isin(["", "NA", "UNKNOWN", "UNMEASURED"]).any():
        raise ValueError("true-dup weight requires an explicit later-divergence marginalization model")
    return ledger


def classify_presence(structure, qc_state):
    """Classify one upstream row as present / absent_noncarrier / technical_missing.

    QC_State is authoritative when it declares a technical error, so a Structure spelling this
    script has not seen before can never be silently absorbed into biological absence.
    """
    s = str(structure).strip()
    q = str(qc_state).strip()
    if q == TECHNICAL_ERROR or s in TECHNICAL_ERROR_FLAGS:
        return "technical_missing"
    if s in BIOLOGICAL_ABSENCE_FLAGS:
        return "absent_noncarrier"
    return "present"


def technical_reason(structure, qc_reason):
    """The explicit cause of a technical failure. Never dropped, never left blank."""
    reason = str(qc_reason).strip()
    if reason and reason.lower() not in ("nan", "na", "none"):
        return reason
    s = str(structure).strip()
    if s and s.lower() not in ("nan", "na", "none") and s != TECHNICAL_ERROR:
        return s.upper()
    return "UNSPECIFIED_TECHNICAL_ERROR"


def atomic_write_dataframe(df, out_path, sep=",", **to_csv_kwargs):
    """Write a DataFrame to out_path ATOMICALLY.

    The combined output must never be observed half-written or interleaved: a reader (or a restart
    of a downstream step) must see either the complete previous file or the complete new one, never a
    truncated mix. We therefore serialize to a uniquely-named temp file in the SAME directory (so the
    final rename stays on one filesystem), flush+fsync it to durable storage, then os.replace() it
    onto the destination -- an atomic rename on POSIX. This also makes the step restartable: a crash
    mid-write leaves the prior good file intact and only an orphan temp file (cleaned up on failure).
    """
    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    os.makedirs(out_dir, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(prefix="." + os.path.basename(out_path) + ".", suffix=".tmp", dir=out_dir)
    try:
        with os.fdopen(fd, "w", newline="") as fh:
            df.to_csv(fh, sep=sep, **to_csv_kwargs)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp_path, out_path)   # atomic on POSIX (same filesystem)
    except BaseException:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise


TSD_PUBLISHED_COLUMNS = [
    "5'_TSD", "3'_TSD", "oriented_tsd_pair", "tsd_length",
    "present_day_tsd_mismatch_count", "tsd_observable_state",
    "tsd_call_status", "tsd_pair_admission", "tsd_call_provenance",
    "tsd_source_interval", "tsd_evidence_role",
    "latent_insertion_tsd_relationship", "genomic_forward_left_tsd",
    "genomic_forward_right_tsd", "TSD_Orientation",
    "integrated_reference_build", "integrated_reference_interval",
    "integrated_sequence_authority", "orthologous_empty_site_authority",
]

DEFINITIVE_NO_TSD_3Q12_3_FRAGMENT_IDS = frozenset({
    "HG02451_pat_hprc_r2_v1.0.1_HML-2_3q12.3",
    "HG02559_pat_hprc_r2_v1.0.1_HML-2_3q12.3",
    "HG02583_hap1_hprc_r2_v1.0.1_HML-2_3q12.3",
    "HG03579_pat_hprc_r2_v1.0.1_HML-2_3q12.3",
    "NA18879_hap2_hprc_r2_v1.0.1_HML-2_3q12.3",
    "NA19338_hap1_hprc_r2_v1.0.1_HML-2_3q12.3",
})
TSD_MERGE_COLUMNS = TSD_PUBLISHED_COLUMNS[:15]

TERMINAL_STATE_MUTABLE_COLUMNS = frozenset({
    "Structure", "Strand", "arm_resolved", "presence_call",
    *TSD_MERGE_COLUMNS,
})
TERMINAL_STATE_SCAN_COLUMNS = (
    "Structure", "Strand", "arm_resolved", "presence_call",
    "tsd_observable_state", "tsd_call_status", "tsd_pair_admission",
    "tsd_call_provenance", "tsd_evidence_role",
    "latent_insertion_tsd_relationship", "TSD_Orientation",
)
NONTERMINAL_STATE_PATTERN = re.compile(
    r"NO_MATCH|CENSORED|NO_DIRECT|UNRESOLVED|UNKNOWN|FALLBACK|"
    r"NO_[A-Z0-9_]*EVIDENCE",
    re.IGNORECASE,
)
NONTERMINAL_EXACT_VALUES = frozenset({
    "UNAVAILABLE",
    "TECHNICAL_MISSING",
    "TECHNICAL_FAILURE",
    "NO_COPY_OR_UNRESOLVED",
    "NOT_APPLICABLE_NO_COPY_OR_UNRESOLVED",
})
AGGREGATE_TSD_LOCI = frozenset({"HML-2_acro_type1", "HML-2_acro_type2"})
ASM_GAP_EXACT_JUNCTION_REQUIRED_ID = (
    "HG00658_pat_hprc_r2_v1.0.1_HML-2_7p22.1_asmdup"
)


def _required_authority_columns(table, columns, label):
    missing = sorted(set(columns) - set(table.columns))
    if missing:
        raise ValueError(f"{label} lacks required columns: {missing}")


def _authority_scalar(value):
    return str(value or "").strip()


def _nonterminal_value_counts(table):
    counts = {}
    for field in TERMINAL_STATE_SCAN_COLUMNS:
        if field not in table.columns:
            continue
        for value, count in table[field].astype(str).value_counts().items():
            if value.upper() in NONTERMINAL_EXACT_VALUES or NONTERMINAL_STATE_PATTERN.search(value):
                counts[f"{field}={value}"] = int(count)
    return dict(sorted(counts.items()))


def terminalize_published_tsd_states(
    rows,
    paired_authority,
    callability_authority,
    evidence_review,
):
    """Replace every generic TSD/QC fallthrough with an authority-backed state.

    The paired table authenticates the biological observation class for every
    row.  Non-carriers and missing observations then require an exact
    callability decision.  Sequence-present no-calls are separated by retained
    junction review: unsupported outward boundary, one resolved terminal
    boundary only, or a concrete retained-junction observation gap.  Any row
    outside these closed classes fails the complete materialization.
    """
    result = rows.copy().reset_index(drop=True)
    paired = paired_authority.copy().reset_index(drop=True)
    callability = callability_authority.copy().reset_index(drop=True)
    review = evidence_review.copy().reset_index(drop=True)
    _required_authority_columns(
        result,
        {
            "Locus", "orig_Locus", "ID_Full", "Structure", "Strand",
            "arm_resolved", "presence_call", *TSD_MERGE_COLUMNS,
        },
        "combined ORF table",
    )
    _required_authority_columns(
        paired,
        {"Locus", "ID_Full", "presence_call", "observation_state"},
        "paired TSD authority",
    )
    _required_authority_columns(
        callability,
        {
            "Locus", "ID_Full", "decision", "match_scope",
            "n_matching_raw_records", "raw_classifications", "raw_paths",
        },
        "ORF callability authority",
    )
    _required_authority_columns(
        review,
        {
            "candidate_locus", "ID_Full", "evidence_resolution_status",
            "evidence_5prime_tsd", "evidence_3prime_tsd",
            "left_terminal_support", "right_terminal_support",
            "review_flag",
        },
        "direct junction review authority",
    )
    if len(result) != len(paired):
        raise ValueError(
            "combined/paired authority row counts differ: "
            f"{len(result)} != {len(paired)}"
        )

    projected_loci = result["Locus"].where(
        ~result["Locus"].isin(AGGREGATE_TSD_LOCI),
        result["orig_Locus"],
    )
    combined_keys = list(zip(projected_loci.astype(str), result["ID_Full"].astype(str)))
    paired_keys = list(zip(paired["Locus"].astype(str), paired["ID_Full"].astype(str)))
    for position, (combined_key, paired_key) in enumerate(
        zip(combined_keys, paired_keys), start=2
    ):
        if combined_key != paired_key:
            raise ValueError(
                "combined/paired authority order differs at row "
                f"{position}: {combined_key} != {paired_key}"
            )
    paired_presence = paired["presence_call"].astype(str).str.strip().str.lower()
    expected_presence = {"present", "absent_noncarrier", "technical_missing"}
    unexpected_presence = sorted(set(paired_presence) - expected_presence)
    if unexpected_presence:
        raise ValueError(
            f"paired authority contains unknown presence calls: {unexpected_presence}"
        )

    callability_by_key = {}
    callability_fields = (
        "decision", "match_scope", "n_matching_raw_records",
        "raw_classifications", "raw_paths",
    )
    for record in callability.to_dict("records"):
        key = (_authority_scalar(record["Locus"]), _authority_scalar(record["ID_Full"]))
        authority = tuple(_authority_scalar(record[field]) for field in callability_fields)
        prior = callability_by_key.get(key)
        if prior is not None and prior != authority:
            raise ValueError(f"conflicting callability authority for {key}")
        callability_by_key[key] = authority

    review_by_key = {}
    for record in review.to_dict("records"):
        key = (
            _authority_scalar(record["candidate_locus"]),
            _authority_scalar(record["ID_Full"]),
        )
        if key in review_by_key:
            raise ValueError(f"duplicate direct junction review authority for {key}")
        review_by_key[key] = record

    before_nonterminal = _nonterminal_value_counts(result)
    before = result.copy()
    classifications = {
        "authenticated_noncarrier": 0,
        "recoverable_technical_zero_raw_orf_records": 0,
        "sequence_present_outward_terminal_unsupported": 0,
        "sequence_present_single_terminal_boundary_only": 0,
        "sequence_present_retained_junction_observation_absent": 0,
        "sequence_present_assembly_gap_exact_junction_absent": 0,
        "sequence_present_established_length_pair_rejected": 0,
        "sequence_present_gt2_pair_rejected": 0,
        "sequence_present_inverted_ltr_geometry": 0,
        "aggregate_arm_not_assigned": 0,
    }

    def set_no_call(
        index,
        status,
        admission,
        provenance,
        evidence_role,
        observable_state="UNOBSERVABLE",
        orientation="INSERTION",
    ):
        result.at[index, "5'_TSD"] = "NONE"
        result.at[index, "3'_TSD"] = "NONE"
        result.at[index, "oriented_tsd_pair"] = "NONE/NONE"
        result.at[index, "tsd_length"] = ""
        result.at[index, "present_day_tsd_mismatch_count"] = ""
        result.at[index, "tsd_observable_state"] = observable_state
        result.at[index, "tsd_call_status"] = status
        result.at[index, "tsd_pair_admission"] = admission
        result.at[index, "tsd_call_provenance"] = provenance
        result.at[index, "tsd_evidence_role"] = evidence_role
        result.at[index, "genomic_forward_left_tsd"] = "NONE"
        result.at[index, "genomic_forward_right_tsd"] = "NONE"
        result.at[index, "TSD_Orientation"] = orientation

    for index, key in enumerate(combined_keys):
        presence = paired_presence.iloc[index]
        result.at[index, "presence_call"] = presence
        if presence == "absent_noncarrier":
            authority = callability_by_key.get(key)
            if authority is None:
                raise ValueError(f"non-carrier lacks exact callability authority: {key}")
            decision, scope, record_count, raw_classification, raw_path = authority
            if (
                decision != "pipeline_noncarrier"
                or scope != "exact_stem"
                or record_count != "1"
                or raw_classification != "pipeline_noncarrier_sentinel"
                or not raw_path
            ):
                raise ValueError(
                    f"non-carrier authority is not terminal for {key}: {authority}"
                )
            set_no_call(
                index,
                "AUTHENTICATED_NONCARRIER_EMPTY_SITE",
                "NOT_APPLICABLE_AUTHENTICATED_NONCARRIER",
                "AUTHENTICATED_PIPELINE_NONCARRIER_SENTINEL_EXACT_STEM",
                "AUTHENTICATED_EMPTY_SITE_NONCARRIER",
                "NOT_APPLICABLE_AUTHENTICATED_NONCARRIER",
                "NOT_APPLICABLE",
            )
            result.at[index, "latent_insertion_tsd_relationship"] = (
                "NOT_APPLICABLE_AUTHENTICATED_NONCARRIER"
            )
            result.at[index, "tsd_source_interval"] = raw_path
            result.at[index, "Strand"] = "NOT_APPLICABLE"
            result.at[index, "presence_call"] = "absent_noncarrier"
            classifications["authenticated_noncarrier"] += 1
            continue
        if presence == "technical_missing":
            authority = callability_by_key.get(key)
            if authority is None:
                raise ValueError(
                    f"technical observation lacks exact callability authority: {key}"
                )
            decision, scope, record_count, raw_classification, raw_path = authority
            if (
                decision != "synthetic_missing_no_raw_evidence"
                or scope != "none"
                or record_count != "0"
                or raw_classification
                or raw_path
            ):
                raise ValueError(
                    f"technical callability cause is not exact for {key}: {authority}"
                )
            set_no_call(
                index,
                "RECOVERABLE_TECHNICAL_ZERO_RAW_ORF_RECORDS_TSD_NOT_ASSESSED",
                "NOT_ASSESSED_TECHNICAL_ZERO_RAW_ORF_RECORDS",
                (
                    "AUTHENTICATED_CALLABILITY_ZERO_RAW_ORF_RECORDS;"
                    "RECOVERY_REQUIRES_SOURCE_SEQUENCE_OR_ORF_RERUN"
                ),
                "TECHNICAL_ZERO_RAW_ORF_RECORDS_TSD_NOT_ASSESSED",
                "NOT_ASSESSED_TECHNICAL_ZERO_RAW_ORF_RECORDS",
                "NOT_APPLICABLE",
            )
            result.at[index, "latent_insertion_tsd_relationship"] = (
                "NOT_ASSESSED_TECHNICAL_ZERO_RAW_ORF_RECORDS"
            )
            result.at[index, "tsd_source_interval"] = (
                "CALLABILITY_AUTHORITY_EXACT_LOCUS_AND_ID"
            )
            result.at[index, "Structure"] = "Technical_Zero_Raw_ORF_Records"
            result.at[index, "Strand"] = "NOT_ASSESSED_ZERO_RAW_ORF_RECORDS"
            result.at[index, "presence_call"] = (
                "recoverable_technical_zero_raw_orf_records"
            )
            classifications["recoverable_technical_zero_raw_orf_records"] += 1
            continue

        provenance = _authority_scalar(result.at[index, "tsd_call_provenance"])
        status = _authority_scalar(result.at[index, "tsd_call_status"])
        if provenance in {
            "NO_DIRECT_RESOLVED_ADJACENT_JUNCTION_EVIDENCE_CENSORED",
            "NO_DIRECT_EXACT_JUNCTION_EVIDENCE_CENSORED",
        }:
            review_record = review_by_key.get(
                (_authority_scalar(result.at[index, "Locus"]), key[1])
            )
            if review_record is None:
                if key[1] == ASM_GAP_EXACT_JUNCTION_REQUIRED_ID:
                    set_no_call(
                        index,
                        "SEQUENCE_PRESENT_ASSEMBLY_GAP_EXACT_JUNCTION_NOT_ASSESSED",
                        "NOT_ASSESSED_ASSEMBLY_GAP_EXACT_JUNCTION_NOT_RETAINED",
                        (
                            "AUTHENTICATED_SEQUENCE_PRESENT_ASSEMBLY_GAP_REPRESENTATIVE;"
                            "EXACT_INSERTION_JUNCTIONS_NOT_RETAINED;"
                            "BASE_LOCUS_TSD_INHERITANCE_FORBIDDEN"
                        ),
                        "ASSEMBLY_GAP_EXACT_JUNCTION_OBSERVATION_NOT_RETAINED",
                    )
                    classifications[
                        "sequence_present_assembly_gap_exact_junction_absent"
                    ] += 1
                else:
                    set_no_call(
                        index,
                        "SEQUENCE_PRESENT_TSD_BOUNDARY_NOT_ASSESSED",
                        "NOT_ASSESSED_RETAINED_JUNCTION_OBSERVATION_ABSENT",
                        (
                            "AUTHENTICATED_SEQUENCE_PRESENT;"
                            "RETAINED_DIRECT_JUNCTION_OBSERVATION_ABSENT;"
                            "TSD_BOUNDARY_NOT_ASSESSED"
                        ),
                        "SEQUENCE_PRESENT_RETAINED_JUNCTION_OBSERVATION_ABSENT",
                    )
                    classifications[
                        "sequence_present_retained_junction_observation_absent"
                    ] += 1
            else:
                resolution = _authority_scalar(
                    review_record["evidence_resolution_status"]
                )
                pair = (
                    _authority_scalar(review_record["evidence_5prime_tsd"]),
                    _authority_scalar(review_record["evidence_3prime_tsd"]),
                )
                support = (
                    _authority_scalar(review_record["left_terminal_support"]),
                    _authority_scalar(review_record["right_terminal_support"]),
                )
                if (
                    resolution == "OUTWARD_TERMINAL_NOT_SUPPORTED"
                    and pair == ("NONE", "NONE")
                    and support == ("UNSUPPORTED", "UNSUPPORTED")
                ):
                    set_no_call(
                        index,
                        "SEQUENCE_PRESENT_OUTWARD_TERMINAL_BOUNDARY_NOT_ESTABLISHED",
                        "NOT_ADMITTED_OUTWARD_TERMINAL_BOUNDARY_UNSUPPORTED",
                        (
                            "AUTHENTICATED_SEQUENCE_PRESENT;"
                            "OUTWARD_TERMINAL_BOUNDARY_UNSUPPORTED_BY_RETAINED_ALIGNMENT"
                        ),
                        "SEQUENCE_PRESENT_OUTWARD_TERMINAL_BOUNDARY_UNSUPPORTED",
                    )
                    classifications[
                        "sequence_present_outward_terminal_unsupported"
                    ] += 1
                elif (
                    resolution == "RESOLVED"
                    and pair == ("NONE", "NONE")
                    and support == ("DIRECT_PAF_ENDPOINT", "UNSUPPORTED")
                ):
                    set_no_call(
                        index,
                        "SEQUENCE_PRESENT_SINGLE_TERMINAL_BOUNDARY_INSUFFICIENT_FOR_TSD",
                        "NOT_ADMITTED_SINGLE_TERMINAL_BOUNDARY_ONLY",
                        (
                            "AUTHENTICATED_SEQUENCE_PRESENT;"
                            "SINGLE_TERMINAL_HOST_BOUNDARY_RESOLVED;"
                            "OPPOSITE_BOUNDARY_UNSUPPORTED;TSD_PAIR_NOT_OBSERVABLE"
                        ),
                        "ONE_TERMINAL_BOUNDARY_ONLY_TSD_PAIR_NOT_OBSERVABLE",
                    )
                    classifications[
                        "sequence_present_single_terminal_boundary_only"
                    ] += 1
                else:
                    raise ValueError(
                        "unclassified direct junction no-call for "
                        f"{key}: resolution={resolution} pair={pair} support={support}"
                    )
            result.at[index, "latent_insertion_tsd_relationship"] = (
                "IDENTICAL_COPIES_AT_INSERTION"
            )
            continue

        if "SHORTER_FALLBACK_NOT_ACCEPTED" in provenance:
            if not provenance.startswith(
                "DIRECT_RESOLVED_BOUNDARIES_FAIL_ESTABLISHED_"
            ):
                raise ValueError(f"unclassified shorter-length rejection for {key}")
            set_no_call(
                index,
                "SEQUENCE_PRESENT_TSD_PAIR_REJECTED_ESTABLISHED_LOCUS_LENGTH_GT2_MISMATCH",
                "REJECTED_ESTABLISHED_LOCUS_LENGTH_GT2_MISMATCH",
                provenance.replace(
                    "SHORTER_FALLBACK_NOT_ACCEPTED",
                    "SHORTER_LENGTH_ADMISSION_FORBIDDEN",
                ),
                "SEQUENCE_PRESENT_BOUNDARY_PAIR_REJECTED_BIOLOGICAL_CONSTRAINT",
            )
            classifications[
                "sequence_present_established_length_pair_rejected"
            ] += 1
            continue
        if provenance == (
            "RAW_FULL_GRAPH_BOUNDARIES_GGAAA_GCTAT;PAIRED_REJECTED_GT2;"
            "NO_ACCEPTED_PAIRED_TSD"
        ):
            set_no_call(
                index,
                "SEQUENCE_PRESENT_TSD_PAIR_REJECTED_GT2_MISMATCH",
                "REJECTED_PRESENT_DAY_GT2_MISMATCH",
                (
                    "RAW_FULL_GRAPH_BOUNDARIES_GGAAA_GCTAT;"
                    "PAIR_REJECTED_GT2_MISMATCH;EXCLUDED_FROM_ACCEPTED_TSD_SET"
                ),
                "SEQUENCE_PRESENT_BOUNDARY_PAIR_REJECTED_BIOLOGICAL_CONSTRAINT",
            )
            classifications["sequence_present_gt2_pair_rejected"] += 1
            continue
        if provenance == "INVERTED_LTR_OUTER_FACES_NO_PAIRED_INSERTION_BOUNDARY":
            set_no_call(
                index,
                "SEQUENCE_PRESENT_INVERTED_LTR_GEOMETRY_TSD_NOT_OBSERVABLE",
                "NOT_APPLICABLE_INVERTED_LTR_GEOMETRY",
                "INVERTED_LTR_OUTER_FACES;PAIRING_GEOMETRY_NOT_APPLICABLE",
                "SEQUENCE_PRESENT_INVERTED_LTR_GEOMETRY",
            )
            classifications["sequence_present_inverted_ltr_geometry"] += 1
            continue
        if status == "ONE_SIDED_NO_PAIR":
            result.at[index, "tsd_evidence_role"] = (
                "ONE_SIDED_TERMINAL_BOUNDARY_OBSERVATION"
            )

    unresolved_arm = result["arm_resolved"].astype(str).str.lower().eq("unresolved")
    invalid_unresolved_arm = unresolved_arm & ~result["Locus"].isin(AGGREGATE_TSD_LOCI)
    if invalid_unresolved_arm.any():
        keys = result.loc[
            invalid_unresolved_arm, ["Locus", "ID_Full"]
        ].head(10).to_dict("records")
        raise ValueError(f"non-aggregate rows retain unresolved arm state: {keys}")
    result.loc[unresolved_arm, "arm_resolved"] = "AGGREGATE_ARM_NOT_ASSIGNED"
    classifications["aggregate_arm_not_assigned"] = int(unresolved_arm.sum())

    after_nonterminal = _nonterminal_value_counts(result)
    if after_nonterminal:
        raise ValueError(
            f"published table retains nonterminal states: {after_nonterminal}"
        )
    changed_fields = {}
    for field in result.columns:
        count = int((result[field].astype(str) != before[field].astype(str)).sum())
        if count:
            if field not in TERMINAL_STATE_MUTABLE_COLUMNS:
                raise ValueError(
                    f"terminal-state repair changed unauthorized field {field}: {count}"
                )
            changed_fields[field] = count
    report = {
        "rows": len(result),
        "classifications": dict(sorted(classifications.items())),
        "changed_fields": dict(sorted(changed_fields.items())),
        "nonterminal_value_counts_before": before_nonterminal,
        "nonterminal_value_counts_after": after_nonterminal,
        "unauthorized_field_changes": 0,
    }
    return result, report


def _evidence_identity(value):
    return str(value or "").lstrip(">").split("|", 1)[0]


def _evidence_fasta_stem(value):
    name = os.path.basename(str(value or ""))
    for suffix in (".fasta.gz", ".fa.gz", ".fasta", ".fa"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return os.path.splitext(name)[0]


EVIDENCE_SUFFIX_CHAIN = re.compile(
    r"(?:_alt\d+|_(?:MULTI|DOUBLE)|_part\d+|_asmdup)+$"
)


def _evidence_base_identity(value):
    return EVIDENCE_SUFFIX_CHAIN.sub("", _evidence_identity(value))


def _direct_tsd_is_adjacent(payload, side, tsd):
    if not tsd:
        return True
    context_field = (
        "left_junction_context_in_insertion_orientation"
        if side == "5" else "right_junction_context_in_insertion_orientation"
    )
    offset_field = "left_context_boundary_offset" if side == "5" else "right_context_boundary_offset"
    context = str(payload.get(context_field, "")).upper()
    try:
        boundary = int(payload.get(offset_field, 18))
    except (TypeError, ValueError):
        return False
    if side == "5":
        return boundary >= len(tsd) and context[boundary - len(tsd):boundary] == tsd
    return len(context) >= boundary + len(tsd) and context[boundary:boundary + len(tsd)] == tsd


def apply_direct_tsd_evidence(current, evidence_root):
    """Replace every present-row TSD from resolved adjacent junction evidence.

    Exact observation identity is required. Tandem ``_part#`` rows alone may
    use their array-base identity because all parts share the array's outer
    insertion boundary. A missing, ambiguous, unresolved, non-LTR-owned,
    nonadjacent, malformed, or discordant observation is an explicit no-call.
    """
    root = Path(evidence_root)
    paths = sorted(root.glob("HML-2_*/junction_evidence.jsonl"))
    current_loci = set(current["Locus"].astype(str))
    evidence_loci = {path.parent.name for path in paths}
    unknown_loci = evidence_loci - current_loci
    if unknown_loci:
        raise ValueError(
            "direct TSD evidence contains loci outside the current 101-locus table: "
            f"{sorted(unknown_loci)}"
        )
    evidence = {}
    for path in paths:
        locus = path.parent.name
        with open(path) as handle:
            for line in handle:
                if not line.strip():
                    continue
                payload = json.loads(line)
                identities = {
                    _evidence_identity(payload.get("source_sequence_id", "")),
                    _evidence_identity(_evidence_fasta_stem(payload.get("source_fasta", ""))),
                }
                identities.discard("")
                for identity in identities:
                    evidence.setdefault((locus, identity), []).append(payload)

    repaired = current.copy()
    present_mask = repaired["presence_call"].astype(str).str.lower().eq("present")
    counts = {
        "evidence_files_present": len(paths),
        "evidence_loci_without_file_no_signal": len(current_loci - evidence_loci),
        "present_rows": int(present_mask.sum()),
        "exact_supported": 0,
        "tandem_array_base_supported": 0,
        "exact_paired_supported": 0,
        "exact_one_sided_supported": 0,
        "grch38_exact_coordinate_supported": 0,
        "no_exact_signal": 0,
        "ambiguous_exact_signal": 0,
        "rejected_diagnostic_pair": 0,
    }
    exact_supported = []
    for index, row in repaired.loc[present_mask].iterrows():
        key = (str(row["Locus"]), _evidence_identity(row["ID_Full"]))
        matches = evidence.get(key, [])
        used_tandem_base = False
        if not matches and re.search(r"_part\d+", str(row["ID_Full"])):
            matches = evidence.get(
                (str(row["Locus"]), _evidence_base_identity(row["ID_Full"])),
                [],
            )
            used_tandem_base = bool(matches)
        is_grch38_reference = (
            str(row["Locus"]) == "HML-2_8q24.3c"
            and re.search(
                r"GCA_000001405\.15|GRCh38",
                str(row["ID_Full"]),
                re.IGNORECASE,
            )
        )
        if is_grch38_reference:
            matches = [
                payload for payload in matches
                if payload.get("evidence_authority")
                == "GRCH38_EXACT_INTEGRATED_COORDINATES"
            ]
        pairs = set()
        for payload in matches:
            if payload.get("resolution_status", "") != "RESOLVED":
                continue
            grch38_coordinate_evidence = (
                is_grch38_reference
                and payload.get("boundary_owner")
                == "GRCH38_EXACT_INTEGRATED_ALLELE_COORDINATES"
                and payload.get("reference_build") == "GRCh38"
                and payload.get("left_reference_interval")
                == "chr8:145021238-145021243"
                and payload.get("right_reference_interval")
                == "chr8:145028836-145028841"
                and payload.get("tsd_source_interval")
                == (
                    "chr8:145021238-145021243;"
                    "chr8:145028836-145028841"
                )
            )
            if (
                payload.get("boundary_owner") != "TERMINAL_LTR_ALIGNMENT_ONLY"
                and not grch38_coordinate_evidence
            ):
                continue
            candidate_structure = str(row.get("Structure", "")).upper()
            evidence_structure = str(payload.get("structure", "")).upper()
            expected_structure = (
                "SOLO_LTR" if candidate_structure == "SOLO-LTR"
                else "PROVIRUS_OR_ARRAY"
            )
            if evidence_structure != expected_structure:
                continue
            invalid_motif_shift = False
            for side in ("left", "right"):
                if payload.get(f"{side}_terminal_support") != "TERMINAL_MOTIF":
                    continue
                try:
                    shift = int(payload.get(f"{side}_boundary_shift"))
                except (TypeError, ValueError):
                    invalid_motif_shift = True
                    break
                if shift not in {-2, -1, 0, 1, 2}:
                    invalid_motif_shift = True
                    break
            if invalid_motif_shift:
                continue
            left = str(payload.get("tsd_5prime_in_insertion_orientation", "")).upper()
            right = str(payload.get("tsd_3prime_in_insertion_orientation", "")).upper()
            left = "" if left in {"", "NA", "NONE"} else left
            right = "" if right in {"", "NA", "NONE"} else right
            if any(
                len(value) not in {4, 5, 6} or set(value) - set("ACGT")
                for value in (left, right) if value
            ):
                continue
            if left and right and (
                len(left) != len(right)
                or sum(a != b for a, b in zip(left, right)) > 2
            ):
                counts["rejected_diagnostic_pair"] += 1
                continue
            terminal_support = {"DIRECT_PAF_ENDPOINT", "TERMINAL_MOTIF"}
            if (
                left and payload.get("left_terminal_support") not in terminal_support
                and not grch38_coordinate_evidence
            ):
                continue
            if (
                right and payload.get("right_terminal_support") not in terminal_support
                and not grch38_coordinate_evidence
            ):
                continue
            if (
                left and right and evidence_structure == "PROVIRUS_OR_ARRAY"
                and payload.get("viral_body_boundary_owner") != "SUPPORTED"
                and not grch38_coordinate_evidence
            ):
                continue
            if not _direct_tsd_is_adjacent(payload, "5", left):
                continue
            if not _direct_tsd_is_adjacent(payload, "3", right):
                continue
            if left or right:
                if grch38_coordinate_evidence and (left, right) != (
                    "GATTGT", "GATTAT"
                ):
                    continue
                pairs.add((left, right))
        if len(pairs) == 1:
            left, right = pairs.pop()
            repaired.at[index, "5'_TSD"] = left or "NONE"
            repaired.at[index, "3'_TSD"] = right or "NONE"
            repaired.at[index, "TSD_Orientation"] = "INSERTION"
            exact_supported.append(index)
            counts["exact_supported"] += 1
            if used_tandem_base:
                counts["tandem_array_base_supported"] += 1
            if is_grch38_reference:
                counts["grch38_exact_coordinate_supported"] += 1
            counts[
                "exact_paired_supported"
                if left and right else "exact_one_sided_supported"
            ] += 1
        else:
            repaired.at[index, "5'_TSD"] = "NONE"
            repaired.at[index, "3'_TSD"] = "NONE"
            repaired.at[index, "TSD_Orientation"] = "INSERTION"
            counts["ambiguous_exact_signal" if len(pairs) > 1 else "no_exact_signal"] += 1

    normalized_present = normalize_present_day_tsd_semantics(
        repaired.loc[present_mask].copy()
    )
    repaired.loc[present_mask, TSD_PUBLISHED_COLUMNS] = normalized_present[
        TSD_PUBLISHED_COLUMNS
    ]
    if exact_supported:
        repaired.loc[exact_supported, "tsd_call_provenance"] += (
            ";RESOLVED_ADJACENT_DIRECT_JUNCTION"
        )
    no_signal = present_mask & ~repaired.index.isin(exact_supported)
    repaired.loc[no_signal, "tsd_call_provenance"] = (
        "NO_DIRECT_RESOLVED_ADJACENT_JUNCTION_EVIDENCE_CENSORED"
    )
    logging.info("Direct TSD evidence census: %s", counts)
    return repaired, counts


def update_existing_master_tsd(
    loci_df,
    results_dir,
    master_path,
    output_path=None,
    current_master_path=None,
    evidence_root=None,
    paired_authority_path=None,
    callability_authority_path=None,
    evidence_review_path=None,
):
    """Overlay current TSD fields while preserving every existing biological row.

    The acrocentric aggregate rows retain their aggregate ``Locus`` identity but
    obtain TSD fields from the exact current observation keyed by
    ``(orig_Locus, ID_Full)``.  Source and output paths may differ so a refreshed
    candidate never destroys the canonical input being reviewed.
    """
    output_path = output_path or master_path
    master = pd.read_csv(master_path, sep="\t", dtype=str, keep_default_na=False)
    source_columns = list(master.columns)
    required_master = {"Locus", "ID_Full", "ID", "Haplotype", "Structure"}
    missing_master = sorted(required_master - set(master.columns))
    if missing_master:
        raise ValueError(f"existing combined ORF table lacks identity fields: {missing_master}")
    if master["Locus"].eq("HML-2_8q24.3b").any():
        raise ValueError("existing master still contains retired 8q24.3b rows")
    master = master.copy()
    master["_tsd_source_row_order"] = range(len(master))

    if current_master_path:
        current = pd.read_csv(
            current_master_path, sep="\t", dtype=str, keep_default_na=False
        )
        required_current = {"Locus", "ID_Full", "ID", "Haplotype", "Structure", *TSD_PUBLISHED_COLUMNS}
        missing_current = sorted(required_current - set(current.columns))
        if missing_current:
            raise ValueError(
                f"current TSD master lacks required fields: {missing_current}"
            )
        current = current[~current["Locus"].eq("HML-2_8q24.3b")].copy()
        current_primary = set(current["Locus"])
        if len(current_primary) != 101 or "HML-2_8q24.3c" not in current_primary:
            raise ValueError(
                f"current TSD master is not the 101-locus domain: {len(current_primary)}"
            )
        no_orientation = current["tsd_call_status"].isin(
            {"NO_COPY_OR_UNRESOLVED", "NO_CURRENT_TSD_OBSERVATION"}
        )
        current.loc[no_orientation, "TSD_Orientation"] = "NOT_APPLICABLE"
        current.loc[~no_orientation, "TSD_Orientation"] = "INSERTION"
    else:
        current_parts = []
        for _, locus in loci_df.iterrows():
            locus_name = str(locus["locus_name"])
            type_str = VALID_LOCUS_TYPES[str(locus["type"]).strip().lower()]
            result_path = os.path.join(
                results_dir,
                locus_name,
                f"{locus_name}_orf_integrity_results_{type_str}_KCON.csv",
            )
            if not os.path.isfile(result_path):
                raise ValueError(f"missing current ORF/TSD result: {result_path}")
            current_part = pd.read_csv(result_path, dtype=str, keep_default_na=False)
            if (
                current_part.empty
                or "ID_Full" not in current_part.columns
                or "Structure" not in current_part.columns
            ):
                raise ValueError(f"unusable current ORF/TSD result: {result_path}")
            qc_state = current_part.get(
                "QC_State", pd.Series("", index=current_part.index)
            )
            current_part["presence_call"] = [
                classify_presence(structure, qc)
                for structure, qc in zip(current_part["Structure"], qc_state)
            ]
            current_part = normalize_present_day_tsd_semantics(current_part)
            current_part.insert(0, "Locus", locus_name)
            current_parts.append(current_part)
        current = pd.concat(current_parts, ignore_index=True, sort=False)
    if evidence_root:
        current, evidence_counts = apply_direct_tsd_evidence(current, evidence_root)
        # This assembly-gap representative requires its own exact junctions and
        # may never borrow the locus-level 7p22.1 pair.
        uncovered_7p_asmdup = (
            current["Locus"].eq("HML-2_7p22.1")
            & current["ID_Full"].eq(
                "HG00658_pat_hprc_r2_v1.0.1_HML-2_7p22.1_asmdup"
            )
        )
        current.loc[uncovered_7p_asmdup, "5'_TSD"] = "NONE"
        current.loc[uncovered_7p_asmdup, "3'_TSD"] = "NONE"
        current.loc[uncovered_7p_asmdup, "oriented_tsd_pair"] = "NONE/NONE"
        current.loc[uncovered_7p_asmdup, "tsd_length"] = ""
        current.loc[uncovered_7p_asmdup, "present_day_tsd_mismatch_count"] = ""
        current.loc[uncovered_7p_asmdup, "tsd_observable_state"] = "UNOBSERVABLE"
        current.loc[uncovered_7p_asmdup, "tsd_call_status"] = "NO_OBSERVABLE_TSD"
        current.loc[uncovered_7p_asmdup, "tsd_pair_admission"] = "UNAVAILABLE"
        current.loc[uncovered_7p_asmdup, "tsd_call_provenance"] = (
            "NO_DIRECT_EXACT_JUNCTION_EVIDENCE_CENSORED"
        )
        current.loc[uncovered_7p_asmdup, "tsd_evidence_role"] = (
            "NO_PAIRED_TSD_OBSERVATION"
        )
        current.loc[uncovered_7p_asmdup, "latent_insertion_tsd_relationship"] = (
            "IDENTICAL_COPIES_AT_INSERTION"
        )
        current.loc[uncovered_7p_asmdup, "genomic_forward_left_tsd"] = "NONE"
        current.loc[uncovered_7p_asmdup, "genomic_forward_right_tsd"] = "NONE"
        current.loc[uncovered_7p_asmdup, "TSD_Orientation"] = "INSERTION"
    c_current = current["Locus"].eq("HML-2_8q24.3c")
    c_reference_current = c_current & current["ID_Full"].str.contains(
        r"GCA_000001405\.15|GRCh38", case=False, regex=True, na=False
    )
    if int(c_reference_current.sum()) != 1:
        raise ValueError(
            "expected exactly one GRCh38 8q24.3c reference observation; "
            f"observed {int(c_reference_current.sum())}"
        )
    authority_fields = [
        "integrated_reference_build",
        "integrated_reference_interval",
        "integrated_sequence_authority",
        "orthologous_empty_site_authority",
    ]
    current.loc[c_current, authority_fields] = ""
    current.loc[c_reference_current, "integrated_reference_build"] = "GRCh38"
    current.loc[c_reference_current, "integrated_reference_interval"] = (
        "chr8:145021244-145028834"
    )
    current.loc[c_reference_current, "integrated_sequence_authority"] = (
        "GRCh38_PRESENT_INTEGRATED_ALLELE"
    )
    current.loc[c_reference_current, "orthologous_empty_site_authority"] = (
        "T2T_CHM13_PRE_INTEGRATION_EMPTY_SITE_ONLY"
    )
    reference_pair = current.loc[
        c_reference_current, ["5'_TSD", "3'_TSD", "tsd_pair_admission"]
    ].iloc[0].to_dict()
    if reference_pair != {
        "5'_TSD": "GATTGT",
        "3'_TSD": "GATTAT",
        "tsd_pair_admission": "PAIRED_ACCEPTED",
    }:
        raise ValueError(
            "8q24.3c GRCh38 TSD lacks authenticated exact-coordinate evidence: "
            f"{reference_pair}"
        )
    current_keys = ["Locus", "ID_Full"]
    duplicated = current.duplicated(current_keys, keep=False)
    if duplicated.any():
        rows = current.loc[duplicated, current_keys].head(10).to_dict("records")
        raise ValueError(f"duplicate current ORF/TSD identities: {rows}")

    master_base = master.drop(
        columns=[column for column in TSD_MERGE_COLUMNS if column in master.columns],
        errors="ignore",
    )
    # Current observations use their physical 101-locus identity.  Historical
    # acrocentric aggregate rows are deliberate aliases of those same retained
    # biological observations, and every one carries the exact source locus in
    # orig_Locus.  Match those aliases to the current call without rewriting
    # their aggregate identity or any ORF/sequence field.
    aggregate = master_base["Locus"].isin({"HML-2_acro_type1", "HML-2_acro_type2"})
    if aggregate.any() and "orig_Locus" not in master_base.columns:
        raise ValueError("acrocentric aggregate rows require orig_Locus for exact TSD mapping")
    master_base["_tsd_lookup_locus"] = master_base["Locus"]
    master_base.loc[aggregate, "_tsd_lookup_locus"] = master_base.loc[aggregate, "orig_Locus"]
    if master_base.loc[aggregate, "_tsd_lookup_locus"].str.strip().eq("").any():
        raise ValueError("acrocentric aggregate row has empty orig_Locus")

    updates = current[current_keys + TSD_MERGE_COLUMNS].copy().rename(
        columns={"Locus": "_tsd_lookup_locus"}
    )
    refreshed = master_base.merge(
        updates,
        on=["_tsd_lookup_locus", "ID_Full"],
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    unmatched = refreshed["_merge"].eq("left_only")
    refreshed.loc[unmatched, "5'_TSD"] = "NONE"
    refreshed.loc[unmatched, "3'_TSD"] = "NONE"
    refreshed.loc[unmatched, "oriented_tsd_pair"] = "NONE/NONE"
    refreshed.loc[unmatched, "tsd_observable_state"] = "UNOBSERVABLE"
    refreshed.loc[unmatched, "tsd_call_status"] = "NO_CURRENT_TSD_OBSERVATION"
    refreshed.loc[unmatched, "tsd_pair_admission"] = "UNAVAILABLE"
    refreshed.loc[unmatched, "tsd_call_provenance"] = "CURRENT_101_LOCUS_RERUN_NO_MATCH"
    refreshed.loc[unmatched, "tsd_evidence_role"] = "NO_PAIRED_TSD_OBSERVATION"
    refreshed.loc[unmatched, "latent_insertion_tsd_relationship"] = "NO_TSD_OBSERVATION"
    refreshed.loc[unmatched, "genomic_forward_left_tsd"] = "NONE"
    refreshed.loc[unmatched, "genomic_forward_right_tsd"] = "NONE"
    refreshed.loc[unmatched, "TSD_Orientation"] = "NOT_APPLICABLE"
    refreshed.drop(columns=["_merge", "_tsd_lookup_locus"], inplace=True)

    # These six legacy master rows were mislabeled Solo-LTR/unknown-strand. The
    # exact current ORF rows classify each retained source as Fragment on the
    # reverse genomic strand. Their two terminal alignments face in opposite
    # directions: AGGT and CTGT sit outside two inverted LTR 5' faces, not at
    # the 5'/3' boundaries of one insertion. Publish a definitive no-call while
    # repairing the structural fields; this is not a quarantine.
    resolved_current = (
        current["Locus"].eq("HML-2_3q12.3")
        & current["ID_Full"].isin(DEFINITIVE_NO_TSD_3Q12_3_FRAGMENT_IDS)
    )
    resolved_master = (
        refreshed["Locus"].eq("HML-2_3q12.3")
        & refreshed["ID_Full"].isin(DEFINITIVE_NO_TSD_3Q12_3_FRAGMENT_IDS)
    )
    if int(resolved_current.sum()) != 6 or int(resolved_master.sum()) != 6:
        raise ValueError("expected all six resolved 3q12.3 fragment observations")
    current_resolution = current.loc[resolved_current]
    if not (
        current_resolution["Structure"].eq("Fragment").all()
        # Present-day normalization represents the current ORF row's raw '-'
        # strand as REVERSE before the TSD merge.
        and current_resolution["Strand"].eq("REVERSE").all()
    ):
        raise ValueError("current ORF source does not resolve all six 3q12.3 fragments")
    refreshed.loc[resolved_master, "Structure"] = "Fragment"
    refreshed.loc[resolved_master, "Strand"] = "-"
    refreshed.loc[resolved_master, "5'_TSD"] = "NONE"
    refreshed.loc[resolved_master, "3'_TSD"] = "NONE"
    refreshed.loc[resolved_master, "oriented_tsd_pair"] = "NONE/NONE"
    refreshed.loc[resolved_master, "tsd_length"] = ""
    refreshed.loc[resolved_master, "present_day_tsd_mismatch_count"] = ""
    refreshed.loc[resolved_master, "tsd_observable_state"] = "UNOBSERVABLE"
    refreshed.loc[resolved_master, "tsd_call_status"] = "NO_OBSERVABLE_TSD"
    refreshed.loc[resolved_master, "tsd_pair_admission"] = "UNAVAILABLE"
    refreshed.loc[resolved_master, "tsd_call_provenance"] = (
        "INVERTED_LTR_OUTER_FACES_NO_PAIRED_INSERTION_BOUNDARY"
    )
    refreshed.loc[resolved_master, "tsd_evidence_role"] = (
        "NO_PAIRED_TSD_OBSERVATION"
    )
    refreshed.loc[resolved_master, "genomic_forward_left_tsd"] = "NONE"
    refreshed.loc[resolved_master, "genomic_forward_right_tsd"] = "NONE"

    terminal_authorities = (
        paired_authority_path,
        callability_authority_path,
        evidence_review_path,
    )
    if not all(terminal_authorities):
        raise ValueError(
            "terminal TSD materialization requires paired, ORF-callability, "
            "and direct-junction review authorities"
        )
    paired_authority = pd.read_csv(
        paired_authority_path, sep="\t", dtype=str, keep_default_na=False
    )
    callability_authority = pd.read_csv(
        callability_authority_path, sep="\t", dtype=str, keep_default_na=False
    )
    evidence_review = pd.read_csv(
        evidence_review_path, sep="\t", dtype=str, keep_default_na=False
    )
    refreshed, terminal_state_report = terminalize_published_tsd_states(
        refreshed,
        paired_authority,
        callability_authority,
        evidence_review,
    )
    logging.info("Terminal TSD/QC state census: %s", terminal_state_report)

    existing_keys = set(zip(master_base["Locus"], master_base["ID_Full"]))
    is_new_current = pd.Series(
        [key not in existing_keys for key in zip(current["Locus"], current["ID_Full"])],
        index=current.index,
    )
    new_current = current.loc[is_new_current].copy()
    if not new_current.empty:
        logging.info(
            "TSD-only refresh ignored %d current observations absent from the "
            "fixed 63,750-row master domain",
            len(new_current),
        )

    loci = set(refreshed["Locus"])
    primary_loci = loci - {"HML-2_acro_type1", "HML-2_acro_type2"}
    if (
        len(primary_loci) != 101
        or "HML-2_8q24.3b" in loci
        or "HML-2_8q24.3c" not in primary_loci
    ):
        raise ValueError(
            "wrong TSD-refreshed domain: "
            f"primary_loci={len(primary_loci)} all_locus_labels={len(loci)}"
        )
    source_aggregate_rows = int(aggregate.sum())
    refreshed_aggregate_rows = int(
        refreshed["Locus"].isin({"HML-2_acro_type1", "HML-2_acro_type2"}).sum()
    )
    if refreshed_aggregate_rows != source_aggregate_rows:
        raise ValueError(
            "acrocentric aggregate rows were not preserved exactly: "
            f"source={source_aggregate_rows} refreshed={refreshed_aggregate_rows}"
        )
    reference = refreshed[
        refreshed["Locus"].eq("HML-2_8q24.3c")
        & refreshed["ID_Full"].str.contains(
            r"GCA_000001405\.15|GRCh38", case=False, regex=True
        )
    ]
    if len(reference) != 1:
        raise ValueError(f"expected one 8q24.3c GRCh38 row, observed {len(reference)}")
    reference = reference.iloc[0]
    expected = {
        "5'_TSD": "GATTGT",
        "3'_TSD": "GATTAT",
        "TSD_Orientation": "INSERTION",
        "present_day_tsd_mismatch_count": "1",
        "tsd_observable_state": "TWO_SIDED",
        "tsd_pair_admission": "PAIRED_ACCEPTED",
        "integrated_reference_interval": "chr8:145021244-145028834",
    }
    wrong = {
        field: (reference.get(field, ""), wanted)
        for field, wanted in expected.items()
        if reference.get(field, "") != wanted
    }
    if wrong:
        raise ValueError(f"8q24.3c reference TSD semantics incorrect: {wrong}")

    if len(refreshed) != len(master):
        raise ValueError(
            f"TSD-only refresh changed row count: {len(master)} -> {len(refreshed)}"
        )
    refreshed.sort_values("_tsd_source_row_order", inplace=True)
    refreshed.drop(columns=["_tsd_source_row_order"], inplace=True)
    refreshed = refreshed.reindex(columns=source_columns)
    atomic_write_dataframe(refreshed, output_path, sep="\t", index=False, na_rep="NA")
    return refreshed, int((~unmatched).sum()), len(new_current), refreshed_aggregate_rows


def setup_logging(log_file):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )

def map_generic_hap(h):
    h = str(h).lower()
    if h in ['h1', 'hap1', '1', 'mat', 'maternal']: return '1'
    if h in ['h2', 'hap2', '2', 'pat', 'paternal']: return '2'
    return h

def normalize_id(val):
    s = str(val).strip()
    if 'GCA_' in s: return 'GCA'
    if 'chm13' in s.lower(): return 'chm13v2.0'
    m = re.match(r'([hH][gG]|[nN][aA])([0-9]+)', s)
    if m: return f"{m.group(1).upper()}{m.group(2)}"
    return s.split('_')[0].split('#')[0]

def load_universe(file_path):
    universe = set()
    try:
        df = pd.read_csv(file_path, sep='\t')
        for _, row in df.iterrows():
            clean_id = normalize_id(str(row['ID']).strip())
            clean_hap = str(row['Haplotype']).strip()
            universe.add((clean_id, clean_hap))
        logging.info(f"Loaded {len(universe)} unique pairs from Universe file.")
        return universe
    except Exception as e:
        logging.error(f"FATAL: Could not load universe file: {e}")
        sys.exit(1)

def load_hgsvc3_haplotypes(dir_path):
    haplotypes = set()
    try:
        if not os.path.exists(dir_path): return set()
        for filename in os.listdir(dir_path):
            if filename.endswith((".fa", ".fasta")):
                parts = filename.split('_')
                if len(parts) >= 2 and parts[0].upper().startswith(("HG", "NA")):
                    sample_id = parts[0].upper()
                    hap = parts[1].replace('.fa', '').replace('.fasta', '')
                    if hap in ['h1', 'h2']:
                        haplotypes.add((sample_id, hap))
        return haplotypes
    except:
        return set()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--locus-file', required=True)
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--universe-file', required=True)
    parser.add_argument('--hgsvc3-sample-dir', required=False, default="")
    parser.add_argument('--copy-state-ledger', required=False, default="",
                        help='Exact observation-level latent-copy weights and marginalized ORF expectations')
    parser.add_argument(
        '--tsd-update-existing-master', default="",
        help='Update only TSD fields in this existing combined ORF TSV; no copy-state rebuild',
    )
    parser.add_argument(
        '--tsd-output-master', default="",
        help='Write the TSD-refreshed candidate here; defaults to the source path',
    )
    parser.add_argument(
        '--tsd-current-master', default="",
        help='Use this already-refreshed 101-locus table as the exact TSD/current-row overlay source',
    )
    parser.add_argument(
        '--tsd-direct-evidence-root', default="",
        help='Require resolved boundary-adjacent evidence for every present-row TSD call',
    )
    parser.add_argument(
        '--tsd-paired-authority', default="",
        help='Complete paired TSD observation authority used for terminal presence classes',
    )
    parser.add_argument(
        '--tsd-callability-authority', default="",
        help='Exact ORF callability authority for noncarrier and recoverable technical rows',
    )
    parser.add_argument(
        '--tsd-evidence-review', default="",
        help='Exact retained direct-junction review authority for sequence-present no-calls',
    )
    
    args = parser.parse_args()

    log_file = os.path.join(args.results_dir, "validation_and_completion.log")
    setup_logging(log_file)
    
    logging.info("--- Starting Validation and Completion ---")
    
    try:
        loci_df = pd.read_csv(args.locus_file, sep='\t', header=None, names=['locus_name', 'chr', 'start', 'end', 'type', 'padding'])
    except FileNotFoundError:
        logging.error(f"FATAL: Locus file not found")
        sys.exit(1)

    # --- Fail closed on an unknown locus type ---------------------------------------------------
    # A typo such as "Typel" for "TypeI" must never be absorbed by an else-branch and silently
    # analysed against the Type II reference. The WHOLE table is validated here, before any output
    # is written, so an invalid generation cannot replace the previously published master.
    invalid_types = [
        f"{row['locus_name']} (type={row['type']!r})"
        for _, row in loci_df.iterrows()
        if str(row['type']).strip().lower() not in VALID_LOCUS_TYPES
    ]
    if invalid_types:
        logging.error(
            "FATAL: invalid locus type(s) in %s -- refusing to run. Valid types are TypeI/TypeII; "
            "any previously published master is left untouched: %s",
            args.locus_file, "; ".join(invalid_types),
        )
        sys.exit(1)

    if args.tsd_update_existing_master:
        try:
            refreshed, matched_rows, added_rows, aggregate_rows = update_existing_master_tsd(
                loci_df,
                args.results_dir,
                args.tsd_update_existing_master,
                args.tsd_output_master or args.tsd_update_existing_master,
                args.tsd_current_master or None,
                args.tsd_direct_evidence_root or None,
                args.tsd_paired_authority or None,
                args.tsd_callability_authority or None,
                args.tsd_evidence_review or None,
            )
        except Exception as exc:
            logging.error(f"FATAL: TSD update failed: {exc}")
            sys.exit(1)
        logging.info(
            "TSD update complete: rows=%d locus_labels=%d matched_current_rows=%d "
            "ignored_out_of_master_current_rows=%d "
            "preserved_acrocentric_aggregate_rows=%d output=%s",
            len(refreshed), refreshed["Locus"].nunique(), matched_rows, added_rows,
            aggregate_rows, args.tsd_output_master or args.tsd_update_existing_master,
        )
        return

    if not args.copy_state_ledger:
        logging.error("FATAL: --copy-state-ledger is required for a full biological combine")
        sys.exit(1)
    master_haplotypes = load_universe(args.universe_file)
    try:
        copy_state_ledger = load_copy_state_ledger(args.copy_state_ledger)
    except Exception as exc:
        logging.error(f"FATAL: invalid copy-state ledger: {exc}")
        sys.exit(1)
    used_copy_state_keys = set()

    hgsvc3_haplotypes = set()
    if args.hgsvc3_sample_dir:
        hgsvc3_haplotypes = load_hgsvc3_haplotypes(args.hgsvc3_sample_dir)

    generic_hgsvc3 = {(uid, map_generic_hap(hap)) for uid, hap in hgsvc3_haplotypes}

    total_placeholders_added = 0
    all_completed_dfs = []

    for index, locus in loci_df.iterrows():
        original_locus_name = locus['locus_name']
        locus_name = (
            original_locus_name[:-5]
            if original_locus_name.endswith("_hg38")
            else original_locus_name
        )
        
        if original_locus_name == 'HML2_7p22.1b': continue
        if original_locus_name == 'HML2_7p22.1a': locus_name = 'HML2_7p22.1'
        
        type_str = VALID_LOCUS_TYPES[str(locus['type']).strip().lower()]   # validated above
        kcon_ref_id = f"{type_str}_KCON"
        accessory_orf = 'np9' if type_str == 'type1' else 'rec'
        
        core_orfs = ["gag", "pro", "pol", "env", accessory_orf]
        # presence_call (biological non-carrier vs technical missing vs present) and the
        # array-aware copy-number fields are surfaced up front; the physical-Pol-region vs
        # translated-Pol-compatibility fields are surfaced immediately after the physical `pol`
        # block. All are additive and must survive into the master combined TSV.
        ordered_cols = ["ID_Full", "ID", "Haplotype", "Source_Identifier", "Structure", "presence_call",
                        "observation_state",
                        "QC_State", "QC_Reason", "Technical_Error_Reason",
                        "copy_class", "independent_copy", "copy_observation_class", "observation_role",
                        "observation_sequence_bytes", "observation_sequence_sha256", "artifact_evidence_retained",
                        "copy_state_model_status", "observation_weight", *COPY_STATE_WEIGHT_COLS,
                        "expected_biological_copy_count", *COPY_EXPECTATION_COLS,
                        "copy_state_weight_provenance", "artifact_evidence_sha256", "artifact_evidence_bytes",
                        "duplication_birth_constraint", "divergence_marginalization_model",
                        "consequence_marginalization_provenance", "5'_TSD", "3'_TSD",
                        "oriented_tsd_pair", "tsd_length",
                        "present_day_tsd_mismatch_count", "tsd_observable_state",
                        "tsd_call_status", "tsd_pair_admission",
                        "tsd_call_provenance", "tsd_source_interval",
                        "tsd_evidence_role", "latent_insertion_tsd_relationship",
                        "genomic_forward_left_tsd", "genomic_forward_right_tsd",
                        "integrated_reference_build", "integrated_reference_interval",
                        "integrated_sequence_authority",
                        "orthologous_empty_site_authority"]
        for orf in core_orfs:
            ordered_cols.extend([
                orf, f"missense_{orf}", f"{orf}_dna_identity",
                f"{orf}_protein_identity", f"{orf}_coverage"
            ])
            if orf == "pol":
                ordered_cols.extend([
                    "pol_region_present", "pol_translation_compatible", "pol_translation_basis",
                    "pol_qc_flag"
                ])
        ordered_cols.append("all_orfs_intact")
        
        results_csv_path = os.path.join(args.results_dir, original_locus_name, f"{original_locus_name}_orf_integrity_results_{kcon_ref_id}.csv")
        completed_csv_path = os.path.join(args.results_dir, original_locus_name, f"{locus_name}_orf_integrity_results_{kcon_ref_id}_completed.csv")
        
        logging.info(f"--- Processing Locus: {original_locus_name} ---")
        
        # --- ABSENCE SEMANTICS: biological non-carrier  vs  technical missingness --------------------
        # An expected sample/haplotype with no provirus row is NOT the same as one whose extraction
        # failed, and NEITHER is a reason to drop the locus:
        #   * biological absence  -> the haplotype was assayed and carries NO insertion here (a true
        #                            non-carrier). This is real, informative data and is RETAINED --
        #                            an entirely absent locus is a result, not a missing one.
        #   * technical missing   -> the extraction could not be evaluated (missing/empty/unreadable
        #                            result, BAM/index error, graph or ORF failure). Absence here is
        #                            UNKNOWN, is never counted as a confident absence, and is emitted
        #                            explicitly with its cause so the gap is visible in the master.
        # No locus is ever silently excluded, and no result is ever carried over from an older run:
        # if this generation has no readable result for a locus, the locus is published as technical
        # errors for this generation, overwriting any previous *_completed.csv.
        results_df = pd.DataFrame()
        locus_technical_reason = None
        if not os.path.exists(results_csv_path):
            locus_technical_reason = "MISSING_UPSTREAM_RESULT"
        else:
            try:
                results_df = pd.read_csv(results_csv_path)
            except Exception as e:
                logging.error(f"Could not read {results_csv_path}: {e}")
                locus_technical_reason = "UNREADABLE_UPSTREAM_RESULT"
            else:
                if results_df.empty:
                    # A results file with no rows is not evidence of absence: the ORF stage emits a
                    # row for every haplotype it interrogated, including non-carriers.
                    locus_technical_reason = "EMPTY_UPSTREAM_RESULT"
                elif 'Structure' not in results_df.columns:
                    locus_technical_reason = "MALFORMED_UPSTREAM_RESULT"
        if locus_technical_reason:
            results_df = pd.DataFrame()
            logging.error(f"[{original_locus_name}] {locus_technical_reason}: emitting an explicit "
                          f"{TECHNICAL_ERROR} row for every expected sample/haplotype. This locus was "
                          f"NOT interrogated; these rows are not biological absence.")

        if not results_df.empty and 'ID' in results_df.columns:
            results_df['ID_Clean'] = results_df['ID'].apply(normalize_id)

            def normalize_haplotype(row):
                # Haplotype-resolution fallback chain, most authoritative source first. Each step is
                # a POSITIVE identification; the terminal default is reached only when nothing at all
                # identifies the haplotype, and it is logged rather than silently applied so a parsing
                # regression cannot masquerade as real h1 data (and collide with a genuine h1 row).
                raw_id = str(row['ID']).lower()
                id_full = str(row.get('ID_Full', '')).lower()
                if '.mat' in raw_id or '_mat' in raw_id or '_mat' in id_full: return 'mat'
                if '.pat' in raw_id or '_pat' in raw_id or '_pat' in id_full: return 'pat'
                if 'chm13' in raw_id or 'chm13' in id_full: return 'h1'   # haploid reference -> single hap
                if 'gca' in raw_id or 'gca' in id_full: return 'h1'       # GRCh38 primary -> single hap

                existing = str(row.get('Haplotype', '')).lower()
                if existing in ['pat', 'mat']: return existing
                if existing in ['h1', '1', 'hap1']: return 'h1'
                if existing in ['h2', '2', 'hap2']: return 'h2'

                # Last resort: the Source_Identifier often carries the graph path haplotype (…#1#… / …#2#…).
                src = str(row.get('Source_Identifier', '')).lower()
                m = re.search(r'#(\d)#', src)
                if m:
                    return 'h1' if m.group(1) == '1' else ('h2' if m.group(1) == '2' else 'UNKNOWN')
                for tok in ('_h1_', '_h2_', '.h1', '.h2'):
                    if tok in id_full:
                        return 'h1' if '1' in tok else 'h2'

                # Genuinely unidentifiable: preserve UNKNOWN. Defaulting to h1
                # corrupts exact sample-copy identity and can collide with a real h1 row.
                logging.warning(f"Haplotype unresolved for ID='{row.get('ID')}' "
                                f"ID_Full='{row.get('ID_Full')}' Source='{row.get('Source_Identifier')}'; "
                                f"preserving explicit UNKNOWN copy state.")
                return 'UNKNOWN'

            results_df['Haplotype_Clean'] = results_df.apply(normalize_haplotype, axis=1)
            results_df['Haplotype'] = results_df['Haplotype_Clean']

            jargon_map = {
                'ref_map_error': 'Mapping_Error',
                'invalid_chars': 'Invalid_Characters',
                'too_short': 'Sequence_Too_Short',
                'too_short_after_trim': 'Sequence_Too_Short',
                'translation_error': 'Translation_Error',
                'ComparisonError': 'Reference_Missing',
                'Deletion': 'Protein_Missing', 
                'AlignmentFailed': 'Alignment_Failed',
                'nonsense': 'Nonsense',
                'intact': 'Intact',
                'no_stop': 'No_Stop_Codon'
            }

            cols_to_clean = [c for c in results_df.columns if c in ['gag', 'pro', 'pol', 'env', 'np9', 'rec'] or c.startswith('missense_')]
            for col in cols_to_clean:
                results_df[col] = results_df[col].replace(jargon_map)
            
        else:
            results_df['ID_Clean'] = []
            results_df['Haplotype_Clean'] = []

        # Tag every REAL row with an explicit presence_call so downstream never has to guess whether a
        # non-"present" Structure is a real non-carrier or a technical failure, then canonicalize the
        # two non-carriage states onto the single spelling each. Upstream stages report absence and
        # failure under several historical names; they are collapsed HERE, once, and only within
        # their own class -- a technical failure keeps its cause and never becomes absence.
        if not results_df.empty:
            _blank = pd.Series([""] * len(results_df), index=results_df.index)
            _qc_state = results_df['QC_State'] if 'QC_State' in results_df.columns else _blank
            _qc_reason = results_df['QC_Reason'] if 'QC_Reason' in results_df.columns else _blank
            results_df['presence_call'] = [
                classify_presence(s, q) for s, q in zip(results_df['Structure'], _qc_state)
            ]
            _is_technical = results_df['presence_call'].eq("technical_missing")
            results_df['Technical_Error_Reason'] = [
                technical_reason(s, r) if t else ""
                for s, r, t in zip(results_df['Structure'], _qc_reason, _is_technical)
            ]
            results_df.loc[_is_technical, 'Structure'] = TECHNICAL_ERROR
            results_df.loc[results_df['presence_call'].eq("absent_noncarrier"), 'Structure'] = BIOLOGICAL_ABSENCE
            results_df['observation_state'] = results_df['presence_call'].map({
                'present': 'PRESENT',
                'absent_noncarrier': 'EMPTY_SITE_SUPPORTED',
                'technical_missing': 'TECHNICAL_MISSING',
            })

            # Candidate rows are observations, not automatic copies. Replace the
            # ORF-stage UNMEASURED placeholders only from the external latent model.
            # Every currently observed row must have weights.  Extra ledger rows are
            # not a biological error: a corrected extraction can turn a formerly
            # false-positive candidate into an empty-site call (8q24.3c previously
            # carried one such T2T-routed row for nearly every haplotype).  Ignore
            # those obsolete rows instead of using them or rejecting current data.
            # Artifact-weighted observations
            # retain and verify their original sequence bytes/digest. True-dup
            # states require identical-at-birth coupling and an explicit later-
            # divergence marginalization model (validated at ledger load).
            present_mask = results_df['presence_call'].eq("present")
            present_ids = set(results_df.loc[present_mask, 'ID_Full'].astype(str))
            locus_ledger = copy_state_ledger[copy_state_ledger['locus'] == str(original_locus_name)].copy()
            ledger_ids = set(locus_ledger['ID_Full'].astype(str))
            missing_ids = present_ids - ledger_ids
            if missing_ids:
                logging.error("FATAL [%s]: currently observed candidates lack copy-state weights: missing=%s",
                              original_locus_name, sorted(missing_ids))
                sys.exit(1)
            locus_ledger = locus_ledger[
                locus_ledger['ID_Full'].astype(str).isin(present_ids)
            ].copy()
            model_cols = [c for c in COPY_STATE_LEDGER_COLUMNS if c not in ("locus", "ID_Full")]
            results_df.drop(columns=[c for c in model_cols if c in results_df.columns], inplace=True, errors='ignore')
            results_df = results_df.merge(locus_ledger.drop(columns=['locus']), on='ID_Full', how='left', validate='one_to_one')
            observed = results_df['presence_call'].eq("present")
            if results_df.loc[observed, 'observation_weight'].isna().any():
                logging.error("FATAL [%s]: a present observation lacks latent-copy weights", original_locus_name)
                sys.exit(1)
            group_sums = results_df.loc[observed].groupby(['ID_Clean', 'Haplotype_Clean'])['observation_weight'].sum()
            if (group_sums - 1.0).abs().gt(1e-9).any():
                logging.error("FATAL [%s]: observation weights do not normalize to 1 within sample/haplotype: %s",
                              original_locus_name, group_sums.to_dict())
                sys.exit(1)
            results_df.loc[observed, 'copy_state_model_status'] = "WEIGHTED_MARGINALIZED"
            results_df.loc[observed, 'copy_class'] = "latent_state_mixture"
            results_df.loc[observed, 'independent_copy'] = "MARGINALIZED_NOT_ROW_COUNTED"
            used_copy_state_keys.update((str(original_locus_name), ident) for ident in present_ids)

            # A real call for a key makes that key's absence row redundant: the site cannot be both
            # occupied and empty, so drop the placeholder. Only the redundant ABSENCE row goes --
            # genuine _part1.._partN tandem copies are separate PHYSICAL rows and are all kept, and a
            # technical-failure row is never dropped because its cause must stay visible.
            _present_keys = set(zip(
                results_df.loc[results_df['presence_call'].eq("present"), 'ID_Clean'],
                results_df.loc[results_df['presence_call'].eq("present"), 'Haplotype_Clean'],
            ))
            if _present_keys:
                _redundant = results_df['presence_call'].eq("absent_noncarrier") & pd.Series(
                    [(i, h) in _present_keys
                     for i, h in zip(results_df['ID_Clean'], results_df['Haplotype_Clean'])],
                    index=results_df.index,
                )
                if _redundant.any():
                    logging.info(f"[{original_locus_name}] Dropped {int(_redundant.sum())} redundant "
                                 f"absence row(s) for sample/haplotypes that have a real call.")
                    results_df = results_df[~_redundant]

        # Samples that had a TECHNICAL failure at this locus: their missing haplotypes must be labelled
        # technical_missing (unknown), NOT confidently imputed as non-carriers.
        technical_samples = set()
        if not results_df.empty:
            technical_samples = set(
                results_df.loc[results_df['presence_call'].eq("technical_missing"), 'ID_Clean']
            )

        if not results_df.empty:
            samples_with_assembly = set(
                results_df[results_df['Haplotype_Clean'].isin(['pat', 'mat'])]['ID_Clean']
            )
            if samples_with_assembly:
                condition = ~ (
                    (results_df['ID_Clean'].isin(samples_with_assembly)) &
                    (results_df['Haplotype_Clean'].isin(['h1', 'h2']))
                )
                results_df = results_df[condition]

        processed_set = set(zip(results_df['ID_Clean'], results_df['Haplotype_Clean'])) if not results_df.empty else set()
        
        generic_master = {}
        for uid, hap in master_haplotypes:
            generic_master[(uid, map_generic_hap(hap))] = hap
            
        generic_processed = set()
        for uid, hap in processed_set:
            generic_processed.add((uid, map_generic_hap(hap)))

        missing_generic = set(generic_master.keys()) - generic_processed
        missing_haplotypes = [(uid, generic_master[(uid, ghap)]) for uid, ghap in missing_generic]

        placeholders = []
        if missing_haplotypes:
            logging.info(f"Adding {len(missing_haplotypes)} placeholders.")
            for sample_id, haplotype in missing_haplotypes:
                is_hgsvc3 = (sample_id, map_generic_hap(haplotype)) in generic_hgsvc3

                # A haplotype with no row is a confident biological non-carrier ONLY if it was really
                # interrogated. If this locus produced no readable result at all, or this sample failed
                # technically at this locus, its absence is UNKNOWN -- never impute it as a non-carrier.
                if locus_technical_reason:
                    ph_structure, ph_presence, ph_reason = TECHNICAL_ERROR, "technical_missing", locus_technical_reason
                elif sample_id in technical_samples:
                    ph_structure, ph_presence, ph_reason = TECHNICAL_ERROR, "technical_missing", "UPSTREAM_TECHNICAL_FAILURE_AT_LOCUS"
                else:
                    ph_structure, ph_presence, ph_reason = BIOLOGICAL_ABSENCE, "absent_noncarrier", ""

                if sample_id == 'GCA':
                     id_full = "GCA_GRCh38"
                     src_id = "Validator_Script"
                elif is_hgsvc3:
                     hap_num = '1' if map_generic_hap(haplotype) == '1' else '2'
                     id_full = f"{sample_id}_{haplotype}_{locus_name}_graph_hgsvc3"
                     src_id = f"{sample_id}#{hap_num}#haplotype{hap_num}"
                else:
                     suffix = "technical_error" if ph_presence == "technical_missing" else "placeholder"
                     id_full = f"{sample_id}_{haplotype}_{locus_name}_{suffix}"
                     src_id = "Validator_Script"

                placeholders.append({
                    'ID': sample_id,
                    'ID_Full': id_full,
                    'Haplotype': haplotype,
                    'Source_Identifier': src_id,
                    'Structure': ph_structure,
                    'presence_call': ph_presence,
                    'observation_state': (
                        "EMPTY_SITE_SUPPORTED"
                        if ph_presence == "absent_noncarrier"
                        else "TECHNICAL_MISSING"
                    ),
                    'QC_State': TECHNICAL_ERROR if ph_presence == "technical_missing" else "PASS_EMPTY_SITE",
                    'QC_Reason': ph_reason,
                    'Technical_Error_Reason': ph_reason,
                })

        if placeholders:
            combined_df = pd.concat([results_df, pd.DataFrame(placeholders)], ignore_index=True)
            total_placeholders_added += len(missing_haplotypes)
        else:
            combined_df = results_df.copy()

        if not combined_df.empty:
            for col in ordered_cols:
                if col not in combined_df.columns:
                    combined_df[col] = "NA"
                    
            combined_df.drop(columns=['ID_Clean', 'Haplotype_Clean'], inplace=True, errors='ignore')
            
            final_cols = ordered_cols + [c for c in combined_df.columns if c not in ordered_cols]
            combined_df = combined_df[final_cols]
            
            # Add Locus identifier for the master combined TSV
            combined_df.insert(0, 'Locus', locus_name)
            
            atomic_write_dataframe(combined_df, completed_csv_path, sep=",", index=False, na_rep='NA')
            all_completed_dfs.append(combined_df)

    # --- FINAL COMBINATION ---
    logging.info("--- Combining all loci into master TSV ---")
    if all_completed_dfs:
        master_df = pd.concat(all_completed_dfs, ignore_index=True)

        # 8q24.3c is a reverse-orientation integrated allele in GRCh38.  CHM13
        # is the orthologous pre-integration empty site and must never replace
        # the integrated sequence authority.  The reference row is a genuine
        # GRCh38 extraction/ORF observation; its two present-day TSD copies were
        # read directly from the genomic flanks at the coordinates below.  They
        # differ by one current base, while the insertion-time latent state is
        # still identical copies.  Do not infer these values for sample rows.
        c_mask = master_df['Locus'].eq('HML-2_8q24.3c')
        if c_mask.any():
            c_reference = (
                c_mask
                & master_df['ID_Full'].astype(str).str.contains(
                    r'GCA_000001405\.15|GRCh38', case=False, regex=True, na=False
                )
                & master_df['presence_call'].eq('present')
            )
            if int(c_reference.sum()) != 1:
                raise ValueError(
                    '8q24.3c requires exactly one present GRCh38 reference '
                    f'ORF row; observed {int(c_reference.sum())}'
                )
            authority_fields = [
                'integrated_reference_build',
                'integrated_reference_interval',
                'integrated_sequence_authority',
                'orthologous_empty_site_authority',
            ]
            master_df.loc[c_mask, authority_fields] = ''
            master_df.loc[c_reference, 'integrated_reference_build'] = 'GRCh38'
            master_df.loc[c_reference, 'integrated_reference_interval'] = (
                'chr8:145021244-145028834'
            )
            master_df.loc[c_reference, 'integrated_sequence_authority'] = (
                'GRCh38_PRESENT_INTEGRATED_ALLELE'
            )
            master_df.loc[c_reference, 'orthologous_empty_site_authority'] = (
                'T2T_CHM13_PRE_INTEGRATION_EMPTY_SITE_ONLY'
            )
            master_df.loc[c_reference, 'Strand'] = 'REVERSE'
            master_df.loc[c_reference, "5'_TSD"] = 'ATAATC'
            master_df.loc[c_reference, "3'_TSD"] = 'ACAATC'
            master_df.loc[c_reference, 'TSD_Orientation'] = 'GENOMIC_FORWARD'
            master_df.loc[c_reference, 'tsd_source_interval'] = (
                'chr8:145021238-145021243;chr8:145028836-145028841'
            )
        master_df = normalize_present_day_tsd_semantics(master_df)
        master_tsv_path = os.path.join(args.results_dir, "combined_hml2_orf_analysis.tsv")

        # --- Bodyless tandem-unit QC FLAG -----------------------------------------------------
        # At a dead solo-LTR locus that happens to carry a tandem ARRAY of solo-LTRs (e.g. 1p31.1b),
        # the shared-LTR splitter can emit "Provirus_from_Multi" units BETWEEN the solo-LTRs that
        # have no provirus body -- pol Protein_Missing at ~0 coverage and only a pulled-in env/np9
        # fragment. That is worth recording: a genuine tandem-array unit has a real pol (the 539
        # real _part units all score pol Intact), so such a unit is suspect.
        #
        # It is NOT grounds for rewriting Structure. Pol damage is a QC observation about ONE
        # segment of a copy, not evidence that the physical tandem unit is absent -- a damaged Pol
        # cannot erase or reclassify a unit that the splitter physically resolved. Rewriting the
        # Structure here would silently delete a real physical copy from the array and corrupt the
        # per-locus copy number. So the row keeps its Structure and identity, and is FLAGGED for
        # downstream filtering to apply (or not) with the biology in view.
        _polcov = pd.to_numeric(master_df.get("pol_coverage"), errors="coerce")
        _polbad = ~master_df.get("pol", "").astype(str).isin(["Intact", "Intact_FS_End"])
        _bodyless = (master_df["Structure"] == "Provirus_from_Multi") & _polbad & (_polcov < 0.30)
        if _bodyless.any():
            master_df.loc[_bodyless, "pol_qc_flag"] = "BODYLESS_UNIT_POL_ABSENT"
            logging.info(f"Bodyless tandem-unit QC flag: flagged {int(_bodyless.sum())} "
                         f"Provirus_from_Multi row(s) whose pol body is absent "
                         f"(loci: {sorted(master_df.loc[_bodyless, 'Locus'].unique())}); "
                         f"Structure and physical copy identity preserved.")

        # Sort master dataframe nicely
        master_df.sort_values(by=["Locus", "ID", "Haplotype"], inplace=True)

        # One state-weighted/marginalized consequence row per biological
        # sample-haplotype-locus. Raw observation rows remain in the master as
        # evidence but are never duplicate-counted and no observed copy is chosen
        # as founder truth.
        modeled = master_df[master_df.get('presence_call', '').eq('present')].copy()
        summary_rows = []
        for (locus, sample, hap), group in modeled.groupby(['Locus', 'ID', 'Haplotype'], dropna=False):
            weights = pd.to_numeric(group['observation_weight'], errors='raise')
            row = {'Locus': locus, 'ID': sample, 'Haplotype': hap,
                   'observation_weight_sum': weights.sum(),
                   'copy_state_single_weight': (weights * pd.to_numeric(group['copy_state_single_weight'])).sum(),
                   'copy_state_artifact_weight': (weights * pd.to_numeric(group['copy_state_artifact_weight'])).sum(),
                   'copy_state_true_dup_weight': (weights * pd.to_numeric(group['copy_state_true_dup_weight'])).sum(),
                   'expected_biological_copy_count': (weights * pd.to_numeric(group['expected_biological_copy_count'])).sum(),
                   'copy_state_weight_provenance': ';'.join(sorted(set(group['copy_state_weight_provenance'].astype(str)))),
                   'consequence_marginalization_provenance': ';'.join(sorted(set(group['consequence_marginalization_provenance'].astype(str))))}
            for col in COPY_EXPECTATION_COLS:
                row[col] = (weights * pd.to_numeric(group[col])).sum()
            summary_rows.append(row)
        summary_cols = ['Locus', 'ID', 'Haplotype', 'observation_weight_sum', *COPY_STATE_WEIGHT_COLS,
                        'expected_biological_copy_count', *COPY_EXPECTATION_COLS,
                        'copy_state_weight_provenance', 'consequence_marginalization_provenance']
        summary_df = pd.DataFrame(summary_rows, columns=summary_cols)
        summary_path = os.path.join(args.results_dir, "combined_hml2_orf_state_weighted_summary.tsv")
        atomic_write_dataframe(summary_df, summary_path, sep='\t', index=False, na_rep='NA')
        logging.info("Saved state-weighted/marginalized ORF summary with %d sample-locus-haplotype rows to %s",
                     len(summary_df), summary_path)

        # Atomic write: the master combined TSV is the pipeline's primary product and is consumed by
        # downstream steps -- it must be produced whole, never half-written, and be restartable.
        atomic_write_dataframe(master_df, master_tsv_path, sep='\t', index=False, na_rep='NA')
        logging.info(f"Saved master combined TSV with {len(master_df)} rows to {master_tsv_path}")
    else:
        logging.warning("No data found to combine into a master TSV.")

    logging.info("--- Finished ---")
    logging.info(f"Total placeholders added across all loci: {total_placeholders_added}")

if __name__ == "__main__":
    main()
