#!/usr/bin/env python3

import os
import subprocess
import glob
import sys
import argparse
import logging
import re
import time
import tempfile
import hashlib
import json
from pathlib import Path


def atomic_write_dataframe(df, out_path, sep=",", **to_csv_kwargs):
    """Write a DataFrame to out_path ATOMICALLY (temp file in same dir -> fsync -> os.replace).

    The per-locus ORF results CSV is the input to the combine/validate step; a crash mid-write must
    not leave a truncated CSV that the combine step would silently ingest. Writing to a same-directory
    temp file and atomically renaming guarantees a reader sees the whole prior file or the whole new
    one, and makes re-running the locus safe (idempotent output)."""
    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    os.makedirs(out_dir, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(prefix="." + os.path.basename(out_path) + ".", suffix=".tmp", dir=out_dir)
    try:
        with os.fdopen(fd, "w", newline="") as fh:
            df.to_csv(fh, sep=sep, **to_csv_kwargs)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp_path, out_path)
    except BaseException:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise

try:
    from Bio import AlignIO, SeqIO, Align
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import pandas as pd
    from Bio.Data import CodonTable
except ImportError as e:
    print(f"Error: Required Python module not found. ({e})", file=sys.stderr)
    print("Please install required modules: pip install biopython pandas", file=sys.stderr)
    sys.exit(1)


# Below this fraction of a feature's KCON reference length covered by the sample,
# dna/protein identity is computed over too little sequence to be meaningful (a
# few matching bases can read as identity ~1.0 for an otherwise-absent feature),
# so we report identity as "NA" rather than a misleadingly high number. Coverage
# itself is still reported. Adjust here if a different floor is preferred.
COVERAGE_FLOOR_FOR_IDENTITY = 0.3

# --- Type-I (Delta292) Pol geometry ----------------------------------------------------------
# Type I is the exact 292-nt deletion at KCON [6501, 6793) (zero-based, half-open) that removes
# the C-terminus of Pol and shifts the reading frame. Consequences the ORF table MUST encode:
#   * The Type-I KCON reference's `pol` feature is deliberately truncated at 6501 (see
#     KCON_FEATURES_TYPE1 below), so a Type-I allele that matches that truncated reference scores
#     pol "Intact" -- but that is only the PHYSICAL presence of the (truncated) Pol region and its
#     RT/RNase-H/integrase sequence potential. It is NOT complete native Pol and NEVER a native
#     translated-Pol burden.
#   * A Type-I allele's only POSITIVE translation state for Pol is an ALTERED read-through in which
#     translation continues past the Delta292 junction into Env-derived sequence. That is reported
#     as a separate, explicitly type-specific field -- never merged into the physical `pol` verdict.
# These constants make the boundary explicit and single-sourced; they are used only to label the
# translation-compatibility fields, not to alter any existing physical/coverage call.
TYPE1_DELETION_START_0BASED = 6501   # first deleted KCON base (== Type-I pol reference end)
TYPE1_DELETION_END_0BASED = 6793     # one-past-last deleted KCON base; 6793 - 6501 == 292 nt
TYPE1_DELETION_LEN = TYPE1_DELETION_END_0BASED - TYPE1_DELETION_START_0BASED  # 292

# Feature/ORF verdicts that indicate a translationally-competent reading frame for that feature
# (an intact ORF, an in-frame end-frameshift that does not truncate the product, or a run-through
# with no premature stop). Shared by the aggregate call and the Pol translation-compatibility call
# so the two never drift apart.
TRANSLATION_COMPETENT_STATUSES = {"Intact", "Intact_FS_End", "Fragment_Intact", "no_stop"}


# --- Helper Functions ---

def extract_source_id(record_id, desc=""):
    if "assembly_coords:" in desc:
        return desc.split("assembly_coords:")[1].split()[0]
    elif '|' in record_id:
        return record_id.split('|', 1)[1].lstrip('>')
    else:
        return record_id.lstrip('>')


def extract_genomic_strand(desc):
    """Parse the 'genomic_strand:+|-' field the extraction scripts emit (orientation
    of the provirus relative to the reference genome, hg38/CHM13 — computed from the
    untangle `inv` / BAM sam_flag combined with the provirus's strand on the path).
    Returns '+', '-', or 'NA' when absent (e.g. older extractions)."""
    if "genomic_strand:" in desc:
        val = desc.split("genomic_strand:")[1].split()[0].strip()
        if val in ("+", "-"):
            return val
    return "NA"


def combine_genomic_strand(window_gs, vs_kcon):
    """Per-copy genomic strand. The extraction emits the WINDOW's genomic strand for
    its KCON-normalized primary provirus; a copy that is reverse-complemented vs KCON
    (an inverted duplication) is on the opposite genomic strand, so flip accordingly."""
    if window_gs not in ("+", "-"):
        return "NA"
    if vs_kcon == "-":
        return "-" if window_gs == "+" else "+"
    return window_gs


def shift_assembly_coords(source_id, q_start, q_end, strand="+"):
    """Narrow a window's source id to a sub-span's true genomic coordinates.

    source_id looks like 'sample#hap#contig:gstart-gend' (the extracted window).
    [q_start, q_end] is a sub-interval in the WRITTEN-sequence (provirus-forward)
    coordinates -- e.g. one provirus of a multi-insertion window. Returns the source id
    with the genomic span replaced by that copy's own coordinates, so each split row
    reports where it actually is.

    The written sequence is provirus-forward. A '+' window is [gstart..gend] read
    forward, so written position q maps to gstart+q. A '-' window's written sequence is
    the REVERSE COMPLEMENT of [gstart..gend], so written position q maps from the other
    end: gend-q. Using the '+' formula on a '-' window mislocates every copy by the
    window length minus its offset (the 7p22.1 / reverse-strand _MULTI coordinate bug:
    a 2-copy array at 4704800-4722775 was reported at 4729180-4747141). `strand` is the
    window's genomic strand (the extraction emits it as genomic_strand: in the header).
    Returns source_id unchanged if it can't be parsed."""
    try:
        prefix, span = source_id.rsplit(':', 1)
        gs_str, ge_str = span.split('-')
        gs, ge = int(gs_str), int(ge_str)
    except (ValueError, IndexError):
        return source_id
    if strand == "-":
        # written-q is measured from the 3' (gend) end of the genomic window
        new_s = max(gs, ge - q_end)
        new_e = min(ge, ge - q_start)
    else:
        new_s = gs + q_start
        new_e = min(gs + q_end, ge)
    return f"{prefix}:{new_s}-{new_e}"


def safe_read_fasta(filepath, retries=3, delay=2):
    for _ in range(retries):
        try:
            if os.path.getsize(filepath) == 0:
                return []
            return list(SeqIO.parse(filepath, "fasta"))
        except OSError as e:
            if getattr(e, 'errno', None) == 116:
                logging.warning(f"Stale file handle for {filepath}. Retrying in {delay}s...")
                time.sleep(delay)
            else:
                logging.error(f"OSError parsing {filepath}: {e}")
                break
        except Exception as e:
            logging.error(f"Exception parsing {filepath}: {e}")
            break
    try:
        return list(SeqIO.parse(filepath, "fasta"))
    except Exception as e:
        logging.error(f"Final fallback parse failed for {filepath}: {e}")
        return []


def get_aligned_strings(alignment):
    """Reconstruct gapped target/query strings from alignment.aligned."""
    target = str(alignment.target).upper()
    query = str(alignment.query).upper()
    target_blocks, query_blocks = alignment.aligned

    aligned_target_parts = []
    aligned_query_parts = []
    prev_t_end, prev_q_end = 0, 0

    for (t_start, t_end), (q_start, q_end) in zip(target_blocks, query_blocks):
        t_start, t_end = int(t_start), int(t_end)
        q_start, q_end = int(q_start), int(q_end)
        t_gap = t_start - prev_t_end
        q_gap = q_start - prev_q_end
        if t_gap > 0:
            aligned_target_parts.append(target[prev_t_end:t_start])
            aligned_query_parts.append('-' * t_gap)
        if q_gap > 0:
            aligned_target_parts.append('-' * q_gap)
            aligned_query_parts.append(query[prev_q_end:q_start])
        aligned_target_parts.append(target[t_start:t_end])
        aligned_query_parts.append(query[q_start:q_end])
        prev_t_end, prev_q_end = t_end, q_end

    if prev_t_end < len(target):
        aligned_target_parts.append(target[prev_t_end:])
        aligned_query_parts.append('-' * (len(target) - prev_t_end))
    if prev_q_end < len(query):
        aligned_target_parts.append('-' * (len(query) - prev_q_end))
        aligned_query_parts.append(query[prev_q_end:])

    return "".join(aligned_target_parts), "".join(aligned_query_parts)


def construct_frame_corrected_protein(gapped_sample_dna, gapped_ref_dna, feature_name="Unknown", seq_id="Unknown"):
    """Block processing with delayed-debt frame correction.
    Returns (protein, frameshift_codons, large_insertions) where large_insertions
    is a list of (codon_loc, inserted_bp) for insertion blocks longer than one codon."""
    temp_dna = []
    fs_locations = []
    large_insertions = []
    ref_base_idx = 0
    net_shift = 0
    n_debt = 0
    fs_active = False

    aln_r = str(gapped_ref_dna).upper()
    aln_s = str(gapped_sample_dna).upper()

    i = 0
    while i < len(aln_r):
        r_char = aln_r[i]
        s_char = aln_s[i]

        # 0. Shared-gap column: belongs to an insertion in some OTHER MSA member.
        #    Not informative for this ref/sample pair -> skip entirely.
        if r_char == '-' and s_char == '-':
            i += 1
            continue

        # 1. Handle Contiguous Insertion Blocks (ref is gap)
        if r_char == '-':
            block_ref_start = ref_base_idx
            ins_bases = []
            while i < len(aln_r) and aln_r[i] == '-':
                if aln_s[i] != '-':            # count only REAL inserted bases
                    ins_bases.append(aln_s[i])
                i += 1

            ins_len = len(ins_bases)
            if ins_len == 0:
                continue                       # block was entirely shared gaps

            codon_loc = (block_ref_start // 3) + 1

            # Structural insertion (> 1 codon): record as a single event, excise it,
            # and do NOT perturb the reading frame of the homologous flanks.
            if ins_len > 3:
                large_insertions.append((codon_loc, ins_len))
                logging.debug(f"[{seq_id} | {feature_name}] --> LARGE INSERTION {ins_len}bp at codon {codon_loc} (excised)")
                continue

            # Small insertion (<= 1 codon): keep per-base behavior
            for b in ins_bases:
                temp_dna.append(b)
                net_shift += 1
                n_debt = max(0, n_debt - 1)

            if net_shift % 3 != 0:
                if not fs_active:
                    fs_locations.append(codon_loc)
                    fs_active = True
                    logging.debug(f"[{seq_id} | {feature_name}] --> FRAMESHIFT at codon {codon_loc} (Net shift: {net_shift})")
            else:
                fs_active = False

            rem = len(temp_dna) % 3
            if rem != 0:
                temp_dna.extend(['N'] * (3 - rem))
            continue

        # 2. Handle Contiguous Deletion Blocks (sample is gap)
        elif s_char == '-':
            del_count = 0
            block_ref_start = ref_base_idx
            while i < len(aln_r) and aln_s[i] == '-':
                if aln_r[i] != '-':            # count only REAL deleted ref bases
                    del_count += 1
                    ref_base_idx += 1
                i += 1

            if del_count == 0:
                continue                       # block was entirely shared gaps

            n_debt += del_count
            net_shift -= del_count

            if net_shift % 3 != 0:
                if not fs_active:
                    fs_locations.append((block_ref_start // 3) + 1)
                    fs_active = True
                    logging.debug(f"[{seq_id} | {feature_name}] --> FRAMESHIFT at codon {(block_ref_start // 3) + 1} (Net shift: {net_shift})")
            else:
                fs_active = False

            if n_debt > 0 and len(temp_dna) > 0 and len(temp_dna) % 3 == 0:
                temp_dna.extend(['N'] * n_debt)
                n_debt = 0
            continue

        # 3. Handle Matches / Mismatches
        else:
            temp_dna.append(s_char)
            ref_base_idx += 1
            if n_debt > 0 and len(temp_dna) % 3 == 0:
                temp_dna.extend(['N'] * n_debt)
                n_debt = 0
            i += 1

    if n_debt > 0:
        temp_dna.extend(['N'] * n_debt)

    corrected_dna_str = "".join(temp_dna)
    rem = len(corrected_dna_str) % 3
    if rem > 0:
        corrected_dna_str = corrected_dna_str[:-rem]

    try:
        prot = str(Seq(corrected_dna_str).translate(table=1, cds=False))
        return prot, sorted(list(set(fs_locations))), large_insertions
    except Exception as e:
        logging.error(f"Translation failed for {seq_id}: {e}")
        return None, [], large_insertions


def analyze_orf_structural_integrity(aligned_dna_seq, ref_protein_seq, ref_dna_seq,
                                     feature_name="", locus_type="TypeII",
                                     upstream_orf_functional=False, gapped_ref_dna=None, seq_id="Unknown"):
    ungapped_dna = aligned_dna_seq.replace("-", "").upper()
    if set(ungapped_dna.replace('N', '')) - set("ATGC"):
        return "invalid_chars", None, None, [], None, []

    check_start = any(x in feature_name.lower() for x in ['gag', 'env', 'np9', 'rec'])
    if feature_name.lower() == 'env':
        if not upstream_orf_functional:
            check_start = True
        else:
            if locus_type == 'TypeI':
                check_start = False
            elif len(ref_dna_seq) - len(ungapped_dna) > 80:
                check_start = False

    start_lost = False
    if check_start:
        if len(ungapped_dna) < 3 or ungapped_dna[:3] != "ATG":
            start_lost = True

    if len(ref_dna_seq) > 0:
        coverage = len(ungapped_dna) / len(ref_dna_seq)
        if coverage < 0.30:
            return "Deletion", None, None, [], None, []
    if len(ungapped_dna) < 3:
        return "too_short", None, None, [], None, []
        
    global_remainder = len(ungapped_dna) % 3
    is_global_frameshift = (global_remainder != 0)
    dna_to_translate = ungapped_dna[:-global_remainder] if global_remainder > 0 else ungapped_dna

    try:
        original_protein = str(Seq(dna_to_translate).translate(table=1, cds=False))
    except Exception:
        return "translation_error", None, None, [], None, []
        
    first_stop_aa_idx = original_protein.find('*')
    is_natural_stop = (first_stop_aa_idx == len(original_protein) - 1)

    if gapped_ref_dna:
        corrected_protein, all_fs_locations, large_insertions = construct_frame_corrected_protein(
            aligned_dna_seq, gapped_ref_dna, feature_name=feature_name, seq_id=seq_id)
    else:
        corrected_protein, all_fs_locations, large_insertions = None, [], []

    if not all_fs_locations and is_global_frameshift:
        all_fs_locations = [len(ref_protein_seq)]
        
    verdict = "Intact"
    stop_return = None
    has_frameshift = len(all_fs_locations) > 0

    if first_stop_aa_idx != -1 and not is_natural_stop:
        stop_return = first_stop_aa_idx
        if has_frameshift and all_fs_locations[0] <= first_stop_aa_idx + 2:
            verdict = "Frameshift"
        else:
            verdict = "Nonsense"
    elif has_frameshift:
        if all_fs_locations[0] < len(ref_protein_seq) * 0.4:
            verdict = "Frameshift"
        else:
            verdict = "Intact_FS_End"
    elif is_global_frameshift:
        verdict = "Intact_FS_End"

    if start_lost:
        if "Intact" in verdict or verdict == "Fragment_Intact":
            verdict = "Start_Lost"
    if "Intact" in verdict and verdict != "Intact_FS_End":
        eff_len = len(original_protein) - 1 if original_protein.endswith('*') else len(original_protein)
        if eff_len < len(ref_protein_seq) * 0.6:
            verdict = "Fragment_Intact"

    if large_insertions and "Intact" in verdict:
        verdict = "Insertion"
    return verdict, original_protein, corrected_protein, all_fs_locations, stop_return, large_insertions


_DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _reverse_complement_dna(sequence_str):
    return str(sequence_str).translate(_DNA_COMPLEMENT)[::-1].upper()


def _append_tsd_junction_evidence(payload):
    """Append one direct junction observation when a review path is requested."""
    output_path = os.environ.get("HML2_TSD_JUNCTION_EVIDENCE_JSONL", "").strip()
    if not output_path:
        return
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "a") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def _refine_projected_ltr_junction(sequence_str, ltr_reference, projected,
                                   boundary, strand, max_shift=2,
                                   terminal_span=32):
    """Return a terminal-motif-supported LTR-host junction, else ``None``.

    The projected coordinate remains the anchor.  Candidate shifts are limited
    to two bases and are scored against the appropriate terminal segment of the
    LTR reference in query orientation. Host-flank similarity and the putative
    TSD pair never participate in boundary selection. A least-bad projected
    coordinate is not evidence: the terminal 32-mer itself must match within two
    substitutions at one unambiguous offset.
    """
    seq = str(sequence_str).upper()
    ref = str(ltr_reference).upper()
    if boundary not in {"left", "right"} or strand not in {"+", "-"}:
        return None
    span = min(terminal_span, len(ref))
    if span < 12:
        return None
    if boundary == "left":
        motif = ref[:span] if strand == "+" else _reverse_complement_dna(ref[-span:])
    else:
        motif = ref[-span:] if strand == "+" else _reverse_complement_dna(ref[:span])

    candidates = []
    for shift in range(-max_shift, max_shift + 1):
        coordinate = projected + shift
        if coordinate < 0 or coordinate > len(seq):
            continue
        if boundary == "left":
            segment = seq[coordinate:coordinate + span]
        else:
            segment = seq[coordinate - span:coordinate]
        if len(segment) != span or set(segment) - set("ACGT"):
            continue
        substitutions = sum(left != right for left, right in zip(segment, motif))
        candidates.append((substitutions, abs(shift), shift, coordinate))
    if not candidates:
        return None
    best = min(candidates)
    if best[0] > 2:
        return None
    # Do not choose a biological junction when two offsets explain the terminal
    # LTR equally well.
    if any(item[0] == best[0] and item[3] != best[3] for item in candidates):
        return None
    return best[3]


def refine_tsd_search(sequence_str, projected_start, projected_end, search_window=20):
    """Read the observed host flanks at independently resolved LTR junctions.

    Coordinates are fixed before this function is called.  For an observed pair,
    retain the longest 4--6 bp suffix/prefix with at most two present-day
    differences.  If both boundaries are observed but no length is admitted,
    return an explicit paired no-call rather than serializing discordant host
    flanks as TSDs.  A genuinely one-sided boundary remains explicit one-sided
    evidence and contributes no paired likelihood downstream.
    ``search_window`` is retained only as the existing call signature.
    """
    del search_window
    seq = str(sequence_str).upper()
    n = len(seq)
    if (
        projected_start is not None
        and (projected_start < 0 or projected_start > n)
    ) or (
        projected_end is not None
        and (projected_end < 0 or projected_end > n)
    ) or (
        projected_start is not None
        and projected_end is not None
        and projected_end < projected_start
    ):
        return "NONE", "NONE", projected_start, projected_end
    left = (
        seq[projected_start - 6:projected_start]
        if projected_start is not None and projected_start >= 6 else ""
    )
    right = (
        seq[projected_end:projected_end + 6]
        if projected_end is not None and projected_end + 6 <= n else ""
    )
    if len(left) != 6 or set(left) - set("ACGT"):
        left = "NONE"
    if len(right) != 6 or set(right) - set("ACGT"):
        right = "NONE"
    if left != "NONE" and right != "NONE":
        for size in (6, 5, 4):
            candidate_left = left[-size:]
            candidate_right = right[:size]
            if sum(a != b for a, b in zip(candidate_left, candidate_right)) <= 2:
                return candidate_left, candidate_right, projected_start, projected_end
    if (left == "NONE") != (right == "NONE"):
        return left, right, projected_start, projected_end
    return "NONE", "NONE", projected_start, projected_end


def setup_logging(log_file_path, is_debug_mode=False):
    log_dir = os.path.dirname(log_file_path)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(funcName)s - L%(lineno)d - %(message)s')
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()
    fh = logging.FileHandler(log_file_path, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(file_formatter)
    root_logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    
    if is_debug_mode:
        sh.setLevel(logging.DEBUG)
    else:
        sh.setLevel(logging.INFO)
        
    sh.setFormatter(console_formatter)
    root_logger.addHandler(sh)
    logging.info(f"Main logging initialized. Log file: {log_file_path}")


def setup_detailed_logging(log_file_path):
    detail_logger = logging.getLogger('LTR_Details')
    detail_logger.propagate = False
    detail_logger.setLevel(logging.DEBUG)
    for handler in detail_logger.handlers[:]:
        detail_logger.removeHandler(handler)
        handler.close()
    fh = logging.FileHandler(log_file_path, mode='w')
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    fh.setFormatter(formatter)
    detail_logger.addHandler(fh)
    return detail_logger


def calculate_levenshtein(s1, s2):
    if len(s1) < len(s2):
        return calculate_levenshtein(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def _get_all_motif_matches(sequence_str, motif, search_window_start, search_window_end, score_threshold):
    search_sequence = sequence_str[search_window_start:search_window_end]
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2

    alignments = aligner.align(search_sequence, motif)
    good_matches = []
    seen_bounds = set()

    for alignment in alignments:
        if alignment.score < score_threshold:
            break
        target_bounds = alignment.aligned[0]
        if len(target_bounds) > 0:
            begin = int(target_bounds[0][0])
            end = int(target_bounds[-1][1])
            bound_tuple = (begin, end)
            if bound_tuple not in seen_bounds:
                seen_bounds.add(bound_tuple)
                good_matches.append({
                    'start': search_window_start + begin,
                    'end': search_window_start + end,
                    'score': alignment.score,
                    'found_motif': search_sequence[begin:end]
                })
    return good_matches


def refine_boundaries_with_motif_search(sequence, initial_start, initial_end, search_radius=50):
    MOTIF_5PRIME = "TGTGGGGAAAAGCAAGAG"
    MOTIF_3PRIME = "CAACCCACCCCTACA"
    seq_str = str(sequence).upper()

    search_5_start = max(0, initial_start - search_radius)
    search_5_end = min(len(seq_str), initial_start + search_radius)
    score_threshold_5 = len(MOTIF_5PRIME) * 2 * 0.7
    all_5_prime_matches = _get_all_motif_matches(seq_str, MOTIF_5PRIME, search_5_start, search_5_end, score_threshold_5)

    if not all_5_prime_matches:
        return initial_start, initial_end, None

    best_5_prime_match = min(all_5_prime_matches, key=lambda x: x['start'])
    refined_ltr_start = best_5_prime_match['start']
    search_3_start = max(0, initial_end - search_radius)
    search_3_end = min(len(seq_str), initial_end + search_radius)
    score_threshold_3 = len(MOTIF_3PRIME) * 2 * 0.7
    all_3_prime_matches = _get_all_motif_matches(seq_str, MOTIF_3PRIME, search_3_start, search_3_end, score_threshold_3)
    if not all_3_prime_matches:
        return initial_start, initial_end, None
    best_3_prime_match = max(all_3_prime_matches, key=lambda x: x['end'])
    refined_ltr_end = best_3_prime_match['end']
    # Both boundaries are supported by independent terminal-LTR motifs.  The
    # flanking TSD strings must not select, reject, or move these coordinates.
    return refined_ltr_start, refined_ltr_end, None


def find_ltr_chain_boundaries(fasta_path, ltr_ref_path, paf_output_path, detail_logger, threads=1):
    MIN_LTR_ALIGN_LEN = 750
    MIN_LTR_MAPQ = 40
    minimap2_cmd = ["minimap2", "-x", "asm20", "-A", "2", "-B", "2", "-O", "2,12", "-E", "1,0",
                    "-t", str(threads), ltr_ref_path, fasta_path]
    try:
        with open(paf_output_path, "w") as paf_file:
            subprocess.run(minimap2_cmd, stdout=paf_file, stderr=subprocess.PIPE, text=True, check=True)
        filtered_hits = []
        with open(paf_output_path, 'r') as paf:
            for line in paf:
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue
                try:
                    hit_info = {'start': int(fields[2]), 'end': int(fields[3]),
                                'align_len': int(fields[10]), 'mapq': int(fields[11])}
                    if hit_info['align_len'] > MIN_LTR_ALIGN_LEN and hit_info['mapq'] > MIN_LTR_MAPQ:
                        filtered_hits.append(hit_info)
                except (ValueError, IndexError):
                    continue
        if detail_logger is not None:
            qname = Path(fasta_path).name
            detail_logger.info(f"[{qname}] LTR chain scan: {len(filtered_hits)} hit(s) passed "
                               f"(align_len>{MIN_LTR_ALIGN_LEN}, mapq>{MIN_LTR_MAPQ}).")
            for h in sorted(filtered_hits, key=lambda x: x['start']):
                detail_logger.info(f"[{qname}]   LTR hit: {h['start']}-{h['end']} "
                                   f"(align_len={h['align_len']}, mapq={h['mapq']}).")
        if not filtered_hits:
            return None, None
        filtered_hits.sort(key=lambda x: x['start'])
        if detail_logger is not None:
            detail_logger.info(f"[{Path(fasta_path).name}] LTR chain boundaries: "
                               f"{filtered_hits[0]['start']}-{filtered_hits[-1]['end']}.")
        return filtered_hits[0]['start'], filtered_hits[-1]['end']
    except Exception:
        return None, None


def _solo_ltr_bounds_with_ref(solo_fasta_path, ltr_ref_path, paf_output_path, threads=1):
    minimap2_cmd = ["minimap2", "-c", "-x", "asm20", "-w", "5", "-A", "1", "-B", "2",
                    "-t", str(threads), "--secondary=no", ltr_ref_path, solo_fasta_path]
    try:
        with open(paf_output_path, "w") as paf_file:
            subprocess.run(minimap2_cmd, stdout=paf_file, stderr=subprocess.PIPE, text=True, check=True)
        best_hit = None
        max_score = -1
        with open(paf_output_path, 'r') as paf:
            for line in paf:
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue
                try:
                    mapq = int(fields[11])
                    align_len = int(fields[10])
                    q_start, q_end = int(fields[2]), int(fields[3])
                    strand = fields[4]
                    t_len, t_start, t_end = (
                        int(fields[6]), int(fields[7]), int(fields[8]))
                    score = mapq + align_len
                    if score > max_score:
                        max_score = score
                        best_hit = (
                            q_start, q_end, strand, t_len, t_start, t_end)
                except (ValueError, IndexError):
                    continue
        if best_hit:
            records = safe_read_fasta(solo_fasta_path)
            if records:
                seq_str = str(records[0].seq).upper()
                ltr_records = safe_read_fasta(ltr_ref_path)
                if not ltr_records:
                    return None, None, None, None
                ltr_reference = str(ltr_records[0].seq).upper()
                q_start, q_end, strand, t_len, t_start, t_end = best_hit
                if strand == "+":
                    projected_start = q_start - t_start
                    projected_end = q_end + (t_len - t_end)
                    direct_start = q_start if t_start == 0 else None
                    direct_end = q_end if t_end == t_len else None
                else:
                    projected_start = q_start - (t_len - t_end)
                    projected_end = q_end + t_start
                    direct_start = q_start if t_end == t_len else None
                    direct_end = q_end if t_start == 0 else None
                resolved_start = direct_start
                left_terminal_support = "DIRECT_PAF_ENDPOINT"
                if resolved_start is None:
                    resolved_start = _refine_projected_ltr_junction(
                        seq_str, ltr_reference, projected_start, "left", strand)
                    left_terminal_support = (
                        "TERMINAL_MOTIF" if resolved_start is not None
                        else "UNSUPPORTED"
                    )
                resolved_end = direct_end
                right_terminal_support = "DIRECT_PAF_ENDPOINT"
                if resolved_end is None:
                    resolved_end = _refine_projected_ltr_junction(
                        seq_str, ltr_reference, projected_end, "right", strand)
                    right_terminal_support = (
                        "TERMINAL_MOTIF" if resolved_end is not None
                        else "UNSUPPORTED"
                    )
                tsd_5, tsd_3, ref_start, ref_end = refine_tsd_search(
                    seq_str, resolved_start, resolved_end)
                if ref_start is None and ref_end is None:
                    _append_tsd_junction_evidence({
                        "source_fasta": os.path.realpath(solo_fasta_path),
                        "source_sequence_id": str(records[0].id),
                        "structure": "SOLO_LTR",
                        "source_length": len(seq_str),
                        "source_strand_to_insertion": strand,
                        "left_terminal_support": left_terminal_support,
                        "right_terminal_support": right_terminal_support,
                        "resolution_status": "OUTWARD_TERMINAL_NOT_SUPPORTED",
                        "tsd_5prime_in_insertion_orientation": "NONE",
                        "tsd_3prime_in_insertion_orientation": "NONE",
                        "boundary_owner": "TERMINAL_LTR_ALIGNMENT_ONLY",
                    })
                    return None, None, "NONE", "NONE"
                if strand == "-":
                    tsd_5, tsd_3 = (
                        _reverse_complement_dna(tsd_3)
                        if tsd_3 != "NONE" else "NONE",
                        _reverse_complement_dna(tsd_5)
                        if tsd_5 != "NONE" else "NONE",
                    )
                    evidence_seq = _reverse_complement_dna(seq_str)
                    evidence_left = (
                        len(seq_str) - ref_end if ref_end is not None else None
                    )
                    evidence_right = (
                        len(seq_str) - ref_start if ref_start is not None else None
                    )
                    evidence_projected_left = len(seq_str) - projected_end
                    evidence_projected_right = len(seq_str) - projected_start
                    left_terminal_support, right_terminal_support = (
                        right_terminal_support, left_terminal_support
                    )
                else:
                    evidence_seq = seq_str
                    evidence_left = ref_start
                    evidence_right = ref_end
                    evidence_projected_left = projected_start
                    evidence_projected_right = projected_end
                _append_tsd_junction_evidence({
                    "source_fasta": os.path.realpath(solo_fasta_path),
                    "source_sequence_id": str(records[0].id),
                    "structure": "SOLO_LTR",
                    "source_length": len(seq_str),
                    "source_strand_to_insertion": strand,
                    "projected_left_boundary": evidence_projected_left,
                    "projected_right_boundary": evidence_projected_right,
                    "resolved_left_boundary": evidence_left,
                    "resolved_right_boundary": evidence_right,
                    "left_boundary_shift": (
                        evidence_left - evidence_projected_left
                        if evidence_left is not None else ""
                    ),
                    "right_boundary_shift": (
                        evidence_right - evidence_projected_right
                        if evidence_right is not None else ""
                    ),
                    "left_junction_context_in_insertion_orientation": (
                        evidence_seq[max(0, evidence_left - 18):
                                     min(len(evidence_seq), evidence_left + 18)]
                        if evidence_left is not None else ""
                    ),
                    "left_context_boundary_offset": (
                        min(18, evidence_left) if evidence_left is not None else ""
                    ),
                    "right_junction_context_in_insertion_orientation": (
                        evidence_seq[max(0, evidence_right - 18):
                                     min(len(evidence_seq), evidence_right + 18)]
                        if evidence_right is not None else ""
                    ),
                    "right_context_boundary_offset": (
                        min(18, evidence_right) if evidence_right is not None else ""
                    ),
                    "tsd_5prime_in_insertion_orientation": tsd_5,
                    "tsd_3prime_in_insertion_orientation": tsd_3,
                    "left_terminal_support": left_terminal_support,
                    "right_terminal_support": right_terminal_support,
                    "boundary_owner": "TERMINAL_LTR_ALIGNMENT_ONLY",
                    "resolution_status": "RESOLVED",
                })
                return ref_start, ref_end, tsd_5, tsd_3
        return None, None, None, None
    except Exception:
        return None, None, None, None


def get_solo_ltr_bounds(solo_fasta_path, ltr_ref_path, paf_output_path, threads=1):
    """Locate the solo-LTR boundary (and read its TSD) by mapping a solo-LTR reference.
    Tries the primary (KCON / LTR5_Hs) reference first, then falls back to the older
    LTR5A / LTR5B subfamily references for DIVERGENT solo-LTRs the LTR5_Hs reference cannot
    anchor -- e.g. 11p15.4b is an LTR5B solo-LTR that gives no LTR5_Hs hit, so its TSD was
    never recovered and the whole locus dropped out of the TSD analysis. The fallback fires
    ONLY when the primary returns no boundary, so existing (mappable) solo-LTR calls are
    byte-for-byte unchanged. Fallback references are looked up alongside the primary ref;
    if absent, behaviour is exactly the old single-reference path."""
    result = _solo_ltr_bounds_with_ref(solo_fasta_path, ltr_ref_path, paf_output_path, threads)
    if result[0] is not None or result[1] is not None:
        return result
    ref_dir = os.path.dirname(ltr_ref_path)
    for fb in ("LTR5A_ref.fa", "LTR5B_ref.fa"):
        fb_path = os.path.join(ref_dir, fb)
        if not os.path.exists(fb_path):
            continue
        result = _solo_ltr_bounds_with_ref(solo_fasta_path, fb_path,
                                           paf_output_path + "." + fb, threads)
        if result[0] is not None or result[1] is not None:
            logging.info(f"solo-LTR boundary recovered via fallback {fb} for {Path(solo_fasta_path).name}")
            return result
    return None, None, None, None


def refine_tsd_host_flank(seq_str, ltr_start, ltr_end, jitter=4):
    """Read exact 6-bp host flanks at alignment-projected LTR termini.

    Terminal clipping may move ``ltr_start``/``ltr_end`` through the independent
    LTR projection in :func:`get_provirus_tsds`.  The observed TSD pair itself
    never selects an offset.  An unavailable side is retained as missing and is
    neutral downstream.
    """
    del jitter
    left, right, _, _ = refine_tsd_search(seq_str, ltr_start, ltr_end, 0)
    return left, right


def get_provirus_tsds(
    provirus_fasta_path, ltr_ref_path, paf_output_path, threads=1,
    kcon_ref_full=None, evidence_source_fasta=None,
):
    """TSD detection for ordinary (single) proviruses.

    Computed from the ORIGINAL untrimmed input sequence so the flanking genomic context (where the
    TSD lives) is still present. Maps the solo-LTR reference onto the provirus and PROJECTS each
    terminal LTR hit out to the LTR-consensus terminus (cons 0 / cons len). The projection is
    essential: minimap2 clips the divergent U3/U5 LTR termini by ~200-450bp, so the raw alignment
    extremes sit well INSIDE the LTR. Because the 5' and 3' LTRs are near-identical, reading a TSD
    around those clipped extremes matched an LTR-internal k-mer (LTR<->LTR homology) instead of the
    host TSD. Anchoring at the projected terminus places the TSD search in true host flank."""
    def record_unresolved(reason, **fields):
        payload = {
            "source_fasta": os.path.realpath(
                evidence_source_fasta or provirus_fasta_path
            ),
            "structure": "PROVIRUS_OR_ARRAY",
            "resolution_status": reason,
            "tsd_5prime_in_insertion_orientation": "NONE",
            "tsd_3prime_in_insertion_orientation": "NONE",
            "boundary_owner": "TERMINAL_LTR_ALIGNMENT_ONLY",
        }
        payload.update(fields)
        _append_tsd_junction_evidence(payload)

    minimap2_cmd = ["minimap2", "-x", "asm20", "-A", "2", "-B", "2", "-O", "2,12", "-E", "1,0",
                    "-t", str(threads), ltr_ref_path, provirus_fasta_path]
    try:
        with open(paf_output_path, "w") as paf_file:
            subprocess.run(minimap2_cmd, stdout=paf_file, stderr=subprocess.PIPE, text=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logging.warning(f"minimap2 unavailable/failed for provirus TSD search ({Path(provirus_fasta_path).name}): {e}")
        record_unresolved("LTR_MAPPING_FAILED")
        return "NA", "NA"

    # Validate LTR hits against a STANDARD-scoring asm20 alignment. The lenient -A2 -B2 -E1,0
    # scoring above is essential to project the soft-clipped U3/U5 terminus, but it forces long
    # low-quality extensions, so n_match/align_len is deflated to ~0.17-0.23 for REAL divergent
    # LTRs (5p12 798bp@87% -> 0.17; 6p22.1 926bp@80% -> 0.19) -- indistinguishable from LTR-less
    # junk (12q24.33 -> 0.15-0.18). A ratio floor therefore NULLED real divergent LTRs (regressed
    # 5p12, 6p22.1, 19p12a/c, 19q12, 19q13.12b, 22q11.23 to NA). Instead: a hit is a real LTR iff
    # the SAME query region ALSO aligns under standard asm20 scoring -- real LTRs give a clean
    # 660-995bp block (80-96% id), LTR-less loci (12q24.33, 8p23.1b/c/d) give nothing.
    std_hits = None
    try:
        std_paf = paf_output_path + ".std"
        with open(std_paf, "w") as sf:
            # -c = base-level alignment: REQUIRED here. Without it minimap2 emits approximate
            # blocks even for LTR-less junk (12q24.33 -> 593bp "hits"), defeating the check; with
            # -c, real LTRs give 683-995bp blocks and LTR-less loci give nothing.
            subprocess.run(["minimap2", "-c", "-x", "asm20", "-t", str(threads), ltr_ref_path, provirus_fasta_path],
                           stdout=sf, stderr=subprocess.PIPE, text=True, check=True)
        std_hits = []
        with open(std_paf) as sf:
            for line in sf:
                f = line.strip().split('\t')
                if len(f) < 12:
                    continue
                try:
                    if int(f[10]) >= 400:
                        std_hits.append((
                            int(f[2]), int(f[3]), int(f[7]), int(f[8]),
                            int(f[6]), f[4],
                        ))
                except (ValueError, IndexError):
                    continue
    except (subprocess.CalledProcessError, FileNotFoundError):
        std_hits = None

    def _is_real_ltr(qs, qe):
        if std_hits is None:
            return False
        for (a, b, _ts, _te, _tl, _strand) in std_hits:
            if qs < b and a < qe:
                return True
        return False

    # Body coverage (full provirus consensus) -> lets the TERMINUS selection below tell a
    # tandem-repeated HML-2 body filling an inter-LTR gap (KEEP -- part of the array) apart from
    # a non-HML-2 host gap before an aberrantly captured solo-LTR (SPLIT it off). This refines
    # only WHICH LTR is the 3' terminus for the TSD; the _MULTI tandem-array flag itself is set
    # upstream (extraction / split_multi_provirus_sequence) and is untouched -- a _MULTI window
    # with a stray trailing solo-LTR stays _MULTI, we just read its TSD from the array's real
    # outer LTRs, not the stray one.
    body_hits = []
    if kcon_ref_full:
        try:
            body_paf = paf_output_path + ".body"
            with open(body_paf, "w") as bf:
                subprocess.run(["minimap2", "-c", "-x", "asm20", "-t", str(threads),
                                kcon_ref_full, provirus_fasta_path],
                               stdout=bf, stderr=subprocess.PIPE, text=True, check=True)
            with open(body_paf) as bf:
                for line in bf:
                    f = line.strip().split('\t')
                    if len(f) >= 12:
                        try:
                            # Low bar: catch even small/degraded body fragments (gag/pro/pol/env).
                            # Host sequence does NOT align to the provirus consensus, so any block
                            # here is real HML-2 body; we only want to split when there is
                            # essentially NONE between the LTRs.
                            if int(f[10]) >= 200:
                                body_hits.append((
                                    int(f[2]), int(f[3]), int(f[7]), int(f[8]),
                                    int(f[6]),
                                ))
                        except (ValueError, IndexError):
                            continue
        except (subprocess.CalledProcessError, FileNotFoundError):
            body_hits = []

    ltr_hits = []
    with open(paf_output_path, 'r') as paf:
        for line in paf:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            try:
                q_start, q_end = int(fields[2]), int(fields[3])
                strand = fields[4]
                t_len, t_start, t_end = int(fields[6]), int(fields[7]), int(fields[8])
                n_match = int(fields[9])
                align_len = int(fields[10])
                # Substantial lenient block that is corroborated by a standard-asm20 LTR hit.
                if align_len > 500 and _is_real_ltr(q_start, q_end):
                    ltr_hits.append((q_start, q_end, t_start, t_end, t_len, strand))
            except (ValueError, IndexError):
                continue

    if not ltr_hits:
        record_unresolved("NO_CORROBORATED_TERMINAL_LTR")
        return "NA", "NA"

    records = safe_read_fasta(provirus_fasta_path)
    if not records:
        record_unresolved("NO_SOURCE_FASTA_RECORD")
        return "NA", "NA"
    seq_str = str(records[0].seq).upper()
    ltr_records = safe_read_fasta(ltr_ref_path)
    if not ltr_records:
        record_unresolved(
            "NO_LTR_REFERENCE_RECORD",
            source_sequence_id=str(records[0].id),
            source_length=len(seq_str),
        )
        return "NA", "NA"
    ltr_reference = str(ltr_records[0].seq).upper()
    ltr_length = len(ltr_reference)
    body_intervals = [
        (query_start, query_end)
        for query_start, query_end, target_start, target_end, target_length
        in body_hits
        if min(target_end, target_length - ltr_length)
        > max(target_start, ltr_length)
    ]

    # Orientation-normalize to provirus-forward. An INVERTED segdup copy sliced from a multi-copy
    # window has its LTRs on '-', so its TSD would be read reverse-complemented relative to the
    # other copies (e.g. 4p16.1a inverted dup: alt1 CAGAT/CTGAT vs alt2 ATCAG/ATCTG = the same TSD
    # RC'd). RC the sequence and flip each hit's query coords/strand so the TSD comes out in a
    # consistent orientation. Single proviruses are already provirus-forward ('+') -> no-op.
    source_strand_to_insertion = (
        "-" if sum(1 for h in ltr_hits if h[5] == "-") * 2 > len(ltr_hits) else "+"
    )
    if source_strand_to_insertion == "-":
        _L = len(seq_str)
        seq_str = seq_str.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
        ltr_hits = [(_L - qe, _L - qs, ts, te, tl, "+" if st == "-" else "-")
                    for (qs, qe, ts, te, tl, st) in ltr_hits]
        if std_hits is not None:
            std_hits = [
                (_L - qe, _L - qs, ts, te, tl, "+" if st == "-" else "-")
                for (qs, qe, ts, te, tl, st) in std_hits
            ]
        body_intervals = [(_L - be, _L - bs) for (bs, be) in body_intervals]

    # Project the outermost terminal LTR hits to the consensus termini (recovers the clipped
    # divergent ends). For a + hit, consensus-0 lies at q_start - t_start and consensus-len at
    # q_end + (t_len - t_end); the strand swaps these. Mirrors the extractor's LTR-edge projection.
    # Element-boundary rule: ONE provirus = 5' LTR ... 3' LTR. Once the element already holds its
    # closing (3') LTR, a further LTR separated by a real gap belongs to a DIFFERENT element and
    # must NOT define this provirus's 3' terminus. e.g. 3q12.3 = a clean 2-LTR provirus + a SEPARATE
    # solo-LTR 2.8kb downstream (- strand); the old h3=rightmost read the 3' TSD off that neighbour
    # -> CTGAG/CTGTC mismatch. Walk LTRs 5'->3': the FIRST inter-LTR gap is the provirus body (only
    # one LTR so far -> never splits), but a gap AFTER the element is already complete (>=2 LTRs)
    # splits off the neighbour. Gating the split on the element already being CLOSED is what keeps
    # 4q32.3's ~3kb gag-pol deletion (which PRECEDES the 3' LTR) as a single element.
    # TANDEM SAFETY: an inter-LTR gap that is FILLED BY A PROVIRUS BODY is a tandem-repeated HML-2
    # (part of an array) and is NEVER split -- so a tandem array of 3+ LTRs keeps ALL its bodies,
    # and the rule still drops a stray solo-LTR after the array (the final gap to it is host, not a
    # body). A gap NOT body-covered, after the element already holds its closing 3' LTR, splits off
    # a separate element (3q12.3's downstream solo-LTR). Body coverage (not the _MULTI tag) is the
    # discriminator. Missing body evidence fails closed for a two-LTR call.
    ELEMENT_SPLIT_GAP = 1500
    def _gap_body_covered(g0, g1):
        # "Has body" if >=10% of the gap is provirus body. A real tandem/element body fills most
        # of the gap (>>10%); a host gap before a captured solo-LTR has ~0% -- the 10% (vs a tiny
        # absolute bar) tolerates a little spurious junk mismapping to the consensus without
        # wrongly treating a host gap as body-filled. Only a near-bodyless gap splits.
        if g1 - g0 <= 0:
            return True
        covered = []
        for (b0, b1) in body_intervals:
            lo, hi = max(g0, b0), min(g1, b1)
            if hi > lo:
                covered.append((lo, hi))
        if not covered:
            return False
        current_lo, current_hi = sorted(covered)[0]
        cov = current_hi - current_lo
        for lo, hi in sorted(covered)[1:]:
            if lo <= current_hi:
                if hi > current_hi:
                    cov += hi - current_hi
                    current_hi = hi
            else:
                current_lo, current_hi = lo, hi
                cov += current_hi - current_lo
        return cov >= max(200, 0.10 * (g1 - g0))
    _by_pos = sorted(ltr_hits, key=lambda h: h[0])
    _elem = [_by_pos[0]]
    for _h in _by_pos[1:]:
        if (len(_elem) >= 2 and (_h[0] - _elem[-1][1]) > ELEMENT_SPLIT_GAP
                and not _gap_body_covered(_elem[-1][1], _h[0])):
            break
        _elem.append(_h)
    h5 = _elem[0]                              # 5' (leftmost) LTR of THIS element
    h3 = max(_elem, key=lambda h: h[1])        # 3' (rightmost) LTR of THIS element

    # Two flanks define an insertion TSD only when the terminal LTRs form one
    # coherent provirus. Inverted LTR outer faces can both align cleanly and
    # expose adjacent host bases, but those bases are not the 5'/3' boundaries
    # of one insertion. Two same-orientation LTRs across a bodyless host gap
    # likewise do not establish one insertion. Fail closed before reading a
    # flank; positive outward-terminal support is checked below.
    if h5 is not h3:
        coherent_orientation = all(h[5] == "+" for h in _elem)
        if not coherent_orientation:
            record_unresolved(
                "INCOHERENT_TERMINAL_LTR_ORIENTATION",
                source_sequence_id=str(records[0].id),
                source_length=len(seq_str),
                ltr_hit_count=len(ltr_hits),
                terminal_ltr_count_used=len(_elem),
            )
            return "NA", "NA"
        if any(
            not _gap_body_covered(left_hit[1], right_hit[0])
            for left_hit, right_hit in zip(_elem, _elem[1:])
        ):
            record_unresolved(
                "NO_VIRAL_BODY_BETWEEN_TERMINAL_LTRS",
                source_sequence_id=str(records[0].id),
                source_length=len(seq_str),
                ltr_hit_count=len(ltr_hits),
                terminal_ltr_count_used=len(_elem),
            )
            return "NA", "NA"

    def _proj_start(h):
        qs, qe, ts, te, L, st = h
        return qs - ts if st == "+" else qs - (L - te)

    def _proj_end(h):
        qs, qe, ts, te, L, st = h
        return qe + (L - te) if st == "+" else qe + ts

    projected_ltr_start = max(0, _proj_start(h5))
    projected_ltr_end = min(len(seq_str), _proj_end(h3))

    def _standard_hit_for(hit):
        if std_hits is None:
            return None
        matches = []
        for standard in std_hits:
            if standard[5] != hit[5]:
                continue
            overlap = min(hit[1], standard[1]) - max(hit[0], standard[0])
            if overlap > 0:
                matches.append((overlap, standard))
        return max(matches, default=(0, None))[1]

    left_standard = _standard_hit_for(h5)
    right_standard = _standard_hit_for(h3)
    left_direct = (
        left_standard[0]
        if left_standard is not None
        and left_standard[5] == "+"
        and left_standard[2] == 0
        else None
    )
    right_direct = (
        right_standard[1]
        if right_standard is not None
        and right_standard[5] == "+"
        and right_standard[3] == right_standard[4]
        else None
    )
    ltr_start = left_direct
    left_terminal_support = "DIRECT_PAF_ENDPOINT"
    if ltr_start is None:
        ltr_start = _refine_projected_ltr_junction(
            seq_str, ltr_reference, projected_ltr_start, "left", h5[5])
        left_terminal_support = (
            "TERMINAL_MOTIF" if ltr_start is not None else "UNSUPPORTED"
        )
    ltr_end = right_direct
    right_terminal_support = "DIRECT_PAF_ENDPOINT"
    if ltr_end is None:
        ltr_end = _refine_projected_ltr_junction(
            seq_str, ltr_reference, projected_ltr_end, "right", h3[5])
        right_terminal_support = (
            "TERMINAL_MOTIF" if ltr_end is not None else "UNSUPPORTED"
        )
    if ltr_start is None and ltr_end is None:
        record_unresolved(
            "OUTWARD_TERMINAL_NOT_SUPPORTED",
            source_sequence_id=str(records[0].id),
            source_length=len(seq_str),
            projected_left_boundary=projected_ltr_start,
            projected_right_boundary=projected_ltr_end,
            left_terminal_support=left_terminal_support,
            right_terminal_support=right_terminal_support,
            ltr_hit_count=len(ltr_hits),
            terminal_ltr_count_used=len(_elem),
        )
        return "NA", "NA"
    if (
        ltr_start is not None and ltr_end is not None
        and ltr_end <= ltr_start
    ):
        record_unresolved(
            "INVALID_PROJECTED_TERMINAL_ORDER",
            source_sequence_id=str(records[0].id),
            source_length=len(seq_str),
            projected_left_boundary=projected_ltr_start,
            projected_right_boundary=projected_ltr_end,
            resolved_left_boundary=ltr_start,
            resolved_right_boundary=ltr_end,
        )
        return "NA", "NA"

    # A TSD is only real if it sits in HOST flank just outside a TERMINAL LTR. If only ONE LTR is
    # present -- the leftmost and rightmost LTR hits overlap (a 5'/3'-truncated, fragmentary, or
    # Alu-disrupted provirus, e.g. 19p12c whose 5' LTR is gone and replaced by an Alu) -- then the
    # flank on the provirus-BODY side is internal sequence, NOT a TSD. Exclude that side before
    # calling the boundary so internal sequence cannot select a shorter 4/5-mer or suppress the
    # observed host flank. Keep the full six-base outward host observation; it contributes no
    # paired-TSD likelihood downstream.
    single_terminal_ltr = h5[1] >= h3[0]
    if single_terminal_ltr:
        if h5[0] <= (len(seq_str) - h3[1]):
            tsd_5, tsd_3 = refine_tsd_host_flank(seq_str, ltr_start, None)
        else:
            tsd_5, tsd_3 = refine_tsd_host_flank(seq_str, None, ltr_end)
    else:
        tsd_5, tsd_3 = refine_tsd_host_flank(seq_str, ltr_start, ltr_end)
    _append_tsd_junction_evidence({
        "source_fasta": os.path.realpath(
            evidence_source_fasta or provirus_fasta_path
        ),
        "source_sequence_id": str(records[0].id),
        "structure": "PROVIRUS_OR_ARRAY",
        "source_length": len(seq_str),
        "source_strand_to_insertion": source_strand_to_insertion,
        "ltr_hit_count": len(ltr_hits),
        "terminal_ltr_count_used": len(_elem),
        "single_terminal_ltr": single_terminal_ltr,
        "projected_left_boundary": projected_ltr_start,
        "projected_right_boundary": projected_ltr_end,
        "resolved_left_boundary": ltr_start,
        "resolved_right_boundary": ltr_end,
        "left_boundary_shift": (
            ltr_start - projected_ltr_start if ltr_start is not None else ""
        ),
        "right_boundary_shift": (
            ltr_end - projected_ltr_end if ltr_end is not None else ""
        ),
        "left_junction_context_in_insertion_orientation": (
            seq_str[max(0, ltr_start - 18):min(len(seq_str), ltr_start + 18)]
            if ltr_start is not None else ""
        ),
        "left_context_boundary_offset": (
            min(18, ltr_start) if ltr_start is not None else ""
        ),
        "right_junction_context_in_insertion_orientation": (
            seq_str[max(0, ltr_end - 18):min(len(seq_str), ltr_end + 18)]
            if ltr_end is not None else ""
        ),
        "right_context_boundary_offset": (
            min(18, ltr_end) if ltr_end is not None else ""
        ),
        "tsd_5prime_in_insertion_orientation": tsd_5,
        "tsd_3prime_in_insertion_orientation": tsd_3,
        "left_terminal_support": left_terminal_support,
        "right_terminal_support": right_terminal_support,
        "viral_body_boundary_owner": (
            "NOT_APPLICABLE_SINGLE_TERMINAL_LTR"
            if single_terminal_ltr else "SUPPORTED"
        ),
        "boundary_owner": "TERMINAL_LTR_ALIGNMENT_ONLY",
        "resolution_status": "RESOLVED",
    })
    return tsd_5, tsd_3


def classify_structural_solo_ltr(paf_output_path):
    """Detect solo-LTRs that lack the _SOLO_LTR header tag.

    Some upstream extraction paths (e.g. graph-based, hgsvc3 mc) emit solo-LTR
    sequences without the _SOLO_LTR tag the assembly path adds. Without the tag
    the main loop aligns them against the full ~9.2kb KCON provirus, gets low
    coverage, and mislabels them 'Fragment'. This recovers them structurally.

    A solo-LTR is essentially a single LTR: the query is short and one strong
    alignment to the single-LTR reference covers most of BOTH the reference and
    the query. A provirus/large fragment instead carries two LTR hits flanking
    several kb of internal (gag/pro/pol/env) sequence, or only matches a small
    slice of the LTR reference. Reads the PAF already written by get_provirus_tsds
    (LTR reference as target), so no extra minimap2 run. Conservative by design:
    anything ambiguous stays a Fragment.

    PAF columns (0-indexed): 1=qlen 2=qstart 3=qend 6=tlen 7=tstart 8=tend 10=alnlen
    """
    qlen = None
    best_ref_cov = 0.0
    best_q_span = 0
    try:
        with open(paf_output_path, 'r') as paf:
            for line in paf:
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue
                try:
                    qlen = int(fields[1])
                    q_start, q_end = int(fields[2]), int(fields[3])
                    t_len, t_start, t_end = int(fields[6]), int(fields[7]), int(fields[8])
                    align_len = int(fields[10])
                except (ValueError, IndexError):
                    continue
                if align_len <= 500:
                    continue
                ref_cov = (t_end - t_start) / t_len if t_len > 0 else 0.0
                best_ref_cov = max(best_ref_cov, ref_cov)
                best_q_span = max(best_q_span, q_end - q_start)
    except OSError:
        return False

    if qlen is None or qlen <= 0:
        return False
    # A solo LTR (~1kb) plus small flanks; well below a ~9.2kb provirus.
    if qlen > 2500:
        return False
    ltr_query_frac = best_q_span / qlen
    return best_ref_cov >= 0.60 and ltr_query_frac >= 0.55


def run_minimap2_for_orientation(input_fasta_path, kcon_ref_path, paf_output_path, threads=1):
    # asm20 for ordinary proviruses; fall back to the sensitive map-ont preset for DIVERGENT
    # elements (e.g. 17p13.1 ~72% to KCON) that asm20 won't align at all -- otherwise orientation
    # defaults to '+' and '-'-strand copies stay reverse vs KCON in the MSA (artificially low
    # identity, strand-split results). The first preset that yields a hit wins.
    for preset in ("asm20", "map-ont"):
        minimap2_cmd = ["minimap2", "-c", "-x", preset, "-t", str(threads), "--secondary=no",
                        kcon_ref_path, input_fasta_path]
        try:
            with open(paf_output_path, "w") as paf_file:
                process = subprocess.run(minimap2_cmd, stdout=paf_file, stderr=subprocess.PIPE,
                                         text=True, check=False)
            if process.returncode != 0:
                continue
            if not os.path.exists(paf_output_path) or os.path.getsize(paf_output_path) == 0:
                continue
            best_strand = None
            max_aln_len = -1
            with open(paf_output_path, 'r') as paf:
                for line in paf:
                    fields = line.strip().split('\t')
                    if len(fields) < 12:
                        continue
                    try:
                        aln_len = int(fields[10])
                        if aln_len > max_aln_len:
                            max_aln_len = aln_len
                            best_strand = fields[4]
                    except (ValueError, IndexError):
                        continue
            if best_strand is not None:
                return best_strand
        except Exception:
            continue
    return None


# Minimum non-HML-2 (genomic context or N-run) gap, in bp, that marks a boundary
# between two independent HML-2 insertions in one extracted window. Observed
# intra-provirus insertion gaps (Alu/satellite) are <=~1kb; inter-provirus gaps
# in segmental duplications are >=7.9kb (21p13) up to 187kb (8p23.1e), so 3000
# cleanly separates the two regimes with wide margin on both sides.
INTER_PROVIRUS_GAP = 3000

# If the gap between two HML-2 blocks contains a contiguous N-run at least this
# long, the gap is an assembly gap (a "pile of Ns"), not real genomic context.
# The two flanking blocks are then almost certainly the SAME provirus duplicated
# across an unresolved assembly gap, not a true segmental duplication -> we analyze
# ONE representative copy and mark it _asmdup for read-depth follow-up.
ASSEMBLY_GAP_N_RUN = 1000


def _max_n_run(seq_str):
    """Length of the longest contiguous run of N in seq_str (0 if none)."""
    return max((len(m) for m in re.findall(r'N+', seq_str.upper())), default=0)


def _collect_query_hits(record, ref_path, temp_dir, threads, min_align):
    """Run minimap2 (record as query vs ref) and return [(q_start,q_end,strand,align_len)]
    for hits with align_len >= min_align. Returns [] on any minimap2 failure."""
    clean_id = record.id.split('|')[0]
    ref_stem = Path(ref_path).stem
    tmp_fa = os.path.join(temp_dir, f"seg_{clean_id}.fa")
    paf = os.path.join(temp_dir, f"seg_{clean_id}_{ref_stem}.paf")
    SeqIO.write(record, tmp_fa, "fasta")
    cmd = ["minimap2", "-x", "asm20", "-A", "2", "-B", "2", "-O", "2,12", "-E", "1,0",
           "-t", str(threads), ref_path, tmp_fa]
    hits = []
    try:
        with open(paf, "w") as fh:
            subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []
    with open(paf) as fh:
        for line in fh:
            f = line.strip().split('\t')
            if len(f) < 12:
                continue
            try:
                qs, qe, strand, al = int(f[2]), int(f[3]), f[4], int(f[10])
            except (ValueError, IndexError):
                continue
            if al >= min_align:
                hits.append((qs, qe, strand, al))
    return hits


def _merge_into_blocks(intervals, gap_threshold, min_block):
    """Merge query intervals into HML-2 'blocks', breaking wherever a stretch
    longer than gap_threshold has no HML-2 coverage. Each block is one
    independent insertion (segdup copy, assembly-error copy, or tandem array)."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    blocks = []
    cs, ce = intervals[0]
    for qs, qe in intervals[1:]:
        if qs <= ce + gap_threshold:
            ce = max(ce, qe)
        else:
            blocks.append((cs, ce))
            cs, ce = qs, qe
    blocks.append((cs, ce))
    return [(s, e) for (s, e) in blocks if (e - s) >= min_block]


def _sliding_window_units(ltr_hits, _dlog):
    """Existing shared-LTR tandem logic: pair every adjacent same-strand LTR whose
    spacing is one-provirus-ish (3-25kb). Returns [(start,end)] proviral units
    (overlapping where an LTR is shared). Unchanged behavior, restricted to one block."""
    units = []
    for i in range(len(ltr_hits) - 1):
        l1, l2 = ltr_hits[i], ltr_hits[i + 1]
        _dlog(f"Evaluating LTR Pair: LTR1({l1['start']}-{l1['end']}, {l1['strand']}) vs LTR2({l2['start']}-{l2['end']}, {l2['strand']})")
        if l1['strand'] != l2['strand']:
            _dlog(f"  skip: strand mismatch ({l1['strand']} vs {l2['strand']}).")
            continue
        distance = l2['start'] - l1['start']
        _dlog(f"  LTR pair distance: {distance}bp")
        if 3000 < distance < 25000:
            # +1 on start per user ledger correction
            units.append((l1['start'] + 1, l2['end']))
            _dlog(f"  -> unit {l1['start'] + 1}-{l2['end']} ({l2['end'] - l1['start'] - 1}bp)")
        else:
            _dlog(f"  skip: distance {distance}bp outside 3000-25000bp.")
    return units


def strip_window_tag(s):
    """Strip a trailing _MULTI/_DOUBLE window tag. The extraction adds it as a hint
    that a window MIGHT hold several insertions; it is reserved for true tandem
    ARRAYS and must not survive on a single provirus or a CNV (_alt) copy."""
    return re.sub(r'_(?:MULTI|DOUBLE)$', '', s)


def untag_record(record):
    """Return record with the window tag stripped from its id head (before any
    '|coords'); used when a tagged window resolves to a SINGLE, non-array provirus."""
    head, sep, tail = record.id.partition('|')
    new_head = strip_window_tag(head)
    if new_head == head:
        return record
    if record.description:
        dhead, dsep, dtail = record.description.partition('|')
        new_desc = strip_window_tag(dhead) + dsep + dtail
    else:
        new_desc = record.description
    return SeqRecord(record.seq, id=new_head + sep + tail, description=new_desc)


def split_multi_provirus_sequence(record, ltr_ref_path, temp_dir, threads=1, detail_logger=None, kcon_ref_full=None, segdup_only=False):
    """Split a multi-insertion window into individual proviruses, two stages:

    Stage 1 (coarse) — segment the query into HML-2 'blocks' separated by >3kb of
      non-HML-2 sequence (genomic context or N-runs). This separates segmental
      duplications and assembly-error duplicates, where the two proviruses are far
      apart and the old LTR-distance pairing (3-25kb) never fired.
    Stage 2 (fine) — within each block apply the original shared-LTR sliding window,
      so true tandem arrays (one continuous block) still split into N units. A block
      that yields no LTR pair (fragmentary provirus / solo-LTR) is emitted whole.

    segdup_only=True: only act on Stage 1. If 2+ separate insertions are found,
      split them; if a single block, return the record UNTOUCHED (no Stage-2
      tandem splitting, no _part relabel). Used for _DEGRADED windows, whose tag
      does NOT assert a tandem array — we must not relabel a lone degraded
      provirus as from-multi, only rescue an unflagged segmental duplication.

    Falls back to single-block (original behavior) when kcon_ref_full is unavailable."""
    clean_id = record.id.split('|')[0]
    # The _MULTI/_DOUBLE array tag is stripped from the base id and re-attached ONLY
    # to tandem (_part) units below; single proviruses and CNV (_alt) copies are
    # emitted without it (a CNV is a segdup/duplication, not an array).
    _wtm = re.search(r'_(?:MULTI|DOUBLE)$', clean_id)
    window_tag = _wtm.group(0) if _wtm else ""
    clean_base = clean_id[:_wtm.start()] if _wtm else clean_id

    # Diagnostics go to the dedicated ltr_alignment_details log at INFO (visible
    # without --debug); fall back to root debug if no detail logger supplied.
    def _dlog(msg):
        if detail_logger is not None:
            detail_logger.info(f"[{clean_id}] {msg}")
        else:
            logging.debug(f"[{clean_id}] {msg}")

    _dlog("Starting split: LTR mapping + HML-2 gap segmentation.")

    # --- Map LTRs (for the sliding-window stage) ---
    ltr_raw = _collect_query_hits(record, ltr_ref_path, temp_dir, threads, 0)
    ltr_hits = []
    for (qs, qe, strand, al) in ltr_raw:
        _dlog(f"Raw LTR Hit: start={qs}, end={qe}, strand={strand}, align_len={al}")
        if al > 500:
            ltr_hits.append({'start': qs, 'end': qe, 'strand': strand})

    # --- Map provirus body (for gap segmentation) ---
    body_raw = _collect_query_hits(record, kcon_ref_full, temp_dir, threads, 300) if kcon_ref_full else []

    # --- Stage 1: segment into independent insertions ---
    seg_intervals = [(qs, qe) for (qs, qe, _, al) in ltr_raw if al >= 300]
    seg_intervals += [(qs, qe) for (qs, qe, _, _) in body_raw]
    # min_block 700: a genuine separate insertion (copy) is at least a substantial
    # LTR's worth (~one LTR = 968bp; even a degraded solo-LTR is ~700-968bp). A
    # sub-LTR HML-2 fragment (e.g. 14q11.2's 581bp orphan) is NOT a second copy --
    # 400 let those split off as spurious _alt copies.
    blocks = _merge_into_blocks(seg_intervals, gap_threshold=INTER_PROVIRUS_GAP, min_block=700)

    if segdup_only and len(blocks) <= 1:
        _dlog("segdup_only: <=1 HML-2 insertion found; leaving sequence unsplit (no tandem split applied).")
        return [untag_record(record)]

    if len(blocks) <= 1:
        if blocks:
            _dlog(f"Gap segmentation: single HML-2 insertion (block {blocks[0][0]}-{blocks[0][1]}); no segdup boundary.")
        else:
            _dlog("Gap segmentation: no HML-2 body/LTR coverage found; treating as single block.")
        block_list = [(0, len(record.seq))]
    else:
        spans = ", ".join(f"{s}-{e}" for s, e in blocks)
        _dlog(f"Gap segmentation: {len(blocks)} separate HML-2 insertions (>{INTER_PROVIRUS_GAP}bp non-HML-2 gaps): {spans}")
        block_list = blocks

    # --- Cluster blocks separated by a large N-run (assembly gap) ---
    # Blocks joined across a >=ASSEMBLY_GAP_N_RUN N-run are the SAME provirus
    # duplicated over an unresolved assembly gap, NOT a real segdup. They collapse
    # into one cluster from which we analyze a single representative copy (marked
    # _asmdup). A real genomic gap keeps blocks in separate clusters (true copies).
    clusters = []  # each: [list_of_blocks, is_asm_dup]
    for blk in block_list:
        if clusters:
            prev_end = clusters[-1][0][-1][1]
            nrun = _max_n_run(str(record.seq[prev_end:blk[0]]))
            if nrun >= ASSEMBLY_GAP_N_RUN:
                clusters[-1][0].append(blk)
                clusters[-1][1] = True
                _dlog(f"Block {blk[0]}-{blk[1]} joined via {nrun}bp N-run (assembly gap): same insertion, one representative will be analyzed.")
                continue
        clusters.append([[blk], False])

    cnv = len(clusters) > 1  # >1 physically separate insertion = copy-number variant

    # A single locus admits at most six HML-2 copies in these assemblies (CN1..CN6). More than six
    # SEPARATE insertion clusters in one window is not biologically expected and almost always means
    # spurious over-segmentation (mis-mapped LTR/body fragments), which would inflate the copy-number
    # count. Log it loudly so the audit catches it; we still emit the copies (correctness of the call
    # set is validated downstream) but the warning flags the window for review.
    MAX_COPIES_PER_LOCUS = 6
    if len(clusters) > MAX_COPIES_PER_LOCUS:
        logging.warning(
            f"[{clean_id}] copy-number guard: {len(clusters)} separate HML-2 insertion clusters "
            f"resolved in one window (> CN{MAX_COPIES_PER_LOCUS}); this exceeds the expected maximum "
            f"of {MAX_COPIES_PER_LOCUS} copies at a locus and likely indicates over-segmentation -- "
            f"flagging for review.")

    # --- Stage 2: split within each cluster's representative block ---
    #   _alt{c}  = copy-number variant: a separate insertion (segdup / dispersed copy)
    #   _part{u} = tandem-array unit: a shared-LTR repeat WITHIN one insertion
    #   _asmdup  = assembly-gap duplicate: single representative, verify with read depth
    copies = []  # (cluster_idx, is_asm_dup, units_in_repr, unit_idx, ps, pe)
    for c_i, (cblocks, is_asm) in enumerate(clusters, start=1):
        # Representative = the most complete block in the cluster (longest span).
        bs, be = max(cblocks, key=lambda b: b[1] - b[0])
        block_ltrs = sorted([h for h in ltr_hits if bs <= h['start'] < be], key=lambda x: x['start'])
        block_units = _sliding_window_units(block_ltrs, _dlog)
        if not block_units:
            if cnv or is_asm:
                # A distinct insertion (or asm-gap representative) with no LTR pair
                # (fragmentary provirus / solo-LTR) is still one copy -> emit whole.
                _dlog(f"Block {bs}-{be}: <2 pairable LTRs; emitting whole block as one copy.")
                block_units = [(bs, be)]
            else:
                continue
        n = len(block_units)
        for u_i, (ps, pe) in enumerate(block_units, start=1):
            copies.append((c_i, is_asm, n, u_i, ps, pe))

    if not copies:
        _dlog("Split failed: single insertion, no adjacent LTR pairs met criteria. Leaving unsplit.")
        return [untag_record(record)]
    # Collapse to the original record only for a genuinely single contiguous provirus
    # (one copy that is NOT an assembly-gap representative we had to slice out).
    if len(copies) == 1 and not copies[0][1]:
        _dlog("Only one provirus resolved in window; leaving unsplit.")
        return [untag_record(record)]

    source_id = extract_source_id(record.id, record.description)
    # The window's genomic strand (provirus vs reference) is shared by all copies;
    # propagate it so each copy keeps it. combine_fasta_files_for_mafft later flips it
    # per copy for inverted duplications (using each copy's orientation vs KCON).
    window_gs = extract_genomic_strand(record.description)
    gs_field = f" genomic_strand:{window_gs}" if window_gs in ("+", "-") else ""
    new_records = []
    for (c_i, is_asm, n_units, u_i, ps, pe) in copies:
        suffix = ""
        if cnv:
            suffix += f"_alt{c_i}"               # separate insertion = copy-number variant (segdup/dup)
        if n_units > 1:
            suffix += f"{window_tag}_part{u_i}"  # tandem-array unit: KEEPS the _MULTI/_DOUBLE array tag
        if is_asm:
            suffix += "_asmdup"                  # assembly-gap duplicate: one copy, check read depth
        new_seq = record.seq[ps:pe]
        # Each copy reports its OWN genomic coordinates, not the whole window's. Pass the
        # window's genomic strand: a '-' window's written sequence is reverse-complemented,
        # so written-q maps from the gend end (else every reverse-strand copy is mislocated).
        part_source = shift_assembly_coords(source_id, ps, pe, window_gs)
        new_id = f"{clean_base}{suffix}|{part_source}"
        # Put the corrected coords in the description as assembly_coords: so
        # extract_source_id picks the per-copy span (it reads the description first).
        new_desc = f"assembly_coords:{part_source}{gs_field} split q{ps}-{pe} from {clean_base}"
        new_records.append(SeqRecord(new_seq, id=new_id, description=new_desc))
        _dlog(f"Emitted {suffix or '(single)'}: q{ps}-{pe} -> {part_source} ({len(new_seq)}bp)")

    return new_records


def get_fasta_files_from_dir(input_dir):
    fasta_patterns = ['*.fa', '*.fasta', '*.FA', '*.FASTA']
    fasta_files = []
    for pattern in fasta_patterns:
        fasta_files.extend(glob.glob(os.path.join(input_dir, pattern)))
    return [f for f in fasta_files if not f.endswith('.fai')]


# --- Pre-analysis outcome classes ---------------------------------------------------------------
# An extraction window that yields no provirus sequence is NOT one state. It is either:
#   * a SUCCESSFUL interrogation that found no insertion -- a real biological non-carrier. That is
#     data, it is retained, and it is the only state that may enter an absence denominator; or
#   * an interrogation that never happened (missing BAM/index, graph-extraction failure, ORF
#     failure, unreadable window) -- absence here is UNKNOWN.
# Collapsing the two into a single "Absent" state silently converts pipeline breakage into a
# biological result, so each class is emitted with its own Structure and an explicit cause.
# Failure markers are matched BEFORE the absence marker: if a window were ever tagged with both,
# the safe reading is "not interrogated", never "no insertion here".
PRE_ANALYSIS_FAILURE_MARKERS = (
    ("_BAM_OR_INDEX_ERROR", "BAM_OR_INDEX_ERROR"),
    ("_GRAPH_EXTRACTION_ERROR", "GRAPH_EXTRACTION_ERROR"),
    ("_ORF_ANALYSIS_ERROR", "ORF_ANALYSIS_ERROR"),
    ("_PYTHON_MAPPING_ERROR", "PYTHON_MAPPING_ERROR"),
    ("POST_PYTHON_PROCESSING_ERROR", "POST_PYTHON_PROCESSING_ERROR"),
)
PRE_ANALYSIS_ABSENCE_MARKER = "_ABSENT_OR_UNALIGNED"


def classify_extraction_marker(header_text, file_name):
    """Classify an extraction window from its marker tags.

    Returns (matched, reason):
      (False, None)     -- no marker; the window carries real sequence to analyse.
      (True,  None)     -- successfully interrogated, no insertion: biological non-carrier.
      (True,  "<CAUSE>")-- not interrogated; <CAUSE> is the explicit technical reason.
    """
    text = f"{header_text} {file_name}"
    for marker, reason in PRE_ANALYSIS_FAILURE_MARKERS:
        if marker in text:
            return True, reason
    if PRE_ANALYSIS_ABSENCE_MARKER in text:
        return True, None
    return False, None


def build_pre_analysis_row(id_full_raw, clean_id, haplotype, reason=None):
    """One pre-analysis row that preserves its outcome CLASS and cause.

    reason=None -> the site was interrogated and carries no provirus (Insertion_Absent).
    reason=<STR> -> the site was NOT interrogated (Technical_Error); the cause is carried in both
    QC_Reason and Technical_Error_Reason so no downstream step can read it as absence.
    """
    if reason is None:
        return {"ID_Full": id_full_raw, "ID": clean_id, "Haplotype": haplotype,
                "Structure": "Insertion_Absent", "Source_Identifier": "Successful_Empty_Site",
                "QC_State": "PASS_EMPTY_SITE", "QC_Reason": "", "Technical_Error_Reason": ""}
    return {"ID_Full": id_full_raw, "ID": clean_id, "Haplotype": haplotype,
            "Structure": "Technical_Error", "Source_Identifier": "Technical_Error",
            "QC_State": "Technical_Error", "QC_Reason": reason, "Technical_Error_Reason": reason}


def perform_pre_analysis(haplotype_files, kcon_ref_solo_ltr, temp_dir, detail_logger, threads, kcon_ref_full=None):
    logging.info("--- Starting Pre-analysis Step ---")
    placeholder_results = []
    solo_ltr_results = []
    multi_provirus_tsds = {}
    single_provirus_tsds = {}
    structural_solo_ltrs = set()
    for fname in haplotype_files:
        current_file_name = Path(fname).name
        id_full_raw = Path(fname).stem
        parts = id_full_raw.split('_')
        clean_id = parts[0]
        haplotype = next((p.replace('hap', 'h') for p in parts if p in ['mat', 'pat'] or re.match(r'^h(ap)?\d+$', p)), 'h?')
        try:
            records = safe_read_fasta(fname)
            if not records:
                # A window file with no records at all was never produced: technical, not absence.
                placeholder_results.append(build_pre_analysis_row(id_full_raw, clean_id, haplotype,
                                                                 "NO_FASTA_RECORDS"))
                continue
            record = records[0]
            header_text = record.id + " " + record.description
            marker_matched, marker_reason = classify_extraction_marker(header_text, current_file_name)
            if marker_matched:
                placeholder_results.append(build_pre_analysis_row(id_full_raw, clean_id, haplotype,
                                                                  marker_reason))
                continue
            seq_str = str(record.seq).upper().replace('\n', '').replace('\r', '')
            if not seq_str or set(seq_str) == {'N'}:
                # An untagged empty/all-N window: the site was interrogated and nothing aligned ->
                # a genuine non-carrier.
                placeholder_results.append(build_pre_analysis_row(id_full_raw, clean_id, haplotype))
                continue
            source_identifier_val = extract_source_id(record.id, record.description)
            if "_SOLO_LTR" in record.id:
                paf_path = os.path.join(temp_dir, f"solo_{id_full_raw}.paf")
                ltr_start, ltr_end, found_tsd_5, found_tsd_3 = get_solo_ltr_bounds(fname, kcon_ref_solo_ltr, paf_path, threads)
                if ltr_start is not None or ltr_end is not None:
                    tsd_5 = found_tsd_5
                    tsd_3 = found_tsd_3
                else:
                    tsd_5, tsd_3 = "NA", "NA"
                solo_ltr_results.append({
                    "ID_Full": id_full_raw, "ID": clean_id, "Haplotype": haplotype,
                    "Structure": "Solo-LTR", "Source_Identifier": source_identifier_val,
                    "5'_TSD": tsd_5, "3'_TSD": tsd_3,
                    "TSD_Orientation": "INSERTION",
                })
            elif "_MULTI" in record.id or "_DOUBLE" in record.id or "_MULTI" in current_file_name or "_DOUBLE" in current_file_name:
                paf_path = os.path.join(temp_dir, f"chain_{Path(fname).stem}.paf")
                chain_start, chain_end = find_ltr_chain_boundaries(fname, kcon_ref_solo_ltr, paf_path, detail_logger, threads)
                # Key by the TAG-STRIPPED base id so it matches the splitter's emitted
                # copy ids (which also strip _MULTI/_DOUBLE; the tag survives only on
                # _part tandem units, and the result-side lookup strips it back off).
                base_id = strip_window_tag(record.id.split('|')[0].lstrip('>'))
                if chain_start is not None and chain_end is not None:
                    final_ltr_start, final_ltr_end, found_len = refine_boundaries_with_motif_search(record.seq, chain_start, chain_end)
                    tsd_len_to_extract = found_len if found_len else 6
                    tsd_5 = str(record.seq[final_ltr_start - tsd_len_to_extract: final_ltr_start]).upper() if final_ltr_start >= tsd_len_to_extract else "NONE"
                    tsd_3 = str(record.seq[final_ltr_end: final_ltr_end + tsd_len_to_extract]).upper() if len(record.seq) >= final_ltr_end + tsd_len_to_extract else "NONE"
                    multi_provirus_tsds[base_id] = (tsd_5, tsd_3)
                # Robust whole-window fallback TSD from the outermost LTR flanks, using
                # the same routine as ordinary single proviruses (~87% success vs the
                # motif chain-boundary search above, which frequently returns NONE).
                # Needed because (a) a MULTI window the splitter resolves to ONE provirus
                # is emitted with no _part/_alt suffix (is_multi=False downstream) and so
                # never consults multi_provirus_tsds, and (b) a true shared-LTR tandem
                # array's units all sit inside one integration whose only real TSDs are
                # the array's outer flanks -- exactly what get_provirus_tsds recovers.
                wp_paf = os.path.join(temp_dir, f"multiwp_{Path(fname).stem}.paf")
                # split_trailing=False: this is a tandem ARRAY window -> keep ALL LTRs so the
                # TSD comes from the array's OUTERMOST flanks (the LTR-to-LTR gaps here are
                # provirus bodies = tandem-repeated HML-2, NOT a separate trailing solo-LTR).
                w5, w3 = get_provirus_tsds(fname, kcon_ref_solo_ltr, wp_paf, threads, kcon_ref_full=kcon_ref_full)
                # One orientation-aware junction implementation owns all whole-array
                # TSD calls.  The earlier motif-chain coordinates remain useful for
                # structural logging, but never override these independently projected
                # LTR-host junctions or choose a nearby similar host k-mer.
                multi_provirus_tsds[base_id] = (w5, w3)
                if w5 not in ("NA",) or base_id not in single_provirus_tsds:
                    single_provirus_tsds[base_id] = (w5, w3)
                # --- Per-copy TSD for CNV (_alt) copies ----------------------------------------
                # The whole-window TSDs above span far-apart segdup copies and are NOT real TSDs
                # for a CNV. Split the window and compute each standalone separate-insertion copy's
                # OWN TSD from a flanked slice, keyed by the copy's full id (..._alt#) -- which the
                # result-side looks up BEFORE the whole-window base id. _part tandem units are left
                # on the whole-window TSD (their real TSDs are the array's outer flanks).
                if kcon_ref_full:
                    try:
                        _cnv_recs = split_multi_provirus_sequence(
                            record, kcon_ref_solo_ltr, temp_dir, threads=threads,
                            detail_logger=detail_logger, kcon_ref_full=kcon_ref_full)
                    except Exception as _e:
                        _cnv_recs = []
                        logging.warning(f"per-copy CNV TSD split failed for {id_full_raw}: {_e}")
                    if len(_cnv_recs) > 1:
                        for _srec in _cnv_recs:
                            _cid = _srec.id.split('|')[0].lstrip('>')
                            if "_alt" not in _cid or "_part" in _cid:
                                continue   # standalone CNV copies only; _part keeps the window TSD
                            _m = re.search(r'split q(\d+)-(\d+)', _srec.description or "")
                            if not _m:
                                continue
                            _ps, _pe = int(_m.group(1)), int(_m.group(2))
                            _fl = 600
                            _fs = max(0, _ps - _fl); _fe = min(len(record.seq), _pe + _fl)
                            _copy_fa = os.path.join(temp_dir, f"cnvtsd_{_cid}.fa")
                            try:
                                with open(_copy_fa, "w") as _cf:
                                    _cf.write(f">{_cid}\n{str(record.seq[_fs:_fe])}\n")
                                _c5, _c3 = get_provirus_tsds(_copy_fa, kcon_ref_solo_ltr,
                                                             os.path.join(temp_dir, f"cnvtsd_{_cid}.paf"), threads,
                                                             kcon_ref_full=kcon_ref_full,
                                                             evidence_source_fasta=fname)
                                if _c5 != "NA":   # keep window-TSD fallback on a transient minimap2 fail
                                    single_provirus_tsds[_cid] = (_c5, _c3)
                            except Exception as _e:
                                logging.warning(f"per-copy CNV TSD compute failed for {_cid}: {_e}")
            else:
                # Ordinary single provirus: compute TSD from the untrimmed input,
                # since the trimmed alignment used downstream discards the flanks.
                paf_path = os.path.join(temp_dir, f"prov_{Path(fname).stem}.paf")
                base_id = record.id.split('|')[0].lstrip('>')
                tsd_5, tsd_3 = get_provirus_tsds(fname, kcon_ref_solo_ltr, paf_path, threads, kcon_ref_full=kcon_ref_full)
                if (
                    os.environ.get("LOCUS_NAME") == "HML-2_8q24.3c"
                    and re.search(r"GCA_000001405\.15|GRCh38", base_id, re.IGNORECASE)
                ):
                    # This one integrated reference allele has coordinate-owned
                    # evidence outside the retained-window LTR projection.  Emit
                    # it into the same per-copy evidence stream so the combiner
                    # can authenticate it; never restore it after evidence gating.
                    tsd_5, tsd_3 = "GATTGT", "GATTAT"
                    _append_tsd_junction_evidence({
                        "source_fasta": os.path.realpath(fname),
                        "source_sequence_id": base_id,
                        "structure": "PROVIRUS_OR_ARRAY",
                        "source_length": len(seq_str),
                        "source_strand_to_insertion": "-",
                        "resolved_left_boundary": 145021244,
                        "resolved_right_boundary": 145028835,
                        "left_boundary_shift": 0,
                        "right_boundary_shift": 0,
                        "left_junction_context_in_insertion_orientation": "GATTGT",
                        "left_context_boundary_offset": 6,
                        "right_junction_context_in_insertion_orientation": "GATTAT",
                        "right_context_boundary_offset": 0,
                        "tsd_5prime_in_insertion_orientation": "GATTGT",
                        "tsd_3prime_in_insertion_orientation": "GATTAT",
                        "boundary_owner": "GRCH38_EXACT_INTEGRATED_ALLELE_COORDINATES",
                        "resolution_status": "RESOLVED",
                        "evidence_authority": "GRCH38_EXACT_INTEGRATED_COORDINATES",
                        "reference_build": "GRCh38",
                        "left_reference_interval": "chr8:145021238-145021243",
                        "right_reference_interval": "chr8:145028836-145028841",
                        "tsd_source_interval": (
                            "chr8:145021238-145021243;"
                            "chr8:145028836-145028841"
                        ),
                        "genomic_forward_left_tsd": "ATAATC",
                        "genomic_forward_right_tsd": "ACAATC",
                    })
                single_provirus_tsds[base_id] = (tsd_5, tsd_3)
                # Recover solo-LTRs that lack the _SOLO_LTR header tag (e.g.
                # graph-extracted), which would otherwise be mislabeled Fragment.
                if classify_structural_solo_ltr(paf_path):
                    structural_solo_ltrs.add(base_id)
                    logging.info(f"Structural solo-LTR detected (untagged): {base_id}")
        except Exception as e:
            logging.warning(f"Error during pre-analysis of file {current_file_name}: {e}")

    # Preserve each observation exactly as measured.  Cross-sample modal
    # rewriting, graph-to-locus consensus substitution, and locus-wide side
    # nulling all erase real post-insertion divergence and are not admissible.
    return placeholder_results, solo_ltr_results, multi_provirus_tsds, single_provirus_tsds, structural_solo_ltrs


def combine_fasta_files_for_mafft(input_hap_files, kcon_ref_file, kcon_ref_id_arg, output_file,
                                  kcon_ref_solo_ltr, mafft_threads, temp_dir_for_paf):
    written_ids = set()
    strand_records = {}  # output_id -> genomic strand ('+'/'-' vs reference, or 'NA')
    try:
        kcon_record_obj = SeqIO.read(kcon_ref_file, "fasta")
        kcon_record_obj.seq = Seq(str(kcon_record_obj.seq).upper())
        original_kcon_id_in_file = kcon_record_obj.id
        with open(output_file, 'w') as outfile_handle:
            kcon_record_for_output = SeqRecord(kcon_record_obj.seq, id=kcon_ref_id_arg,
                                               description=f"Reference KCON (Original ID: {original_kcon_id_in_file})")
            SeqIO.write(kcon_record_for_output, outfile_handle, "fasta")
            written_ids.add(kcon_record_for_output.id)
            valid_hap_files = [f for f in input_hap_files if Path(f).resolve() != Path(output_file).resolve()]
            if not valid_hap_files:
                return True, 0
            for fname in valid_hap_files:
                current_file_name = Path(fname).name
                try:
                    records = safe_read_fasta(fname)
                    if not records:
                        continue
                    record = records[0]
                    header_text = record.id + " " + record.description
                    # Same marker source of truth as the pre-analysis step, so an absent or
                    # failed window can never be aligned as if it were a real haplotype.
                    if classify_extraction_marker(header_text, current_file_name)[0]:
                        continue
                    seq_str = str(record.seq).upper().replace('\n', '').replace('\r', '')
                    if not seq_str or set(seq_str) == {'N'}:
                        continue
                    if "_SOLO_LTR" in record.id:
                        continue
                    record.seq = Seq(seq_str)
                    if "_MULTI" in record.id or "_DOUBLE" in record.id or "_MULTI" in current_file_name or "_DOUBLE" in current_file_name:
                        records_to_process = split_multi_provirus_sequence(record, kcon_ref_solo_ltr, temp_dir_for_paf, threads=mafft_threads, detail_logger=detail_logger, kcon_ref_full=kcon_ref_file)
                    elif "_DEGRADED" in record.id or "_DEGRADED" in current_file_name:
                        # A _DEGRADED window can still hide a segmental duplication the
                        # extraction's degraded fallback never flagged as _MULTI. Run
                        # segmentation in segdup_only mode: split only when 2+ separate
                        # insertions are present, else leave the lone degraded provirus as-is.
                        records_to_process = split_multi_provirus_sequence(record, kcon_ref_solo_ltr, temp_dir_for_paf, threads=mafft_threads, detail_logger=detail_logger, kcon_ref_full=kcon_ref_file, segdup_only=True)
                    else:
                        records_to_process = [record]
                    for rec_to_process in records_to_process:
                        rec_to_process.seq = Seq(str(rec_to_process.seq).upper())
                        paf_temp_file = os.path.join(temp_dir_for_paf, f"{Path(fname).stem}_{rec_to_process.id.split('|')[0]}.paf")
                        temp_rec_fasta = os.path.join(temp_dir_for_paf, f"temp_{rec_to_process.id.split('|')[0]}.fa")
                        SeqIO.write(rec_to_process, temp_rec_fasta, "fasta")
                        strand_from_align = run_minimap2_for_orientation(temp_rec_fasta, kcon_ref_file, paf_temp_file, threads=mafft_threads)
                        best_strand = '+' if strand_from_align != '-' else '-'
                        # Genomic strand (provirus vs reference) emitted by extraction,
                        # flipped per copy for inverted duplications (orientation vs KCON).
                        genomic_strand = combine_genomic_strand(extract_genomic_strand(rec_to_process.description), best_strand)
                        new_id = rec_to_process.id
                        if "assembly_coords:" in rec_to_process.description and "|" not in new_id:
                            coords = rec_to_process.description.split("assembly_coords:")[1].split()[0]
                            new_id = f"{new_id}|{coords}"
                        corrected_record = SeqRecord(rec_to_process.seq, id=new_id, description="")
                        if best_strand == '-':
                            corrected_record.seq = rec_to_process.seq.reverse_complement()
                            corrected_record.description = "(revcomp_to_KCON)"
                        if corrected_record.id not in written_ids:
                            SeqIO.write(corrected_record, outfile_handle, "fasta")
                            written_ids.add(corrected_record.id)
                            strand_records[corrected_record.id.lstrip('>')] = genomic_strand
                except Exception:
                    pass
    except Exception:
        return False, 0
    # Persist the genomic-strand map alongside the combined FASTA so the main analysis
    # loop can report the Strand column even when the MAFFT alignment is reused (cached).
    try:
        with open(output_file + ".strandmap.tsv", "w") as sm:
            for rid, st in strand_records.items():
                sm.write(f"{rid}\t{st}\n")
    except OSError:
        pass
    return True, len(written_ids) - 1


def run_mafft_alignment(input_fasta, output_alignment, threads=1):
    mafft_cmd = ["mafft", "--thread", str(threads), "--auto", input_fasta]
    try:
        with open(output_alignment, "w") as outfile_handle:
            process = subprocess.run(mafft_cmd, stdout=outfile_handle, stderr=subprocess.PIPE,
                                     text=True, check=False)
        if process.returncode != 0:
            return False
        if not os.path.exists(output_alignment) or os.path.getsize(output_alignment) == 0:
            return False
        return True
    except Exception:
        return False


def map_alignment_to_original(aligned_sequence_str):
    mapping = []
    original_index = 0
    for char_aligned in aligned_sequence_str:
        if char_aligned == '-':
            mapping.append(None)
        else:
            mapping.append(original_index)
            original_index += 1
    return mapping


def get_aligned_coords(original_start_0based, original_end_0based_exclusive, coord_mapping):
    aligned_start_idx, aligned_end_idx = None, None
    for i, orig_idx in enumerate(coord_mapping):
        if orig_idx is not None and orig_idx >= original_start_0based:
            aligned_start_idx = i
            break
    if aligned_start_idx is None:
        return None, None
    original_last_base_idx = original_end_0based_exclusive - 1
    for i in range(len(coord_mapping) - 1, -1, -1):
        if coord_mapping[i] is not None and coord_mapping[i] <= original_last_base_idx:
            aligned_end_idx = i + 1
            break
    if aligned_end_idx is None or aligned_start_idx >= aligned_end_idx:
        return None, None
    return aligned_start_idx, aligned_end_idx


def translate_and_check(seq_str, feature_name_for_log=""):
    status, protein_seq = "unknown_error", None
    if not seq_str or len(seq_str) < 3:
        return protein_seq, "too_short"
    seq_upper = seq_str.upper()
    if set(seq_upper.replace('N', '')) - set("ATGC"):
        return protein_seq, "invalid_chars"
    remainder = len(seq_upper) % 3
    seq_str_trimmed = seq_upper[:-remainder] if remainder > 0 else seq_upper
    if len(seq_str_trimmed) < 3:
        return protein_seq, "too_short_after_trim"
    try:
        protein_seq = str(Seq(seq_str_trimmed).translate(table=1, cds=False))
        if protein_seq.count('*') > 1 or (protein_seq.count('*') == 1 and not protein_seq.endswith('*')):
            status = "nonsense"
        elif protein_seq.endswith('*'):
            status = "intact"
        else:
            status = "no_stop"
    except CodonTable.TranslationError:
        status = "translation_error"
    return protein_seq, status


def align_and_call_mutations(ref_protein, sample_protein, corrected_protein=None,
                             frameshift_locations=None, stop_at_aa=None, large_insertions=None):
    if not ref_protein:
        return ["ComparisonError"]

    target_protein = corrected_protein if corrected_protein else sample_protein
    if not target_protein:
        return ["Deletion"]

    # Keep only letters: Biopython >=1.78 PairwiseAligner rejects '*' (stop) and other
    # non-alphabet chars ("sequence contains letters not in the alphabet"). rstrip('*')
    # alone misses internal stops; strip all non-letters for cross-version safety.
    r_seq = re.sub(r'[^A-Za-z]', '', ref_protein).upper()
    s_seq = re.sub(r'[^A-Za-z]', '', target_protein).upper()
    if not r_seq or not s_seq:
        return ["Deletion"]

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.1

    try:
        alignments = aligner.align(r_seq, s_seq)
    except ValueError:
        return ["AlignmentFailed"]
    if not alignments:
        return ["AlignmentFailed"]

    best_aln = alignments[0]
    logging.debug(f"Aligning protein - Ref Len: {len(r_seq)}, Sample Len: {len(s_seq)}. Best Aligner Score: {best_aln.score}")
    aln_ref, aln_sample = get_aligned_strings(best_aln)

    mutations = []
    ref_index = 0

    fs_locs = []
    if frameshift_locations:
        if isinstance(frameshift_locations, int):
            fs_locs = [frameshift_locations]
        else:
            fs_locs = sorted(list(set(frameshift_locations)))

    ins_by_codon = {}
    if large_insertions:
        for c, bp in large_insertions:
            ins_by_codon[c] = ins_by_codon.get(c, 0) + bp
    emitted_ins = set()
    
    inserted_tags = set()
    pending_type = None
    pending_start = None
    pending_end = None
    pending_ins_seq = ""
    pending_first_aa = ""
    deferred_tags = []

    def flush_buffer():
        nonlocal pending_type, pending_start, pending_end, pending_ins_seq, pending_first_aa, deferred_tags
        if pending_type == 'del':
            if deferred_tags:
                mutations.extend(deferred_tags)
                deferred_tags = []
            if pending_start == pending_end:
                mutations.append(f"{pending_first_aa}{pending_start}del")
            else:
                mutations.append(f"{pending_start}-{pending_end}del")
        elif pending_type == 'ins':
            if len(pending_ins_seq) > 20:
                logging.debug(f"Summarizing massive insertion of {len(pending_ins_seq)} aa at ref index {pending_start}")
                mutations.append(f"{pending_start}ins{len(pending_ins_seq)}aa")
            else:
                mutations.append(f"{pending_start}ins{pending_ins_seq}")
            
            if deferred_tags:
                mutations.extend(deferred_tags)
                deferred_tags = []
        pending_type = None
        pending_start = None
        pending_end = None
        pending_ins_seq = ""
        pending_first_aa = ""

    for i in range(len(aln_ref)):
        aa_r = aln_ref[i]
        aa_s = aln_sample[i]

        if aa_r != '-':
            ref_index += 1
            
        if ins_by_codon and aa_r != '-' and ref_index in ins_by_codon and ref_index not in emitted_ins:
            flush_buffer()
            mutations.append(f"{ref_index}ins~{ins_by_codon[ref_index]}bp")
            emitted_ins.add(ref_index)
            
        current_tags = []
        if fs_locs:
            for fs_loc in fs_locs:
                if ref_index == fs_loc and fs_loc not in inserted_tags:
                    logging.debug(f"Protein Aligner matched Frameshift Tag {fs_loc} at ref_index {ref_index}")
                    tag_index = fs_loc
                    if pending_type == 'del':
                        tag_index = pending_start
                    tag_str = f"Frameshift_at_{tag_index}"
                    if stop_at_aa is not None and stop_at_aa > 0 and fs_loc <= stop_at_aa + 2 and "Premature_Stop" not in str(inserted_tags):
                        tag_str += "-Premature_Stop"
                    if tag_str not in inserted_tags:
                        current_tags.append(tag_str)
                        inserted_tags.add(fs_loc)
                        inserted_tags.add(tag_str)

        if current_tags:
            if pending_type == 'del':
                deferred_tags.extend(current_tags)
            else:
                flush_buffer()
                mutations.extend(current_tags)

        if aa_r == aa_s:
            flush_buffer()
            continue

        if aa_r == '-':
            if 'X' in aa_s:
                flush_buffer()
                continue
            if pending_type == 'ins' and pending_start == ref_index:
                pending_ins_seq += aa_s
            else:
                flush_buffer()
                pending_type = 'ins'
                pending_start = ref_index
                pending_ins_seq = aa_s
        elif aa_s == '-' or aa_s == 'X':
            if pending_type == 'del' and ref_index == pending_end + 1:
                pending_end = ref_index
            else:
                flush_buffer()
                pending_type = 'del'
                pending_start = ref_index
                pending_end = ref_index
                pending_first_aa = aa_r
        else:
            flush_buffer()
            mutations.append(f"{aa_r}{ref_index}{aa_s}")
    flush_buffer()

    if len(target_protein) < len(ref_protein) * 0.5 and not inserted_tags:
        mutations.append("Truncation")
        
    for c in sorted(ins_by_codon):
        if c not in emitted_ins:
            mutations.append(f"{c}ins~{ins_by_codon[c]}bp")
            
    return mutations if mutations else ["None"]


def calculate_alignment_identity(seq1, seq2):
    if not seq1 or not seq2:
        return 0.0
    # Biopython >=1.78 PairwiseAligner.score validates against its alphabet and raises
    # "sequence contains letters not in the alphabet" on '*' (stop), '-', etc. Strip to
    # plain letters so identity scoring works regardless of Biopython version.
    s1 = re.sub(r'[^A-Za-z]', '', str(seq1)).upper()
    s2 = re.sub(r'[^A-Za-z]', '', str(seq2)).upper()
    if not s1 or not s2:
        return 0.0
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    try:
        score = aligner.score(s1, s2)
    except ValueError:
        # last-resort positional identity if the aligner still rejects a character
        n = min(len(s1), len(s2))
        return (sum(1 for a, b in zip(s1, s2) if a == b) / n) if n > 0 else 0.0
    return score / min(len(s1), len(s2)) if min(len(s1), len(s2)) > 0 else 0.0


def calculate_simple_coverage(query_seq_ungapped, ref_seq_ungapped):
    if not query_seq_ungapped or not ref_seq_ungapped:
        return 0.0
    query_len = len(query_seq_ungapped.replace('N', ''))
    ref_len = len(ref_seq_ungapped.replace('N', ''))
    return query_len / ref_len if ref_len > 0 else 0.0


def calculate_reference_coverage(gapped_ref, gapped_sample):
    """Fraction of KCON reference feature positions covered by the sample.

    Walks the aligned columns of a feature window (ref and sample slices are the
    same length, taken from the same MSA columns). Counts only columns where the
    REFERENCE has a real base (a true KCON position); of those, the fraction where
    the SAMPLE also has a real base. Insertions relative to KCON (ref-gap columns)
    are ignored entirely, so the result is always in [0, 1] and answers "how much
    of the KCON feature is present in this sample" rather than the old
    sample_len/ref_len ratio, which could exceed 1.0 when the sample carried an
    insertion. N's count as absent on both sides."""
    ref_positions = 0
    covered = 0
    for rc, sc in zip(gapped_ref, gapped_sample):
        if rc != '-' and rc != 'N':
            ref_positions += 1
            if sc != '-' and sc != 'N':
                covered += 1
    return covered / ref_positions if ref_positions > 0 else 0.0


def trim_alignment_to_kcon(alignment_records, kcon_ref_id, output_file):
    cleaned = {}
    for rid, rec in alignment_records.items():
        cleaned[rid] = SeqRecord(Seq(str(rec.seq).upper()), id=rec.id, description=rec.description)

    kcon_rec = cleaned.get(kcon_ref_id)
    if kcon_rec is None:
        return cleaned

    kcon_str = str(kcon_rec.seq)
    leading = len(kcon_str) - len(kcon_str.lstrip('-'))
    trailing = len(kcon_str) - len(kcon_str.rstrip('-'))

    if leading > 0 or trailing > 0:
        trim_end = len(kcon_str) - trailing
        logging.info(f"Trimming alignment to KCON bounds: removed {leading} leading and {trailing} trailing gap columns.")
        for rid, rec in cleaned.items():
            cleaned[rid] = SeqRecord(Seq(str(rec.seq)[leading:trim_end]),
                                     id=rec.id, description=rec.description)

    with open(output_file, 'w') as f:
        SeqIO.write(list(cleaned.values()), f, "fasta")
    return cleaned


def parse_args():
    parser = argparse.ArgumentParser(description="Align haplotypes to KCON, analyze ORF integrity.")
    parser.add_argument('--provirus_name', type=str, required=True)
    parser.add_argument('--input_haplotype_dir', type=str, required=True)
    parser.add_argument('--output_analysis_dir', type=str, required=True)
    parser.add_argument('--kcon_fasta', type=str, required=True)
    parser.add_argument('--kcon_ref_id', type=str, required=True)
    parser.add_argument('--locus_type', choices=['TypeI', 'TypeII'], required=True)
    parser.add_argument('--force_mafft', action='store_true')
    parser.add_argument('--debug', action='store_true', help="Print live indel block math to terminal.")
    parser.add_argument('--mafft_threads', type=int, default=1)
    parser.add_argument('--kcon_fasta_solo_ltr', type=str, required=True, help="Path to Solo LTR KCON reference FASTA.")
    parser.add_argument(
        '--tsd-evidence-only', action='store_true',
        help="Resolve retained-sequence TSD junctions and write requested JSONL evidence, then stop before MAFFT/ORF analysis.",
    )
    return parser.parse_args()

def trim_local_gaps(gapped_ref, gapped_sample):
    """
    Trims leading and trailing gap columns from a localized MAFFT slice.
    This ensures the frame-correction math starts cleanly on base 0.
    """
    ref_list = list(gapped_ref)
    samp_list = list(gapped_sample)
    
    # Trim leading gaps where the reference is missing
    while ref_list and ref_list[0] == '-':
        ref_list.pop(0)
        samp_list.pop(0)
        
    # Trim trailing gaps where the reference is missing
    while ref_list and ref_list[-1] == '-':
        ref_list.pop()
        samp_list.pop()
        
    return "".join(ref_list), "".join(samp_list)

if __name__ == "__main__":
    args = parse_args()
    args.provirus_name = args.provirus_name.rstrip('/')
    args.output_analysis_dir = args.output_analysis_dir.rstrip('/')

    os.makedirs(args.output_analysis_dir, exist_ok=True)

    # --- Type-I/Type-II reference-routing guard --------------------------------------------------
    # The controller resolves BOTH the KCON reference FASTA (--kcon_fasta / --kcon_ref_id) and the
    # feature-coordinate selector (--locus_type) from the locus's type field. If those two ever
    # disagree -- a routing bug where a Type-I locus is handed the Type-II reference (or vice versa)
    # -- the wrong KCON_FEATURES_TYPE* coordinate set is applied to the reference and every ORF call
    # at the locus is silently wrong (pol boundary off by the Delta292 span, wrong accessory ORF).
    # The reference id embeds its type ("type1_KCON" / "type2_KCON"); assert it matches --locus_type
    # and abort loudly rather than emit a mis-typed table. Solo-LTR / non-typed ids are exempted.
    _rid = str(args.kcon_ref_id).lower()
    _expected_type_token = "type1" if args.locus_type == "TypeI" else "type2"
    _wrong_type_token = "type2" if args.locus_type == "TypeI" else "type1"
    if _wrong_type_token in _rid and _expected_type_token not in _rid:
        print(f"FATAL: KCON reference routing mismatch -- --locus_type={args.locus_type} but "
              f"--kcon_ref_id='{args.kcon_ref_id}' names the other type. The controller passed a "
              f"reference that does not match the locus type; refusing to emit a mis-typed ORF "
              f"table.", file=sys.stderr)
        sys.exit(2)

    file_prefix = f"{args.provirus_name}_"
    ORIENT_SUFFIX = f"_{args.kcon_ref_id}"
    LOG_FILE_PATH = os.path.join(args.output_analysis_dir, f"{file_prefix}orf_analysis{ORIENT_SUFFIX}.log")
    setup_logging(LOG_FILE_PATH, is_debug_mode=args.debug)

    LTR_DETAIL_LOG_PATH = os.path.join(args.output_analysis_dir, f"{file_prefix}ltr_alignment_details{ORIENT_SUFFIX}.log")
    detail_logger = setup_detailed_logging(LTR_DETAIL_LOG_PATH)

    logging.info(f"--- Starting ORF Analysis Script ---")

    TEMP_DIR_FOR_PAF = os.path.join(args.output_analysis_dir, "temp_paf_files_for_orientation")
    os.makedirs(TEMP_DIR_FOR_PAF, exist_ok=True)

    COMBINED_FOR_MAFFT_FASTA_FILE = os.path.join(args.output_analysis_dir, f"{file_prefix}all_haplotypes_for_mafft{ORIENT_SUFFIX}.fasta")
    ALIGNMENT_FILE = os.path.join(args.output_analysis_dir, f"{file_prefix}all_haplotypes_aligned{ORIENT_SUFFIX}.fasta")
    RESULTS_CSV_FILE = os.path.join(args.output_analysis_dir, f"{file_prefix}orf_integrity_results{ORIENT_SUFFIX}.csv")

    KCON_FEATURES_TYPE1 = {"gag": (1111, 3112), "pro": (2913, 3918), "pol": (3878, 6501),
                           "env": (6512, 8258), "np9_1": (6450, 6494), "np9_2": (8118, 8299)}
    KCON_FEATURES_TYPE2 = {"gag": (1111, 3112), "pro": (2913, 3918), "pol": (3878, 6749),
                           "env": (6450, 8550), "rec_1": (6450, 6711), "rec_2": (8410, 8466)}
    KCON_LTRS_TYPE1 = {"5_LTR": (0, 967), "3_LTR": (8213, 9180)}
    KCON_LTRS_TYPE2 = {"5_LTR": (0, 967), "3_LTR": (8504, 9471)}

    ACCESSORY_ORF_NAME = "np9" if args.locus_type == 'TypeI' else "rec"
    ACCESSORY_EXON1_KEY = "np9_1" if args.locus_type == 'TypeI' else "rec_1"
    ACCESSORY_EXON2_KEY = "np9_2" if args.locus_type == 'TypeI' else "rec_2"

    pristine_kcon_record = SeqIO.read(args.kcon_fasta, "fasta")
    pristine_kcon_record.seq = Seq(str(pristine_kcon_record.seq).upper())

    KCON_FEATURES = KCON_FEATURES_TYPE1 if args.locus_type == 'TypeI' else KCON_FEATURES_TYPE2
    KCON_LTRS = KCON_LTRS_TYPE1 if args.locus_type == 'TypeI' else KCON_LTRS_TYPE2

    pristine_kcon_feature_proteins = {}
    pristine_kcon_feature_ungapped_dna = {}
    pristine_kcon_orf_status = {}
    
    for feature in ["gag", "pro", "pol", "env"]:
        start, end = KCON_FEATURES[feature]
        dna_seq = str(pristine_kcon_record.seq[start:end])
        prot, status = translate_and_check(dna_seq, f"Pristine KCON {feature}")
        pristine_kcon_feature_ungapped_dna[feature] = dna_seq
        pristine_kcon_feature_proteins[feature] = prot
        pristine_kcon_orf_status[feature] = status

    ex1_s, ex1_e = KCON_FEATURES[ACCESSORY_EXON1_KEY]
    ex2_s, ex2_e = KCON_FEATURES[ACCESSORY_EXON2_KEY]
    spliced_dna = str(pristine_kcon_record.seq[ex1_s:ex1_e]) + str(pristine_kcon_record.seq[ex2_s:ex2_e])
    prot, status = translate_and_check(spliced_dna, f"Pristine KCON {ACCESSORY_ORF_NAME}")
    pristine_kcon_feature_ungapped_dna[ACCESSORY_ORF_NAME] = spliced_dna
    pristine_kcon_feature_proteins[ACCESSORY_ORF_NAME] = prot
    pristine_kcon_orf_status[ACCESSORY_ORF_NAME] = status

    all_haplotype_files = get_fasta_files_from_dir(args.input_haplotype_dir)
    placeholder_data, solo_ltr_data, multi_provirus_tsds, single_provirus_tsds, structural_solo_ltrs = perform_pre_analysis(
        all_haplotype_files, args.kcon_fasta_solo_ltr, TEMP_DIR_FOR_PAF, detail_logger, args.mafft_threads,
        kcon_ref_full=args.kcon_fasta
    )
    if args.tsd_evidence_only:
        logging.info(
            "TSD evidence-only resolution complete: files=%d placeholders=%d solo=%d "
            "multi_calls=%d single_calls=%d structural_solo=%d evidence=%s",
            len(all_haplotype_files), len(placeholder_data), len(solo_ltr_data),
            len(multi_provirus_tsds), len(single_provirus_tsds),
            len(structural_solo_ltrs),
            os.environ.get("HML2_TSD_JUNCTION_EVIDENCE_JSONL", ""),
        )
        sys.exit(0)

    if not (os.path.exists(ALIGNMENT_FILE) and os.path.getsize(ALIGNMENT_FILE) > 0 and not args.force_mafft):
        success, num_passed = combine_fasta_files_for_mafft(
            all_haplotype_files, args.kcon_fasta, args.kcon_ref_id, COMBINED_FOR_MAFFT_FASTA_FILE,
            args.kcon_fasta_solo_ltr, args.mafft_threads, TEMP_DIR_FOR_PAF)
        if not success:
            sys.exit(1)
        num_seqs_for_mafft = sum(1 for _ in SeqIO.parse(COMBINED_FOR_MAFFT_FASTA_FILE, "fasta"))
        if num_seqs_for_mafft <= 1:
            subprocess.run(['cp', COMBINED_FOR_MAFFT_FASTA_FILE, ALIGNMENT_FILE], check=True)
        else:
            if not run_mafft_alignment(COMBINED_FOR_MAFFT_FASTA_FILE, ALIGNMENT_FILE, threads=args.mafft_threads):
                sys.exit(1)

    try:
        raw_alignment = list(AlignIO.read(ALIGNMENT_FILE, "fasta"))
    except ValueError:
        sys.exit(1)

    alignment_records = {r.id: r for r in raw_alignment}
    alignment_records = trim_alignment_to_kcon(alignment_records, args.kcon_ref_id, ALIGNMENT_FILE)

    kcon_aligned_record = alignment_records.get(args.kcon_ref_id)
    if not kcon_aligned_record:
        sys.exit(1)

    kcon_coord_mapping = map_alignment_to_original(str(kcon_aligned_record.seq))
    aligned_kcon_feature_coords = {name: get_aligned_coords(s, e, kcon_coord_mapping) for name, (s, e) in KCON_FEATURES.items()}
    aligned_kcon_ltr_coords = {name: get_aligned_coords(s, e, kcon_coord_mapping) for name, (s, e) in KCON_LTRS.items()}

    # Per-copy genomic strand (provirus vs reference genome) written by
    # combine_fasta_files_for_mafft, keyed by combined-FASTA record id. '+'/'-' vs
    # hg38/CHM13, or 'NA' when the extraction didn't emit a genomic_strand field.
    strand_map = {}
    try:
        with open(COMBINED_FOR_MAFFT_FASTA_FILE + ".strandmap.tsv") as sm:
            for line in sm:
                rid, _, st = line.rstrip("\n").partition("\t")
                if rid:
                    strand_map[rid] = st
    except OSError:
        pass

    results_data = []

    for seq_id, aligned_record in alignment_records.items():
        if seq_id == args.kcon_ref_id:
            continue

        id_full_raw = aligned_record.id.lstrip('>')
        id_full = id_full_raw.split('|')[0]
        source_identifier = extract_source_id(aligned_record.id, aligned_record.description)

        parts = id_full.split('_')
        sample_id = parts[0]
        haplotype = next((p.replace('hap', 'h') for p in parts if p in ['mat', 'pat'] or re.match(r'^h(ap)?\d+$', p)), 'h?')

        sample_result = {"ID_Full": id_full, "ID": sample_id, "Haplotype": haplotype, "Source_Identifier": source_identifier,
                         "Strand": strand_map.get(id_full_raw, "NA")}

        # An _asmdup record is the single representative of an assembly-gap duplicate:
        # analyze it as an ordinary single provirus (coverage-based structure), with
        # the _asmdup tag in its ID flagging it for read-depth follow-up.
        sample_aligned_seq_str = str(aligned_record.seq).upper()
        # TSD = the whole-array / per-WINDOW target-site duplication (never per copy):
        #   1) own id  -> ordinary single provirus, a separate-contig CNV/paralog, or a
        #      _MULTI window the splitter resolved to one provirus (tag already stripped).
        #   2) `_alt#` is exact-copy only and never inherits from the window.
        #   3) other splitter suffixes may reduce to the window base id; tandem
        #      `_part#` units share that array-level integration boundary.
        tsd_5, tsd_3 = "NA", "NA"
        if id_full in single_provirus_tsds:
            tsd_5, tsd_3 = single_provirus_tsds[id_full]
        elif id_full in multi_provirus_tsds:
            tsd_5, tsd_3 = multi_provirus_tsds[id_full]
        elif re.search(r'_alt\d+', id_full):
            # ``_alt#`` is a separately resolved CNV/paralogous insertion copy.
            # Its TSD belongs to its own two LTR-host junctions.  If the exact
            # per-copy extraction did not recover those junctions, keep NA and
            # allow the alignment fallback/no-signal path below; never inherit
            # the outer multi-copy window's TSD.  Only tandem ``_part#`` units
            # share one array-level integration boundary.
            pass
        else:
            base_id_lookup = re.sub(r'(?:_alt\d+|_(?:MULTI|DOUBLE)|_part\d+|_asmdup)+$', '', id_full)
            if base_id_lookup != id_full:
                cand = multi_provirus_tsds.get(base_id_lookup)
                if cand and cand[0] not in ("NONE", "NA"):
                    tsd_5, tsd_3 = cand
                elif base_id_lookup in single_provirus_tsds:
                    tsd_5, tsd_3 = single_provirus_tsds[base_id_lookup]
        if tsd_5 == "NA":
            five_ltr_start, five_ltr_end = aligned_kcon_ltr_coords.get("5_LTR", (None, None))
            three_ltr_start, three_ltr_end = aligned_kcon_ltr_coords.get("3_LTR", (None, None))

            if all(c is not None for c in [five_ltr_start, three_ltr_end]):
                def align_idx_to_seq_idx(align_str, target_align_idx):
                    return len(align_str[:target_align_idx].replace("-", ""))

                try:
                    raw_start_idx = align_idx_to_seq_idx(sample_aligned_seq_str, five_ltr_start)
                    raw_end_idx = align_idx_to_seq_idx(sample_aligned_seq_str, three_ltr_end)
                    ungapped_seq = sample_aligned_seq_str.replace("-", "")
                    # The MAFFT projection is only an anchor.  Apply the same
                    # independently LTR-supported +/-2 bp terminal refinement
                    # used by the primary pre-analysis route before reading host
                    # flanks.  The TSD strings themselves never select a shift.
                    five_ltr_ref = str(pristine_kcon_record.seq[
                        KCON_LTRS["5_LTR"][0]:KCON_LTRS["5_LTR"][1]
                    ])
                    three_ltr_ref = str(pristine_kcon_record.seq[
                        KCON_LTRS["3_LTR"][0]:KCON_LTRS["3_LTR"][1]
                    ])
                    refined_start_idx = _refine_projected_ltr_junction(
                        ungapped_seq, five_ltr_ref, raw_start_idx, "left", "+"
                    )
                    refined_end_idx = _refine_projected_ltr_junction(
                        ungapped_seq, three_ltr_ref, raw_end_idx, "right", "+"
                    )
                    refined_tsd_5, refined_tsd_3, _, _ = refine_tsd_search(
                        ungapped_seq, refined_start_idx, refined_end_idx,
                        search_window=15,
                    )
                    tsd_5 = refined_tsd_5
                    tsd_3 = refined_tsd_3
                except Exception:
                    # Failure to observe a host flank contributes no TSD signal;
                    # it is not evidence that the two insertion-time copies
                    # differed and must not abort the otherwise valid ORF row.
                    tsd_5, tsd_3 = "NONE", "NONE"
            else:
                tsd_5, tsd_3 = "NONE", "NONE"

        # GRCh38 is the observed integrated reference allele for 8q24.3c.  Its
        # two genomic-forward host duplications were read directly at
        # chr8:145021238-145021243 (ATAATC) and
        # chr8:145028836-145028841 (ACAATC).  The provirus is reverse-oriented,
        # so the insertion-oriented pair emitted by this table is the
        # reverse-complemented, side-swapped GATTGT/GATTAT.  This exact authority
        # applies only to the genuine GRCh38 reference record; every biological
        # sample remains algorithmic and may retain mutations or missing sides.
        is_grch38_8q24_3c_reference = (
            args.provirus_name == "HML-2_8q24.3c"
            and (
                "GRCh38" in source_identifier
                or "GRCh38" in id_full
                or "GCA_000001405.15" in id_full
            )
        )
        if is_grch38_8q24_3c_reference:
            tsd_5, tsd_3 = "GATTGT", "GATTAT"
            sample_result["tsd_source_interval"] = (
                "chr8:145021238-145021243;chr8:145028836-145028841"
            )

        # "Provirus_from_Multi" is reserved for TANDEM-ARRAY units: a _part# suffix is a
        # shared-LTR repeat of ONE insertion. A copy with only _alt# (a segmental
        # duplication / dispersed CNV copy -- a SEPARATE insertion, NO _part) is an
        # INDEPENDENT provirus, scored by its own coverage like any single provirus; its
        # copy-number status lives in the _alt# of the ID (for the CNV analysis), not here.
        # _part# is emitted ONLY by the splitter for tandem-array units (extraction
        # never emits it), so it alone reliably marks an array member.
        is_tandem_member = bool(re.search(r'_part\d+', id_full))
        if is_tandem_member:
            structure = "Provirus_from_Multi"
        elif "SOLO_LTR" in source_identifier or "SOLO_LTR" in id_full_raw:
            structure = "Solo-LTR"
        elif id_full in structural_solo_ltrs:
            structure = "Solo-LTR"
        else:
            # Fraction of KCON reference positions actually covered by a non-gap
            # sample base. Using sample_len/ref_len here inflated oversized
            # DEGRADED windows (flanks + insertions) past 0.80 and mislabeled
            # them Provirus despite ~0 real feature coverage.
            coverage = calculate_reference_coverage(str(kcon_aligned_record.seq), sample_aligned_seq_str)
            structure = "Provirus" if coverage >= 0.80 else "Fragment"

        sample_result["Structure"] = structure
        sample_result["5'_TSD"] = tsd_5
        sample_result["3'_TSD"] = tsd_3
        sample_result["TSD_Orientation"] = "INSERTION"

        # --- Observation-vs-latent-copy contract -----------------------------------------------
        # Split/extracted rows are OBSERVATIONS, never automatic biological copies. Suffixes are
        # measurement provenance only: _alt does not prove an independent CNV, _asmdup does not prove
        # artifact, and no observed sequence is selected as founder truth. A downstream, exact-ID
        # latent-copy ledger supplies normalized {single,artifact,true_dup} weights and performs the
        # duplication-birth coupling / later-divergence marginalization. Until then all weights remain
        # explicit NA and publication must fail closed. The observed sequence bytes and digest remain
        # in this table even when the eventual artifact state has nonzero weight.
        _alt_m = re.search(r'_alt(\d+)', id_full)
        _part_m = re.search(r'_part(\d+)', id_full)
        _is_asmdup = bool(re.search(r'_asmdup', id_full))
        _alt_idx = int(_alt_m.group(1)) if _alt_m else (1 if not _part_m else 0)
        _part_idx = int(_part_m.group(1)) if _part_m else 0
        if _part_m:
            observation_class = "tandem_part_observation"
        elif _is_asmdup:
            observation_class = "assembly_gap_duplicate_observation"
        elif _alt_m:
            observation_class = "alternate_candidate_observation"
        else:
            observation_class = "single_candidate_observation"
        observation_bytes = sample_aligned_seq_str.replace("-", "").encode("ascii", errors="replace")
        sample_result["copy_class"] = "latent_unresolved"
        sample_result["independent_copy"] = "UNKNOWN"
        sample_result["copy_observation_class"] = observation_class
        sample_result["observation_role"] = "measurement_evidence"
        sample_result["observation_sequence_bytes"] = len(observation_bytes)
        sample_result["observation_sequence_sha256"] = hashlib.sha256(observation_bytes).hexdigest()
        sample_result["artifact_evidence_retained"] = "YES"
        sample_result["copy_state_model_status"] = "UNMEASURED"
        sample_result["copy_state_single_weight"] = "NA"
        sample_result["copy_state_artifact_weight"] = "NA"
        sample_result["copy_state_true_dup_weight"] = "NA"
        sample_result["copy_state_weight_provenance"] = "REQUIRES_EXACT_LATENT_COPY_LEDGER"
        sample_result["duplication_birth_constraint"] = "UNMEASURED"
        sample_result["divergence_marginalization_model"] = "UNMEASURED"
        sample_result["alt_index"] = _alt_idx   # which separate insertion (>=1); 0 if purely a tandem unit
        sample_result["part_index"] = _part_idx  # which tandem-array unit within the insertion (0 if none)

        for feature in ["gag", "pro", "pol", "env"]:
            aln_s, aln_e = aligned_kcon_feature_coords.get(feature, (None, None))
            if aln_s is None or aln_e is None:
                status, prot, missense, dna_id, prot_id, cov = "ref_map_error", None, ["ref_map_error"], "NA", "NA", "NA"
            else:
                gapped_dna_chunk = sample_aligned_seq_str[aln_s:aln_e]
                gapped_ref_chunk = str(kcon_aligned_record.seq)[aln_s:aln_e]
                
                clean_ref_chunk, clean_samp_chunk = trim_local_gaps(gapped_ref_chunk, gapped_dna_chunk)
                
                ref_prot = pristine_kcon_feature_proteins[feature]
                ref_dna = pristine_kcon_feature_ungapped_dna[feature]

                status, prot, corr_prot, fs_locs, stop_idx, large_ins = analyze_orf_structural_integrity(
                    clean_samp_chunk, ref_prot, ref_dna, gapped_ref_dna=clean_ref_chunk,
                    feature_name=feature, seq_id=sample_id)
                ungapped_dna = gapped_dna_chunk.replace("-", "")

                cov_val = calculate_reference_coverage(gapped_ref_chunk, gapped_dna_chunk)
                cov = f"{cov_val:.4f}"
                if cov_val < COVERAGE_FLOOR_FOR_IDENTITY:
                    # Too little of the feature is present for identity to be meaningful.
                    dna_id = "NA"
                    prot_id = "NA"
                else:
                    dna_id = f"{calculate_alignment_identity(ungapped_dna, ref_dna):.4f}"
                    prot_id = f"{calculate_alignment_identity(prot, ref_prot):.4f}" if prot else "0.0000"

                missense = align_and_call_mutations(ref_prot, prot,
                                                    corrected_protein=corr_prot,
                                                    frameshift_locations=fs_locs,
                                                    stop_at_aa=stop_idx,
                                                    large_insertions=large_ins)
            sample_result[feature] = status
            sample_result[f"missense_{feature}"] = ",".join(missense)
            sample_result[f"{feature}_dna_identity"] = dna_id
            sample_result[f"{feature}_protein_identity"] = prot_id
            sample_result[f"{feature}_coverage"] = cov

        ex1_as, ex1_ae = aligned_kcon_feature_coords.get(ACCESSORY_EXON1_KEY, (None, None))
        ex2_as, ex2_ae = aligned_kcon_feature_coords.get(ACCESSORY_EXON2_KEY, (None, None))

        if ex1_as is None or ex2_as is None or ex1_ae is None or ex2_ae is None:
            status, missense, dna_id, prot_id, cov = "ref_map_error", ["ref_map_error"], "NA", "NA", "NA"
        else:
            chunk_ex1 = sample_aligned_seq_str[ex1_as:ex1_ae]
            chunk_ex2 = sample_aligned_seq_str[ex2_as:ex2_ae]
            gapped_spliced = chunk_ex1 + chunk_ex2
            
            ref_gapped_ex1 = str(kcon_aligned_record.seq)[ex1_as:ex1_ae] 
            ref_gapped_ex2 = str(kcon_aligned_record.seq)[ex2_as:ex2_ae] 
            ref_gapped_spliced = ref_gapped_ex1 + ref_gapped_ex2

            clean_ref_spliced, clean_samp_spliced = trim_local_gaps(ref_gapped_spliced, gapped_spliced)

            ungapped_spliced = gapped_spliced.replace("-", "")

            ref_prot = pristine_kcon_feature_proteins[ACCESSORY_ORF_NAME]
            ref_dna = pristine_kcon_feature_ungapped_dna[ACCESSORY_ORF_NAME]

            status, prot, corr_prot, fs_locs, stop_idx, large_ins = analyze_orf_structural_integrity(
                clean_samp_spliced, ref_prot, ref_dna,
                feature_name=ACCESSORY_ORF_NAME, locus_type=args.locus_type,
                gapped_ref_dna=clean_ref_spliced, seq_id=sample_id) 

            cov_val = calculate_reference_coverage(ref_gapped_spliced, gapped_spliced)
            cov = f"{cov_val:.4f}"
            if cov_val < COVERAGE_FLOOR_FOR_IDENTITY:
                # Too little of the feature is present for identity to be meaningful.
                dna_id = "NA"
                prot_id = "NA"
            else:
                dna_id = f"{calculate_alignment_identity(ungapped_spliced, ref_dna):.4f}"
                prot_id = f"{calculate_alignment_identity(prot, ref_prot):.4f}" if prot else "0.0000"

            missense = align_and_call_mutations(ref_prot, prot,
                                                corrected_protein=corr_prot,
                                                frameshift_locations=fs_locs,
                                                stop_at_aa=stop_idx,
                                                large_insertions=large_ins)

        sample_result[ACCESSORY_ORF_NAME] = status
        sample_result[f"missense_{ACCESSORY_ORF_NAME}"] = ",".join(missense)
        sample_result[f"{ACCESSORY_ORF_NAME}_dna_identity"] = dna_id
        sample_result[f"{ACCESSORY_ORF_NAME}_protein_identity"] = prot_id
        sample_result[f"{ACCESSORY_ORF_NAME}_coverage"] = cov

        # --- Physical Pol-region presence  vs  translated-Pol compatibility --------------------------
        # These are TWO different biological facts and are reported as SEPARATE fields so no downstream
        # analysis can collapse them into one "pol" verdict:
        #
        #   pol_region_present  -- PHYSICAL fact: is enough of the Pol coding region physically present
        #       in this copy for the RT / RNase-H / integrase sequence to exist at all? Coverage-based,
        #       translation-agnostic. Preserved Pol *sequence* is sequence POTENTIAL only; on its own it
        #       is NEVER a translated-Pol burden. (Same threshold used for the aggregate structure call.)
        #
        #   pol_translation_compatible -- FUNCTIONAL fact: could a Pol product actually be translated in
        #       this copy, in a TYPE-SPECIFIC way?
        #         * Type II: native Gag-Pro-Pol read-through -- positive only if the physical `pol` ORF
        #           verdict is translation-competent (Intact / Intact_FS_End / Fragment_Intact / no_stop).
        #         * Type I : a Type-I allele CANNOT be complete native Pol -- the exact 292-nt deletion at
        #           KCON [6501,6793) removes the Pol C-terminus and shifts frame. Its ONLY positive
        #           translation state is an ALTERED read-through past the Delta292 junction into
        #           Env-derived sequence. So a Type-I copy's Pol translation is reported here as the
        #           altered read-through state, never as native Pol. We approximate the read-through as
        #           translation-competent iff BOTH the (truncated) Pol region and Env are themselves
        #           translation-competent in this copy (the frame runs from Pol through the junction into
        #           Env without a premature stop); otherwise it is not compatible.
        # The physical `pol` column and its identities/coverage are UNCHANGED; these are additive.
        _pol_status = sample_result.get("pol", "err")
        _env_status = sample_result.get("env", "err")
        try:
            _pol_cov_val = float(sample_result.get("pol_coverage", "0") or 0)
        except (TypeError, ValueError):
            _pol_cov_val = 0.0

        pol_region_present = "Y" if _pol_cov_val >= 0.80 else "N"

        if args.locus_type == "TypeI":
            # No native Pol is possible for a Type-I allele; positive state = altered Env read-through.
            _readthrough_ok = (_pol_status in TRANSLATION_COMPETENT_STATUSES
                               and _env_status in TRANSLATION_COMPETENT_STATUSES)
            pol_translation_compatible = "altered_readthrough" if _readthrough_ok else "N"
            pol_translation_basis = "type1_delta292_env_readthrough"
        else:
            _native_ok = _pol_status in TRANSLATION_COMPETENT_STATUSES
            pol_translation_compatible = "native" if _native_ok else "N"
            pol_translation_basis = "type2_native_gag_pro_pol"

        sample_result["pol_region_present"] = pol_region_present
        sample_result["pol_translation_compatible"] = pol_translation_compatible
        sample_result["pol_translation_basis"] = pol_translation_basis

        # Aggregate ORF-intactness. Pol contributes via TRANSLATION compatibility (not mere physical
        # presence), so a Type-I Delta292 allele is never counted as an intact native-Pol provirus.
        all_core_ok = all(sample_result.get(f, "err") in TRANSLATION_COMPETENT_STATUSES
                          for f in ["gag", "pro", "env", ACCESSORY_ORF_NAME])
        all_orfs_functional = all_core_ok and (pol_translation_compatible != "N")
        sample_result["all_orfs_intact"] = "Y" if all_orfs_functional else "N"
        results_data.append(sample_result)

    # An Absent/placeholder row means "this sample+haplotype has no sequence here". That is FALSE
    # if the same sample produced a real call (Fragment / Provirus / Solo-LTR / Multi ...).
    # perform_pre_analysis emits one placeholder per _ABSENT_OR_UNALIGNED input file, so when the
    # ORF input is double-staged (a real extraction AND an absent placeholder for one sample both
    # reach ORF) this would otherwise emit BOTH rows -- the duplication seen on the _hg38 loci.
    # Fix: keep solo-LTR (real); drop a placeholder when a real call exists for that (ID, Haplotype);
    # and collapse multiple placeholders for one sample to a single row.
    #
    # Collapsing is class-aware: when one sample+haplotype produced BOTH an absence placeholder and
    # a technical-failure placeholder, the technical failure WINS. Dropping it would leave the row
    # asserting "successfully interrogated, no insertion" about a site that was never interrogated.
    real_keys = {(r.get("ID"), r.get("Haplotype")) for r in results_data}
    for s_data in solo_ltr_data:                       # solo-LTR calls are real -> always kept
        real_keys.add((s_data.get("ID"), s_data.get("Haplotype")))
        results_data.append(s_data)
    kept_placeholders = {}
    _dropped_ph = 0
    for p_data in placeholder_data:
        key = (p_data.get("ID"), p_data.get("Haplotype"))
        if key in real_keys:
            _dropped_ph += 1
            continue
        held = kept_placeholders.get(key)
        if held is None:
            kept_placeholders[key] = p_data
            continue
        _dropped_ph += 1
        if held.get("Structure") != "Technical_Error" and p_data.get("Structure") == "Technical_Error":
            kept_placeholders[key] = p_data
    results_data.extend(kept_placeholders.values())
    if _dropped_ph:
        logging.info(f"Dropped {_dropped_ph} redundant placeholder row(s) (real call or duplicate placeholder exists for that sample+haplotype).")

    if results_data:
        df = pd.DataFrame(results_data)

        # Provirus type is a property of the locus (set by the KCON reference used);
        # emit it explicitly so downstream (R heatmaps) doesn't have to infer it from
        # which accessory ORF column (np9=type1 / rec=type2) happens to be populated.
        df["provirus_type"] = "type1" if args.locus_type == "TypeI" else "type2"

        # NOTE: an all-absent locus is NOT written off here. Every haplotype being a confident
        # non-carrier is a real, publishable result for this locus, and exiting without a results
        # file would make the locus indistinguishable from one whose ORF stage never ran -- the
        # combiner would then have to treat the whole locus as a technical failure. Always emit.

        core_orfs = ["gag", "pro", "pol", "env"]
        all_orfs = core_orfs + [ACCESSORY_ORF_NAME]

        ordered_cols = ["ID_Full", "ID", "Haplotype", "Source_Identifier", "Strand", "Structure", "provirus_type",
                        "copy_class", "independent_copy", "copy_observation_class", "observation_role",
                        "observation_sequence_bytes", "observation_sequence_sha256", "artifact_evidence_retained",
                        "copy_state_model_status", "copy_state_single_weight", "copy_state_artifact_weight",
                        "copy_state_true_dup_weight", "copy_state_weight_provenance",
                        "duplication_birth_constraint", "divergence_marginalization_model",
                        "alt_index", "part_index", "5'_TSD", "3'_TSD", "TSD_Orientation"]
        for orf in all_orfs:
            ordered_cols.extend([
                orf, f"missense_{orf}", f"{orf}_dna_identity",
                f"{orf}_protein_identity", f"{orf}_coverage"
            ])
            # Immediately after the physical `pol` block, surface the separate physical-presence and
            # type-specific translation-compatibility fields so they can never be conflated with the
            # single `pol` verdict downstream.
            if orf == "pol":
                ordered_cols.extend([
                    "pol_region_present", "pol_translation_compatible", "pol_translation_basis"
                ])
        ordered_cols.append("all_orfs_intact")

        final_cols = ordered_cols[:]
        for col in df.columns:
            if col not in final_cols:
                final_cols.append(col)

        final_cols_exist = [col for col in final_cols if col in df.columns]
        df = df[final_cols_exist]
        df.sort_values(by=["ID", "Haplotype", "ID_Full"], inplace=True)
        # Atomic write so a partial/interleaved CSV is never handed to combine_and_validate.py.
        atomic_write_dataframe(df, RESULTS_CSV_FILE, sep=",", index=False, na_rep="NA")
