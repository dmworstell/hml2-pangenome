#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates padded query coordinates (1-based) corresponding to a specific
target *reference* locus, given an alignment's CIGAR string and reference
start position. Uses pysam for robust CIGAR parsing and coordinate mapping.

Version 6 (orientation_fixed):
- Correctly calculates the locus orientation on the query strand.
- Finds the query coordinates mapping to the min/max reference positions
  within the target locus to determine if the mapping is inverted.
"""
import argparse
import sys
import traceback

try:
    import pysam
except ImportError:
    print("FATAL: Required Python module 'pysam' not found.", file=sys.stderr)
    sys.exit(1)

def calculate_coords_robust(cigar_string, ref_start_0based, target_ref_start_1based, target_ref_end_1based, padding):
    """
    Robustly calculates 1-based query coordinates and orientation for a target reference locus.
    """
    if padding < 0:
        raise ValueError("Padding must be non-negative.")
    if not cigar_string or cigar_string == '*':
        raise ValueError("Cannot calculate coordinates from empty or unmapped CIGAR string.")

    # Convert to 0-based for internal calculations
    target_ref_start_0based = target_ref_start_1based - 1
    target_ref_end_0based = target_ref_end_1based - 1

    # Use a dummy AlignedSegment to leverage pysam's CIGAR parser
    segment = pysam.AlignedSegment()
    segment.cigarstring = cigar_string
    segment.reference_start = ref_start_0based

    initial_hard_clip_len = 0
    if segment.cigartuples and segment.cigartuples[0][0] == pysam.CHARD_CLIP:
        initial_hard_clip_len = segment.cigartuples[0][1]

    min_q_idx_rel = float('inf')
    max_q_idx_rel = float('-inf')
    
    # Store the query index that corresponds to the first and last aligned base of the locus
    q_idx_at_min_ref = -1
    q_idx_at_max_ref = -1
    min_ref_idx_found = float('inf')
    max_ref_idx_found = float('-inf')
    found_overlap = False

    for q_idx, r_idx in segment.get_aligned_pairs(matches_only=True):
        if r_idx >= target_ref_start_0based and r_idx <= target_ref_end_0based:
            found_overlap = True
            
            # Update overall query span
            if q_idx < min_q_idx_rel: min_q_idx_rel = q_idx
            if q_idx > max_q_idx_rel: max_q_idx_rel = q_idx

            # Find query base mapping to the start-most part of the locus
            if r_idx < min_ref_idx_found:
                min_ref_idx_found = r_idx
                q_idx_at_min_ref = q_idx
            
            # Find query base mapping to the end-most part of the locus
            if r_idx > max_ref_idx_found:
                max_ref_idx_found = r_idx
                q_idx_at_max_ref = q_idx

    if not found_overlap:
        raise ValueError(f"Alignment CIGAR ({cigar_string[:60]}...) does not have any bases mapping inside the target reference locus {target_ref_start_1based}-{target_ref_end_1based}.")

    # *** FIX: Determine orientation ***
    # If the query position at the start of the locus is greater than the query position
    # at the end of the locus, the feature is inverted on the query contig.
    if q_idx_at_min_ref > q_idx_at_max_ref:
        locus_orientation_flag = "1" # Reversed
    else:
        locus_orientation_flag = "0" # Forward

    # Unpadded, absolute coordinates on the original query contig
    raw_query_locus_start_1based = initial_hard_clip_len + min_q_idx_rel + 1
    raw_query_locus_end_1based   = initial_hard_clip_len + max_q_idx_rel + 1

    # Apply padding
    padded_qstart = raw_query_locus_start_1based - padding
    padded_qend = raw_query_locus_end_1based + padding

    if padded_qstart < 1:
        padded_qstart = 1
        
    return padded_qstart, padded_qend, locus_orientation_flag, raw_query_locus_start_1based, raw_query_locus_end_1based

def main():
    parser = argparse.ArgumentParser(description="Calculate query coordinates from CIGAR.")
    parser.add_argument("--ref_start", type=int, required=True, help="1-based start position of the alignment (BAM POS field).")
    parser.add_argument("--target_ref_start", type=int, required=True, help="1-based start of target locus.")
    parser.add_argument("--target_ref_end", type=int, required=True, help="1-based end of target locus.")
    parser.add_argument("-p", "--padding", type=int, default=5000, help="Padding to add.")
    
    # This script will now read the CIGAR from stdin for simplicity with the calling script
    args = parser.parse_args()

    try:
        cigar_string = sys.stdin.readline().strip()
        if not cigar_string:
            raise IOError("CIGAR string input from stdin is empty.")
            
        ref_start_0based = args.ref_start - 1

        padded_start, padded_end, locus_orientation_flag, raw_start, raw_end = calculate_coords_robust(
            cigar_string, ref_start_0based, args.target_ref_start, args.target_ref_end, args.padding
        )
        print(f"{padded_start} {padded_end} {locus_orientation_flag} {raw_start} {raw_end}")

    except Exception as e:
        print(f"FATAL ERROR in calculate_padded_coords_pysam.py: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()