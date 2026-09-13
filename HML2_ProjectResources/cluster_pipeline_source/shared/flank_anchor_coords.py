#!/usr/bin/env python3
"""Production flank-anchored coordinate helper for process_locus_assemblies_batch.sh.

Replaces the BAM-CIGAR + calculate_padded_coords_pysam.py locating step. Given the
reference locus window and ONE candidate assembly contig (already identified by the
global BAM), it anchors on the UNIQUE, HML-2-free genome flanking the locus and reports
the contig coordinates of the element that sits between the anchors.

Why this fixes the paralog/phantom drift: the provirus body is paralogous (or absent, in
phantom loci) so the old window-overlap mapping landed on the wrong copy or on empty
genome. The flanks are unique per orthologue, so they pin the correct insertion site;
whatever lies between them is the real allele (provirus / solo-LTR / pre-integration).

Maps anchors against the single candidate contig only (chromosome-scale) -> the paralog
discrimination happens within the chromosome where the paralogs actually are, and memory
stays tiny (never indexes the whole assembly).

Output (one tab line to stdout) consumed by the bash:
    STATUS  inner_lo  inner_hi  strand  contig_len
STATUS is "OK" or "FAIL_<reason>". inner_lo/inner_hi are 0-based contig coords of the
expected element span (between the inner anchor edges); strand is + / - (locus orientation
on the contig vs the reference). On FAIL the numeric fields are 0.
"""
import argparse, sys, os, subprocess, tempfile
import pysam
from flank_anchored_locate import clean_anchor, select_pair, anchor_inner_in_asm


def map_pair(anchor_fa, contig_fa, minimap2):
    """Map a 2-record anchor fasta against one contig in a single minimap2 call.
    Returns {anchor_name: [hit dicts]}."""
    p = subprocess.run([minimap2, "-x", "asm20", "-t", "4", "--secondary=yes",
                        "-N", "5", contig_fa, anchor_fa],
                       capture_output=True, text=True)
    out = {}
    for ln in p.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 12:
            continue
        out.setdefault(f[0], []).append(
            dict(qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
                 ctg=f[5], ts=int(f[7]), te=int(f[8]),
                 nmatch=int(f[9]), blen=int(f[10]), mapq=int(f[11])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="reference fasta the locus coords are in (T2T or hg38)")
    ap.add_argument("--kcon", required=True, help="type2_KCON.fa for HML-2-free anchor screening")
    ap.add_argument("--contig-fa", required=True, help="single candidate contig fasta (from BAM qname)")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--inner", type=int, default=500)
    ap.add_argument("--flanks", default="5000,12000,25000")
    ap.add_argument("--expect-pos", type=int, default=-1,
                    help="approx contig position of the locus from the global BAM; when "
                         "set, ambiguous (segdup) anchors are resolved by proximity to it")
    ap.add_argument("--expect-tol", type=int, default=40000)
    args = ap.parse_args()
    expect_pos = None if args.expect_pos < 0 else args.expect_pos

    flanks = [int(x) for x in args.flanks.split(",")]
    ref = pysam.FastaFile(args.ref)
    if args.chrom not in set(ref.references):
        print(f"FAIL_chrom_not_in_ref\t0\t0\t.\t0"); return
    reflen = ref.get_reference_length(args.chrom)
    cfa = pysam.FastaFile(args.contig_fa)
    contig_name = cfa.references[0]
    contig_len = cfa.get_reference_length(contig_name)

    tmp = tempfile.mkdtemp(prefix="fac_")
    pair = None
    why = "no_pair"
    try:
        for fl in flanks:                       # widen anchors until a pair resolves
            a5s, a5e = clean_anchor(ref, args.chrom, max(0, args.start - args.inner),
                                    "5p", fl, args.kcon, args.minimap2, reflen)
            a3s, a3e = clean_anchor(ref, args.chrom, min(reflen, args.end + args.inner),
                                    "3p", fl, args.kcon, args.minimap2, reflen)
            # reference gap between the two anchors' inner edges (element-facing): a5e is
            # the 5' block's inner end, a3s the 3' block's inner start.
            ref_gap = a3s - a5e
            apath = os.path.join(tmp, "anchors.fa")
            with open(apath, "w") as fh:
                fh.write(f">a5\n{ref.fetch(args.chrom, a5s, a5e)}\n")
                fh.write(f">a3\n{ref.fetch(args.chrom, a3s, a3e)}\n")
            hits = map_pair(apath, args.contig_fa, args.minimap2)
            pair = select_pair(hits.get("a5", []), hits.get("a3", []), ref_gap,
                               expect_pos=expect_pos, expect_tol=args.expect_tol)
            if pair:
                break
            why = (f"5p_hits={len(hits.get('a5', []))} 3p_hits={len(hits.get('a3', []))} "
                   f"ref_gap={ref_gap}")
    finally:
        try:
            os.remove(os.path.join(tmp, "anchors.fa"))
        except OSError:
            pass
        try:
            os.rmdir(tmp)
        except OSError:
            pass

    if not pair:
        print(f"FAIL_no_pair\t0\t0\t.\t{contig_len}")
        sys.stderr.write(f"[flank_anchor] {args.chrom}:{args.start}-{args.end} on "
                         f"{contig_name}: {why}\n")
        return
    _h5, _h3, lo, hi, strand = pair
    print(f"OK\t{lo}\t{hi}\t{strand}\t{contig_len}")


if __name__ == "__main__":
    main()
