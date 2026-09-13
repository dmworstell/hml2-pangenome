#!/usr/bin/env python3
"""Flank-anchored locator for HML-2 extraction.

Instead of mapping the (often paralogous / phantom) provirus window into an assembly,
anchor on the UNIQUE genome flanking the locus: pull a block of reference sequence just
5' and just 3' of [START,END] (held off the element by an inner buffer so the anchor
never includes the LTR), map each block to the assembly, and take whatever sits BETWEEN
the two anchors. That span is the element (provirus / solo-LTR) when present, or an empty
pre-integration junction when absent -- which is the correct answer for insertionally
polymorphic / phantom loci where the reference window itself holds no element.

Self-contained: shells out to minimap2, parses PAF. No samtools needed (uses pysam).
"""
import argparse, subprocess, sys, tempfile, os
import pysam

def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]

def map_anchor(anchor_fa, asm_target, minimap2):
    """Map one anchor block to the assembly; return list of PAF hit dicts.
    asm_target may be a FASTA or a prebuilt minimap2 .mmi index (much faster)."""
    p = subprocess.run([minimap2, "-x", "asm20", "-t", "4", "--secondary=yes",
                        "-N", "5", asm_target, anchor_fa],
                       capture_output=True, text=True)
    hits = []
    for ln in p.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 12:
            continue
        hits.append(dict(qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
                         ctg=f[5], ts=int(f[7]), te=int(f[8]),
                         nmatch=int(f[9]), blen=int(f[10]), mapq=int(f[11])))
    return hits

def kcon_hits_in(ref, chrom, a, b, kcon_ref, minimap2):
    """Return list of (hs,he) KCON-homologous intervals (ref coords) within [a,b]."""
    if not kcon_ref:
        return []
    import tempfile, os
    seq = ref.fetch(chrom, max(0, a), b)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as t:
        t.write(f">blk\n{seq}\n"); fn = t.name
    out = subprocess.run([minimap2, "-x", "asm20", "-p", "0.05", "-N", "20", "-t", "4",
                          kcon_ref, fn], capture_output=True, text=True).stdout
    os.remove(fn)
    iv = []
    for ln in out.splitlines():
        f = ln.split("\t")
        if int(f[10]) >= 200:                       # ignore tiny incidental hits
            iv.append((a + int(f[2]), a + int(f[3])))
    return iv

def clean_anchor(ref, chrom, inner_edge, side, fl, kcon_ref, minimap2,
                 reflen, max_walk=40000, min_gap_anchor=3500):
    """Return (a,b) for a KCON-free anchor block of length ~fl on the given side
    ('5p' -> block ends at inner_edge going left; '3p' -> starts at inner_edge going
    right), sliding OUTWARD past any HML-2/LTR homology so the anchor is unique genome.

    Neighbour handling depends on the unique gap between that HML-2 and the locus:
      * WIDE gap (>= min_gap_anchor): the neighbour is a genuinely SEPARATE element
        (e.g. the solo-LTR 4.3kb 5' of the 4q32.3 degraded provirus; measured flank gap
        3835bp). Anchor INSIDE the clean gap, keeping the inner edge at the locus, so the
        neighbour is NOT pulled into the extracted span. Threshold 3500 is set ABOVE the
        largest real single-element internal HML-2-free gap seen across extractions_longseqs
        (4q32.3 gag-pol deletion ~3.1kb, 3q12.3 ~2.8kb) so a genuine degraded provirus is
        never mistaken for two elements, and BELOW the 4q32.3 neighbour gap (3835bp).
      * NARROW gap (< min_gap_anchor): the neighbour abuts the locus -- typically just
        slightly-off coordinates -- so slide the whole anchor outward past it, which
        also pulls that adjacent LTR into the span so a real element is not missed."""
    walked = 0
    while walked <= max_walk:
        if side == "5p":
            a, b = max(0, inner_edge - fl), inner_edge
        else:
            a, b = inner_edge, min(reflen, inner_edge + fl)
        iv = kcon_hits_in(ref, chrom, a, b, kcon_ref, minimap2)
        if not iv:
            return a, b
        # push the inner edge outward, just past the HML-2 hit nearest the element
        if side == "5p":
            # element-facing edge of the HML-2 nearest the locus = highest inner end.
            nearest_inner = max(h[1] for h in iv)
            if inner_edge - nearest_inner >= min_gap_anchor:
                return nearest_inner + 200, inner_edge   # clean gap abutting the locus
            new_inner = min(h[0] for h in iv) - 200
            if new_inner >= inner_edge:           # no progress -> give up, return as-is
                return a, b
            walked += inner_edge - new_inner; inner_edge = new_inner
        else:
            nearest_inner = min(h[0] for h in iv)
            if nearest_inner - inner_edge >= min_gap_anchor:
                return inner_edge, nearest_inner - 200   # clean gap abutting the locus
            new_inner = max(h[1] for h in iv) + 200
            if new_inner <= inner_edge:
                return a, b
            walked += new_inner - inner_edge; inner_edge = new_inner
    if side == "5p":
        return max(0, inner_edge - fl), inner_edge
    return inner_edge, min(reflen, inner_edge + fl)

def best_unique(hits, min_mapq=20, min_frac=0.5, expect_pos=None, expect_tol=40000):
    """Pick the single best anchor hit, else None.

    Two modes:
    * Global (expect_pos=None): the anchor must be globally unique -- top hit mapq>=
      min_mapq, covers >=min_frac of the anchor, and clearly beats the runner-up. Used
      when there is no positional prior (standalone loci, or local testing without a BAM).
    * Localized (expect_pos set): the global BAM has already pinned the locus's neighbour-
      hood, so restrict hits to within expect_tol of expect_pos and take the best-covered
      one there. mapq is NOT required (in a segmental duplication the correct flank gets
      mapq 0 precisely because identical paralogs exist elsewhere -- the BAM prior is what
      disambiguates, not mapq). This is what makes the dense subtelomeric clusters
      (8p23.1 / 1p36.21 / Xq28) resolvable."""
    if not hits:
        return None, "no_hit"
    hits = sorted(hits, key=lambda h: -h["nmatch"])
    if expect_pos is not None:
        near = [h for h in hits
                if abs(((h["ts"] + h["te"]) // 2) - expect_pos) <= expect_tol]
        if not near:
            return None, "none_near_expect"
        near.sort(key=lambda h: -h["nmatch"])
        top = near[0]
        if top["nmatch"] < min_frac * top["qlen"]:
            return None, f"low_cov({top['nmatch']}/{top['qlen']})"
        return top, "ok_localized"
    top = hits[0]
    if top["mapq"] < min_mapq:
        return None, f"low_mapq({top['mapq']})"
    if top["nmatch"] < min_frac * top["qlen"]:
        return None, f"low_cov({top['nmatch']}/{top['qlen']})"
    if len(hits) > 1 and hits[1]["nmatch"] > 0.66 * top["nmatch"]:
        return None, "ambiguous(2nd_hit_close)"
    return top, "ok"

def select_pair(h5s, h3s, ref_gap, expect_pos=None, expect_tol=40000,
                min_ins=-12000, max_ins=40000, min_frac=0.5):
    """Joint 5'/3' anchor-pair selection for dense / closely-spaced duplicate clusters.

    Independently picking the best 5' and best 3' anchor fails when the locus sits in a
    tandem segmental duplication: the 5' anchor may land on copy A while the 3' lands on
    copy B, bracketing a chimaeric span (or mismatched strands). Instead, score every
    (5',3') hit PAIR and keep the one that:
      - shares a strand,
      - brackets an assembly span close to the reference gap between the anchors
        (ref_gap + insertion / - small deletion) -> both anchors are from the SAME copy,
      - lies nearest the BAM positional prior (expect_pos).
    ref_gap = reference distance between the 5' anchor's inner edge and the 3' anchor's
    inner edge. A cross-copy pair inflates the span by the copy spacing; proximity to a
    precise expect_pos is what ultimately resolves near-identical adjacent copies.
    Returns (h5, h3, lo, hi, strand) or None."""
    def good(h):
        return h["nmatch"] >= min_frac * h["qlen"]
    h5s = [h for h in h5s if good(h)]
    h3s = [h for h in h3s if good(h)]
    best = None
    best_key = None
    for h5 in h5s:
        for h3 in h3s:
            if h5["strand"] != h3["strand"]:
                continue
            p5 = anchor_inner_in_asm(h5, "5p")
            p3 = anchor_inner_in_asm(h3, "3p")
            lo, hi = sorted((p5, p3))
            ins = (hi - lo) - ref_gap
            if ins < min_ins or ins > max_ins:
                continue
            mid = (lo + hi) // 2
            prox = abs(mid - expect_pos) if expect_pos is not None else 0
            if expect_pos is not None and prox > expect_tol:
                continue
            # primary: nearest the prior; secondary: smallest span anomaly (favours the
            # tight within-copy bracket over a sprawling cross-copy one on ties).
            key = (prox, abs(ins))
            if best_key is None or key < best_key:
                best_key = key
                best = (h5, h3, lo, hi, h5["strand"])
    return best


def anchor_inner_in_asm(hit, which):
    """Assembly coordinate of the anchor edge that faces the element.
    The anchor block is oriented 5'->3' along the reference. 'which' is "5p" or "3p"
    (the reference side). On a + strand mapping the inner edge of the 5' anchor is its
    te; on - strand the reference 5' anchor maps reversed so its inner edge is ts."""
    if which == "5p":
        return hit["te"] if hit["strand"] == "+" else hit["ts"]
    else:  # 3p anchor inner edge
        return hit["ts"] if hit["strand"] == "+" else hit["te"]

def locate(ref, asm_fa, chrom, start, end, minimap2,
           inner=500, flanks=(5000, 12000, 25000), pad=1500,
           max_span=80000, asm_idx=None, kcon_ref=None, expect_pos=None, expect_tol=40000):
    """Return dict(contig, qstart, qend, strand, span, status, detail).
    asm_idx = optional prebuilt minimap2 .mmi (mapping target); asm_fa is still used
    for pysam length/fetch. Falls back to asm_fa for mapping if asm_idx is None.
    kcon_ref = type2_KCON.fa; if given, anchors are screened HML-2-free (slid outward
    past neighboring LTRs in dense clusters)."""
    asm_target = asm_idx or asm_fa
    reflen = ref.get_reference_length(chrom)
    tmp = tempfile.mkdtemp(prefix="flank_")
    a5_hit = a3_hit = None
    used5 = used3 = None
    try:
        for fl in flanks:  # widen until both anchors are unique
            a5s, a5e = clean_anchor(ref, chrom, max(0, start - inner), "5p", fl,
                                    kcon_ref, minimap2, reflen)
            a3s, a3e = clean_anchor(ref, chrom, min(reflen, end + inner), "3p", fl,
                                    kcon_ref, minimap2, reflen)
            f5 = os.path.join(tmp, "a5.fa"); f3 = os.path.join(tmp, "a3.fa")
            open(f5, "w").write(f">a5\n{ref.fetch(chrom, a5s, a5e)}\n")
            open(f3, "w").write(f">a3\n{ref.fetch(chrom, a3s, a3e)}\n")
            if a5_hit is None:
                h, why5 = best_unique(map_anchor(f5, asm_target, minimap2),
                                      expect_pos=expect_pos, expect_tol=expect_tol)
                if h: a5_hit, used5 = h, fl
            if a3_hit is None:
                h, why3 = best_unique(map_anchor(f3, asm_target, minimap2),
                                      expect_pos=expect_pos, expect_tol=expect_tol)
                if h: a3_hit, used3 = h, fl
            if a5_hit and a3_hit:
                break
    finally:
        for fn in ("a5.fa", "a3.fa"):
            try: os.remove(os.path.join(tmp, fn))
            except OSError: pass
        try: os.rmdir(tmp)
        except OSError: pass

    if not a5_hit and not a3_hit:
        return dict(status="FAIL_both_anchors", detail=f"5p:{why5} 3p:{why3}")
    if not a5_hit or not a3_hit:
        miss = "5p" if not a5_hit else "3p"
        return dict(status="FAIL_one_anchor", detail=f"missing {miss} (5p:{why5} 3p:{why3})")
    if a5_hit["ctg"] != a3_hit["ctg"]:
        return dict(status="FAIL_split_contigs",
                    detail=f"5p->{a5_hit['ctg']} 3p->{a3_hit['ctg']}")
    if a5_hit["strand"] != a3_hit["strand"]:
        return dict(status="FAIL_strand_mismatch", detail="")
    ctg = a5_hit["ctg"]; strand = a5_hit["strand"]
    p5 = anchor_inner_in_asm(a5_hit, "5p")
    p3 = anchor_inner_in_asm(a3_hit, "3p")
    lo, hi = sorted((p5, p3))
    span = hi - lo
    if span > max_span:
        return dict(status="FAIL_span_too_big", detail=f"span={span}")
    ctglen = ref_asm_len(asm_fa, ctg)
    qstart = max(0, lo - pad); qend = min(ctglen, hi + pad)
    return dict(status="ok", contig=ctg, qstart=qstart, qend=qend, strand=strand,
                span=span, inner_lo=lo, inner_hi=hi,
                flank5=used5, flank3=used3, detail="")

_ASMIDX = {}
def ref_asm_len(asm_fa, ctg):
    if asm_fa not in _ASMIDX:
        _ASMIDX[asm_fa] = pysam.FastaFile(asm_fa)
    return _ASMIDX[asm_fa].get_reference_length(ctg)

def extract(asm_fa, ctg, qstart, qend, strand):
    fa = pysam.FastaFile(asm_fa)
    seq = fa.fetch(ctg, qstart, qend)
    return revcomp(seq) if strand == "-" else seq

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--asm", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--inner", type=int, default=500)
    ap.add_argument("--pad", type=int, default=1500)
    ap.add_argument("--out", default="")
    ap.add_argument("--name", default="element")
    args = ap.parse_args()
    ref = pysam.FastaFile(args.ref)
    res = locate(ref, args.asm, args.chrom, args.start, args.end, args.minimap2,
                 inner=args.inner, pad=args.pad)
    if res["status"] != "ok":
        sys.stderr.write(f"{args.name}\t{res['status']}\t{res.get('detail','')}\n")
        sys.exit(2)
    seq = extract(args.asm, res["contig"], res["qstart"], res["qend"], res["strand"])
    hdr = (f">{args.name} {res['contig']}:{res['qstart']}-{res['qend']} strand={res['strand']} "
           f"span={res['span']} flanks={res['flank5']}/{res['flank3']}")
    if args.out:
        open(args.out, "w").write(hdr + "\n" + seq + "\n")
    else:
        sys.stdout.write(hdr + "\n" + seq + "\n")
    sys.stderr.write(f"{args.name}\tok\tspan={res['span']}\t{res['contig']}:"
                     f"{res['qstart']}-{res['qend']}\tstrand={res['strand']}\n")

if __name__ == "__main__":
    main()
