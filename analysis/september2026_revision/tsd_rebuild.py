import os
from pathlib import Path
INPUT = Path(__file__).resolve().parent / "inputs"
WORK = Path(os.environ["HML2_REVISION_OUTPUT"])
"""Reconcile manuscript record membership with retained paired-boundary evidence."""
from pathlib import Path
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = INPUT
OUT = WORK
CATALOG = INPUT/'tsd/catalog_columns.tsv'
BOUNDARIES = INPUT/'tsd/paired_tsd_boundary_observations_v1.tsv'
BOUNDARY_META = INPUT/'tsd/paired_tsd_boundary_observations_v1.json'
DIRECT = INPUT/'tsd/candidate_tsd_direct_evidence_join.tsv'
DISCREPANCY = INPUT/'tsd/orf_tsd_cross_sample_discrepancy_audit_v1.tsv'
POP = INPUT/'tsd/igsr_samples.tsv'
FIBER = INPUT/'tsd/fiberseq_donor_ids.tsv'

def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write(name, rows, fields=None):
    fields = fields or list(rows[0])
    with (OUT / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def label(locus):
    return locus.removeprefix("HML-2_").removesuffix("_hg38").removesuffix("_new")

def length_comparisons(a, b):
    out = []
    for n in (4, 5, 6):
        left = a[-n:] if n <= len(a) else "?" * (n-len(a)) + a
        right = b[:n] if n <= len(b) else b + "?" * (n-len(b))
        known_differences = sum(x != y for x,y in zip(left,right) if "?" not in (x,y))
        out.append((n, left, right, known_differences, "?" not in left+right))
    return out

catalog = read(CATALOG)
assert len({(r["Locus"], r["ID_Full"]) for r in catalog}) == len(catalog)
boundary_rows = read(BOUNDARIES)
boundary = {}
for r in boundary_rows:
    key = r["Locus"], r["ID_Full"]
    # The source contains replicated alias projections. Collapse only byte-equivalent rows.
    assert key not in boundary or boundary[key] == r, key
    boundary[key] = r
digest = hashlib.sha256(BOUNDARIES.read_bytes()).hexdigest()
assert digest == json.loads(BOUNDARY_META.read_text())["data_file"]["sha256"]
pop = {r["Sample name"]: r["Superpopulation code"] for r in read(POP)}
direct = defaultdict(list)
for r in read(DIRECT):
    direct[(r["candidate_locus"], r["ID_Full"])].append(r)

overrides, unmatched, included = [], [], []
for r in catalog:
    key = r["Locus"], r["ID_Full"]
    original = r["orig_Locus"], r["ID_Full"]
    choices = [boundary[k] for k in dict.fromkeys((key, original)) if k in boundary]
    if len(choices) > 1:
        assert all({k: v for k, v in x.items() if k != "Locus"} ==
                   {k: v for k, v in choices[0].items() if k != "Locus"}
                   for x in choices[1:]), key
    q = choices[0] if choices else None
    in_panel = bool(re.fullmatch(r"(?:HG|NA)\d+", r["ID"])) and r["analysis_include"] == "1"
    o = {
        "Locus": r["Locus"], "ID_Full": r["ID_Full"],
        "orig_Locus": r["orig_Locus"], "ID": r["ID"], "Haplotype": r["Haplotype"],
        "Source_Identifier": r["Source_Identifier"],
        "stable_source_identity_v3": r["stable_source_identity_v3"],
        "analysis_include": r["analysis_include"], "manuscript_public_panel": int(in_panel),
        "catalog_5prime_TSD": r["5'_TSD"], "catalog_3prime_TSD": r["3'_TSD"],
        "boundary_join": ("EXACT_LOCUS_AND_ID" if key in boundary else "EXACT_ORIGINAL_LOCUS_AND_ID") if q else "NO_EXACT_RETAINED_RECORD",
        "boundary_source_locus": q["Locus"] if q else "",
        "boundary_5prime_sequence": q["5'_TSD"] if q else "",
        "boundary_3prime_sequence": q["3'_TSD"] if q else "",
        "pair_admission": q["pair_admission"] if q else "NOT_ASSESSED_NO_EXACT_RETAINED_RECORD",
        "mismatch_count": q["mismatch_count"] if q else "",
        "boundary_presence_call": q["presence_call"] if q else "",
        "boundary_observation_state": q["observation_state"] if q else "",
        "catalog_observation_state": r["observation_state"],
        "superpopulation": pop.get(r["ID"], "Unassigned") or "Unassigned",
        "length_robust_candidate_difference": "",
        "length_ambiguous_difference": "",
        "candidate_length_comparisons": "",
    }
    if q and q["pair_admission"] == "PAIRED_ACCEPTED":
        a, b = q["5'_TSD"], q["3'_TSD"]
        assert len(a) == len(b) and len(a) in (4, 5, 6)
        assert set(a + b) <= set("ACGT")
        assert sum(x != y for x, y in zip(a, b)) == int(q["mismatch_count"])
        candidates = length_comparisons(a, b)
        o["candidate_length_comparisons"] = ";".join(f"{n}:{x}/{y}:{m}" for n,x,y,m,known in candidates)
        o["length_robust_candidate_difference"] = int(all(c[3] > 0 for c in candidates))
        o["length_ambiguous_difference"] = int(int(q["mismatch_count"]) > 0 and not o["length_robust_candidate_difference"])
    overrides.append(o)
    if not q:
        unmatched.append(o)
    if in_panel:
        included.append(o)

assert len(included) == 59656
assert len({(r["ID"], r["Haplotype"]) for r in included}) == 584
write("manuscript_terminal_boundary_measurements.tsv", overrides)
write("manuscript_TSD_unmatched_records.tsv", unmatched)
write("Figure_S7_record_level_boundary_calls.tsv", included)

per_locus, per_pop = defaultdict(list), defaultdict(list)
for r in included:
    per_locus[r["Locus"]].append(r)
    per_pop[r["superpopulation"]].append(r)

def summarize(group, name, rows):
    paired = [r for r in rows if r["pair_admission"] == "PAIRED_ACCEPTED"]
    c = Counter(r["pair_admission"] for r in rows)
    mm = Counter(int(r["mismatch_count"]) for r in paired)
    return {
        group: name, "catalog_records": len(rows),
        "donors": len({r["ID"] for r in rows}),
        "donor_haplotypes": len({(r["ID"], r["Haplotype"]) for r in rows}),
        "paired_accepted": len(paired), "identical_pairs": mm[0],
        "one_mismatch_pairs": mm[1], "two_mismatch_pairs": mm[2],
        "different_pairs": mm[1] + mm[2],
        "fraction_different": (mm[1] + mm[2]) / len(paired) if paired else "",
        "length_robust_candidate_differences": sum(r["length_robust_candidate_difference"] for r in paired),
        "length_ambiguous_differences": sum(r["length_ambiguous_difference"] for r in paired),
        "one_sided": c["ONE_SIDED"], "unavailable": c["UNAVAILABLE"],
        "no_exact_retained_record": c["NOT_ASSESSED_NO_EXACT_RETAINED_RECORD"],
    }

locus_summary = [summarize("locus", l, r) for l, r in sorted(per_locus.items())]
pop_summary = [summarize("superpopulation", p, per_pop[p]) for p in ("AFR", "AMR", "EAS", "EUR", "SAS", "Unassigned")]
write("Figure_S7_TSD_per_locus.tsv", locus_summary)
write("Figure_S7_TSD_population_denominators.tsv", pop_summary)
pair_counts = Counter((r["Locus"], r["boundary_5prime_sequence"], r["boundary_3prime_sequence"], r["mismatch_count"]) for r in included if r["pair_admission"] == "PAIRED_ACCEPTED")
write("Figure_S7_TSD_sequence_pairs.tsv", [dict(locus=l, five_prime=a, three_prime=b, retained_window_differences=int(m), retained_window_length=len(a), count=n,
    length_robust_candidate_difference=int(all(c[3]>0 for c in length_comparisons(a,b))),
    candidate_length_comparisons=";".join(f"{k}:{x}/{y}:{d}" for k,x,y,d,known in length_comparisons(a,b))) for (l,a,b,m),n in sorted(pair_counts.items())])
focus = {"HML-2_3q12.3", "HML-2_3q21.2", "HML-2_4q32.3", "HML-2_5p12", "HML-2_8p23.1a_hg38"}
write("comment90_current_boundary_disposition.tsv", [r for r in locus_summary if r["locus"] in focus])
selected_direct = []
for (l, identity), records in direct.items():
    if l in focus:
        for r in records:
            selected_direct.append({k: r[k] for k in ("candidate_locus", "ID_Full", "tsd_call_status", "evidence_5prime_tsd", "evidence_3prime_tsd", "left_terminal_support", "right_terminal_support", "viral_body_boundary_owner", "review_flag")})
write("comment90_retained_terminal_evidence.tsv", selected_direct)

plt.rcParams.update({"font.family": "Arial", "font.size": 9, "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(7.25, 4.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1.3], wspace=.48)
ax = fig.add_subplot(gs[0])
plot = sorted((r for r in locus_summary if r["paired_accepted"] >= 20 and r["length_robust_candidate_differences"]), key=lambda r: (r["length_robust_candidate_differences"]/r["paired_accepted"], r["locus"]))
ax.barh(range(len(plot)), [r["length_robust_candidate_differences"]/r["paired_accepted"] for r in plot], color="#D55E00", height=.58)
ax.set_yticks(range(len(plot)), [label(r["locus"]) for r in plot])
ax.set_xlim(0, .032)
ax.set_xticks([0, .01, .02, .03])
ax.xaxis.set_major_formatter(PercentFormatter(1))
ax.set_xlabel("Pairs differing at every 4–6 bp length")
ax.set_title("Candidate TSD differences", fontsize=11, pad=14)
for i, r in enumerate(plot):
    ax.text(r["length_robust_candidate_differences"]/r["paired_accepted"] + .0007, i, f'{r["length_robust_candidate_differences"]}/{r["paired_accepted"]}', va="center", fontsize=9)
ax.spines[["right", "top"]].set_visible(False)
ax.text(-.20, 1.04, "A", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax = fig.add_subplot(gs[1]); ax.axis("off")
ax.set_title("8p23.1a depends on length", fontsize=11, pad=14)
ax.text(-.05, 1.04, "B", transform=ax.transAxes, fontweight="bold", fontsize=12)
ax.text(.5, .91, "6-base boundary comparison", ha="center", fontsize=9)
for y, seq, side in ((.79, "ACCTTT", "5′"), (.64, "CCTTTT", "3′")):
    ax.text(.04, y, side, ha="right", va="center", fontsize=10)
    for i, nt in enumerate(seq):
        different = i in (0, 2)
        ax.text(.17 + i*.14, y, nt, ha="center", va="center", family="monospace", fontsize=14,
                color="#D55E00" if different else "#242424", fontweight="bold" if different else "normal")
ax.text(.5, .51, "2 differences", ha="center", fontsize=9)
ax.text(.5, .33, "5-base candidate TSD", ha="center", fontsize=9)
ax.text(.5, .20, "CCTTT / CCTTT", ha="center", va="center", fontsize=12, family="monospace", color="#0072B2")
ax.text(.5, .07, "Identical in all 42 pairs", ha="center", fontsize=9)
fig.subplots_adjust(left=.13, right=.98, bottom=.19, top=.83)
for ext in ("png", "pdf", "svg"):
    fig.savefig(OUT / f"Figure_S6_corrected.{ext}", dpi=450, facecolor="white")
fig.canvas.draw()
for letter, panel in zip(("A", "B"), fig.axes):
    bbox = panel.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted()).expanded(1.05, 1.06)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"Figure_S6{letter}_corrected.{ext}", dpi=450, facecolor="white", bbox_inches=bbox)
plt.close(fig)

ids = sorted({r["Individual_ID"] for r in read(FIBER)})
assert len(ids) == 39 and "HG002" in ids
write("Fiberseq_cell_type_roster.tsv", [dict(Individual_ID=i, cell_type="lymphoblastoid cell line", panel_role="benchmark" if i=="HG002" else "HPRC", source="Lucas et al. 2026, HPRC2; paired long-read chromatin epigenome and Cell Line Production sections") for i in ids])
summary = summarize("scope", "included manuscript public-donor records", included)
summary.update({
    "catalog_all_records": len(catalog), "all_exact_joined_records": len(catalog)-len(unmatched),
    "all_unmatched_records": len(unmatched),
    "identical_duplicate_boundary_rows_collapsed": len(boundary_rows)-len(boundary),
    "join_policy": "Exact current Locus + ID_Full, or exact retained orig_Locus + unchanged ID_Full. Conflicting identities rejected. No fuzzy or donor-only joins.",
    "population_metadata_source": str(POP.relative_to(ROOT)),
    "cell_type_source_url": "https://www.biorxiv.org/content/10.64898/2026.07.21.739710v1.full",
    "cell_type_evidence": "HPRC2 describes Fiber-seq on 38 paired HPRC LCLs and HG002 LCLs, yielding the 39-sample consensus-peak panel. Local retained summary contains HG002 and 38 HPRC IDs.",
    "interpretation": "Boundary measurements are not replacements for catalog TSD annotations. Of 62 differing retained windows, 20 differ at every 4–6 bp candidate length; 42 at 8p23.1a have identical 5-bp candidate TSDs and are length-ambiguous, not demonstrated TSD substitutions.",
    "inputs": [{"path":str(p.relative_to(ROOT)),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in (CATALOG,BOUNDARIES,BOUNDARY_META,DIRECT,DISCREPANCY,POP,FIBER)],
})
(OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
obsolete = OUT / "manuscript_TSD_authoritative_overlay.tsv"
if obsolete.exists():
    obsolete.unlink()
print(json.dumps(summary, indent=2))
