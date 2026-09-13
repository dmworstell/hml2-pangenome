#!/usr/bin/env bash
# Merge one manifest-declared target and produce coverage diagnostics.
# The immutable, plan-hashed static renderer is src/plot_cnv_depth_v2.py.
# Provisional workflow: coverage is diagnostic/assembly-support evidence only.

#SBATCH --partition=batch
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G

set -Eeuo pipefail
IFS=$'\n\t'
PYTHON=${HML2_PYTHON:-python3.11}

die() { printf 'CNV-v2 analysis: %s\n' "$*" >&2; exit 2; }
usage() {
    cat <<'EOF'
Usage: 03_cnv_analysis_v2.sh --plan PLAN.json [--verify-only]

The plan may instead be supplied in CNV_V2_PLAN. Existing state is reused only
when the exact ANALYSIS_DONE receipt and every bound artifact still validate.
EOF
}

plan=${CNV_V2_PLAN:-}
verify_only=false
while (($#)); do
    case "$1" in
        --plan) (($# >= 2)) || die "--plan requires a value"; plan=$2; shift 2 ;;
        --verify-only) verify_only=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done
[[ -n "$plan" && -s "$plan" ]] || die "a nonempty plan is required"
command -v "$PYTHON" >/dev/null 2>&1 || die "CPython 3.11 runtime is required: $PYTHON"
command -v samtools >/dev/null 2>&1 || die "samtools is required"
command -v minimap2 >/dev/null 2>&1 || die "minimap2 is required"
plan=$("$PYTHON" - "$plan" <<'PY'
import os, sys
print(os.path.realpath(sys.argv[1]))
PY
)

# Validate the immutable plan, the current seven-artifact PREPARED receipt, and
# every ordered part receipt before emitting shell-safe assignments. A part is
# accepted only if it binds the same current PREPARED digest and fingerprints.
eval "$("$PYTHON" - "$plan" <<'PY'
import hashlib, json, os, re, shlex, sys

def canonical_digest(value):
    raw=(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=True)+"\n").encode()
    return hashlib.sha256(raw).hexdigest()
def file_sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()
def fp(path, allow_empty=False):
    real=os.path.realpath(path)
    if not os.path.isfile(real) or (not allow_empty and os.path.getsize(real)<=0): raise SystemExit(f"missing/empty artifact: {real}")
    return {"path":real,"size":os.path.getsize(real),"sha256":file_sha(real)}
def qname_fp(path):
    previous=None; count=0
    with open(path,"rb") as handle:
        for number,raw in enumerate(handle,1):
            if not raw.endswith(b"\n"): raise SystemExit(f"QNAME inventory lacks terminal LF: {path}:{number}")
            value=raw[:-1]
            if not value or any(byte<=32 or byte>=127 for byte in value): raise SystemExit(f"invalid QNAME inventory: {path}:{number}")
            if previous is not None and value<=previous: raise SystemExit(f"QNAME inventory is not byte-sorted unique: {path}:{number}")
            previous=value; count+=1
    return {**fp(path,allow_empty=True),"qname_count":count,
            "encoding":"SAM-QNAME ASCII bytes, LF terminated",
            "order":"LC_ALL=C bytewise ascending, unique"}
def receipt(path, schema):
    with open(path,encoding="utf-8") as f: value=json.load(f)
    if value.get("schema_version")!=schema: raise SystemExit(f"receipt schema mismatch: {path}")
    claimed=value.get("receipt_digest")
    material=dict(value); material.pop("receipt_digest",None)
    if not isinstance(claimed,str) or canonical_digest(material)!=claimed:
        raise SystemExit(f"receipt digest mismatch: {path}")
    return value

plan_path=os.path.realpath(sys.argv[1])
with open(plan_path,encoding="utf-8") as f: p=json.load(f)
if p.get("schema_version")!="cnv_v2_plan_1" or p.get("workflow_status")!="provisional" or p.get("production_authorized") is not False:
    raise SystemExit("unsupported or production-authorized plan")
material=dict(p); material.pop("plan_digest",None)
if canonical_digest(material)!=p.get("plan_digest"): raise SystemExit("plan digest mismatch")
plan_file_sha=file_sha(plan_path)
parts=p.get("parts"); expected=p.get("target",{}).get("expected_part_count")
if not isinstance(expected,int) or expected<1 or not isinstance(parts,list) or [x.get("part_id") for x in parts]!=list(range(1,expected+1)):
    raise SystemExit("plan parts are not exactly ordered 1..expected_part_count")
for record in [p.get("manifest"), *p.get("code",{}).values()]:
    if not isinstance(record,dict) or fp(record.get("path",""))["sha256"]!=record.get("sha256"):
        raise SystemExit("plan-bound source/manifest fingerprint mismatch")
ref=p.get("reference",{})
for key in ("full_window_bed","core_body_union_bed","outer_flanks_bed"):
    observed=fp(ref.get(key,""))
    if observed["sha256"]!=ref.get(key+"_sha256"): raise SystemExit(f"plan BED mismatch: {key}")
if fp(ref.get("kcon_fasta",""))["sha256"]!=ref.get("kcon_sha256"): raise SystemExit("KCON digest mismatch")
if fp(p.get("plotter",{}).get("path",""))["sha256"]!=p.get("plotter",{}).get("sha256"):
    raise SystemExit("static plotter source digest mismatch")

prepared=receipt(p["prepared_receipt"],"cnv_v2_prepared_1")
for key,want in {"status":"PREPARED","run_id":p["run_id"],"target_id":p["target_id"],
                  "target_digest":p["target_digest"],"plan_digest":p["plan_digest"],
                  "manifest_sha256":p["manifest"]["sha256"]}.items():
    if prepared.get(key)!=want: raise SystemExit(f"PREPARED binding mismatch: {key}")
seven={name:ref[name] for name in ("combined_fasta","combined_fai","combined_mmi","bait_fasta",
                                    "full_window_bed","core_body_union_bed","outer_flanks_bed")}
if set(prepared.get("artifacts",{}))!=set(seven): raise SystemExit("PREPARED artifact set is not the exact seven")
for name,path in seven.items():
    if prepared["artifacts"].get(name)!=fp(path): raise SystemExit(f"PREPARED fingerprint mismatch: {name}")

part_receipt_digests=[]; bams=[]; bais=[]; qnames=[]; part_claim_digests=[]
for index,part in enumerate(parts):
    done=receipt(part["receipt"],"cnv_v2_part_done_1")
    bindings={"status":"DONE","run_id":p["run_id"],"target_id":p["target_id"],
              "target_digest":p["target_digest"],"plan_digest":p["plan_digest"],
              "manifest_sha256":p["manifest"]["sha256"],"plan_file_sha256":plan_file_sha,
              "part_id":part["part_id"],"part_index":index,"read_uri":part["read_uri"],
              "read_version":part["read_version"],"read_etag":part["read_etag"],
              "read_size":part["read_size"],"read_sha256":part["read_sha256"],
              "prepared_receipt_digest":prepared["receipt_digest"],
              "prepared_artifacts":prepared["artifacts"]}
    for key,want in bindings.items():
        if done.get(key)!=want: raise SystemExit(f"part {part['part_id']} binding mismatch: {key}")
    if set(done.get("outputs",{}))!={"bam","bai","primary_qnames"}: raise SystemExit("part output set mismatch")
    if done["outputs"]["bam"]!=fp(part["bam"]) or done["outputs"]["bai"]!=fp(part["bai"]):
        raise SystemExit(f"part {part['part_id']} output fingerprint mismatch")
    if done["outputs"]["primary_qnames"]!=qname_fp(part["primary_qnames"]):
        raise SystemExit(f"part {part['part_id']} QNAME inventory fingerprint mismatch")
    claim_digest=done.get("claim_digest")
    if not isinstance(claim_digest,str) or not re.fullmatch(r"[0-9a-f]{64}",claim_digest):
        raise SystemExit(f"part {part['part_id']} claim digest is missing")
    part_receipt_digests.append(done["receipt_digest"]); bams.append(part["bam"]); bais.append(part["bai"])
    qnames.append(part["primary_qnames"]); part_claim_digests.append(claim_digest)

out=p.get("outputs",{}); required=("merged_bam","merged_bai","raw_depth_tsv","mapq10_depth_tsv",
                                    "summary_json","raw_pdf","mapq10_pdf","kcon_paf","analysis_receipt")
if set(out)!=set(required): raise SystemExit("analysis output schema mismatch")
for name in required:
    if not os.path.isabs(out[name]): raise SystemExit(f"non-absolute output: {name}")
analysis_dir=os.path.dirname(out["analysis_receipt"])
if any(os.path.commonpath([os.path.realpath(p["target_dir"]),os.path.realpath(v)])!=os.path.realpath(p["target_dir"]) for v in out.values()):
    raise SystemExit("analysis output escapes target directory")
values={"RUN_ID":p["run_id"],"TARGET_ID":p["target_id"],"TARGET_DIGEST":p["target_digest"],
        "PLAN_DIGEST":p["plan_digest"],"PLAN_FILE_SHA256":plan_file_sha,
        "MANIFEST_SHA256":p["manifest"]["sha256"],"PREPARED_RECEIPT_DIGEST":prepared["receipt_digest"],
        "PREPARED_RECEIPT":p["prepared_receipt"],"COMBINED_FASTA":ref["combined_fasta"],
        "FULL_WINDOW_BED":ref["full_window_bed"],"CORE_BODY_BED":ref["core_body_union_bed"],
        "OUTER_FLANKS_BED":ref["outer_flanks_bed"],"KCON_FASTA":ref["kcon_fasta"],
        "PLOTTER":p["plotter"]["path"],"PLOTTER_SHA256":p["plotter"]["sha256"],
        "PLANNER":p["code"]["cnv_plan.py"]["path"],"ANALYSIS_CLAIM":p["analysis_claim"],
        "ANALYSIS_DIR":analysis_dir,**{k.upper():v for k,v in out.items()}}
for key,value in values.items(): print(key+"="+shlex.quote(str(value)))
print("PART_BAMS=("+" ".join(shlex.quote(x) for x in bams)+")")
print("PART_BAIS=("+" ".join(shlex.quote(x) for x in bais)+")")
print("PART_QNAME_INVENTORIES=("+" ".join(shlex.quote(x) for x in qnames)+")")
print("PART_CLAIM_DIGESTS=("+" ".join(shlex.quote(x) for x in part_claim_digests)+")")
print("PART_RECEIPT_DIGESTS=("+" ".join(shlex.quote(x) for x in part_receipt_digests)+")")
print("EXPECTED_PART_COUNT="+str(expected))
PY
)"

for ((i=0;i<EXPECTED_PART_COUNT;i++)); do
    "$PYTHON" "$PLANNER" validate-part --plan "$plan" --part-id "$((i+1))" --samtools samtools >/dev/null \
        || die "part $((i+1)) receipt, claim, BAM, index, or QNAME inventory failed fresh validation"
done
[[ "$MERGED_BAI" == "${MERGED_BAM}.bai" ]] || die "merged BAI must be BAM.bai"

output_paths=("$MERGED_BAM" "$MERGED_BAI" "$RAW_DEPTH_TSV" "$MAPQ10_DEPTH_TSV" "$SUMMARY_JSON" "$RAW_PDF" "$MAPQ10_PDF" "$KCON_PAF" "$ANALYSIS_RECEIPT")
state=0; for path in "${output_paths[@]}"; do [[ -e "$path" || -L "$path" ]] && state=$((state+1)); done

verify_analysis_receipt() {
    "$PYTHON" - "$plan" "$ANALYSIS_RECEIPT" "$PREPARED_RECEIPT_DIGEST" "$ANALYSIS_CLAIM" "${PART_RECEIPT_DIGESTS[@]}" <<'PY'
import hashlib,json,os,sys
def digest(v): return hashlib.sha256((json.dumps(v,sort_keys=True,separators=(",",":"),ensure_ascii=True)+"\n").encode()).hexdigest()
def fp(path):
 h=hashlib.sha256()
 with open(path,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""): h.update(b)
 return {"path":os.path.realpath(path),"size":os.path.getsize(path),"sha256":h.hexdigest()}
p=json.load(open(sys.argv[1],encoding="utf-8")); r=json.load(open(sys.argv[2],encoding="utf-8"))
m=dict(r); claimed=m.pop("receipt_digest",None)
if r.get("schema_version")!="cnv_v2_analysis_done_1" or digest(m)!=claimed: raise SystemExit("analysis receipt digest mismatch")
for k,w in {"status":"DONE","run_id":p["run_id"],"target_id":p["target_id"],"target_digest":p["target_digest"],
            "plan_digest":p["plan_digest"],"manifest_sha256":p["manifest"]["sha256"],
            "prepared_receipt_digest":sys.argv[3]}.items():
 if r.get(k)!=w: raise SystemExit(f"analysis receipt binding mismatch: {k}")
claim=json.load(open(sys.argv[4],encoding="utf-8")); cm=dict(claim); claim_digest=cm.pop("claim_digest",None)
if claim.get("schema_version")!="cnv_v2_claim_1" or digest(cm)!=claim_digest: raise SystemExit("analysis claim digest mismatch")
for k,w in {"status":"CLAIMED","kind":"analysis","run_id":p["run_id"],"target_id":p["target_id"],
            "target_digest":p["target_digest"],"plan_digest":p["plan_digest"],
            "manifest_sha256":p["manifest"]["sha256"]}.items():
 if claim.get(k)!=w: raise SystemExit(f"analysis claim binding mismatch: {k}")
if r.get("claim_digest")!=claim_digest: raise SystemExit("analysis receipt claim binding mismatch")
if r.get("part_receipt_digests")!=sys.argv[5:]: raise SystemExit("analysis receipt part binding mismatch")
expected={k:v for k,v in p["outputs"].items() if k!="analysis_receipt"}
if set(r.get("outputs",{}))!=set(expected): raise SystemExit("analysis receipt output set mismatch")
for k,path in expected.items():
 if r["outputs"].get(k)!=fp(path): raise SystemExit(f"analysis output fingerprint mismatch: {k}")
PY
    local verify_union
    verify_union=$(mktemp "${ANALYSIS_DIR}/.${TARGET_ID}.verify_union.XXXXXX")
    rm -f -- "$verify_union"
    local verify_args=()
    local inventory
    for inventory in "${PART_QNAME_INVENTORIES[@]}"; do verify_args+=(--inventory "$inventory"); done
    "$PYTHON" "$PLANNER" verify-qname-union "${verify_args[@]}" --output-json "$verify_union"
    "$PYTHON" - "$ANALYSIS_RECEIPT" "$verify_union" <<'PY'
import json,sys
receipt=json.load(open(sys.argv[1],encoding="utf-8"))
fresh=json.load(open(sys.argv[2],encoding="utf-8"))
if receipt.get("qname_union")!=fresh: raise SystemExit("analysis receipt QNAME-union binding mismatch")
PY
    rm -f -- "$verify_union"
    samtools quickcheck -v "$MERGED_BAM"
    samtools idxstats "$MERGED_BAM" >/dev/null
    "$PYTHON" - "$RAW_PDF" "$MAPQ10_PDF" <<'PY'
import sys
for name in sys.argv[1:]:
 data=open(name,"rb").read()
 if len(data)<256 or not data.startswith(b"%PDF-") or not data.rstrip().endswith(b"%%EOF"):
  raise SystemExit(f"invalid PDF: {name}")
PY
}

if ((state)); then
    ((state==${#output_paths[@]})) || die "partial analysis state exists; refusing to overwrite it"
    verify_analysis_receipt || die "existing analysis state is stale or invalid"
    printf 'Exact CNV-v2 analysis DONE receipt validated: %s / %s\n' "$RUN_ID" "$TARGET_ID"
    exit 0
fi
[[ "$verify_only" != true ]] || die "--verify-only requested but no completed analysis exists"

mkdir -p -- "$ANALYSIS_DIR"
claimant_token="analysis:${SLURM_JOB_ID:-local}:${HOSTNAME:-unknown}:$$"
ANALYSIS_CLAIM_DIGEST=$("$PYTHON" "$PLANNER" claim --plan "$plan" --kind analysis --claimant-token "$claimant_token") \
    || die "exact analysis claim already exists, is stale, or could not be acquired"
[[ "$ANALYSIS_CLAIM_DIGEST" =~ ^[0-9a-f]{64}$ ]] || die "analysis claim did not return an exact digest"
scratch=$(mktemp -d "${ANALYSIS_DIR}/.${TARGET_ID}.analysis.XXXXXX")
trap 'rm -rf -- "$scratch"' EXIT INT TERM

# Freshly validate the ordered canonical inventories and prove their union is
# pairwise disjoint with a streaming merge before BAM merge or depth work.
qname_union="$scratch/qname_union.json"
union_args=()
for inventory in "${PART_QNAME_INVENTORIES[@]}"; do union_args+=(--inventory "$inventory"); done
"$PYTHON" "$PLANNER" verify-qname-union "${union_args[@]}" --output-json "$qname_union" \
    || die "cross-part primary-QNAME closure failed"
threads=${SLURM_CPUS_PER_TASK:-1}; [[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "invalid CPU count"
((threads>10)) && threads=10
stage_bam="$scratch/merged.bam"
samtools merge -@ "$threads" -o "$stage_bam" -- "${PART_BAMS[@]}"
samtools index -@ "$threads" "$stage_bam"
samtools quickcheck -v "$stage_bam"; samtools idxstats "$stage_bam" >/dev/null
[[ $(samtools view -c -f 0x100 "$stage_bam") == 0 && $(samtools view -c -f 0x800 "$stage_bam") == 0 ]] || die "merged BAM contains non-primary alignments"

raw="$scratch/raw.tsv"; mapq="$scratch/mapq10.tsv"; kcon="$scratch/kcon.paf"; summary="$scratch/summary.json"
samtools depth -a -b "$FULL_WINDOW_BED" "$stage_bam" >"$raw"
samtools depth -a -Q 10 -b "$FULL_WINDOW_BED" "$stage_bam" >"$mapq"
minimap2 -c -t "$threads" "$COMBINED_FASTA" "$KCON_FASTA" >"$kcon"
[[ -s "$raw" && -s "$mapq" ]] || die "a depth output is empty"
[[ -f "$kcon" ]] || die "assembly-resolved KCON output was not created"

"$PYTHON" - "$plan" "$raw" "$mapq" "$kcon" "$summary" <<'PY'
import json,math,statistics,sys
p=json.load(open(sys.argv[1],encoding="utf-8"))
def bed(path):
 out=[]
 for line in open(path,encoding="utf-8"):
  if line.strip():
   f=line.rstrip().split("\t"); out.append((f[0],int(f[1]),int(f[2]),f[3] if len(f)>3 else ""))
 return out
full=bed(p["reference"]["full_window_bed"]); bodies=bed(p["reference"]["core_body_union_bed"]); flanks=bed(p["reference"]["outer_flanks_bed"])
expected={(c,pos) for c,s,e,_ in full for pos in range(s,e)}
def depths(path):
 d={}
 for n,line in enumerate(open(path,encoding="utf-8"),1):
  f=line.rstrip().split("\t")
  if len(f)!=3: raise SystemExit(f"invalid depth row {path}:{n}")
  key=(f[0],int(f[1])-1); value=int(f[2])
  if key in d or value<0: raise SystemExit(f"duplicate/negative depth row {path}:{n}")
  d[key]=value
 if set(d)!=expected: raise SystemExit(f"depth coordinate closure mismatch: {path}")
 return d
raw=depths(sys.argv[2]); mq=depths(sys.argv[3])
if any(mq[k]>raw[k] for k in expected): raise SystemExit("MAPQ>=10 depth exceeds raw depth")
def vals(d, intervals): return [d[(c,pos)] for c,s,e,_ in intervals for pos in range(s,e)]
assemblies=[]; warnings=[]
for a in p["target"]["assemblies"]:
 c=a["combined_contig"]; bi=[x for x in bodies if x[0]==c]; fi=[x for x in flanks if x[0]==c]
 planned_flanks=[x for x in a["flanks"] if x["available_bp"]>0]
 if len(bi)!=1 or len(fi)!=len(planned_flanks): raise SystemExit(f"interval ownership mismatch for {a['assembly_id']}")
 flank_intervals={}
 for planned in planned_flanks:
  expected_name=f"{a['assembly_id']}_{planned['side']}_outer_flank"
  owned=[x for x in fi if x[3]==expected_name]
  if len(owned)!=1 or owned[0][1]!=planned["start"] or owned[0][2]!=planned["end"]:
   raise SystemExit(f"{planned['side']} flank interval mismatch for {a['assembly_id']}")
  flank_intervals[planned["side"]]=owned[0]
 record={"assembly_id":a["assembly_id"],"combined_contig":c,
         "body_interval_0based_half_open":[bi[0][1],bi[0][2]],
         "flank_completeness":a["flank_completeness"],
         "boundary_truncated_sides":a["boundary_truncated_sides"],
         "flank_geometry":a["flanks"],
         "placement_depth_informative":True}
 for name,d in (("raw",raw),("mapq10",mq)):
  bv=vals(d,bi); fv=vals(d,list(flank_intervals.values())); bm=statistics.median(bv); fm=statistics.median(fv)
  if name=="raw" and sum(bv)<=0: raise SystemExit(f"nonzero body depth gate failed for {a['assembly_id']}")
  record[name]={"body_zero_inclusive_median_depth":bm,"own_outer_flank_zero_inclusive_median_depth":fm,
                "body_own_flank_normalized_depth":(bm/fm if fm>0 else None),
                "flank_sides":{
                 planned["side"]:{
                  "requested_bp":planned["requested_bp"],
                  "available_bp":planned["available_bp"],
                  "callable_bp":planned["callable_bp"],
                  "boundary_truncated":planned["boundary_truncated"],
                  "interval_0based_half_open":[planned["start"],planned["end"]],
                  "zero_inclusive_median_depth":statistics.median(vals(d,[flank_intervals[planned["side"]]])),
                  "positive_depth_bp":sum(value>0 for value in vals(d,[flank_intervals[planned["side"]]]))}
                 for planned in planned_flanks}}
  if fm<=0: warnings.append(f"{a['assembly_id']} has zero {name} own-flank depth; pooled flanks do not rescue this assembly-specific scale.")
 record["own_flank_depth_scale"]=record["raw"]["own_outer_flank_zero_inclusive_median_depth"]
 record["body_own_flank_normalized_depth"]=record["raw"]["body_own_flank_normalized_depth"]
 assemblies.append(record)
pooled_raw=statistics.median(vals(raw,flanks)); pooled_mq=statistics.median(vals(mq,flanks))
raw_sum=sum(raw.values()); mq_sum=sum(mq.values()); lost=(raw_sum-mq_sum)/raw_sum if raw_sum else 1.0
warnings.append(f"MAPQ filtering removed {lost:.1%} of summed target-window depth; ambiguous/low-MAPQ alignment remains an interpretation warning.")
# Keep the best KCON alignment for each assembly, explicitly assembly-resolved.
hits={a["assembly_id"]:[] for a in p["target"]["assemblies"]}
for line in open(sys.argv[4],encoding="utf-8"):
 f=line.rstrip().split("\t")
 if len(f)<12: raise SystemExit("invalid KCON PAF row")
 for a in p["target"]["assemblies"]:
  target_start,target_end=int(f[7]),int(f[8])
  # A hit belongs to an assembly only on its exact planned combined contig and
  # only when its target interval overlaps that assembly's planned body.
  if (f[5]==a["combined_contig"] and target_start<a["body_end"]
      and target_end>a["body_start"]):
   hits[a["assembly_id"]].append({"query":f[0],"target":f[5],"target_start":target_start,"target_end":target_end,
                                  "matching_bases":int(f[9]),"alignment_block_length":int(f[10]),"mapq":int(f[11])})
for rec in assemblies:
 hs=sorted(hits[rec["assembly_id"]],key=lambda x:(x["matching_bases"],x["mapq"]),reverse=True)
 rec["assembly_resolved_kcon_best_alignment"]=hs[0] if hs else None
 if not hs: warnings.append(f"No KCON alignment was reported for {rec['assembly_id']}.")
 elif len(hs)>1 and hs[1]["matching_bases"]>=0.95*hs[0]["matching_bases"]: warnings.append(f"{rec['assembly_id']} has near-tied KCON alignments; locus assignment may be ambiguous.")
summary={"schema_version":"cnv_v2_coverage_summary_1","run_id":p["run_id"],"target_id":p["target_id"],
 "sample_id":p["target"]["sample_id"],"locus_id":p["target"]["locus_id"],
 "coverage_role":"diagnostic/assembly-support only","read_inferred_copy_number":"not_estimated",
 "normalization_statement":"Each combined-reference haplotype uses every available base from its own named outer flank side or sides as a one-copy depth scale; a missing boundary side is neither padded nor imputed, and ratios are not divided by two.",
 "boundary_flank_rule":"ONE_SIDED_BOUNDARY_FLANK is informative when exactly one non-empty callable flank remains. Shortened boundary flanks contribute only their actual callable bases. Body depth alone never rescues a no-flank observation.",
 "pooled_outer_flank_baseline":pooled_raw,
 "pooled_outer_flank_baselines":{"raw_zero_inclusive_median_depth":pooled_raw,"mapq10_zero_inclusive_median_depth":pooled_mq,"descriptive_only":True,"may_rescue_failed_haplotype":False},
 "assembly_summaries":assemblies,"warnings":warnings}
with open(sys.argv[5],"x",encoding="utf-8") as f: json.dump(summary,f,sort_keys=True,indent=2); f.write("\n")
PY

raw_pdf="$scratch/raw.pdf"; mapq_pdf="$scratch/mapq10.pdf"
[[ $("$PYTHON" - "$PLOTTER" <<'PY'
import hashlib,sys
h=hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest(); print(h)
PY
) == "$PLOTTER_SHA256" ]] || die "static plotter changed before rendering"
"$PYTHON" "$PLOTTER" --depth-tsv "$raw" --output-pdf "$raw_pdf" --metric raw --summary-json "$summary" --title "$TARGET_ID coverage diagnostic"
"$PYTHON" "$PLOTTER" --depth-tsv "$mapq" --output-pdf "$mapq_pdf" --metric mapq10 --summary-json "$summary" --title "$TARGET_ID coverage diagnostic"

for pdf in "$raw_pdf" "$mapq_pdf"; do
 "$PYTHON" - "$pdf" <<'PY'
import sys
d=open(sys.argv[1],"rb").read()
if len(d)<256 or not d.startswith(b"%PDF-") or not d.rstrip().endswith(b"%%EOF"): raise SystemExit("invalid PDF")
PY
done

# Publish canonical outputs with same-filesystem create-only links, then the
# self-digested DONE receipt last.  No destination can be overwritten.
"$PYTHON" "$PLANNER" publish --staged "$stage_bam" --destination "$MERGED_BAM" >/dev/null
"$PYTHON" "$PLANNER" publish --staged "${stage_bam}.bai" --destination "$MERGED_BAI" >/dev/null
"$PYTHON" "$PLANNER" publish --staged "$raw" --destination "$RAW_DEPTH_TSV" >/dev/null
"$PYTHON" "$PLANNER" publish --staged "$mapq" --destination "$MAPQ10_DEPTH_TSV" >/dev/null
"$PYTHON" "$PLANNER" publish --staged "$summary" --destination "$SUMMARY_JSON" >/dev/null
"$PYTHON" "$PLANNER" publish --staged "$raw_pdf" --destination "$RAW_PDF" >/dev/null
"$PYTHON" "$PLANNER" publish --staged "$mapq_pdf" --destination "$MAPQ10_PDF" >/dev/null
"$PYTHON" "$PLANNER" publish --staged "$kcon" --destination "$KCON_PAF" >/dev/null
receipt_stage="$scratch/ANALYSIS_DONE.json"
"$PYTHON" - "$plan" "$PREPARED_RECEIPT_DIGEST" "$ANALYSIS_CLAIM_DIGEST" "$qname_union" "$receipt_stage" "${PART_RECEIPT_DIGESTS[@]}" <<'PY'
import hashlib,json,os,sys
def digest(v): return hashlib.sha256((json.dumps(v,sort_keys=True,separators=(",",":"),ensure_ascii=True)+"\n").encode()).hexdigest()
def fp(path):
 h=hashlib.sha256()
 with open(path,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""): h.update(b)
 return {"path":os.path.realpath(path),"size":os.path.getsize(path),"sha256":h.hexdigest()}
p=json.load(open(sys.argv[1],encoding="utf-8")); expected=p["target"]["expected_part_count"]
qname_union=json.load(open(sys.argv[4],encoding="utf-8")); part_digests=sys.argv[6:]
if len(part_digests)!=expected: raise SystemExit("part receipt digest count mismatch")
r={"schema_version":"cnv_v2_analysis_done_1","status":"DONE","run_id":p["run_id"],"target_id":p["target_id"],
   "target_digest":p["target_digest"],"plan_digest":p["plan_digest"],"manifest_sha256":p["manifest"]["sha256"],
   "prepared_receipt_digest":sys.argv[2],"claim_digest":sys.argv[3],"part_receipt_digests":part_digests,
   "qname_union":qname_union,
   "outputs":{k:fp(v) for k,v in p["outputs"].items() if k!="analysis_receipt"}}
r["receipt_digest"]=digest(r)
with open(sys.argv[5],"x",encoding="utf-8") as f: json.dump(r,f,sort_keys=True,indent=2); f.write("\n")
PY
"$PYTHON" "$PLANNER" publish --staged "$receipt_stage" --destination "$ANALYSIS_RECEIPT" >/dev/null
verify_analysis_receipt
printf 'CNV-v2 coverage diagnostics complete: %s / %s\n' "$RUN_ID" "$TARGET_ID"
printf 'Read-inferred copy number: not estimated. Nucleotide-frequency analysis was not invoked.\n'
