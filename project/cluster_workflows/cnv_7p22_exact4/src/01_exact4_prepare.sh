#!/usr/bin/env bash
# Exact-four preparation entry point.  It has no submission or retry client.
set -euo pipefail
export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1
umask 027

PLAN=""
RELEASE=""
SAMPLE_INDEX=""
VERIFY_ONLY=false
while (($#)); do
    case "$1" in
        --plan) PLAN=${2:?}; shift 2 ;;
        --sample-index) SAMPLE_INDEX=${2:?}; shift 2 ;;
        --execution-release) RELEASE=${2:?}; shift 2 ;;
        --verify-only) VERIFY_ONLY=true; shift ;;
        *) printf 'unknown argument: %s\n' "$1" >&2; exit 2 ;;
    esac
done
[[ -n "$PLAN" && -n "$RELEASE" && "$SAMPLE_INDEX" =~ ^[0-3]$ ]] || {
    printf 'usage: %s --plan PLAN.json --sample-index 0..3 --execution-release RELEASE.json [--verify-only]\n' "$0" >&2
    exit 2
}
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
PLANNER="$HERE/exact4_plan.py"
PYTHON=${PYTHON:?PYTHON must be the absolute interpreter path pinned by the target-only plan}
[[ "$PYTHON" = /* ]] || { printf 'PYTHON must be an absolute pinned path\n' >&2; exit 2; }
"$PYTHON" "$PLANNER" validate-plan --plan "$PLAN" >/dev/null
mapfile -t RUNTIME_PATHS < <("$PYTHON" - "$PLAN" <<'PY'
import json,sys
p=json.load(open(sys.argv[1],encoding="utf-8")); tools={x["name"]:x["path"] for x in p["immutable_run_material"]["resident_identity"]["tools"]}
print(sys.executable); print(tools["python"]); print(tools["bash"])
PY
)
[[ "${RUNTIME_PATHS[0]}" == "${RUNTIME_PATHS[1]}" && "$BASH" == "${RUNTIME_PATHS[2]}" ]] || {
    printf 'active Python/Bash differ from immutable tool receipt\n' >&2; exit 2
}

"$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
if "$VERIFY_ONLY"; then
    "$PYTHON" "$PLANNER" validate-prepared --plan "$PLAN" --sample-index "$SAMPLE_INDEX" >/dev/null
    exit 0
fi

"$PYTHON" "$PLANNER" claim --plan "$PLAN" --execution-release "$RELEASE" \
    --kind prepare --index "$SAMPLE_INDEX" --token "exact4-prepare-${SAMPLE_INDEX}" >/dev/null

mapfile -t FIELDS < <("$PYTHON" - "$PLAN" "$SAMPLE_INDEX" <<'PY'
import json,sys
p=json.load(open(sys.argv[1],encoding="utf-8")); i=int(sys.argv[2]); m=p["preparation_map"][i]
tools={x["name"]:x["path"] for x in p["immutable_run_material"]["resident_identity"]["tools"]}
rows=[x for x in p["immutable_run_material"]["resident_identity"]["assemblies"] if x["sample_id"]==m["sample_id"]]
print(m["staging_root"]+"/bundle")
print(m["bundle_destination"])
print(rows[0]["fasta"]["path"]); print(rows[1]["fasta"]["path"])
print(tools["minimap2"]); print(tools["samtools"])
PY
)
STAGE=${FIELDS[0]}; DESTINATION=${FIELDS[1]}; FASTA1=${FIELDS[2]}; FASTA2=${FIELDS[3]}
MINIMAP2=${FIELDS[4]}; SAMTOOLS=${FIELDS[5]}

"$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
[[ ! -e "$STAGE" && ! -L "$STAGE" ]] || { printf 'staging path already exists: %s\n' "$STAGE" >&2; exit 2; }
mkdir -m 0750 -- "$STAGE"

HEADERS="$STAGE/.headers"
grep '^>' "$FASTA1" "$FASTA2" | sed 's/^[^:]*:>//; s/^>//; s/[[:space:]].*$//' > "$HEADERS"
[[ -s "$HEADERS" ]] || { printf 'combined reference has no FASTA identifiers\n' >&2; exit 2; }
[[ $(wc -l < "$HEADERS") -eq $(sort -u "$HEADERS" | wc -l) ]] || {
    printf 'duplicate contig identifier across combined haplotypes\n' >&2; exit 2
}
cat -- "$FASTA1" "$FASTA2" > "$STAGE/combined_reference.fa"
"$SAMTOOLS" faidx "$STAGE/combined_reference.fa"
"$MINIMAP2" -d "$STAGE/combined_reference.mmi" "$STAGE/combined_reference.fa"

"$PYTHON" - "$PLAN" "$SAMPLE_INDEX" "$STAGE" <<'PY'
import json,sys
p=json.load(open(sys.argv[1],encoding="utf-8")); i=int(sys.argv[2]); out=sys.argv[3]; sample=p["preparation_map"][i]["sample_id"]
beds=[x for x in p["immutable_run_material"]["resident_identity"]["bed_authorities"] if x["sample_id"]==sample]
core=[]; full=[]
for b in beds:
    for s,e in b["core_intervals"]: core.append((b["contig"],s,e,b["haplotype"]+"_core"))
    for s,e in b["full_intervals"]: full.append((b["contig"],s,e,b["haplotype"]+"_full"))
for name,rows in (("core.bed",core),("full_window.bed",full)):
    with open(out+"/"+name,"x",encoding="ascii",newline="\n") as h:
        for row in sorted(rows): h.write("\t".join(map(str,row))+"\n")
PY
awk 'BEGIN{OFS=""} {print $1,":",$2+1,"-",$3}' "$STAGE/full_window.bed" > "$STAGE/.regions"
"$SAMTOOLS" faidx -r "$STAGE/.regions" "$STAGE/combined_reference.fa" > "$STAGE/bait.fa"
rm -- "$STAGE/.regions"
rm -- "$HEADERS"

PLAN_PATH="$PLAN" SAMPLE_INDEX_VALUE="$SAMPLE_INDEX" STAGE_PATH="$STAGE" FINAL_PATH="$DESTINATION" \
"$PYTHON" - <<'PY'
import hashlib,json,os,pathlib
p=json.load(open(os.environ["PLAN_PATH"],encoding="utf-8")); m=p["preparation_map"][int(os.environ["SAMPLE_INDEX_VALUE"])]
stage=pathlib.Path(os.environ["STAGE_PATH"]); final=pathlib.Path(os.environ["FINAL_PATH"])
names=["combined_reference.fa","combined_reference.fa.fai","combined_reference.mmi","core.bed","full_window.bed","bait.fa"]
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
arts=[{"path":str(final/n),"size_bytes":(stage/n).stat().st_size,"sha256":sha(stage/n)} for n in names]
details={"haplotypes":m["haplotypes"],"core_bed_sha256":sha(stage/"core.bed"),"full_window_bed_sha256":sha(stage/"full_window.bed"),"combined_reference_sha256":sha(stage/"combined_reference.fa"),"combined_fai_sha256":sha(stage/"combined_reference.fa.fai"),"combined_mmi_sha256":sha(stage/"combined_reference.mmi"),"bait_sha256":sha(stage/"bait.fa"),"publication_mode":"same_filesystem_directory_noreplace"}
r={"schema_version":"hml2_7p22_exact4_prepared_1","kind":"prepare","run_id":p["run_id"],"plan_digest":p["plan_digest"],"bundle_id":p["bundle_id"],"index":m["index"],"sample_id":m["sample_id"],"locus":m["locus"],"claim_token":"exact4-prepare-"+str(m["index"]),"artifacts":arts,"scientific_assertions":{"diagnostic_only":True,"complete_unfiltered_alignment_required":True,"competitive_combined_reference_required":True,"primary_only_depth_required":True,"read_cap_or_downsampling_allowed":False},"details":details,"read_inferred_copy_number":"not_estimated","receipt_digest":""}
material=dict(r); material.pop("receipt_digest"); enc=lambda v:(json.dumps(v,sort_keys=True,separators=(",",":"),ensure_ascii=False)+"\n").encode()
r["receipt_digest"]=hashlib.sha256(enc(material)).hexdigest(); (stage/"PREPARED.json").write_bytes(enc(r))
PY
chmod -R a-w -- "$STAGE"
"$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
"$PYTHON" "$PLANNER" publish --plan "$PLAN" --execution-release "$RELEASE" \
    --staged "$STAGE" --destination "$DESTINATION" >/dev/null
"$PYTHON" "$PLANNER" validate-prepared --plan "$PLAN" --sample-index "$SAMPLE_INDEX" >/dev/null
