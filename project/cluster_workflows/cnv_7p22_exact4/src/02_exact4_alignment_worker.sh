#!/usr/bin/env bash
# Exact-four complete-unfiltered competitive alignment worker.
set -euo pipefail
export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1
umask 027

PLAN=""; RELEASE=""; PART_INDEX=""; VERIFY_ONLY=false
while (($#)); do
    case "$1" in
        --plan) PLAN=${2:?}; shift 2 ;;
        --part-index) PART_INDEX=${2:?}; shift 2 ;;
        --execution-release) RELEASE=${2:?}; shift 2 ;;
        --verify-only) VERIFY_ONLY=true; shift ;;
        *) printf 'unknown argument: %s\n' "$1" >&2; exit 2 ;;
    esac
done
[[ -n "$PLAN" && -n "$RELEASE" && "$PART_INDEX" =~ ^([0-9]|1[0-2])$ ]] || {
    printf 'usage: %s --plan PLAN.json --part-index 0..12 --execution-release RELEASE.json [--verify-only]\n' "$0" >&2
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
    "$PYTHON" "$PLANNER" validate-part --plan "$PLAN" --part-index "$PART_INDEX" >/dev/null
    exit 0
fi

"$PYTHON" "$PLANNER" claim --plan "$PLAN" --execution-release "$RELEASE" \
    --kind part --index "$PART_INDEX" --token "exact4-part-${PART_INDEX}" >/dev/null

mapfile -t FIELDS < <("$PYTHON" - "$PLAN" "$PART_INDEX" <<'PY'
import json,sys,urllib.parse
p=json.load(open(sys.argv[1],encoding="utf-8")); m=p["alignment_map"][int(sys.argv[2])]
tools={x["name"]:x["path"] for x in p["immutable_run_material"]["resident_identity"]["tools"]}
o=m["provider_object"]; u=urllib.parse.urlsplit(o["uri"])
prep=p["preparation_map"][m["sample_index"]]
for v in (m["staging_root"]+"/work",m["node_scratch_basename"],u.netloc,urllib.parse.unquote(u.path.lstrip("/")),o["version_id"],o["etag"],str(o["content_length_bytes"]),tools["aws"],tools["minimap2"],tools["samtools"],prep["destinations"][0],prep["destinations"][2],prep["destinations"][4],m["unfiltered_bam"],m["primary_bam"],m["primary_bai"],m["source_qname_inventory"],m["decoded_qname_inventory"],m["aligned_primary_qname_inventory"],m["supplementary_qname_inventory"],m["receipt"]): print(v)
PY
)
STAGE=${FIELDS[0]}; SCRATCH_NAME=${FIELDS[1]}; BUCKET=${FIELDS[2]}; KEY=${FIELDS[3]}
VERSION_ID=${FIELDS[4]}; EXPECTED_ETAG=${FIELDS[5]}; EXPECTED_LENGTH=${FIELDS[6]}
AWS=${FIELDS[7]}; MINIMAP2=${FIELDS[8]}; SAMTOOLS=${FIELDS[9]}
COMBINED_FA=${FIELDS[10]}; COMBINED_MMI=${FIELDS[11]}; FULL_BED=${FIELDS[12]}
FINAL_UNFILTERED=${FIELDS[13]}; FINAL_PRIMARY=${FIELDS[14]}; FINAL_BAI=${FIELDS[15]}
FINAL_SOURCE_QNAMES=${FIELDS[16]}; FINAL_DECODED_QNAMES=${FIELDS[17]}
FINAL_ALIGNED_QNAMES=${FIELDS[18]}; FINAL_SUPPLEMENTARY_QNAMES=${FIELDS[19]}
FINAL_RECEIPT=${FIELDS[20]}

"$PYTHON" "$PLANNER" validate-prepared --plan "$PLAN" --sample-index "$("$PYTHON" -c 'import json,sys; p=json.load(open(sys.argv[1])); print(p["alignment_map"][int(sys.argv[2])]["sample_index"])' "$PLAN" "$PART_INDEX")" >/dev/null
"$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
[[ ! -e "$STAGE" && ! -L "$STAGE" ]] || { printf 'staging path already exists: %s\n' "$STAGE" >&2; exit 2; }
mkdir -m 0750 -- "$STAGE"

WORKER_USER=${USER:-$(id -un)}
NODE_TMP=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}
[[ "$NODE_TMP" = /* && -d "$NODE_TMP" && -w "$NODE_TMP" ]] || {
    printf 'scratch parent must be an absolute writable directory: %s\n' "$NODE_TMP" >&2
    exit 2
}
SCRATCH_ROOT=${NODE_TMP%/}/$WORKER_USER
install -d -m 0700 -- "$SCRATCH_ROOT"
SCRATCH="$SCRATCH_ROOT/$SCRATCH_NAME"
[[ ! -e "$SCRATCH" && ! -L "$SCRATCH" ]] || { printf 'scratch source path already exists: %s\n' "$SCRATCH" >&2; exit 2; }
cleanup_scratch() {
    rm -f -- "$SCRATCH"
}
trap cleanup_scratch EXIT INT TERM
"$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
"$AWS" s3api head-object --no-sign-request --bucket "$BUCKET" --key "$KEY" --version-id "$VERSION_ID" > "$STAGE/provider_head.json"
HEAD_PATH="$STAGE/provider_head.json" EXPECTED_ETAG_VALUE="$EXPECTED_ETAG" EXPECTED_LENGTH_VALUE="$EXPECTED_LENGTH" VERSION_VALUE="$VERSION_ID" \
"$PYTHON" - <<'PY'
import json,os
v=json.load(open(os.environ["HEAD_PATH"],encoding="utf-8"))
if v.get("ETag") != os.environ["EXPECTED_ETAG_VALUE"] or v.get("ContentLength") != int(os.environ["EXPECTED_LENGTH_VALUE"]) or v.get("VersionId") != os.environ["VERSION_VALUE"]:
    raise SystemExit("provider identity drift immediately before download")
PY
"$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
"$AWS" s3api get-object --no-sign-request --bucket "$BUCKET" --key "$KEY" --version-id "$VERSION_ID" "$SCRATCH" > "$STAGE/download_response.json"
[[ $("$PYTHON" -c 'import os,sys; print(os.stat(sys.argv[1]).st_size)' "$SCRATCH") -eq "$EXPECTED_LENGTH" ]] || {
    printf 'downloaded object length differs from provider receipt\n' >&2; exit 2
}
"$SAMTOOLS" quickcheck "$SCRATCH"
DOWNLOAD_SHA=$(shasum -a 256 "$SCRATCH" | awk '{print $1}')

"$SAMTOOLS" view -F 0x900 "$SCRATCH" | "$PYTHON" -c '
import sys
out=open(sys.argv[1],"x",encoding="ascii",newline="\n"); count=bases=0
for n,line in enumerate(sys.stdin,1):
    f=line.rstrip("\n").split("\t")
    if len(f)<11 or not f[0]: raise SystemExit(f"malformed primary source row {n}")
    flag=int(f[1]); cigar=f[5]; seq=f[9]
    if flag & 0x1 or seq=="*" or "H" in cigar: raise SystemExit(f"source row {n} is not one complete ONT molecule")
    out.write(f[0]+"\n"); count+=1; bases+=len(seq)
out.close()
if count < 1: raise SystemExit("source object contains no primary molecules")
' "$STAGE/source_qnames.unsorted"
LC_ALL=C sort "$STAGE/source_qnames.unsorted" > "$STAGE/source_primary_qnames.txt"
[[ $(wc -l < "$STAGE/source_qnames.unsorted") -eq $(uniq "$STAGE/source_primary_qnames.txt" | wc -l) ]] || {
    printf 'duplicate primary QNAME within source object\n' >&2; exit 2
}
DECODED_MOLECULES=$(wc -l < "$STAGE/source_primary_qnames.txt")
"$SAMTOOLS" fasta -F 0x900 "$SCRATCH" > "$STAGE/decoded_molecules.fa"
read -r DECODED_FASTA_COUNT DECODED_BASES < <("$PYTHON" -c '
import sys
src,out=sys.argv[1:]
names=[]; bases=0; current=None
with open(src,"rb") as handle:
    for number,raw in enumerate(handle,1):
        line=raw.rstrip(b"\r\n")
        if raw[-1:]!=b"\n": raise SystemExit("decoded FASTA is not LF terminated")
        if line.startswith(b">"):
            name=line[1:]
            if not name or any(b<33 or b>126 for b in name) or b" " in name or b"\t" in name: raise SystemExit(f"malformed decoded QNAME at line {number}")
            names.append(name); current=name
        else:
            if current is None or not line: raise SystemExit(f"malformed decoded sequence at line {number}")
            bases+=len(line)
with open(out,"xb") as target:
    for name in names: target.write(name+b"\n")
print(len(names),bases)
' "$STAGE/decoded_molecules.fa" "$STAGE/decoded_qnames.unsorted")
LC_ALL=C sort "$STAGE/decoded_qnames.unsorted" > "$STAGE/decoded_qnames.txt"
[[ "$DECODED_FASTA_COUNT" -eq "$DECODED_MOLECULES" && $(uniq "$STAGE/decoded_qnames.txt" | wc -l) -eq "$DECODED_MOLECULES" ]] || { printf 'decoded FASTA multiplicity differs from source inventory\n' >&2; exit 2; }
cmp -s "$STAGE/source_primary_qnames.txt" "$STAGE/decoded_qnames.txt" || { printf 'decoded FASTA QNAME set differs from source inventory\n' >&2; exit 2; }

"$MINIMAP2" -t 4 -ax map-ont --secondary=no "$COMBINED_MMI" "$STAGE/decoded_molecules.fa" \
    | "$SAMTOOLS" view -@ 1 -b -o "$STAGE/complete_unfiltered_alignment.bam" -
"$SAMTOOLS" quickcheck "$STAGE/complete_unfiltered_alignment.bam"
read -r PRIMARY_RECORDS SUPPLEMENTARY_RECORDS < <("$SAMTOOLS" view "$STAGE/complete_unfiltered_alignment.bam" | "$PYTHON" -c '
import sys
primary={}; supplementary=[]
for line in sys.stdin:
    f=line.split("\t"); q=f[0]; flag=int(f[1])
    if not q or any(ord(c)<33 or ord(c)>126 for c in q): raise SystemExit("malformed aligned QNAME")
    if flag & 0x100: raise SystemExit("secondary alignment present despite --secondary=no")
    if flag & 0x800: supplementary.append(q)
    else: primary[q]=primary.get(q,0)+1
if any(v!=1 for v in primary.values()): raise SystemExit("aligned molecule lacks exactly one primary")
with open(sys.argv[1],"x",encoding="ascii",newline="\n") as out:
    for q in primary: out.write(q+"\n")
with open(sys.argv[2],"x",encoding="ascii",newline="\n") as out:
    for q in supplementary: out.write(q+"\n")
print(len(primary),len(supplementary))
' "$STAGE/aligned_primary_qnames.unsorted" "$STAGE/supplementary_qnames.unsorted")
[[ "$PRIMARY_RECORDS" -eq "$DECODED_MOLECULES" ]] || { printf 'primary alignment inventory differs from decoded molecules\n' >&2; exit 2; }
LC_ALL=C sort "$STAGE/aligned_primary_qnames.unsorted" > "$STAGE/aligned_primary_qnames.txt"
[[ $(uniq "$STAGE/aligned_primary_qnames.txt" | wc -l) -eq "$PRIMARY_RECORDS" ]] || { printf 'aligned primary QNAME multiplicity is not one\n' >&2; exit 2; }
cmp -s "$STAGE/source_primary_qnames.txt" "$STAGE/aligned_primary_qnames.txt" || { printf 'aligned primary QNAME set differs from decoded/source set\n' >&2; exit 2; }
LC_ALL=C sort -u "$STAGE/supplementary_qnames.unsorted" > "$STAGE/supplementary_qnames.txt"
[[ -z "$(comm -23 "$STAGE/supplementary_qnames.txt" "$STAGE/source_primary_qnames.txt")" ]] || { printf 'supplementary alignment has unknown QNAME\n' >&2; exit 2; }
"$PYTHON" "$PLANNER" validate-qname-ownership \
    --source "$STAGE/source_primary_qnames.txt" \
    --decoded "$STAGE/decoded_qnames.txt" \
    --aligned-primary "$STAGE/aligned_primary_qnames.txt" \
    --supplementary "$STAGE/supplementary_qnames.txt" >/dev/null

"$SAMTOOLS" view -u -F 0x900 -L "$FULL_BED" "$STAGE/complete_unfiltered_alignment.bam" \
    | "$SAMTOOLS" sort -@ 4 -o "$STAGE/primary_full_window.bam" -
"$SAMTOOLS" index -@ 4 "$STAGE/primary_full_window.bam"
"$SAMTOOLS" quickcheck "$STAGE/primary_full_window.bam"
"$SAMTOOLS" view "$STAGE/primary_full_window.bam" | "$PYTHON" -c '
import sys
for n,line in enumerate(sys.stdin,1):
    if int(line.split("\t",2)[1]) & 0x900: raise SystemExit(f"excluded 0x900 record retained at row {n}")
'

PLAN_PATH="$PLAN" PART_INDEX_VALUE="$PART_INDEX" STAGE_PATH="$STAGE" DOWNLOAD_SHA_VALUE="$DOWNLOAD_SHA" DOWNLOAD_SIZE_VALUE="$EXPECTED_LENGTH" DECODED_MOLECULES_VALUE="$DECODED_MOLECULES" DECODED_BASES_VALUE="$DECODED_BASES" PRIMARY_RECORDS_VALUE="$PRIMARY_RECORDS" SUPPLEMENTARY_RECORDS_VALUE="$SUPPLEMENTARY_RECORDS" \
"$PYTHON" - <<'PY'
import hashlib,json,os,pathlib
p=json.load(open(os.environ["PLAN_PATH"],encoding="utf-8")); m=p["alignment_map"][int(os.environ["PART_INDEX_VALUE"])]; stage=pathlib.Path(os.environ["STAGE_PATH"])
pairs=[("complete_unfiltered_alignment.bam",m["unfiltered_bam"]),("primary_full_window.bam",m["primary_bam"]),("primary_full_window.bam.bai",m["primary_bai"]),("source_primary_qnames.txt",m["source_qname_inventory"]),("decoded_qnames.txt",m["decoded_qname_inventory"]),("aligned_primary_qnames.txt",m["aligned_primary_qname_inventory"]),("supplementary_qnames.txt",m["supplementary_qname_inventory"])]
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest(); arts=[{"path":dst,"size_bytes":(stage/src).stat().st_size,"sha256":sha(stage/src)} for src,dst in pairs]
d={"provider_object":m["provider_object"],"download_size_bytes":int(os.environ["DOWNLOAD_SIZE_VALUE"]),"download_sha256":os.environ["DOWNLOAD_SHA_VALUE"],"decoded_molecules":int(os.environ["DECODED_MOLECULES_VALUE"]),"decoded_bases":int(os.environ["DECODED_BASES_VALUE"]),"primary_records":int(os.environ["PRIMARY_RECORDS_VALUE"]),"supplementary_records":int(os.environ["SUPPLEMENTARY_RECORDS_VALUE"]),"qname_inventory_sha256":sha(stage/"aligned_primary_qnames.txt"),"source_qname_inventory_sha256":sha(stage/"source_primary_qnames.txt"),"decoded_qname_inventory_sha256":sha(stage/"decoded_qnames.txt"),"aligned_primary_qname_inventory_sha256":sha(stage/"aligned_primary_qnames.txt"),"supplementary_qname_inventory_sha256":sha(stage/"supplementary_qnames.txt"),"qname_set_cardinality":int(os.environ["DECODED_MOLECULES_VALUE"]),"exact_qname_set_and_multiplicity_equal":True,"supplementary_qnames_subset":True,"exclusion_mask_hex":"0x900","exclusion_mask_decimal":2304,"decode_argv":["samtools","fasta","-F","0x900"],"alignment_argv":["minimap2","-ax","map-ont","--secondary=no"],"primary_view_argv":["samtools","view","-u","-F","0x900","-L","full_window.bed"]}
r={"schema_version":"hml2_7p22_exact4_part_1","kind":"part","run_id":p["run_id"],"plan_digest":p["plan_digest"],"bundle_id":p["bundle_id"],"index":m["index"],"sample_id":m["sample_id"],"locus":m["locus"],"claim_token":"exact4-part-"+str(m["index"]),"artifacts":arts,"scientific_assertions":{"diagnostic_only":True,"complete_unfiltered_alignment_required":True,"competitive_combined_reference_required":True,"primary_only_depth_required":True,"read_cap_or_downsampling_allowed":False},"details":d,"read_inferred_copy_number":"not_estimated","receipt_digest":""}
enc=lambda v:(json.dumps(v,sort_keys=True,separators=(",",":"),ensure_ascii=False)+"\n").encode(); material=dict(r); material.pop("receipt_digest"); r["receipt_digest"]=hashlib.sha256(enc(material)).hexdigest(); (stage/"PART.json").write_bytes(enc(r))
PY

for pair in \
    "complete_unfiltered_alignment.bam|$FINAL_UNFILTERED" \
    "primary_full_window.bam|$FINAL_PRIMARY" \
    "primary_full_window.bam.bai|$FINAL_BAI" \
    "source_primary_qnames.txt|$FINAL_SOURCE_QNAMES" \
    "decoded_qnames.txt|$FINAL_DECODED_QNAMES" \
    "aligned_primary_qnames.txt|$FINAL_ALIGNED_QNAMES" \
    "supplementary_qnames.txt|$FINAL_SUPPLEMENTARY_QNAMES" \
    "PART.json|$FINAL_RECEIPT"
do
    source=${pair%%|*}; destination=${pair#*|}
    "$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
    "$PYTHON" "$PLANNER" publish --plan "$PLAN" --execution-release "$RELEASE" --staged "$STAGE/$source" --destination "$destination" >/dev/null
done
"$PYTHON" "$PLANNER" validate-part --plan "$PLAN" --part-index "$PART_INDEX" >/dev/null
