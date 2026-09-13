#!/usr/bin/env bash
# Exact-four primary-only depth and diagnostic analysis entry point.
set -euo pipefail
export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1
umask 027

PLAN=""; RELEASE=""; SAMPLE_INDEX=""; VERIFY_ONLY=false
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
PLOTTER="$HERE/plot_exact4_depth.py"
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
    "$PYTHON" "$PLANNER" validate-analysis --plan "$PLAN" --sample-index "$SAMPLE_INDEX" >/dev/null
    exit 0
fi

mapfile -t PART_INDICES < <("$PYTHON" -c 'import json,sys; p=json.load(open(sys.argv[1])); print(*p["analysis_map"][int(sys.argv[2])]["expected_global_part_indices"],sep="\n")' "$PLAN" "$SAMPLE_INDEX")
for index in "${PART_INDICES[@]}"; do
    "$PYTHON" "$PLANNER" validate-part --plan "$PLAN" --part-index "$index" >/dev/null
done
"$PYTHON" "$PLANNER" claim --plan "$PLAN" --execution-release "$RELEASE" \
    --kind analysis --index "$SAMPLE_INDEX" --token "exact4-analysis-${SAMPLE_INDEX}" >/dev/null

mapfile -t FIELDS < <("$PYTHON" - "$PLAN" "$SAMPLE_INDEX" <<'PY'
import json,sys
p=json.load(open(sys.argv[1],encoding="utf-8")); m=p["analysis_map"][int(sys.argv[2])]
tools={x["name"]:x["path"] for x in p["immutable_run_material"]["resident_identity"]["tools"]}
kcon=next(x["path"] for x in p["immutable_run_material"]["resident_identity"]["resident_files"] if x["role"]=="type2_kcon")
prep=p["preparation_map"][m["index"]]
for v in (m["staging_root"]+"/work",tools["samtools"],tools["minimap2"],kcon,prep["destinations"][4],prep["destinations"][5],m["merged_primary_bam"],m["merged_primary_bai"],m["raw_depth"],m["mapq10_depth"],m["summary"],m["kcon_paf"],m["raw_pdf"],m["mapq10_pdf"],m["receipt"],m["done"]): print(v)
PY
)
STAGE=${FIELDS[0]}; SAMTOOLS=${FIELDS[1]}; MINIMAP2=${FIELDS[2]}; KCON=${FIELDS[3]}
FULL_BED=${FIELDS[4]}; BAIT=${FIELDS[5]}; FINAL_MERGED=${FIELDS[6]}; FINAL_BAI=${FIELDS[7]}
FINAL_RAW=${FIELDS[8]}; FINAL_MAPQ=${FIELDS[9]}; FINAL_SUMMARY=${FIELDS[10]}; FINAL_PAF=${FIELDS[11]}
FINAL_RAW_PDF=${FIELDS[12]}; FINAL_MAPQ_PDF=${FIELDS[13]}; FINAL_RECEIPT=${FIELDS[14]}; FINAL_DONE=${FIELDS[15]}

"$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
[[ ! -e "$STAGE" && ! -L "$STAGE" ]] || { printf 'staging path already exists: %s\n' "$STAGE" >&2; exit 2; }
mkdir -m 0750 -- "$STAGE"

mapfile -t PRIMARY_BAMS < <("$PYTHON" -c 'import json,sys; p=json.load(open(sys.argv[1])); m=p["analysis_map"][int(sys.argv[2])]; print(*(p["alignment_map"][i]["primary_bam"] for i in m["expected_global_part_indices"]),sep="\n")' "$PLAN" "$SAMPLE_INDEX")
"$SAMTOOLS" merge -@ 4 -o "$STAGE/merged_primary_full_window.bam" "${PRIMARY_BAMS[@]}"
"$SAMTOOLS" index -@ 4 "$STAGE/merged_primary_full_window.bam"
"$SAMTOOLS" quickcheck "$STAGE/merged_primary_full_window.bam"
"$SAMTOOLS" view "$STAGE/merged_primary_full_window.bam" | "$PYTHON" -c '
import sys
for n,line in enumerate(sys.stdin,1):
    if int(line.split("\t",2)[1]) & 0x900: raise SystemExit(f"merged depth BAM contains excluded record at {n}")
'
"$SAMTOOLS" depth -a -g 0x600 -b "$FULL_BED" "$STAGE/merged_primary_full_window.bam" > "$STAGE/raw.depth.tsv"
"$SAMTOOLS" depth -a -g 0x600 -Q 10 -b "$FULL_BED" "$STAGE/merged_primary_full_window.bam" > "$STAGE/mapq10.depth.tsv"
"$MINIMAP2" -c -p 0.1 -N 10 "$BAIT" "$KCON" > "$STAGE/type2_KCON.paf"

PYTHONPATH="$HERE" PLAN_PATH="$PLAN" SAMPLE_INDEX_VALUE="$SAMPLE_INDEX" STAGE_PATH="$STAGE" \
"$PYTHON" - <<'PY'
import json,os,pathlib
import exact4_plan,plot_exact4_depth
p=exact4_plan.validate_target_only_plan(json.load(open(os.environ["PLAN_PATH"],encoding="utf-8")))
s=pathlib.Path(os.environ["STAGE_PATH"])
summary=plot_exact4_depth._build_summary_from_plan(p,int(os.environ["SAMPLE_INDEX_VALUE"]),s/"raw.depth.tsv",s/"mapq10.depth.tsv",s/"type2_KCON.paf")
(s/"depth_summary.json").write_bytes(exact4_plan.canonical_json_bytes(summary))
PY
"$PYTHON" "$PLOTTER" --depth "$STAGE/raw.depth.tsv" --summary "$STAGE/depth_summary.json" --output "$STAGE/raw.depth.pdf" --metric raw --target-only-plan "$PLAN" --sample-index "$SAMPLE_INDEX" --kcon-paf "$STAGE/type2_KCON.paf"
"$PYTHON" "$PLOTTER" --depth "$STAGE/mapq10.depth.tsv" --summary "$STAGE/depth_summary.json" --output "$STAGE/mapq10.depth.pdf" --metric mapq10 --target-only-plan "$PLAN" --sample-index "$SAMPLE_INDEX" --kcon-paf "$STAGE/type2_KCON.paf"

PLAN_PATH="$PLAN" SAMPLE_INDEX_VALUE="$SAMPLE_INDEX" STAGE_PATH="$STAGE" PLANNER_DIR="$HERE" \
"$PYTHON" - <<'PY'
import hashlib,json,os,pathlib,sys
sys.path.insert(0,os.environ["PLANNER_DIR"]); import exact4_plan,plot_exact4_depth
p=json.load(open(os.environ["PLAN_PATH"],encoding="utf-8")); m=p["analysis_map"][int(os.environ["SAMPLE_INDEX_VALUE"])]; stage=pathlib.Path(os.environ["STAGE_PATH"])
qpaths=[p["alignment_map"][i]["qname_inventory"] for i in m["expected_global_part_indices"]]
qunion=exact4_plan.verify_qname_inventory_union(qpaths)
pairs=[("merged_primary_full_window.bam",m["merged_primary_bam"]),("merged_primary_full_window.bam.bai",m["merged_primary_bai"]),("raw.depth.tsv",m["raw_depth"]),("mapq10.depth.tsv",m["mapq10_depth"]),("depth_summary.json",m["summary"]),("type2_KCON.paf",m["kcon_paf"]),("raw.depth.pdf",m["raw_pdf"]),("mapq10.depth.pdf",m["mapq10_pdf"])]
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest(); arts=[{"path":dst,"size_bytes":(stage/src).stat().st_size,"sha256":sha(stage/src)} for src,dst in pairs]
d={"expected_global_part_indices":m["expected_global_part_indices"],"qname_union":qunion,"merged_primary_bam_sha256":sha(stage/"merged_primary_full_window.bam"),"raw_depth_sha256":sha(stage/"raw.depth.tsv"),"mapq10_depth_sha256":sha(stage/"mapq10.depth.tsv"),"summary_sha256":sha(stage/"depth_summary.json"),"kcon_paf_sha256":sha(stage/"type2_KCON.paf"),"raw_pdf_sha256":sha(stage/"raw.depth.pdf"),"mapq10_pdf_sha256":sha(stage/"mapq10.depth.pdf"),"raw_depth_argv":["samtools","depth","-a","-g","0x600","-b","full_window.bed","merged_primary_full_window.bam"],"mapq10_depth_argv":["samtools","depth","-a","-g","0x600","-Q","10","-b","full_window.bed","merged_primary_full_window.bam"],"kcon_argv":["minimap2","-c","-p","0.1","-N","10","full-window.fa","type2_KCON.fa"],"kcon_selected_annotations":plot_exact4_depth._select_kcon_annotations_from_plan(p,int(os.environ["SAMPLE_INDEX_VALUE"]),stage/"type2_KCON.paf"),"deployment_manifest_sha256":p["deployment_manifest_sha256"],"entrypoint_hashes":p["entrypoint_hashes"],"coverage_ratio_schema_identity":exact4_plan.COVERAGE_RATIO_SCHEMA_IDENTITY,"coverage_ratio_schema_sha256":exact4_plan.COVERAGE_RATIO_SCHEMA_SHA256,"done_written_last":True}
d.update(exact4_plan.CONTRACT_BINDINGS)
r={"schema_version":"hml2_7p22_exact4_analysis_1","kind":"analysis","run_id":p["run_id"],"plan_digest":p["plan_digest"],"bundle_id":p["bundle_id"],"index":m["index"],"sample_id":m["sample_id"],"locus":m["locus"],"claim_token":"exact4-analysis-"+str(m["index"]),"artifacts":arts,"scientific_assertions":{"diagnostic_only":True,"complete_unfiltered_alignment_required":True,"competitive_combined_reference_required":True,"primary_only_depth_required":True,"read_cap_or_downsampling_allowed":False},"details":d,"read_inferred_copy_number":"not_estimated","receipt_digest":""}
enc=exact4_plan.canonical_json_bytes; material=dict(r); material.pop("receipt_digest"); r["receipt_digest"]=hashlib.sha256(enc(material)).hexdigest(); (stage/"ANALYSIS.json").write_bytes(enc(r))
PY

for pair in \
    "merged_primary_full_window.bam|$FINAL_MERGED" \
    "merged_primary_full_window.bam.bai|$FINAL_BAI" \
    "raw.depth.tsv|$FINAL_RAW" \
    "mapq10.depth.tsv|$FINAL_MAPQ" \
    "depth_summary.json|$FINAL_SUMMARY" \
    "type2_KCON.paf|$FINAL_PAF" \
    "raw.depth.pdf|$FINAL_RAW_PDF" \
    "mapq10.depth.pdf|$FINAL_MAPQ_PDF" \
    "ANALYSIS.json|$FINAL_RECEIPT"
do
    source=${pair%%|*}; destination=${pair#*|}
    "$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
    "$PYTHON" "$PLANNER" publish --plan "$PLAN" --execution-release "$RELEASE" --staged "$STAGE/$source" --destination "$destination" >/dev/null
done

"$PYTHON" "$PLANNER" validate-release --release "$RELEASE" --plan "$PLAN" >/dev/null
PLAN_PATH="$PLAN" SAMPLE_INDEX_VALUE="$SAMPLE_INDEX" RECEIPT_PATH="$FINAL_RECEIPT" DONE_STAGE="$STAGE/DONE.json" PLANNER_DIR="$HERE" \
"$PYTHON" - <<'PY'
import hashlib,json,os,sys
sys.path.insert(0,os.environ["PLANNER_DIR"]); import exact4_plan
p=json.load(open(os.environ["PLAN_PATH"],encoding="utf-8")); m=p["analysis_map"][int(os.environ["SAMPLE_INDEX_VALUE"])]
d={"schema_version":"hml2_7p22_exact4_done_1","run_id":p["run_id"],"plan_digest":p["plan_digest"],"bundle_id":p["bundle_id"],"deployment_manifest_sha256":p["deployment_manifest_sha256"],"entrypoint_hashes":p["entrypoint_hashes"],"sample_id":m["sample_id"],"analysis_receipt_sha256":exact4_plan.sha256_file(os.environ["RECEIPT_PATH"]),"coverage_ratio_schema_identity":exact4_plan.COVERAGE_RATIO_SCHEMA_IDENTITY,"coverage_ratio_schema_sha256":exact4_plan.COVERAGE_RATIO_SCHEMA_SHA256,"read_inferred_copy_number":"not_estimated","done_digest":""}
d.update(exact4_plan.CONTRACT_BINDINGS)
material=dict(d); material.pop("done_digest"); d["done_digest"]=hashlib.sha256(exact4_plan.canonical_json_bytes(material)).hexdigest(); open(os.environ["DONE_STAGE"],"xb").write(exact4_plan.canonical_json_bytes(d))
PY
"$PYTHON" "$PLANNER" publish --plan "$PLAN" --execution-release "$RELEASE" --staged "$STAGE/DONE.json" --destination "$FINAL_DONE" >/dev/null
"$PYTHON" "$PLANNER" validate-analysis --plan "$PLAN" --sample-index "$SAMPLE_INDEX" >/dev/null
