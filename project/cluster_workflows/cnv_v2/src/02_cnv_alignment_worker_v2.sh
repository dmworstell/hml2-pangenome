#!/usr/bin/env bash
#SBATCH --partition=batch
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G

# One immutable manifest part -> one target-bound combined-reference BAM.

set -Eeuo pipefail
IFS=$'\n\t'
PYTHON=${HML2_PYTHON:-python3.11}

usage() {
    cat <<'EOF'
Usage: 02_cnv_alignment_worker_v2.sh [--plan PLAN.json] [--part-index N] [--verify-only]

When submitted by the controller, PLAN comes from CNV_V2_PLAN and N comes from
SLURM_ARRAY_TASK_ID.  --verify-only validates an existing exact DONE receipt,
BAM, and index without downloading or changing anything.
EOF
}

die() {
    printf 'CNV-v2 alignment worker: %s\n' "$*" >&2
    exit 2
}

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum -- "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 -- "$1" | awk '{print $1}'
    else
        die "sha256sum or shasum is required"
    fi
}

plan=${CNV_V2_PLAN:-}
part_index=${SLURM_ARRAY_TASK_ID:-}
verify_only=false

while (($#)); do
    case "$1" in
        --plan)
            (($# >= 2)) || die "--plan requires a value"
            plan=$2
            shift 2
            ;;
        --part-index)
            (($# >= 2)) || die "--part-index requires a value"
            part_index=$2
            shift 2
            ;;
        --verify-only)
            verify_only=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown argument: $1"
            ;;
    esac
done

[[ -n "$plan" ]] || die "a plan is required (--plan or CNV_V2_PLAN)"
[[ -s "$plan" ]] || die "plan is missing or empty: $plan"
[[ "$part_index" =~ ^[0-9]+$ ]] || die "part index must be a zero-based integer"
command -v "$PYTHON" >/dev/null 2>&1 || die "CPython 3.11 runtime is required: $PYTHON"

plan=$("$PYTHON" - "$plan" <<'PY'
import os, sys
print(os.path.realpath(sys.argv[1]))
PY
)
plan_file_sha256=$(sha256_file "$plan")

# Extract a fixed, strictly validated record from the planner-generated JSON.
eval "$("$PYTHON" - "$plan" "$part_index" <<'PY'
import json, os, re, shlex, sys

path, index_text = sys.argv[1:]
with open(path, encoding="utf-8") as handle:
    plan = json.load(handle)
if plan.get("schema_version") != "cnv_v2_plan_1":
    raise SystemExit("unsupported CNV-v2 plan schema")
parts = plan.get("parts")
if not isinstance(parts, list) or not parts:
    raise SystemExit("plan has no declared parts")
index = int(index_text)
if index < 0 or index >= len(parts):
    raise SystemExit("part index is outside the declared plan")
part = parts[index]
manifest = plan.get("manifest")
reference = plan.get("reference")
if not isinstance(manifest, dict) or not isinstance(reference, dict):
    raise SystemExit("plan is missing manifest/reference records")

values = {
    "RUN_ID": plan.get("run_id"),
    "TARGET_ID": plan.get("target_id"),
    "TARGET_DIGEST": plan.get("target_digest"),
    "MANIFEST_SHA256": manifest.get("sha256"),
    "PLAN_DIGEST": plan.get("plan_digest"),
    "PREPARED_RECEIPT": plan.get("prepared_receipt"),
    "COMBINED_FAI": reference.get("combined_fai"),
    "COMBINED_MMI": reference.get("combined_mmi"),
    "BAIT_FASTA": reference.get("bait_fasta"),
    "FULL_WINDOW_BED": reference.get("full_window_bed"),
    "READ_URI": part.get("read_uri"),
    "READ_VERSION": part.get("read_version"),
    "READ_ETAG": part.get("read_etag"),
    "READ_SHA256": part.get("read_sha256"),
    "READ_FORMAT": part.get("read_format"),
    "PRESET": part.get("minimap2_preset"),
    "PART_BAM": part.get("bam"),
    "PART_BAI": part.get("bai"),
    "PART_PRIMARY_QNAMES": part.get("primary_qnames"),
    "PART_CLAIM": part.get("claim"),
    "PART_RECEIPT": part.get("receipt"),
    "PLANNER": plan.get("code", {}).get("cnv_plan.py", {}).get("path"),
    "TARGET_DIR": plan.get("target_dir", plan.get("run_dir")),
}
for name, value in values.items():
    if not isinstance(value, str) or not value or "\n" in value or "\x00" in value:
        raise SystemExit(f"invalid or missing {name} in plan")
for name in ("MANIFEST_SHA256", "PLAN_DIGEST", "TARGET_DIGEST", "READ_SHA256"):
    if not re.fullmatch(r"[0-9a-f]{64}", values[name]):
        raise SystemExit(f"invalid {name} in plan")
read_size = part.get("read_size")
if not isinstance(read_size, int) or read_size <= 0:
    raise SystemExit("read_size must be a positive integer")
part_id = part.get("part_id")
if not isinstance(part_id, int) or part_id != index + 1:
    raise SystemExit("part_id must be the exact positive plan index plus one")
bait_min = plan.get("bait_min_aligned_bp")
if not isinstance(bait_min, int) or bait_min <= 0:
    raise SystemExit("bait_min_aligned_bp must be a positive integer")
if values["READ_FORMAT"] != "bam":
    raise SystemExit("CNV-v2 currently requires BAM input parts")
if values["PRESET"] not in {"map-ont", "map-hifi"}:
    raise SystemExit("unsupported minimap2 preset")
for name in ("PREPARED_RECEIPT", "COMBINED_FAI", "COMBINED_MMI", "BAIT_FASTA", "FULL_WINDOW_BED", "PART_BAM", "PART_BAI", "PART_PRIMARY_QNAMES", "PART_CLAIM", "PART_RECEIPT", "PLANNER", "TARGET_DIR"):
    if not os.path.isabs(values[name]):
        raise SystemExit(f"{name} must be absolute")
if os.path.commonpath([os.path.realpath(values["TARGET_DIR"]), os.path.realpath(values["PART_BAM"])]) != os.path.realpath(values["TARGET_DIR"]):
    raise SystemExit("part BAM is outside the canonical target directory")
output_dirs = {os.path.dirname(os.path.realpath(values[name])) for name in ("PART_BAM", "PART_BAI", "PART_PRIMARY_QNAMES", "PART_CLAIM", "PART_RECEIPT")}
if len(output_dirs) != 1:
    raise SystemExit("part BAM, BAI, QNAME inventory, claim, and DONE receipt must share one canonical directory")
if part.get("part_index") not in (None, index):
    raise SystemExit("part_index does not match the array index")
for name, value in values.items():
    print(name + "=" + shlex.quote(value))
print("READ_SIZE=" + str(read_size))
print("PART_ID=" + str(part_id))
print("BAIT_MIN_ALIGNED_BP=" + str(bait_min))
print("DECLARED_PART_INDEX=" + str(index))
PY
)"

[[ "$PART_BAI" == "${PART_BAM}.bai" ]] || die "plan must bind the BAI to BAM.bai"

verify_prepared_receipt() {
    "$PYTHON" - "$plan" "$PREPARED_RECEIPT" <<'PY'
import hashlib, json, os, sys

plan_path, receipt_path = sys.argv[1:]
with open(plan_path, encoding="utf-8") as handle:
    plan = json.load(handle)
with open(receipt_path, encoding="utf-8") as handle:
    receipt = json.load(handle)
claimed_digest = receipt.get("receipt_digest")
material = dict(receipt)
material.pop("receipt_digest", None)
canonical = (json.dumps(material, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode()
if not isinstance(claimed_digest, str) or hashlib.sha256(canonical).hexdigest() != claimed_digest:
    raise SystemExit("PREPARED receipt digest mismatch")
expected = {
    "schema_version": "cnv_v2_prepared_1",
    "status": "PREPARED",
    "run_id": plan["run_id"],
    "target_id": plan["target_id"],
    "target_digest": plan["target_digest"],
    "manifest_sha256": plan["manifest"]["sha256"],
    "plan_digest": plan["plan_digest"],
}
for key, value in expected.items():
    if receipt.get(key) != value:
        raise SystemExit(f"PREPARED receipt mismatch: {key}")
artifacts = receipt.get("artifacts")
if not isinstance(artifacts, dict):
    raise SystemExit("PREPARED receipt has no artifacts map")
required = {
    "combined_fasta": plan["reference"]["combined_fasta"],
    "combined_fai": plan["reference"]["combined_fai"],
    "combined_mmi": plan["reference"]["combined_mmi"],
    "bait_fasta": plan["reference"]["bait_fasta"],
    "full_window_bed": plan["reference"]["full_window_bed"],
    "core_body_union_bed": plan["reference"]["core_body_union_bed"],
    "outer_flanks_bed": plan["reference"]["outer_flanks_bed"],
}
if set(artifacts) != set(required):
    raise SystemExit("PREPARED receipt artifact set is not the exact seven-artifact closure")
for label, path in required.items():
    record = artifacts.get(label)
    if not os.path.isfile(path):
        raise SystemExit(f"PREPARED artifact missing: {label}")
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    observed = {"path": path, "size": os.path.getsize(path), "sha256": digest.hexdigest()}
    if record != observed:
        raise SystemExit(f"PREPARED artifact digest mismatch: {label}")
PY
}

verify_part_receipt() {
    "$PYTHON" "$PLANNER" validate-part --plan "$plan" --part-id "$PART_ID" --samtools samtools >/dev/null
}

verify_prepared_receipt
command -v samtools >/dev/null 2>&1 || die "samtools is required"

state_count=0
[[ -e "$PART_BAM" || -L "$PART_BAM" ]] && state_count=$((state_count + 1))
[[ -e "$PART_BAI" || -L "$PART_BAI" ]] && state_count=$((state_count + 1))
[[ -e "$PART_PRIMARY_QNAMES" || -L "$PART_PRIMARY_QNAMES" ]] && state_count=$((state_count + 1))
[[ -e "$PART_RECEIPT" || -L "$PART_RECEIPT" ]] && state_count=$((state_count + 1))
if ((state_count > 0)); then
    ((state_count == 4)) || die "partial prior state exists; refusing to overwrite or rebuild it"
    verify_part_receipt || die "existing part state is stale or invalid; refusing reuse"
    printf 'Exact DONE receipt validated for run=%s target=%s part=%s\n' "$RUN_ID" "$TARGET_ID" "$PART_ID"
    exit 0
fi
[[ "$verify_only" != true ]] || die "no completed part exists to verify"

for tool in minimap2 sort awk; do
    command -v "$tool" >/dev/null 2>&1 || die "$tool is required"
done

part_output_dir=$(dirname -- "$PART_BAM")
mkdir -p -- "$part_output_dir"
claimant_token="worker:${SLURM_JOB_ID:-local}:${SLURM_ARRAY_JOB_ID:-none}:${SLURM_ARRAY_TASK_ID:-$part_index}:${HOSTNAME:-unknown}:$$"
CLAIM_DIGEST=$("$PYTHON" "$PLANNER" claim --plan "$plan" --kind part --part-id "$PART_ID" --claimant-token "$claimant_token") \
    || die "exact part claim already exists, is stale, or could not be acquired"
[[ "$CLAIM_DIGEST" =~ ^[0-9a-f]{64}$ ]] || die "part claim did not return an exact digest"
ncpus=${SLURM_CPUS_PER_TASK:-1}
[[ "$ncpus" =~ ^[1-9][0-9]*$ ]] || die "SLURM_CPUS_PER_TASK must be a positive integer"
decode_threads=$((ncpus >= 4 ? 2 : 1))
scan_threads=$((ncpus - decode_threads))
((scan_threads >= 1)) || scan_threads=1
realign_threads=$((ncpus >= 4 ? ncpus - 2 : 1))
sort_threads=$((ncpus - realign_threads))
((sort_threads >= 1)) || sort_threads=1

scratch=$(mktemp -d "${part_output_dir}/.${TARGET_ID}.part_${PART_ID}.work.XXXXXX")
cleanup() {
    rc=$?
    rm -rf -- "$scratch"
    exit "$rc"
}
trap cleanup EXIT INT TERM

local_read="${scratch}/input.bam"
if [[ "$READ_URI" == file://* ]]; then
    eval "$("$PYTHON" - "$READ_URI" <<'PY'
import os, shlex, sys, urllib.parse
parsed = urllib.parse.urlsplit(sys.argv[1])
if parsed.scheme != "file" or parsed.netloc not in {"", "localhost"}:
    raise SystemExit("file URI must be local")
path = urllib.parse.unquote(parsed.path)
if not os.path.isabs(path):
    raise SystemExit("file URI must contain an absolute path")
print("local_source=" + shlex.quote(path))
PY
)"
    [[ -f "$local_source" ]] || die "local read object is missing: $local_source"
    cp -- "$local_source" "$local_read"
elif [[ "$READ_URI" == s3://* ]]; then
    command -v aws >/dev/null 2>&1 || die "aws CLI is required for s3 read objects"
    [[ -n "$READ_VERSION" && -n "$READ_ETAG" ]] || die "s3 objects require expected immutable version and ETag metadata"
    eval "$("$PYTHON" - "$READ_URI" <<'PY'
import shlex, sys, urllib.parse
parsed = urllib.parse.urlsplit(sys.argv[1])
if parsed.scheme != "s3" or not parsed.netloc or not parsed.path.lstrip("/"):
    raise SystemExit("invalid s3 URI")
print("S3_BUCKET=" + shlex.quote(parsed.netloc))
print("S3_KEY=" + shlex.quote(urllib.parse.unquote(parsed.path.lstrip("/"))))
PY
)"
    head_json="${scratch}/head.json"
    aws s3api head-object \
        --bucket "$S3_BUCKET" \
        --key "$S3_KEY" \
        --version-id "$READ_VERSION" >"$head_json"
    "$PYTHON" - "$head_json" "$READ_VERSION" "$READ_ETAG" "$READ_SIZE" <<'PY'
import json, sys
path, version, etag, size = sys.argv[1:]
with open(path, encoding="utf-8") as handle:
    head = json.load(handle)
if head.get("VersionId") != version:
    raise SystemExit("S3 VersionId differs from the manifest")
if str(head.get("ETag", "")).strip('"') != etag.strip('"'):
    raise SystemExit("S3 ETag differs from the manifest")
if head.get("ContentLength") != int(size):
    raise SystemExit("S3 object size differs from the manifest")
PY
    aws s3api get-object \
        --bucket "$S3_BUCKET" \
        --key "$S3_KEY" \
        --version-id "$READ_VERSION" \
        "$local_read" >/dev/null
else
    die "read URI must use s3:// or file://"
fi

actual_size=$("$PYTHON" - "$local_read" <<'PY'
import os, sys
print(os.path.getsize(sys.argv[1]))
PY
)
[[ "$actual_size" == "$READ_SIZE" ]] || die "downloaded object size differs from the manifest"
actual_sha256=$(sha256_file "$local_read")
[[ "$actual_sha256" == "$READ_SHA256" ]] || die "downloaded object SHA-256 differs from the manifest"
samtools quickcheck -v "$local_read"

# A source BAM may carry secondary/supplementary records for a molecule.  They
# are not independent reads and must not multiply its weight.  Restrict both
# decode passes to primary records, and fail closed if the source nevertheless
# contains more than one primary record with the same query name.
primary_names="${scratch}/primary.names"
duplicate_primary_names="${scratch}/duplicate_primary.names"
samtools view -F 0x900 "$local_read" | awk '{print $1}' | LC_ALL=C sort >"$primary_names"
uniq -d "$primary_names" >"$duplicate_primary_names"
[[ ! -s "$duplicate_primary_names" ]] || die "source BAM contains multiple primary records for at least one query name"
rm -f -- "$primary_names" "$duplicate_primary_names"

[[ -s "$BAIT_FASTA" ]] || die "prepared bait FASTA is missing or empty"
[[ -s "$COMBINED_MMI" ]] || die "prepared combined-reference minimap2 index is missing or empty"
[[ -s "$FULL_WINDOW_BED" ]] || die "prepared full analysis-window BED is missing or empty"

names="${scratch}/captured.names"
captured="${scratch}/captured.fastq"

# Scan every input read against the bait.  The summed PAF block threshold is
# fixed by the validated manifest.  There is no read-count ceiling, sampling,
# truncation, or fallback to another computation.
samtools view -u -F 0x900 "$local_read" \
    | samtools fastq -@ "$decode_threads" - \
    | minimap2 -x "$PRESET" -t "$scan_threads" "$BAIT_FASTA" - \
    | awk -v minimum="$BAIT_MIN_ALIGNED_BP" '{total[$1]+=$11} END {for (name in total) if (total[name] >= minimum) print name}' \
    | sort -u >"$names"

if [[ -s "$names" ]]; then
    samtools view -N "$names" -u -F 0x900 "$local_read" \
        | samtools fastq -@ "$decode_threads" - >"$captured"
else
    : >"$captured"
fi

staged_bam="${scratch}/part.sorted.bam"

# Each captured molecule enters exactly one base-level alignment, against the
# single combined two-assembly reference.  Filtering uses the complete clamped
# analysis window, not the insertion body interval.
minimap2 -ax "$PRESET" --secondary=no -t "$realign_threads" "$COMBINED_MMI" "$captured" \
    | samtools view -u -F 0x900 -L "$FULL_WINDOW_BED" - \
    | samtools sort -@ "$sort_threads" -o "$staged_bam" -
samtools index "$staged_bam"
samtools quickcheck -v "$staged_bam"
samtools idxstats "$staged_bam" >/dev/null
[[ -s "$staged_bam" && -s "${staged_bam}.bai" ]] || die "alignment did not produce a BAM and index"

# The final primary-only BAM owns the exact canonical molecule inventory used
# for cross-part closure.  Bytewise sorting is deterministic and has no cap.
staged_qnames="${scratch}/primary_qnames.txt"
samtools view -F 0x900 "$staged_bam" | awk '{print $1}' | LC_ALL=C sort -u >"$staged_qnames"
"$PYTHON" "$PLANNER" verify-qname-union --inventory "$staged_qnames" >/dev/null \
    || die "final BAM did not produce a canonical primary-QNAME inventory"

# Canonical outputs are create-only publications from the destination
# filesystem.  A crash in between intentionally leaves detectable partial
# state; no repeat can overwrite it.  DONE remains last.
"$PYTHON" "$PLANNER" publish --staged "$staged_bam" --destination "$PART_BAM" >/dev/null
"$PYTHON" "$PLANNER" publish --staged "${staged_bam}.bai" --destination "$PART_BAI" >/dev/null
"$PYTHON" "$PLANNER" publish --staged "$staged_qnames" --destination "$PART_PRIMARY_QNAMES" >/dev/null

receipt_stage="${scratch}/DONE.json"
"$PYTHON" - "$plan" "$DECLARED_PART_INDEX" "$PART_BAM" "$PART_BAI" "$PART_PRIMARY_QNAMES" "$plan_file_sha256" "$PREPARED_RECEIPT" "$CLAIM_DIGEST" "$receipt_stage" <<'PY'
import hashlib, json, os, sys

plan_path, index_text, bam, bai, qnames, plan_file_sha, prepared_path, claim_digest, output = sys.argv[1:]
with open(plan_path, encoding="utf-8") as handle:
    plan = json.load(handle)
index = int(index_text)
part = plan["parts"][index]
with open(prepared_path, encoding="utf-8") as handle:
    prepared = json.load(handle)

def artifact(path, allow_empty=False):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    size=os.path.getsize(path)
    if not allow_empty and size <= 0: raise SystemExit(f"empty output artifact: {path}")
    return {"path": os.path.realpath(path), "size": size, "sha256": digest.hexdigest()}

def qname_artifact(path):
    previous=None; count=0
    with open(path,"rb") as handle:
        for number,raw in enumerate(handle,1):
            if not raw.endswith(b"\n"): raise SystemExit(f"QNAME inventory lacks terminal LF at line {number}")
            value=raw[:-1]
            if not value or any(byte<=32 or byte>=127 for byte in value): raise SystemExit("invalid QNAME inventory byte")
            if previous is not None and value<=previous: raise SystemExit("QNAME inventory is not byte-sorted unique")
            previous=value; count+=1
    return {**artifact(path,allow_empty=True),"qname_count":count,
            "encoding":"SAM-QNAME ASCII bytes, LF terminated",
            "order":"LC_ALL=C bytewise ascending, unique"}

receipt = {
    "schema_version": "cnv_v2_part_done_1",
    "status": "DONE",
    "run_id": plan["run_id"],
    "target_id": plan["target_id"],
    "target_digest": plan["target_digest"],
    "manifest_sha256": plan["manifest"]["sha256"],
    "plan_digest": plan["plan_digest"],
    "plan_file_sha256": plan_file_sha,
    "part_id": part["part_id"],
    "part_index": index,
    "read_uri": part["read_uri"],
    "read_version": part["read_version"],
    "read_etag": part["read_etag"],
    "read_size": part["read_size"],
    "read_sha256": part["read_sha256"],
    "prepared_receipt_digest": prepared["receipt_digest"],
    "claim_digest": claim_digest,
    "prepared_artifacts": prepared["artifacts"],
    "outputs": {"bam": artifact(bam), "bai": artifact(bai), "primary_qnames": qname_artifact(qnames)},
}
canonical = (json.dumps(receipt, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode()
receipt["receipt_digest"] = hashlib.sha256(canonical).hexdigest()
with open(output, "x", encoding="utf-8") as handle:
    json.dump(receipt, handle, sort_keys=True, indent=2)
    handle.write("\n")
PY
"$PYTHON" "$PLANNER" publish --staged "$receipt_stage" --destination "$PART_RECEIPT" >/dev/null

verify_part_receipt
printf 'CNV-v2 part complete: run=%s target=%s part=%s\n' "$RUN_ID" "$TARGET_ID" "$PART_ID"
