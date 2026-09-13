#!/usr/bin/env bash
# Generic, manifest-driven CNV-v2 controller.
#
# This entry point is deliberately a local submission driver, not a Slurm job.
# Rendering and validation are the default.  Nothing enters the scheduler unless
# --submit is present and the immutable run preparation has succeeded.

set -Eeuo pipefail
IFS=$'\n\t'
PYTHON=${HML2_PYTHON:-python3.11}

usage() {
    cat <<'EOF'
Usage:
  01_cnv_controller_v2.sh \
    --manifest MANIFEST.tsv \
    --expected-manifest-sha256 HEX \
    --run-root DIRECTORY \
    [--submit]

Default: validate the manifest and render the canonical run plan only.
--submit: prepare the run-scoped combined reference and submit the exact part
          array plus a receipt-conditioned analysis release controller.
EOF
}

die() {
    printf 'CNV-v2 controller: %s\n' "$*" >&2
    exit 2
}

manifest=""
expected_manifest_sha256=""
run_root=""
submit=false

while (($#)); do
    case "$1" in
        --manifest)
            (($# >= 2)) || die "--manifest requires a value"
            manifest=$2
            shift 2
            ;;
        --expected-manifest-sha256)
            (($# >= 2)) || die "--expected-manifest-sha256 requires a value"
            expected_manifest_sha256=$2
            shift 2
            ;;
        --run-root)
            (($# >= 2)) || die "--run-root requires a value"
            run_root=$2
            shift 2
            ;;
        --submit)
            submit=true
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

[[ -n "$manifest" ]] || die "--manifest is required"
[[ -n "$expected_manifest_sha256" ]] || die "--expected-manifest-sha256 is required"
[[ -n "$run_root" ]] || die "--run-root is required"
[[ -f "$manifest" ]] || die "manifest does not exist: $manifest"
[[ "$expected_manifest_sha256" =~ ^[0-9a-fA-F]{64}$ ]] || die "expected manifest SHA-256 must be 64 hexadecimal characters"

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
planner="${script_dir}/cnv_plan.py"
worker="${script_dir}/02_cnv_alignment_worker_v2.sh"
analysis="${script_dir}/03_cnv_analysis_v2.sh"
receipt_controller="${script_dir}/04_cnv_receipt_controller_v2.sh"

[[ -f "$planner" ]] || die "planner is missing: $planner"
[[ -f "$worker" ]] || die "worker is missing: $worker"
[[ -f "$analysis" ]] || die "analysis entry point is missing: $analysis"
[[ -f "$receipt_controller" ]] || die "receipt controller is missing: $receipt_controller"
command -v "$PYTHON" >/dev/null 2>&1 || die "CPython 3.11 runtime is required: $PYTHON"

# The planner performs strict manifest/header/part/assembly/hash validation and
# emits one absolute plan path per sample/locus target.  A multi-target manifest
# is therefore rendered in one pass without mixing output or resume identities.
render_output=$(mktemp "${TMPDIR:-/tmp}/cnv_v2.render.XXXXXX")
trap 'rm -f -- "$render_output"' EXIT INT TERM
"$PYTHON" "$planner" render \
    --manifest "$manifest" \
    --expected-manifest-sha256 "$expected_manifest_sha256" \
    --run-root "$run_root" \
    --print-plan-path >"$render_output"

plan_paths=()
while IFS= read -r plan_path || [[ -n "$plan_path" ]]; do
    [[ -n "$plan_path" ]] || continue
    [[ "$plan_path" = /* ]] || die "planner returned a non-absolute plan path: $plan_path"
    [[ -s "$plan_path" ]] || die "rendered plan is missing or empty: $plan_path"
    plan_paths[${#plan_paths[@]}]=$plan_path
done <"$render_output"
rm -f -- "$render_output"
trap - EXIT INT TERM
((${#plan_paths[@]} > 0)) || die "planner did not return any plan paths"

# Read only the small set of scheduler values needed by the controller.  The
# planner-generated plan is treated as data and validated before assignments are
# emitted with shell-safe quoting.
read_plan_info() {
    local selected_plan=$1
    eval "$("$PYTHON" - "$selected_plan" <<'PY'
import json, os, shlex, sys

path = os.path.realpath(sys.argv[1])
with open(path, encoding="utf-8") as handle:
    plan = json.load(handle)
if plan.get("schema_version") != "cnv_v2_plan_1":
    raise SystemExit("unsupported CNV-v2 plan schema")
parts = plan.get("parts")
if not isinstance(parts, list) or not parts:
    raise SystemExit("plan has no declared parts")
run_id = plan.get("run_id")
target_id = plan.get("target_id")
run_dir = plan.get("run_dir")
for label, value in (("run_id", run_id), ("target_id", target_id), ("run_dir", run_dir)):
    if not isinstance(value, str) or not value or "\n" in value:
        raise SystemExit(f"invalid {label} in plan")
if not os.path.isabs(run_dir):
    raise SystemExit("run_dir must be absolute")
print("RUN_ID=" + shlex.quote(run_id))
print("TARGET_ID=" + shlex.quote(target_id))
print("RUN_DIR=" + shlex.quote(run_dir))
print("PART_COUNT=" + str(len(parts)))
PY
)"
}

printf 'CNV-v2 rendered %s canonical target plan(s)\n' "${#plan_paths[@]}"
for plan_path in "${plan_paths[@]}"; do
    read_plan_info "$plan_path"
    printf '  target_id=%s run_id=%s parts=%s\n' "$TARGET_ID" "$RUN_ID" "$PART_COUNT"
    printf '    plan: %s\n' "$plan_path"
done

if [[ "$submit" != true ]]; then
    printf 'DRY RUN ONLY: no reference preparation, scheduler submission, download, alignment, or analysis was started.\n'
    printf 'Re-run the same command with --submit only after reviewing this canonical plan.\n'
    exit 0
fi

command -v sbatch >/dev/null 2>&1 || die "--submit requested but sbatch is not available"

# Preparation must create and validate the one run-scoped combined diploid
# reference/index, bait, clamped interval files, and an immutable PREPARED receipt.
# A stale or partial preparation is a hard failure; it is never silently rebuilt.
for plan_path in "${plan_paths[@]}"; do
    "$PYTHON" "$planner" prepare --plan "$plan_path"
done

for plan_path in "${plan_paths[@]}"; do
    read_plan_info "$plan_path"
    log_dir="${RUN_DIR}/slurm_logs"
    mkdir -p -- "$log_dir"
    array_max=$((PART_COUNT - 1))
    worker_job=$(
        sbatch --parsable \
            --job-name="cnv2_${TARGET_ID}_parts" \
            --array="0-${array_max}" \
            --output="${log_dir}/part_%A_%a.out" \
            --error="${log_dir}/part_%A_%a.err" \
            "$worker" --plan "$plan_path"
    )
    worker_job=${worker_job%%;*}
    [[ "$worker_job" =~ ^[0-9]+$ ]] || die "unexpected worker sbatch response: $worker_job"

    controller_job=$(
        sbatch --parsable \
            --job-name="cnv2_${TARGET_ID}_receipts" \
            --output="${log_dir}/receipt_controller_%j.out" \
            --error="${log_dir}/receipt_controller_%j.err" \
            "$receipt_controller" --plan "$plan_path" --analysis "$analysis"
    )
    controller_job=${controller_job%%;*}
    [[ "$controller_job" =~ ^[0-9]+$ ]] || die "unexpected receipt-controller sbatch response: $controller_job"

    printf 'Submitted exact CNV-v2 run %s for target %s\n' "$RUN_ID" "$TARGET_ID"
    printf '  part array job: %s (0-%s)\n' "$worker_job" "$array_max"
    printf '  receipt controller: %s (no Slurm dependency; analysis releases from DONE receipts)\n' "$controller_job"
done
printf 'Nucleotide-frequency analysis is not a coverage-completion condition and was not submitted.\n'
