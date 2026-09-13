#!/usr/bin/env bash
# Receipt-conditioned closure release.  This script is short-lived and
# resubmits itself without a Slurm parent dependency while unit receipts remain
# incomplete.
set -euo pipefail
export LC_ALL=C PYTHONDONTWRITEBYTECODE=1
umask 027

if [[ $# -ne 2 ]]; then
    printf 'usage: 03_cnv_array_receipt_controller.sh MANIFEST.json CLOSURE.sh\n' >&2
    exit 2
fi
MANIFEST=$1
CLOSURE=$2
PYTHON=${HML2_CLUSTER_PYTHON:-/path/to/hml2_workspace/condaenv/repeatmaskerenv/bin/python3.11}
[[ -x $PYTHON ]] || {
    printf 'selected HML2_CLUSTER_PYTHON is not executable: %s\n' "$PYTHON" >&2
    exit 2
}
PLANNER=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["execution_paths"]["planner"])' "$MANIFEST")
[[ $MANIFEST == /* && $PLANNER == /* && $CLOSURE == /* ]] || {
    printf 'manifest, planner and closure paths must be absolute\n' >&2
    exit 2
}

ready=0
for index in 0 1 2 3; do
    if "$PYTHON" "$PLANNER" validate-unit \
            --manifest "$MANIFEST" --unit-index "$index" >/dev/null 2>&1; then
        ready=$((ready + 1))
    fi
done

if [[ $ready -eq 4 ]]; then
    job=$(/usr/bin/sbatch --parsable --partition=batch --cpus-per-task=1 \
        --mem=1G --time=00:05:00 --job-name=cnv7p22_close \
        --export="HML2_CLUSTER_PYTHON=$PYTHON" "$CLOSURE" "$MANIFEST")
    printf 'CNV_EXACT4_CLOSURE_RELEASED_FROM_RECEIPTS job=%s\n' "${job%%;*}"
    exit 0
fi

next=$(/usr/bin/sbatch --parsable --partition=batch --cpus-per-task=1 \
    --mem=1G --time=00:05:00 --job-name=cnv7p22_receipts \
    --begin=now+5minutes --export="HML2_CLUSTER_PYTHON=$PYTHON" \
    "$0" "$MANIFEST" "$CLOSURE")
printf 'CNV_EXACT4_RECEIPTS_PENDING ready=%s expected=4 next=%s\n' \
    "$ready" "${next%%;*}"
