#!/usr/bin/env bash
# Create-only closure; release only after the controller validates all four unit
# checkpoint receipts.
set -euo pipefail
export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1
export APPTAINER_BIND=/cluster
umask 027

DRY_RUN=false
if [[ ${1:-} == --dry-run ]]; then
    DRY_RUN=true
    shift
fi
if [[ $# -ne 1 ]]; then
    printf 'usage: 02_cnv_array_closure.sh [--dry-run] MANIFEST.json\n' >&2
    exit 2
fi
MANIFEST=$1
PYTHON=${HML2_CLUSTER_PYTHON:-/path/to/hml2_workspace/condaenv/repeatmaskerenv/bin/python3.11}
if [[ $DRY_RUN == true && ! -x $PYTHON ]]; then
    PYTHON=$(command -v python3)
fi
[[ -x $PYTHON ]] || {
    printf 'selected HML2_CLUSTER_PYTHON is not executable: %s\n' "$PYTHON" >&2
    exit 2
}
PLANNER=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["execution_paths"]["planner"])' "$MANIFEST")
[[ $MANIFEST == /* && $PLANNER == /* ]] || {
    printf 'manifest and manifest-bound planner paths must be absolute\n' >&2
    exit 2
}
if [[ $DRY_RUN == true ]]; then
    printf '%q %q close-run --manifest %q\n' "$PYTHON" "$PLANNER" "$MANIFEST"
    exit 0
fi
exec "$PYTHON" "$PLANNER" close-run --manifest "$MANIFEST"
