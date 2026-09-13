#!/usr/bin/env bash
# Rocky9 exact-four assembly structural validation unit.
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
if [[ $# -ne 1 || ( $DRY_RUN == false && ! ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-3]$ ) ]]; then
    printf 'usage: 01_cnv_array_unit.sh [--dry-run] MANIFEST.json\n' >&2
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
    printf '%q %q run-unit --manifest %q --unit-index %q\n' "$PYTHON" "$PLANNER" "$MANIFEST" '${SLURM_ARRAY_TASK_ID}'
    exit 0
fi

module purge
module load samtools/1.21
[[ $(command -v samtools) == /path/to/software/samtools/1.21/bin/samtools ]] || {
    printf 'unexpected samtools path\n' >&2
    exit 2
}
exec "$PYTHON" "$PLANNER" run-unit --manifest "$MANIFEST" --unit-index "$SLURM_ARRAY_TASK_ID"
