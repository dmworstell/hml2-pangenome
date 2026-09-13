#!/usr/bin/env bash
# Short receipt-census controller.  It allocates no long-lived polling worker
# and never encodes a parent job status as scientific readiness.
#SBATCH --partition=batch,preempt
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -Eeuo pipefail
IFS=$'\n\t'

plan=""
analysis=""
while (($#)); do
    case "$1" in
        --plan) plan=$2; shift 2 ;;
        --analysis) analysis=$2; shift 2 ;;
        *) printf 'CNV-v2 receipt controller: unknown argument: %s\n' "$1" >&2; exit 2 ;;
    esac
done
[[ -s "$plan" && -s "$analysis" ]] || {
    printf 'CNV-v2 receipt controller: plan or analysis is unavailable\n' >&2
    exit 2
}

eval "$(python3 - "$plan" <<'PY'
import json, os, shlex, sys
p = json.load(open(sys.argv[1], encoding="utf-8"))
receipts = [part.get("receipt", "") for part in p.get("parts", [])]
if not receipts or any(not isinstance(x, str) or not os.path.isabs(x) for x in receipts):
    raise SystemExit("plan has no absolute part receipt census")
ready = sum(os.path.isfile(x) and os.path.getsize(x) > 0 for x in receipts)
print("EXPECTED=" + str(len(receipts)))
print("READY=" + str(ready))
print("TARGET_ID=" + shlex.quote(str(p.get("target_id", "unknown"))))
PY
)"

if [[ "$READY" -eq "$EXPECTED" ]]; then
    submitted=$(sbatch --parsable --job-name="cnv2_${TARGET_ID}_analysis" \
        "$analysis" --plan "$plan")
    submitted=${submitted%%;*}
    [[ "$submitted" =~ ^[0-9]+$ ]]
    printf 'CNV_V2_ANALYSIS_RELEASED_FROM_RECEIPTS job=%s ready=%s\n' \
        "$submitted" "$READY"
    exit 0
fi

next=$(sbatch --parsable --begin=now+5minutes \
    --job-name="cnv2_${TARGET_ID}_receipts" "$0" \
    --plan "$plan" --analysis "$analysis")
printf 'CNV_V2_RECEIPTS_PENDING ready=%s expected=%s next=%s\n' \
    "$READY" "$EXPECTED" "${next%%;*}"
