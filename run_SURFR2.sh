#!/usr/bin/env bash
# =============================================================================
# SURFR2 - launcher
#
# Validates the config, then runs or submits:
#   1. sample  - per-sample FASTQ -> miRTrace -> KMC     (SLURM array, 1 task/sample)
#   2. matrix  - candidate k-mers + per-sample extraction (1 job, after all samples)
#   3. filter  - normalisation, per-cohort filters, mergeTags (1 job, after matrix)
#   4. report  - static + interactive (HTML) scatterplots      (1 job, after filter)
#
# Resumable: re-running the launcher skips any sample/step whose outputs exist AND
# were produced with identical parameters and inputs (fingerprinted .done markers).
# After a partial failure, fix the cause and simply re-launch.
#
# Run on the login node (not via sbatch):
#   bash run_SURFR2.sh -c config.yaml [options]
#
# Options:
#   -c, --config FILE     SURFR2 YAML config (required)
#   --container PATH      override execution.container (also used to run the
#                         validator when the host python lacks PyYAML)
#   --from STEP           start at sample (default) | matrix | filter | report
#   --dry-run             validate and print the commands; submit/run nothing
#   --no-file-check       skip existence checks of input files (dry runs off-cluster)
# =============================================================================

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/bin"

CONFIG="" CONTAINER_OVERRIDE="__unset__" FROM=sample DRY_RUN=false FILE_CHECK=true
while [ "$#" -gt 0 ]; do
    case "$1" in
        -c|--config)     CONFIG=${2:?}; shift 2 ;;
        --container)     CONTAINER_OVERRIDE=${2?}; shift 2 ;;
        --from)          FROM=${2:?}; shift 2 ;;
        --dry-run)       DRY_RUN=true; shift ;;
        --no-file-check) FILE_CHECK=false; shift ;;
        -h|--help)       sed -n "2,28p" "$0"; exit 0 ;;
        *) die "unknown argument: $1 (see --help)" ;;
    esac
done
[ -n "${CONFIG}" ] || die "missing -c/--config (see --help)"
[ -f "${CONFIG}" ] || die "config not found: ${CONFIG}"
case "${FROM}" in sample|matrix|filter|report) ;; *) die "--from must be sample, matrix, filter or report" ;; esac
for s in surfr2_config.py surfr2_sample.sh surfr2_matrix.sh surfr2_filter.R surfr2_report.R surfr2_plotlib.R; do
    [ -f "${BIN}/${s}" ] || die "missing pipeline script ${BIN}/${s}"
done

# -----------------------------------------------------------------------------
# 1. Validate and resolve the configuration
# -----------------------------------------------------------------------------
if python3 -c 'import yaml' 2>/dev/null; then
    VALIDATE=(python3 "${BIN}/surfr2_config.py")
elif [ "${CONTAINER_OVERRIDE}" != "__unset__" ] && command -v singularity > /dev/null; then
    VALIDATE=(singularity exec "${CONTAINER_OVERRIDE}" python3 "${BIN}/surfr2_config.py")
else
    die "host python3 lacks PyYAML: pass --container <sandbox> (after 'ml singularity') or 'pip install --user pyyaml'"
fi

validate() {   # validate <out_dir>
    local extra=()
    [ "${FILE_CHECK}" = true ] || extra+=(--no-file-check)
    [ "${CONTAINER_OVERRIDE}" = "__unset__" ] || extra+=(--container "${CONTAINER_OVERRIDE}")
    "${VALIDATE[@]}" "${CONFIG}" --out "$1" "${extra[@]}"
}

PROBE=$(mktemp -d)
trap 'rm -rf "${PROBE}"' EXIT
validate "${PROBE}"
OUTDIR=$(. "${PROBE}/params.env" && echo "${OUTDIR}")

if [ "${DRY_RUN}" = true ]; then
    RUN_DIR="${PROBE}"
else
    RUN_DIR="${OUTDIR}/runs/$(date +%Y%m%d-%H%M%S)"
    mkdir -p "${RUN_DIR}/logs"
    validate "${RUN_DIR}" 2>/dev/null      # re-resolve so params.env records RUN_DIR
    cp "${CONFIG}" "${RUN_DIR}/config.submitted.yaml"
    ln -sfn "${RUN_DIR}" "${OUTDIR}/runs/latest"
fi
# shellcheck source=/dev/null
source "${RUN_DIR}/params.env"
LOGS="${RUN_DIR}/logs"

# Command prefix that runs a step inside the container (if one is configured)
EXEC=""
if [ -n "${CONTAINER}" ]; then
    EXEC="singularity exec"
    [ -z "${BIND}" ] || EXEC="${EXEC} -B ${BIND}"
    EXEC="${EXEC} ${CONTAINER}"
fi

echo "============================================================"
echo " SURFR2  project=${PROJECT}  samples=${N_SAMPLES}  executor=${EXECUTOR}"
echo " outdir:  ${OUTDIR}"
echo " run dir: ${RUN_DIR}$([ "${DRY_RUN}" = true ] && echo '  (dry run: temporary)')"
echo " from:    ${FROM}"
echo "============================================================"

# -----------------------------------------------------------------------------
# 2a. Local execution
# -----------------------------------------------------------------------------
if [ "${EXECUTOR}" = "local" ]; then
    case "${FROM}" in
        sample) steps=(sample matrix filter report) ;;
        matrix) steps=(matrix filter report) ;;
        filter) steps=(filter report) ;;
        report) steps=(report) ;;
    esac
    for step in "${steps[@]}"; do
        case "${step}" in
            sample) cmd="for i in \$(seq 1 ${N_SAMPLES}); do ${EXEC} bash ${BIN}/surfr2_sample.sh ${RUN_DIR} \$i; done" ;;
            matrix) cmd="${EXEC} bash ${BIN}/surfr2_matrix.sh ${RUN_DIR}" ;;
            filter) cmd="${EXEC} ${RSCRIPT} ${BIN}/surfr2_filter.R ${RUN_DIR}" ;;
            report) cmd="${EXEC} ${RSCRIPT} ${BIN}/surfr2_report.R ${RUN_DIR}" ;;
        esac
        if [ "${DRY_RUN}" = true ]; then
            echo "[dry-run] ${step}: ${cmd}"
        else
            echo "[$(date '+%F %T')] ${step} (log: ${LOGS}/${step}.log)"
            bash -c "set -euo pipefail; ${cmd}" > "${LOGS}/${step}.log" 2>&1 \
                || die "${step} failed - see ${LOGS}/${step}.log"
        fi
    done
    [ "${DRY_RUN}" = true ] || echo "Done. Results: ${OUTDIR}/results (summary: run_summary.txt)"
    exit 0
fi

# -----------------------------------------------------------------------------
# 2b. SLURM execution
# -----------------------------------------------------------------------------
PRE=""
[ -z "${SLURM_MODULES}" ] || PRE="ml ${SLURM_MODULES} && "   # load modules on the HOST node

submit() {   # submit <sbatch args...>; prints job id
    if [ "${DRY_RUN}" = true ]; then
        printf '[dry-run] sbatch' >&2; printf ' %q' "$@" >&2; echo >&2
        echo "DRYRUN$((RANDOM))"
    else
        sbatch --parsable "$@" | cut -d';' -f1
    fi
}
common=(--account="${SLURM_ACCOUNT}" --kill-on-invalid-dep=yes)

dep=""
if [ "${FROM}" = sample ]; then
    # Array indices must be < MaxArraySize: split large sample sets over several arrays
    max_array=$( (scontrol show config 2>/dev/null | awk '/^MaxArraySize/ { print $3 }') || true)
    chunk=$(( ${max_array:-1001} - 1 ))
    array_ids=()
    for (( offset = 0; offset < N_SAMPLES; offset += chunk )); do
        n=$(( N_SAMPLES - offset < chunk ? N_SAMPLES - offset : chunk ))
        jid=$(submit "${common[@]}" \
            --partition="${SLURM_SAMPLE_PARTITION}" --cpus-per-task="${SLURM_SAMPLE_CPUS}" \
            --time="${SLURM_SAMPLE_TIME}" --array="1-${n}%${SLURM_ARRAY_THROTTLE}" \
            --job-name="SURFR2_${PROJECT}_sample" --output="${LOGS}/sample_$((offset))+%a.log" \
            --export=ALL,SURFR2_INDEX_OFFSET="${offset}" \
            --wrap="${PRE}${EXEC} bash ${BIN}/surfr2_sample.sh ${RUN_DIR}")
        echo "sample array (samples $((offset + 1))-$((offset + n))): ${jid}"
        array_ids+=("${jid}")
    done
    dep="--dependency=afterok:$(IFS=:; echo "${array_ids[*]}")"
fi

if [ "${FROM}" = sample ] || [ "${FROM}" = matrix ]; then
    jid=$(submit "${common[@]}" ${dep:+"${dep}"} \
        --partition="${SLURM_MATRIX_PARTITION}" --cpus-per-task="${SLURM_MATRIX_CPUS}" \
        --time="${SLURM_MATRIX_TIME}" --job-name="SURFR2_${PROJECT}_matrix" \
        --output="${LOGS}/matrix.log" \
        --wrap="${PRE}${EXEC} bash ${BIN}/surfr2_matrix.sh ${RUN_DIR}")
    echo "matrix: ${jid}"
    dep="--dependency=afterok:${jid}"
fi

if [ "${FROM}" != report ]; then
    jid=$(submit "${common[@]}" ${dep:+"${dep}"} \
        --partition="${SLURM_FILTER_PARTITION}" --cpus-per-task="${SLURM_FILTER_CPUS}" \
        --time="${SLURM_FILTER_TIME}" --job-name="SURFR2_${PROJECT}_filter" \
        --output="${LOGS}/filter.log" \
        --wrap="${PRE}${EXEC} ${RSCRIPT} ${BIN}/surfr2_filter.R ${RUN_DIR}")
    echo "filter: ${jid}"
    dep="--dependency=afterok:${jid}"
fi

# Plots are cheap: short job on the filter partition
jid=$(submit "${common[@]}" ${dep:+"${dep}"} \
    --partition="${SLURM_FILTER_PARTITION}" --cpus-per-task=1 \
    --time="01:00:00" --job-name="SURFR2_${PROJECT}_report" \
    --output="${LOGS}/report.log" \
    --wrap="${PRE}${EXEC} ${RSCRIPT} ${BIN}/surfr2_report.R ${RUN_DIR}")
echo "report: ${jid}"
echo
echo "Logs: ${LOGS}   Monitor: squeue -u \$USER -n SURFR2_${PROJECT}_sample,SURFR2_${PROJECT}_matrix,SURFR2_${PROJECT}_filter,SURFR2_${PROJECT}_report"
echo "If a sample task fails, downstream jobs are cancelled; fix the cause and re-launch (finished work is reused)."
