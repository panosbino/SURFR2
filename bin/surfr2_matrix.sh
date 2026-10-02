#!/usr/bin/env bash
# =============================================================================
# SURFR2 - candidate k-mer set and per-sample count extraction
#
# A dense k-mer x sample matrix is infeasible (tens of millions of distinct
# 17-mers x hundreds of samples), so set operations stay in KMC's binary format:
#
#  A. Candidates: k-mers with >= CAND_MIN_COUNT reads in >= CAND_MIN_SAMPLES case OR
#     control samples. Selection is CONDITION-BLIND on purpose: selecting on case
#     counts and then testing case vs control on the same data would bias the
#     p-values (Bourgon et al. 2010, PNAS). Each DB is reduced to an indicator
#     (count -> 1), indicators are summed by a pairwise union tree (= number of
#     samples carrying the k-mer), and the sum is thresholded.
#  B. For EVERY sample (case, control, external), the counts of the candidate
#     k-mers are extracted by intersection, keeping the sample's own count.
#     Controls were counted with -ci1, so any trace of a candidate is retained.
#  C. median_ratio only: reference k-mers present in ALL samples (intersection
#     tree), with per-sample counts, for size-factor estimation in R.
#
# Usage: surfr2_matrix.sh <run_config_dir>
#
# Outputs in <outdir>/matrix/:
#   library_sizes.tsv                 all samples (from per-sample step)
#   candidates.tsv                    kmer <TAB> n_samples_with_>=min_count (case+control)
#   counts/candidates/<sample>.tsv.gz kmer <TAB> count  (candidate k-mers present)
#   counts/reference/<sample>.tsv.gz  (median_ratio only)
# =============================================================================

set -euo pipefail
export LC_ALL=C

die() { echo "[$(date '+%F %T')] ERROR [matrix]: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] [matrix] $*" >&2; }

[ "$#" -eq 1 ] || die "usage: $0 <run_config_dir>"
RUN_CONFIG_DIR=$1
# shellcheck source=/dev/null
source "${RUN_CONFIG_DIR}/params.env"

THREADS=${SLURM_CPUS_PER_TASK:-$(nproc)}
# kmc_tools set operations are largely I/O-bound: run several in parallel with
# a couple of threads each rather than one job with many threads.
PAR=$(( THREADS >= 4 ? THREADS / 2 : 1 ))
T_EACH=$(( THREADS / PAR )); [ "${T_EACH}" -ge 1 ] || T_EACH=1
# Counter ceiling for indicator sums: cannot exceed the number of samples.
# (kmc_tools parses -cs with atoi, so values >= 2^31 would overflow.)
CS_IND=$(( N_SAMPLES + 1 ))
# Ceiling for real per-sample counts when copying them through an intersection.
CS_COUNTS=2147483647

MDIR="${OUTDIR}/matrix"
WORK="${MDIR}/work"
SAMPLES_DIR="${OUTDIR}/samples"

mapfile -t ALL_IDS  < <(awk -F'\t' 'NR > 1 { print $2 }' "${RUN_CONFIG_DIR}/samples.tsv")
mapfile -t TEST_IDS < <(awk -F'\t' 'NR > 1 && ($5 == "case" || $5 == "control") { print $2 }' "${RUN_CONFIG_DIR}/samples.tsv")
[ "${#ALL_IDS[@]}" -eq "${N_SAMPLES}" ] || die "samples.tsv / params.env mismatch"

# Fingerprint: pre-filter/normalisation parameters, sample roles, and every sample's own
# fingerprint (which changes whenever that sample is recomputed).
FINGERPRINT=$( {
    echo "surfr2-matrix-v2 min_count=${CAND_MIN_COUNT} min_samples=${CAND_MIN_SAMPLES} norm=${NORM_METHOD}"
    awk -F'\t' 'NR > 1 { print $2, $5 }' "${RUN_CONFIG_DIR}/samples.tsv"
    for id in "${ALL_IDS[@]}"; do cat "${SAMPLES_DIR}/${id}/.done" 2>/dev/null || echo "${id} missing"; done
} | sha1sum | cut -d' ' -f1 )

if [ -f "${MDIR}/.done" ] && [ "$(cat "${MDIR}/.done")" = "${FINGERPRINT}" ]; then
    log "already complete with identical inputs and parameters - skipping"
    exit 0
fi
rm -f "${MDIR}/.done"
rm -rf "${WORK}" "${MDIR}/counts"
mkdir -p "${WORK}" "${MDIR}/counts/candidates"

db_ok() { [ -s "$1.kmc_pre" ] && [ -s "$1.kmc_suf" ]; }
db_rm() { rm -f "$1.kmc_pre" "$1.kmc_suf"; }

# Run each line of a job file as an independent bash command, PAR at a time.
# xargs exits non-zero if any job fails, which aborts this script (set -e).
run_parallel() {
    [ -s "$1" ] || return 0
    xargs -d '\n' -P "${PAR}" -I CMD bash -c 'set -euo pipefail; CMD' < "$1"
}

# tree_reduce <union|intersect> <-oc mode> <-cs value> <tag> <out_db> <in_db>...
# Balanced pairwise reduction: log2(n) levels, pairs within a level in parallel.
tree_reduce() {
    local op=$1 mode=$2 cs=$3 tag=$4 out=$5; shift 5
    local -a cur=("$@") next
    local level=0 i jobs tdir="${WORK}/tree_${tag}"
    [ "${#cur[@]}" -ge 1 ] || die "tree_reduce(${tag}): no inputs"
    mkdir -p "${tdir}"
    while [ "${#cur[@]}" -gt 1 ]; do
        next=(); jobs="${tdir}/L${level}.jobs"; : > "${jobs}"
        for (( i = 0; i < ${#cur[@]}; i += 2 )); do
            if (( i + 1 < ${#cur[@]} )); then
                local o="${tdir}/L${level}_${i}"
                echo "'${KMC_TOOLS}' -hp -t${T_EACH} simple '${cur[i]}' '${cur[i+1]}' ${op} '${o}' ${mode} -cs${cs} >/dev/null" >> "${jobs}"
                next+=("${o}")
            else
                next+=("${cur[i]}")
            fi
        done
        log "  ${tag}: level ${level}, ${#cur[@]} -> ${#next[@]} databases"
        run_parallel "${jobs}"
        # Free the previous level's intermediates - but never original inputs, and never
        # the odd database carried unchanged into the next level (it is still needed).
        if [ "${level}" -gt 0 ]; then
            local carried="${next[${#next[@]}-1]}"
            for db in "${cur[@]}"; do
                [ "${db}" = "${carried}" ] && continue
                case "${db}" in "${tdir}"/*) db_rm "${db}" ;; esac
            done
        fi
        cur=("${next[@]}"); level=$(( level + 1 ))
    done
    if [ "${level}" -eq 0 ]; then       # single input: copy, keep original intact
        cp "${cur[0]}.kmc_pre" "${out}.kmc_pre"; cp "${cur[0]}.kmc_suf" "${out}.kmc_suf"
    else
        mv "${cur[0]}.kmc_pre" "${out}.kmc_pre"; mv "${cur[0]}.kmc_suf" "${out}.kmc_suf"
    fi
    rm -rf "${tdir}"
}

# extract_counts <filter_db> <outsubdir>: per-sample counts of the k-mers in filter_db
extract_counts() {
    local filt=$1 sub=$2 jobs="${WORK}/extract_${2}.jobs" id
    mkdir -p "${MDIR}/counts/${sub}" "${WORK}/x_${sub}"
    : > "${jobs}"
    for id in "${ALL_IDS[@]}"; do
        local x="${WORK}/x_${sub}/${id}" out="${MDIR}/counts/${sub}/${id}.tsv.gz"
        echo "'${KMC_TOOLS}' -hp -t${T_EACH} simple '${SAMPLES_DIR}/${id}/kmc/${id}' '${filt}' intersect '${x}' -ocleft -cs${CS_COUNTS} >/dev/null && '${KMC_TOOLS}' -hp -t${T_EACH} transform '${x}' dump -s '${x}.txt' >/dev/null && '${PIGZ}' -p1 -c '${x}.txt' > '${out}.part' && mv '${out}.part' '${out}' && rm -f '${x}.txt' '${x}.kmc_pre' '${x}.kmc_suf'" >> "${jobs}"
    done
    run_parallel "${jobs}"
    rm -rf "${WORK}/x_${sub}"
}

# -----------------------------------------------------------------------------
# 0. Check per-sample outputs and collect library sizes
# -----------------------------------------------------------------------------
missing=0
for id in "${ALL_IDS[@]}"; do
    if [ ! -f "${SAMPLES_DIR}/${id}/.done" ] || ! db_ok "${SAMPLES_DIR}/${id}/kmc/${id}"; then
        echo "  incomplete: ${id}" >&2; missing=$(( missing + 1 ))
    fi
done
[ "${missing}" -eq 0 ] || die "${missing} sample(s) incomplete - re-run the sample step"

# Merge by column NAME: samples processed by different SURFR2 versions can have different
# columns (e.g. without umi_molecules); concatenating under one header would misalign them.
python3 - "${MDIR}/library_sizes.tsv" "${ALL_IDS[@]/#/${SAMPLES_DIR}/}" <<'PY'
import csv, sys
out, dirs = sys.argv[1], sys.argv[2:]
rows, cols = [], []
for d in dirs:
    sid = d.rstrip("/").split("/")[-1]
    with open(f"{d}/library_size.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
            cols += [c for c in r if c not in cols]
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="NA", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
PY
log "${N_SAMPLES} samples (${#TEST_IDS[@]} case/control); threads=${THREADS} parallel=${PAR}x${T_EACH}"

# -----------------------------------------------------------------------------
# A. Candidate set
# -----------------------------------------------------------------------------
log "A. candidates (condition-blind): >= ${CAND_MIN_COUNT} reads in >= ${CAND_MIN_SAMPLES} case/control samples"
mkdir -p "${WORK}/ind"
jobs="${WORK}/ind.jobs"; : > "${jobs}"
ind_dbs=()
for id in "${TEST_IDS[@]}"; do
    # input -ci drops k-mers below min_count BEFORE set_counts turns counts into 1
    echo "'${KMC_TOOLS}' -hp -t${T_EACH} transform '${SAMPLES_DIR}/${id}/kmc/${id}' -ci${CAND_MIN_COUNT} set_counts 1 '${WORK}/ind/${id}' >/dev/null" >> "${jobs}"
    ind_dbs+=("${WORK}/ind/${id}")
done
run_parallel "${jobs}"

tree_reduce union -ocsum "${CS_IND}" prevalence "${WORK}/prevalence" "${ind_dbs[@]}"
rm -rf "${WORK}/ind"

"${KMC_TOOLS}" -hp -t"${THREADS}" transform "${WORK}/prevalence" -ci"${CAND_MIN_SAMPLES}" \
    reduce "${WORK}/candidates" dump -s "${MDIR}/candidates.tsv" > /dev/null
db_rm "${WORK}/prevalence"
N_CAND=$(wc -l < "${MDIR}/candidates.tsv")
[ "${N_CAND}" -gt 0 ] || die "no candidate k-mers - candidates.* thresholds are too strict for these data"
log "   ${N_CAND} candidate k-mers"

# -----------------------------------------------------------------------------
# B. Candidate counts for every sample
# -----------------------------------------------------------------------------
log "B. extracting candidate counts for ${N_SAMPLES} samples"
extract_counts "${WORK}/candidates" candidates

# -----------------------------------------------------------------------------
# C. Reference set for median-of-ratios normalisation
# -----------------------------------------------------------------------------
if [ "${NORM_METHOD}" = "median_ratio" ]; then
    log "C. reference k-mers present in all ${N_SAMPLES} samples"
    all_dbs=()
    for id in "${ALL_IDS[@]}"; do all_dbs+=("${SAMPLES_DIR}/${id}/kmc/${id}"); done
    tree_reduce intersect -ocmin "${CS_COUNTS}" reference "${WORK}/reference" "${all_dbs[@]}"
    "${KMC_TOOLS}" -hp -t"${THREADS}" transform "${WORK}/reference" dump -s "${WORK}/reference.txt" > /dev/null
    N_REF=$(wc -l < "${WORK}/reference.txt")
    [ "${N_REF}" -gt 0 ] || die "no k-mer is present in all samples; median_ratio is impossible (a failed library?) - check library_sizes.tsv or use method: cpm"
    log "   ${N_REF} reference k-mers"
    extract_counts "${WORK}/reference" reference
fi

rm -rf "${WORK}"
echo "${FINGERPRINT}" > "${MDIR}/.done"
log "done"
