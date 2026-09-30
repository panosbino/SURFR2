#!/usr/bin/env bash
# =============================================================================
# SURFR2 - per-sample processing (one SLURM array task per sample)
#
#   input (BAM or FASTQ) -> FASTQ -> miRTrace QC -> uncollapsed FASTA
#                        -> KMC k-mer database + library-size record
#
# Usage:
#   surfr2_sample.sh <run_config_dir> [sample_index]
#   (sample_index defaults to $SLURM_ARRAY_TASK_ID; 1-based row of samples.tsv)
#
# Outputs in <outdir>/samples/<sample_id>/:
#   kmc/<sample_id>.kmc_pre|.kmc_suf  per-sample k-mer counts (-ci1: nothing hidden)
#   kmc_stats.json                    KMC summary (total k-mers = library size)
#   library_size.tsv                  one-row table consumed by the R step
#   mirtrace/                         miRTrace QC report (FASTA removed unless kept)
#   .done                             parameter fingerprint, written last; the task is
#                                     skipped only if it matches the current parameters
#
# Runs inside the container; does not load modules or call singularity itself.
# =============================================================================

set -euo pipefail
export LC_ALL=C

die() { echo "[$(date '+%F %T')] ERROR [${SAMPLE:-?}]: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] [${SAMPLE:-?}] $*" >&2; }

[ "$#" -ge 1 ] || die "usage: $0 <run_config_dir> [sample_index]"
RUN_CONFIG_DIR=$1
IDX=${2:-}
if [ -z "${IDX}" ]; then
    [ -n "${SLURM_ARRAY_TASK_ID:-}" ] || die "no sample index given and SLURM_ARRAY_TASK_ID is unset"
    # Offset lets the launcher split >MaxArraySize samples over several arrays
    IDX=$(( SLURM_ARRAY_TASK_ID + ${SURFR2_INDEX_OFFSET:-0} ))
fi
[ -f "${RUN_CONFIG_DIR}/params.env" ] || die "missing ${RUN_CONFIG_DIR}/params.env"
# shellcheck source=/dev/null
source "${RUN_CONFIG_DIR}/params.env"

THREADS=${SLURM_CPUS_PER_TASK:-$(nproc)}

row=$(awk -F'\t' -v i="${IDX}" 'NR > 1 && $1 == i' "${RUN_CONFIG_DIR}/samples.tsv")
[ -n "${row}" ] || die "sample index ${IDX} not found in samples.tsv"
IFS=$'\t' read -r _ SAMPLE COHORT CONDITION ROLE FTYPE INPUT <<< "${row}"

SDIR="${OUTDIR}/samples/${SAMPLE}"
[ -f "${INPUT}" ] || die "input not found: ${INPUT}"

# Everything that determines this sample's outputs. A finished sample is reused only
# if its fingerprint matches, so changing k, QC settings or the input file forces a redo.
FINGERPRINT="surfr2-sample-v1 input=${INPUT} type=${FTYPE} size=$(stat -L -c %s "${INPUT}") mtime=$(stat -L -c %Y "${INPUT}") k=${K} canonical=${CANONICAL} species=${MIRTRACE_SPECIES} adapter=${MIRTRACE_ADAPTER} phred=${MIRTRACE_PHRED}"
if [ -f "${SDIR}/.done" ]; then
    if [ "$(cat "${SDIR}/.done")" = "${FINGERPRINT}" ]; then
        log "already complete with identical parameters - skipping"
        exit 0
    fi
    log "previous run used different parameters or input - recomputing"
    rm -f "${SDIR}/.done"
fi
mkdir -p "${SDIR}"
TMP="${SDIR}/tmp"          # on the project filesystem: Dardel nodes have no local scratch
rm -rf "${TMP}" "${SDIR}/kmc" "${SDIR}/mirtrace"
mkdir -p "${TMP}/kmc_tmp" "${SDIR}/kmc"
trap 'rm -rf "${TMP}"' EXIT

log "cohort=${COHORT} condition=${CONDITION} role=${ROLE} type=${FTYPE} threads=${THREADS}"

# -----------------------------------------------------------------------------
# 1. Obtain FASTQ
# -----------------------------------------------------------------------------
case "${FTYPE}" in
    bam)
        # Small-RNA libraries are single-end. Refuse paired data rather than
        # silently counting both mates (check the first 10k records only).
        n_paired=$("${SAMTOOLS}" view "${INPUT}" 2>/dev/null | head -n 10000 \
                   | awk '{ if (and($2, 1)) n++ } END { print n + 0 }') || true
        [ "${n_paired:-0}" -eq 0 ] || die "BAM contains paired-end reads; SURFR2 expects single-end small-RNA data"

        FASTQ="${TMP}/${SAMPLE}.fastq.gz"
        log "BAM -> FASTQ"
        # -F 0x900: drop secondary/supplementary alignments so each read is counted once.
        # samtools fastq reverse-complements reads aligned to '-', restoring the
        # original read orientation - essential because k-mers are counted stranded.
        "${SAMTOOLS}" fastq -F 0x900 -@ "${THREADS}" "${INPUT}" 2> "${SDIR}/samtools_fastq.log" \
            | "${PIGZ}" -p "${THREADS}" > "${FASTQ}"
        ;;
    fastq)
        case "${INPUT}" in
            *.gz) FASTQ="${TMP}/${SAMPLE}.fastq.gz" ;;
            *)    FASTQ="${TMP}/${SAMPLE}.fastq"    ;;
        esac
        ln -s "${INPUT}" "${FASTQ}"
        ;;
    *) die "unsupported file_type '${FTYPE}'" ;;
esac

# -----------------------------------------------------------------------------
# 2. miRTrace QC -> uncollapsed FASTA of QC-passed reads
# -----------------------------------------------------------------------------
log "miRTrace QC"
mt_args=(qc --species "${MIRTRACE_SPECIES}" --output-dir "${SDIR}/mirtrace"
         --write-fasta --uncollapse-fasta --num-threads "${THREADS}" --force)
if [ -n "${MIRTRACE_PHRED}" ]; then
    # A PHRED offset can only be given through miRTrace's per-sample config CSV
    # (path,name,adapter,phred); in that mode the phred field is mandatory.
    printf '%s,%s,%s,%s\n' "${FASTQ}" "${SAMPLE}" "${MIRTRACE_ADAPTER}" "${MIRTRACE_PHRED}" \
        > "${TMP}/mirtrace_config.csv"
    mt_args+=(--config "${TMP}/mirtrace_config.csv")
else
    if [ -n "${MIRTRACE_ADAPTER}" ]; then mt_args+=(--adapter "${MIRTRACE_ADAPTER}"); fi
    mt_args+=("${FASTQ}")
fi
"${MIRTRACE}" "${mt_args[@]}" > "${SDIR}/mirtrace.log" 2>&1 \
    || die "miRTrace failed; see ${SDIR}/mirtrace.log"
# miRTrace 1.0.1 can exit 0 after aborting (e.g. PHRED auto-detection failure)
if grep -q -E '^ERROR|Error parsing|aborting' "${SDIR}/mirtrace.log"; then
    die "miRTrace reported an error (exit status was 0): $(grep -m1 -E '^ERROR|Error parsing|Could not' "${SDIR}/mirtrace.log")"
fi

# One input per run, so exactly one FASTA is expected. Locate it rather than
# predicting miRTrace's naming convention.
mapfile -t fastas < <(find "${SDIR}/mirtrace/qc_passed_reads.all.uncollapsed" \
                           -maxdepth 1 -type f -name '*.fasta*' 2>/dev/null)
[ "${#fastas[@]}" -eq 1 ] || die "expected 1 miRTrace FASTA, found ${#fastas[@]}"

FASTA="${TMP}/${SAMPLE}.qc.fasta.gz"
case "${fastas[0]}" in
    *.gz) mv "${fastas[0]}" "${FASTA}" ;;
    *)    "${PIGZ}" -p "${THREADS}" -c "${fastas[0]}" > "${FASTA}" && rm -f "${fastas[0]}" ;;
esac

QC_READS=$("${PIGZ}" -dc "${FASTA}" | awk 'substr($0, 1, 1) == ">" { n++ } END { print n + 0 }')
[ "${QC_READS}" -gt 0 ] || die "no reads passed miRTrace QC"
log "QC-passed reads: ${QC_READS}"

# -----------------------------------------------------------------------------
# 3. KMC k-mer counting
#   -ci1  keep every k-mer. SURFR1 used -ci30 on pooled cancer data; per sample,
#         any cut-off here would silently turn low control counts into zeros and
#         inflate 'specificity'. Filtering happens later, on normalised values.
#   -b    no canonical form (stranded library)          [SURFR1 comment was wrong]
#   -cs   32-bit counter ceiling (2^32 - 1)
#   -j    must be glued to its value; with a space KMC mis-parses the arguments
# -----------------------------------------------------------------------------
log "KMC k=${K}"
kmc_args=(-hp -fa "-k${K}" -ci1 -cs4294967295 "-m${KMC_MEMORY_GB}" "-t${THREADS}"
          "-j${SDIR}/kmc_stats.json")
if [ "${CANONICAL}" = "false" ]; then kmc_args+=(-b); fi
"${KMC}" "${kmc_args[@]}" "${FASTA}" "${SDIR}/kmc/${SAMPLE}" "${TMP}/kmc_tmp" \
    > "${SDIR}/kmc.log" 2>&1 || die "KMC failed; see ${SDIR}/kmc.log"

read -r TOTAL_KMERS UNIQUE_KMERS KMC_READS < <(python3 - "${SDIR}/kmc_stats.json" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))["Stats"]
print(s["#Total no. of k-mers"], s["#Unique_counted_k-mers"], s["#Total_reads"])
PY
)
[ "${TOTAL_KMERS}" -gt 0 ] || die "KMC counted zero k-mers"
# KMC must see exactly the reads that passed QC; a mismatch means a parsing problem
# (e.g. multi-line FASTA) and would corrupt the library size.
[ "${KMC_READS}" -eq "${QC_READS}" ] \
    || die "KMC read count (${KMC_READS}) != QC-passed reads (${QC_READS})"

# -----------------------------------------------------------------------------
# 4. Library-size record, intermediates, completion marker
# -----------------------------------------------------------------------------
{
    printf 'sample_id\tcohort\tcondition\trole\tqc_reads\ttotal_kmers\tunique_kmers\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${SAMPLE}" "${COHORT}" "${CONDITION}" "${ROLE}" \
        "${QC_READS}" "${TOTAL_KMERS}" "${UNIQUE_KMERS}"
} > "${SDIR}/library_size.tsv.part"
mv "${SDIR}/library_size.tsv.part" "${SDIR}/library_size.tsv"

if [ "${KEEP_INTERMEDIATES}" = "true" ]; then
    mv "${FASTA}" "${SDIR}/"
    if [ "${FTYPE}" = "bam" ]; then mv "${FASTQ}" "${SDIR}/"; fi
fi

printf '%s\n' "${FINGERPRINT}" > "${SDIR}/.done.part" && mv "${SDIR}/.done.part" "${SDIR}/.done"
log "done: total_kmers=${TOTAL_KMERS} unique_kmers=${UNIQUE_KMERS}"
