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
IFS=$'\t' read -r _ SAMPLE COHORT CONDITION ROLE FTYPE INPUT _BLOCK <<< "${row}"

SDIR="${OUTDIR}/samples/${SAMPLE}"
[ -f "${INPUT}" ] || die "input not found: ${INPUT}"

# Everything that determines this sample's outputs. A finished sample is reused only
# if its fingerprint matches, so changing k, QC settings or the input file forces a redo.
FINGERPRINT="surfr2-sample-v1 input=${INPUT} type=${FTYPE} size=$(stat -L -c %s "${INPUT}") mtime=$(stat -L -c %Y "${INPUT}") k=${K} canonical=${CANONICAL} species=${MIRTRACE_SPECIES} adapter=${MIRTRACE_ADAPTER} phred=${MIRTRACE_PHRED}"
# appended only when used, so samples processed before UMI support keep their fingerprint
if [ "${UMI_LENGTH}" -gt 0 ]; then FINGERPRINT="${FINGERPRINT} umi=${UMI_LENGTH}"; fi

# Artifact references: built-in files beside the installed tools, plus extra FASTA files.
# Their content (checksums) and the adapter enter the fingerprint, so changing any of
# them recomputes the sample.
ART_FILES=()
if [ "${ARTIFACTS_ENABLED}" = "true" ]; then
    if [ "${ARTIFACTS_BUILTIN}" = "true" ]; then
        kmc_path=$(command -v "${KMC}") || die "kmc not found"
        art_dir="$(cd "$(dirname "$(readlink -f "${kmc_path}")")/.." && pwd)/share/surfr2/artifacts"
        for f in illumina_adapters.fa phix174.fa.gz; do
            [ -s "${art_dir}/${f}" ] || die "built-in artifact reference missing: ${art_dir}/${f} (rerun container/install_tools.sh, or set artifacts.builtin: false)"
            ART_FILES+=("${art_dir}/${f}")
        done
    fi
    if [ -n "${ARTIFACTS_EXTRA}" ]; then
        IFS=':' read -r -a extra <<< "${ARTIFACTS_EXTRA}"
        for f in "${extra[@]}"; do [ -s "${f}" ] || die "artifact FASTA missing: ${f}"; ART_FILES+=("${f}"); done
    fi
    ART_ID=$( { for f in "${ART_FILES[@]}"; do sha256sum "${f}" | cut -c1-64; done; echo "adapter=${MIRTRACE_ADAPTER}"; } \
              | sha256sum | cut -c1-16 )
    FINGERPRINT="${FINGERPRINT} artifacts=${ART_ID}"
fi
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
# 1b. UMI libraries (qc.umi_length > 0, e.g. QIAseq miRNA: 12):
#     read = <insert><adapter><UMI>...  ->  one read per (insert, UMI) molecule
# -----------------------------------------------------------------------------
# PCR copies of one molecule share insert AND UMI; counting reads instead of molecules
# inflates and distorts counts, most of all in low-input libraries such as EVs.
# Deduplication sorts on disk (low memory for ~100 M reads), keeping the first read of
# each (insert, UMI). Reads are then already trimmed, so miRTrace gets no adapter.
MT_ADAPTER="${MIRTRACE_ADAPTER}"
INPUT_READS=NA; UMI_MOLECULES=NA
if [ "${UMI_LENGTH}" -gt 0 ]; then
    log "UMI deduplication: adapter ${MIRTRACE_ADAPTER}, ${UMI_LENGTH}-nt UMI"
    DEDUP="${TMP}/${SAMPLE}.umi_dedup.fastq.gz"
    mkdir -p "${TMP}/sort"
    case "${FASTQ}" in *.gz) cat_fq=("${PIGZ}" -dc "${FASTQ}") ;; *) cat_fq=(cat "${FASTQ}") ;; esac
    "${cat_fq[@]}" \
      | awk -v A="${MIRTRACE_ADAPTER}" -v U="${UMI_LENGTH}" -v S="${TMP}/umi_counts.tsv" '
          NR % 4 == 2 { s = $0 }
          NR % 4 == 0 {
              n++
              p = index(s, A)
              if (p == 0)                                  { noad++;  next }
              if (p == 1)                                  { dimer++; next }   # no insert
              if (p - 1 + length(A) + U > length(s))       { short++; next }   # UMI cut off
              umi = substr(s, p + length(A), U)
              if (umi ~ /[^ACGT]/)                         { badumi++; next }
              print substr(s, 1, p - 1) "\t" umi "\t" substr($0, 1, p - 1)
              ok++
          }
          END {
              printf "input_reads\t%d\nno_adapter\t%d\nadapter_dimer\t%d\numi_incomplete\t%d\numi_with_N\t%d\nwith_umi\t%d\n", \
                     n, noad, dimer, short, badumi, ok > S
          }' \
      | LC_ALL=C sort -t$'\t' -k1,1 -k2,2 -u -S 2G -T "${TMP}/sort" --parallel="${THREADS}" \
      | awk -F'\t' -v id="${SAMPLE}" '{ printf "@%s_%d\n%s\n+\n%s\n", id, NR, $1, $3 }' \
      | "${PIGZ}" -p "${THREADS}" > "${DEDUP}"
    INPUT_READS=$(awk -F'\t' '$1 == "input_reads" { print $2 }' "${TMP}/umi_counts.tsv")
    WITH_UMI=$(awk -F'\t' '$1 == "with_umi" { print $2 }' "${TMP}/umi_counts.tsv")
    UMI_MOLECULES=$("${PIGZ}" -dc "${DEDUP}" | awk 'END { print NR / 4 }')
    [ "${UMI_MOLECULES}" -gt 0 ] || die "no reads with adapter and complete UMI - check qc.adapter and qc.umi_length"
    { cat "${TMP}/umi_counts.tsv"; printf 'unique_molecules\t%s\n' "${UMI_MOLECULES}"; } > "${SDIR}/umi_stats.tsv"
    log "UMI: ${INPUT_READS} reads, ${WITH_UMI} with adapter+UMI, ${UMI_MOLECULES} unique molecules ($(( 100 * UMI_MOLECULES / (WITH_UMI > 0 ? WITH_UMI : 1) ))% of them)"
    if [ "${WITH_UMI}" -lt $(( INPUT_READS / 2 )) ]; then
        log "WARNING: only ${WITH_UMI} of ${INPUT_READS} reads have the adapter and a complete UMI; see ${SDIR}/umi_stats.tsv"
    fi
    if [ "${FTYPE}" = "bam" ]; then rm -f "${FASTQ}"; fi
    FASTQ="${DEDUP}"
    MT_ADAPTER=""
fi

# -----------------------------------------------------------------------------
# 2. miRTrace QC -> uncollapsed FASTA of QC-passed reads
# -----------------------------------------------------------------------------
log "miRTrace QC"
mt_args=(qc --species "${MIRTRACE_SPECIES}" --output-dir "${SDIR}/mirtrace"
         --write-fasta --uncollapse-fasta --num-threads "${THREADS}" --force)
if [ -n "${MIRTRACE_PHRED}" ]; then
    # A PHRED offset can only be given through miRTrace's per-sample config CSV
    # (path,name,adapter,phred); in that mode the phred field is mandatory.
    printf '%s,%s,%s,%s\n' "${FASTQ}" "${SAMPLE}" "${MT_ADAPTER}" "${MIRTRACE_PHRED}" \
        > "${TMP}/mirtrace_config.csv"
    mt_args+=(--config "${TMP}/mirtrace_config.csv")
else
    if [ -n "${MT_ADAPTER}" ]; then mt_args+=(--adapter "${MT_ADAPTER}"); fi
    mt_args+=("${FASTQ}")
fi
# Java heap for the SURFR2 miRTrace launcher (container/install_tools.sh). The job's
# memory must exceed it by ~1 GB of JVM overhead; the validator checks slurm.sample.mem.
export MIRTRACE_HEAP_GB="${MIRTRACE_MEMORY_GB}"
"${MIRTRACE}" "${mt_args[@]}" > "${SDIR}/mirtrace.log" 2>&1 \
    || die "miRTrace failed (exit $?); see ${SDIR}/mirtrace.log"
if grep -q 'OutOfMemoryError' "${SDIR}/mirtrace.log"; then
    die "miRTrace ran out of Java heap (${MIRTRACE_HEAP_GB} GB): raise qc.mirtrace_memory_gb and slurm.sample.mem"
fi
if grep -q -E '^ERROR|Error parsing|aborting' "${SDIR}/mirtrace.log"; then
    die "miRTrace reported an error: $(grep -m1 -E '^ERROR|Error parsing|Could not' "${SDIR}/mirtrace.log")"
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

# Completion check. miRTrace writes its per-category read counts only after a run
# finishes; every input read is in exactly one category. A FASTA truncated by a killed
# run cannot match "total minus failing categories". (Comparing with KMC's read count
# alone cannot catch this: both read the same, possibly truncated, file.)
QCSTAT="${SDIR}/mirtrace/mirtrace-stats-qcstatus.tsv"
[ -s "${QCSTAT}" ] || die "miRTrace did not finish: ${QCSTAT} missing (killed? see ${SDIR}/mirtrace.log)"
read -r MT_TOTAL MT_FAILED < <(awk -F'\t' 'NR > 1 { t += $2; if ($1 ~ /^LOW_|SHORTER/) f += $2 } END { print t + 0, f + 0 }' "${QCSTAT}")
[ "${QC_READS}" -eq $(( MT_TOTAL - MT_FAILED )) ] \
    || die "QC-passed FASTA has ${QC_READS} reads but miRTrace reports $(( MT_TOTAL - MT_FAILED )) passing of ${MT_TOTAL}: incomplete output"
log "QC-passed reads: ${QC_READS} of ${MT_TOTAL} ($(( 100 * QC_READS / MT_TOTAL ))%)"

# -----------------------------------------------------------------------------
# 2b. Artifact removal: drop QC-passed reads containing ANY k-mer of an artifact
#     sequence (Illumina adapters/primers/indexes, PhiX, the configured adapter,
#     extra FASTA files), on both strands.
# -----------------------------------------------------------------------------
# Whole reads are removed, not only artifact k-mers: a partly-artifact read also yields
# junction k-mers that match no reference, and removing the read keeps the library size
# consistent with the counts. Matching uses the analysis k, so fragments shorter than k
# cannot be detected; a shorter k would delete genuine small RNAs that happen to share
# a short word with PhiX or an adapter, in every sample.
ARTIFACT_READS=NA
if [ "${ARTIFACTS_ENABLED}" = "true" ]; then
    REF="${TMP}/artifact_refs.fa"
    {
        for f in "${ART_FILES[@]}"; do
            case "${f}" in *.gz) "${PIGZ}" -dc "${f}" ;; *) cat "${f}" ;; esac
            echo
        done
        if [ -n "${MIRTRACE_ADAPTER}" ]; then printf '>configured_adapter\n%s\n' "${MIRTRACE_ADAPTER}"; fi
    } | awk '
        function emit() {
            if (seq == "") return
            rc = ""
            for (i = length(seq); i > 0; i--) {
                c = substr(seq, i, 1)
                rc = rc (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : "N")
            }
            n++; print ">ref" n "_fwd"; print seq
            print ">ref" n "_rc"; print rc
        }
        /^>/ { emit(); seq = ""; next }
        { l = toupper($0); gsub(/[ \t\r]/, "", l); gsub(/U/, "T", l); seq = seq l }   # RNA -> DNA first
        END { emit() }' > "${REF}"
    N_REF=$(( $(grep -c '^>' "${REF}") / 2 ))
    mkdir -p "${TMP}/kmc_art_tmp"
    art_args=(-hp -fm "-k${K}" -ci1 -cs255 -m2 "-t${THREADS}")
    if [ "${CANONICAL}" = "false" ]; then art_args+=(-b); fi
    "${KMC}" "${art_args[@]}" "${REF}" "${TMP}/artifacts" "${TMP}/kmc_art_tmp" > "${SDIR}/kmc_artifacts.log" 2>&1 \
        || die "building the artifact k-mer database failed; see ${SDIR}/kmc_artifacts.log"
    # keep reads with 0 artifact k-mers (-ci0 -cx0 on the read set)
    "${KMC_TOOLS}" -hp -t"${THREADS}" filter "${TMP}/artifacts" "${FASTA}" -fa -ci0 -cx0 "${TMP}/clean.fa" \
        >> "${SDIR}/kmc_artifacts.log" 2>&1 || die "artifact filtering failed; see ${SDIR}/kmc_artifacts.log"
    CLEAN_READS=$(awk 'substr($0, 1, 1) == ">" { n++ } END { print n + 0 }' "${TMP}/clean.fa")
    [ "${CLEAN_READS}" -gt 0 ] || die "every QC-passed read matched an artifact sequence - check the references"
    ARTIFACT_READS=$(( QC_READS - CLEAN_READS ))
    "${PIGZ}" -p "${THREADS}" -c "${TMP}/clean.fa" > "${FASTA}"
    rm -f "${TMP}/clean.fa"
    printf 'qc_passed_reads\t%s\nartifact_reads\t%s\nclean_reads\t%s\nreference_sequences\t%s\nreference_id\t%s\n' \
        "${QC_READS}" "${ARTIFACT_READS}" "${CLEAN_READS}" "${N_REF}" "${ART_ID}" > "${SDIR}/artifact_stats.tsv"
    log "artifacts: removed ${ARTIFACT_READS} of ${QC_READS} QC-passed reads ($(awk -v a="${ARTIFACT_READS}" -v q="${QC_READS}" 'BEGIN { printf "%.2f", 100 * a / q }')%) matching ${N_REF} reference sequences"
    QC_READS=${CLEAN_READS}
fi

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
    printf 'sample_id\tcohort\tcondition\trole\tqc_reads\ttotal_kmers\tunique_kmers\tinput_reads\tumi_molecules\tartifact_reads\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${SAMPLE}" "${COHORT}" "${CONDITION}" "${ROLE}" \
        "${QC_READS}" "${TOTAL_KMERS}" "${UNIQUE_KMERS}" "${INPUT_READS}" "${UMI_MOLECULES}" "${ARTIFACT_READS}"
} > "${SDIR}/library_size.tsv.part"
mv "${SDIR}/library_size.tsv.part" "${SDIR}/library_size.tsv"

if [ "${KEEP_INTERMEDIATES}" = "true" ]; then
    mv "${FASTA}" "${SDIR}/"
    if [ "${FTYPE}" = "bam" ]; then mv "${FASTQ}" "${SDIR}/"; fi
fi

printf '%s\n' "${FINGERPRINT}" > "${SDIR}/.done.part" && mv "${SDIR}/.done.part" "${SDIR}/.done"
log "done: total_kmers=${TOTAL_KMERS} unique_kmers=${UNIQUE_KMERS}"
