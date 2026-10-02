#!/usr/bin/env bash
# =============================================================================
# SURFR2 - identify the small-RNA library structure from raw reads
#
# Usage:  surfr2_detect_adapter.sh <reads.fastq[.gz]> [n_reads (default 200000)]
#         ... | surfr2_detect_adapter.sh - [n_reads]
#
# Reports, from the first n reads:
#  - read length (untrimmed reads are much longer than ~22-nt miRNAs)
#  - how many reads contain each known 3' adapter, and where
#  - the insert length (sequence before the adapter): miRNAs give a peak at ~22 nt;
#    a peak at ~30 nt means 4 random bases on each side (NEXTflex)
#  - base diversity at the first read positions: random bases (NEXTflex 4N) are
#    near-uniform, real RNA 5' ends are not
#  - for QIAseq: whether the 12 nt after the adapter look random (UMI)
# and a recommendation for SURFR2's qc.adapter setting.
# =============================================================================

set -euo pipefail
[ "$#" -ge 1 ] || { echo "usage: $0 <reads.fastq[.gz] | -> [n_reads]" >&2; exit 2; }
IN=$1; N=${2:-200000}

reader() {
    case "${IN}" in
        -) cat ;;
        *.gz) gzip -dc "${IN}" ;;
        *) cat "${IN}" ;;
    esac
}

# head closes the pipe early; that is expected, so don't let pipefail abort on it
{ reader 2>/dev/null || true; } | head -n $(( N * 4 )) | awk -v N="${N}" '
BEGIN {
    # name, 3-prime adapter, structure note
    ad["TruSeq/NEXTflex 3prime"] = "TGGAATTCTCGGGTGCCAAGG"
    ad["QIAseq miRNA"]           = "AACTGTAGGCACCATCAAT"
    ad["NEBNext/Illumina univ."] = "AGATCGGAAGAGCACACGTCT"
    PL = 12                       # prefix length used for matching
}
NR % 4 == 2 {
    n++; len = length($0); lensum += len; if (len > lenmax) lenmax = len
    for (a in ad) {
        p = index($0, substr(ad[a], 1, PL))
        if (p > 0) {
            hit[a]++; ins[a, p - 1]++
            if (a == "QIAseq miRNA" && p - 1 + length(ad[a]) + 12 <= len) {
                u = substr($0, p + length(ad[a]), 12); umi[u]++; numi++
            }
        }
    }
    for (i = 1; i <= 6 && i <= len; i++) base[i, substr($0, i, 1)]++
}
END {
    if (n == 0) { print "no reads read"; exit 1 }
    printf "reads examined: %d   mean length: %.1f   max: %d\n\n", n, lensum / n, lenmax
    best = ""; bestn = 0
    printf "%-24s %9s   %s\n", "3prime adapter", "reads", "most common insert lengths (nt:% of hits)"
    for (a in ad) {
        h = hit[a] + 0
        line = ""
        if (h > 0) {
            # top 4 insert lengths
            for (k = 1; k <= 4; k++) {
                m = -1; ml = -1
                for (L = 0; L <= lenmax; L++) if ((a, L) in ins && ins[a, L] > m && !((a, L) in used)) { m = ins[a, L]; ml = L }
                if (ml < 0) break
                used[a, ml] = 1
                line = line sprintf("%d:%.0f%%  ", ml, 100 * m / h)
                if (k == 1) mode[a] = ml
            }
        }
        printf "%-24s %8.1f%%   %s\n", a, 100 * h / n, line
        if (h > bestn) { bestn = h; best = a }
    }
    printf "\nbase composition at read positions 1-6 (A/C/G/T %%):\n"
    for (i = 1; i <= 6; i++) {
        t = base[i,"A"] + base[i,"C"] + base[i,"G"] + base[i,"T"]
        if (t == 0) continue
        mx = 0; for (b in arr) delete arr[b]
        split("A C G T", B, " ")
        s = ""; for (j = 1; j <= 4; j++) { f = 100 * base[i, B[j]] / t; s = s sprintf("%5.1f", f); if (f > mx) mx = f }
        printf "  pos %d: %s   (max %.0f%%)\n", i, s, mx
        if (i <= 4) maxsum += mx
    }
    random4 = (maxsum / 4 < 35)      # all four bases ~25%: random adapter bases

    if (numi > 100) {
        distinct = 0; for (u in umi) distinct++
        printf "\nQIAseq UMI check: %d distinct 12-mers after the adapter in %d reads (%.0f%%)\n", distinct, numi, 100 * distinct / numi
    }

    printf "\n==> "
    if (bestn < 0.3 * n && lensum / n <= 35) {
        printf "reads are already adapter-trimmed (mean length %.1f nt): set qc.adapter: null\n", lensum / n
    } else if (bestn < 0.3 * n) {
        printf "no known adapter in most reads (best: %s, %.0f%%) although reads are long (mean %.0f nt).\n", (best == "" ? "none" : best), 100 * bestn / n, lensum / n
        printf "    Another library kit: the adapter must be identified before running SURFR2.\n"
    } else if (best == "QIAseq miRNA") {
        printf "QIAseq miRNA library: set qc.adapter: %s\n", ad[best]
        printf "    Reads carry 12-nt UMIs after the adapter. SURFR2 does not yet deduplicate UMIs,\n"
        printf "    so counts include PCR duplicates.\n"
    } else if (random4 || (mode[best] >= 28 && mode[best] <= 32)) {
        printf "NEXTflex-type library with 4 random bases on each side of the insert (insert peak %d nt%s).\n", mode[best], random4 ? ", random first bases" : ""
        printf "    SURFR2 cannot handle this correctly yet: the random bases must be removed before k-mer\n"
        printf "    counting, which miRTrace does not do. Do not run SURFR2 on these reads as they are.\n"
    } else {
        printf "%s adapter, insert peak %d nt: set qc.adapter: %s\n", best, mode[best], ad[best]
    }
}'
