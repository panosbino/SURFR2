#!/usr/bin/env bash
# SURFR2 regression test: bash tests/run_test.sh [--container SANDBOX]
# Runs the full pipeline locally on synthetic data and checks the planted truth.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/surfr2_test.XXXXXX")
cp "${HERE}/make_synthetic.py" "${HERE}/test_config.yaml" "${WORK}/"
cd "${WORK}"
python3 make_synthetic.py
bash "${HERE}/../run_SURFR2.sh" -c test_config.yaml "$@"
fail=0
while read -r name seq; do
    found=$(awk -F'\t' -v s="${seq}" 'NR > 1 && $2 == s' out/results/case_specific_sequences.tsv | wc -l)
    case "${name}" in S1_*|S5_*) want=1 ;; *) want=0 ;; esac
    if [ "${found}" -eq "${want}" ]; then echo "PASS ${name}"; else echo "FAIL ${name} (found ${found}, want ${want})"; fail=1; fi
done < truth.tsv
n=$(( $(wc -l < out/results/case_specific_sequences.tsv) - 1 ))
if [ "${n}" -eq 2 ]; then echo "PASS exactly 2 sequences"; else echo "FAIL ${n} sequences, want 2"; fail=1; fi
[ "${fail}" -eq 0 ] && echo "ALL TESTS PASSED (work dir: ${WORK})" || { echo "TESTS FAILED (work dir: ${WORK})"; exit 1; }
