#!/usr/bin/env bash
# SURFR2 regression tests on synthetic data with planted truth.
#   bash tests/run_test.sh [scenario ...] [-- launcher options, e.g. --container SANDBOX]
# Scenarios: multi multi_descriptive small single null (default: all). Needs samtools + python3-numpy.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
scenarios=() launcher_opts=()
while [ "$#" -gt 0 ]; do
    if [ "$1" = "--" ]; then shift; launcher_opts=("$@"); break; fi
    scenarios+=("$1"); shift
done
[ "${#scenarios[@]}" -gt 0 ] || scenarios=(multi multi_descriptive small single null)

WORK=$(mktemp -d "${TMPDIR:-/tmp}/surfr2_test.XXXXXX")
status=0
for s in "${scenarios[@]}"; do
    echo "=== scenario: ${s}"
    mkdir -p "${WORK}/${s}" && cd "${WORK}/${s}"
    python3 "${HERE}/make_synthetic.py" "${s%%_*}"    # multi_descriptive reuses the multi data
    cp "${HERE}/configs/${s}.yaml" config.yaml
    if bash "${HERE}/../run_SURFR2.sh" -c config.yaml "${launcher_opts[@]}" > launcher.log 2>&1; then
        python3 "${HERE}/check_results.py" "${s}" out || status=1
    else
        echo "FAIL pipeline did not complete - see ${WORK}/${s}/launcher.log"; status=1
    fi
done
if [ "${status}" -eq 0 ]; then echo "ALL TESTS PASSED (work dir: ${WORK})"; else echo "TESTS FAILED (work dir: ${WORK})"; fi
exit "${status}"
