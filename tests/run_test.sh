#!/usr/bin/env bash
# SURFR2 regression tests on synthetic data with planted truth.
#   bash tests/run_test.sh [scenario ...] [-- launcher options, e.g. --container SANDBOX]
# Scenarios: multi multi_descriptive small single null small_dea single_dea null_dea umi (default: all). Needs python3-numpy.
#
# Environment for the pipeline steps (default: modules, i.e. tools from the login PATH):
#   SURFR2_TEST_ENV=container SURFR2_TEST_CONTAINER=/path/surfr2.sif   use the image
#   SURFR2_TEST_MODULES="PDC samtools R"    modules to load      (either environment)
#   SURFR2_TEST_TOOLS=/path/prefix/bin      added to PATH        (modules environment)
#   SURFR2_TEST_RLIBS=/path/R/library       personal R library   (modules environment)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
scenarios=() launcher_opts=()
while [ "$#" -gt 0 ]; do
    if [ "$1" = "--" ]; then shift; launcher_opts=("$@"); break; fi
    scenarios+=("$1"); shift
done
[ "${#scenarios[@]}" -gt 0 ] || scenarios=(multi multi_descriptive small single null small_dea single_dea null_dea umi)

WORK=$(mktemp -d "${TMPDIR:-/tmp}/surfr2_test.XXXXXX")
status=0
for s in "${scenarios[@]}"; do
    echo "=== scenario: ${s}"
    mkdir -p "${WORK}/${s}" && cd "${WORK}/${s}"
    python3 "${HERE}/make_synthetic.py" "${s%%_*}"    # multi_descriptive reuses the multi data
    # replace the scenario's execution line with the requested test environment
    exe="executor: local, environment: ${SURFR2_TEST_ENV:-modules}"
    [ -z "${SURFR2_TEST_MODULES:-}" ] || exe="${exe}, modules: [$(echo "${SURFR2_TEST_MODULES}" | sed 's/ \+/, /g')]"
    if [ "${SURFR2_TEST_ENV:-modules}" = container ]; then
        exe="${exe}, container: ${SURFR2_TEST_CONTAINER:?set SURFR2_TEST_CONTAINER}"
    else
        [ -z "${SURFR2_TEST_TOOLS:-}" ]   || exe="${exe}, path_prepend: [${SURFR2_TEST_TOOLS}]"
        [ -z "${SURFR2_TEST_RLIBS:-}" ]   || exe="${exe}, r_libs: ${SURFR2_TEST_RLIBS}"
    fi
    sed "s|^execution: .*|execution: {${exe}}|" "${HERE}/configs/${s}.yaml" > config.yaml
    if bash "${HERE}/../run_SURFR2.sh" -c config.yaml "${launcher_opts[@]}" > launcher.log 2>&1; then
        python3 "${HERE}/check_results.py" "${s}" out || status=1
    else
        echo "FAIL pipeline did not complete - see ${WORK}/${s}/launcher.log"; status=1
    fi
done
if [ "${status}" -eq 0 ]; then echo "ALL TESTS PASSED (work dir: ${WORK})"; else echo "TESTS FAILED (work dir: ${WORK})"; fi
exit "${status}"
