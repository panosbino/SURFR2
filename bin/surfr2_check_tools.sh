#!/usr/bin/env bash
# =============================================================================
# SURFR2 - check that every tool and R package resolves in the configured
# environment, and record exactly which versions a run used.
#
# Run by the launcher before anything is submitted, inside the same environment
# as the steps (container, or modules + path_prepend + r_libs).
#
# Usage:  surfr2_check_tools.sh <run_config_dir>
# Writes: <run_config_dir>/tool_versions.txt   (copied to results/ by the filter step)
# Exit 1 if a required tool or package is missing; optional ones only warn.
# =============================================================================

set -uo pipefail          # no -e: collect every problem, then report them together
export LC_ALL=C

[ "$#" -eq 1 ] || { echo "usage: $0 <run_config_dir>" >&2; exit 2; }
RUN_CONFIG_DIR=$1
# shellcheck source=/dev/null
source "${RUN_CONFIG_DIR}/params.env"
OUT="${RUN_CONFIG_DIR}/tool_versions.txt"

missing=() warnings=() rows=()

resolve() { command -v "$1" 2>/dev/null || true; }
sha() { sha256sum "$1" 2>/dev/null | cut -c1-16; }

# check <label> <command> <required:yes|no> <version command...>
check() {
    local label=$1 cmd=$2 required=$3; shift 3
    local path version
    path=$(resolve "${cmd}")
    if [ -z "${path}" ]; then
        if [ "${required}" = yes ]; then missing+=("${label} (${cmd})"); else warnings+=("${label} (${cmd}) not found"); fi
        rows+=("${label}|MISSING|${cmd}|-|-")
        return
    fi
    version=$("$@" 2>&1 | grep -v '^[[:space:]]*$' | head -n 1 || true)
    rows+=("${label}|ok|${path}|${version:-unknown}|$(sha "${path}")")
}

# A container image must have a complete root filesystem. A .sif converted with its temp
# directory on Lustre came out with only /etc/resolv.conf: tools that live outside /etc
# still resolve, so check for it explicitly before checking tools.
if [ "${ENVIRONMENT}" = "container" ] && { [ ! -e /etc/passwd ] || [ ! -d /etc/alternatives ]; }; then
    missing+=("complete image: /etc in the container has $(find /etc -mindepth 1 -maxdepth 1 2>/dev/null | wc -l) entries (expected ~100) - the .sif conversion lost files; rebuild it (README: Building the container)")
fi

check samtools  "${SAMTOOLS}"  yes "${SAMTOOLS}" --version
check pigz      "${PIGZ}"      yes "${PIGZ}" --version
check python3   python3        yes python3 --version
check kmc       "${KMC}"       yes "${KMC}"
check kmc_tools "${KMC_TOOLS}" yes "${KMC_TOOLS}"
check mirtrace  "${MIRTRACE}"  yes "${MIRTRACE}" --version
check java      java           yes java -version

# miRTrace must be SURFR2's launcher, not the upstream Python wrapper: the wrapper sizes
# Java's heap from the whole node's RAM (OOM kills on shared nodes) and discards Java's
# exit status (a killed run looks successful).
mt_path=$(resolve "${MIRTRACE}")
if [ -n "${mt_path}" ] && ! grep -q 'surfr2-mirtrace-launcher' "${mt_path}" 2>/dev/null; then
    missing+=("SURFR2 miRTrace launcher: '${mt_path}' is not it (upstream wrapper?) - rerun container/install_tools.sh, or rebuild the image")
fi
check mergeTags "${MERGETAGS}"  yes echo "no version flag; identified by sha256"
check Rscript   "${RSCRIPT}"   yes "${RSCRIPT}" --version
check pandoc    pandoc         no  pandoc --version

# KMC and kmc_tools have no version flag; their usage header (first line) carries it.

# ---- artifact references (artifacts.builtin) --------------------------------------
if [ "${ARTIFACTS_ENABLED}" = "true" ] && [ "${ARTIFACTS_BUILTIN}" = "true" ]; then
    kp=$(resolve "${KMC}")
    if [ -n "${kp}" ]; then
        ad="$(cd "$(dirname "$(readlink -f "${kp}")")/.." && pwd)/share/surfr2/artifacts"
        for f in illumina_adapters.fa phix174.fa.gz; do
            if [ -s "${ad}/${f}" ]; then rows+=("artifacts|ok|${ad}/${f}|reference|$(sha "${ad}/${f}")")
            else missing+=("artifact reference ${ad}/${f} - rerun container/install_tools.sh, or set artifacts.builtin: false"); fi
        done
    fi
fi

# ---- R packages ---------------------------------------------------------------
required_r="jsonlite readr dplyr ggplot2"
[ "${NEEDS_EDGER}" = "true" ] && required_r="${required_r} edgeR"
optional_r="plotly htmlwidgets ggvenn"
r_report=""
if [ -n "$(resolve "${RSCRIPT}")" ]; then
    r_report=$("${RSCRIPT}" --vanilla -e "
      req <- strsplit('${required_r}', ' ')[[1]]; opt <- strsplit('${optional_r}', ' ')[[1]]
      cat('R', as.character(getRversion()), '| libraries:', paste(.libPaths(), collapse = ' : '), '\n')
      for (p in c(req, opt)) {
        ok <- requireNamespace(p, quietly = TRUE)
        cat(sprintf('%s|%s|%s\n', p, if (ok) as.character(packageVersion(p)) else 'MISSING',
                    if (p %in% req) 'required' else 'optional'))
      }" 2>&1)
    # R must actually run: if it crashes (e.g. a missing shared library) there are no
    # package lines, and nothing would otherwise be flagged.
    if ! grep -q '^R [0-9]' <<< "${r_report}"; then
        missing+=("R does not run: $(grep -m1 -i -E 'error|cannot' <<< "${r_report}" || echo "${r_report}" | head -n 1)")
    fi
    while IFS='|' read -r pkg ver kind; do
        [ "${ver:-}" = "MISSING" ] || continue
        if [ "${kind}" = required ]; then missing+=("R package ${pkg}"); else warnings+=("R package ${pkg} not installed (${kind})"); fi
    done <<< "$(grep '|' <<< "${r_report}")"
fi

# ---- write the record ---------------------------------------------------------
{
    echo "# SURFR2 tool versions - $(date '+%F %T') on $(hostname)"
    echo "# environment: ${ENVIRONMENT}${CONTAINER:+ (${CONTAINER}, $(stat -c '%s bytes, modified %y' "${CONTAINER}" 2>/dev/null || echo 'not readable here'))}"
    [ -n "${MODULES}" ] && echo "# modules requested: ${MODULES}"
    # LOADEDMODULES is set by Lmod and passed into containers; calling 'module list' inside
    # a container fails, because the host's Lmod is not there.
    if [ -n "${LOADEDMODULES:-}" ]; then
        echo "# modules loaded: ${LOADEDMODULES//:/ }"
    fi
    [ -n "${PATH_PREPEND}" ] && echo "# path_prepend: ${PATH_PREPEND}"
    [ -n "${R_LIBS_DIR}" ] && echo "# r_libs: ${R_LIBS_DIR}"
    echo
    printf '%-10s %-8s %-48s %-16s %s\n' tool status path sha256[:16] version
    for r in "${rows[@]}"; do
        IFS='|' read -r t st p v h <<< "${r}"
        printf '%-10s %-8s %-48s %-16s %s\n' "${t}" "${st}" "${p}" "${h}" "${v}"
    done
    echo
    [ -n "${r_report}" ] && echo "${r_report}" | sed 's/|/  /g'
} > "${OUT}"
cat "${OUT}"

if [ "${#warnings[@]}" -gt 0 ]; then
    printf 'WARNING: %s\n' "${warnings[@]}" >&2
fi
if [ "${#missing[@]}" -gt 0 ]; then
    echo >&2
    echo "ERROR: missing in the '${ENVIRONMENT}' environment:" >&2
    printf '  - %s\n' "${missing[@]}" >&2
    if [ "${ENVIRONMENT}" = "modules" ]; then
        cat >&2 <<'EOF'
Fix one of three ways (see README, "Development environment"):
  - load a module that provides it:  module spider <name>, then add it to execution.modules
  - install kmc, kmc_tools, mirtrace, mergeTags, pigz into a prefix:
        bash container/install_tools.sh <prefix> [--with-pigz]
    and add <prefix>/bin to execution.path_prepend
  - R packages: Rscript container/install_r_packages.R <dir>, then set execution.r_libs: <dir>
  - or give an absolute binary path in tools.<name>
EOF
    else
        echo "The container image lacks these; rebuild it from container/Dockerfile." >&2
    fi
    exit 1
fi
