#!/usr/bin/env bash
# =============================================================================
# SURFR2 - install the tools that are not available as system packages or modules
#
# The single source of truth for these tool versions: the Dockerfile runs this
# script, and so does a development setup on a cluster. Both therefore get
# identical binaries.
#
#   bash install_tools.sh <prefix> [--with-pigz]
#
# Installs into <prefix>/bin:
#   kmc, kmc_tools   KMC 3.2.4, official statically linked release binaries
#   mirtrace         miRTrace 1.0.1 launcher (jar in <prefix>/share/mirtrace); see below
#   mergeTags        dekupl-mergeTags, built from a pinned commit
#   pigz             only with --with-pigz (most systems provide it)
#
# Every download is verified against a SHA-256 checksum; a mismatch aborts.
# Requirements: bash, curl, tar, unzip, gcc, make, zlib headers; Java >= 8 at
# run time for miRTrace.
# =============================================================================

set -euo pipefail

# ---- pinned versions (change here, and only here) ---------------------------
KMC_VERSION=3.2.4
KMC_URL="https://github.com/refresh-bio/KMC/releases/download/v${KMC_VERSION}/KMC${KMC_VERSION}.linux.x64.tar.gz"
KMC_SHA256=158f2084f8d928b3f33b8aaf7d1220fee4183bf46837787e5e6b16bbdf54d31d

MIRTRACE_VERSION=1.0.1
MIRTRACE_URL="https://github.com/friedlanderlab/mirtrace/releases/download/v${MIRTRACE_VERSION}/mirtrace-v${MIRTRACE_VERSION}.zip"
MIRTRACE_SHA256=952e9b07d7a16ee475652683e780ce3e6b9f8261a75c2154893f25ae5177b8cd

# dekupl-mergeTags has no releases; pin the commit (SURFR1 cloned 'master' unpinned)
MERGETAGS_COMMIT=4cdad2c5ce45c3a30458aa73ce970e31c7646699
MERGETAGS_URL="https://codeload.github.com/Transipedia/dekupl-mergeTags/tar.gz/${MERGETAGS_COMMIT}"

PIGZ_VERSION=2.8
PIGZ_URL="https://github.com/madler/pigz/archive/refs/tags/v${PIGZ_VERSION}.tar.gz"
PIGZ_SHA256=2f7f6a6986996d21cb8658535fff95f1c7107ddce22b5324f4b41890e2904706
# -----------------------------------------------------------------------------

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[install_tools] $*" >&2; }

[ "$#" -ge 1 ] || die "usage: $0 <prefix> [--with-pigz]"
PREFIX=$(mkdir -p "$1" && cd "$1" && pwd)
WITH_PIGZ=false
[ "${2:-}" = "--with-pigz" ] && WITH_PIGZ=true

for c in curl tar unzip gcc make sha256sum; do
    command -v "$c" > /dev/null || die "'$c' is required"
done

WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT
mkdir -p "${PREFIX}/bin" "${PREFIX}/share"

fetch() {   # fetch <url> <dest> [sha256]
    curl -fsSL --retry 3 -o "$2" "$1" || die "download failed: $1"
    if [ -n "${3:-}" ]; then
        echo "$3  $2" | sha256sum -c --quiet - || die "checksum mismatch for $1"
    fi
}

# ---- KMC ---------------------------------------------------------------------
log "KMC ${KMC_VERSION}"
fetch "${KMC_URL}" "${WORK}/kmc.tgz" "${KMC_SHA256}"
mkdir -p "${WORK}/kmc" && tar -xzf "${WORK}/kmc.tgz" -C "${WORK}/kmc"
install -m 755 "${WORK}/kmc/bin/kmc" "${WORK}/kmc/bin/kmc_tools" "${PREFIX}/bin/"

# ---- miRTrace ----------------------------------------------------------------
log "miRTrace ${MIRTRACE_VERSION}"
fetch "${MIRTRACE_URL}" "${WORK}/mirtrace.zip" "${MIRTRACE_SHA256}"
unzip -q "${WORK}/mirtrace.zip" -d "${WORK}/mirtrace"
src=$(find "${WORK}/mirtrace" -name mirtrace.jar -printf '%h\n' | head -n 1)
[ -n "${src}" ] || die "mirtrace.jar not found in the release archive"
rm -rf "${PREFIX}/share/mirtrace" && mkdir -p "${PREFIX}/share/mirtrace"
cp "${src}/mirtrace.jar" "${PREFIX}/share/mirtrace/"
# SURFR2 launches the jar itself instead of miRTrace's Python wrapper. The wrapper sizes
# the Java heap at half of the WHOLE NODE's physical RAM (not the job's allocation),
# which gets jobs OOM-killed on shared nodes, and it discards Java's exit status, so a
# killed run looks successful. The heap comes from MIRTRACE_HEAP_GB (set by SURFR2 from
# qc.mirtrace_memory_gb). -Xms must equal -Xmx: miRTrace sizes its read hash table from
# the heap committed at start-up (Runtime.totalMemory()) and grows into all of it.
# 'exec' passes Java's exit status through unchanged.
cat > "${PREFIX}/bin/mirtrace" <<EOS
#!/bin/sh
# surfr2-mirtrace-launcher (marker checked by surfr2_check_tools.sh)
exec java -Xms"\${MIRTRACE_HEAP_GB:-4}g" -Xmx"\${MIRTRACE_HEAP_GB:-4}g" -jar "${PREFIX}/share/mirtrace/mirtrace.jar" "\$@"
EOS
chmod 755 "${PREFIX}/bin/mirtrace"

# ---- dekupl-mergeTags ----------------------------------------------------------
log "dekupl-mergeTags ${MERGETAGS_COMMIT:0:8}"
fetch "${MERGETAGS_URL}" "${WORK}/mergetags.tgz"
mkdir -p "${WORK}/mergetags" && tar -xzf "${WORK}/mergetags.tgz" -C "${WORK}/mergetags" --strip-components=1
make -C "${WORK}/mergetags" -s > "${WORK}/mergetags.log" 2>&1 || { cat "${WORK}/mergetags.log" >&2; die "mergeTags build failed (needs gcc and zlib headers)"; }
install -m 755 "${WORK}/mergetags/mergeTags" "${PREFIX}/bin/"

# ---- pigz (optional) ---------------------------------------------------------
if [ "${WITH_PIGZ}" = true ]; then
    log "pigz ${PIGZ_VERSION}"
    fetch "${PIGZ_URL}" "${WORK}/pigz.tgz" "${PIGZ_SHA256}"
    mkdir -p "${WORK}/pigz" && tar -xzf "${WORK}/pigz.tgz" -C "${WORK}/pigz" --strip-components=1
    make -C "${WORK}/pigz" -s > "${WORK}/pigz.log" 2>&1 || { cat "${WORK}/pigz.log" >&2; die "pigz build failed (needs zlib headers)"; }
    install -m 755 "${WORK}/pigz/pigz" "${PREFIX}/bin/"
fi

# ---- smoke test ----------------------------------------------------------------
# Both print usage and exit non-zero without arguments: capture first, then match
# (a pipe into grep would fail under 'pipefail' even when the text matches).
out=$("${PREFIX}/bin/kmc" 2>&1 || true);       [[ "${out}" == *"ver. ${KMC_VERSION}"* ]] || die "kmc smoke test failed"
out=$("${PREFIX}/bin/mergeTags" 2>&1 || true); [[ "${out}" == *"Usage"* ]] || die "mergeTags smoke test failed"
out=$("${PREFIX}/bin/mirtrace" --help 2>&1 || true); [[ "${out}" == *"miRTrace"* ]] || die "mirtrace smoke test failed (java on PATH?)"
if [ "${WITH_PIGZ}" = true ]; then "${PREFIX}/bin/pigz" --version > /dev/null 2>&1 || die "pigz smoke test failed"; fi
log "installed into ${PREFIX}/bin"
ls -1 "${PREFIX}/bin" >&2
