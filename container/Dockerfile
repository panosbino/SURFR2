# =============================================================================
# SURFR2 container (linux/amd64)
#
# Build (from the SURFR2 repository root):
#   docker buildx build --platform linux/amd64 -f container/Dockerfile \
#     -t surfr2:<version> --load container/
# Convert to Singularity on the build machine (NOT with a temp/cache dir on Lustre: on
# Dardel that silently produced an image with an almost empty /etc):
#   docker save surfr2:<version> -o surfr2_<version>.tar
#   docker run --rm --privileged -v "$(pwd)":/work quay.io/singularity/singularity:v4.1.0 \
#     build /work/surfr2_<version>.sif docker-archive:///work/surfr2_<version>.tar
#
# Only the TOOLS live in this image. SURFR2's own scripts (bin/, app/) are run
# from the repository through the image, so changing SURFR2 code never requires
# a rebuild - only adding or upgrading a tool or R package does.
#
# Tool versions are defined once, in install_tools.sh and install_r_packages.R,
# which a module-based development setup runs as well.
#
#   R           4.4.3 (rocker/r-ver: CRAN packages from a dated snapshot)
#   Bioconductor matching R 4.4 (edgeR)
#   samtools    1.23.1 (source)
#   KMC         3.2.4  (official static binaries, checksum-verified)
#   miRTrace    1.0.1  (release zip, checksum-verified) on OpenJDK 17
#   mergeTags   pinned commit (source)
#   pigz, pandoc, python3 + PyYAML from Ubuntu
# =============================================================================

# --platform: Dardel is x86-64, and many builds happen on Apple Silicon laptops
# hadolint ignore=DL3029
FROM --platform=linux/amd64 rocker/r-ver:4.4.3

# Every RUN below uses bash with pipefail, so a failing command inside a pipe
# (e.g. 'samtools --version | head') fails the build instead of being masked.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

LABEL org.opencontainers.image.title="SURFR2" \
      org.opencontainers.image.description="Tools for the SURFR2 small-RNA k-mer pipeline" \
      surfr2.r="4.4.3" surfr2.samtools="1.23.1" surfr2.kmc="3.2.4" surfr2.mirtrace="1.0.1"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    JAVA_TOOL_OPTIONS="-Djava.awt.headless=true" \
    PATH=/opt/surfr2/bin:/usr/local/bin:$PATH

# ---------------------------------------------------------------------------
# 1. System packages: build tools and libraries for samtools, mergeTags and
#    R packages; runtime tools (Java for miRTrace, pandoc for HTML reports).
#    Not pinned individually: Ubuntu withdraws superseded versions, which would
#    break rebuilds; the base image pins the distribution release.
# ---------------------------------------------------------------------------
# hadolint ignore=DL3008
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        ca-certificates curl unzip bzip2 make gcc g++ gfortran \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
        libxml2-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev \
        libfribidi-dev libpng-dev libjpeg-dev libtiff-dev \
        pigz pandoc python3 python3-yaml openjdk-17-jre-headless && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 2. samtools from source
# ---------------------------------------------------------------------------
ARG SAMTOOLS_VERSION=1.23.1
# hadolint ignore=DL3003
RUN curl -fsSL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
        -o /tmp/samtools.tar.bz2 && \
    tar -xjf /tmp/samtools.tar.bz2 -C /tmp && \
    cd /tmp/samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --without-curses && \
    make -j"$(nproc)" && make install && \
    cd / && rm -rf /tmp/samtools*

# ---------------------------------------------------------------------------
# 3. KMC, miRTrace, mergeTags (shared installer) and R packages (shared installer)
# ---------------------------------------------------------------------------
COPY install_tools.sh install_r_packages.R /tmp/surfr2/
RUN bash /tmp/surfr2/install_tools.sh /opt/surfr2 && \
    Rscript /tmp/surfr2/install_r_packages.R && \
    rm -rf /tmp/surfr2

# ---------------------------------------------------------------------------
# 4. No dependence on /etc/alternatives. Ubuntu reaches java, BLAS and LAPACK through
#    symlinks into /etc/alternatives. Under Singularity on some HPC systems (seen on
#    Dardel) the image's own /etc/alternatives is not what the container sees, so java
#    was "not found" and R failed to load libblas.so.3. Point every such link directly
#    at its final target.
# ---------------------------------------------------------------------------
RUN find / -xdev -type l -lname '/etc/alternatives/*' -print 2>/dev/null | \
    while read -r l; do t=$(readlink -f "$l") && [ -e "$t" ] && ln -sfn "$t" "$l"; done

# ---------------------------------------------------------------------------
# 5. Smoke tests: the build fails here if any tool is broken. They run with
#    /etc/alternatives moved away, reproducing the HPC situation above.
#    Every check fails the build on error: no "echo $(cmd)", which would hide failures.
# ---------------------------------------------------------------------------
RUN set -e; mv /etc/alternatives /etc/alternatives.smoketest; \
    samtools --version | head -n 1; \
    pigz --version; \
    pandoc --version | head -n 1; \
    python3 -c "import yaml; print('PyYAML', yaml.__version__)"; \
    java -version; \
    out=$(kmc 2>&1 || true);       case "$out" in *"ver. 3.2.4"*) echo "kmc OK";; *) exit 1;; esac; \
    out=$(kmc_tools 2>&1 || true); case "$out" in *"ver. 3.2.4"*) echo "kmc_tools OK";; *) exit 1;; esac; \
    out=$(mergeTags 2>&1 || true); case "$out" in *Usage*) echo "mergeTags OK";; *) exit 1;; esac; \
    v=$(mirtrace --version); echo "miRTrace ${v}"; \
    test -s /opt/surfr2/share/surfr2/artifacts/illumina_adapters.fa; \
    test -s /opt/surfr2/share/surfr2/artifacts/phix174.fa.gz; echo "artifact references OK"; \
    Rscript -e "for (p in c('jsonlite','readr','dplyr','ggplot2','edgeR','plotly','htmlwidgets','shiny')) suppressPackageStartupMessages(library(p, character.only = TRUE)); x <- edgeR::DGEList(matrix(1:8, 2)); cat('R packages OK (BLAS:', La_library(), ')\n')"; \
    mv /etc/alternatives.smoketest /etc/alternatives

CMD ["bash", "-c", "echo 'SURFR2 tool container. Run the pipeline with run_SURFR2.sh; see README_SURFR2.md.'"]
