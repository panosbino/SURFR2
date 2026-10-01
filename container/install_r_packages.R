# =============================================================================
# SURFR2 - install the R packages SURFR2 needs
#
# Used by the Dockerfile and by a development setup on a cluster:
#   Rscript install_r_packages.R [library_dir]
# library_dir defaults to the first entry of .libPaths(); for a personal library
# on a cluster pass e.g. ~/R/surfr2-4.4 and set execution.r_libs to the same path.
#
# install.packages() only WARNS when a package fails to build, so every package is
# loaded afterwards and the script exits non-zero if any is missing.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else .libPaths()[1]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

if (is.null(getOption("repos")) || identical(unname(getOption("repos")["CRAN"]), "@CRAN@"))
  options(repos = c(CRAN = "https://cloud.r-project.org"))

cran <- c(
  "jsonlite", "readr", "dplyr", "ggplot2",   # filter step (required)
  "plotly", "htmlwidgets",                    # interactive report
  "shiny",                                    # explorer app
  "ggvenn",                                   # cohort Venn diagram (optional)
  "BiocManager"
)
bioc <- c("edgeR")                            # statistical test

ncpus <- max(1L, parallel::detectCores() - 1L)
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran))
  install.packages(missing_cran, lib = lib, Ncpus = ncpus)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc))
  BiocManager::install(missing_bioc, lib = lib, update = FALSE, ask = FALSE, Ncpus = ncpus)

ok <- vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)
for (p in names(ok))
  cat(sprintf("%-12s %s\n", p, if (ok[[p]]) as.character(packageVersion(p)) else "MISSING"))
bioc_version <- if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else "unknown"
cat(sprintf("R %s, Bioconductor %s, library %s\n", getRversion(), bioc_version, lib))
if (!all(ok)) {
  message("ERROR: failed to install: ", paste(names(ok)[!ok], collapse = ", "))
  quit(status = 1)
}
