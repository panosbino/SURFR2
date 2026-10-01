#!/usr/bin/env Rscript
# =============================================================================
# SURFR2 - report: scatterplots of case vs control per replicate cohort
#
#   results/plots/scatter_<cohort>.pdf               static, for publication
#   results/interactive/scatter_<cohort>.html        self-contained, hover on
#                                                    case-specific k-mers shows
#                                                    per-condition (and, for small
#                                                    designs, per-sample) counts
#
# Runs after surfr2_filter.R. Separate from it so plots can be regenerated
# without re-filtering (run_SURFR2.sh --from report).
#
# Usage: Rscript surfr2_report.R <run_config_dir>
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)))
source(file.path(dirname(self), "surfr2_plotlib.R"))

run_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(run_dir)) stop("usage: Rscript surfr2_report.R <run_config_dir>")
msg <- function(...) message(sprintf("[%s] [report] %s", format(Sys.time(), "%F %T"), sprintf(...)))

cfg <- fromJSON(file.path(run_dir, "resolved_config.json"), simplifyVector = TRUE)
res_dir <- file.path(cfg$outdir, "results")
if (!file.exists(file.path(res_dir, "run_summary.txt"))) stop("filter step has not completed: ", res_dir)

# Make results/ self-contained for the Shiny app and for sharing
invisible(file.copy(file.path(run_dir, "resolved_config.json"), file.path(res_dir, "resolved_config.json"),
          overwrite = TRUE))

res <- surfr2_load_results(res_dir)
plots_dir <- file.path(res_dir, "plots")
html_dir  <- file.path(res_dir, "interactive")
unlink(html_dir, recursive = TRUE)
dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

have_plotly <- requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)
have_pandoc <- nzchar(Sys.which("pandoc"))
if (!have_plotly) message("WARNING: R packages plotly/htmlwidgets missing - interactive plots skipped")
if (have_plotly && !have_pandoc)
  message("WARNING: pandoc not found - HTML files need their *_files/ folder alongside them")

for (co in cfg$replicate_cohorts) {
  g <- surfr2_scatter_ggplot(res, co)
  ggsave(file.path(plots_dir, sprintf("scatter_%s.pdf", co)), g, width = 7, height = 7)
  msg("%s: static scatter written", co)

  if (have_plotly) {
    p <- surfr2_scatter_plotly(res, co)
    out <- file.path(normalizePath(html_dir), sprintf("scatter_%s.html", co))
    htmlwidgets::saveWidget(p, out, selfcontained = have_pandoc,
                            title = sprintf("SURFR2 %s - %s", cfg$project, co))
    # saveWidget leaves an empty/unused lib folder when self-contained
    if (have_pandoc) unlink(sub("\\.html$", "_files", out), recursive = TRUE)
    msg("%s: interactive scatter -> %s", co, out)
  }
}
msg("done")
