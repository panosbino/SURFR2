#!/usr/bin/env Rscript
# =============================================================================
# SURFR2 - report: scatter and volcano plots per compared cohort
#
#   results/plots/{scatter,volcano}_<cohort>.pdf          static, for publication
#   results/interactive/{scatter,volcano}_<cohort>.html   self-contained; hovering a
#                                                         passing k-mer shows counts per
#                                                         condition (and per sample for
#                                                         small designs)
# Volcano plots need p-values, so they are made only for cohorts with a statistical test.
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
  types <- c("scatter", if (surfr2_has_test(res, co)) "volcano")
  if (!surfr2_has_test(res, co)) msg("%s: no statistical test (descriptive mode) - no volcano plot", co)
  for (type in types) {
    g <- surfr2_ggplot(res, co, type)
    ggsave(file.path(plots_dir, sprintf("%s_%s.pdf", type, co)), g, width = 7, height = 7.5)
    if (have_plotly) {
      p <- surfr2_plotly(res, co, type)
      out <- file.path(normalizePath(html_dir), sprintf("%s_%s.html", type, co))
      htmlwidgets::saveWidget(p, out, selfcontained = have_pandoc,
                              title = sprintf("SURFR2 %s - %s %s", cfg$project, type, co))
      # saveWidget leaves an unused lib folder when self-contained
      if (have_pandoc) unlink(sub("\\.html$", "_files", out), recursive = TRUE)
    }
    msg("%s: %s plot written%s", co, type, if (have_plotly) " (PDF + interactive HTML)" else " (PDF)")
  }
}
msg("done")
