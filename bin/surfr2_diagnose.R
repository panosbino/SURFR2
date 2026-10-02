#!/usr/bin/env Rscript
# =============================================================================
# SURFR2 - diagnose a differential-expression result with few or no significant calls
#
# Usage: Rscript surfr2_diagnose.R <run_config_dir> [cohort]
#   e.g.  Rscript bin/surfr2_diagnose.R <outdir>/runs/latest
# Run in the same environment as the pipeline (modules/r_libs or the container).
#
# Prints, and plots to <outdir>/results/plots/diagnostics_<cohort>.pdf:
#  1. p-value distribution  - is there ANY signal (excess of small p-values)?
#  2. sample similarity     - do samples group by condition, or by something else
#                             (pairs, species, batch)? Which samples pair up?
#  3. variability (BCV)     - how much do replicates disagree?
#  4. replicate agreement   - are k-mers found in all replicates of a condition, or
#                             mostly in a single sample?
#  5. positive controls     - statistics of k-mers of known sequences (optional,
#                             SURFR2_CONTROLS="name=SEQ,name=SEQ"; default: miR-122)
# =============================================================================

suppressPackageStartupMessages({ library(jsonlite); library(readr); library(dplyr); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript surfr2_diagnose.R <run_config_dir> [cohort]")
run_dir <- args[1]
cfg <- fromJSON(file.path(run_dir, "resolved_config.json"), simplifyVector = TRUE)
co  <- if (length(args) >= 2) args[2] else cfg$replicate_cohorts[1]
CASE <- cfg$comparison$case; CTRL <- cfg$comparison$control
res_dir <- file.path(cfg$outdir, "results"); mdir <- file.path(cfg$outdir, "matrix")
out <- function(...) cat(sprintf(...), "\n", sep = "")
hr <- function(t) cat("\n==== ", t, " ", strrep("=", max(0, 60 - nchar(t))), "\n", sep = "")

samples <- read_tsv(file.path(run_dir, "samples.tsv"), show_col_types = FALSE, col_types = cols(.default = "c"))
norm <- read_tsv(file.path(res_dir, "normalisation.tsv"), show_col_types = FALSE)
ss <- samples |> filter(cohort == co, role %in% c("case", "control")) |>
  left_join(select(norm, sample_id, qc_reads, eff_libsize), by = "sample_id") |>
  arrange(factor(role, levels = c("case", "control")), sample_id)
stats <- read_tsv(file.path(res_dir, "kmer_stats.tsv.gz"), show_col_types = FALSE, progress = FALSE)
pdf_file <- file.path(res_dir, "plots", sprintf("diagnostics_%s.pdf", co))
pdf(pdf_file, width = 8, height = 6)

# ---- 1. p-values ---------------------------------------------------------------------
hr("1. p-value distribution")
pcol <- paste0(co, "_pvalue")
if (pcol %in% names(stats)) {
  key <- paste(stats[[pcol]], stats[[paste0(co, "_edger_log2fc_", CASE, "_vs_", CTRL)]])
  p <- stats[[pcol]][!duplicated(key) & !is.na(stats[[pcol]])]   # one per distinct profile
  out("tested distinct profiles: %s", format(length(p), big.mark = ","))
  out("smallest p-value: %.3g   smallest FDR: %.3g", min(p), min(stats[[paste0(co, "_fdr")]], na.rm = TRUE))
  out("fraction p < 0.001: %.4f   < 0.01: %.4f   < 0.05: %.4f   (uniform/no signal: 0.001, 0.01, 0.05)",
      mean(p < 0.001), mean(p < 0.01), mean(p < 0.05))
  out("smallest p needed for FDR 0.05 with this many tests: about %.2g", 0.05 / length(p))
  hist(p, breaks = 50, main = sprintf("%s: p-values of %s tested profiles", co, format(length(p), big.mark = ",")),
       xlab = "p-value", col = "grey80", border = "white")
} else out("no p-values (descriptive mode)")

# ---- 2. sample similarity ----------------------------------------------------------------
hr("2. sample similarity (log CPM of abundant k-mers)")
cand <- read_tsv(file.path(mdir, "candidates.tsv"), col_names = c("kmer", "n"), show_col_types = FALSE)
M <- sapply(ss$sample_id, function(id) {
  d <- read_tsv(file.path(mdir, "counts", "candidates", paste0(id, ".tsv.gz")),
                col_names = c("kmer", "count"), show_col_types = FALSE, progress = FALSE)
  v <- numeric(nrow(cand)); v[match(d$kmer, cand$kmer)] <- d$count; v })
rownames(M) <- cand$kmer
cpm <- sweep(M, 2, ss$eff_libsize, "/") * 1e6
prof <- !duplicated(apply(M, 1, paste, collapse = ","))
top <- head(order(rowMeans(cpm) * prof, decreasing = TRUE), 50000)
lc <- log2(cpm[top, , drop = FALSE] + 1)
C <- cor(lc, method = "spearman"); dimnames(C) <- list(ss$sample_id, ss$sample_id)
out("conditions: %s", paste(sprintf("%s=%s", ss$sample_id, ss$condition), collapse = ", "))
options(width = 250); print(round(C, 2))
within <- mean(c(C[ss$role == "case", ss$role == "case"][upper.tri(diag(sum(ss$role == "case")))],
                 C[ss$role == "control", ss$role == "control"][upper.tri(diag(sum(ss$role == "control")))]))
between <- mean(C[ss$role == "case", ss$role == "control"])
out("mean correlation within %s/%s: %.3f   between them: %.3f", CASE, CTRL, within, between)
if (within <= between)
  out("!! samples are NOT more similar within a condition than between conditions: replicates disagree")
for (i in which(ss$role == "case")) {
  j <- which(ss$role == "control")[which.max(C[i, ss$role == "control"])]
  out("  %s is most similar to %s %s (r = %.3f)", ss$sample_id[i], CTRL, ss$sample_id[j], C[i, j])
}
hc <- hclust(as.dist(1 - C), method = "average")
plot(hc, main = sprintf("%s: sample clustering (1 - Spearman r)", co), xlab = "", sub = "")
mds <- cmdscale(as.dist(1 - C), k = 2)
plot(mds, pch = 19, col = ifelse(ss$role == "case", "#e6550d", "#3182bd"),
     main = sprintf("%s: MDS of samples", co), xlab = "dim 1", ylab = "dim 2")
text(mds, labels = ss$sample_id, pos = 3, cex = 0.7)
legend("topright", legend = c(CASE, CTRL), col = c("#e6550d", "#3182bd"), pch = 19)

# ---- 3. variability --------------------------------------------------------------------
hr("3. biological coefficient of variation (edgeR)")
if (requireNamespace("edgeR", quietly = TRUE) && min(table(ss$role)) >= 2) {
  cond <- factor(ss$role, levels = c("control", "case")); X <- model.matrix(~ cond)
  y <- edgeR::DGEList(M[top, , drop = FALSE], lib.size = ss$eff_libsize)
  y <- edgeR::estimateDisp(y, X)
  out("common BCV: %.2f  (well-controlled replicates: ~0.1-0.4; > 1: replicates differ wildly)",
      sqrt(y$common.dispersion))
  edgeR::plotBCV(y, main = sprintf("%s: BCV (%d abundant profiles)", co, length(top)))
}

# ---- 4. replicate agreement ---------------------------------------------------------------
hr("4. replicate agreement: in how many samples is each candidate k-mer found?")
for (cn in c(CASE, CTRL)) {
  nf <- stats[[paste0("n_found_", cn)]]; n <- sum(ss$condition == cn)
  tab <- table(factor(nf, levels = 0:max(n, max(nf))))
  out("%-10s %s", cn, paste(sprintf("%s:%.1f%%", names(tab), 100 * tab / sum(tab)), collapse = "  "))
}
out("(a high share found in a single sample means replicates contain different sequences,")
out(" e.g. different species, strains or contamination - typical for mixed-species designs)")

# ---- 5. positive controls ----------------------------------------------------------------
hr("5. positive controls")
ctrl_spec <- Sys.getenv("SURFR2_CONTROLS", "miR-122=TGGAGTGTGACAATGGTGTTTG")
for (item in strsplit(ctrl_spec, ",")[[1]]) {
  nm <- sub("=.*", "", item); sq <- toupper(chartr("U", "T", sub(".*=", "", item)))
  k <- cfg$kmer$k; n <- nchar(sq) - k + 1
  km <- substring(sq, seq_len(n), seq_len(n) + k - 1)
  r <- stats[stats$kmer %in% km, , drop = FALSE]
  out("%s (%s): %d of %d k-mers among candidates", nm, sq, nrow(r), n)
  if (nrow(r)) {
    cols <- c("kmer", paste0("n_found_", c(CASE, CTRL)), paste0(co, "_mean_cpm_", c(CASE, CTRL)),
              paste0(co, "_log2fc_", CASE, "_vs_", CTRL), pcol, paste0(co, "_fdr"))
    print(as.data.frame(r[, intersect(cols, names(r))]), digits = 3, row.names = FALSE)
  }
}
invisible(dev.off())
out("\nplots: %s", pdf_file)
