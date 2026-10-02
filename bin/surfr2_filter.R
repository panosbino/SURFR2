#!/usr/bin/env Rscript
# =============================================================================
# SURFR2 - normalisation, statistics and calling of k-mers
#
# Two analyses (config: analysis), both comparing case vs control per cohort:
#   dea       differential expression: k-mers significantly higher in EITHER condition
#   specific  k-mers present in the case condition and absent from the control
#
#  1. Normalisation -> effective library size L_i; CPM = count / L_i * 1e6
#       cpm          : L_i = QC-passed reads of sample i (counts per million reads)
#       median_ratio : DESeq2 median-of-ratios size factors on k-mers present in all
#                      samples, rescaled to L_i = s_i * geomean(QC-passed reads)
#  2. One pass over the samples: per-group sums/maxima of CPM, detection counts, and a
#     bit-packed record of which samples contain each k-mer (for found_in_* columns).
#  3. Per compared cohort:
#       edgeR quasi-likelihood test when the design allows ('test' mode)
#       dea:      FDR <= statistics.fdr (test mode), |log2FC| >= dea.min_abs_log2fc,
#                 mean CPM in the higher condition >= dea.min_mean_cpm
#       specific: detection, abundance, absence and fold-change filters, plus FDR
#     A k-mer must pass in every compared cohort (dea: in the same direction) and, for
#     specific, in every external control.
#  4. Matrices of passing k-mers, dekupl-mergeTags, QC plots, run summary.
#
# Results are labelled with the CONDITION NAMES from the config, never case/control.
# Usage: Rscript surfr2_filter.R <run_config_dir>
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

msg <- function(...) message(sprintf("[%s] [filter] %s", format(Sys.time(), "%F %T"), sprintf(...)))
warnings_log <- character()
warn <- function(...) {
  w <- sprintf(...)
  warnings_log <<- c(warnings_log, w)
  message("WARNING: ", w)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript surfr2_filter.R <run_config_dir>")
run_dir <- args[1]

cfg      <- fromJSON(file.path(run_dir, "resolved_config.json"), simplifyVector = TRUE)
analysis <- cfg$analysis
f        <- cfg$filters
sp       <- cfg$specific
dea_cfg  <- cfg$dea
st       <- cfg$statistics
design   <- cfg$design
CASE     <- cfg$comparison$case
CTRL     <- cfg$comparison$control
VS       <- paste0(CASE, "_vs_", CTRL)        # log2FC = log2(CASE / CTRL)
pc       <- f$pseudocount_cpm
prefix   <- if (analysis == "dea") "de" else "specific"
label    <- if (analysis == "dea") "differentially expressed" else paste0(CASE, "-specific")

mdir <- file.path(cfg$outdir, "matrix")
if (!file.exists(file.path(mdir, ".done"))) stop("matrix step incomplete: ", mdir)

res_dir   <- file.path(cfg$outdir, "results")
plots_dir <- file.path(res_dir, "plots")
# Start from an empty results directory: a re-run that produces fewer outputs (e.g. no
# k-mer passes) must not leave files from an earlier run looking current.
unlink(res_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Provenance: the tools this run used (written by the launcher's tool check) and R details
tv <- file.path(run_dir, "tool_versions.txt")
if (file.exists(tv)) invisible(file.copy(tv, file.path(res_dir, "tool_versions.txt")))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo.txt"))

samples <- read_tsv(file.path(run_dir, "samples.tsv"), show_col_types = FALSE,
                    col_types = cols(.default = col_character(), index = col_integer()))
lib <- read_tsv(file.path(mdir, "library_sizes.tsv"), show_col_types = FALSE,
                col_types = cols(.default = col_character(), qc_reads = col_double(),
                                 total_kmers = col_double(), unique_kmers = col_double()))
stopifnot(setequal(samples$sample_id, lib$sample_id), !anyDuplicated(lib$sample_id))
extra_cols <- intersect(c("input_reads", "umi_molecules", "artifact_reads"), names(lib))
lib <- lib |> mutate(across(all_of(extra_cols), ~ suppressWarnings(as.numeric(.x))))
samples <- samples |> left_join(select(lib, sample_id, qc_reads, total_kmers, unique_kmers, all_of(extra_cols)),
                                by = "sample_id")
samples$block[is.na(samples$block)] <- ""

count_file <- function(sub, id) file.path(mdir, "counts", sub, paste0(id, ".tsv.gz"))
read_counts <- function(sub, id) {
  read_tsv(count_file(sub, id), col_names = c("kmer", "count"), show_col_types = FALSE,
           col_types = cols(kmer = col_character(), count = col_double()), progress = FALSE)
}
geomean <- function(x) exp(mean(log(x)))
pmin_cols <- function(df, cols) do.call(pmin, unname(as.list(df[cols])))
pmax_cols <- function(df, cols) do.call(pmax, unname(as.list(df[cols])))

# Cauchy combination (ACAT) of possibly dependent p-values
acat <- function(p) {
  p <- p[!is.na(p)]
  if (!length(p)) return(NA_real_)
  p <- pmin(p, 1 - 1e-15)
  t <- ifelse(p < 1e-15, 1 / (pmax(p, 1e-300) * pi), tan((0.5 - p) * pi))
  tt <- mean(t)
  if (tt > 1e15) 1 / (tt * pi) else 0.5 - atan(tt) / pi
}

# =============================================================================
# 1. NORMALISATION (library size = QC-passed reads)
# =============================================================================
# A k-mer occurs at most about once per read, so count / reads is the fraction of reads
# containing it: CPM = counts per million QC-passed reads. QC-passed reads (not raw input
# reads) are the population the k-mers were counted from; raw reads would make the scale
# depend on each library's QC failure rate (e.g. adapter dimers).
msg("normalisation: %s (library size: QC-passed reads)", cfg$normalization$method)

if (cfg$normalization$method == "cpm") {
  samples$size_factor <- samples$qc_reads / geomean(samples$qc_reads)
  samples$eff_libsize <- samples$qc_reads
} else {
  first <- read_counts("reference", samples$sample_id[1])
  n_ref_all <- nrow(first)
  step <- max(1L, ceiling(n_ref_all / cfg$normalization$max_reference_kmers))
  ref_kmers <- first$kmer[seq(1, n_ref_all, by = step)]   # deterministic subsample
  M <- matrix(NA_real_, nrow = length(ref_kmers), ncol = nrow(samples),
              dimnames = list(NULL, samples$sample_id))
  for (j in seq_len(nrow(samples))) {
    d <- read_counts("reference", samples$sample_id[j])
    M[, j] <- d$count[match(ref_kmers, d$kmer)]
  }
  if (anyNA(M) || any(M <= 0)) stop("reference k-mer missing/zero in some sample - matrix step inconsistent")
  logM <- log(M)
  sf <- exp(apply(logM - rowMeans(logM), 2, median))
  sf <- sf / geomean(sf)
  samples$size_factor <- unname(sf[samples$sample_id])
  samples$eff_libsize <- samples$size_factor * geomean(samples$qc_reads)
  msg("median-of-ratios on %d of %d reference k-mers", length(ref_kmers), n_ref_all)
  r <- cor(log(samples$size_factor), log(samples$qc_reads))
  msg("cor(log size factor, log QC-passed reads) = %.3f", r)
  if (r < 0.5) warn(paste0("size factors correlate weakly with depth (r=%.2f): strong composition ",
                           "differences between samples; inspect plots/normalisation.pdf"), r)
  rm(M, logM)
}

write_tsv(samples |> select(sample_id, cohort, condition, any_of(c("input_reads", "umi_molecules", "artifact_reads")),
                            qc_reads, total_kmers, unique_kmers, size_factor, eff_libsize),
          file.path(res_dir, "normalisation.tsv"))

# =============================================================================
# 2. GROUPS
# =============================================================================
rep_cohorts <- cfg$replicate_cohorts
samples$group <- NA_character_
in_rep <- samples$cohort %in% rep_cohorts & samples$role %in% c("case", "control")
samples$group[in_rep] <- paste(samples$cohort[in_rep], samples$role[in_rep], sep = "|")
is_ext <- samples$role == "external"
samples$group[is_ext] <- paste("external", samples$condition[is_ext], sep = "|")

# Detection threshold per sample (specific analysis): case uses min_cpm_case, others max_cpm_control
samples$expr_thr <- if (analysis == "specific") {
  ifelse(samples$role == "case", sp$min_cpm_case, sp$max_cpm_control)
} else 0

groups <- sort(unique(na.omit(samples$group)))
n_in_group <- table(samples$group)[groups]
msg("groups: %s", paste(sprintf("%s(n=%d)", groups, as.integer(n_in_group)), collapse = ", "))

# Specific analysis: the raw-count pre-filter (>= min_count reads) is only guaranteed to be
# looser than the CPM detection filter if min_count <= min_cpm_case * L_i / 1e6 in every
# case sample.
if (analysis == "specific") {
  case_rep <- samples$role == "case" & samples$cohort %in% rep_cohorts
  count_equiv <- sp$min_cpm_case * samples$eff_libsize[case_rep] / 1e6
  n_viol <- sum(cfg$candidates$min_count > count_equiv + 1e-9)
  if (n_viol > 0) {
    warn(paste0("candidates.min_count=%d exceeds min_cpm_case*libsize/1e6 in %d/%d %s samples ",
                "(smallest equivalent count %.1f): detection in %s may be underestimated near ",
                "threshold. Lower candidates.min_count to <= %d."),
         cfg$candidates$min_count, n_viol, sum(case_rep), CASE, min(count_equiv), CASE,
         max(1L, floor(min(count_equiv))))
  }
}

# =============================================================================
# 3. ONE PASS OVER ALL SAMPLES
# =============================================================================
cand <- read_tsv(file.path(mdir, "candidates.tsv"), col_names = c("kmer", "n_samples_prefilter"),
                 show_col_types = FALSE,
                 col_types = cols(kmer = col_character(), n_samples_prefilter = col_integer()))
n_cand <- nrow(cand)
msg("%d candidate k-mers; reading %d samples", n_cand, nrow(samples))

G <- length(groups)
sum_cpm <- matrix(0,  nrow = n_cand, ncol = G, dimnames = list(NULL, groups))
max_cpm <- matrix(0,  nrow = n_cand, ncol = G, dimnames = list(NULL, groups))
n_expr  <- matrix(0L, nrow = n_cand, ncol = G, dimnames = list(NULL, groups))

# Which samples contain each k-mer (count >= 1): 30 samples per integer column, so the
# memory cost is n_candidates x ceiling(n_samples / 30) x 4 bytes.
BITS <- 30L
found_bits <- matrix(0L, nrow = n_cand, ncol = ceiling(nrow(samples) / BITS))

for (j in seq_len(nrow(samples))) {
  s <- samples[j, ]
  d <- read_counts("candidates", s$sample_id)
  idx <- match(d$kmer, cand$kmer)
  if (anyNA(idx)) stop("sample ", s$sample_id, " has k-mers outside the candidate set")
  col <- (j - 1L) %/% BITS + 1L
  found_bits[idx, col] <- bitwOr(found_bits[idx, col], bitwShiftL(1L, (j - 1L) %% BITS))
  if (!is.na(s$group)) {
    cpm <- d$count / s$eff_libsize * 1e6
    g <- s$group
    sum_cpm[idx, g] <- sum_cpm[idx, g] + cpm
    max_cpm[idx, g] <- pmax(max_cpm[idx, g], cpm)
    n_expr[idx, g]  <- n_expr[idx, g] + as.integer(cpm >= s$expr_thr - 1e-12 & d$count > 0)
  }
  if (j %% 50 == 0) msg("  %d / %d samples", j, nrow(samples))
}
is_found <- function(j) bitwAnd(found_bits[, (j - 1L) %/% BITS + 1L], bitwShiftL(1L, (j - 1L) %% BITS)) != 0L

stats <- cand
# found_in_<condition>: comma-separated samples (all cohorts) with >= 1 read of the k-mer
conds <- unique(c(CASE, CTRL, samples$condition[samples$role == "external"]))
for (cn in conds) {
  js <- which(samples$condition == cn)
  lst <- rep("", n_cand); n <- integer(n_cand)
  for (j in js) {
    h <- is_found(j)
    lst[h] <- ifelse(n[h] == 0L, samples$sample_id[j], paste0(lst[h], ",", samples$sample_id[j]))
    n <- n + h
  }
  stats[[paste0("n_found_", cn)]] <- n
  stats[[paste0("found_in_", cn)]] <- lst
}

# =============================================================================
# 4. PER-COHORT STATISTICS, TEST AND CALLS
# =============================================================================
pass_all <- rep(TRUE, n_cand)
direction <- rep(NA_character_, n_cand)     # dea: condition with higher expression

# edgeR quasi-likelihood test of case vs control for all candidates in one cohort.
# Normalisation is SURFR2's own (effective library sizes, norm.factors = 1): TMM on the
# candidate set would be biased, because candidates are enriched for differential k-mers.
run_edger <- function(ss) {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("R package 'edgeR' is required for statistics.test")
  M <- matrix(0, nrow = n_cand, ncol = nrow(ss))
  msg("  edgeR: %d candidates x %d samples (~%.1f GB for the count matrix)",
      n_cand, nrow(ss), n_cand * nrow(ss) * 8 / 1e9)
  for (j in seq_len(nrow(ss))) {
    d <- read_counts("candidates", ss$sample_id[j])
    M[match(d$kmer, cand$kmer), j] <- d$count
  }
  keep <- which(rowSums(M) > 0)    # absent from this cohort: untestable (condition-blind)
  # Overlapping k-mers of one molecule have IDENTICAL count profiles (a 22-nt read yields
  # six 17-mers). Fed to edgeR as separate features they are pseudo-replicates: they inflate
  # the empirical-Bayes prior df (over-trusting the dispersion trend) and multiply the
  # number of BH tests. Identical profiles necessarily get identical results, so each
  # distinct profile is tested once and the result mapped back - lossless for p-values.
  Mk  <- M[keep, , drop = FALSE]
  key <- do.call(paste, c(asplit(Mk, 2), sep = ","))
  uniq <- !duplicated(key)
  map  <- match(key, key[uniq])
  msg("  edgeR: %d k-mers collapse into %d distinct count profiles", length(keep), sum(uniq))
  cond <- factor(ss$role, levels = c("control", "case"))
  X <- if (isTRUE(design$cohorts[[ss$cohort[1]]]$blocked)) {
    model.matrix(~ factor(ss$block) + cond)
  } else model.matrix(~ cond)
  y <- edgeR::DGEList(counts = Mk[uniq, , drop = FALSE], lib.size = ss$eff_libsize)
  y <- edgeR::estimateDisp(y, X)
  fit <- edgeR::glmQLFit(y, X, robust = TRUE)
  res <- edgeR::glmQLFTest(fit, coef = ncol(X))$table
  res$FDR <- p.adjust(res$PValue, method = "BH")      # over distinct profiles
  out <- data.frame(pvalue = rep(NA_real_, n_cand), fdr = NA_real_, logFC = NA_real_)
  out$pvalue[keep] <- res$PValue[map]
  out$logFC[keep]  <- res$logFC[map]
  out$fdr[keep]    <- res$FDR[map]
  msg("  edgeR: %d of %d profiles at FDR <= %g (either direction; QL prior df %.1f)",
      sum(res$FDR <= st$fdr), nrow(res), st$fdr, median(fit$df.prior))
  out
}

mode_of <- character()
for (co in rep_cohorts) {
  dz <- design$cohorts[[co]]
  gc <- paste(co, "case", sep = "|"); gn <- paste(co, "control", sep = "|")
  nc <- dz$n_case; nn <- dz$n_control
  stopifnot(nc == as.integer(n_in_group[gc]), nn == as.integer(n_in_group[gn]))
  mean_c <- sum_cpm[, gc] / nc
  mean_n <- sum_cpm[, gn] / nn
  lfc <- log2((mean_c + pc) / (mean_n + pc))
  col <- function(x) paste0(co, "_", x)

  stats[[col(paste0("mean_cpm_", CASE))]] <- mean_c
  stats[[col(paste0("mean_cpm_", CTRL))]] <- mean_n
  stats[[col(paste0("log2fc_", VS))]]     <- lfc
  if (analysis == "specific") {
    stats[[col(paste0("n_detected_", CASE))]] <- n_expr[, gc]
    stats[[col(paste0("n_detected_", CTRL))]] <- n_expr[, gn]
    stats[[col(paste0("max_cpm_", CTRL))]]    <- max_cpm[, gn]
  }

  et <- NULL
  if (dz$mode == "test") {
    ss <- samples[samples$cohort == co & samples$role %in% c("case", "control"), ]
    et <- run_edger(ss)
    stats[[col(paste0("edger_log2fc_", VS))]] <- et$logFC
    stats[[col("pvalue")]] <- et$pvalue
    stats[[col("fdr")]]    <- et$fdr
  } else {
    warn("cohort '%s' ran in DESCRIPTIVE mode (%s n=%d, %s n=%d): no statistical test; results are exploratory",
         co, CASE, nc, CTRL, nn)
  }

  if (analysis == "dea") {
    # fold change: edgeR's estimate when tested, otherwise the descriptive CPM ratio
    lfc_used <- if (is.null(et)) lfc else et$logFC
    higher <- ifelse(is.na(lfc_used), NA_character_, ifelse(lfc_used > 0, CASE, CTRL))
    mean_higher <- ifelse(lfc_used > 0, mean_c, mean_n)
    pass <- !is.na(lfc_used) & abs(lfc_used) >= dea_cfg$min_abs_log2fc &
            mean_higher >= dea_cfg$min_mean_cpm
    if (!is.null(et)) pass <- pass & !is.na(et$fdr) & et$fdr <= st$fdr
    stats[[col("higher_in")]] <- higher
    # replication: the same direction in every compared cohort
    agree <- is.na(direction) | direction == higher
    direction <- ifelse(is.na(direction), higher, direction)
    pass_all <- pass_all & pass & !is.na(higher) & agree
    msg("%s [%s]: %d k-mers pass (higher in %s: %d, higher in %s: %d)", co, dz$mode, sum(pass),
        CASE, sum(pass & higher == CASE, na.rm = TRUE), CTRL, sum(pass & higher == CTRL, na.rm = TRUE))
  } else {
    pass <- n_expr[, gc] >= dz$need_case &
            mean_c >= sp$min_mean_cpm_case &
            n_expr[, gn] <= dz$allow_control &
            lfc >= sp$min_log2fc
    if (!is.null(et)) pass <- pass & !is.na(et$fdr) & et$fdr <= st$fdr
    msg("%s [%s]: %s n=%d (need >= %d detected), %s n=%d (allow <= %d detected): %d pass",
        co, dz$mode, CASE, nc, dz$need_case, CTRL, nn, dz$allow_control, sum(pass))
    pass_all <- pass_all & pass
  }
  mode_of[co] <- dz$mode
  stats[[col("pass")]] <- pass
}

ext_groups <- grep("^external\\|", groups, value = TRUE)
for (ge in ext_groups) {            # specific analysis only (validator)
  cond <- sub("^external\\|", "", ge)
  ne <- as.integer(n_in_group[ge])
  allow_ext <- design$external[[cond]]$allow
  pass <- n_expr[, ge] <= allow_ext
  msg("external %s: n=%d (allow <= %d detected): %d candidates pass", cond, ne, allow_ext, sum(pass))
  stats[[paste0("external_", cond, "_mean_cpm")]]   <- sum_cpm[, ge] / ne
  stats[[paste0("external_", cond, "_max_cpm")]]    <- max_cpm[, ge]
  stats[[paste0("external_", cond, "_n_detected")]] <- n_expr[, ge]
  stats[[paste0("external_", cond, "_pass")]]       <- pass
  pass_all <- pass_all & pass
}
rm(sum_cpm, max_cpm, n_expr, found_bits); invisible(gc())

# ---- summary columns across cohorts ------------------------------------------------
fdr_cols <- intersect(paste0(rep_cohorts, "_fdr"), names(stats))
p_cols   <- intersect(paste0(rep_cohorts, "_pvalue"), names(stats))
if (length(fdr_cols)) {
  stats$max_fdr    <- pmax_cols(stats, fdr_cols)
  stats$max_pvalue <- pmax_cols(stats, p_cols)
}
if (analysis == "dea") {
  stats$higher_in <- direction
  # signed fold change with the smallest magnitude across cohorts (conservative)
  fc_cols <- vapply(rep_cohorts, function(co) {
    e <- paste0(co, "_edger_log2fc_", VS)
    if (e %in% names(stats)) e else paste0(co, "_log2fc_", VS)
  }, character(1))
  fc <- as.matrix(stats[fc_cols])
  stats[[paste0("log2fc_", VS)]] <- fc[cbind(seq_len(n_cand), max.col(-abs(replace(fc, is.na(fc), Inf)), ties.method = "first"))]
} else {
  stats[[paste0("min_log2fc_", VS)]]       <- pmin_cols(stats, paste0(rep_cohorts, "_log2fc_", VS))
  stats[[paste0("min_mean_cpm_", CASE)]]   <- pmin_cols(stats, paste0(rep_cohorts, "_mean_cpm_", CASE))
  stats[[paste0("max_mean_cpm_", CTRL)]]   <- pmax_cols(stats, paste0(rep_cohorts, "_mean_cpm_", CTRL))
}
stats$pass_all <- pass_all

write_tsv(stats, file.path(res_dir, "kmer_stats.tsv.gz"))
final <- stats |> filter(pass_all)
fc_name <- if (analysis == "dea") paste0("log2fc_", VS) else paste0("min_log2fc_", VS)
ord <- order(if (length(fdr_cols)) final$max_fdr else rep(0, nrow(final)), -abs(final[[fc_name]]))
final <- final[ord, , drop = FALSE]
write_tsv(final, file.path(res_dir, paste0(prefix, "_kmers.tsv")))
if (analysis == "dea") {
  msg("%d differentially expressed k-mers (higher in %s: %d, higher in %s: %d) in all of: %s",
      nrow(final), CASE, sum(final$higher_in == CASE), CTRL, sum(final$higher_in == CTRL),
      paste(rep_cohorts, collapse = ", "))
} else {
  msg("%d %s k-mers pass in all of: %s%s", nrow(final), label, paste(rep_cohorts, collapse = ", "),
      if (length(ext_groups)) paste0(" + external ", paste(sub("^external\\|", "", ext_groups), collapse = ", ")) else "")
}

# =============================================================================
# 5. MATRICES OF PASSING K-MERS (all samples)
# =============================================================================
if (nrow(final) > 0) {
  cnt <- matrix(0, nrow = nrow(final), ncol = nrow(samples),
                dimnames = list(final$kmer, samples$sample_id))
  for (j in seq_len(nrow(samples))) {
    d <- read_counts("candidates", samples$sample_id[j])
    hit <- match(final$kmer, d$kmer)
    cnt[, j] <- ifelse(is.na(hit), 0, d$count[hit])
  }
  cpm_mat <- sweep(cnt, 2, samples$eff_libsize, "/") * 1e6
  write_tsv(as_tibble(cnt, rownames = "kmer"),     file.path(res_dir, paste0(prefix, "_kmers_counts.tsv")))
  write_tsv(as_tibble(cpm_mat, rownames = "kmer"), file.path(res_dir, paste0(prefix, "_kmers_cpm.tsv")))
} else {
  warn("no k-mer passed; matrices, merging and k-mer plots are skipped")
}

# =============================================================================
# 6. dekupl-mergeTags: assemble overlapping / offset k-mers into sequences
# =============================================================================
if (isTRUE(cfg$mergetags$enabled) && nrow(final) > 0) {
  mt_in  <- file.path(res_dir, "mergetags_input.tsv")
  mt_out <- file.path(res_dir, paste0(prefix, "_sequences.tsv"))
  # mergeTags hard-codes DE-kupl's layout (see mergeTags.c): value column 1 is read as a
  # p-value (the LOWEST becomes the contig representative) and value column 4 as log2FC
  # (> 0 'up', otherwise 'down'; the two sets are merged separately). It reads column 4
  # unconditionally, so >= 4 value columns are mandatory. Without p-values (descriptive
  # mode, or the specific analysis' ranking) a rank (1 = best) fills the p-value slot.
  use_p <- analysis == "dea" && length(p_cols) > 0
  rank_or_p <- if (use_p) final$max_pvalue else seq_len(nrow(final))
  mean_case <- pmin_cols(final, paste0(rep_cohorts, "_mean_cpm_", CASE))
  mean_ctrl <- pmax_cols(final, paste0(rep_cohorts, "_mean_cpm_", CTRL))
  mt_tab <- tibble(tag = final$kmer, !!(if (use_p) "pvalue" else "rank") := rank_or_p,
                   !!paste0("mean_cpm_", CASE) := mean_case,
                   !!paste0("mean_cpm_", CTRL) := mean_ctrl,
                   !!paste0("log2fc_", VS) := final[[fc_name]])
  write_tsv(mt_tab, mt_in)
  mt_args <- c("-k", cfg$kmer$k, "-m", cfg$mergetags$min_overlap)
  # '-n' means UNSTRANDED merging (mergeTags.c: case 'n': stranded = 0)
  if (isTRUE(cfg$kmer$canonical)) mt_args <- c(mt_args, "-n")
  status <- system2(cfg$tools$mergetags, c(mt_args, mt_in), stdout = mt_out)
  if (!identical(as.integer(status), 0L)) stop("dekupl-mergeTags failed with status ", status)
  contigs <- read_tsv(mt_out, show_col_types = FALSE)
  msg("mergeTags: %d k-mers -> %d sequences", nrow(final), nrow(contigs))

  if (nrow(contigs)) {
    k <- cfg$kmer$k
    members <- lapply(contigs$contig, function(sq) {
      n <- nchar(sq) - k + 1
      intersect(substring(sq, seq_len(n), seq_len(n) + k - 1), final$kmer)
    })
    if (analysis == "dea")
      contigs$higher_in <- ifelse(contigs[[paste0("log2fc_", VS)]] > 0, CASE, CTRL)
    # samples in which ANY k-mer of the sequence is found
    for (cn in conds) {
      fl <- setNames(final[[paste0("found_in_", cn)]], final$kmer)
      contigs[[paste0("found_in_", cn)]] <- vapply(members, function(m) {
        ids <- unique(unlist(strsplit(fl[m][nzchar(fl[m])], ",", fixed = TRUE)))
        paste(samples$sample_id[samples$sample_id %in% ids], collapse = ",")
      }, character(1))
    }
    # FDR is controlled at the k-mer level; overlapping k-mers of one molecule are strongly
    # correlated, so per-sequence p-values use the Cauchy combination (ACAT; Liu & Xie
    # 2020, JASA), valid under dependence.
    for (co in sub("_fdr$", "", fdr_cols)) {
      pv <- setNames(final[[paste0(co, "_pvalue")]], final$kmer)
      fd <- setNames(final[[paste0(co, "_fdr")]], final$kmer)
      contigs[[paste0(co, "_acat_pvalue")]]  <- vapply(members, function(m) acat(pv[m]), numeric(1))
      contigs[[paste0(co, "_min_kmer_fdr")]] <- vapply(members, function(m) min(fd[m]), numeric(1))
    }
    write_tsv(contigs, mt_out)
  }
}

# =============================================================================
# 7. QC PLOTS
# =============================================================================
p_lib <- ggplot(samples, aes(x = interaction(condition, cohort, sep = "\n"), y = qc_reads)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  scale_y_log10() + theme_bw() +
  labs(x = NULL, y = "QC-passed reads (log10)", title = "Library sizes")
ggsave(file.path(plots_dir, "library_sizes.pdf"), p_lib, width = 8, height = 5)

if (cfg$normalization$method == "median_ratio") {
  p_sf <- ggplot(samples, aes(qc_reads / geomean(qc_reads), size_factor, colour = condition)) +
    geom_abline(linetype = "dashed") + geom_point() + scale_x_log10() + scale_y_log10() + theme_bw() +
    labs(x = "relative library size (QC-passed reads)", y = "median-of-ratios size factor",
         title = "Deviation from the diagonal = composition effect")
  ggsave(file.path(plots_dir, "normalisation.pdf"), p_sf, width = 7, height = 5)
}

if (nrow(final) > 0) {
  # One point per k-mer per (cohort, condition): mean CPM across that group's samples
  grp <- paste(samples$cohort, samples$condition, sep = " ")
  per_group <- lapply(split(seq_len(nrow(samples)), grp), function(j)
    rowMeans(cpm_mat[, j, drop = FALSE]))
  vdf <- data.frame(group = rep(names(per_group), each = nrow(final)),
                    mean_cpm = unlist(per_group, use.names = FALSE))
  vdf$condition <- setNames(samples$condition, grp)[vdf$group]
  if (analysis == "dea")
    vdf$higher_in <- rep(final$higher_in, times = length(per_group))
  p_v <- ggplot(vdf, aes(group, log10(mean_cpm + pc), fill = condition)) +
    geom_violin(scale = "width") +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.5) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = sprintf("log10(mean CPM + %g)", pc),
         title = sprintf("%d %s k-mers", nrow(final), label))
  if (analysis == "dea") p_v <- p_v + facet_wrap(~ paste("higher in", higher_in))
  ggsave(file.path(plots_dir, paste0(prefix, "_kmers_violin.pdf")), p_v, width = 10, height = 5)
}

if (length(rep_cohorts) %in% 2:4 && requireNamespace("ggvenn", quietly = TRUE)) {
  sets <- lapply(setNames(rep_cohorts, rep_cohorts),
                 function(co) stats$kmer[stats[[paste0(co, "_pass")]]])
  ggsave(file.path(plots_dir, "cohort_overlap_venn.pdf"), ggvenn::ggvenn(sets), width = 7, height = 6)
}

# =============================================================================
# 8. RUN SUMMARY
# =============================================================================
n_seq <- if (file.exists(file.path(res_dir, paste0(prefix, "_sequences.tsv"))))
  nrow(read_tsv(file.path(res_dir, paste0(prefix, "_sequences.tsv")), show_col_types = FALSE)) else 0L
summary_lines <- c(
  sprintf("SURFR2 run summary - %s", format(Sys.time(), "%F %T")),
  if (any(mode_of == "descriptive"))
    "*** EXPLORATORY: at least one cohort had no statistical test (see modes below) ***",
  sprintf("project: %s", cfg$project),
  sprintf("analysis: %s", if (analysis == "dea") "differential expression" else paste0(CASE, "-specific k-mers")),
  sprintf("comparison: %s vs %s (log2FC = log2(%s / %s))", CASE, CTRL, CASE, CTRL),
  sprintf("normalisation: %s; library size = QC-passed reads (CPM = counts per million reads)",
          cfg$normalization$method),
  sprintf("cohorts (mode): %s", paste(sprintf("%s (%s)", rep_cohorts, mode_of[rep_cohorts]), collapse = ", ")),
  if (length(fdr_cols)) sprintf("test: edgeR quasi-likelihood F-test, BH FDR <= %g per cohort", st$fdr)
  else "test: none",
  if (analysis == "dea")
    sprintf("thresholds: |log2FC| >= %g, mean CPM in the higher condition >= %g",
            dea_cfg$min_abs_log2fc, dea_cfg$min_mean_cpm),
  if (analysis == "specific")
    sprintf("external controls: %s", if (length(ext_groups)) paste(sub("^external\\|", "", ext_groups), collapse = ", ") else "none"),
  if ("artifact_reads" %in% names(samples) && any(!is.na(samples$artifact_reads)))
    sprintf("artifact removal: %s reads removed (%.2f%% of QC-passed reads; per sample %.2f-%.2f%%)",
            format(sum(samples$artifact_reads, na.rm = TRUE), big.mark = ","),
            100 * sum(samples$artifact_reads, na.rm = TRUE) / sum(samples$qc_reads + samples$artifact_reads, na.rm = TRUE),
            min(100 * samples$artifact_reads / (samples$qc_reads + samples$artifact_reads), na.rm = TRUE),
            max(100 * samples$artifact_reads / (samples$qc_reads + samples$artifact_reads), na.rm = TRUE))
  else "artifact removal: off",
  sprintf("candidate k-mers (condition-blind pre-filter: >= %d reads in >= %d samples): %d",
          cfg$candidates$min_count, design$candidates_min_samples, n_cand),
  if (analysis == "dea")
    sprintf("differentially expressed k-mers: %d (higher in %s: %d, higher in %s: %d)",
            nrow(final), CASE, sum(final$higher_in == CASE), CTRL, sum(final$higher_in == CTRL))
  else sprintf("%s k-mers: %d", label, nrow(final)),
  sprintf("merged sequences: %d", n_seq),
  "", "warnings:", if (length(warnings_log)) paste0("  - ", warnings_log) else "  none"
)
writeLines(summary_lines, file.path(res_dir, "run_summary.txt"))
msg("done - see %s", file.path(res_dir, "run_summary.txt"))
