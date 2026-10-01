#!/usr/bin/env Rscript
# =============================================================================
# SURFR2 - normalisation and case-specific k-mer filtering
#
#  1. Normalisation -> effective library size L_i per sample; CPM = count / L_i * 1e6
#       cpm          : L_i = total k-mers counted by KMC in sample i
#       median_ratio : DESeq2 median-of-ratios size factors on k-mers present in all
#                      samples, rescaled to L_i = s_i * geomean(total k-mers) so that
#                      CPM thresholds keep a comparable meaning under both methods
#  2. Streams over samples, accumulating per-k-mer, per-group statistics
#     (sum of CPM, number of samples detecting the k-mer, max CPM). Groups are
#     <cohort>|case and <cohort>|control, and external|<condition> for external
#     controls pooled over cohorts.
#  3. Per compared cohort, a k-mer passes if
#       detected (CPM >= min_cpm_case) in >= need_case case samples
#       mean_cpm_case >= min_mean_cpm_case
#       detected (CPM >= max_cpm_control) in <= allow_control control samples
#       log2((mean_case + pc) / (mean_ctrl + pc)) >= min_log2fc
#     and, in 'test' mode, edgeR quasi-likelihood FDR <= statistics.fdr.
#     Cohorts with a single sample per group (or no residual df) run in 'descriptive'
#     mode: same filters, no test, results flagged as exploratory.
#     External controls: detected in <= allow samples. Thresholds are resolved per
#     cohort by the config validator (resolved_config.json -> design).
#     Final set = passes in ALL compared cohorts and ALL external controls.
#  4. Count/CPM matrices of final k-mers, dekupl-mergeTags, QC plots.
#
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

cfg <- fromJSON(file.path(run_dir, "resolved_config.json"), simplifyVector = TRUE)
f   <- cfg$filters
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
samples <- samples |> left_join(select(lib, sample_id, qc_reads, total_kmers, unique_kmers),
                                by = "sample_id")
samples$block[is.na(samples$block)] <- ""

count_file <- function(sub, id) file.path(mdir, "counts", sub, paste0(id, ".tsv.gz"))
read_counts <- function(sub, id) {
  read_tsv(count_file(sub, id), col_names = c("kmer", "count"), show_col_types = FALSE,
           col_types = cols(kmer = col_character(), count = col_double()), progress = FALSE)
}
geomean <- function(x) exp(mean(log(x)))

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
# 1. NORMALISATION
# =============================================================================
msg("normalisation: %s", cfg$normalization$method)

if (cfg$normalization$method == "cpm") {
  samples$size_factor <- samples$total_kmers / geomean(samples$total_kmers)
  samples$eff_libsize <- samples$total_kmers
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
  samples$eff_libsize <- samples$size_factor * geomean(samples$total_kmers)
  msg("median-of-ratios on %d of %d reference k-mers", length(ref_kmers), n_ref_all)
  r <- cor(log(samples$size_factor), log(samples$total_kmers))
  msg("cor(log size factor, log total k-mers) = %.3f", r)
  if (r < 0.5) warn(paste0("size factors correlate weakly with depth (r=%.2f): strong composition ",
                           "differences between samples; inspect plots/normalisation.pdf"), r)
  rm(M, logM)
}

write_tsv(samples |> select(sample_id, cohort, condition, role, qc_reads, total_kmers,
                            unique_kmers, size_factor, eff_libsize),
          file.path(res_dir, "normalisation.tsv"))

# =============================================================================
# 2. GROUPS AND PRE-FILTER CONSISTENCY
# =============================================================================
rep_cohorts <- cfg$replicate_cohorts
samples$group <- NA_character_
in_rep <- samples$cohort %in% rep_cohorts & samples$role %in% c("case", "control")
samples$group[in_rep] <- paste(samples$cohort[in_rep], samples$role[in_rep], sep = "|")
is_ext <- samples$role == "external"
samples$group[is_ext] <- paste("external", samples$condition[is_ext], sep = "|")

# Threshold defining "expressed" per sample: case uses min_cpm_case, others max_cpm_control
samples$expr_thr <- ifelse(samples$role == "case", f$min_cpm_case, f$max_cpm_control)

groups <- sort(unique(na.omit(samples$group)))
n_in_group <- table(samples$group)[groups]
msg("groups: %s", paste(sprintf("%s(n=%d)", groups, as.integer(n_in_group)), collapse = ", "))

# The raw-count pre-filter (>= min_count reads) is only guaranteed to be looser than
# the final CPM filter if min_count <= min_cpm_case * L_i / 1e6 for every case sample.
case_rep <- samples$role == "case" & samples$cohort %in% rep_cohorts
count_equiv <- f$min_cpm_case * samples$eff_libsize[case_rep] / 1e6
n_viol <- sum(cfg$candidates$min_count > count_equiv + 1e-9)
if (n_viol > 0) {
  warn(paste0("candidates.min_count=%d exceeds min_cpm_case*libsize/1e6 in %d/%d case samples ",
              "(smallest equivalent count %.1f). In those samples, k-mers at CPM >= min_cpm_case ",
              "but < %d reads were not counted by the pre-filter, so case detection may be ",
              "underestimated for k-mers near threshold. Lower candidates.min_count to <= %d."),
       cfg$candidates$min_count, n_viol, sum(case_rep), min(count_equiv),
       cfg$candidates$min_count, max(1L, floor(min(count_equiv))))
}

# =============================================================================
# 3. STREAMING ACCUMULATION
# =============================================================================
cand <- read_tsv(file.path(mdir, "candidates.tsv"), col_names = c("kmer", "n_samples_prefilter"),
                 show_col_types = FALSE, col_types = cols(kmer = col_character(), n_samples_prefilter = col_integer()))
n_cand <- nrow(cand)
msg("%d candidate k-mers; accumulating over %d grouped samples", n_cand, sum(!is.na(samples$group)))

G <- length(groups)
sum_cpm <- matrix(0,  nrow = n_cand, ncol = G, dimnames = list(NULL, groups))
max_cpm <- matrix(0,  nrow = n_cand, ncol = G, dimnames = list(NULL, groups))
n_expr  <- matrix(0L, nrow = n_cand, ncol = G, dimnames = list(NULL, groups))

grouped <- which(!is.na(samples$group))
for (k in seq_along(grouped)) {
  s <- samples[grouped[k], ]
  d <- read_counts("candidates", s$sample_id)
  idx <- match(d$kmer, cand$kmer)
  if (anyNA(idx)) stop("sample ", s$sample_id, " has k-mers outside the candidate set")
  cpm <- d$count / s$eff_libsize * 1e6
  g <- s$group
  sum_cpm[idx, g] <- sum_cpm[idx, g] + cpm
  max_cpm[idx, g] <- pmax(max_cpm[idx, g], cpm)
  n_expr[idx, g]  <- n_expr[idx, g] + as.integer(cpm >= s$expr_thr - 1e-12)
  if (k %% 50 == 0) msg("  %d / %d samples", k, length(grouped))
}

# =============================================================================
# 4. PER-COHORT FILTERS AND STATISTICAL TEST
# =============================================================================
design <- cfg$design
st <- cfg$statistics
stats <- cand
pass_all <- rep(TRUE, n_cand)
pc <- f$pseudocount_cpm

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
  msg("  edgeR: %d of %d profiles at FDR <= %g (any direction; QL prior df %.1f)",
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
  pass <- n_expr[, gc] >= dz$need_case &
          mean_c >= f$min_mean_cpm_case &
          n_expr[, gn] <= dz$allow_control &
          lfc >= f$min_log2fc
  msg("%s [%s]: case n=%d (need >= %d detected), control n=%d (allow <= %d detected): %d pass filters",
      co, dz$mode, nc, dz$need_case, nn, dz$allow_control, sum(pass))

  stats[[paste0(co, "_n_detected_case")]]    <- n_expr[, gc]
  stats[[paste0(co, "_n_detected_control")]] <- n_expr[, gn]
  stats[[paste0(co, "_mean_cpm_case")]]      <- mean_c
  stats[[paste0(co, "_mean_cpm_control")]]   <- mean_n
  stats[[paste0(co, "_max_cpm_control")]]    <- max_cpm[, gn]
  stats[[paste0(co, "_log2fc")]]             <- lfc

  if (dz$mode == "test") {
    ss <- samples[samples$cohort == co & samples$role %in% c("case", "control"), ]
    et <- run_edger(ss)
    stats[[paste0(co, "_edger_logFC")]] <- et$logFC
    stats[[paste0(co, "_pvalue")]]      <- et$pvalue
    stats[[paste0(co, "_fdr")]]         <- et$fdr
    pass <- pass & !is.na(et$fdr) & et$fdr <= st$fdr
    msg("%s: %d pass filters AND FDR <= %g", co, sum(pass), st$fdr)
  } else {
    warn("cohort '%s' ran in DESCRIPTIVE mode (case n=%d, control n=%d): no statistical test; results are exploratory",
         co, nc, nn)
  }
  mode_of[co] <- dz$mode
  stats[[paste0(co, "_pass")]] <- pass
  pass_all <- pass_all & pass
}

ext_groups <- grep("^external\\|", groups, value = TRUE)
for (ge in ext_groups) {
  cond <- sub("^external\\|", "", ge)
  ne <- as.integer(n_in_group[ge])
  allow_ext <- design$external[[cond]]$allow
  pass <- n_expr[, ge] <= allow_ext
  lab <- paste0("ext_", cond)
  msg("%s: n=%d (allow <= %d detected): %d candidates pass", ge, ne, allow_ext, sum(pass))
  stats[[paste0(lab, "_mean_cpm")]]   <- sum_cpm[, ge] / ne
  stats[[paste0(lab, "_max_cpm")]]    <- max_cpm[, ge]
  stats[[paste0(lab, "_n_detected")]] <- n_expr[, ge]
  stats[[paste0(lab, "_pass")]]       <- pass
  pass_all <- pass_all & pass
}
stats$pass_all <- pass_all
rm(sum_cpm, max_cpm, n_expr); invisible(gc())

lfc_cols  <- paste0(rep_cohorts, "_log2fc")
mean_cols <- paste0(rep_cohorts, "_mean_cpm_case")
fdr_cols  <- intersect(paste0(rep_cohorts, "_fdr"), names(stats))
stats$min_log2fc           <- do.call(pmin, unname(as.list(stats[lfc_cols])))
stats$min_mean_cpm_case    <- do.call(pmin, unname(as.list(stats[mean_cols])))
stats$max_mean_cpm_control <- do.call(pmax, unname(as.list(stats[paste0(rep_cohorts, "_mean_cpm_control")])))
if (length(fdr_cols)) stats$max_fdr <- do.call(pmax, unname(as.list(stats[fdr_cols])))

write_tsv(stats, file.path(res_dir, "candidate_kmer_stats.tsv.gz"))
final <- stats |> filter(pass_all) |>
  arrange(if (length(fdr_cols)) max_fdr else 0, desc(min_log2fc), desc(min_mean_cpm_case))
write_tsv(final, file.path(res_dir, "case_specific_kmers.tsv"))
msg("%d case-specific k-mers pass in all of: %s%s", nrow(final), paste(rep_cohorts, collapse = ", "),
    if (length(ext_groups)) paste0(" + ", paste(ext_groups, collapse = ", ")) else "")

# =============================================================================
# 5. MATRICES FOR FINAL K-MERS (all samples, incl. non-replicate cohorts)
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
  write_tsv(as_tibble(cnt, rownames = "kmer"),     file.path(res_dir, "case_specific_kmers_counts.tsv"))
  write_tsv(as_tibble(cpm_mat, rownames = "kmer"), file.path(res_dir, "case_specific_kmers_cpm.tsv"))
} else {
  warn("no k-mer passed all filters; matrices, mergeTags and k-mer plots are skipped")
}

# =============================================================================
# 6. dekupl-mergeTags: assemble overlapping / offset k-mers into sequences
# =============================================================================
if (isTRUE(cfg$mergetags$enabled) && nrow(final) > 0) {
  mt_in  <- file.path(res_dir, "mergetags_input.tsv")
  mt_out <- file.path(res_dir, "case_specific_sequences.tsv")
  # mergeTags hard-codes DE-kupl's layout (see mergeTags.c): value column 1 is read as a
  # p-value (the LOWEST becomes the contig representative) and value column 4 as log2FC
  # (> 0 = 'up'); it reads column 4 unconditionally, so >= 4 value columns are mandatory.
  # We supply a rank (1 = best k-mer) in the p-value slot. SURFR1 passed raw TCGA cancer
  # counts there, so each contig was represented by its LEAST abundant k-mer.
  mt_tab <- final |>
    mutate(rank = row_number()) |>   # 'final' is sorted best-first (max FDR, then effect size)
    select(tag = kmer, rank, min_mean_cpm_case, max_mean_cpm_control, min_log2fc)
  write_tsv(mt_tab, mt_in)
  mt_args <- c("-k", cfg$kmer$k, "-m", cfg$mergetags$min_overlap)
  # '-n' means UNSTRANDED merging (mergeTags.c: case 'n': stranded = 0). SURFR1 used it on
  # stranded small-RNA data, allowing k-mers to merge with reverse complements.
  if (isTRUE(cfg$kmer$canonical)) mt_args <- c(mt_args, "-n")
  status <- system2(cfg$tools$mergetags, c(mt_args, mt_in), stdout = mt_out)
  if (!identical(as.integer(status), 0L)) stop("dekupl-mergeTags failed with status ", status)
  contigs <- read_tsv(mt_out, show_col_types = FALSE)
  msg("mergeTags: %d k-mers -> %d sequences", nrow(final), nrow(contigs))

  # Per-sequence summary of the test. FDR is controlled at the k-mer level; overlapping
  # k-mers of one molecule are strongly correlated, so their p-values are combined with
  # the Cauchy combination test (ACAT; Liu & Xie 2020, JASA), valid under dependence.
  if (length(fdr_cols) && nrow(contigs)) {
    k <- cfg$kmer$k
    members <- lapply(contigs$contig, function(sq) {
      n <- nchar(sq) - k + 1
      intersect(substring(sq, seq_len(n), seq_len(n) + k - 1), final$kmer)
    })
    for (co in sub("_fdr$", "", fdr_cols)) {
      pv <- setNames(final[[paste0(co, "_pvalue")]], final$kmer)
      fd <- setNames(final[[paste0(co, "_fdr")]], final$kmer)
      contigs[[paste0(co, "_acat_pvalue")]] <- vapply(members, function(m) acat(pv[m]), numeric(1))
      contigs[[paste0(co, "_min_kmer_fdr")]] <- vapply(members, function(m) min(fd[m]), numeric(1))
    }
    write_tsv(contigs, mt_out)
  }
}

# =============================================================================
# 7. QC PLOTS
# =============================================================================
p_lib <- ggplot(samples, aes(x = interaction(condition, cohort, sep = "\n"), y = total_kmers)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  scale_y_log10() + theme_bw() +
  labs(x = NULL, y = "total k-mers (log10)", title = "Library sizes")
ggsave(file.path(plots_dir, "library_sizes.pdf"), p_lib, width = 8, height = 5)

if (cfg$normalization$method == "median_ratio") {
  p_sf <- ggplot(samples, aes(total_kmers / geomean(total_kmers), size_factor, colour = condition)) +
    geom_abline(linetype = "dashed") + geom_point() + scale_x_log10() + scale_y_log10() + theme_bw() +
    labs(x = "relative library size", y = "median-of-ratios size factor",
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
  role_of <- setNames(samples$role, grp)
  vdf$role <- role_of[vdf$group]
  p_v <- ggplot(vdf, aes(group, log10(mean_cpm + pc), fill = role)) +
    geom_violin(scale = "width") +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.5) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = sprintf("log10(mean CPM + %g)", pc),
         title = sprintf("%d case-specific k-mers", nrow(final)))
  ggsave(file.path(plots_dir, "case_specific_kmers_violin.pdf"), p_v, width = 9, height = 5)
}

if (length(rep_cohorts) %in% 2:4 && requireNamespace("ggvenn", quietly = TRUE)) {
  sets <- lapply(setNames(rep_cohorts, rep_cohorts),
                 function(co) stats$kmer[stats[[paste0(co, "_pass")]]])
  ggsave(file.path(plots_dir, "cohort_overlap_venn.pdf"), ggvenn::ggvenn(sets), width = 7, height = 6)
}

# =============================================================================
# 8. RUN SUMMARY
# =============================================================================
summary_lines <- c(
  sprintf("SURFR2 run summary - %s", format(Sys.time(), "%F %T")),
  sprintf("project: %s   case: %s   control: %s", cfg$project, cfg$comparison$case, cfg$comparison$control),
  sprintf("normalisation: %s", cfg$normalization$method),
  if (any(mode_of == "descriptive"))
    "*** EXPLORATORY: at least one cohort had no statistical test (see modes below) ***",
  sprintf("cohorts (mode): %s", paste(sprintf("%s (%s)", rep_cohorts, mode_of[rep_cohorts]), collapse = ", ")),
  if (length(fdr_cols)) sprintf("test: edgeR quasi-likelihood F-test, BH FDR <= %g per cohort", st$fdr)
  else "test: none",
  sprintf("external controls: %s", if (length(ext_groups)) paste(ext_groups, collapse = ", ") else "none"),
  sprintf("candidate k-mers (condition-blind pre-filter: >= %d reads in >= %d samples): %d",
          cfg$candidates$min_count, design$candidates_min_samples, n_cand),
  sprintf("case-specific k-mers (all filters): %d", nrow(final)),
  "", "warnings:", if (length(warnings_log)) paste0("  - ", warnings_log) else "  none"
)
writeLines(summary_lines, file.path(res_dir, "run_summary.txt"))
msg("done - see %s", file.path(res_dir, "run_summary.txt"))
