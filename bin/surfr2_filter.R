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
#     (sum of CPM, number of samples above threshold, max CPM). Groups are
#     <cohort>|case and <cohort>|control for replicate cohorts, and
#     external|<condition> for external controls pooled over cohorts.
#  3. Per replicate cohort, a k-mer passes if
#       prevalence_case >= min_prevalence_case   (CPM >= min_cpm_case)
#       mean_cpm_case   >= min_mean_cpm_case
#       prevalence_ctrl <= max_prevalence_control (CPM >= max_cpm_control)
#       log2((mean_case + pc) / (mean_ctrl + pc)) >= min_log2fc
#     and for every external control: prevalence <= max_prevalence_external.
#     Final set = passes in ALL replicate cohorts and ALL external controls.
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

samples <- read_tsv(file.path(run_dir, "samples.tsv"), show_col_types = FALSE,
                    col_types = cols(.default = col_character(), index = col_integer()))
lib <- read_tsv(file.path(mdir, "library_sizes.tsv"), show_col_types = FALSE,
                col_types = cols(.default = col_character(), qc_reads = col_double(),
                                 total_kmers = col_double(), unique_kmers = col_double()))
stopifnot(setequal(samples$sample_id, lib$sample_id), !anyDuplicated(lib$sample_id))
samples <- samples |> left_join(select(lib, sample_id, qc_reads, total_kmers, unique_kmers),
                                by = "sample_id")

count_file <- function(sub, id) file.path(mdir, "counts", sub, paste0(id, ".tsv.gz"))
read_counts <- function(sub, id) {
  read_tsv(count_file(sub, id), col_names = c("kmer", "count"), show_col_types = FALSE,
           col_types = cols(kmer = col_character(), count = col_double()), progress = FALSE)
}
geomean <- function(x) exp(mean(log(x)))

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
              "but < %d reads were not counted by the pre-filter, so prevalence_case may be ",
              "underestimated for k-mers near threshold. Lower candidates.min_count to <= %d."),
       cfg$candidates$min_count, n_viol, sum(case_rep), min(count_equiv),
       cfg$candidates$min_count, max(1L, floor(min(count_equiv))))
}

# =============================================================================
# 3. STREAMING ACCUMULATION
# =============================================================================
cand <- read_tsv(file.path(mdir, "candidates.tsv"), col_names = c("kmer", "n_case_prefilter"),
                 show_col_types = FALSE, col_types = cols(kmer = col_character(), n_case_prefilter = col_integer()))
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
# 4. PER-COHORT STATISTICS AND FILTERS
#    Prevalence thresholds are compared as integer sample counts, using the same
#    rounding as the config validator, to avoid floating-point edge cases.
# =============================================================================
stats <- cand
pass_all <- rep(TRUE, n_cand)
pc <- f$pseudocount_cpm

for (co in rep_cohorts) {
  gc <- paste(co, "case", sep = "|"); gn <- paste(co, "control", sep = "|")
  nc <- as.integer(n_in_group[gc]); nn <- as.integer(n_in_group[gn])
  mean_c <- sum_cpm[, gc] / nc
  mean_n <- sum_cpm[, gn] / nn
  lfc <- log2((mean_c + pc) / (mean_n + pc))
  need_case <- max(1L, ceiling(f$min_prevalence_case * nc - 1e-9))
  allow_ctrl <- floor(f$max_prevalence_control * nn + 1e-9)
  pass <- n_expr[, gc] >= need_case &
          mean_c >= f$min_mean_cpm_case &
          n_expr[, gn] <= allow_ctrl &
          lfc >= f$min_log2fc
  msg("%s: case n=%d (need >= %d expressing), control n=%d (allow <= %d): %d pass",
      co, nc, need_case, nn, allow_ctrl, sum(pass))
  stats[[paste0(co, "_mean_cpm_case")]]   <- mean_c
  stats[[paste0(co, "_prev_case")]]       <- n_expr[, gc] / nc
  stats[[paste0(co, "_mean_cpm_control")]] <- mean_n
  stats[[paste0(co, "_max_cpm_control")]]  <- max_cpm[, gn]
  stats[[paste0(co, "_prev_control")]]    <- n_expr[, gn] / nn
  stats[[paste0(co, "_log2fc")]]          <- lfc
  stats[[paste0(co, "_pass")]]            <- pass
  pass_all <- pass_all & pass
}

ext_groups <- grep("^external\\|", groups, value = TRUE)
for (ge in ext_groups) {
  ne <- as.integer(n_in_group[ge])
  allow_ext <- floor(f$max_prevalence_external * ne + 1e-9)
  pass <- n_expr[, ge] <= allow_ext
  lab <- sub("^external\\|", "ext_", ge)
  msg("%s: n=%d (allow <= %d expressing): %d candidates pass", ge, ne, allow_ext, sum(pass))
  stats[[paste0(lab, "_mean_cpm")]] <- sum_cpm[, ge] / ne
  stats[[paste0(lab, "_max_cpm")]]  <- max_cpm[, ge]
  stats[[paste0(lab, "_prev")]]     <- n_expr[, ge] / ne
  stats[[paste0(lab, "_pass")]]     <- pass
  pass_all <- pass_all & pass
}
stats$pass_all <- pass_all
rm(sum_cpm, max_cpm, n_expr); invisible(gc())

lfc_cols  <- paste0(rep_cohorts, "_log2fc")
mean_cols <- paste0(rep_cohorts, "_mean_cpm_case")
stats$min_log2fc        <- do.call(pmin, unname(as.list(stats[lfc_cols])))
stats$min_mean_cpm_case <- do.call(pmin, unname(as.list(stats[mean_cols])))
stats$max_mean_cpm_control <- do.call(pmax, unname(as.list(stats[paste0(rep_cohorts, "_mean_cpm_control")])))

write_tsv(stats, file.path(res_dir, "candidate_kmer_stats.tsv.gz"))
final <- stats |> filter(pass_all) |> arrange(desc(min_log2fc), desc(min_mean_cpm_case))
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
    mutate(rank = row_number()) |>            # 'final' is sorted best-first above
    select(tag = kmer, rank, min_mean_cpm_case, max_mean_cpm_control, min_log2fc)
  write_tsv(mt_tab, mt_in)
  mt_args <- c("-k", cfg$kmer$k, "-m", cfg$mergetags$min_overlap)
  # '-n' means UNSTRANDED merging (mergeTags.c: case 'n': stranded = 0). SURFR1 used it on
  # stranded small-RNA data, allowing k-mers to merge with reverse complements.
  if (isTRUE(cfg$kmer$canonical)) mt_args <- c(mt_args, "-n")
  status <- system2(cfg$tools$mergetags, c(mt_args, mt_in), stdout = mt_out)
  if (!identical(as.integer(status), 0L)) stop("dekupl-mergeTags failed with status ", status)
  n_seq <- nrow(read_tsv(mt_out, show_col_types = FALSE))
  msg("mergeTags: %d k-mers -> %d sequences", nrow(final), n_seq)
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
  sprintf("replicate cohorts: %s", paste(rep_cohorts, collapse = ", ")),
  sprintf("external controls: %s", if (length(ext_groups)) paste(ext_groups, collapse = ", ") else "none"),
  sprintf("candidate k-mers (pre-filter): %d", n_cand),
  sprintf("case-specific k-mers (all filters): %d", nrow(final)),
  "", "warnings:", if (length(warnings_log)) paste0("  - ", warnings_log) else "  none"
)
writeLines(summary_lines, file.path(res_dir, "run_summary.txt"))
msg("done - see %s", file.path(res_dir, "run_summary.txt"))
