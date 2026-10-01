# =============================================================================
# SURFR2 - shared plotting library
#
# Used by bin/surfr2_report.R (static PDF + self-contained interactive HTML) and by
# app/app.R (Shiny). One implementation, so the pipeline output and the app always
# show the same thing.
#
# Everything is read from <outdir>/results/, which is therefore self-contained:
# copy that folder to a laptop and the Shiny app works without the rest of the run.
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Loading
# -----------------------------------------------------------------------------
surfr2_load_results <- function(res_dir) {
  need <- c("resolved_config.json", "candidate_kmer_stats.tsv.gz", "normalisation.tsv",
            "case_specific_kmers.tsv")
  miss <- need[!file.exists(file.path(res_dir, need))]
  if (length(miss)) stop("not a SURFR2 results directory (missing: ", paste(miss, collapse = ", "), "): ", res_dir)

  cfg  <- fromJSON(file.path(res_dir, "resolved_config.json"), simplifyVector = TRUE)
  rd   <- function(f, ...) read_tsv(file.path(res_dir, f), show_col_types = FALSE, progress = FALSE, ...)
  stats <- rd("candidate_kmer_stats.tsv.gz")
  final <- rd("case_specific_kmers.tsv")
  norm  <- rd("normalisation.tsv", col_types = cols(.default = col_character(),
               qc_reads = col_double(), total_kmers = col_double(), unique_kmers = col_double(),
               size_factor = col_double(), eff_libsize = col_double()))

  cnt <- cpm <- NULL
  if (nrow(final) > 0) {
    to_mat <- function(f) {
      d <- rd(f)
      m <- as.matrix(d[, -1]); rownames(m) <- d$kmer
      m[, norm$sample_id, drop = FALSE]               # fixed sample order
    }
    cnt <- to_mat("case_specific_kmers_counts.tsv")
    cpm <- to_mat("case_specific_kmers_cpm.tsv")
  }
  seq_file <- file.path(res_dir, "case_specific_sequences.tsv")
  contigs <- if (file.exists(seq_file)) rd("case_specific_sequences.tsv") else NULL

  list(cfg = cfg, stats = stats, final = final, norm = norm, cnt = cnt, cpm = cpm,
       contigs = contigs, contig_of = surfr2_kmer_to_contig(final$kmer, contigs, cfg$kmer$k))
}

# mergeTags reports each contig with one representative k-mer only; recover full
# membership by decomposing every contig into its k-mers (stranded: no rev-comp).
surfr2_kmer_to_contig <- function(kmers, contigs, k) {
  out <- setNames(rep(NA_character_, length(kmers)), kmers)
  if (is.null(contigs) || !nrow(contigs) || !length(kmers)) return(out)
  for (s in contigs$contig) {
    n <- nchar(s) - k + 1
    if (n < 1) next
    km <- substring(s, seq_len(n), seq_len(n) + k - 1)
    hit <- intersect(km, kmers)
    out[hit] <- s
  }
  out
}

# -----------------------------------------------------------------------------
# Hover text for passing k-mers: per-condition counts (+ per-sample if few samples)
# -----------------------------------------------------------------------------
surfr2_hover_text <- function(res, cohort, max_samples_listed = 12) {
  f <- res$final
  if (!nrow(f)) return(character())
  norm <- res$norm
  grp  <- paste(norm$cohort, norm$condition)
  lfc  <- f[[paste0(cohort, "_log2fc")]]
  contig <- res$contig_of[f$kmer]

  # one line per (cohort, condition) group, for each passing k-mer
  lines <- lapply(split(seq_len(nrow(norm)), grp), function(j) {
    cn <- res$cnt[f$kmer, j, drop = FALSE]
    cp <- res$cpm[f$kmer, j, drop = FALSE]
    list(cn = rowSums(cn), cp = rowMeans(cp), det = rowSums(cn > 0), n = length(j))
  })
  cond_block <- do.call(paste, c(lapply(names(lines), function(g) {
    l <- lines[[g]]
    sprintf("%s: counts=%s | mean CPM=%.2f | detected %d/%d",
            g, format(l$cn, big.mark = ",", trim = TRUE), l$cp, l$det, l$n)
  }), sep = "<br>"))

  sample_block <- rep("", nrow(f))
  if (nrow(norm) <= max_samples_listed) {
    cn <- res$cnt[f$kmer, , drop = FALSE]
    lab <- sprintf("%s (%s)", norm$sample_id, norm$condition)
    sample_block <- apply(cn, 1, function(v)
      paste0("<br><br><b>per sample (raw counts)</b><br>",
             paste(sprintf("%s: %s", lab, format(v, big.mark = ",", trim = TRUE)), collapse = "<br>")))
  }

  paste0("<b>", f$kmer, "</b>",
         "<br>sequence: ", ifelse(is.na(contig), "(not merged)", contig),
         sprintf("<br>log2FC (%s): %.2f", cohort, lfc),
         if (paste0(cohort, "_fdr") %in% names(f))
           sprintf("<br>FDR (%s, edgeR QL): %.3g", cohort, f[[paste0(cohort, "_fdr")]])
         else "<br>no statistical test (descriptive mode)",
         "<br><br><b>per condition</b><br>", cond_block,
         sample_block)
}

# -----------------------------------------------------------------------------
# Scatter data for one cohort
# -----------------------------------------------------------------------------
surfr2_scatter_data <- function(res, cohort, nbins = 150) {
  pc <- res$cfg$filters$pseudocount_cpm
  mc <- paste0(cohort, "_mean_cpm_case"); mn <- paste0(cohort, "_mean_cpm_control")
  if (!all(c(mc, mn) %in% names(res$stats))) stop("cohort not in results: ", cohort)
  x <- log10(res$stats[[mn]] + pc); y <- log10(res$stats[[mc]] + pc)

  # Exact 2-D histogram of ALL candidates (compact background; no subsampling)
  lim <- range(c(x, y), finite = TRUE)
  br  <- seq(lim[1], lim[2] + 1e-9, length.out = nbins + 1)
  mid <- (head(br, -1) + tail(br, -1)) / 2
  ix <- findInterval(x, br, rightmost.closed = TRUE)
  iy <- findInterval(y, br, rightmost.closed = TRUE)
  z  <- matrix(0, nbins, nbins)
  tab <- table(factor(iy, levels = 1:nbins), factor(ix, levels = 1:nbins))
  z[] <- as.numeric(tab)                         # rows = y bins, cols = x bins
  z[z == 0] <- NA

  pass <- res$stats$kmer %in% res$final$kmer
  fin  <- match(res$final$kmer, res$stats$kmer)

  list(cohort = cohort, pc = pc, lim = lim, mid = mid, z = z,
       n_all = length(x), n_pass = sum(pass),
       fx = x[fin], fy = y[fin], fkmer = res$final$kmer,
       y_min = log10(res$cfg$filters$min_mean_cpm_case + pc),
       lfc_offset = res$cfg$filters$min_log2fc * log10(2))
}

# -----------------------------------------------------------------------------
# Interactive plot (plotly). 'highlight' = k-mers to mark (e.g. a sequence search).
# -----------------------------------------------------------------------------
surfr2_scatter_plotly <- function(res, cohort, highlight = character(), source = "scatter") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("R package 'plotly' is required")
  d <- surfr2_scatter_data(res, cohort)
  cfg <- res$cfg
  hover <- surfr2_hover_text(res, cohort)

  p <- plotly::plot_ly(source = source) |>
    plotly::add_heatmap(x = d$mid, y = d$mid, z = log10(d$z), hoverinfo = "skip",
                        colorscale = list(c(0, "#c4c4c4"), c(1, "#3d3d3d")),
                        showscale = FALSE, name = "all candidates")
  if (d$n_pass > 0) {
    p <- p |> plotly::add_trace(
      type = "scattergl", mode = "markers", x = d$fx, y = d$fy,
      text = hover, hoverinfo = "text", customdata = d$fkmer,
      marker = list(color = "#febf38", size = 7, line = list(color = "#7a5a00", width = 0.6)),
      name = sprintf("case-specific (%d)", d$n_pass))
    hl <- which(d$fkmer %in% highlight)
    if (length(hl)) {
      p <- p |> plotly::add_trace(
        type = "scattergl", mode = "markers", x = d$fx[hl], y = d$fy[hl],
        text = hover[hl], hoverinfo = "text", customdata = d$fkmer[hl],
        marker = list(color = "#d62728", size = 11, symbol = "diamond",
                      line = list(color = "black", width = 1)),
        name = sprintf("search hits (%d)", length(hl)))
    }
  }
  lo <- d$lim[1]; hi <- d$lim[2]
  line <- function(x0, y0, x1, y1) list(type = "line", x0 = x0, y0 = y0, x1 = x1, y1 = y1,
                                        line = list(dash = "dash", width = 1, color = "#b22222"))
  plotly::layout(p,
    title = list(text = sprintf("%s: %s vs %s  (%s candidates, %d case-specific)",
                                cohort, cfg$comparison$case, cfg$comparison$control,
                                format(d$n_all, big.mark = ","), d$n_pass), font = list(size = 14)),
    xaxis = list(title = sprintf("log10(mean CPM %s + %g)", cfg$comparison$control, d$pc),
                 range = c(lo, hi), zeroline = FALSE),
    yaxis = list(title = sprintf("log10(mean CPM %s + %g)", cfg$comparison$case, d$pc),
                 range = c(lo, hi), zeroline = FALSE, scaleanchor = "x"),
    shapes = list(line(lo, d$y_min, hi, d$y_min),                       # min_mean_cpm_case
                  line(lo, lo + d$lfc_offset, hi - d$lfc_offset, hi)),  # min_log2fc (exact)
    hoverlabel = list(align = "left", font = list(family = "monospace", size = 11)),
    legend = list(orientation = "h", y = -0.15))
}

# -----------------------------------------------------------------------------
# Static plot (ggplot2) for publication
# -----------------------------------------------------------------------------
surfr2_scatter_ggplot <- function(res, cohort) {
  d <- surfr2_scatter_data(res, cohort)
  cfg <- res$cfg
  bg <- expand.grid(y = d$mid, x = d$mid)          # column-major matches z[] layout
  bg$n <- as.vector(d$z)
  bg <- bg[!is.na(bg$n), ]
  step <- diff(d$mid[1:2])
  g <- ggplot() +
    geom_tile(data = bg, aes(x, y, fill = log10(n)), width = step, height = step) +
    scale_fill_gradient(low = "#c4c4c4", high = "#3d3d3d", guide = "none") +
    geom_hline(yintercept = d$y_min, linetype = "dashed", colour = "#b22222") +
    geom_abline(intercept = d$lfc_offset, slope = 1, linetype = "dashed", colour = "#b22222")
  if (d$n_pass > 0)
    g <- g + geom_point(aes(x = d$fx, y = d$fy), shape = 21, size = 1.8,
                        fill = "#febf38", colour = "#7a5a00", stroke = 0.3)
  g + coord_equal(xlim = d$lim, ylim = d$lim) + theme_classic(base_size = 14) +
    labs(x = sprintf("log10(mean CPM %s + %g)", cfg$comparison$control, d$pc),
         y = sprintf("log10(mean CPM %s + %g)", cfg$comparison$case, d$pc),
         title = sprintf("%s: %d case-specific of %s candidate k-mers",
                         cohort, d$n_pass, format(d$n_all, big.mark = ",")))
}
