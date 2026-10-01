# =============================================================================
# SURFR2 - shared plotting library
#
# Used by bin/surfr2_report.R (static PDF + self-contained interactive HTML) and by
# app/app.R (Shiny). One implementation, so the pipeline output and the app always
# show the same thing.
#
# Plots per compared cohort:
#   scatter  mean CPM of the control condition (x) vs the case condition (y)
#   volcano  log2 fold change (x) vs -log10 p-value (y); cohorts with a test only
# All candidate k-mers form an exact binned background (not hoverable); passing k-mers
# are hoverable points, coloured by the condition they are higher in (DEA).
#
# Everything is read from <outdir>/results/, which is self-contained: copy that folder
# to a laptop and the Shiny app works without the rest of the run.
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
  if (!file.exists(file.path(res_dir, "resolved_config.json")))
    stop("not a SURFR2 results directory (no resolved_config.json): ", res_dir)
  cfg <- fromJSON(file.path(res_dir, "resolved_config.json"), simplifyVector = TRUE)
  prefix <- if (identical(cfg$analysis, "specific")) "specific" else "de"
  need <- c("kmer_stats.tsv.gz", "normalisation.tsv", paste0(prefix, "_kmers.tsv"))
  miss <- need[!file.exists(file.path(res_dir, need))]
  if (length(miss)) stop("incomplete SURFR2 results (missing: ", paste(miss, collapse = ", "), "): ", res_dir)

  rd <- function(f, ...) read_tsv(file.path(res_dir, f), show_col_types = FALSE, progress = FALSE, ...)
  stats <- rd("kmer_stats.tsv.gz")
  final <- rd(paste0(prefix, "_kmers.tsv"))
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
    cnt <- to_mat(paste0(prefix, "_kmers_counts.tsv"))
    cpm <- to_mat(paste0(prefix, "_kmers_cpm.tsv"))
  }
  seq_file <- paste0(prefix, "_sequences.tsv")
  contigs <- if (file.exists(file.path(res_dir, seq_file))) rd(seq_file) else NULL

  CASE <- cfg$comparison$case; CTRL <- cfg$comparison$control
  list(cfg = cfg, stats = stats, final = final, norm = norm, cnt = cnt, cpm = cpm,
       contigs = contigs, contig_of = surfr2_kmer_to_contig(final$kmer, contigs, cfg$kmer$k),
       analysis = if (identical(cfg$analysis, "specific")) "specific" else "dea",
       CASE = CASE, CTRL = CTRL, VS = paste0(CASE, "_vs_", CTRL),
       label = if (identical(cfg$analysis, "specific")) paste0(CASE, "-specific") else "differentially expressed")
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
    out[intersect(km, kmers)] <- s
  }
  out
}

surfr2_has_test <- function(res, cohort) paste0(cohort, "_pvalue") %in% names(res$stats)

# Fold-change column of a cohort: edgeR's estimate when tested, else the descriptive one
surfr2_fc_col <- function(res, cohort) {
  e <- paste0(cohort, "_edger_log2fc_", res$VS)
  if (e %in% names(res$stats)) e else paste0(cohort, "_log2fc_", res$VS)
}

# Group of each passing k-mer for colouring: the higher condition (DEA) or the label
surfr2_pass_group <- function(res) {
  if (res$analysis == "dea") paste("higher in", res$final$higher_in) else rep(res$label, nrow(res$final))
}

# -----------------------------------------------------------------------------
# Hover text for passing k-mers: per-condition counts (+ per-sample if few samples)
# -----------------------------------------------------------------------------
surfr2_hover_text <- function(res, cohort, max_samples_listed = 12) {
  f <- res$final
  if (!nrow(f)) return(character())
  norm <- res$norm
  grp  <- paste(norm$cohort, norm$condition)
  fc_col <- surfr2_fc_col(res, cohort)
  fc_lab <- if (grepl("_edger_", fc_col)) "log2FC, edgeR" else "log2FC"
  contig <- res$contig_of[f$kmer]

  lines <- lapply(split(seq_len(nrow(norm)), grp), function(j) {
    cn <- res$cnt[f$kmer, j, drop = FALSE]
    cp <- res$cpm[f$kmer, j, drop = FALSE]
    list(cn = rowSums(cn), cp = rowMeans(cp), det = rowSums(cn > 0), n = length(j))
  })
  cond_block <- do.call(paste, c(lapply(names(lines), function(g) {
    l <- lines[[g]]
    sprintf("%s: counts=%s | mean CPM=%.2f | found in %d/%d",
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
         if (res$analysis == "dea") paste0("<br>higher in: ", f$higher_in) else "",
         sprintf("<br>%s (%s / %s, %s): %.2f", fc_lab, res$CASE, res$CTRL, cohort, f[[fc_col]]),
         if (paste0(cohort, "_fdr") %in% names(f))
           sprintf("<br>FDR (%s, edgeR QL): %.3g", cohort, f[[paste0(cohort, "_fdr")]])
         else "<br>no statistical test (descriptive mode)",
         "<br><br><b>per condition</b><br>", cond_block,
         sample_block)
}

# -----------------------------------------------------------------------------
# Plot data: background histogram, passing points, threshold lines
# -----------------------------------------------------------------------------
bin2d <- function(x, y, xlim, ylim, nbins) {
  bx <- seq(xlim[1], xlim[2] + 1e-9, length.out = nbins + 1)
  by <- seq(ylim[1], ylim[2] + 1e-9, length.out = nbins + 1)
  tab <- table(factor(findInterval(y, by, rightmost.closed = TRUE), levels = 1:nbins),
               factor(findInterval(x, bx, rightmost.closed = TRUE), levels = 1:nbins))
  z <- matrix(as.numeric(tab), nbins, nbins)          # rows = y bins, cols = x bins
  z[z == 0] <- NA
  list(xmid = (head(bx, -1) + tail(bx, -1)) / 2, ymid = (head(by, -1) + tail(by, -1)) / 2, z = z,
       xstep = diff(bx[1:2]), ystep = diff(by[1:2]))
}

surfr2_plot_data <- function(res, cohort, type = c("scatter", "volcano"), nbins = 150) {
  type <- match.arg(type)
  cfg <- res$cfg; s <- res$stats
  pc <- cfg$filters$pseudocount_cpm
  fin <- match(res$final$kmer, s$kmer)
  seg <- function(x0, y0, x1, y1, kind = "threshold") data.frame(x0, y0, x1, y1, kind)

  if (type == "scatter") {
    mc <- paste0(cohort, "_mean_cpm_", res$CASE); mn <- paste0(cohort, "_mean_cpm_", res$CTRL)
    if (!all(c(mc, mn) %in% names(s))) stop("cohort not in results: ", cohort)
    x <- log10(s[[mn]] + pc); y <- log10(s[[mc]] + pc)
    lim <- range(c(x, y), finite = TRUE); xlim <- ylim <- lim
    if (res$analysis == "specific") {
      y_min <- log10(cfg$specific$min_mean_cpm_case + pc); off <- cfg$specific$min_log2fc * log10(2)
      lines <- rbind(seg(lim[1], y_min, lim[2], y_min),                    # min_mean_cpm_case
                     seg(lim[1], lim[1] + off, lim[2] - off, lim[2]))       # min_log2fc (exact)
    } else {
      lines <- seg(lim[1], lim[1], lim[2], lim[2], "reference")            # equal expression
    }
    xlab <- sprintf("log10(mean CPM %s + %g)", res$CTRL, pc)
    ylab <- sprintf("log10(mean CPM %s + %g)", res$CASE, pc)
    keep <- rep(TRUE, length(x))
  } else {
    if (!surfr2_has_test(res, cohort)) stop("no statistical test in cohort ", cohort, " (descriptive mode)")
    p <- s[[paste0(cohort, "_pvalue")]]; fdr <- s[[paste0(cohort, "_fdr")]]
    keep <- !is.na(p)
    x <- s[[surfr2_fc_col(res, cohort)]]
    y <- -log10(pmax(p, .Machine$double.xmin))
    xlim <- range(x[keep], finite = TRUE); xlim <- c(-1, 1) * max(abs(xlim))
    ylim <- c(0, max(y[keep], finite = TRUE) * 1.03)
    ok <- keep & fdr <= cfg$statistics$fdr
    lines <- NULL
    if (any(ok, na.rm = TRUE)) {     # p-value at the FDR cut-off: largest p still passing
      yt <- -log10(max(p[which(ok)]))
      lines <- seg(xlim[1], yt, xlim[2], yt)
    }
    if (res$analysis == "dea") {
      t <- cfg$dea$min_abs_log2fc
      lines <- rbind(lines, seg(-t, ylim[1], -t, ylim[2]), seg(t, ylim[1], t, ylim[2]))
    }
    xlab <- sprintf("log2 fold change (%s / %s)", res$CASE, res$CTRL)
    ylab <- "-log10 p-value (edgeR QL)"
  }
  bg <- bin2d(x[keep], y[keep], xlim, ylim, nbins)
  list(type = type, cohort = cohort, xlim = xlim, ylim = ylim, bg = bg, lines = lines,
       n_all = sum(keep), fx = x[fin], fy = y[fin], fkmer = res$final$kmer,
       fgroup = surfr2_pass_group(res), xlab = xlab, ylab = ylab,
       title = sprintf("%s: %s vs %s (%s %s k-mers, %d %s)", cohort, res$CASE, res$CTRL,
                       format(sum(keep), big.mark = ","), if (type == "volcano") "tested" else "candidate",
                       length(fin), res$label))
}

surfr2_group_colours <- function(res) {
  if (res$analysis == "dea")
    setNames(c("#e6550d", "#3182bd"), paste("higher in", c(res$CASE, res$CTRL)))
  else setNames("#febf38", res$label)
}

# -----------------------------------------------------------------------------
# Interactive plot (plotly). 'highlight' = k-mers to mark (e.g. a sequence search).
# -----------------------------------------------------------------------------
surfr2_plotly <- function(res, cohort, type = "scatter", highlight = character(), source = "scatter") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("R package 'plotly' is required")
  d <- surfr2_plot_data(res, cohort, type)
  hover <- surfr2_hover_text(res, cohort)
  cols <- surfr2_group_colours(res)

  p <- plotly::plot_ly(source = source) |>
    plotly::add_heatmap(x = d$bg$xmid, y = d$bg$ymid, z = log10(d$bg$z), hoverinfo = "skip",
                        colorscale = list(c(0, "#c4c4c4"), c(1, "#3d3d3d")),
                        showscale = FALSE, name = "all k-mers")
  for (g in names(cols)) {
    i <- which(d$fgroup == g & is.finite(d$fx) & is.finite(d$fy))
    if (!length(i)) next
    p <- p |> plotly::add_trace(
      type = "scattergl", mode = "markers", x = d$fx[i], y = d$fy[i],
      text = hover[i], hoverinfo = "text", customdata = d$fkmer[i],
      marker = list(color = cols[[g]], size = 7, line = list(color = "#333333", width = 0.5)),
      name = sprintf("%s (%d)", g, length(i)))
  }
  hl <- which(d$fkmer %in% highlight)
  if (length(hl)) {
    p <- p |> plotly::add_trace(
      type = "scattergl", mode = "markers", x = d$fx[hl], y = d$fy[hl],
      text = hover[hl], hoverinfo = "text", customdata = d$fkmer[hl],
      marker = list(color = "#d62728", size = 11, symbol = "diamond", line = list(color = "black", width = 1)),
      name = sprintf("search hits (%d)", length(hl)))
  }
  shapes <- if (is.null(d$lines)) list() else lapply(seq_len(nrow(d$lines)), function(i) {
    l <- d$lines[i, ]
    list(type = "line", x0 = l$x0, y0 = l$y0, x1 = l$x1, y1 = l$y1,
         line = list(dash = "dash", width = 1, color = if (l$kind == "reference") "#888888" else "#b22222"))
  })
  plotly::layout(p,
    title = list(text = d$title, font = list(size = 14)),
    xaxis = list(title = d$xlab, range = d$xlim, zeroline = FALSE),
    yaxis = c(list(title = d$ylab, range = d$ylim, zeroline = FALSE),
              if (type == "scatter") list(scaleanchor = "x")),
    shapes = shapes,
    hoverlabel = list(align = "left", font = list(family = "monospace", size = 11)),
    legend = list(orientation = "h", y = -0.15))
}

# backward-compatible name used by the app
surfr2_scatter_plotly <- function(res, cohort, highlight = character(), source = "scatter")
  surfr2_plotly(res, cohort, "scatter", highlight, source)

# -----------------------------------------------------------------------------
# Static plot (ggplot2) for publication
# -----------------------------------------------------------------------------
surfr2_ggplot <- function(res, cohort, type = "scatter") {
  d <- surfr2_plot_data(res, cohort, type)
  bg <- expand.grid(y = d$bg$ymid, x = d$bg$xmid)     # column-major matches z[] layout
  bg$n <- as.vector(d$bg$z)
  bg <- bg[!is.na(bg$n), ]
  g <- ggplot() +
    geom_tile(data = bg, aes(x, y, fill = log10(n)), width = d$bg$xstep, height = d$bg$ystep) +
    scale_fill_gradient(low = "#c4c4c4", high = "#3d3d3d", guide = "none")
  if (!is.null(d$lines))
    g <- g + geom_segment(data = d$lines, aes(x = x0, y = y0, xend = x1, yend = y1, colour = kind),
                          linetype = "dashed", show.legend = FALSE) +
      scale_colour_manual(values = c(threshold = "#b22222", reference = "#888888"))
  if (length(d$fkmer)) {
    pts <- data.frame(x = d$fx, y = d$fy, group = d$fgroup)
    g <- g + ggnewscale_free_points(pts, surfr2_group_colours(res))
  }
  g <- g + theme_classic(base_size = 14) + labs(x = d$xlab, y = d$ylab, title = d$title) +
    theme(plot.title = element_text(size = 11), legend.position = "bottom", legend.title = element_blank())
  if (type == "scatter") g + coord_equal(xlim = d$xlim, ylim = d$ylim)
  else g + coord_cartesian(xlim = d$xlim, ylim = d$ylim)
}

# Passing points with their own fill scale (the background already uses 'fill' for density,
# so the points use shape 21 with a fixed fill per group drawn group by group).
ggnewscale_free_points <- function(pts, cols) {
  present <- names(cols)[names(cols) %in% pts$group]
  labs <- vapply(present, function(g) sprintf("%s (%d)", g, sum(pts$group == g)), character(1))
  layers <- lapply(seq_along(present), function(i) {
    sub <- pts[pts$group == present[i], , drop = FALSE]
    sub$legend <- labs[[i]]
    geom_point(data = sub, aes(x, y, shape = legend), fill = cols[[present[i]]],
               colour = "#333333", size = 1.8, stroke = 0.3)
  })
  # explicit breaks keep legend keys in the same order as the override fills
  c(layers, list(scale_shape_manual(values = setNames(rep(21, length(labs)), unname(labs)),
                                    breaks = unname(labs),
                                    guide = guide_legend(override.aes = list(fill = unname(cols[present]), size = 3)))))
}

# backward-compatible name
surfr2_scatter_ggplot <- function(res, cohort) surfr2_ggplot(res, cohort, "scatter")
