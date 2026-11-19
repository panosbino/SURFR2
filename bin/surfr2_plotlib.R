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

# Genome names in config order (empty without genome mapping)
surfr2_genomes <- function(res) {
  g <- res$cfg$genomes
  if (is.null(g) || !length(g) || is.null(res$contigs) || !"genome_category" %in% names(res$contigs))
    return(character())
  if (is.data.frame(g)) g$name else vapply(g, function(x) x$name, character(1))
}

# Genome category of each passing k-mer, from the sequence it was merged into
surfr2_kmer_genome <- function(res) {
  cat_of <- setNames(res$contigs$genome_category, res$contigs$contig)
  out <- unname(cat_of[res$contig_of[res$final$kmer]])
  out[is.na(out)] <- "unmapped"
  out
}

surfr2_category_label <- function(cat, genomes) {
  ifelse(cat == "both", "both genomes", ifelse(cat == "multiple", "several genomes", cat))
}

# Group of each passing k-mer for colouring: the genome category when sequences were
# mapped to genomes; otherwise the higher condition (DEA) or the analysis label
surfr2_pass_group <- function(res) {
  gn <- surfr2_genomes(res)
  if (length(gn)) return(surfr2_category_label(surfr2_kmer_genome(res), gn))
  if (res$analysis == "dea") paste("higher in", res$final$higher_in) else rep(res$label, nrow(res$final))
}

surfr2_html_escape <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)

# Pseudocount for plot axes: log10(CPM + this). Separate from filters.pseudocount_cpm,
# which enters fold changes and thresholds. Older results without 'plots' use 1.
surfr2_plot_pc <- function(cfg) {
  pp <- cfg$plots$pseudocount_cpm
  if (is.null(pp)) 1 else pp
}

# Hover block with the best hit per genome for each of the given sequences ("" if no mapping)
surfr2_genome_block <- function(res, contig) {
  gn <- surfr2_genomes(res)
  if (!length(gn) || !length(contig)) return(rep("", length(contig)))
  ct <- res$contigs[match(contig, res$contigs$contig), , drop = FALSE]
  per <- vapply(gn, function(g) {
    mm <- as.character(ct[[paste0(g, "_mismatches")]]); mm[is.na(mm)] <- ""
    hits <- as.character(ct[[paste0(g, "_hits")]]); locus <- as.character(ct[[paste0(g, "_locus")]])
    line <- ifelse(mm == "", sprintf("%s: not mapped", g),
                   sprintf("%s: %s, %s hit(s), %s", g, ifelse(mm == "0", "exact", paste(mm, "mismatch(es)")), hits, locus))
    ori_col <- paste0(g, "_gene_orientation")
    if (ori_col %in% names(ct)) {                 # GTF annotation configured for this genome
      ori <- as.character(ct[[ori_col]]); ori[is.na(ori)] <- ""
      gname <- as.character(ct[[paste0(g, "_gene_name")]]); gname[is.na(gname)] <- ""
      btype <- as.character(ct[[paste0(g, "_gene_biotype")]]); btype[is.na(btype)] <- ""
      nh <- as.character(ct[[paste0(g, "_annotated_hits")]]); nh[is.na(nh)] <- ""
      ann <- ifelse(ori == "", "",
             ifelse(ori == "intergenic", "<br>&nbsp;&nbsp;genes: none (intergenic)",
                    sprintf("<br>&nbsp;&nbsp;genes (%s, %s hits): %s%s", ori, nh, surfr2_html_escape(gname),
                            ifelse(btype == "", "", paste0(" [", surfr2_html_escape(btype), "]")))))
      line <- paste0(line, ann)
    }
    line
  }, character(length(contig)))
  if (is.null(dim(per))) per <- matrix(per, nrow = 1)
  paste0("<br><br><b>genome mapping</b><br>", apply(per, 1, paste, collapse = "<br>"))
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

  genome_block <- surfr2_genome_block(res, contig)

  paste0("<b>", f$kmer, "</b>",
         "<br>sequence: ", ifelse(is.na(contig), "(not merged)", contig),
         if (res$analysis == "dea") paste0("<br>higher in: ", f$higher_in) else "",
         sprintf("<br>%s (%s / %s, %s): %.2f", fc_lab, res$CASE, res$CTRL, cohort, f[[fc_col]]),
         if (paste0(cohort, "_fdr") %in% names(f))
           sprintf("<br>FDR (%s, edgeR QL): %.3g", cohort, f[[paste0(cohort, "_fdr")]])
         else "<br>no statistical test (descriptive mode)",
         genome_block,
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

# Sequence-level data: one point per merged sequence. Position = median over its k-mers
# of their mean CPM per condition (robust to one odd k-mer); colour = genome category when
# mapped, else the higher condition (DEA) or the analysis label.
surfr2_sequence_points <- function(res, cohort) {
  ct <- res$contigs
  if (is.null(ct) || !nrow(ct)) return(NULL)
  pc <- surfr2_plot_pc(res$cfg)
  mc <- paste0(cohort, "_mean_cpm_", res$CASE); mn <- paste0(cohort, "_mean_cpm_", res$CTRL)
  st <- res$stats[match(names(res$contig_of), res$stats$kmer), , drop = FALSE]
  # members = every called k-mer the sequence contains. Not res$contig_of, which assigns a
  # k-mer to one sequence only: mergeTags can return a sequence nested inside another, and
  # the nested one would be left without k-mers (NA position, silently not drawn).
  k <- res$cfg$kmer$k; kn <- names(res$contig_of)
  mem <- lapply(ct$contig, function(s) {
    n <- nchar(s) - k + 1
    if (n < 1) return(integer())
    m <- match(substring(s, seq_len(n), seq_len(n) + k - 1), kn)
    unique(m[!is.na(m)])
  })
  med <- function(col) vapply(mem, function(i) if (length(i)) stats::median(st[[col]][i]) else NA_real_, numeric(1))
  m_case <- med(mc); m_ctrl <- med(mn)
  gn <- surfr2_genomes(res)
  group <- if (length(gn)) surfr2_category_label(ct$genome_category, gn)
           else if (res$analysis == "dea") paste("higher in", ct$higher_in) else rep(res$label, nrow(ct))
  fdr_col <- paste0(cohort, "_min_kmer_fdr"); fc_col <- grep("^log2fc_", names(ct), value = TRUE)[1]
  p_col <- paste0(cohort, "_acat_pvalue")
  num <- function(col) if (!is.na(col) && col %in% names(ct)) as.numeric(ct[[col]]) else rep(NA_real_, nrow(ct))
  fc <- num(fc_col); fdr <- num(fdr_col); pval <- num(p_col)
  hover <- paste0("<b>", ct$contig, "</b>",
                  sprintf("<br>%d nt, %d k-mers merged", nchar(ct$contig), as.integer(ct$nb_merged_kmers)),
                  if (res$analysis == "dea") paste0("<br>higher in: ", ct$higher_in) else "",
                  if (!is.na(fc_col)) sprintf("<br>log2FC (%s / %s, representative k-mer): %.2f", res$CASE, res$CTRL, ct[[fc_col]]) else "",
                  if (fdr_col %in% names(ct)) sprintf("<br>smallest k-mer FDR (%s): %.3g", cohort, ct[[fdr_col]]) else "",
                  if (p_col %in% names(ct)) sprintf("<br>sequence p-value (ACAT of its k-mers): %.3g", pval) else "",
                  sprintf("<br>median mean CPM: %s %.2f, %s %.2f", res$CASE, m_case, res$CTRL, m_ctrl),
                  surfr2_genome_block(res, ct$contig))
  list(x = log10(m_ctrl + pc), y = log10(m_case + pc), key = ct$tag, group = group, hover = hover,
       copy = ct$contig, fc = fc, fdr = fdr, pval = pval, cpm_case = unname(m_case), cpm_ctrl = unname(m_ctrl),
       len = nchar(ct$contig), ann = surfr2_annotation_text(res, ct))
}

# One-line gene annotation per sequence over the annotated genomes ("" without a GTF),
# e.g. "mouse: Mir124a-1;Mir124a-2 [miRNA] (sense)". Used in the filtered-table download.
surfr2_annotation_text <- function(res, ct) {
  gn <- surfr2_genomes(res)
  parts <- lapply(gn, function(g) {
    oc <- paste0(g, "_gene_orientation")
    if (!oc %in% names(ct)) return(NULL)
    o <- as.character(ct[[oc]]); o[is.na(o)] <- ""
    nm <- as.character(ct[[paste0(g, "_gene_name")]]); nm[is.na(nm)] <- ""
    bt <- as.character(ct[[paste0(g, "_gene_biotype")]]); bt[is.na(bt)] <- ""
    ifelse(o == "", "", ifelse(o == "intergenic", paste0(g, ": intergenic"),
                               sprintf("%s: %s%s (%s)", g, nm, ifelse(bt == "", "", paste0(" [", bt, "]")), o)))
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(rep("", nrow(ct)))
  apply(do.call(cbind, parts), 1, function(x) paste(x[nzchar(x)], collapse = "; "))
}

# Text copied from a point of an interactive k-mer plot: the merged sequence the k-mer
# belongs to, or the k-mer itself when it was not merged.
surfr2_copy_text <- function(res, kmers) {
  ct <- if (is.null(res$contig_of)) rep(NA_character_, length(kmers)) else unname(res$contig_of[kmers])
  ifelse(is.na(ct), kmers, ct)
}

# Types: scatter / volcano (one point per called k-mer), sequences / sequences_volcano
# (one point per merged sequence; the volcano uses the sequence's ACAT p-value).
surfr2_plot_data <- function(res, cohort, type = c("scatter", "volcano", "sequences", "sequences_volcano"),
                             nbins = 150) {
  type <- match.arg(type)
  cfg <- res$cfg; s <- res$stats
  pc <- surfr2_plot_pc(cfg)                  # axes: log10(CPM + pc)
  pcf <- cfg$filters$pseudocount_cpm         # fold changes / thresholds were computed with this
  fin <- match(res$final$kmer, s$kmer)
  seg <- function(x0, y0, x1, y1, kind = "threshold") data.frame(x0, y0, x1, y1, kind)

  if (type %in% c("scatter", "sequences")) {
    mc <- paste0(cohort, "_mean_cpm_", res$CASE); mn <- paste0(cohort, "_mean_cpm_", res$CTRL)
    if (!all(c(mc, mn) %in% names(s))) stop("cohort not in results: ", cohort)
    x <- log10(s[[mn]] + pc); y <- log10(s[[mc]] + pc)
    lim <- range(c(x, y), finite = TRUE); xlim <- ylim <- lim
    if (res$analysis == "specific" && type == "scatter") {
      y_min <- log10(cfg$specific$min_mean_cpm_case + pc)
      # min_log2fc: log2((case + pcf) / (ctrl + pcf)) >= fc, drawn in log10(x + pc) axes.
      # A straight diagonal when pc == pcf; otherwise a curve, drawn as short segments.
      xs <- seq(lim[1], lim[2], length.out = 200)
      ctrl <- pmax(10^xs - pc, 0)
      ys <- log10(2^cfg$specific$min_log2fc * (ctrl + pcf) - pcf + pc)
      ok <- is.finite(ys) & ys <= lim[2]
      xs <- xs[ok]; ys <- ys[ok]
      lines <- rbind(seg(lim[1], y_min, lim[2], y_min),                    # min_mean_cpm_case
                     if (length(xs) > 1) seg(head(xs, -1), head(ys, -1), xs[-1], ys[-1]))   # min_log2fc
    } else {
      lines <- seg(lim[1], lim[1], lim[2], lim[2], "reference")            # equal expression
    }
    xlab <- sprintf("log10(mean CPM %s + %g)", res$CTRL, pc)
    ylab <- sprintf("log10(mean CPM %s + %g)", res$CASE, pc)
    keep <- rep(TRUE, length(x))
  } else {          # volcano and sequences_volcano: background = all tested k-mers
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
  out <- list(type = type, cohort = cohort, xlim = xlim, ylim = ylim, bg = bg, lines = lines,
              n_all = sum(keep), fx = x[fin], fy = y[fin], fkmer = res$final$kmer,
              fcopy = surfr2_copy_text(res, res$final$kmer),
              fgroup = surfr2_pass_group(res), hover = NULL, xlab = xlab, ylab = ylab,
              title = sprintf("%s: %s vs %s (%s %s k-mers, %d %s)", cohort, res$CASE, res$CTRL,
                              format(sum(keep), big.mark = ","), if (type == "volcano") "tested" else "candidate",
                              length(fin), res$label))
  if (type %in% c("sequences", "sequences_volcano")) {
    sp <- surfr2_sequence_points(res, cohort)
    if (is.null(sp)) stop("no merged sequences to plot")
    out$fkmer <- sp$key; out$fgroup <- sp$group; out$hover <- sp$hover; out$fcopy <- sp$copy
    out$seq <- sp                            # per-sequence values for the interactive filter panel
    if (type == "sequences") {
      out$fx <- sp$x; out$fy <- sp$y
      out$title <- sprintf("%s: %d %s sequences (%s vs %s)", cohort, length(sp$x), res$label, res$CASE, res$CTRL)
      out$xlab <- sprintf("log10(median k-mer mean CPM %s + %g)", res$CTRL, pc)
      out$ylab <- sprintf("log10(median k-mer mean CPM %s + %g)", res$CASE, pc)
    } else {
      if (all(is.na(sp$pval))) stop("no sequence p-values (ACAT) in these results: re-run the filter step")
      out$fx <- sp$fc; out$fy <- -log10(pmax(sp$pval, .Machine$double.xmin))
      ok <- is.finite(out$fx) & is.finite(out$fy)
      out$xlim <- c(-1, 1) * max(abs(c(out$xlim, out$fx[ok])))
      out$ylim <- c(0, max(c(out$ylim[2], out$fy[ok] * 1.03)))
      # the k-mer FDR line does not apply to sequence p-values; keep the fold-change cut-offs
      out$lines <- NULL
      if (res$analysis == "dea") {
        t <- cfg$dea$min_abs_log2fc
        out$lines <- rbind(seg(-t, out$ylim[1], -t, out$ylim[2]), seg(t, out$ylim[1], t, out$ylim[2]))
      }
      out$title <- sprintf("%s: %d %s sequences (%s vs %s); background: tested k-mers", cohort, length(sp$x),
                           res$label, res$CASE, res$CTRL)
      out$xlab <- sprintf("log2 fold change (%s / %s, representative k-mer)", res$CASE, res$CTRL)
      out$ylab <- "-log10 sequence p-value (ACAT of its k-mers)"
    }
  }
  out
}

surfr2_group_colours <- function(res) {
  gn <- surfr2_genomes(res)
  if (length(gn)) {
    # Okabe-Ito (colour-blind safe) for genomes; purple = several, black = unmapped
    pal <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#E69F00", "#F0E442")
    cols <- setNames(rep_len(pal, length(gn)), gn)
    multi <- if (length(gn) == 2) "both genomes" else "several genomes"
    return(c(cols, setNames("#7B3294", multi), unmapped = "#000000"))
  }
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
  hover <- if (!is.null(d$hover)) d$hover else surfr2_hover_text(res, cohort)
  is_seq <- type %in% c("sequences", "sequences_volcano")
  hover <- paste0(hover, "<br><i>click or press C to copy the ",
                  if (is_seq) "sequence" else "sequence (k-mer if not merged)", "</i>")
  cols <- surfr2_group_colours(res)
  # per-point values travel in the trace's meta: copy text, and for sequence plots the
  # values the filter panel filters on (as.list keeps one-point traces as arrays)
  copy_meta <- function(i) {
    m <- list(surfr2_copy = as.list(d$fcopy[i]))
    if (is_seq) {
      s <- d$seq
      na <- function(v) as.list(v[i])          # NA -> JSON null (htmlwidgets na = "null")
      m$surfr2_filt <- list(fdr = na(s$fdr), fc = na(s$fc), cpm_case = na(s$cpm_case),
                            cpm_ctrl = na(s$cpm_ctrl), len = as.list(s$len[i]), ann = as.list(s$ann[i]),
                            conds = list(res$CASE, res$CTRL))
    }
    m
  }

  p <- plotly::plot_ly(source = source) |>
    plotly::add_heatmap(x = d$bg$xmid, y = d$bg$ymid, z = log10(d$bg$z), hoverinfo = "skip",
                        colorscale = list(c(0, "#c4c4c4"), c(1, "#3d3d3d")),
                        showscale = FALSE, name = "all k-mers")
  for (g in names(cols)) {
    i <- which(d$fgroup == g & is.finite(d$fx) & is.finite(d$fy))
    if (!length(i)) next
    p <- p |> plotly::add_trace(
      type = "scattergl", mode = "markers", x = d$fx[i], y = d$fy[i],
      text = hover[i], hoverinfo = "text", customdata = d$fkmer[i], meta = copy_meta(i),
      marker = list(color = cols[[g]], size = if (is_seq) 10 else 7,
                    line = list(color = "#333333", width = 0.5)),
      name = sprintf("%s (%d)", g, length(i)))
  }
  hl <- which(d$fkmer %in% highlight)
  if (length(hl)) {
    p <- p |> plotly::add_trace(
      type = "scattergl", mode = "markers", x = d$fx[hl], y = d$fy[hl],
      text = hover[hl], hoverinfo = "text", customdata = d$fkmer[hl], meta = copy_meta(hl),
      marker = list(color = "#d62728", size = 11, symbol = "diamond", line = list(color = "black", width = 1)),
      name = sprintf("search hits (%d)", length(hl)))
  }
  shapes <- if (is.null(d$lines)) list() else lapply(seq_len(nrow(d$lines)), function(i) {
    l <- d$lines[i, ]
    list(type = "line", x0 = l$x0, y0 = l$y0, x1 = l$x1, y1 = l$y1,
         line = list(dash = "dash", width = 1, color = if (l$kind == "reference") "#888888" else "#b22222"))
  })
  p <- plotly::layout(p,
    title = list(text = d$title, font = list(size = 14)),
    xaxis = list(title = d$xlab, range = d$xlim, zeroline = FALSE),
    yaxis = c(list(title = d$ylab, range = d$ylim, zeroline = FALSE),
              if (!type %in% c("volcano", "sequences_volcano")) list(scaleanchor = "x")),
    shapes = shapes,
    hoverlabel = list(align = "left", font = list(family = "monospace", size = 11)),
    legend = list(orientation = "h", y = -0.15))
  p <- htmlwidgets::onRender(p, SURFR2_COPY_JS)
  htmlwidgets::onRender(p, SURFR2_FILTER_JS)
}

# Filter panel for the sequence plots (sequences, sequences_volcano): number boxes above
# the plot narrow the shown sequences by FDR (smallest k-mer FDR), |log2FC|, CPM in the
# higher condition (median k-mer mean CPM, the plotted value) and length; empty = no
# limit. Shows "N of M" and downloads the shown sequences as TSV. The values come from
# each trace's meta.surfr2_filt. Filters only narrow the called set: a looser cut-off than
# the run's needs a re-run. Removes itself on plots without filter data (Shiny re-renders).
SURFR2_FILTER_JS <- "
function(el) {
  var gd = el;
  var has = (gd.data || []).some(function(t) { return t.meta && t.meta.surfr2_filt; });
  if (!has) { if (el._surfr2Panel) { el._surfr2Panel.remove(); el._surfr2Panel = null; } return; }
  var arr = function(v, n) { return Array.isArray(v) ? v.slice() : (v == null ? [] : (n > 1 ? Array(n).fill(v) : [v])); };
  var orig = [];
  gd.data.forEach(function(t, i) {
    if (!(t.meta && t.meta.surfr2_filt)) return;
    var x = arr(t.x), n = x.length, f = t.meta.surfr2_filt, ff = {};
    Object.keys(f).forEach(function(k) { ff[k] = k === 'conds' ? f[k] : arr(f[k], n); });
    orig.push({i: i, x: x, y: arr(t.y, n), text: arr(t.text, n), cd: arr(t.customdata, n),
               copy: arr(t.meta.surfr2_copy, n), f: ff, name: String(t.name).replace(/ [(][0-9]+[)]$/, '')});
  });
  var conds = orig[0].f.conds || ['case', 'control'];
  var total = orig.reduce(function(s, o) { return s + o.x.length; }, 0);
  var all = function(k) { return [].concat.apply([], orig.map(function(o) { return o.f[k]; })).filter(function(v) { return v != null; }); };
  var hiCpm = function(o, j) { return Math.max(o.f.cpm_case[j] || 0, o.f.cpm_ctrl[j] || 0); };
  var cpms = [].concat.apply([], orig.map(function(o) { return o.x.map(function(_, j) { return hiCpm(o, j); }); }));
  var lens = all('len'), hasFdr = all('fdr').length > 0;
  var fmt = function(v) { return v >= 100 ? Math.round(v) : +v.toPrecision(3); };
  var panel = el._surfr2Panel;
  if (!panel) {
    panel = document.createElement('div'); el._surfr2Panel = panel;
    panel.style.cssText = 'font:12px sans-serif;display:flex;flex-wrap:wrap;gap:6px 14px;align-items:center;' +
      'padding:6px 8px;margin:0 0 4px 0;border:1px solid #ccc;border-radius:4px;background:#fafafa;color:#222';
    var box = function(id, label, title) {
      return '<label title=\"' + title + '\">' + label + ' <input data-k=\"' + id + '\" type=\"number\" step=\"any\" min=\"0\" ' +
             'style=\"width:72px;font:12px sans-serif\"></label>';
    };
    panel.innerHTML = '<b>Filter sequences</b>' +
      '<span data-g=\"fdr\">' + box('fdr', 'FDR &le;', 'smallest k-mer FDR of the sequence') + '</span>' +
      box('fc', '|log2FC| &ge;', 'absolute log2 fold change (representative k-mer)') +
      '<span>CPM ' + box('cpm_lo', 'from', 'median k-mer mean CPM in the higher condition') + box('cpm_hi', 'to', '') + '</span>' +
      '<span>length ' + box('len_lo', 'from', 'sequence length (nt)') + box('len_hi', 'to', '') + '</span>' +
      '<button data-a=\"reset\" style=\"font:12px sans-serif\">reset</button>' +
      '<button data-a=\"tsv\" style=\"font:12px sans-serif\">download shown (TSV)</button>' +
      '<span data-o=\"count\" style=\"color:#555\"></span>';
    el.parentNode.insertBefore(panel, el);
    var timer = null;
    panel.addEventListener('input', function() { clearTimeout(timer); timer = setTimeout(function() { el._surfr2Apply(); }, 150); });
    panel.querySelector('[data-a=reset]').addEventListener('click', function() {
      panel.querySelectorAll('input').forEach(function(i) { i.value = ''; }); el._surfr2Apply();
    });
    panel.querySelector('[data-a=tsv]').addEventListener('click', function() { el._surfr2Tsv(); });
  }
  panel.querySelector('[data-g=fdr]').style.display = hasFdr ? '' : 'none';
  var ph = {fdr: 'any', fc: 'any', cpm_lo: fmt(Math.min.apply(null, cpms)), cpm_hi: fmt(Math.max.apply(null, cpms)),
            len_lo: Math.min.apply(null, lens), len_hi: Math.max.apply(null, lens)};
  Object.keys(ph).forEach(function(k) { panel.querySelector('[data-k=' + k + ']').placeholder = ph[k]; });
  var val = function(k) { var v = parseFloat(panel.querySelector('[data-k=' + k + ']').value); return isNaN(v) ? null : v; };
  var keepIdx = function(o) {
    var fdr = val('fdr'), fc = val('fc'), clo = val('cpm_lo'), chi = val('cpm_hi'), llo = val('len_lo'), lhi = val('len_hi');
    var out = [];
    for (var j = 0; j < o.x.length; j++) {
      var F = o.f, c = hiCpm(o, j);
      if (fdr !== null && hasFdr && !(F.fdr[j] != null && F.fdr[j] <= fdr)) continue;
      if (fc !== null && !(F.fc[j] != null && Math.abs(F.fc[j]) >= fc)) continue;
      if (clo !== null && c < clo) continue;
      if (chi !== null && c > chi) continue;
      if (llo !== null && F.len[j] < llo) continue;
      if (lhi !== null && F.len[j] > lhi) continue;
      out.push(j);
    }
    return out;
  };
  el._surfr2Apply = function() {
    var upd = {x: [], y: [], text: [], customdata: [], name: [], meta: []}, idx = [], shown = 0;
    orig.forEach(function(o) {
      var k = keepIdx(o), pick = function(a) { return k.map(function(j) { return a[j]; }); }, f = {};
      Object.keys(o.f).forEach(function(key) { f[key] = key === 'conds' ? o.f[key] : pick(o.f[key]); });
      upd.x.push(pick(o.x)); upd.y.push(pick(o.y)); upd.text.push(pick(o.text)); upd.customdata.push(pick(o.cd));
      upd.name.push(o.name + ' (' + k.length + ')'); upd.meta.push({surfr2_copy: pick(o.copy), surfr2_filt: f});
      idx.push(o.i); shown += k.length;
    });
    Plotly.restyle(gd, upd, idx);
    panel.querySelector('[data-o=count]').textContent = 'showing ' + shown + ' of ' + total + ' sequences';
  };
  el._surfr2Tsv = function() {
    var rows = [['sequence', 'group', 'log2fc', 'min_kmer_fdr', 'median_cpm_' + conds[0], 'median_cpm_' + conds[1],
                 'length', 'annotation'].join('\\t')];
    orig.forEach(function(o) {
      keepIdx(o).forEach(function(j) {
        var F = o.f, s = function(v) { return v == null ? '' : String(v); };
        rows.push([o.copy[j], o.name, s(F.fc[j]), s(F.fdr[j]), s(F.cpm_case[j]), s(F.cpm_ctrl[j]), s(F.len[j]),
                   s(F.ann && F.ann[j])].join('\\t'));
      });
    });
    var blob = new Blob([rows.join('\\n') + '\\n'], {type: 'text/tab-separated-values'});
    var a = document.createElement('a'); a.href = URL.createObjectURL(blob);
    a.download = 'surfr2_sequences_filtered.tsv'; document.body.appendChild(a); a.click();
    setTimeout(function() { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  };
  el._surfr2Apply();
}"

# Copy a point's sequence to the clipboard: click it, or press C while hovering it
# (Ctrl/Cmd+C also works when no page text is selected). The text comes from the
# trace's meta.surfr2_copy, indexed like the points. Works in the saved HTML (file://)
# and in the Shiny app; falls back to execCommand where the Clipboard API is unavailable.
SURFR2_COPY_JS <- "
function(el) {
  if (el._surfr2Copy) return;            // Shiny re-renders call this again: register once
  el._surfr2Copy = true;
  var gd = el, current = null;
  function textOf(pt) {
    var tr = gd.data[pt.curveNumber], m = tr && tr.meta && tr.meta.surfr2_copy;
    if (m == null) return null;
    var v = Array.isArray(m) ? m[pt.pointNumber] : m;
    return v == null ? null : String(v);
  }
  function toast(msg) {
    var t = el.querySelector('.surfr2-toast');
    if (!t) {
      t = document.createElement('div'); t.className = 'surfr2-toast';
      t.style.cssText = 'position:absolute;left:50%;top:8px;transform:translateX(-50%);z-index:1000;' +
        'background:rgba(30,30,30,.9);color:#fff;font:12px monospace;padding:6px 10px;border-radius:4px;' +
        'max-width:90%;overflow-wrap:anywhere;pointer-events:none;transition:opacity .3s';
      if (getComputedStyle(el).position === 'static') el.style.position = 'relative';
      el.appendChild(t);
    }
    t.textContent = msg; t.style.opacity = 1;
    clearTimeout(t._h); t._h = setTimeout(function() { t.style.opacity = 0; }, 1800);
  }
  function fallback(s) {
    var ta = document.createElement('textarea');
    ta.value = s; ta.style.cssText = 'position:fixed;top:0;left:0;opacity:0';
    document.body.appendChild(ta); ta.select();
    var ok = false; try { ok = document.execCommand('copy'); } catch (e) {}
    document.body.removeChild(ta); return ok;
  }
  function copy(s) {
    if (!s) return;
    var done = function() { toast('Copied ' + s); };
    var fail = function() { if (fallback(s)) done(); else toast('Copy failed - sequence: ' + s); };
    if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(s).then(done, fail);
    else fail();
  }
  gd.on('plotly_hover', function(e) { current = e && e.points && e.points[0] ? textOf(e.points[0]) : null; });
  gd.on('plotly_unhover', function() { current = null; });
  gd.on('plotly_click', function(e) { if (e && e.points && e.points[0]) copy(textOf(e.points[0])); });
  document.addEventListener('keydown', function(e) {
    if (!current || e.altKey) return;
    var a = document.activeElement, tag = a && a.tagName;
    if (tag === 'INPUT' || tag === 'TEXTAREA' || (a && a.isContentEditable)) return;
    var k = (e.key || '').toLowerCase();
    var plain = k === 'c' && !e.ctrlKey && !e.metaKey;
    var sysCopy = k === 'c' && (e.ctrlKey || e.metaKey) && !String(window.getSelection()).length;
    if (plain || sysCopy) { e.preventDefault(); copy(current); }
  });
}"

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
    g <- g + ggnewscale_free_points(pts, surfr2_group_colours(res),
                                    size = if (type %in% c("sequences", "sequences_volcano")) 3 else 1.8)
  }
  g <- g + theme_classic(base_size = 14) + labs(x = d$xlab, y = d$ylab, title = d$title) +
    theme(plot.title = element_text(size = 11), legend.position = "bottom", legend.title = element_blank())
  if (!type %in% c("volcano", "sequences_volcano")) g + coord_equal(xlim = d$xlim, ylim = d$ylim)
  else g + coord_cartesian(xlim = d$xlim, ylim = d$ylim)
}

# Passing points with their own fill scale (the background already uses 'fill' for density,
# so the points use shape 21 with a fixed fill per group drawn group by group).
ggnewscale_free_points <- function(pts, cols, size = 1.8) {
  present <- names(cols)[names(cols) %in% pts$group]
  labs <- vapply(present, function(g) sprintf("%s (%d)", g, sum(pts$group == g)), character(1))
  layers <- lapply(seq_along(present), function(i) {
    sub <- pts[pts$group == present[i], , drop = FALSE]
    sub$legend <- labs[[i]]
    geom_point(data = sub, aes(x, y, shape = legend), fill = cols[[present[i]]],
               colour = "#333333", size = size, stroke = 0.3)
  })
  # explicit breaks keep legend keys in the same order as the override fills
  c(layers, list(scale_shape_manual(values = setNames(rep(21, length(labs)), unname(labs)),
                                    breaks = unname(labs),
                                    guide = guide_legend(override.aes = list(fill = unname(cols[present]), size = 3)))))
}

# backward-compatible name
surfr2_scatter_ggplot <- function(res, cohort) surfr2_ggplot(res, cohort, "scatter")
