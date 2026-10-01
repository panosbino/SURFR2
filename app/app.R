# =============================================================================
# SURFR2 - interactive explorer (Shiny)
#
# Adds to the static HTML report:
#   * click a case-specific k-mer -> per-sample counts/CPM (plot + downloadable table)
#   * search a known sequence (DNA or RNA, e.g. a mature miRNA) -> highlight its k-mers
#   * switch between replicate cohorts
#
# Needs only the <outdir>/results folder (copy it to a laptop).
#
# Run:
#   SURFR2_RESULTS=/path/to/outdir/results \
#     Rscript -e 'shiny::runApp("SURFR2/app", launch.browser = TRUE)'
# R packages: shiny, plotly, dplyr, readr, ggplot2, jsonlite
# =============================================================================

library(shiny)
source(file.path("..", "bin", "surfr2_plotlib.R"), local = TRUE)

# Which case-specific k-mers does a query sequence hit? Query >= k: k-mers contained in
# the query, plus k-mers whose merged sequence contains it. Query < k: k-mers containing it.
search_hits <- function(res, query) {
  q <- toupper(gsub("[^ACGTUacgtu]", "", query))
  q <- chartr("U", "T", q)
  k <- res$cfg$kmer$k
  km <- res$final$kmer
  if (!nzchar(q) || !length(km)) return(character())
  if (nchar(q) >= k) {
    in_query  <- vapply(km, function(x) grepl(x, q, fixed = TRUE), logical(1))
    contig    <- res$contig_of[km]
    in_contig <- !is.na(contig) & vapply(contig, function(s) !is.na(s) && grepl(q, s, fixed = TRUE), logical(1))
    km[in_query | in_contig]
  } else {
    km[grepl(q, km, fixed = TRUE)]
  }
}

ui <- fluidPage(
  titlePanel("SURFR2 explorer"),
  sidebarLayout(
    sidebarPanel(width = 3,
      textInput("res_dir", "Results directory", value = Sys.getenv("SURFR2_RESULTS", "")),
      actionButton("load", "Load", class = "btn-primary"),
      hr(),
      selectInput("cohort", "Cohort", choices = character()),
      textInput("query", "Highlight sequence (DNA/RNA)", placeholder = "e.g. UGGAGUGUGACAAUGGUGUUUG"),
      htmlOutput("search_info"),
      hr(),
      htmlOutput("summary"),
      helpText("Hover a gold dot for per-condition counts; click it for per-sample counts below.")
    ),
    mainPanel(width = 9,
      plotly::plotlyOutput("scatter", height = "650px"),
      hr(),
      h4(textOutput("sel_title")),
      plotOutput("sel_plot", height = "320px"),
      downloadButton("download", "Download per-sample table (CSV)"),
      tableOutput("sel_table")
    )
  )
)

server <- function(input, output, session) {
  res <- reactiveVal(NULL)
  selected <- reactiveVal(NULL)

  load_dir <- function(d) {
    r <- tryCatch(surfr2_load_results(d), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL); NULL
    })
    if (is.null(r)) return()
    res(r); selected(NULL)
    updateSelectInput(session, "cohort", choices = r$cfg$replicate_cohorts,
                      selected = r$cfg$replicate_cohorts[1])
  }
  observeEvent(input$load, load_dir(input$res_dir))
  observe({                                   # auto-load when started with SURFR2_RESULTS
    isolate(if (is.null(res()) && nzchar(input$res_dir)) load_dir(input$res_dir))
  })

  hits <- reactive({ req(res()); search_hits(res(), input$query) })

  output$summary <- renderUI({
    r <- req(res())
    HTML(sprintf("<b>%s</b><br>%s vs %s<br>%d samples<br>%s candidate k-mers<br>%d case-specific k-mers<br>%d merged sequences",
                 r$cfg$project, r$cfg$comparison$case, r$cfg$comparison$control, nrow(r$norm),
                 format(nrow(r$stats), big.mark = ","), nrow(r$final),
                 if (is.null(r$contigs)) 0L else nrow(r$contigs)))
  })
  output$search_info <- renderUI({
    if (!nzchar(input$query)) return(NULL)
    n <- length(hits())
    HTML(if (n) sprintf("<span style='color:#d62728'>%d case-specific k-mer(s) match</span>", n)
         else "No case-specific k-mer matches (the sequence may not have passed the filters).")
  })

  output$scatter <- plotly::renderPlotly({
    r <- req(res()); req(input$cohort %in% r$cfg$replicate_cohorts)
    surfr2_scatter_plotly(r, input$cohort, highlight = hits(), source = "scatter") |>
      plotly::event_register("plotly_click")
  })

  observeEvent(plotly::event_data("plotly_click", source = "scatter"), {
    ev <- plotly::event_data("plotly_click", source = "scatter")
    if (!is.null(ev$customdata)) selected(as.character(ev$customdata[[1]]))
  })

  sel_table <- reactive({
    r <- req(res()); km <- req(selected())
    req(km %in% rownames(r$cnt))
    data.frame(sample_id = r$norm$sample_id, cohort = r$norm$cohort, condition = r$norm$condition,
               count = as.numeric(r$cnt[km, ]), cpm = round(as.numeric(r$cpm[km, ]), 3),
               eff_libsize = r$norm$eff_libsize, check.names = FALSE)
  })

  output$sel_title <- renderText({
    km <- selected()
    if (is.null(km)) return("Click a case-specific k-mer to see per-sample counts")
    contig <- res()$contig_of[km]
    sprintf("%s  (sequence: %s)", km, ifelse(is.na(contig), "not merged", contig))
  })
  output$sel_plot <- renderPlot({
    d <- sel_table()
    d$sample_id <- factor(d$sample_id, levels = d$sample_id[order(d$condition, d$cohort, -d$cpm)])
    ggplot(d, aes(sample_id, cpm, fill = condition)) +
      geom_col() + facet_grid(~ cohort, scales = "free_x", space = "free_x") +
      theme_bw(base_size = 13) + labs(x = NULL, y = "CPM") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = if (nrow(d) > 60) 5 else 9))
  })
  output$sel_table <- renderTable(sel_table(), digits = 3)
  output$download <- downloadHandler(
    filename = function() sprintf("SURFR2_%s_%s.csv", res()$cfg$project, selected()),
    content  = function(f) write.csv(sel_table(), f, row.names = FALSE)
  )
}

shinyApp(ui, server)
