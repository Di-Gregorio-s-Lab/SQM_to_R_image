source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

project_dir <- Sys.getenv("SQM_CS8_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) {
  stop("SQM_CS8_PROJECT_DIR is required.", call. = FALSE)
}
if (!dir.exists(project_dir)) {
  stop("SQM_CS8_PROJECT_DIR does not exist: ", project_dir, call. = FALSE)
}

script_env <- t1_source_sqm_plots_without_main()
sqm <- script_env$load_sqm_project(project_dir, "prokfilter")
selected_samples <- c("CS8T0", "CS8T2", "CS8T3", "CS8T4", "CS8T6")
expected_k01563 <- c(
  CS8T0 = 20.461747379587,
  CS8T2 = 64.241150618838,
  CS8T3 = 48.990246839864,
  CS8T4 = 111.090699101642,
  CS8T6 = 38.836884026958
)

fixture <- t1_read_k01563_fixture()
fixture_nodes <- fixture[
  as.character(fixture$pathway_id) == "00361" &
    as.character(fixture$type) == "ortholog",
  ,
  drop = FALSE
]
fixture_kos <- unique(unlist(lapply(
  as.character(fixture_nodes$kegg_names),
  function(value) {
    matches <- gregexpr("K[0-9]{5}", value, perl = TRUE)
    regmatches(value, matches)[[1L]]
  }
)))
t1_expect_identical(
  fixture_kos,
  "K01563",
  "Fixture 00361 no longer identifies only K01563"
)

analysis <- script_env$build_pathway_analysis(
  list(
    pathway_name = "Chlorocyclohexane and chlorobenzene degradation",
    pathway_id = "00361",
    pathway_selection = "defined",
    pathway_sqm = sqm,
    context_sqm = sqm,
    pathway_ko_ids = fixture_kos
  ),
  selected_samples
)
flow_rank <- script_env$build_flow_table_for_rank(
  orf_long = analysis$orf_long_result$data,
  rank = "family",
  selected_samples = selected_samples,
  top_n_taxa = 10L,
  top_n_ko = 20L,
  ko_lookup = analysis$ko_lookup
)

render_dir <- tempfile("t5_cs8_simple_flow_")
dir.create(render_dir, recursive = TRUE)
on.exit(unlink(render_dir, recursive = TRUE, force = TRUE), add = TRUE)

for (sample_name in selected_samples) {
  flow_tbl <- script_env$build_flow_table_for_sample(
    flow_rank,
    analysis$pathway_name,
    "family",
    sample_name
  )
  flow_before_render <- flow_tbl
  t1_expect_equal(
    sum(flow_tbl$TPM[as.character(flow_tbl$KO) == "K01563"]),
    unname(expected_k01563[[sample_name]]),
    paste0("CS8 FLOW K01563 TPM changed for ", sample_name),
    tolerance = 1e-8
  )

  plot_object <- script_env$make_flow_plot(
    flow_tbl,
    analysis$pathway_name,
    "family",
    sample_name
  )
  built_plot <- ggplot2::ggplot_build(plot_object)
  flow_intervals <- built_plot$data[[1L]] |>
    dplyr::filter(.data$x == 1) |>
    dplyr::group_by(.data$stratum) |>
    dplyr::summarise(
      flow_ymin = min(.data$ymin),
      flow_ymax = max(.data$ymax),
      .groups = "drop"
    ) |>
    dplyr::mutate(stratum = as.character(.data$stratum))
  stratum_intervals <- built_plot$data[[2L]] |>
    dplyr::filter(.data$x == 1) |>
    dplyr::transmute(
      stratum = as.character(.data$stratum),
      stratum_ymin = .data$ymin,
      stratum_ymax = .data$ymax
    )
  geometry <- dplyr::left_join(
    flow_intervals,
    stratum_intervals,
    by = "stratum",
    relationship = "one-to-one"
  ) |>
    dplyr::mutate(
      delta = pmax(
        abs(.data$flow_ymin - .data$stratum_ymin),
        abs(.data$flow_ymax - .data$stratum_ymax)
      )
    )
  if (anyNA(geometry$delta) || any(geometry$delta > 1e-8)) {
    stop("CS8 FLOW ribbons and taxonomy strata are displaced for ", sample_name, ".", call. = FALSE)
  }
  if (!identical(plot_object$theme$legend.position, "none")) {
    stop("CS8 simple FLOW PNG rendered a legend for ", sample_name, ".", call. = FALSE)
  }
  if (!identical(unique(toupper(as.character(built_plot$data[[2L]]$fill))), "GREY95")) {
    stop("CS8 simple FLOW columns are not neutral for ", sample_name, ".", call. = FALSE)
  }

  sankey <- script_env$make_flow_sankey(
    flow_tbl,
    analysis$pathway_name,
    "family",
    sample_name
  )
  trace <- plotly::plotly_build(sankey)$x$data[[1L]]
  if (!identical(as.character(trace$arrangement), "snap") ||
      length(trace$node$x) > 0L || length(trace$node$y) > 0L) {
    stop("CS8 FLOW HTML still forces node positions for ", sample_name, ".", call. = FALSE)
  }

  png_path <- file.path(render_dir, paste0(sample_name, ".png"))
  html_path <- file.path(render_dir, paste0(sample_name, ".html"))
  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = 16,
    height = 9,
    dpi = 75,
    units = "in"
  )
  script_env$save_html_widget(sankey, html_path)
  dependency_dir <- paste0(tools::file_path_sans_ext(html_path), "_files")
  if (!file.exists(png_path) || file.info(png_path)$size <= 0L ||
      !file.exists(html_path) || file.info(html_path)$size <= 0L ||
      !dir.exists(dependency_dir) ||
      length(list.files(dependency_dir, recursive = TRUE)) == 0L) {
    stop("CS8 simple FLOW rendering is incomplete for ", sample_name, ".", call. = FALSE)
  }
  if (!identical(flow_tbl, flow_before_render)) {
    stop("CS8 FLOW rendering mutated scientific data for ", sample_name, ".", call. = FALSE)
  }
}

message(
  "PASS: CS8 K01563 simple FLOW has aligned PNG geometry and automatic HTML layout"
)
