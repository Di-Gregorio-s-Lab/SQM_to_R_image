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
pathway_name <- "Xylene degradation"
pathway_id <- "00622"
pathway_sqm <- script_env$subset_pathway(sqm, pathway_name)
native_tpm <- as.matrix(
  pathway_sqm$functions$KEGG$tpm[, selected_samples, drop = FALSE]
)

analysis <- script_env$build_pathway_analysis(
  list(
    pathway_name = pathway_name,
    pathway_id = pathway_id,
    pathway_selection = "defined",
    pathway_sqm = pathway_sqm,
    context_sqm = sqm,
    pathway_ko_ids = rownames(native_tpm)
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

funz <- script_env$build_ko_plot_table(
  orf_long = analysis$orf_long_result$data,
  selected_samples = selected_samples,
  top_n_ko = 20L,
  ko_lookup = analysis$ko_lookup,
  pathway_name = pathway_name
)
funz_tpm <- xtabs(tpm ~ ko_id + sample, data = funz)
t1_expect_equal(
  unname(funz_tpm[rownames(native_tpm), selected_samples, drop = FALSE]),
  unname(native_tpm),
  "CS8 FUNZ margins differ from subsetFun()",
  tolerance = 1e-8
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
  flow_tpm <- tapply(
    flow_tbl$TPM,
    factor(as.character(flow_tbl$KO), levels = rownames(native_tpm)),
    sum,
    default = 0
  )
  t1_expect_equal(
    unname(flow_tpm),
    unname(native_tpm[, sample_name]),
    paste0("CS8 FLOW margins differ from subsetFun() for ", sample_name),
    tolerance = 1e-8
  )
  pie_tpm <- vapply(rownames(native_tpm), function(ko_id) {
    sum(script_env$build_pie_chart_table(
      orf_long = analysis$orf_long_result$data,
      sample_name = sample_name,
      ko_id_filter = ko_id,
      rank_name = "family",
      top_n_taxa = 10L,
      pathway_name = pathway_name,
      pathway_id = pathway_id,
      pathway_selection = "defined"
    )$tpm)
  }, numeric(1))
  t1_expect_equal(
    unname(pie_tpm),
    unname(native_tpm[, sample_name]),
    paste0("CS8 PIE margins differ from subsetFun() for ", sample_name),
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
  if (!identical(plot_object$theme$legend.position, "right")) {
    stop("CS8 FLOW PNG did not render its legends for ", sample_name, ".", call. = FALSE)
  }
  flow_palette <- script_env$build_flow_color_map(
    levels(flow_tbl$taxon),
    levels(flow_tbl$KO)
  )
  strata <- built_plot$data[[2L]]
  expected_strata_colors <- toupper(unname(
    flow_palette[as.character(strata$stratum)]
  ))
  if (anyNA(expected_strata_colors) || !identical(
    toupper(as.character(strata$fill)),
    expected_strata_colors
  )) {
    stop("CS8 FLOW columns do not use shared colors for ", sample_name, ".", call. = FALSE)
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
  expected_node_colors <- toupper(unname(flow_palette[c(
    levels(flow_tbl$taxon),
    levels(flow_tbl$KO)
  )]))
  if (!identical(
    toupper(as.character(unlist(trace$node$color, use.names = FALSE))),
    expected_node_colors
  )) {
    stop("CS8 FLOW HTML and PNG colors differ for ", sample_name, ".", call. = FALSE)
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
  "PASS: CS8 FUNZ/FLOW/PIE match subsetFun(); FLOW rendering is aligned"
)
