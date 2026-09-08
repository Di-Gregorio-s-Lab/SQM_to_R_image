source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

project_dir <- Sys.getenv("SQM_CS8_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) {
  stop("SQM_CS8_PROJECT_DIR is required.", call. = FALSE)
}
if (!dir.exists(project_dir)) {
  stop("SQM_CS8_PROJECT_DIR does not exist: ", project_dir, call. = FALSE)
}

extract_rendered_colors <- function(layer_data, x_value = NULL) {
  fill_column <- grep("^fill", names(layer_data), value = TRUE)
  if (length(fill_column) != 1L) {
    stop("Expected exactly one rendered fill column.", call. = FALSE)
  }
  if (!is.null(x_value)) {
    layer_data <- layer_data[layer_data$x == x_value, , drop = FALSE]
  }
  result <- data.frame(
    category = as.character(layer_data$stratum),
    color = toupper(as.character(layer_data[[fill_column]])),
    stringsAsFactors = FALSE
  )
  unique(result)
}

find_scale <- function(built_plot, prefix) {
  matches <- Filter(
    function(scale) {
      is.character(scale$name) &&
        length(scale$name) == 1L &&
        grepl(prefix, scale$name, fixed = TRUE)
    },
    built_plot$plot$scales$scales
  )
  if (length(matches) != 1L) {
    stop("Expected one FLOW scale matching ", prefix, ".", call. = FALSE)
  }
  matches[[1L]]
}

expect_same_category_colors <- function(observed, expected, label) {
  comparison <- merge(
    observed,
    expected,
    by = "category",
    suffixes = c("_observed", "_expected")
  )
  if (nrow(comparison) != nrow(expected) ||
      !all(comparison$color_observed == comparison$color_expected)) {
    stop(label, call. = FALSE)
  }
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
  rank = "phylum",
  selected_samples = selected_samples,
  top_n_taxa = 999L,
  top_n_ko = 1L,
  ko_lookup = analysis$ko_lookup
)

render_dir <- tempfile("t5_cs8_flow_colors_")
dir.create(render_dir, recursive = TRUE)
on.exit(unlink(render_dir, recursive = TRUE, force = TRUE), add = TRUE)

for (sample_name in selected_samples) {
  flow_tbl <- script_env$build_flow_table_for_sample(
    flow_rank,
    analysis$pathway_name,
    "phylum",
    sample_name
  )
  flow_before_render <- flow_tbl
  observed_tpm <- sum(flow_tbl$TPM[as.character(flow_tbl$KO) == "K01563"])
  t1_expect_equal(
    observed_tpm,
    unname(expected_k01563[[sample_name]]),
    paste0("CS8 FLOW K01563 TPM changed for ", sample_name),
    tolerance = 1e-8
  )

  plot_object <- script_env$make_flow_plot(
    flow_tbl,
    analysis$pathway_name,
    "phylum",
    sample_name
  )
  built_plot <- ggplot2::ggplot_build(plot_object)
  taxon_categories <- levels(flow_tbl$taxon)
  ko_categories <- levels(flow_tbl$KO)

  taxon_scale <- find_scale(built_plot, "Taxonomy")
  ko_scale <- find_scale(built_plot, "Function")
  expected_taxon <- data.frame(
    category = taxon_categories,
    color = toupper(taxon_scale$map(taxon_categories)),
    stringsAsFactors = FALSE
  )
  expected_ko <- data.frame(
    category = ko_categories,
    color = toupper(ko_scale$map(ko_categories)),
    stringsAsFactors = FALSE
  )

  if (anyNA(expected_taxon$color) ||
      any(expected_taxon$color %in% c("WHITE", "#FFFFFF"))) {
    stop("CS8 FLOW taxonomy scale contains missing or white colors for ", sample_name, ".", call. = FALSE)
  }
  expect_same_category_colors(
    extract_rendered_colors(built_plot$data[[1L]], x_value = 1),
    expected_taxon,
    paste0("CS8 FLOW alluvia diverge from taxonomy colors for ", sample_name)
  )
  expect_same_category_colors(
    extract_rendered_colors(built_plot$data[[2L]]),
    expected_taxon,
    paste0("CS8 FLOW taxonomy strata diverge from legend for ", sample_name)
  )
  expect_same_category_colors(
    extract_rendered_colors(built_plot$data[[3L]]),
    expected_ko,
    paste0("CS8 FLOW KO strata diverge from legend for ", sample_name)
  )

  sankey <- script_env$make_flow_sankey(
    flow_tbl,
    analysis$pathway_name,
    "phylum",
    sample_name
  )
  trace <- plotly::plotly_build(sankey)$x$data[[1L]]
  expected_node_colors <- c(expected_taxon$color, expected_ko$color)
  observed_node_colors <- toupper(as.character(unlist(trace$node$color, use.names = FALSE)))
  if (!identical(observed_node_colors, expected_node_colors)) {
    stop("CS8 FLOW PNG and HTML colors diverge for ", sample_name, ".", call. = FALSE)
  }
  link_sources <- as.integer(unlist(trace$link$source, use.names = FALSE)) + 1L
  expected_link_colors <- toupper(grDevices::adjustcolor(
    expected_node_colors[link_sources],
    alpha.f = 0.65
  ))
  observed_link_colors <- toupper(as.character(unlist(trace$link$color, use.names = FALSE)))
  if (!identical(observed_link_colors, expected_link_colors)) {
    stop("CS8 FLOW Sankey links diverge from source taxa for ", sample_name, ".", call. = FALSE)
  }

  png_path <- file.path(render_dir, paste0(sample_name, ".png"))
  html_path <- file.path(render_dir, paste0(sample_name, ".html"))
  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = 6,
    height = 4,
    dpi = 75,
    units = "in"
  )
  script_env$save_html_widget(sankey, html_path)
  if (!file.exists(png_path) || file.info(png_path)$size <= 0L ||
      !file.exists(html_path) || file.info(html_path)$size <= 0L) {
    stop("CS8 FLOW rendering did not create non-empty outputs for ", sample_name, ".", call. = FALSE)
  }
  if (!identical(flow_tbl, flow_before_render)) {
    stop("CS8 FLOW rendering mutated scientific data for ", sample_name, ".", call. = FALSE)
  }
}

message(
  "PASS: CS8 K01563 FLOW colors match across alluvia, strata, legends, PNG and HTML"
)
