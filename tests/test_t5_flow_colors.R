source("tests/helpers/t1_sqmtools_oracles.R")

script_env <- t1_source_sqm_plots_without_main()

flow_fixture <- tibble::tibble(
  pathway = "Synthetic pathway",
  rank = "phylum",
  sample = "S_flow",
  taxon = factor(
    c("Alpha", "Beta", "Alpha", "Other"),
    levels = c("Beta", "Alpha", "Other")
  ),
  KO = factor(
    c("K00001", "K00001", "K00002", "Other"),
    levels = c("K00002", "K00001", "Other")
  ),
  KO_name = c("Function one", "Function one", "Function two", "Other KOs"),
  ec_codes = c("1.1.1.1", "1.1.1.1", "2.2.2.2", NA_character_),
  TPM = c(35, 20, 35, 10),
  taxon_percent = c(70, 20, 70, 10),
  KO_percent = c(55, 55, 35, 10),
  flow_percent = c(35, 20, 35, 10)
)

extract_layer_colors <- function(layer_data, x_value = NULL) {
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

plot_object <- script_env$make_flow_plot(
  flow_fixture,
  "Synthetic pathway",
  "phylum",
  "S_flow"
)
built_plot <- ggplot2::ggplot_build(plot_object)

flow_colors <- extract_layer_colors(built_plot$data[[1L]], x_value = 1)
taxonomy_stratum_colors <- extract_layer_colors(built_plot$data[[2L]])

taxonomy_scale <- Filter(
  function(scale) {
    is.character(scale$name) &&
      length(scale$name) == 1L &&
      grepl("Taxonomy", scale$name, fixed = TRUE)
  },
  built_plot$plot$scales$scales
)[[1L]]
legend_colors <- data.frame(
  category = levels(flow_fixture$taxon),
  color = toupper(taxonomy_scale$map(levels(flow_fixture$taxon))),
  stringsAsFactors = FALSE
)

observed <- merge(
  merge(flow_colors, taxonomy_stratum_colors, by = "category", suffixes = c("_flow", "_stratum")),
  legend_colors,
  by = "category"
)

if (!all(observed$color_flow == observed$color_stratum &
         observed$color_flow == observed$color)) {
  details <- apply(
    observed,
    1L,
    function(row) paste0(
      row[["category"]],
      ": flow=", row[["color_flow"]],
      ", stratum=", row[["color_stratum"]],
      ", legend=", row[["color"]]
    )
  )
  stop(
    "FLOW taxonomy colors diverge between flow, stratum, and legend: ",
    paste(details, collapse = "; "),
    call. = FALSE
  )
}

message("PASS: FLOW taxonomy stratum uses the flow and legend colors")
