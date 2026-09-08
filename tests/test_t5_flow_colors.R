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

find_named_scale <- function(built_plot, title) {
  title_prefix <- strsplit(title, " ", fixed = TRUE)[[1L]][[1L]]
  matches <- Filter(
    function(scale) {
      is.character(scale$name) &&
        length(scale$name) == 1L &&
        grepl(title_prefix, scale$name, fixed = TRUE)
    },
    built_plot$plot$scales$scales
  )
  if (length(matches) != 1L) {
    stop("Expected one FLOW scale titled ", title, ".", call. = FALSE)
  }
  matches[[1L]]
}

extract_scale_colors <- function(built_plot, title, categories) {
  toupper(find_named_scale(built_plot, title)$map(categories))
}

taxon_categories <- levels(flow_fixture$taxon)
ko_categories <- levels(flow_fixture$KO)
expected_taxon_colors <- c("#5D8AA8", "#E32636", "GREY70")
expected_ko_colors <- c("#EFDECD", "#FFBF00", "GREY70")

observed_taxon_colors <- extract_scale_colors(
  built_plot,
  "Taxonomy | % sample",
  taxon_categories
)
observed_ko_colors <- extract_scale_colors(
  built_plot,
  "Function (KO / EC) | % sample",
  ko_categories
)

if (!identical(observed_taxon_colors, expected_taxon_colors) ||
    !identical(observed_ko_colors, expected_ko_colors)) {
  stop(
    "FLOW category mapping is not shared: taxonomy=",
    paste(observed_taxon_colors, collapse = ","),
    "; function=", paste(observed_ko_colors, collapse = ","),
    "; expected taxonomy=", paste(expected_taxon_colors, collapse = ","),
    "; expected function=", paste(expected_ko_colors, collapse = ","),
    call. = FALSE
  )
}

non_other_colors <- c(
  observed_taxon_colors[taxon_categories != "Other"],
  observed_ko_colors[ko_categories != "Other"]
)
if (anyDuplicated(non_other_colors)) {
  stop("Distinct FLOW categories received duplicate colors.", call. = FALSE)
}

shuffled_plot <- script_env$make_flow_plot(
  flow_fixture[c(4L, 2L, 1L, 3L), ],
  "Synthetic pathway",
  "phylum",
  "S_flow"
)
shuffled_build <- ggplot2::ggplot_build(shuffled_plot)
if (!identical(
      extract_scale_colors(shuffled_build, "Taxonomy | % sample", taxon_categories),
      observed_taxon_colors
    ) ||
    !identical(
      extract_scale_colors(shuffled_build, "Function (KO / EC) | % sample", ko_categories),
      observed_ko_colors
    )) {
  stop("FLOW category colors changed after row reordering.", call. = FALSE)
}

sankey <- script_env$make_flow_sankey(
  flow_fixture,
  "Synthetic pathway",
  "phylum",
  "S_flow"
)
sankey_trace <- plotly::plotly_build(sankey)$x$data[[1L]]
observed_node_colors <- toupper(as.character(unlist(
  sankey_trace$node$color,
  use.names = FALSE
)))
expected_node_colors <- c(expected_taxon_colors, expected_ko_colors)
if (!identical(observed_node_colors, expected_node_colors)) {
  stop(
    "FLOW PNG and HTML node colors diverge: expected ",
    paste(expected_node_colors, collapse = ","),
    "; observed ", paste(observed_node_colors, collapse = ","),
    call. = FALSE
  )
}

link_sources <- as.integer(unlist(sankey_trace$link$source, use.names = FALSE)) + 1L
expected_link_colors <- toupper(grDevices::adjustcolor(
  expected_taxon_colors[link_sources],
  alpha.f = 0.65
))
observed_link_colors <- toupper(as.character(unlist(
  sankey_trace$link$color,
  use.names = FALSE
)))
if (!identical(observed_link_colors, expected_link_colors)) {
  stop("FLOW Sankey links do not reuse their source taxon colors.", call. = FALSE)
}

base_color_count <- length(unique(script_env$colors_hex))
overflow_taxa <- paste0("Taxon_", seq_len(base_color_count + 1L))
overflow_fixture <- tibble::tibble(
  pathway = "Overflow pathway",
  rank = "phylum",
  sample = "S_overflow",
  taxon = factor(overflow_taxa, levels = overflow_taxa),
  KO = factor(rep("K_OVERFLOW", length(overflow_taxa)), levels = "K_OVERFLOW"),
  KO_name = "Overflow function",
  ec_codes = NA_character_,
  TPM = 1,
  taxon_percent = 100 / length(overflow_taxa),
  KO_percent = 100,
  flow_percent = 100 / length(overflow_taxa)
)
overflow_build <- ggplot2::ggplot_build(script_env$make_flow_plot(
  overflow_fixture,
  "Overflow pathway",
  "phylum",
  "S_overflow"
))
overflow_colors <- c(
  extract_scale_colors(
    overflow_build,
    "Taxonomy | % sample",
    overflow_taxa
  ),
  extract_scale_colors(
    overflow_build,
    "Function (KO / EC) | % sample",
    "K_OVERFLOW"
  )
)
if (anyNA(overflow_colors) || anyDuplicated(overflow_colors)) {
  stop(
    "FLOW palette overflow produced missing or recycled colors for ",
    length(overflow_colors), " categories.",
    call. = FALSE
  )
}

message("PASS: FLOW uses one category color mapping across PNG and HTML")
