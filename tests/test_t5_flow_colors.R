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

taxon_levels <- levels(flow_fixture$taxon)
ko_levels <- levels(flow_fixture$KO)
expected_taxon_labels <- c("Beta | 20.0%", "Alpha | 70.0%", "Other | 10.0%")
expected_ko_labels <- c(
  "K00002 / EC 2.2.2.2 | 35.0%",
  "K00001 / EC 1.1.1.1 | 55.0%",
  "Other KOs | 10.0%"
)

flow_palette <- script_env$build_flow_color_map(taxon_levels, ko_levels)
expected_categories <- c("Beta", "Alpha", "K00002", "K00001", "Other")
if (!identical(names(flow_palette), expected_categories)) {
  stop("FLOW palette does not contain taxonomy then function categories.", call. = FALSE)
}
if (!identical(unname(flow_palette[["Other"]]), "grey70")) {
  stop("FLOW palette changed the reserved Other color.", call. = FALSE)
}
if (anyDuplicated(toupper(unname(flow_palette[names(flow_palette) != "Other"])))) {
  stop("FLOW palette assigns one color to multiple displayed categories.", call. = FALSE)
}

plot_object <- script_env$make_flow_plot(
  flow_fixture,
  "Synthetic pathway",
  "phylum",
  "S_flow"
)
built_plot <- ggplot2::ggplot_build(plot_object)

strata <- built_plot$data[[2L]] |>
  dplyr::transmute(
    x = as.numeric(.data$x),
    stratum = as.character(.data$stratum),
    fill = toupper(as.character(.data$fill)),
    ymin = .data$ymin,
    ymax = .data$ymax
  )
expected_strata_fill <- toupper(unname(flow_palette[strata$stratum]))
if (anyNA(expected_strata_fill) || !identical(strata$fill, expected_strata_fill)) {
  stop("FLOW strata do not use the shared category colors.", call. = FALSE)
}
if (any(strata$fill %in% c("WHITE", "#FFFFFF", "GREY95"))) {
  stop("FLOW still contains neutral taxonomy or function strata.", call. = FALSE)
}

alluvia <- built_plot$data[[1L]] |>
  dplyr::transmute(
    x = as.numeric(.data$x),
    stratum = as.character(.data$stratum),
    fill = toupper(as.character(.data$fill)),
    ymin = .data$ymin,
    ymax = .data$ymax
  )
taxon_ribbons <- dplyr::filter(alluvia, .data$x == 1)
if (!identical(
  taxon_ribbons$fill,
  toupper(unname(flow_palette[taxon_ribbons$stratum]))
)) {
  stop("FLOW ribbons do not use their origin taxonomy color.", call. = FALSE)
}

for (axis_position in c(1, 2)) {
  ribbon_bounds <- alluvia |>
    dplyr::filter(.data$x == axis_position) |>
    dplyr::group_by(.data$stratum) |>
    dplyr::summarise(ymin = min(.data$ymin), ymax = max(.data$ymax), .groups = "drop")
  stratum_bounds <- strata |>
    dplyr::filter(.data$x == axis_position) |>
    dplyr::select("stratum", "ymin", "ymax")
  geometry_difference <- dplyr::inner_join(
    ribbon_bounds,
    stratum_bounds,
    by = "stratum",
    suffix = c("_ribbon", "_stratum")
  ) |>
    dplyr::summarise(
      difference = max(
        abs(.data$ymin_ribbon - .data$ymin_stratum),
        abs(.data$ymax_ribbon - .data$ymax_stratum)
      )
    ) |>
    dplyr::pull("difference")
  if (!isTRUE(geometry_difference < 1e-10)) {
    stop("FLOW color layers changed the aligned wide-table geometry.", call. = FALSE)
  }
}

fill_scale <- built_plot$plot$scales$get_scales("fill")
colour_scale <- built_plot$plot$scales$get_scales("colour")
if (!identical(fill_scale$name, "Taxonomy") ||
    !identical(fill_scale$breaks, taxon_levels) ||
    !identical(fill_scale$labels, expected_taxon_labels)) {
  stop("FLOW taxonomy legend title, order, or labels changed.", call. = FALSE)
}
if (!identical(colour_scale$name, "Function (KO / EC)") ||
    !identical(colour_scale$breaks, ko_levels) ||
    !identical(colour_scale$labels, expected_ko_labels)) {
  stop("FLOW function legend title, order, or labels changed.", call. = FALSE)
}
if (!identical(plot_object$theme$legend.position, "right")) {
  stop("FLOW PNG does not display its two legends on the right.", call. = FALSE)
}

sankey <- script_env$make_flow_sankey(
  flow_fixture,
  "Synthetic pathway",
  "phylum",
  "S_flow"
)
trace <- sankey$x$data[[1L]]
expected_node_colors <- toupper(unname(flow_palette[c(taxon_levels, ko_levels)]))
observed_node_colors <- toupper(as.character(unlist(trace$node$color, use.names = FALSE)))
if (!identical(observed_node_colors, expected_node_colors)) {
  stop("FLOW HTML nodes do not share the PNG category colors.", call. = FALSE)
}
expected_node_labels <- c(expected_taxon_labels, expected_ko_labels)
observed_node_labels <- as.character(unlist(trace$node$label, use.names = FALSE))
if (!identical(observed_node_labels, expected_node_labels)) {
  stop("FLOW HTML nodes do not share the PNG legend labels.", call. = FALSE)
}

link_sources <- as.integer(unlist(trace$link$source, use.names = FALSE)) + 1L
expected_link_colors <- toupper(unname(grDevices::adjustcolor(
  flow_palette[taxon_levels[link_sources]],
  alpha.f = 0.65
)))
observed_link_colors <- toupper(as.character(unlist(trace$link$color, use.names = FALSE)))
if (!identical(observed_link_colors, expected_link_colors)) {
  stop("FLOW HTML links do not reuse their origin taxonomy color.", call. = FALSE)
}

base_color_count <- length(unique(script_env$colors_hex))
overflow_categories <- paste0("Category_", seq_len(base_color_count + 1L))
overflow_palette <- script_env$build_flow_color_map(
  overflow_categories,
  character()
)
if (anyNA(overflow_palette) || anyDuplicated(toupper(unname(overflow_palette)))) {
  stop("FLOW palette overflow produced missing or recycled colors.", call. = FALSE)
}

message("PASS: FLOW shares aligned colors and labels across PNG and HTML")
