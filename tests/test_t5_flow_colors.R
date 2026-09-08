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

taxon_categories <- levels(flow_fixture$taxon)
taxon_palette <- script_env$build_flow_taxon_colors(taxon_categories)
if (!identical(names(taxon_palette), taxon_categories)) {
  stop("Simple FLOW palette does not follow taxonomy factor order.", call. = FALSE)
}
if (!identical(unname(taxon_palette[["Other"]]), "grey70")) {
  stop("Simple FLOW palette changed the reserved Other color.", call. = FALSE)
}
if (anyDuplicated(toupper(unname(taxon_palette[names(taxon_palette) != "Other"])))) {
  stop("Simple FLOW palette duplicated taxonomy colors.", call. = FALSE)
}

plot_object <- script_env$make_flow_plot(
  flow_fixture,
  "Synthetic pathway",
  "phylum",
  "S_flow"
)
built_plot <- ggplot2::ggplot_build(plot_object)
alluvium_fill <- grep("^fill", names(built_plot$data[[1L]]), value = TRUE)
rendered_flows <- built_plot$data[[1L]] |>
  dplyr::filter(.data$x == 1) |>
  dplyr::transmute(
    category = as.character(.data$stratum),
    color = toupper(as.character(.data[[alluvium_fill]]))
  ) |>
  dplyr::distinct()
expected_flows <- data.frame(
  category = taxon_categories,
  color = toupper(unname(taxon_palette[taxon_categories])),
  stringsAsFactors = FALSE
)
comparison <- merge(rendered_flows, expected_flows, by = "category")
if (nrow(comparison) != length(taxon_categories) ||
    !all(comparison$color.x == comparison$color.y)) {
  stop("Simple FLOW ribbons do not use the taxonomy palette.", call. = FALSE)
}

stratum_colors <- unique(toupper(as.character(built_plot$data[[2L]]$fill)))
if (!identical(stratum_colors, "GREY95")) {
  stop("Simple FLOW columns must remain neutral grey95.", call. = FALSE)
}
if (!identical(plot_object$theme$legend.position, "none")) {
  stop("Simple FLOW PNG unexpectedly rendered a legend.", call. = FALSE)
}

sankey <- script_env$make_flow_sankey(
  flow_fixture,
  "Synthetic pathway",
  "phylum",
  "S_flow"
)
trace <- plotly::plotly_build(sankey)$x$data[[1L]]
node_colors <- unique(as.character(unlist(trace$node$color, use.names = FALSE)))
if (!identical(node_colors, "rgba(245,245,245,1)")) {
  stop("Simple FLOW HTML nodes must remain neutral.", call. = FALSE)
}
link_sources <- as.integer(unlist(trace$link$source, use.names = FALSE)) + 1L
expected_link_colors <- toupper(grDevices::adjustcolor(
  unname(taxon_palette[taxon_categories])[link_sources],
  alpha.f = 0.65
))
observed_link_colors <- toupper(as.character(unlist(trace$link$color, use.names = FALSE)))
if (!identical(observed_link_colors, expected_link_colors)) {
  stop("Simple FLOW HTML links do not reuse taxonomy colors.", call. = FALSE)
}

base_color_count <- length(unique(script_env$colors_hex))
overflow_categories <- paste0("Taxon_", seq_len(base_color_count + 1L))
overflow_palette <- script_env$build_flow_taxon_colors(overflow_categories)
if (anyNA(overflow_palette) || anyDuplicated(toupper(unname(overflow_palette)))) {
  stop("Simple FLOW palette overflow produced missing or recycled colors.", call. = FALSE)
}

message("PASS: simple FLOW uses taxonomy colors only for ribbons")
