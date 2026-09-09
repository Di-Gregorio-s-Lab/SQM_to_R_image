script_env <- if (identical(Sys.getenv("R_COVR"), "true")) environment() else new.env(parent = globalenv())
if (!identical(Sys.getenv("R_COVR"), "true")) sys.source("sqm_plots_lean.R", envir = script_env)

need <- function(name) {
  if (!exists(name, envir = script_env, mode = "function", inherits = FALSE)) {
    stop("Missing lean renderer helper: ", name, call. = FALSE)
  }
  get(name, envir = script_env, inherits = FALSE)
}

expect_close <- function(actual, expected, label, tolerance = 1e-10) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance))) {
    stop(label, "; expected ", expected, ", observed ", actual, call. = FALSE)
  }
}

artifact_paths <- function(manifest, root) {
  stopifnot(is.data.frame(manifest), "output_file" %in% names(manifest), nrow(manifest) > 0L)
  paths <- as.character(manifest$output_file)
  relative <- !file.exists(paths)
  paths[relative] <- file.path(root, paths[relative])
  stopifnot(all(file.exists(paths)), all(file.info(paths)$size > 0L))
  normalizePath(paths, winslash = "/", mustWork = TRUE)
}

assert_renderer <- function(manifest, root, output_root) {
  paths <- artifact_paths(manifest, root)
  expected_root <- paste0(normalizePath(file.path(root, output_root), winslash = "/", mustWork = TRUE), "/")
  stopifnot(all(startsWith(paste0(dirname(paths), "/"), expected_root)))
  stopifnot(any(grepl("\\.tsv$", paths, ignore.case = TRUE)))
  invisible(paths)
}

allocated <- data.frame(
  orf_id = paste0("orf", 1:4),
  sample = rep("S1", 4L),
  ko_id = rep(c("K00001", "K00002"), each = 2L),
  tpm = c(60, 40, 30, 20),
  taxon = c("Alpha", "Beta", "Alpha", "Beta"),
  phylum = c("Alpha", "Beta", "Alpha", "Beta"),
  kegg_function = rep(c("Primary KO", "Secondary KO"), each = 2L),
  ec_codes = rep(c("1.1.1.1", "2.2.2.2"), each = 2L),
  pathway_id = "00001",
  pathway_name = "Synthetic pathway",
  pathway_root = "Metabolism",
  pathway_category = "Synthetic category",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

official_totals <- data.frame(
  sample = c("S1", "S1"),
  ko_id = c("K00001", "K00002"),
  tpm = c(100, 50),
  stringsAsFactors = FALSE
)

funz <- need("build_funz_table")(
  allocated,
  official_totals,
  samples = "S1",
  top_n = 1L
)
stopifnot(all(c("sample", "ko_id", "tpm", "percent") %in% names(funz)))
expect_close(sum(funz$tpm), 150, "FUNZ changed official KO mass")
expect_close(sum(funz$percent), 100, "FUNZ percentages do not close")
stopifnot("Other" %in% as.character(funz$ko_id))

funz_plot <- need("make_funz_plot")(funz, "Synthetic pathway")
funz_scale <- ggplot2::ggplot_build(funz_plot)$plot$scales$get_scales("fill")
stopifnot(
  identical(funz_plot$labels$title, "KO barplot - Synthetic pathway"),
  identical(funz_plot$labels$x, "Sample"),
  identical(funz_plot$labels$y, "TPM"),
  identical(funz_plot$labels$fill, "KO / EC | % per sample"),
  identical(toupper(unname(funz_scale$map("K00001"))), toupper(colors_hex[[1L]])),
  identical(tolower(unname(funz_scale$map("Other"))), "grey70")
)

enzyme_style_table <- data.frame(
  sample = factor(c("S1", "S2"), levels = c("S1", "S2")),
  ec_code = factor(c("1.1.1.1", "2.2.2.2"), levels = c("1.1.1.1", "2.2.2.2")),
  tpm = c(10, 20)
)
enzyme_bar <- need("make_enzyme_barplot")(enzyme_style_table, "Enzyme barplot - combined")
enzyme_line <- need("make_enzyme_lineplot")(enzyme_style_table, "Enzyme line chart - combined")
stopifnot(
  identical(enzyme_bar$labels$x, "Sample"),
  identical(enzyme_bar$labels$y, "TPM"),
  identical(enzyme_bar$labels$fill, "EC"),
  identical(enzyme_line$labels$x, "Sample"),
  identical(enzyme_line$labels$y, "TPM"),
  identical(enzyme_line$labels$colour, "EC"),
  identical(
    toupper(unname(ggplot2::ggplot_build(enzyme_bar)$plot$scales$get_scales("fill")$map("1.1.1.1"))),
    toupper(colors_hex[[1L]])
  )
)

allocated_many <- rbind(
  allocated,
  transform(allocated[1L, ], orf_id = "orf5", ko_id = "K00003", tpm = 20,
            kegg_function = "KO without EC", ec_codes = NA_character_)
)
official_many <- rbind(
  official_totals,
  data.frame(sample = "S1", ko_id = "K00003", tpm = 20)
)
funz_many <- need("build_funz_table")(allocated_many, official_many, "S1", 1L)
other <- funz_many[funz_many$ko_id == "Other", , drop = FALSE]
stopifnot(nrow(other) == 1L)
expect_close(other$tpm, 70, "FUNZ did not collapse every non-top KO into one Other row")

flow <- need("build_flow_table")(
  allocated,
  samples = "S1",
  rank = "phylum",
  top_n_ko = 1L,
  top_n_taxa = 1L
)
stopifnot(all(c(
  "sample", "ko_id", "taxon", "tpm", "percent", "kegg_function",
  "ec_codes", "taxon_percent", "ko_percent"
) %in% names(flow)))
expect_close(sum(flow$tpm), 150, "FLOW changed allocated mass")
expect_close(sum(flow$percent), 100, "FLOW percentages do not close")
stopifnot(
  identical(unique(flow$kegg_function[flow$ko_id == "K00001"]), "Primary KO"),
  identical(unique(flow$ec_codes[flow$ko_id == "K00001"]), "1.1.1.1"),
  identical(unique(flow$kegg_function[flow$ko_id == "Other"]), "Other KOs"),
  all(is.na(flow$ec_codes[flow$ko_id == "Other"]))
)
expect_close(unique(flow$taxon_percent[flow$taxon == "Alpha"]), 60,
             "FLOW taxonomy share changed")
expect_close(unique(flow$ko_percent[flow$ko_id == "K00001"]), 100 / 1.5,
             "FLOW KO share changed")

taxon_levels <- unique(flow$taxon)
ko_levels <- unique(flow$ko_id)
flow_palette <- need("build_flow_color_map")(taxon_levels, ko_levels)
stopifnot(
  identical(unname(flow_palette["Alpha"]), "#5d8aa8"),
  identical(unname(flow_palette["K00001"]), "#e32636"),
  identical(unname(flow_palette["Other"]), "grey70")
)

flow_plot <- need("make_flow_plot")(flow, "Synthetic pathway", "phylum", "S1")
flow_built <- ggplot2::ggplot_build(flow_plot)
flow_strata <- flow_built$data[[2L]]
expected_strata_colors <- toupper(unname(flow_palette[as.character(flow_strata$stratum)]))
stopifnot(
  !anyNA(expected_strata_colors),
  identical(toupper(as.character(flow_strata$fill)), expected_strata_colors),
  !any(toupper(as.character(flow_strata$fill)) %in% c("WHITE", "#FFFFFF"))
)
fill_scale <- flow_built$plot$scales$get_scales("fill")
colour_scale <- flow_built$plot$scales$get_scales("colour")
stopifnot(
  identical(flow_plot$labels$title, "Flowplot - Synthetic pathway - phylum - S1"),
  identical(flow_plot$labels$subtitle, "Taxonomy-to-KO flow from the ORF x sample x KO table"),
  identical(flow_plot$labels$y, "Relative flow (%)"),
  identical(sort(unique(as.character(flow_built$data[[3L]]$label))),
            sort(c("Alpha", "Other", "K00001"))),
  isTRUE(all.equal(max(flow_built$data[[1L]]$ymax), 100, tolerance = 1e-10)),
  identical(fill_scale$name, "Taxonomy"),
  identical(fill_scale$labels, c("Alpha | 60.0%", "Other | 40.0%")),
  identical(colour_scale$name, "Function (KO / EC)"),
  identical(colour_scale$labels,
            c("K00001 / EC 1.1.1.1 | 66.7%", "Other KOs | 33.3%")),
  identical(flow_plot$theme$legend.position, "right")
)
unclassified <- allocated
unclassified$phylum <- NA_character_
unclassified_flow <- need("build_flow_table")(unclassified, "S1", "phylum", 2L, 2L)
stopifnot(identical(unique(unclassified_flow$taxon), "Unclassified"))

pie <- need("build_pie_table")(
  allocated,
  official_totals,
  samples = "S1",
  ko_id = "K00001",
  rank = "phylum",
  top_n_taxa = 10L
)
stopifnot(all(c("sample", "ko_id", "taxon", "tpm", "pct") %in% names(pie)))
expect_close(sum(pie$tpm), 100, "PIE changed official KO mass")
expect_close(sum(pie$pct), 1, "PIE proportions do not close")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required by the production renderers.", call. = FALSE)
}

output_dir <- tempfile("lean_renderers_")
dir.create(output_dir, recursive = TRUE)
on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
context <- list(
  output_dir = output_dir,
  selection = "defined",
  filtered_taxon = NA_character_,
  filtered_taxon_rank = NA_character_
)
render_args <- list(
  context = context,
  pathway_id = "00001",
  pathway_name = "Synthetic pathway",
  width = 2,
  height = 2,
  dpi = 72
)

funz_manifest <- do.call(need("render_funz"), c(list(table = funz), render_args))
flow_manifest <- do.call(need("render_flow"), c(list(table = flow), render_args))
pie_manifest <- do.call(need("render_pie"), c(list(table = pie), render_args))
assert_renderer(funz_manifest, output_dir, "funz")
flow_paths <- assert_renderer(flow_manifest, output_dir, "flowplot")
assert_renderer(pie_manifest, output_dir, "pie")
flow_tsv <- utils::read.delim(
  flow_paths[grepl("\\.tsv$", flow_paths, ignore.case = TRUE)][[1L]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  all(c("kegg_function", "ec_codes", "taxon_percent", "ko_percent") %in% names(flow_tsv)),
  isTRUE(all.equal(sum(flow_tsv$flow_percent), 100, tolerance = 1e-10))
)

if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
  flow_widget <- need("make_flow_sankey")(flow, "Synthetic pathway", "phylum", "S1")
  flow_trace <- plotly::plotly_build(flow_widget)$x$data[[1L]]
  expected_node_colors <- toupper(unname(flow_palette[c(taxon_levels, ko_levels)]))
  observed_node_colors <- toupper(as.character(unlist(flow_trace$node$color, use.names = FALSE)))
  expected_link_colors <- toupper(unname(grDevices::adjustcolor(
    flow_palette[flow$taxon],
    alpha.f = 0.65
  )))
  observed_link_colors <- toupper(as.character(unlist(flow_trace$link$color, use.names = FALSE)))
  stopifnot(
    identical(observed_node_colors, expected_node_colors),
    identical(observed_link_colors, expected_link_colors),
    identical(as.character(unlist(flow_trace$node$label, use.names = FALSE)), c(
      "Alpha | 60.0%", "Other | 40.0%",
      "K00001 / EC 1.1.1.1 | 66.7%", "Other KOs | 33.3%"
    )),
    isTRUE(all.equal(as.numeric(unlist(flow_trace$link$value)), flow$flow_percent,
                     tolerance = 1e-10)),
    identical(as.character(flow_trace$valueformat), ".2f"),
    identical(as.character(flow_trace$valuesuffix), "%"),
    all(grepl("Taxon share:|Function share:|TPM:|Flow:",
              as.character(unlist(flow_trace$link$customdata)))),
    any(grepl("Function: Primary KO", as.character(unlist(flow_trace$link$customdata)), fixed = TRUE)),
    identical(
      as.character(plotly::plotly_build(flow_widget)$x$layout$title$text),
      "Flowplot - Synthetic pathway - phylum - S1<br><sup>Taxonomy-to-KO flow based on TPM</sup>"
    )
  )

  flow_html <- do.call(
    need("render_flow_html"),
    c(list(table = flow), render_args[c("context", "pathway_id", "pathway_name")])
  )
  html_paths <- artifact_paths(flow_html, output_dir)
  stopifnot(
    all(startsWith(html_paths, paste0(normalizePath(file.path(output_dir, "flowplot"), winslash = "/"), "/"))),
    any(grepl("\\.html$", flow_html$output_file, ignore.case = TRUE))
  )
}

saved_plot <- file.path(output_dir, "direct", "small.png")
need("save_plot")(
  ggplot2::ggplot(data.frame(x = "A", y = 1), ggplot2::aes(x, y)) + ggplot2::geom_col(),
  saved_plot,
  width = 2,
  height = 2,
  dpi = 72
)
stopifnot(file.exists(saved_plot), file.info(saved_plot)$size > 0L)

taxonomy_data <- data.frame(taxon = c("Alpha", "Beta"), tpm = c(90, 60))
taxonomy_call <- NULL
fake_plot_taxonomy <- function(SQM, rank, samples, rescale, ...) {
  taxonomy_call <<- list(SQM = SQM, rank = rank, samples = samples, rescale = rescale)
  ggplot2::ggplot(taxonomy_data, ggplot2::aes(taxon, tpm)) + ggplot2::geom_col()
}
taxonomy_manifest <- need("render_taxonomy")(
  sqm_object = list(marker = "synthetic"),
  samples = "S1",
  rank = "phylum",
  context = context,
  pathway_id = "00001",
  pathway_name = "Synthetic pathway",
  width = 2,
  height = 2,
  dpi = 72,
  plot_fun = fake_plot_taxonomy
)
taxonomy_paths <- assert_renderer(taxonomy_manifest, output_dir, "taxonomy_by_pathway")
stopifnot(identical(taxonomy_call$rescale, FALSE))
taxonomy_tsv <- utils::read.delim(
  taxonomy_paths[grepl("\\.tsv$", taxonomy_paths, ignore.case = TRUE)][[1L]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(isTRUE(all.equal(taxonomy_tsv[names(taxonomy_data)], taxonomy_data, check.attributes = FALSE)))

export_call <- NULL
fake_pathview_export <- function(SQM, pathway_id, count, samples, split_samples,
                                 log_scale, output_dir, output_suffix) {
  export_call <<- list(
    output_dir = output_dir,
    pathway_id = pathway_id,
    samples = samples,
    split_samples = split_samples
  )
  writeLines("synthetic pathview", file.path(output_dir, paste0(output_suffix, ".png")))
}
pathview_dir <- file.path(output_dir, "pathview", "definiti", "00001")
pathview_result <- need("export_pathview_isolated")(
  fake_pathview_export,
  list(marker = "synthetic"),
  "00001",
  "S1",
  pathview_dir
)
pathview_files <- if (is.data.frame(pathview_result)) {
  artifact_paths(pathview_result, output_dir)
} else {
  normalizePath(as.character(pathview_result), winslash = "/", mustWork = TRUE)
}
stopifnot(
  length(pathview_files) > 0L,
  all(file.info(pathview_files)$size > 0L),
  all(startsWith(pathview_files, paste0(normalizePath(pathview_dir, winslash = "/"), "/"))),
  !identical(normalizePath(export_call$output_dir, winslash = "/"), normalizePath(pathview_dir, winslash = "/")),
  !dir.exists(export_call$output_dir),
  identical(export_call$split_samples, FALSE)
)

message("PASS: lean renderer contracts preserve mass and write isolated artifacts")
