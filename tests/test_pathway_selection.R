source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")

stopifnot(identical(
  script_env$normalize_pathway_selection_modes(NULL),
  c("defined", "top20")
))
stopifnot(identical(
  script_env$pie_pathway_selection_modes(c("defined", "top20")),
  "defined"
))
stopifnot(identical(
  script_env$pie_pathway_selection_modes("top20", explicitly_requested = TRUE),
  "top20"
))
stopifnot(identical(script_env$default_pathway_ids, script_env$known_pathways$pathway_id))

invalid_selection_error <- tryCatch(
  {
    script_env$normalize_pathway_selection_modes("all")
    NULL
  },
  error = function(error) error$message
)
stopifnot(grepl("pathway_selection_modes", invalid_selection_error, fixed = TRUE))

fake_sqm <- list(
  orfs = list(
    table = data.frame(
      KEGGPATH = c(
        "Metabolism; Test category; Alpha pathway",
        "Metabolism; Test category; Beta pathway",
        "Metabolism; Test category; Gamma pathway",
        "Metabolism; Test category; Alpha pathway"
      ),
      row.names = c("orf_a", "orf_b", "orf_c", "orf_d"),
      check.names = FALSE
    ),
    tpm = data.frame(
      S0 = c(10, 5, 7, 0),
      S1 = c(0, 5, 3, 10),
      row.names = c("orf_a", "orf_b", "orf_c", "orf_d")
    )
  ),
  misc = list(KEGG_paths = c("Alpha pathway", "Beta pathway", "Gamma pathway"))
)

top_pathways <- script_env$select_top_pathways(fake_sqm, c("S0", "S1"), 2L)
stopifnot(identical(
  vapply(top_pathways, `[[`, character(1), "canonical_pathway_name"),
  c("Alpha pathway", "Beta pathway")
))
stopifnot(identical(
  vapply(top_pathways, `[[`, numeric(1), "total_tpm"),
  c(20, 10)
))
stopifnot(all(is.na(vapply(top_pathways, `[[`, character(1), "pathway_id"))))

defined_groups <- script_env$select_pathway_groups(
  fake_sqm,
  requested_pathways = "00361",
  pathway_selection_modes = "defined",
  selected_samples = c("S0", "S1"),
  pathway_top_n = 2L
)
stopifnot(identical(names(defined_groups), "defined"))
stopifnot(identical(defined_groups$defined[[1]]$pathway_id, "00361"))

default_groups <- script_env$select_pathway_groups(
  fake_sqm,
  requested_pathways = character(),
  pathway_selection_modes = NULL,
  selected_samples = c("S0", "S1"),
  pathway_top_n = 2L
)
stopifnot(identical(names(default_groups), c("defined", "top20")))
stopifnot(identical(
  vapply(default_groups$defined, `[[`, character(1), "pathway_id"),
  script_env$default_pathway_ids
))
stopifnot(!script_env$pathview_is_exportable("top20", NA_character_))
stopifnot(!script_env$pathview_is_exportable("defined", NA_character_))

test_root <- tempfile("pathway_selection_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

calls <- list()
fake_export_pathway <- function(
    SQM,
    pathway_id,
    count,
    samples,
    split_samples,
    log_scale,
    output_dir,
    output_suffix) {
  calls[[length(calls) + 1L]] <<- list(
    output_dir = output_dir,
    split_samples = split_samples,
    log_scale = log_scale
  )
  writeLines("fake pathview output", file.path(output_dir, "map.png"))
}

result <- script_env$run_pathview_mode(
  sqm_object = list(
    fake = TRUE,
    functions = list(KEGG = list(tpm = data.frame(
      S0 = c(1, 2),
      S1 = c(3, 4),
      row.names = c("K00001", "K00002")
    )))
  ),
  output_dir = test_root,
  manifest_base_dir = test_root,
  output_manifests = list(pathview = tibble::tibble()),
  script_name = "sqm_plots.R",
  project_dir = "project",
  tax_mode = "prokfilter",
  pathway_name = "Nitrogen metabolism",
  pathway_id = "00910",
  selected_samples = c("S0", "S1"),
  top_n_taxa = 15L,
  top_n_ko = 20L,
  pathway_selection = "top20",
  pathview_sample_modes = c("insieme", "separato"),
  export_pathway_fn = fake_export_pathway
)

stopifnot(length(calls) == 3L)
stopifnot(all(!vapply(calls, `[[`, logical(1), "split_samples")))
stopifnot(all(!vapply(calls, `[[`, logical(1), "log_scale")))
stopifnot(identical(
  sort(unique(result$pathview$output_scope)),
  c("pathway_top20_insieme", "pathway_top20_separato")
))

fake_pie_sqm <- fake_sqm
fake_pie_sqm$orfs$table$"KEGG ID" <- c("K00001", "K00002", "K00003", "K00004")
fake_pie_sqm$orfs$table$KEGGFUN <- c(
  "Function A [EC:1.1.1.1]", "Function B [EC:1.1.1.2]",
  "Function C [EC:1.1.1.3]", "Function D [EC:1.1.1.4]"
)
fake_pie_sqm$orfs$table$KEGGPATH <- "Top pathway"
fake_pie_sqm$orfs$tax <- data.frame(
  superkingdom = rep("Bacteria", 4),
  phylum = rep("Proteobacteria", 4),
  class = rep("Gammaproteobacteria", 4),
  order = rep("Pseudomonadales", 4),
  family = rep("Pseudomonadaceae", 4),
  genus = rep("Pseudomonas", 4),
  species = rep("Pseudomonas sp.", 4),
  row.names = rownames(fake_pie_sqm$orfs$table),
  check.names = FALSE
)
fake_pie_sqm$misc <- list(
  KEGG_names = c(K00001 = "Function A", K00002 = "Function B", K00003 = "Function C", K00004 = "Function D")
)

pie_result <- script_env$run_pie_mode(
  output_dir = test_root,
  manifest_base_dir = file.path(test_root, "pie"),
  output_manifests = list(pie = tibble::tibble()),
  script_name = "sqm_plots.R",
  project_dir = "project",
  tax_mode = "prokfilter",
  pathway_name = "Top pathway",
  pathway_sqm = fake_pie_sqm,
  selected_samples = c("S0", "S1"),
  taxonomy_ranks = "phylum",
  dimensions = list("2x2" = c(width = 2, height = 2)),
  plot_dpi = 72,
  top_n_taxa = 15L,
  top_n_ko = 20L,
  pathway_selection = "top20"
)

stopifnot(dir.exists(file.path(
  test_root,
  "pie",
  "top20",
  script_env$safe_output_component("Top pathway", max_length = 28L)
)))
stopifnot(nrow(pie_result$pie) > 0L)
stopifnot(identical(unique(pie_result$pie$output_scope), "pathway_top20"))

legacy_pie_row <- script_env$new_manifest_row(
  script_name = "sqm_plots.R",
  project_dir = "project",
  tax_mode = "prokfilter",
  pathway = "Legacy pathway",
  samples = "S0",
  metric = "tpm",
  top_n_taxa = 15L,
  top_n_ko = 20L,
  output_type = "plot_png",
  output_file = "pie/Legacy_pathway/S0/legacy.png",
  mode = "pie",
  format = "png",
  output_scope = "pathway"
)
legacy_pie_path <- file.path(test_root, "pie", "manifest_pie.tsv")
script_env$write_tsv_safe(legacy_pie_row |> select(-.data$ec_code), legacy_pie_path)
merged_pie_path <- script_env$write_section_manifest(
  pie_result$pie,
  test_root,
  "pie",
  "manifest_pie.tsv"
)
merged_pie_manifest <- readr::read_tsv(merged_pie_path, show_col_types = FALSE, na = "NA")
legacy_pie_entry <- merged_pie_manifest |> filter(.data$pathway == "Legacy pathway")
stopifnot(nrow(legacy_pie_entry) == 0L)

script_env$write_section_manifest(pie_result$pie, test_root, "pie", "manifest_pie.tsv")
deduplicated_pie_manifest <- readr::read_tsv(merged_pie_path, show_col_types = FALSE, na = "NA")
stopifnot(nrow(deduplicated_pie_manifest) == nrow(merged_pie_manifest))

message("PASS: Pathway selection ranks top pathways and scopes Pathview output")
