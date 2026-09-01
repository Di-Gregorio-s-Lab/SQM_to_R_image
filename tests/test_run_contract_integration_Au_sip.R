source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)
  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")
project_dir <- file.path("in", "Au_sip")
stopifnot(dir.exists(project_dir))

sqm <- script_env$load_sqm_project(project_dir, "prokfilter")
script_env$validate_sqm_object(sqm)
selected_sample <- colnames(sqm$orfs$tpm)[[1L]]

single_taxon <- script_env$resolve_taxa_filters(
  sqm,
  "Accipitriformes (no class in NCBI)"
)[[1L]]
stopifnot(length(single_taxon$orf_ids) == 1L)
single_subset <- script_env$subset_sqm_by_taxon(sqm, single_taxon$orf_ids)
stopifnot(identical(rownames(single_subset$orfs$table), single_taxon$orf_ids))
stopifnot(nrow(single_subset$orfs$table) == 1L)

empty_taxon <- script_env$resolve_taxa_filters(
  sqm,
  "Archaeoglobi"
)[[1L]]
empty_context <- script_env$subset_sqm_by_taxon(sqm, empty_taxon$orf_ids)
script_env$validate_sqm_object(empty_context)
script_env$validate_sqm_tpm_inputs(
  empty_context,
  selected_sample,
  require_kegg_tpm = FALSE
)

resolved_pathways <- script_env$resolve_pathways(sqm, script_env$default_pathway_ids)
pathway_entries <- lapply(resolved_pathways, function(pathway_info) {
  list(
    pathway_name = pathway_info$canonical_pathway_name,
    pathway_id = pathway_info$pathway_id,
    pathway_selection = "defined"
  )
})
warnings <- character()
prepared <- withCallingHandlers(
  script_env$prepare_context_pathway_subsets(
    context_sqm = empty_context,
    pathway_entries = pathway_entries,
    context_label = paste0(empty_taxon$taxon, "@", empty_taxon$rank)
  ),
  warning = function(condition) {
    warnings <<- c(warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
)
stopifnot(nrow(prepared$skips) > 0L)
stopifnot(length(warnings) == nrow(prepared$skips))

output_root <- tempfile("run_contract_real_taxonomy_")
dir.create(output_root, recursive = TRUE)
on.exit(unlink(output_root, recursive = TRUE, force = TRUE), add = TRUE)
taxonomy_result <- script_env$run_taxonomy_scope(
  sqm_object = empty_context,
  output_dir = output_root,
  manifest_base_dir = output_root,
  output_manifests = list(taxon = tibble::tibble()),
  script_name = "sqm_plots.R",
  project_dir = project_dir,
  tax_mode = "prokfilter",
  pathway_name = NA_character_,
  selected_samples = selected_sample,
  dimensions = list("2x2" = c(width = 2, height = 2)),
  plot_dpi = 72,
  top_n_taxa = 15L,
  top_n_ko = 20L,
  taxonomy_ranks = empty_taxon$rank,
  taxonomy_counts = "abund",
  scope_name = "taxonomy_global",
  ignore_unmapped = TRUE,
  ignore_unclassified = TRUE,
  filtered_taxon = empty_taxon$taxon,
  filtered_taxon_rank = empty_taxon$rank
)
stopifnot(nrow(taxonomy_result$taxon) == 2L)
stopifnot(all(file.exists(file.path(output_root, taxonomy_result$taxon$output_file))))

message(
  "PASS: real run-contract edge cases | single_taxon_orfs=",
  nrow(single_subset$orfs$table),
  " | empty_pathways=", nrow(prepared$skips),
  " | taxonomy_artifacts=", nrow(taxonomy_result$taxon)
)
