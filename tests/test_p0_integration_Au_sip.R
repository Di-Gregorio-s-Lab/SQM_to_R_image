source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")
project_dir <- "in/Au_sip"

if (!dir.exists(project_dir)) {
  stop("Integration fixture not found: ", project_dir, call. = FALSE)
}

sqm <- SQMtools::loadSQM(
  project_path = script_env$normalize_sqm_project_dir(project_dir),
  tax_mode = "prokfilter",
  trusted_functions_only = FALSE,
  load_sequences = FALSE
)
script_env$validate_sqm_object(sqm)

expected_mapping <- c(
  "Carbon fixation in photosynthetic organisms" = "00710",
  "Nitrotoluene degradation" = "00633",
  "Nitrogen metabolism" = "00910"
)
resolved_mapping <- script_env$resolve_pathways(sqm, names(expected_mapping))
actual_mapping <- stats::setNames(
  vapply(resolved_mapping, `[[`, character(1), "pathway_id"),
  vapply(resolved_mapping, `[[`, character(1), "canonical_pathway_name")
)
stopifnot(identical(actual_mapping[names(expected_mapping)], expected_mapping))

bacillota <- script_env$resolve_taxa_filters(sqm, "Bacillota")[[1]]
expected_orf_ids <- rownames(sqm$orfs$tax)[
  tolower(as.character(sqm$orfs$tax[, "phylum"])) == "bacillota"
]
stopifnot(identical(bacillota$orf_ids, expected_orf_ids))

bacillota_sqm <- script_env$subset_sqm_by_taxon(
  sqm = sqm,
  orf_ids = bacillota$orf_ids
)
actual_orf_ids <- rownames(bacillota_sqm$orfs$table)
stopifnot(identical(actual_orf_ids, expected_orf_ids))

pathway_sqm <- script_env$subset_pathway(
  sqm,
  "Chlorocyclohexane and chlorobenzene degradation"
)
selected_samples <- colnames(pathway_sqm$orfs$tpm)
taxonomy_table <- script_env$build_pathway_taxonomy_percent_table(
  sqm_object = pathway_sqm,
  rank = "phylum",
  selected_samples = selected_samples,
  top_n_taxa = 15L,
  pathway_name = "Chlorocyclohexane and chlorobenzene degradation"
)

required_columns <- c(
  "sample", "taxon", "taxon_tpm", "pathway_tpm", "value",
  "denominator", "status", "plotted"
)
stopifnot(all(required_columns %in% colnames(taxonomy_table)))

positive_rows <- taxonomy_table[taxonomy_table$status == "ok", , drop = FALSE]
percent_sums <- tapply(positive_rows$value, positive_rows$sample, sum)
stopifnot(length(percent_sums) > 0L)
stopifnot(all(abs(percent_sums - 100) <= 1e-6))
stopifnot(all(positive_rows$denominator == "pathway_tpm_same_sample"))

message(
  "PASS: P0 integration checks completed | Bacillota ORFs=",
  length(expected_orf_ids),
  " | pathway samples=", paste(names(percent_sums), collapse = ",")
)
