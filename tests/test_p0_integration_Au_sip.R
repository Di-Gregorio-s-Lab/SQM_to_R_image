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
stopifnot(all(positive_rows$denominator == positive_rows$pathway_tpm))

taxonomy_warnings <- character()
global_plot <- withCallingHandlers(
  script_env$make_taxonomy_plot(
    sqm_object = sqm,
    rank = "phylum",
    count = "percent",
    selected_samples = colnames(sqm$orfs$tpm),
    top_n_taxa = 15L,
    ignore_unmapped = TRUE,
    ignore_unclassified = TRUE
  ),
  warning = function(warning_condition) {
    taxonomy_warnings <<- c(taxonomy_warnings, conditionMessage(warning_condition))
    invokeRestart("muffleWarning")
  }
)
direct_plot <- SQMtools::plotTaxonomy(
  SQM = SQMtools::subsetSamples(sqm, samples = colnames(sqm$orfs$tpm)),
  rank = "phylum",
  count = "percent",
  N = 15L,
  ignore_unmapped = TRUE,
  ignore_unclassified = TRUE,
  no_partial_classifications = FALSE
)
stopifnot(identical(global_plot$data, direct_plot$data))
palette_size <- length(unique(global_plot$data$item))
stopifnot(identical(
  global_plot$scales$scales[[1L]]$palette(palette_size),
  direct_plot$scales$scales[[1L]]$palette(palette_size)
))
stopifnot(!any(grepl("colors", taxonomy_warnings, ignore.case = TRUE)))

global_data <- script_env$extract_taxonomy_plot_data(global_plot, "percent")
global_data <- script_env$add_global_taxonomy_percent_metadata(
  global_data,
  sqm_object = sqm,
  rank = "phylum",
  selected_samples = colnames(sqm$orfs$tpm)
)
s13_sum <- unique(global_data$displayed_percent_sum[global_data$sample == "S13_1_8"])
stopifnot(length(s13_sum) == 1L, abs(s13_sum - 21.18628) <= 1e-5)
stopifnot(all(abs(
  global_data$displayed_percent_sum + global_data$excluded_percent -
    global_data$raw_percent_sum
) <= 1e-6))
stopifnot(all(abs(
  global_data$accounted_percent_sum - global_data$raw_percent_sum
) <= 1e-6))
stopifnot(all(global_data$effective_excluded_categories == "Unmapped;Unclassified"))
stopifnot(all(is.na(global_data$retained_requested_categories)))

message(
  "PASS: P0 integration checks completed | Bacillota ORFs=",
  length(expected_orf_ids),
  " | pathway samples=", paste(names(percent_sums), collapse = ",")
)
