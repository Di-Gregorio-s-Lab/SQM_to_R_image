source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)
  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

independent_ko_ids <- function(value) {
  if (is.na(value) || !nzchar(trimws(as.character(value)))) {
    return(character())
  }
  unique(stringr::str_extract_all(as.character(value), "K[0-9]{5}")[[1L]])
}

independent_ec_codes <- function(value) {
  if (is.na(value) || !nzchar(trimws(as.character(value)))) {
    return(character())
  }
  ec_block <- stringr::str_match(as.character(value), "\\[EC:([^]]+)\\]")[, 2]
  if (is.na(ec_block)) {
    return(character())
  }
  codes <- strsplit(trimws(ec_block), "[[:space:]]+")[[1L]]
  sort(unique(codes[nzchar(codes)]))
}

independent_ko_ec_lookup <- function(orf_table, included_orf_ids) {
  matched <- match(included_orf_ids, rownames(orf_table))
  stopifnot(!anyNA(matched))
  selected <- orf_table[matched, , drop = FALSE]
  records <- vector("list", nrow(selected))
  for (row_index in seq_len(nrow(selected))) {
    ko_ids <- independent_ko_ids(selected[["KEGG ID"]][[row_index]])
    ec_codes <- independent_ec_codes(selected[["KEGGFUN"]][[row_index]])
    if (length(ko_ids) == 0L) {
      next
    }
    if (length(ec_codes) == 0L) {
      ec_codes <- NA_character_
    }
    records[[row_index]] <- tidyr::expand_grid(
      ko_id = ko_ids,
      ec_code = ec_codes
    )
  }

  dplyr::bind_rows(records) |>
    dplyr::group_by(.data$ko_id) |>
    dplyr::summarise(
      ec_codes = {
        codes <- sort(unique(stats::na.omit(.data$ec_code)))
        if (length(codes) == 0L) NA_character_ else paste(codes, collapse = ";")
      },
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$ko_id)
}

script_env <- source_without_main("sqm_plots.R")
project_dir <- file.path("in", "Au_sip")
if (!dir.exists(project_dir)) {
  stop("P2 integration fixture missing: ", project_dir, call. = FALSE)
}

# Keep the SQMtools/SqueezeMeta compatibility warning visible.
sqm <- SQMtools::loadSQM(
  project_path = script_env$normalize_sqm_project_dir(project_dir),
  tax_mode = "prokfilter",
  trusted_functions_only = FALSE,
  load_sequences = FALSE
)
selected_samples <- colnames(sqm$orfs$tpm)

# Real FLOW case: Unclassified is fourth by TPM, but must not consume one of
# the three classified Top-N slots.
pathway_00361 <- script_env$subset_pathway(
  sqm,
  "Chlorocyclohexane and chlorobenzene degradation"
)
result_00361 <- script_env$build_orf_long_result(pathway_00361, selected_samples)
orf_00361 <- result_00361$data
manual_taxon_rank <- orf_00361 |>
  dplyr::group_by(.data$phylum) |>
  dplyr::summarise(tpm = sum(.data$tpm), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(.data$tpm), .data$phylum)
stopifnot(match("Unclassified", manual_taxon_rank$phylum) == 4L)

flow_rank <- script_env$build_flow_table_for_rank(
  orf_long = orf_00361,
  rank = "phylum",
  selected_samples = selected_samples,
  top_n_taxa = 3L,
  top_n_ko = max(1L, length(unique(orf_00361$ko_id))),
  ko_lookup = script_env$get_ko_name_lookup(pathway_00361)
)
observed_taxa <- unique(as.character(flow_rank$taxon))
stopifnot(all(c("Unclassified", "Other") %in% observed_taxa))
stopifnot(length(setdiff(observed_taxa, c("Unclassified", "Other"))) == 3L)

expected_mass <- orf_00361 |>
  dplyr::group_by(.data$sample) |>
  dplyr::summarise(TPM = sum(.data$tpm), .groups = "drop") |>
  dplyr::arrange(.data$sample)
observed_mass <- flow_rank |>
  dplyr::group_by(.data$sample) |>
  dplyr::summarise(TPM = sum(.data$TPM), .groups = "drop") |>
  dplyr::arrange(.data$sample)
stopifnot(identical(expected_mass$sample, observed_mass$sample))
stopifnot(isTRUE(all.equal(expected_mass$TPM, observed_mass$TPM, tolerance = 1e-10)))
for (sample_name in selected_samples) {
  sample_flow <- script_env$build_flow_table_for_sample(
    flow_rank,
    "Chlorocyclohexane and chlorobenzene degradation",
    "phylum",
    sample_name
  )
  stopifnot(abs(sum(sample_flow$flow_percent) - 100) <= 1e-6)
}

# Real multi-EC case: compare every KO against an independently tokenized
# reference derived from KEGGFUN, then exercise both FUNZ and PIE consumers.
pathway_00710 <- script_env$subset_pathway(
  sqm,
  "Carbon fixation in photosynthetic organisms"
)
result_00710 <- script_env$build_orf_long_result(pathway_00710, selected_samples)
orf_00710 <- result_00710$data
actual_ec <- script_env$extract_ko_ec_lookup(orf_00710) |>
  dplyr::arrange(.data$ko_id)
reference_ec <- independent_ko_ec_lookup(
  as.data.frame(pathway_00710$orfs$table, check.names = FALSE),
  sort(unique(orf_00710$orf_id))
)
stopifnot(identical(actual_ec$ko_id, reference_ec$ko_id))
stopifnot(identical(actual_ec$ec_codes, reference_ec$ec_codes))
multi_ec_rows <- reference_ec |>
  dplyr::filter(!is.na(.data$ec_codes), grepl(";", .data$ec_codes, fixed = TRUE))
stopifnot(nrow(multi_ec_rows) > 0L)

funz_00710 <- script_env$build_ko_plot_table(
  orf_long = orf_00710,
  selected_samples = selected_samples,
  top_n_ko = 999L,
  ko_lookup = script_env$get_ko_name_lookup(pathway_00710),
  pathway_name = "Carbon fixation in photosynthetic organisms"
)
funz_ec <- funz_00710 |>
  dplyr::filter(.data$plotted, !is.na(.data$ko_id), as.character(.data$ko_id) != "Other") |>
  dplyr::transmute(
    ko_id = as.character(.data$ko_id),
    ec_codes = as.character(.data$ec_codes)
  ) |>
  dplyr::distinct() |>
  dplyr::arrange(.data$ko_id)
stopifnot(identical(funz_ec$ko_id, reference_ec$ko_id))
stopifnot(identical(funz_ec$ec_codes, reference_ec$ec_codes))

pie_ko <- multi_ec_rows$ko_id[[1L]]
pie_samples <- orf_00710 |>
  dplyr::filter(.data$ko_id == pie_ko) |>
  dplyr::pull(.data$sample) |>
  unique() |>
  sort()
pie_sample <- pie_samples[[1L]]
pie_table <- script_env$build_pie_chart_table(
  orf_long = orf_00710,
  sample_name = pie_sample,
  ko_id_filter = pie_ko,
  rank_name = "phylum",
  top_n_taxa = 1L
)
expected_pie_ec <- reference_ec$ec_codes[reference_ec$ko_id == pie_ko][[1L]]
stopifnot(nrow(pie_table) > 0L)
stopifnot(all(as.character(pie_table$ec_codes) == expected_pie_ec))

# Full-project audit is calculated without expanding the 945k-ORF table.
full_orf_table <- as.data.frame(sqm$orfs$table, check.names = FALSE)
actual_audit <- script_env$build_ko_expansion_audit(full_orf_table)
independent_counts <- lengths(lapply(full_orf_table[["KEGG ID"]], independent_ko_ids))
stopifnot(actual_audit$input_orf_count[[1L]] == nrow(full_orf_table))
stopifnot(actual_audit$excluded_orfs_without_ko[[1L]] == sum(independent_counts == 0L))
stopifnot(actual_audit$multi_ko_orf_count[[1L]] == sum(independent_counts > 1L))
stopifnot(actual_audit$orf_ko_association_count[[1L]] == sum(independent_counts))
stopifnot(actual_audit$multi_ko_policy[[1L]] == "full_tpm_per_ko")
stopifnot(actual_audit$ko_denominator_basis[[1L]] == "expanded_orf_sample_ko_tpm")

message(
  "PASS: P2 integration checks completed | SQMtools=",
  as.character(utils::packageVersion("SQMtools")),
  " | multi_ec_00710=", nrow(multi_ec_rows),
  " | excluded_orfs_without_ko=", actual_audit$excluded_orfs_without_ko[[1L]],
  " | multi_ko_orfs=", actual_audit$multi_ko_orf_count[[1L]]
)
