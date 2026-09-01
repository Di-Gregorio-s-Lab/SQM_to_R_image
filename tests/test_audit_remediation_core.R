source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

expect_error_contains <- function(code, text, label) {
  observed <- tryCatch(
    {
      force(code)
      NA_character_
    },
    error = function(error) conditionMessage(error)
  )
  if (is.na(observed) || !grepl(text, observed, fixed = TRUE)) {
    stop(label, "; observed error: ", observed, call. = FALSE)
  }
}

script_env <- source_without_main("sqm_plots.R")

pathway_sqm <- list(
  misc = list(
    KEGG_paths = c(
      "Metabolism; Energy metabolism; Nitrogen metabolism",
      "Brite Hierarchies; Protein families; Transporters",
      "Metabolism; Synthetic category A; Ambiguous pathway",
      "Environmental Information Processing; Synthetic category B; Ambiguous pathway"
    )
  )
)

resolved <- script_env$resolve_pathways(pathway_sqm, "Nitrogen metabolism")
stopifnot(identical(resolved[[1L]]$canonical_pathway_name, "Nitrogen metabolism"))
stopifnot(identical(resolved[[1L]]$pathway_id, "00910"))

expect_error_contains(
  script_env$resolve_pathways(pathway_sqm, "Transporters"),
  "not a KEGG PATHWAY",
  "A BRITE leaf was accepted by defined pathway selection"
)
expect_error_contains(
  script_env$resolve_pathways(pathway_sqm, "Ambiguous pathway"),
  "Ambiguous pathway",
  "A pathway leaf present in multiple hierarchies was not rejected"
)

plot_data <- tibble::tibble(
  sample = c("S13_1_8", "S13_1_8"),
  taxon = c("Bacillota", "Other"),
  value = c(20, 1.18628),
  count = "percent"
)
taxonomy_sqm <- list(
  total_reads = c(S13_1_8 = 1000),
  taxa = list(
    phylum = list(
      percent = matrix(
        c(60, 18.81372, 20, 1.18628),
        ncol = 1L,
        dimnames = list(
          c("Unmapped", "Unclassified", "Bacillota", "Rare taxon"),
          "S13_1_8"
        )
      )
    )
  )
)

annotated <- script_env$add_global_taxonomy_percent_metadata(
  plot_data = plot_data,
  sqm_object = taxonomy_sqm,
  rank = "phylum",
  selected_samples = "S13_1_8",
  excluded_categories = c("Unmapped", "Unclassified")
)
stopifnot(all(annotated$percent_of_total_library == annotated$value))
stopifnot(all(annotated$denominator_type == "total_reads"))
stopifnot(all(annotated$denominator_value == 1000))
stopifnot(all(abs(annotated$displayed_percent_sum - 21.18628) < 1e-8))
stopifnot(all(abs(annotated$excluded_percent - 78.81372) < 1e-8))
stopifnot(all(annotated$excluded_categories == "Unmapped;Unclassified"))
stopifnot(all(abs(
  annotated$displayed_percent_sum + annotated$excluded_percent - 100
) < 1e-8))

manifest_row <- script_env$new_manifest_row(
  script_name = "sqm_plots.R",
  project_dir = "project",
  tax_mode = "prokfilter",
  pathway = NA_character_,
  samples = "S13_1_8",
  metric = "percent",
  top_n_taxa = 3L,
  top_n_ko = NA_integer_,
  output_type = "data_tsv",
  output_file = "taxonomy.tsv",
  mode = "taxon",
  ko_selection_policy = "not_applicable",
  taxonomy_display_policy = "SQMtools_non_rescaled_excluding_unmapped_unclassified",
  denominator_type = "total_reads",
  source_data_file = NA_character_,
  sample_order_basis = NA_character_
)
expected_manifest_columns <- c(
  "ko_selection_policy", "taxonomy_display_policy", "denominator_type",
  "source_data_file", "sample_order_basis"
)
stopifnot(all(expected_manifest_columns %in% colnames(manifest_row)))

expect_error_contains(
  script_env$parse_named_args(c("--project_dir", "--output_dir=out")),
  "Missing value for argument: --project_dir",
  "CLI consumed another option as an argument value"
)
expect_error_contains(
  script_env$validate_tax_mode("invalid"),
  "tax_mode",
  "Invalid tax_mode passed validation"
)
expect_error_contains(
  script_env$validate_tpm_matrix(data.frame(S = "broken"), "synthetic TPM"),
  "numeric",
  "Non-numeric TPM passed validation"
)
for (invalid_value in list(NA_real_, Inf, -1)) {
  expect_error_contains(
    script_env$validate_tpm_matrix(data.frame(S = invalid_value), "synthetic TPM"),
    "finite and non-negative",
    "Invalid numeric TPM passed validation"
  )
}
stopifnot(isTRUE(script_env$validate_tpm_matrix(data.frame(S = c(0, 1)), "valid TPM")))

expect_error_contains(
  script_env$parse_dimensions(list(dimensions = "Infx9")),
  "Invalid dimension",
  "Infinite plot dimension passed validation"
)

safe_a <- script_env$safe_output_component("A/B")
safe_b <- script_env$safe_output_component("A B")
stopifnot(!grepl("[/\\\\]", safe_a))
stopifnot(!identical(safe_a, safe_b))
stopifnot(identical(script_env$safe_output_component("S13_1_8"), "S13_1_8"))
compact_a <- script_env$safe_output_component(
  "Chlorocyclohexane and chlorobenzene degradation",
  max_length = 28L
)
compact_b <- script_env$safe_output_component(
  "Chlorocyclohexane and chlorobenzene degradation variant",
  max_length = 28L
)
stopifnot(nchar(compact_a) <= 28L, nchar(compact_b) <= 28L)
stopifnot(!identical(compact_a, compact_b))

long_pathway_name <- "Chlorocyclohexane and chlorobenzene degradation"
stopifnot(identical(
  script_env$taxonomy_file_stem(
    "taxonomy_by_pathway",
    "abund",
    "phylum",
    long_pathway_name
  ),
  "taxonomy_abund_phylum"
))
stopifnot(identical(
  script_env$taxonomy_file_stem(
    "taxonomy_global",
    "percent",
    "species",
    NA_character_
  ),
  "taxonomy_global_percent_species"
))

remediation_root <- file.path(getwd(), "out", "remediation_20260831_01_all")
taxonomy_rank_dir <- file.path(
  remediation_root,
  "taxonomy_by_pathway",
  "definiti",
  compact_a,
  "abund",
  "phylum"
)
pie_ko_dir <- file.path(
  remediation_root,
  "pie",
  "definiti",
  compact_a,
  "S13_1_8",
  script_env$safe_output_component(
    "K00001_EC1.1.1.1;2.2.2.2;3.3.3.3",
    max_length = 20L
  )
)
minimum_compact_png_name <- "p__0123456789ab_9x16.png"
stopifnot(
  nchar(normalizePath(taxonomy_rank_dir, winslash = "/", mustWork = FALSE)) +
    1L + nchar(minimum_compact_png_name) <= 240L,
  nchar(normalizePath(pie_ko_dir, winslash = "/", mustWork = FALSE)) +
    1L + nchar(minimum_compact_png_name) <= 240L
)

message("PASS: core audit remediation invariants are enforced")
