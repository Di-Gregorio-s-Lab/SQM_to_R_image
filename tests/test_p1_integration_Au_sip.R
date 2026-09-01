source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

manual_top_pathways <- function(sqm_object, selected_samples, top_n, allowed_roots) {
  orf_table <- as.data.frame(sqm_object$orfs$table, check.names = FALSE) |>
    tibble::rownames_to_column("orf_id")
  tpm_table <- as.data.frame(sqm_object$orfs$tpm, check.names = FALSE) |>
    tibble::rownames_to_column("orf_id")

  membership <- orf_table |>
    dplyr::transmute(
      orf_id = .data$orf_id,
      hierarchy = as.character(.data$KEGGPATH)
    ) |>
    dplyr::filter(!is.na(.data$hierarchy) & nzchar(trimws(.data$hierarchy))) |>
    dplyr::mutate(entry = stringr::str_split(.data$hierarchy, "\\s*\\|\\s*")) |>
    tidyr::unnest_longer("entry") |>
    dplyr::mutate(parts = stringr::str_split(.data$entry, "\\s*;\\s*")) |>
    dplyr::filter(lengths(.data$parts) == 3L) |>
    dplyr::transmute(
      orf_id = .data$orf_id,
      pathway_root = purrr::map_chr(.data$parts, 1L),
      pathway_category = purrr::map_chr(.data$parts, 2L),
      canonical_pathway_name = purrr::map_chr(.data$parts, 3L)
    ) |>
    dplyr::filter(.data$pathway_root %in% allowed_roots) |>
    dplyr::distinct(
      .data$orf_id,
      .data$pathway_root,
      .data$pathway_category,
      .data$canonical_pathway_name
    )

  hierarchy_by_name <- membership |>
    dplyr::distinct(
      .data$canonical_pathway_name,
      .data$pathway_root,
      .data$pathway_category
    ) |>
    dplyr::count(.data$canonical_pathway_name, name = "hierarchy_count")
  stopifnot(all(hierarchy_by_name$hierarchy_count == 1L))

  tpm_table |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(selected_samples),
      names_to = "sample",
      values_to = "tpm"
    ) |>
    dplyr::inner_join(membership, by = "orf_id", relationship = "many-to-many") |>
    dplyr::group_by(
      .data$pathway_root,
      .data$pathway_category,
      .data$canonical_pathway_name
    ) |>
    dplyr::summarise(total_tpm = sum(as.numeric(.data$tpm), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(
      dplyr::desc(.data$total_tpm),
      .data$canonical_pathway_name,
      .data$pathway_root,
      .data$pathway_category
    ) |>
    dplyr::slice_head(n = top_n)
}

selected_pathway_table <- function(pathways) {
  dplyr::bind_rows(lapply(pathways, function(pathway) {
    tibble::tibble(
      pathway_root = pathway$pathway_root,
      pathway_category = pathway$pathway_category,
      canonical_pathway_name = pathway$canonical_pathway_name,
      total_tpm = pathway$total_tpm
    )
  }))
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

allowed_roots <- c(
  "Metabolism",
  "Genetic Information Processing",
  "Environmental Information Processing",
  "Cellular Processes",
  "Organismal Systems",
  "Human Diseases"
)
selected_samples <- colnames(sqm$orfs$tpm)

global_expected <- manual_top_pathways(sqm, selected_samples, 20L, allowed_roots)
global_actual <- selected_pathway_table(
  script_env$select_top_pathways(sqm, selected_samples, 20L)
)
stopifnot(identical(global_actual$canonical_pathway_name, global_expected$canonical_pathway_name))
stopifnot(identical(global_actual$pathway_root, global_expected$pathway_root))
stopifnot(isTRUE(all.equal(global_actual$total_tpm, global_expected$total_tpm, tolerance = 1e-10)))

known_non_pathways <- c(
  "Transporters",
  "DNA repair and recombination proteins",
  "Enzymes with EC numbers",
  "Function unknown",
  "Peptidases",
  "Transcription factors"
)
stopifnot(all(global_actual$pathway_root %in% allowed_roots))
stopifnot(!any(global_actual$canonical_pathway_name %in% known_non_pathways))
transporters_error <- tryCatch(
  {
    script_env$resolve_pathways(sqm, "Transporters")
    NA_character_
  },
  error = function(error) conditionMessage(error)
)
stopifnot(
  !is.na(transporters_error),
  grepl("not a KEGG PATHWAY", transporters_error, fixed = TRUE)
)

bacillota <- script_env$resolve_taxa_filters(sqm, "Bacillota")[[1]]
bacillota_sqm <- script_env$subset_sqm_by_taxon(sqm, bacillota$orf_ids)
bacillota_expected <- manual_top_pathways(bacillota_sqm, selected_samples, 20L, allowed_roots)
bacillota_actual <- selected_pathway_table(
  script_env$select_top_pathways(bacillota_sqm, selected_samples, 20L)
)
stopifnot(identical(bacillota_actual$canonical_pathway_name, bacillota_expected$canonical_pathway_name))
stopifnot(identical(bacillota_actual$pathway_root, bacillota_expected$pathway_root))
stopifnot(isTRUE(all.equal(bacillota_actual$total_tpm, bacillota_expected$total_tpm, tolerance = 1e-10)))
stopifnot(!identical(
  global_actual$canonical_pathway_name,
  bacillota_actual$canonical_pathway_name
))
stopifnot(identical(
  utils::head(bacillota_actual$canonical_pathway_name, 3L),
  c("Quorum sensing", "ABC transporters", "Two-component system")
))

pathway_name <- "Chlorocyclohexane and chlorobenzene degradation"
pathway_sqm <- script_env$subset_pathway(sqm, pathway_name)
orf_long <- script_env$build_orf_long_table(pathway_sqm, selected_samples)
flow_rank <- script_env$build_flow_table_for_rank(
  orf_long = orf_long,
  rank = "phylum",
  selected_samples = selected_samples,
  top_n_taxa = length(unique(orf_long$phylum)),
  top_n_ko = length(unique(orf_long$ko_id)),
  ko_lookup = script_env$get_ko_name_lookup(pathway_sqm)
)

expected_mass <- orf_long |>
  dplyr::group_by(.data$sample) |>
  dplyr::summarise(TPM = sum(.data$tpm), .groups = "drop") |>
  dplyr::arrange(.data$sample)
actual_mass <- flow_rank |>
  dplyr::group_by(.data$sample) |>
  dplyr::summarise(TPM = sum(.data$TPM), .groups = "drop") |>
  dplyr::arrange(.data$sample)
stopifnot(identical(actual_mass$sample, expected_mass$sample))
stopifnot(isTRUE(all.equal(actual_mass$TPM, expected_mass$TPM, tolerance = 1e-10)))
stopifnot(!anyDuplicated(flow_rank[c("sample", "taxon", "KO")]))

for (sample_name in selected_samples) {
  sample_flow <- script_env$build_flow_table_for_sample(
    flow_rank,
    pathway_name,
    "phylum",
    sample_name
  )
  stopifnot(abs(sum(sample_flow$flow_percent) - 100) <= 1e-6)
}

message(
  "PASS: P1 integration checks completed | SQMtools=",
  as.character(utils::packageVersion("SQMtools")),
  " | global_top20=", nrow(global_actual),
  " | Bacillota_top20=", nrow(bacillota_actual)
)
