source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

expect_identical <- function(actual, expected, label) {
  if (!identical(actual, expected)) {
    stop(
      label,
      "; expected ", paste(expected, collapse = ", "),
      ", observed ", paste(actual, collapse = ", "),
      call. = FALSE
    )
  }
}

expect_true <- function(condition, label) {
  if (!isTRUE(condition)) {
    stop(label, call. = FALSE)
  }
}

expect_number <- function(actual, expected, label, tolerance = 1e-10) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance))) {
    stop(
      label,
      "; expected ", expected,
      ", observed ", actual,
      call. = FALSE
    )
  }
}

failures <- character()
run_case <- function(name, code) {
  tryCatch(
    {
      force(code)
      message("PASS: ", name)
    },
    error = function(error) {
      failures <<- c(failures, paste0(name, ": ", conditionMessage(error)))
      message("FAIL: ", name, " -- ", conditionMessage(error))
    }
  )
}

require_script_function <- function(script_env, name) {
  if (!exists(name, envir = script_env, mode = "function", inherits = FALSE)) {
    stop("Required P2 provenance helper is absent: ", name, call. = FALSE)
  }
  get(name, envir = script_env, mode = "function", inherits = FALSE)
}

script_env <- source_without_main("sqm_plots.R")

orf_ids <- c("orf_without_ko", "orf_mono_ko", "orf_multi_ko")
fake_pathway_sqm <- list(
  orfs = list(
    table = data.frame(
      "KEGG ID" = c(NA_character_, "K00001", "K00002; K00003*"),
      KEGGFUN = c(
        "No KO annotation",
        "Mono-KO enzyme [EC:1.1.1.1]",
        "Multi-KO enzyme [EC:2.2.2.2 3.3.3.-]"
      ),
      KEGGPATH = rep("Synthetic pathway", 3L),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tax = data.frame(
      superkingdom = rep("Bacteria", 3L),
      phylum = c("Alpha", "Beta", "Gamma"),
      class = rep("Synthetic class", 3L),
      order = rep("Synthetic order", 3L),
      family = rep("Synthetic family", 3L),
      genus = rep("Synthetic genus", 3L),
      species = rep("Synthetic species", 3L),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tpm = data.frame(
      S_positive = c(4, 10, 20),
      S_second = c(6, 0, 30),
      row.names = orf_ids,
      check.names = FALSE
    )
  ),
  functions = list(
    KEGG = list(
      tpm = data.frame(
        S_positive = c(10, 10, 10),
        S_second = c(0, 15, 15),
        row.names = c("K00001", "K00002", "K00003"),
        check.names = FALSE
      )
    )
  ),
  misc = list(
    KEGG_names = c(
      K00001 = "Mono-KO enzyme",
      K00002 = "Multi-KO enzyme A",
      K00003 = "Multi-KO enzyme B"
    )
  )
)

selected_samples <- c("S_positive", "S_second")
required_audit_fields <- c(
  "input_orf_count",
  "excluded_orfs_without_ko",
  "multi_ko_orf_count",
  "orf_ko_association_count",
  "multi_ko_policy",
  "ko_denominator_basis"
)
expected_audit <- tibble::tibble(
  input_orf_count = 3L,
  excluded_orfs_without_ko = 1L,
  multi_ko_orf_count = 1L,
  orf_ko_association_count = 3L,
  multi_ko_policy = "full_tpm_per_ko",
  ko_denominator_basis = "expanded_orf_sample_ko_tpm"
)

expect_manifest_audit <- function(manifest_tbl, label) {
  expect_true(nrow(manifest_tbl) > 0L, paste0(label, " produced no manifest rows"))
  expect_true(
    all(required_audit_fields %in% colnames(manifest_tbl)),
    paste0(label, " manifest is missing KO provenance fields")
  )

  expected_counts <- c(
    input_orf_count = 3L,
    excluded_orfs_without_ko = 1L,
    multi_ko_orf_count = 1L,
    orf_ko_association_count = 3L
  )
  for (field in names(expected_counts)) {
    expect_identical(
      as.integer(manifest_tbl[[field]]),
      rep(unname(expected_counts[[field]]), nrow(manifest_tbl)),
      paste0(label, " manifest did not propagate ", field)
    )
  }
  expect_identical(
    as.character(manifest_tbl$multi_ko_policy),
    rep("full_tpm_per_ko", nrow(manifest_tbl)),
    paste0(label, " manifest did not propagate multi_ko_policy")
  )
  expect_identical(
    as.character(manifest_tbl$ko_denominator_basis),
    rep("expanded_orf_sample_ko_tpm", nrow(manifest_tbl)),
    paste0(label, " manifest did not propagate ko_denominator_basis")
  )
}

run_case("ORF long result exposes stable KO provenance", {
  build_orf_long_result <- require_script_function(
    script_env,
    "build_orf_long_result"
  )
  result <- build_orf_long_result(fake_pathway_sqm, selected_samples)

  expect_true(is.list(result), "The ORF expansion result must be a list")
  expect_identical(
    names(result),
    c("data", "audit"),
    "The ORF expansion result must expose data and audit"
  )
  expect_true(
    is.data.frame(result$data),
    "The data member is not a tabular ORF/sample/KO result"
  )
  expect_true(
    tibble::is_tibble(result$audit) && nrow(result$audit) == 1L,
    "The audit member must be a one-row tibble"
  )

  expect_true(
    all(required_audit_fields %in% colnames(result$audit)),
    "The KO provenance audit is missing required fields"
  )

  observed_counts <- as.integer(unlist(result$audit[1L, c(
    "input_orf_count",
    "excluded_orfs_without_ko",
    "multi_ko_orf_count",
    "orf_ko_association_count"
  )], use.names = FALSE))
  expect_identical(
    observed_counts,
    c(3L, 1L, 1L, 3L),
    "KO provenance counts are not based on unique input ORFs and ORF/KO associations"
  )
  expect_identical(
    as.character(result$audit$multi_ko_policy),
    "full_tpm_per_ko",
    "The multi-KO policy is absent or incorrect"
  )
  expect_identical(
    as.character(result$audit$ko_denominator_basis),
    "expanded_orf_sample_ko_tpm",
    "The KO percentage denominator basis is absent or incorrect"
  )
})

run_case("KO manifest row exposes exact expansion provenance", {
  result <- list(audit = expected_audit)
  manifest_row <- script_env$new_manifest_row(
    script_name = "sqm_plots.R",
    project_dir = "synthetic_project",
    tax_mode = "prokfilter",
    pathway = "Synthetic pathway",
    samples = selected_samples,
    metric = "TPM",
    top_n_taxa = 3L,
    top_n_ko = 3L,
    output_type = "tsv",
    output_file = "funz/synthetic.tsv",
    mode = "funz",
    ko_audit = result$audit
  )

  expect_true(
    all(required_audit_fields %in% colnames(manifest_row)),
    "The KO manifest row is missing provenance fields"
  )
  observed_counts <- as.integer(unlist(manifest_row[1L, c(
    "input_orf_count",
    "excluded_orfs_without_ko",
    "multi_ko_orf_count",
    "orf_ko_association_count"
  )], use.names = FALSE))
  expect_identical(
    observed_counts,
    c(3L, 1L, 1L, 3L),
    "The KO manifest row changed expansion provenance counts"
  )
  expect_identical(
    as.character(manifest_row$multi_ko_policy),
    "full_tpm_per_ko",
    "The KO manifest row changed the multi-KO policy"
  )
  expect_identical(
    as.character(manifest_row$ko_denominator_basis),
    "expanded_orf_sample_ko_tpm",
    "The KO manifest row changed the denominator basis"
  )
})

run_case("non-KO manifest row retains provenance schema with NA", {
  manifest_row <- script_env$new_manifest_row(
    script_name = "sqm_plots.R",
    project_dir = "synthetic_project",
    tax_mode = "prokfilter",
    pathway = "Synthetic pathway",
    samples = selected_samples,
    metric = "TPM",
    top_n_taxa = 3L,
    top_n_ko = 3L,
    output_type = "tsv",
    output_file = "taxonomy/synthetic.tsv",
    mode = "taxon"
  )

  expect_true(
    all(required_audit_fields %in% colnames(manifest_row)),
    "A non-KO manifest row does not retain the append-only provenance schema"
  )
  expect_true(
    all(vapply(manifest_row[required_audit_fields], function(value) {
      length(value) == 1L && is.na(value[[1L]])
    }, logical(1L))),
    "A non-KO manifest row must use NA for every KO provenance field"
  )
})

run_case("KO mode runners propagate expansion provenance to every manifest row", {
  temp_root <- tempfile("p2_ko_manifest_modes_")
  dir.create(temp_root, recursive = TRUE)
  on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)

  common_args <- list(
    output_dir = temp_root,
    manifest_base_dir = temp_root,
    script_name = "sqm_plots.R",
    project_dir = "synthetic_project",
    tax_mode = "prokfilter",
    pathway_name = "Synthetic pathway",
    pathway_sqm = fake_pathway_sqm,
    selected_samples = "S_positive",
    dimensions = list(),
    plot_dpi = 72,
    top_n_taxa = 3L,
    top_n_ko = 3L,
    pathway_id = "00000"
  )

  funz_result <- do.call(
    script_env$run_funz_mode,
    c(common_args, list(output_manifests = list(funz = tibble::tibble())))
  )
  flow_result <- do.call(
    script_env$run_flow_mode,
    c(
      common_args,
      list(
        output_manifests = list(flow = tibble::tibble()),
        taxonomy_ranks = "phylum",
        flowplot_formats = character()
      )
    )
  )
  pie_result <- do.call(
    script_env$run_pie_mode,
    c(
      common_args,
      list(
        output_manifests = list(pie = tibble::tibble()),
        taxonomy_ranks = "phylum"
      )
    )
  )

  expect_manifest_audit(funz_result$funz, "FUNZ")
  expect_manifest_audit(flow_result$flow, "FLOW")
  expect_manifest_audit(pie_result$pie, "PIE")
})

run_case("Pathway analysis is built once and reused by every KO mode", {
  temp_root <- tempfile("shared_pathway_analysis_")
  dir.create(temp_root, recursive = TRUE)
  on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)

  build_calls <- 0L
  original_builder <- script_env$build_orf_long_result
  script_env$build_orf_long_result <- function(pathway_sqm, selected_samples) {
    build_calls <<- build_calls + 1L
    original_builder(pathway_sqm, selected_samples)
  }
  on.exit({
    script_env$build_orf_long_result <- original_builder
  }, add = TRUE)

  pathway_info <- list(
    pathway_name = "Synthetic pathway",
    pathway_id = "00000",
    pathway_selection = "defined",
    pathway_sqm = fake_pathway_sqm
  )
  pathway_analysis <- script_env$build_pathway_analysis(
    pathway_info,
    "S_positive"
  )
  expect_identical(
    build_calls,
    1L,
    "Building one pathway analysis did not perform exactly one ORF expansion"
  )

  common_args <- list(
    output_dir = temp_root,
    manifest_base_dir = temp_root,
    script_name = "sqm_plots.R",
    project_dir = "synthetic_project",
    tax_mode = "prokfilter",
    pathway_name = pathway_info$pathway_name,
    pathway_sqm = pathway_info$pathway_sqm,
    selected_samples = "S_positive",
    dimensions = list(`2x2` = c(width = 2, height = 2)),
    plot_dpi = 75,
    top_n_taxa = 3L,
    top_n_ko = 3L,
    pathway_id = pathway_info$pathway_id,
    pathway_selection = pathway_info$pathway_selection,
    pathway_analysis = pathway_analysis
  )

  funz_result <- do.call(
    script_env$run_funz_mode,
    c(common_args, list(output_manifests = list(funz = tibble::tibble())))
  )
  flow_result <- do.call(
    script_env$run_flow_mode,
    c(
      common_args,
      list(
        output_manifests = list(flow = tibble::tibble()),
        taxonomy_ranks = "phylum",
        flowplot_formats = "png"
      )
    )
  )
  pie_result <- do.call(
    script_env$run_pie_mode,
    c(
      common_args,
      list(
        output_manifests = list(pie = tibble::tibble()),
        taxonomy_ranks = "phylum"
      )
    )
  )

  expect_identical(
    build_calls,
    1L,
    "FUNZ, FLOW, or PIE rebuilt the shared ORF-long result"
  )
  expect_manifest_audit(funz_result$funz, "shared FUNZ")
  expect_manifest_audit(flow_result$flow, "shared FLOW")
  expect_manifest_audit(pie_result$pie, "shared PIE")
  rendered_pngs <- list.files(
    temp_root,
    pattern = "\\.png$",
    recursive = TRUE,
    full.names = TRUE
  )
  expect_true(
    length(rendered_pngs) >= 3L && all(file.info(rendered_pngs)$size > 0L),
    "The shared 75-DPI analysis did not render non-empty PNG artifacts"
  )
  rendered_rows <- dplyr::bind_rows(
    funz_result$funz,
    flow_result$flow,
    pie_result$pie
  ) |>
    dplyr::filter(.data$format == "png")
  expect_true(
    nrow(rendered_rows) == length(rendered_pngs) &&
      all(rendered_rows$dpi == 75),
    "The shared analysis manifests did not preserve 75 DPI for every PNG"
  )
})

run_case("canonical multi-KO allocation splits TPM before pathway filtering", {
  build_pathway_ko_result <- require_script_function(
    script_env,
    "build_pathway_ko_result"
  )
  expanded <- build_pathway_ko_result(
    fake_pathway_sqm,
    selected_samples,
    c("K00001", "K00002", "K00003")
  )$data |>
    dplyr::arrange(.data$orf_id, .data$sample, .data$ko_id)

  expected_keys <- c(
    "orf_mono_ko|S_positive|K00001",
    "orf_multi_ko|S_positive|K00002",
    "orf_multi_ko|S_positive|K00003",
    "orf_multi_ko|S_second|K00002",
    "orf_multi_ko|S_second|K00003"
  )
  observed_keys <- paste(expanded$orf_id, expanded$sample, expanded$ko_id, sep = "|")
  expect_identical(
    observed_keys,
    expected_keys,
    "ORF/sample/KO expansion produced the wrong association keys"
  )
  expect_identical(
    as.numeric(expanded$tpm),
    c(10, 10, 10, 15, 15),
    "A multi-KO ORF was not divided equally among its KO associations"
  )
  expect_number(
    sum(expanded$tpm),
    60,
    "Canonical KO allocations do not conserve the original ORF TPM"
  )
  expect_true(
    !"orf_without_ko" %in% expanded$orf_id,
    "An ORF without KO annotation was retained in KO analysis data"
  )
})

run_case("legacy ORF long wrapper returns the structured data member", {
  build_orf_long_result <- require_script_function(
    script_env,
    "build_orf_long_result"
  )
  result <- build_orf_long_result(fake_pathway_sqm, selected_samples)
  legacy_data <- script_env$build_orf_long_table(fake_pathway_sqm, selected_samples)

  expect_identical(
    colnames(legacy_data),
    colnames(result$data),
    "The compatibility wrapper changed the ORF long schema"
  )
  legacy_keys <- legacy_data |>
    dplyr::arrange(.data$orf_id, .data$sample, .data$ko_id) |>
    dplyr::transmute(key = paste(.data$orf_id, .data$sample, .data$ko_id, sep = "|")) |>
    dplyr::pull(.data$key)
  result_keys <- result$data |>
    dplyr::arrange(.data$orf_id, .data$sample, .data$ko_id) |>
    dplyr::transmute(key = paste(.data$orf_id, .data$sample, .data$ko_id, sep = "|")) |>
    dplyr::pull(.data$key)
  expect_identical(
    legacy_keys,
    result_keys,
    "The compatibility wrapper does not return the structured data member"
  )
  expect_number(
    sum(legacy_data$tpm),
    sum(result$data$tpm),
    "The compatibility wrapper changed expanded TPM"
  )
})

if (length(failures) > 0L) {
  stop(
    paste(
      "P2 KO provenance regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: P2 KO expansion reports provenance and canonical multi-KO splitting")
