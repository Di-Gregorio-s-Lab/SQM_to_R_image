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
    stop("Required FLOW helper is absent: ", name, call. = FALSE)
  }
  get(name, envir = script_env, mode = "function", inherits = FALSE)
}

script_env <- source_without_main("sqm_plots.R")

# Two descriptions belong to one KO. Description A is repeated deliberately:
# metadata must be distinct and deterministic without multiplying abundance.
# Before metadata attachment there are two sample/taxon/KO edges totalling 30.
duplicated_description_orfs <- tibble::tibble(
  orf_id = c("orf_alpha_b", "orf_alpha_a", "orf_beta_a"),
  sample = rep("S_flow", 3L),
  tpm = c(10, 5, 15),
  ko_id = rep("K00001", 3L),
  kegg_function = c("Description B", "Description A", "Description A"),
  ec_codes = c("2.2.2.2;1.1.1.1", "1.1.1.1", "2.2.2.2"),
  phylum = c("Alpha", "Alpha", "Beta")
)

run_case("KO metadata join preserves FLOW edges and TPM", {
  rank_tbl <- script_env$build_flow_table_for_rank(
    orf_long = duplicated_description_orfs,
    rank = "phylum",
    selected_samples = "S_flow",
    top_n_taxa = 2L,
    top_n_ko = 1L,
    ko_lookup = c(K00001 = "Lookup description")
  )

  observed_rows <- nrow(rank_tbl)
  observed_tpm <- sum(rank_tbl$TPM)
  expect_true(
    observed_rows == 2L && isTRUE(all.equal(observed_tpm, 30, tolerance = 1e-10)),
    paste0(
      "Metadata attachment must preserve 2 edges and TPM 30; observed ",
      observed_rows, " edges and TPM ", observed_tpm
    )
  )
  expect_identical(
    anyDuplicated(rank_tbl[c("sample", "taxon", "KO")]),
    0L,
    "Metadata attachment duplicated a sample/taxon/KO key"
  )
  expect_true(
    all(rank_tbl$KO_name == "Description A; Description B"),
    "Distinct KO descriptions were not sorted and concatenated deterministically"
  )
  expect_true(
    all(rank_tbl$ec_codes == "1.1.1.1;2.2.2.2"),
    "FLOW did not retain the complete deterministic KO/EC association"
  )

  sample_tbl <- script_env$build_flow_table_for_sample(
    flow_rank_table = rank_tbl,
    pathway_name = "Synthetic pathway",
    rank = "phylum",
    sample_name = "S_flow"
  )
  expect_identical(nrow(sample_tbl), 2L, "Sample FLOW table changed edge count")
  expect_number(sum(sample_tbl$TPM), 30, "Sample FLOW table changed TPM mass")
  expect_number(
    sum(sample_tbl$flow_percent),
    100,
    "FLOW percentages do not sum to 100 for the sample",
    tolerance = 1e-6
  )
})

run_case("KO metadata uses normalized descriptions and ordered fallbacks", {
  build_flow_ko_metadata <- require_script_function(
    script_env,
    "build_flow_ko_metadata"
  )

  metadata_orfs <- tibble::tibble(
    ko_id = c("K00001", "K00001", "K00001", "K00002", "K00003"),
    kegg_function = c(
      " Description B ", "Description A", "Description A", "   ", NA_character_
    ),
    ec_codes = c(
      "2.2.2.2;1.1.1.1", "1.1.1.1", "2.2.2.2", "", NA_character_
    )
  )
  metadata <- build_flow_ko_metadata(
    orf_long = metadata_orfs,
    ko_lookup = c(K00002 = "Lookup description")
  )

  expect_identical(nrow(metadata), 3L, "KO metadata is not one row per KO")
  expect_identical(
    anyDuplicated(metadata$ko_id),
    0L,
    "KO metadata contains duplicate keys"
  )
  observed_names <- stats::setNames(metadata$KO_name, metadata$ko_id)
  expect_identical(
    unname(observed_names[c("K00001", "K00002", "K00003")]),
    c("Description A; Description B", "Lookup description", "K00003"),
    "KO description normalization or fallback precedence changed"
  )
  observed_ec <- stats::setNames(metadata$ec_codes, metadata$ko_id)
  expect_identical(
    unname(observed_ec[c("K00001", "K00002", "K00003")]),
    c("1.1.1.1;2.2.2.2", NA_character_, NA_character_),
    "FLOW KO/EC metadata did not preserve multiple or missing EC codes"
  )
})

run_case("FLOW labels retain fallback and Other semantics", {
  fallback_orfs <- tibble::tibble(
    orf_id = c("orf_lookup", "orf_id", "orf_other"),
    sample = rep("S_fallback", 3L),
    tpm = c(12, 10, 8),
    ko_id = c("K00002", "K00003", "K00004"),
    kegg_function = c(" ", NA_character_, "Outside top N"),
    ec_codes = c("2.2.2.2", NA_character_, "4.4.4.4"),
    phylum = rep("Alpha", 3L)
  )

  rank_tbl <- script_env$build_flow_table_for_rank(
    orf_long = fallback_orfs,
    rank = "phylum",
    selected_samples = "S_fallback",
    top_n_taxa = 1L,
    top_n_ko = 2L,
    ko_lookup = c(K00002 = "Lookup description")
  )
  observed_names <- stats::setNames(rank_tbl$KO_name, rank_tbl$KO)

  expect_identical(
    unname(observed_names[c("K00002", "K00003", "Other")]),
    c("Lookup description", "K00003", "Other KOs"),
    "FLOW KO fallback or Other label changed"
  )
  observed_ec <- stats::setNames(rank_tbl$ec_codes, rank_tbl$KO)
  expect_true(
    identical(unname(observed_ec[c("K00002", "K00003")]), c("2.2.2.2", NA_character_)) &&
      is.na(unname(observed_ec[["Other"]])),
    "FLOW EC metadata changed missing-EC or Other semantics"
  )

  sample_tbl <- script_env$build_flow_table_for_sample(
    flow_rank_table = rank_tbl,
    pathway_name = "Synthetic pathway",
    rank = "phylum",
    sample_name = "S_fallback"
  )
  expect_number(sum(sample_tbl$TPM), 30, "Fallback labels changed TPM mass")
  expect_number(
    sum(sample_tbl$flow_percent),
    100,
    "FLOW percentages do not sum to 100 after fallback labeling",
    tolerance = 1e-6
  )
})

run_case("non-unique KO metadata fails the join postcondition", {
  join_flow_ko_metadata <- require_script_function(
    script_env,
    "join_flow_ko_metadata"
  )
  summary_tbl <- tibble::tibble(
    sample = "S_flow",
    taxon = "Alpha",
    KO = "K00001",
    TPM = 30
  )
  invalid_metadata <- tibble::tibble(
    ko_id = c("K00001", "K00001"),
    KO_name = c("Description A", "Description B"),
    ec_codes = c("1.1.1.1", "2.2.2.2")
  )

  join_error <- tryCatch(
    {
      join_flow_ko_metadata(summary_tbl, invalid_metadata)
      NULL
    },
    error = function(error) error
  )
  expect_true(
    inherits(join_error, "error"),
    "A many-to-many KO metadata join must fail before returning inflated FLOW data"
  )
  expect_true(
    grepl(
      "metadata.*(unique|one row|one-to-one|many-to-one)|postcondition|duplicate",
      conditionMessage(join_error),
      ignore.case = TRUE
    ),
    "The invalid metadata failure must identify the violated join postcondition"
  )
})

legend_orfs <- tibble::tibble(
  orf_id = c("orf_multi_a", "orf_multi_b", "orf_missing", "orf_other"),
  sample = rep("S_legend", 4L),
  tpm = c(40, 20, 30, 10),
  ko_id = c("K00001", "K00001", "K00002", "K00003"),
  kegg_function = c("Multi EC", "Multi EC", "Missing EC", "Outside Top N"),
  ec_codes = c("2.2.2.2;1.1.1.1", "1.1.1.1", NA_character_, "3.3.3.3"),
  phylum = c("Alpha", "Alpha", "Beta", "Beta")
)

run_case("FLOW legend module exposes taxonomy and KO/EC percentages", {
  build_flow_legend_spec <- require_script_function(
    script_env,
    "build_flow_legend_spec"
  )
  rank_tbl <- script_env$build_flow_table_for_rank(
    orf_long = legend_orfs,
    rank = "phylum",
    selected_samples = "S_legend",
    top_n_taxa = 2L,
    top_n_ko = 2L,
    ko_lookup = character()
  )
  sample_tbl <- script_env$build_flow_table_for_sample(
    rank_tbl,
    "Synthetic pathway",
    "phylum",
    "S_legend"
  )
  expect_true(
    "ec_codes" %in% colnames(sample_tbl),
    "The append-only FLOW TSV schema is missing ec_codes"
  )

  legend_spec <- build_flow_legend_spec(sample_tbl)
  expect_identical(
    names(legend_spec),
    c("taxonomy", "functional"),
    "FLOW legend module changed its two-section interface"
  )
  expect_identical(
    legend_spec$taxonomy$label,
    c("Beta | 40.0%", "Alpha | 60.0%"),
    "Taxonomy legend labels or factor order changed"
  )
  expect_identical(
    legend_spec$functional$label,
    c(
      "K00002 / EC NA | 30.0%",
      "K00001 / EC 1.1.1.1;2.2.2.2 | 60.0%",
      "Other KOs | 10.0%"
    ),
    "Functional KO/EC legend labels or factor order changed"
  )
  expect_number(
    sum(legend_spec$taxonomy$percent),
    100,
    "Taxonomy legend percentages do not sum to 100",
    tolerance = 1e-6
  )
  expect_number(
    sum(legend_spec$functional$percent),
    100,
    "Functional legend percentages do not sum to 100",
    tolerance = 1e-6
  )
})

run_case("FLOW PNG uses the simple legend-free renderer", {
  rank_tbl <- script_env$build_flow_table_for_rank(
    legend_orfs,
    "phylum",
    "S_legend",
    2L,
    2L,
    character()
  )
  sample_tbl <- script_env$build_flow_table_for_sample(
    rank_tbl,
    "Synthetic pathway",
    "phylum",
    "S_legend"
  )
  plot_object <- script_env$make_flow_plot(
    sample_tbl,
    "Synthetic pathway",
    "phylum",
    "S_legend"
  )
  expect_identical(
    plot_object$theme$legend.position,
    "none",
    "FLOW PNG unexpectedly rendered a legend"
  )

  png_path <- tempfile("flow_two_legends_", fileext = ".png")
  on.exit(unlink(png_path, force = TRUE), add = TRUE)
  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = 12,
    height = 9,
    dpi = 75,
    units = "in"
  )
  expect_true(
    file.exists(png_path) && file.info(png_path)$size > 0L,
    "Simple FLOW PNG did not render at 75 DPI"
  )
})

run_case("FLOW Sankey uses plain nodes and detailed hover", {
  rank_tbl <- script_env$build_flow_table_for_rank(
    legend_orfs,
    "phylum",
    "S_legend",
    2L,
    2L,
    character()
  )
  sample_tbl <- script_env$build_flow_table_for_sample(
    rank_tbl,
    "Synthetic pathway",
    "phylum",
    "S_legend"
  )
  widget <- script_env$make_flow_sankey(
    sample_tbl,
    "Synthetic pathway",
    "phylum",
    "S_legend"
  )
  trace <- plotly::plotly_build(widget)$x$data[[1L]]
  node_labels <- as.character(unlist(trace$node$label, use.names = FALSE))
  hover_text <- as.character(unlist(trace$link$customdata, use.names = FALSE))

  expect_true(
    all(c("Alpha", "K00001", "Other KOs") %in% node_labels),
    "FLOW Sankey plain node labels are incomplete"
  )
  expect_true(
    any(grepl("Taxon share: 60.0%", hover_text, fixed = TRUE)) &&
      any(grepl("KO: K00001", hover_text, fixed = TRUE)) &&
      any(grepl("EC: 1.1.1.1;2.2.2.2", hover_text, fixed = TRUE)) &&
      any(grepl("Function share: 60.0%", hover_text, fixed = TRUE)) &&
      any(grepl("TPM:", hover_text, fixed = TRUE)) &&
      any(grepl("Flow:", hover_text, fixed = TRUE)),
    "FLOW Sankey hover does not retain KO/EC, node shares, TPM, and flow percent"
  )
})

if (length(failures) > 0L) {
  stop(
    paste(
      "P1 FLOW KO metadata regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: P1 FLOW KO metadata joins preserve edge keys and TPM mass")
