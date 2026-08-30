source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

expect_true <- function(condition, label) {
  if (!isTRUE(condition)) {
    stop(label, call. = FALSE)
  }
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

expect_named_numbers <- function(actual, expected, label, tolerance = 1e-10) {
  expected_names <- sort(names(expected))
  if (!identical(sort(names(actual)), expected_names)) {
    stop(
      label,
      "; expected categories ", paste(expected_names, collapse = ", "),
      ", observed ", paste(sort(names(actual)), collapse = ", "),
      call. = FALSE
    )
  }
  observed <- unname(actual[expected_names])
  wanted <- unname(expected[expected_names])
  if (!isTRUE(all.equal(observed, wanted, tolerance = tolerance))) {
    stop(
      label,
      "; expected values ", paste(wanted, collapse = ", "),
      ", observed ", paste(observed, collapse = ", "),
      call. = FALSE
    )
  }
}

capture_error <- function(code) {
  tryCatch(
    {
      force(code)
      NULL
    },
    error = function(error) error
  )
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

script_env <- source_without_main("sqm_plots.R")

# Unclassified deliberately has the second-largest TPM. It is informative and
# must remain visible, but it must not consume either of the two classified
# Top-N positions. Gamma is the only classified taxon that belongs in Other.
reserved_taxa_orfs <- tibble::tibble(
  orf_id = c("orf_alpha", "orf_unclassified", "orf_beta", "orf_gamma"),
  sample = rep("S_reserved", 4L),
  tpm = c(50, 40, 30, 20),
  ko_id = rep("K00001", 4L),
  kegg_function = rep("Synthetic KO", 4L),
  phylum = c("Alpha", "Unclassified", "Beta", "Gamma")
)

expected_tpm <- c(
  Alpha = 50,
  Beta = 30,
  Other = 20,
  Unclassified = 40
)

run_case("FLOW keeps Unclassified outside the classified Top N", {
  flow_rank <- script_env$build_flow_table_for_rank(
    orf_long = reserved_taxa_orfs,
    rank = "phylum",
    selected_samples = "S_reserved",
    top_n_taxa = 2L,
    top_n_ko = 1L,
    ko_lookup = c(K00001 = "Synthetic KO")
  )

  flow_tpm <- stats::setNames(flow_rank$TPM, as.character(flow_rank$taxon))
  expect_named_numbers(
    flow_tpm,
    expected_tpm,
    "FLOW did not preserve classified Top N, Unclassified, and Other"
  )
  classified_taxa <- sort(setdiff(names(flow_tpm), c("Unclassified", "Other")))
  expect_identical(
    classified_taxa,
    c("Alpha", "Beta"),
    "FLOW Top N must count classified taxa only"
  )
  expect_number(sum(flow_rank$TPM), 140, "FLOW taxonomy collapse changed TPM mass")

  sample_table <- script_env$build_flow_table_for_sample(
    flow_rank_table = flow_rank,
    pathway_name = "Synthetic pathway",
    rank = "phylum",
    sample_name = "S_reserved"
  )
  expect_number(
    sum(sample_table$flow_percent),
    100,
    "FLOW taxonomy collapse changed sample percentages",
    tolerance = 1e-6
  )
})

run_case("PIE keeps Unclassified outside the classified Top N", {
  pie_table <- script_env$build_pie_chart_table(
    orf_long = reserved_taxa_orfs,
    sample_name = "S_reserved",
    ko_id_filter = "K00001",
    rank_name = "phylum",
    top_n_taxa = 2L
  )

  pie_tpm <- stats::setNames(
    pie_table$tpm,
    as.character(pie_table$taxon_rank)
  )
  expect_named_numbers(
    pie_tpm,
    expected_tpm,
    "PIE did not preserve classified Top N, Unclassified, and Other"
  )
  classified_taxa <- sort(setdiff(names(pie_tpm), c("Unclassified", "Other")))
  expect_identical(
    classified_taxa,
    c("Alpha", "Beta"),
    "PIE Top N must count classified taxa only"
  )
  expect_number(sum(pie_table$tpm), 140, "PIE taxonomy collapse changed TPM mass")
  expect_number(
    sum(pie_table$pct),
    1,
    "PIE taxonomy collapse changed sample percentages",
    tolerance = 1e-10
  )
})

source_other_orfs <- tibble::tibble(
  orf_id = c("orf_alpha", "orf_source_other"),
  sample = rep("S_reserved", 2L),
  tpm = c(10, 5),
  ko_id = rep("K00001", 2L),
  kegg_function = rep("Synthetic KO", 2L),
  phylum = c("Alpha", "Other")
)

run_case("reserved source label Other is rejected by FLOW and PIE", {
  flow_error <- capture_error(script_env$build_flow_table_for_rank(
    orf_long = source_other_orfs,
    rank = "phylum",
    selected_samples = "S_reserved",
    top_n_taxa = 1L,
    top_n_ko = 1L,
    ko_lookup = c(K00001 = "Synthetic KO")
  ))
  pie_error <- capture_error(script_env$build_pie_chart_table(
    orf_long = source_other_orfs,
    sample_name = "S_reserved",
    ko_id_filter = "K00001",
    rank_name = "phylum",
    top_n_taxa = 1L
  ))

  expect_true(
    inherits(flow_error, "error"),
    "FLOW must reject the reserved source taxonomy label Other"
  )
  expect_true(
    inherits(pie_error, "error"),
    "PIE must reject the reserved source taxonomy label Other"
  )
  expect_true(
    grepl("reserved.*Other|Other.*reserved", conditionMessage(flow_error), ignore.case = TRUE),
    "FLOW error must identify Other as a reserved source label"
  )
  expect_true(
    grepl("reserved.*Other|Other.*reserved", conditionMessage(pie_error), ignore.case = TRUE),
    "PIE error must identify Other as a reserved source label"
  )
})

if (length(failures) > 0L) {
  stop(
    paste(
      "P2 reserved taxonomy regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: P2 reserved taxonomy categories remain biologically distinct")
