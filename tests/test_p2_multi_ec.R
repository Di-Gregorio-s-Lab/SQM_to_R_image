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

script_env <- source_without_main("sqm_plots.R")

run_case("EC extraction retains valid incomplete codes", {
  observed <- script_env$extract_ec_codes(c(
    "Synthetic enzyme [EC:2.2.2.2 1.1.1.- 1.1.1.1]",
    "Synthetic enzyme without an EC block"
  ))
  expect_identical(
    observed,
    c("2.2.2.2;1.1.1.-;1.1.1.1", NA_character_),
    "EC extraction lost an incomplete code or invented an absent EC"
  )
})

run_case("KO EC lookup preserves every distinct code deterministically", {
  orf_long <- tibble::tibble(
    orf_id = c("orf_b", "orf_a", "orf_duplicate", "orf_missing", "orf_single"),
    sample = rep("S_multi_ec", 5L),
    tpm = c(7, 11, 3, 5, 13),
    ko_id = c("K00001", "K00001", "K00001", "K00001", "K00003"),
    ec_codes = c(
      "2.2.2.2;1.1.1.-",
      "1.1.1.1;3.3.3.3",
      "1.1.1.1",
      NA_character_,
      "4.4.4.4"
    )
  )

  lookup <- script_env$extract_ko_ec_lookup(orf_long) |>
    dplyr::arrange(.data$ko_id)

  expect_identical(
    colnames(lookup),
    c("ko_id", "ec_codes"),
    "The KO EC lookup schema changed"
  )
  expect_identical(nrow(lookup), 2L, "The KO EC lookup is not one row per KO")
  expect_identical(
    as.character(lookup$ko_id),
    c("K00001", "K00003"),
    "The KO EC lookup changed its KO keys"
  )
  expect_identical(
    as.character(lookup$ec_codes),
    c("1.1.1.-;1.1.1.1;2.2.2.2;3.3.3.3", "4.4.4.4"),
    "Distinct EC tokens were truncated, duplicated, or not sorted"
  )
})

run_case("KO without EC metadata remains one row with NA", {
  no_ec_orfs <- tibble::tibble(
    ko_id = c("K00002", "K00002", "K00002"),
    ec_codes = c(NA_character_, "", NA_character_)
  )
  lookup <- script_env$extract_ko_ec_lookup(no_ec_orfs)

  expect_identical(nrow(lookup), 1L, "A KO without EC metadata was dropped or duplicated")
  expect_identical(as.character(lookup$ko_id), "K00002", "The no-EC KO key changed")
  expect_true(
    length(lookup$ec_codes) == 1L && is.na(lookup$ec_codes[[1L]]),
    "A KO without EC metadata must use NA, not an empty string"
  )
})

run_case("EC metadata attachment cannot duplicate abundance", {
  abundance <- tibble::tibble(
    sample = c("S_join", "S_join", "S_join"),
    taxon = c("Alpha", "Beta", "Gamma"),
    ko_id = c("K00001", "K00001", "K00003"),
    tpm = c(10, 20, 30)
  )
  lookup <- script_env$extract_ko_ec_lookup(tibble::tibble(
    ko_id = c("K00001", "K00001", "K00003"),
    ec_codes = c("2.2.2.2;1.1.1.-", "1.1.1.1", "4.4.4.4")
  ))
  joined <- dplyr::left_join(abundance, lookup, by = "ko_id")

  expect_identical(nrow(joined), 3L, "EC metadata duplicated abundance rows")
  expect_identical(
    anyDuplicated(joined[c("sample", "taxon", "ko_id")]),
    0L,
    "EC metadata duplicated an abundance key"
  )
  expect_number(sum(joined$tpm), 60, "EC metadata changed total TPM")
})

run_case("FUNZ table retains complete EC metadata without changing TPM", {
  orf_long <- tibble::tibble(
    orf_id = c("orf_k1_b", "orf_k1_a", "orf_k2"),
    sample = rep("S_funz", 3L),
    tpm = c(10, 20, 30),
    ko_id = c("K00001", "K00001", "K00002"),
    kegg_function = c("Enzyme one", "Enzyme one", "Enzyme two"),
    ec_codes = c("2.2.2.2;1.1.1.-", "1.1.1.1", "4.4.4.4")
  )

  plot_tbl <- script_env$build_ko_plot_table(
    orf_long = orf_long,
    selected_samples = "S_funz",
    top_n_ko = 2L,
    ko_lookup = c(K00001 = "Enzyme one", K00002 = "Enzyme two")
  ) |>
    dplyr::arrange(.data$ko_id)

  expect_identical(nrow(plot_tbl), 2L, "FUNZ metadata changed the KO row count")
  expect_number(sum(plot_tbl$tpm), 60, "FUNZ metadata changed total TPM")
  expect_number(
    sum(plot_tbl$sample_pathway_percent),
    100,
    "FUNZ metadata changed the sample percentage total",
    tolerance = 1e-6
  )
  observed_ec <- stats::setNames(as.character(plot_tbl$ec_codes), plot_tbl$ko_id)
  expect_identical(
    unname(observed_ec[c("K00001", "K00002")]),
    c("1.1.1.-;1.1.1.1;2.2.2.2", "4.4.4.4"),
    "FUNZ output does not expose the complete deterministic EC lookup"
  )
})

run_case("PIE table exposes complete EC metadata used by its output", {
  orf_long <- tibble::tibble(
    orf_id = c("orf_alpha", "orf_beta", "orf_other_ko"),
    sample = rep("S_pie", 3L),
    tpm = c(10, 20, 30),
    ko_id = c("K00001", "K00001", "K00002"),
    kegg_function = c("Enzyme one", "Enzyme one", "Enzyme two"),
    ec_codes = c("2.2.2.2;1.1.1.-", "1.1.1.1", "4.4.4.4"),
    phylum = c("Alpha", "Beta", "Gamma")
  )

  pie_tbl <- script_env$build_pie_chart_table(
    orf_long = orf_long,
    sample_name = "S_pie",
    ko_id_filter = "K00001",
    rank_name = "phylum",
    top_n_taxa = 2L
  )

  expect_true(
    all(c("sample", "ko_id", "ec_codes") %in% colnames(pie_tbl)),
    "PIE output is missing sample, KO, or complete EC metadata required by its TSV"
  )
  expect_identical(
    unique(as.character(pie_tbl$sample)),
    "S_pie",
    "PIE output does not identify its sample"
  )
  expect_identical(
    unique(as.character(pie_tbl$ko_id)),
    "K00001",
    "PIE output does not identify its KO"
  )
  expect_identical(
    unique(as.character(pie_tbl$ec_codes)),
    "1.1.1.-;1.1.1.1;2.2.2.2",
    "PIE output does not expose the complete deterministic EC lookup"
  )
  expect_number(sum(pie_tbl$tpm), 30, "PIE metadata changed total TPM")
})

run_case("TSV round-trip preserves the complete KO EC lookup", {
  orf_long <- tibble::tibble(
    ko_id = c("K00001", "K00001", "K00002"),
    ec_codes = c("2.2.2.2;1.1.1.-", "1.1.1.1", NA_character_)
  )
  lookup <- script_env$extract_ko_ec_lookup(orf_long) |>
    dplyr::arrange(.data$ko_id)

  temp_root <- tempfile("p2_multi_ec_")
  dir.create(temp_root, recursive = TRUE)
  on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)
  tsv_path <- file.path(temp_root, "ko_ec_lookup.tsv")
  script_env$write_tsv_safe(lookup, tsv_path)
  observed <- readr::read_tsv(
    tsv_path,
    na = "NA",
    col_types = readr::cols(.default = readr::col_character()),
    progress = FALSE
  )

  expect_identical(
    as.character(observed$ko_id),
    c("K00001", "K00002"),
    "TSV KO keys differ from the lookup"
  )
  expect_identical(
    as.character(observed$ec_codes),
    c("1.1.1.-;1.1.1.1;2.2.2.2", NA_character_),
    "TSV output truncated or changed the complete EC metadata"
  )
})

if (length(failures) > 0L) {
  stop(
    paste(
      "P2 multi-EC regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: P2 KO metadata preserves complete deterministic EC associations")
