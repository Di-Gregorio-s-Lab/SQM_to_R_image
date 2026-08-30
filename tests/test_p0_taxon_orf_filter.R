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
    stop(label, call. = FALSE)
  }
}

expect_true <- function(condition, label) {
  if (!isTRUE(condition)) {
    stop(label, call. = FALSE)
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

# The contig taxonomy deliberately disagrees with the ORF taxonomy. Filtering
# Bacillota through contigs would include orf_contig_only and exclude both ORFs
# that the ORF-level resolver selected.
all_orf_ids <- c(
  "orf_bacillota_first",
  "orf_contig_only",
  "orf_bacillota_second",
  "orf_neither"
)
expected_orf_ids <- c("orf_bacillota_first", "orf_bacillota_second")

fake_sqm <- list(
  orfs = list(
    table = data.frame(
      contig_id = c(
        "contig_not_bacillota_1",
        "contig_bacillota",
        "contig_not_bacillota_2",
        "contig_other"
      ),
      row.names = all_orf_ids,
      check.names = FALSE
    ),
    tax = data.frame(
      superkingdom = rep("Bacteria", length(all_orf_ids)),
      phylum = c(
        "Bacillota",
        "Pseudomonadota",
        "Bacillota",
        "Actinomycetota"
      ),
      row.names = all_orf_ids,
      check.names = FALSE
    ),
    tpm = data.frame(
      sample_a = c(10, 20, 30, 40),
      sample_b = c(1, 2, 3, 4),
      row.names = all_orf_ids,
      check.names = FALSE
    )
  ),
  contigs = list(
    tax = data.frame(
      phylum = c(
        "Pseudomonadota",
        "Bacillota",
        "Pseudomonadota",
        "Actinomycetota"
      ),
      row.names = c(
        "contig_not_bacillota_1",
        "contig_bacillota",
        "contig_not_bacillota_2",
        "contig_other"
      ),
      check.names = FALSE
    )
  )
)

run_case("resolver returns exact ordered ORF identifiers", {
  resolved <- script_env$resolve_taxa_filters(fake_sqm, "Bacillota")

  expect_identical(length(resolved), 1L, "Expected one resolved taxon")
  expect_identical(resolved[[1]]$taxon, "Bacillota", "Resolved taxon changed")
  expect_identical(resolved[[1]]$rank, "phylum", "Resolved rank changed")
  expect_identical(
    resolved[[1]]$orf_ids,
    expected_orf_ids,
    "Resolved ORF IDs are missing, reordered, or include contig-level matches"
  )
})

subset_calls <- list()
fake_subset_orfs <- function(
  SQM,
  orfs,
  tax_source = "orfs",
  trusted_functions_only = FALSE,
  ignore_unclassified_functions = FALSE,
  rescale_tpm = FALSE,
  rescale_copy_number = FALSE,
  recalculate_bin_stats = TRUE,
  contigs_override = NULL,
  allow_empty = FALSE
) {
  subset_calls[[length(subset_calls) + 1L]] <<- list(
    orfs = orfs,
    tax_source = tax_source,
    trusted_functions_only = trusted_functions_only,
    ignore_unclassified_functions = ignore_unclassified_functions,
    rescale_tpm = rescale_tpm,
    rescale_copy_number = rescale_copy_number,
    recalculate_bin_stats = recalculate_bin_stats
  )

  result <- SQM
  result$orfs$table <- SQM$orfs$table[orfs, , drop = FALSE]
  result$orfs$tax <- SQM$orfs$tax[orfs, , drop = FALSE]
  result$orfs$tpm <- SQM$orfs$tpm[orfs, , drop = FALSE]
  result
}

run_case("subset uses the ORF boundary without rescaling", {
  subset_calls <- list()
  subset <- script_env$subset_sqm_by_taxon(
    fake_sqm,
    expected_orf_ids,
    subset_orfs_fn = fake_subset_orfs
  )

  expect_identical(length(subset_calls), 1L, "Expected one subsetORFs call")
  call <- subset_calls[[1]]
  expect_identical(call$orfs, expected_orf_ids, "subsetORFs received different ORF IDs")
  expect_identical(call$tax_source, "orfs", "subsetORFs did not use ORF taxonomy")
  expect_identical(
    call$trusted_functions_only,
    FALSE,
    "subsetORFs unexpectedly restricted trusted functions"
  )
  expect_identical(
    call$ignore_unclassified_functions,
    TRUE,
    "subsetORFs must preserve the existing unclassified-function policy"
  )
  expect_identical(call$rescale_tpm, FALSE, "subsetORFs rescaled TPM")
  expect_identical(
    call$rescale_copy_number,
    FALSE,
    "subsetORFs rescaled copy number"
  )
  expect_identical(
    call$recalculate_bin_stats,
    FALSE,
    "subsetORFs recalculated bin statistics"
  )
  expect_identical(
    rownames(subset$orfs$table),
    expected_orf_ids,
    "Returned subset does not contain the exact requested ORF IDs"
  )
})

run_case("subset rejects a backend ORF mismatch", {
  mismatch_backend_called <- FALSE
  mismatching_subset_orfs <- function(
    SQM,
    orfs,
    tax_source = "orfs",
    trusted_functions_only = FALSE,
    ignore_unclassified_functions = FALSE,
    rescale_tpm = FALSE,
    rescale_copy_number = FALSE,
    recalculate_bin_stats = TRUE,
    contigs_override = NULL,
    allow_empty = FALSE
  ) {
    mismatch_backend_called <<- TRUE
    fake_subset_orfs(
      SQM = SQM,
      orfs = orfs[[1]],
      tax_source = tax_source,
      trusted_functions_only = trusted_functions_only,
      ignore_unclassified_functions = ignore_unclassified_functions,
      rescale_tpm = rescale_tpm,
      rescale_copy_number = rescale_copy_number,
      recalculate_bin_stats = recalculate_bin_stats,
      contigs_override = contigs_override,
      allow_empty = allow_empty
    )
  }

  mismatch_error <- tryCatch(
    {
      script_env$subset_sqm_by_taxon(
        fake_sqm,
        expected_orf_ids,
        subset_orfs_fn = mismatching_subset_orfs
      )
      NULL
    },
    error = function(error) error
  )

  expect_true(
    mismatch_backend_called,
    "The injected subsetORFs backend was not called"
  )
  expect_true(
    inherits(mismatch_error, "error"),
    "An incomplete ORF subset must fail its exact-ID postcondition"
  )
})

if (length(failures) > 0L) {
  stop(
    paste(
      "P0 taxon/ORF filter regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: P0 taxon filtering is exact at the ORF level")
