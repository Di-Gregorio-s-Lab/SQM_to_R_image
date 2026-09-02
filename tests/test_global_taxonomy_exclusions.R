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

expect_error_contains <- function(expression, expected, label) {
  message <- tryCatch(
    {
      force(expression)
      NA_character_
    },
    error = function(error) conditionMessage(error)
  )
  expect_true(
    !is.na(message) && grepl(expected, message, fixed = TRUE),
    paste0(label, "\nObserved: ", message)
  )
}

capture_warnings <- function(expression) {
  messages <- character()
  value <- withCallingHandlers(
    force(expression),
    warning = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = messages)
}

make_taxonomy_sqm <- function(percent_matrix, total_reads = NULL) {
  if (is.null(total_reads)) {
    total_reads <- stats::setNames(
      rep(1000, ncol(percent_matrix)),
      colnames(percent_matrix)
    )
  }
  list(
    total_reads = total_reads,
    taxa = list(genus = list(percent = percent_matrix))
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

run_case("CS8 native plot retains Unclassified inside displayed data", {
  samples <- c("CS8T0", "CS8T6")
  raw_percent <- matrix(
    c(
      88.44857, 67.69785,
      0.01071923, 0.88092737,
      8, 20,
      3.54071077, 11.42122263
    ),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(
      c("Unmapped", "Unclassified", "Classified genus", "Rare genus"),
      samples
    )
  )
  plot_data <- tibble::tibble(
    sample = rep(samples, each = 2L),
    taxon = rep(c("Classified genus", "Other"), times = 2L),
    value = c(8, 3.55143, 20, 12.30215),
    count = "percent"
  )
  observed <- capture_warnings(
    script_env$add_global_taxonomy_percent_metadata(
      plot_data = plot_data,
      sqm_object = make_taxonomy_sqm(raw_percent),
      rank = "genus",
      selected_samples = samples,
      context_label = "global"
    )
  )
  annotated <- observed$value
  by_sample <- annotated[!duplicated(annotated$sample), , drop = FALSE]

  expect_true(
    identical(by_sample$raw_percent_sum, c(100, 100)),
    "Raw CS8 percentages were not retained"
  )
  expect_true(
    all(abs(by_sample$excluded_percent - c(88.44857, 67.69785)) < 1e-8),
    "Unmapped was not identified as the only effective exclusion"
  )
  expect_true(
    all(abs(by_sample$accounted_percent_sum - by_sample$raw_percent_sum) < 1e-8),
    "CS8 displayed and effective excluded percentages do not reconcile"
  )
  expect_true(
    all(by_sample$requested_excluded_categories == "Unmapped;Unclassified") &&
      all(by_sample$effective_excluded_categories == "Unmapped") &&
      all(by_sample$excluded_categories == "Unmapped") &&
      all(by_sample$retained_requested_categories == "Unclassified") &&
      all(by_sample$exclusion_resolution_status == "requested_retained"),
    "CS8 exclusion provenance is incomplete"
  )
  expect_true(
    length(observed$warnings) == 1L &&
      grepl("genus", observed$warnings, fixed = TRUE) &&
      grepl("Unclassified", observed$warnings, fixed = TRUE) &&
      grepl("Other", observed$warnings, fixed = TRUE),
    "CS8 retained exclusion did not emit one contextual warning"
  )
})

run_case("non-rescaled taxon context reconciles against its raw percentage", {
  raw_percent <- matrix(
    c(11.348, 0.0847),
    nrow = 1L,
    dimnames = list("Bacillota", c("S1", "S2"))
  )
  plot_data <- tibble::tibble(
    sample = c("S1", "S2"),
    taxon = "Bacillota",
    value = c(11.348, 0.0847),
    count = "percent"
  )
  observed <- capture_warnings(
    script_env$add_global_taxonomy_percent_metadata(
      plot_data = plot_data,
      sqm_object = make_taxonomy_sqm(raw_percent),
      rank = "genus",
      selected_samples = c("S1", "S2"),
      context_label = "Bacillota@phylum"
    )
  )
  annotated <- observed$value
  by_sample <- annotated[!duplicated(annotated$sample), , drop = FALSE]

  expect_true(
    identical(by_sample$raw_percent_sum, c(11.348, 0.0847)),
    "Filtered context was incorrectly forced to 100 percent"
  )
  expect_true(
    all(by_sample$excluded_percent == 0) &&
      all(by_sample$accounted_percent_sum == by_sample$raw_percent_sum),
    "Filtered context did not reconcile against its raw total"
  )
  expect_true(
    all(is.na(by_sample$effective_excluded_categories)) &&
      all(is.na(by_sample$retained_requested_categories)) &&
      all(by_sample$exclusion_resolution_status == "no_positive_requested_exclusions") &&
      length(observed$warnings) == 0L,
    "Absent requested exclusions were not treated as non-observable"
  )
})

run_case("both requested categories can be effectively excluded", {
  raw_percent <- matrix(
    c(60, 18.81372, 20, 1.18628),
    ncol = 1L,
    dimnames = list(
      c("Unmapped", "Unclassified", "Bacillota", "Rare taxon"),
      "S1"
    )
  )
  plot_data <- tibble::tibble(
    sample = c("S1", "S1"),
    taxon = c("Bacillota", "Other"),
    value = c(20, 1.18628),
    count = "percent"
  )
  observed <- capture_warnings(
    script_env$add_global_taxonomy_percent_metadata(
      plot_data = plot_data,
      sqm_object = make_taxonomy_sqm(raw_percent),
      rank = "genus",
      selected_samples = "S1",
      context_label = "global"
    )
  )
  annotated <- observed$value

  expect_true(
    all(annotated$effective_excluded_categories == "Unmapped;Unclassified") &&
      all(annotated$excluded_categories == "Unmapped;Unclassified") &&
      all(is.na(annotated$retained_requested_categories)) &&
      all(annotated$exclusion_resolution_status == "matched_requested") &&
      length(observed$warnings) == 0L,
    "Fully effective exclusions were not recorded"
  )
})

run_case("zero and absent requested categories are non-observable", {
  raw_percent <- matrix(
    c(0, 40),
    ncol = 1L,
    dimnames = list(c("Unmapped", "Alpha"), "S1")
  )
  plot_data <- tibble::tibble(
    sample = "S1",
    taxon = "Alpha",
    value = 40,
    count = "percent"
  )
  observed <- capture_warnings(
    script_env$add_global_taxonomy_percent_metadata(
      plot_data = plot_data,
      sqm_object = make_taxonomy_sqm(raw_percent),
      rank = "genus",
      selected_samples = "S1",
      context_label = "global"
    )
  )
  annotated <- observed$value

  expect_true(
    all(annotated$excluded_percent == 0) &&
      all(is.na(annotated$effective_excluded_categories)) &&
      all(is.na(annotated$retained_requested_categories)) &&
      all(annotated$exclusion_resolution_status == "no_positive_requested_exclusions") &&
      length(observed$warnings) == 0L,
    "Zero or absent exclusions changed the audit"
  )
})

run_case("unexplained percentage loss remains blocking", {
  raw_percent <- matrix(
    c(60, 10, 30),
    ncol = 1L,
    dimnames = list(c("Unmapped", "Unclassified", "Alpha"), "S1")
  )
  plot_data <- tibble::tibble(
    sample = "S1",
    taxon = "Alpha",
    value = 25,
    count = "percent"
  )
  expect_error_contains(
    script_env$add_global_taxonomy_percent_metadata(
      plot_data = plot_data,
      sqm_object = make_taxonomy_sqm(raw_percent),
      rank = "genus",
      selected_samples = "S1",
      context_label = "global"
    ),
    "cannot be reconciled",
    "An unexplained percentage loss was accepted"
  )
})

run_case("ambiguous effective exclusion remains blocking", {
  raw_percent <- matrix(
    c(10, 10, 80),
    ncol = 1L,
    dimnames = list(c("Unmapped", "Unclassified", "Alpha"), "S1")
  )
  plot_data <- tibble::tibble(
    sample = "S1",
    taxon = "Alpha and Other",
    value = 90,
    count = "percent"
  )
  expect_error_contains(
    script_env$add_global_taxonomy_percent_metadata(
      plot_data = plot_data,
      sqm_object = make_taxonomy_sqm(raw_percent),
      rank = "genus",
      selected_samples = "S1",
      context_label = "global"
    ),
    "ambiguous",
    "Two non-equivalent exclusion subsets were accepted"
  )
})

run_case("multi-sample vectors disambiguate equal aggregate exclusions", {
  raw_percent <- matrix(
    c(
      10, 20,
      20, 10,
      70, 70
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("Unmapped", "Unclassified", "Alpha"), c("S1", "S2"))
  )
  observed <- script_env$resolve_effective_taxonomy_exclusions(
    percent_frame = as.data.frame(raw_percent, check.names = FALSE),
    displayed_percent_sum = c(S1 = 90, S2 = 80),
    selected_samples = c("S1", "S2")
  )

  expect_true(
    identical(observed$effective_excluded_categories, "Unmapped") &&
      identical(observed$retained_requested_categories, "Unclassified"),
    "Per-sample vectors did not resolve exclusions with equal aggregate abundance"
  )
})

run_case("identical multi-sample exclusion vectors remain blocking", {
  raw_percent <- matrix(
    c(
      10, 20,
      10, 20,
      80, 60
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("Unmapped", "Unclassified", "Alpha"), c("S1", "S2"))
  )
  expect_error_contains(
    script_env$resolve_effective_taxonomy_exclusions(
      percent_frame = as.data.frame(raw_percent, check.names = FALSE),
      displayed_percent_sum = c(S1 = 90, S2 = 80),
      selected_samples = c("S1", "S2")
    ),
    "ambiguous",
    "Identical multi-sample exclusion vectors were accepted"
  )
})

if (length(failures) > 0L) {
  stop(
    paste(
      "Global taxonomy exclusion regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: global taxonomy exclusions reflect native SQMtools output")
