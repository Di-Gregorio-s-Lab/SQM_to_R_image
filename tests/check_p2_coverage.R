if (!requireNamespace("covr", quietly = TRUE)) {
  stop("Package 'covr' is required for the P2 coverage gate.", call. = FALSE)
}

source_lines <- readLines("sqm_plots.R", warn = FALSE)
main_call <- which(trimws(source_lines) == "main()")
stopifnot(length(main_call) == 1L)

temporary_source <- tempfile("sqm_plots_covr_p2_", fileext = ".R")
writeLines(source_lines[-main_call], temporary_source, useBytes = TRUE)
on.exit(unlink(temporary_source, force = TRUE), add = TRUE)

coverage <- covr::file_coverage(
  source_files = temporary_source,
  test_files = "tests/coverage_p2_exercises.R",
  parent_env = globalenv()
)
coverage_table <- as.data.frame(coverage)

target_functions <- c(
  "select_top_classified_taxa",
  "collapse_taxa_preserving_unclassified",
  "extract_ko_ec_lookup",
  "build_ko_plot_table",
  "build_pie_chart_table",
  "parse_positive_integer_arg",
  "build_ko_expansion_audit",
  "build_orf_long_result",
  "normalize_manifest_ko_audit"
)

function_results <- lapply(target_functions, function(function_name) {
  rows <- coverage_table[
    coverage_table$functions == function_name,
    ,
    drop = FALSE
  ]
  if (nrow(rows) == 0L) {
    return(data.frame(
      function_name = function_name,
      covered = 0L,
      expressions = 0L,
      percent = 0
    ))
  }

  data.frame(
    function_name = function_name,
    covered = sum(rows$value > 0),
    expressions = nrow(rows),
    percent = 100 * sum(rows$value > 0) / nrow(rows)
  )
})
function_results <- do.call(rbind, function_results)

global_percent <- covr::percent_coverage(coverage)
message("covr version: ", as.character(utils::packageVersion("covr")))
message("Instrumented sqm_plots.R P2 baseline: ", sprintf("%.2f%%", global_percent))
for (row_index in seq_len(nrow(function_results))) {
  result <- function_results[row_index, ]
  message(
    result$function_name,
    ": ", sprintf("%.2f%%", result$percent),
    " (", result$covered, "/", result$expressions, ")"
  )
}

failed <- function_results[function_results$percent < 80, , drop = FALSE]
if (nrow(failed) > 0L) {
  stop(
    "P2 coverage below 80% for: ",
    paste(
      paste0(failed$function_name, "=", sprintf("%.2f%%", failed$percent)),
      collapse = ", "
    ),
    call. = FALSE
  )
}

message("PASS: P2 corrected functions meet the 80% coverage gate")
