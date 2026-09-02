if (!requireNamespace("covr", quietly = TRUE)) {
  stop("Package 'covr' is required for the P0 coverage gate.", call. = FALSE)
}

source_lines <- readLines("sqm_plots.R", warn = FALSE)
main_call <- which(trimws(source_lines) == "main()")
stopifnot(length(main_call) == 1L)

temporary_source <- tempfile("sqm_plots_covr_", fileext = ".R")
writeLines(source_lines[-main_call], temporary_source, useBytes = TRUE)
on.exit(unlink(temporary_source, force = TRUE), add = TRUE)

coverage <- covr::file_coverage(
  source_files = temporary_source,
  test_files = "tests/coverage_p0_exercises.R",
  parent_env = globalenv()
)
coverage_table <- as.data.frame(coverage)

target_functions <- c(
  "resolve_pathways",
  "pathview_is_exportable",
  "resolve_taxa_filters",
  "subset_sqm_by_taxon",
  "resolve_effective_taxonomy_exclusions",
  "add_global_taxonomy_percent_metadata",
  "build_pathway_taxonomy_percent_table",
  "make_pathway_taxonomy_percent_plot"
)

function_results <- lapply(target_functions, function(function_name) {
  rows <- coverage_table[coverage_table$functions == function_name, , drop = FALSE]
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
message("Instrumented sqm_plots.R baseline: ", sprintf("%.2f%%", global_percent))
for (row_index in seq_len(nrow(function_results))) {
  result <- function_results[row_index, ]
  message(
    result$function_name, ": ", sprintf("%.2f%%", result$percent),
    " (", result$covered, "/", result$expressions, " expressions)"
  )
}

failed_functions <- function_results$function_name[
  function_results$percent < 80
]
if (length(failed_functions) > 0L) {
  stop(
    "P0 coverage below 80% for: ",
    paste(failed_functions, collapse = ", "),
    call. = FALSE
  )
}

message("PASS: all directly corrected P0 functions meet the 80% coverage gate")
