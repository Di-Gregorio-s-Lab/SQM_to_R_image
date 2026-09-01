if (!requireNamespace("covr", quietly = TRUE)) {
  stop("Package 'covr' is required for the P3 coverage gate.", call. = FALSE)
}

source_lines <- readLines("sqm_plots.R", warn = FALSE)
main_call <- which(trimws(source_lines) == "main()")
stopifnot(length(main_call) == 1L)

temporary_source <- tempfile("sqm_plots_covr_p3_", fileext = ".R")
writeLines(source_lines[-main_call], temporary_source, useBytes = TRUE)
on.exit(unlink(temporary_source, force = TRUE), add = TRUE)

coverage <- covr::file_coverage(
  source_files = temporary_source,
  test_files = "tests/coverage_p3_exercises.R",
  parent_env = globalenv()
)
coverage_table <- as.data.frame(coverage)
target_functions <- c(
  "required_packages_for_mode",
  "check_required_packages",
  "check_flow_html_preflight",
  "stable_path_token",
  "safe_output_component",
  "taxonomy_file_stem",
  "portable_png_output_path",
  "assert_output_artifact",
  "save_png_dimensions",
  "save_html_widget",
  "build_pathview_input_table",
  "export_pathview_isolated",
  "generate_run_id",
  "allocate_run_id",
  "add_run_id_to_path",
  "resolve_sample_selection",
  "subset_pathway",
  "is_empty_pathway_subset",
  "prepare_context_pathway_subsets",
  "manifest_target_status",
  "validate_current_manifest_targets",
  "write_section_manifest",
  "write_combined_manifest",
  "list_run_artifacts",
  "write_run_manifest",
  "write_failed_run_manifests",
  "build_pie_chart_table",
  "make_pie_plot"
)

function_results <- lapply(target_functions, function(function_name) {
  rows <- coverage_table[coverage_table$functions == function_name, , drop = FALSE]
  if (nrow(rows) == 0L) {
    return(data.frame(function_name = function_name, covered = 0L, expressions = 0L, percent = 0))
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
for (row_index in seq_len(nrow(function_results))) {
  result <- function_results[row_index, ]
  message(
    result$function_name, ": ", sprintf("%.2f", result$percent),
    "% (", result$covered, "/", result$expressions, ")"
  )
}
message("P3 global sqm_plots.R coverage baseline: ", sprintf("%.2f", global_percent), "%")

failed <- function_results$function_name[function_results$percent < 80]
if (length(failed) > 0L) {
  stop("P3 coverage below 80% for: ", paste(failed, collapse = ", "), call. = FALSE)
}

message("PASS: P3 corrected functions meet the 80% coverage gate")
