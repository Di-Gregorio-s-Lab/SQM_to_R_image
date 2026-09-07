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

read_bytes <- function(path) {
  readBin(path, what = "raw", n = file.info(path)$size)
}

write_fixture <- function(path, value) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(value, path, useBytes = TRUE)
  path
}

script_env <- source_without_main("sqm_plots.R")
test_root <- tempfile("t4_output_lifecycle_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

run_a <- "20260908T101112_UTCp0200_a111"
run_b <- "20260908T101113_UTCp0200_b222"

# Generated artifacts use semantic, stable names and overwrite in place.
script_env$set_run_context(run_a, cli_args = "--mode=flow")
stable_target <- file.path(test_root, "flowplot", "pathway", "flow.tsv")
expect_true(
  identical(script_env$run_artifact_path(stable_target), stable_target),
  "run_artifact_path() still adds a run-specific suffix"
)
first_path <- script_env$write_tsv_safe(tibble::tibble(value = "first"), stable_target)
script_env$set_run_context(run_b, cli_args = "--mode=flow")
second_path <- script_env$write_tsv_safe(tibble::tibble(value = "second"), stable_target)
expect_true(identical(first_path, second_path), "A rerun selected a different artifact path")
expect_true(
  identical(readr::read_tsv(stable_target, show_col_types = FALSE)$value, "second"),
  "A rerun did not overwrite the stable artifact"
)
expect_true(
  length(list.files(dirname(stable_target), pattern = "^flow.*\\.tsv$")) == 1L,
  "A rerun left multiple run-specific copies of one artifact"
)

# A physical output context owns three section-local manifest registries.
context_root <- file.path(test_root, "taxon_filter", "phylum", "Bacillota")
registry <- script_env$new_context_manifest_registry(context_root)
flow_one <- write_fixture(file.path(context_root, "flowplot", "p1", "one.tsv"), "one")
flow_old <- write_fixture(file.path(context_root, "flowplot", "p2", "old.tsv"), "old")
funz_one <- write_fixture(file.path(context_root, "funz", "p1", "one.tsv"), "one")
pie_one <- write_fixture(file.path(context_root, "pie", "p1", "one.tsv"), "one")
taxon_one <- write_fixture(file.path(context_root, "taxonomy_global", "one.tsv"), "one")
pathview_one <- write_fixture(file.path(context_root, "pathview", "one.tsv"), "one")

registry$flow <- tibble::tibble(
  run_id = run_a,
  output_file = c("p1/one.tsv", "p2/old.tsv")
)
registry$funz <- tibble::tibble(run_id = run_a, output_file = "p1/one.tsv")
registry$pie <- tibble::tibble(run_id = run_a, output_file = "p1/one.tsv")
registry$taxon <- tibble::tibble(
  run_id = run_a,
  output_file = script_env$relative_to_output(taxon_one, context_root)
)
registry$pathview <- tibble::tibble(
  run_id = run_a,
  output_file = script_env$relative_to_output(pathview_one, context_root)
)

written_manifests <- script_env$flush_context_manifests(registry)
expected_manifests <- c(
  flow = file.path(context_root, "flowplot", "manifest_flow.tsv"),
  funz = file.path(context_root, "funz", "manifest_funz.tsv"),
  pie = file.path(context_root, "pie", "manifest_pie.tsv")
)
expect_true(
  identical(unname(normalizePath(written_manifests, winslash = "/")),
            unname(normalizePath(expected_manifests, winslash = "/"))),
  "Manifest files are not local to their physical sections"
)
expect_true(
  !file.exists(file.path(context_root, "manifest_taxon.tsv")) &&
    !file.exists(file.path(context_root, "pathview", "manifest_pathview.tsv")) &&
    !file.exists(file.path(context_root, "manifest_all.tsv")),
  "A forbidden taxonomy, Pathview, or combined manifest was created"
)

# The next FLOW run replaces its manifest but preserves unselected files and sections.
funz_manifest <- expected_manifests[["funz"]]
funz_manifest_before <- read_bytes(funz_manifest)
flow_old_before <- read_bytes(flow_old)
script_env$set_run_context(run_b, cli_args = "--mode=flow")
registry_b <- script_env$new_context_manifest_registry(context_root)
write_fixture(flow_one, "one-rewritten")
registry_b$flow <- tibble::tibble(run_id = run_b, output_file = "p1/one.tsv")
script_env$flush_context_manifests(registry_b)
flow_rows <- readr::read_tsv(expected_manifests[["flow"]], show_col_types = FALSE)
expect_true(
  identical(flow_rows$output_file, "p1/one.tsv") && identical(flow_rows$run_id, run_b),
  "The FLOW manifest does not describe only the current run"
)
expect_true(identical(read_bytes(flow_old), flow_old_before), "An unselected FLOW file changed")
expect_true(identical(read_bytes(funz_manifest), funz_manifest_before), "An unselected section manifest changed")

# Every analytical run writes one <run_id>.log directly in output_dir.
script_env$set_run_context(
  run_a,
  started_at = as.POSIXct("2026-09-08 10:11:12", tz = "UTC"),
  samples = "S1",
  cli_args = c("--project_dir=project", paste0("--output_dir=", test_root), "--mode=flow")
)
log_a <- script_env$initialize_run_log(test_root, "project", "flow", "prokfilter")
script_env$progress_message("test progress")
script_env$record_run_warning("test warning")
script_env$finalize_run_log("SUCCESS")

legacy_manifest <- write_fixture(file.path(test_root, "manifest_all__legacy.tsv"), "legacy")
legacy_before <- read_bytes(legacy_manifest)
script_env$set_run_context(
  run_b,
  started_at = as.POSIXct("2026-09-08 10:11:13", tz = "UTC"),
  samples = "S1",
  cli_args = c("--project_dir=project", paste0("--output_dir=", test_root), "--mode=flow")
)
log_b <- script_env$initialize_run_log(test_root, "project", "flow", "prokfilter")
script_env$record_legacy_outputs(test_root)
script_env$finalize_run_log("FAILED", "controlled failure")

expect_true(identical(log_a, file.path(test_root, paste0(run_a, ".log"))), "Success log has the wrong path")
expect_true(identical(log_b, file.path(test_root, paste0(run_b, ".log"))), "Failure log has the wrong path")
expect_true(length(list.files(test_root, pattern = "^[0-9].*\\.log$")) == 2L, "Runs did not create distinct logs")
expect_true(grepl("STATUS=SUCCESS", paste(readLines(log_a), collapse = "\n"), fixed = TRUE), "Success status is missing from log")
failure_log <- paste(readLines(log_b), collapse = "\n")
expect_true(grepl("STATUS=FAILED", failure_log, fixed = TRUE), "Failure status is missing from log")
expect_true(grepl("controlled failure", failure_log, fixed = TRUE), "Failure reason is missing from log")
expect_true(grepl("LEGACY", failure_log, fixed = TRUE), "Legacy output presence is missing from log")
expect_true(identical(read_bytes(legacy_manifest), legacy_before), "Legacy manifest was modified")

script_env$clear_run_context()
message("PASS: stable overwrite, local manifests, partial preservation, and per-run logs")
