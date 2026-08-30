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

capture_error <- function(code) {
  tryCatch(
    {
      force(code)
      NULL
    },
    error = function(error) conditionMessage(error)
  )
}

write_nonempty_file <- function(path, contents = "fixture") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(contents, path, useBytes = TRUE)
  stopifnot(file.exists(path), file.info(path)$size > 0)
  invisible(path)
}

manifest_row <- function(output_file, marker) {
  tibble::tibble(
    output_file = output_file,
    marker = marker
  )
}

read_file_bytes <- function(path) {
  readBin(path, what = "raw", n = file.info(path)$size)
}

script_env <- source_without_main("sqm_plots.R")

# Existing valid targets survive, stale targets are pruned with a section-specific
# warning, and a new row wins when it has the same output_file as an old row.
section_root <- tempfile("p3_manifest_section_")
dir.create(section_root, recursive = TRUE)
on.exit(unlink(section_root, recursive = TRUE, force = TRUE), add = TRUE)

shared_relative <- "taxon/shared.tsv"
new_relative <- "taxon/new.tsv"
stale_relative <- "taxon/stale.tsv"
write_nonempty_file(file.path(section_root, shared_relative), "shared")
write_nonempty_file(file.path(section_root, new_relative), "new")

section_manifest_path <- file.path(section_root, "manifest_taxon.tsv")
existing_rows <- dplyr::bind_rows(
  manifest_row(shared_relative, "old"),
  manifest_row(stale_relative, "stale")
)
script_env$write_tsv_safe(existing_rows, section_manifest_path)

captured_warnings <- character()
written_path <- withCallingHandlers(
  script_env$write_section_manifest(
    dplyr::bind_rows(
      manifest_row(shared_relative, "new-wins"),
      manifest_row(new_relative, "new-only")
    ),
    section_root,
    "",
    "manifest_taxon.tsv"
  ),
  warning = function(warning_condition) {
    captured_warnings <<- c(captured_warnings, conditionMessage(warning_condition))
    invokeRestart("muffleWarning")
  }
)

expect_true(
  identical(normalizePath(written_path, winslash = "/"), normalizePath(section_manifest_path, winslash = "/")),
  "write_section_manifest() returned a different manifest path"
)
merged_rows <- readr::read_tsv(section_manifest_path, show_col_types = FALSE, na = "NA")
expect_true(
  setequal(merged_rows$output_file, c(shared_relative, new_relative)),
  "The section manifest did not prune the stale target or retain every valid target"
)
expect_true(
  identical(merged_rows$marker[merged_rows$output_file == shared_relative], "new-wins"),
  "An existing manifest row took precedence over the current run"
)
expect_true(
  any(grepl("taxon", captured_warnings, ignore.case = TRUE)),
  "The stale-target warning does not identify the taxon section"
)
expect_true(
  any(grepl("1", captured_warnings, fixed = TRUE)),
  "The stale-target warning does not report the number of pruned rows"
)

# Invalid rows from the current run must fail transactionally: the existing
# manifest must retain both its exact bytes and its original modification time.
baseline_rows <- manifest_row(shared_relative, "baseline")
fixed_mtime <- as.POSIXct("2001-02-03 04:05:06", tz = "UTC")

expect_rejected_without_rewrite <- function(candidate_row, error_pattern, label) {
  script_env$write_tsv_safe(baseline_rows, section_manifest_path)
  Sys.setFileTime(section_manifest_path, fixed_mtime)
  before_bytes <- read_file_bytes(section_manifest_path)
  before_mtime <- file.info(section_manifest_path)$mtime

  error_message <- capture_error(script_env$write_section_manifest(
    candidate_row,
    section_root,
    "",
    "manifest_taxon.tsv"
  ))

  expect_true(!is.null(error_message), paste0(label, " was accepted"))
  expect_true(
    grepl(error_pattern, error_message, ignore.case = TRUE),
    paste0(label, " produced an unclear error: ", error_message)
  )
  expect_true(
    identical(read_file_bytes(section_manifest_path), before_bytes),
    paste0(label, " rewrote the existing manifest contents")
  )
  expect_true(
    identical(file.info(section_manifest_path)$mtime, before_mtime),
    paste0(label, " changed the existing manifest modification time")
  )
}

expect_rejected_without_rewrite(
  manifest_row("taxon/missing.tsv", "missing"),
  "missing|does not exist|target",
  "A missing current-run target"
)

absolute_target <- normalizePath(
  file.path(section_root, shared_relative),
  winslash = "/",
  mustWork = TRUE
)
expect_rejected_without_rewrite(
  manifest_row(absolute_target, "absolute"),
  "relative|absolute|unsafe",
  "An absolute current-run target"
)

outside_target <- file.path(dirname(section_root), "outside_manifest_target.tsv")
write_nonempty_file(outside_target, "outside")
on.exit(unlink(outside_target, force = TRUE), add = TRUE)
expect_rejected_without_rewrite(
  manifest_row("../outside_manifest_target.tsv", "traversal"),
  "relative|\\.\\.|traversal|unsafe|outside",
  "A path-traversal current-run target"
)

# The combined manifest is rebuilt from all section manifests under output_dir.
# Sections reduced to zero valid rows are omitted, and every referenced manifest
# and artifact must exist and be non-empty.
combined_root <- tempfile("p3_manifest_all_")
dir.create(combined_root, recursive = TRUE)
on.exit(unlink(combined_root, recursive = TRUE, force = TRUE), add = TRUE)

flow_relative <- "flowplot/valid_flow.tsv"
pie_relative <- "pie/definiti/valid_pie.tsv"
write_nonempty_file(file.path(combined_root, flow_relative), "flow")
write_nonempty_file(file.path(combined_root, pie_relative), "pie")

flow_manifest_path <- file.path(combined_root, "flowplot", "manifest_flow.tsv")
pie_manifest_path <- file.path(combined_root, "pie", "manifest_pie.tsv")
taxon_manifest_path <- file.path(combined_root, "manifest_taxon.tsv")
dir.create(dirname(flow_manifest_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(pie_manifest_path), recursive = TRUE, showWarnings = FALSE)
script_env$write_tsv_safe(manifest_row(flow_relative, "flow"), flow_manifest_path)
script_env$write_tsv_safe(manifest_row(pie_relative, "pie"), pie_manifest_path)
script_env$write_tsv_safe(manifest_row("taxon/missing.tsv", "stale"), taxon_manifest_path)

# This stale aggregate row must never survive a rebuild-from-zero operation.
script_env$write_tsv_safe(
  tibble::tibble(section = "stale", manifest_file = "missing_manifest.tsv"),
  file.path(combined_root, "manifest_all.tsv")
)

combined_warnings <- character()
combined_manifest_path <- withCallingHandlers(
  script_env$write_combined_manifest(combined_root),
  warning = function(warning_condition) {
    combined_warnings <<- c(combined_warnings, conditionMessage(warning_condition))
    invokeRestart("muffleWarning")
  }
)

expect_true(
  identical(
    normalizePath(combined_manifest_path, winslash = "/", mustWork = TRUE),
    normalizePath(file.path(combined_root, "manifest_all.tsv"), winslash = "/", mustWork = TRUE)
  ),
  "write_combined_manifest() returned a different aggregate manifest path"
)
combined_rows <- readr::read_tsv(combined_manifest_path, show_col_types = FALSE, na = "NA")
expect_true(
  identical(names(combined_rows), c("section", "manifest_file")),
  "manifest_all.tsv changed its public schema"
)
expect_true(
  setequal(combined_rows$section, c("flow", "pie")),
  "manifest_all.tsv included a stale/empty section or omitted a valid section"
)
expect_true(
  any(grepl("taxon", combined_warnings, ignore.case = TRUE)) &&
    any(grepl("1", combined_warnings, fixed = TRUE)),
  "Combined reconciliation did not warn about the one stale taxon row"
)

for (manifest_relative in combined_rows$manifest_file) {
  manifest_path <- file.path(combined_root, manifest_relative)
  expect_true(file.exists(manifest_path), "manifest_all.tsv references a missing manifest")
  expect_true(file.info(manifest_path)$size > 0, "manifest_all.tsv references an empty manifest file")

  section_rows <- readr::read_tsv(manifest_path, show_col_types = FALSE, na = "NA")
  expect_true(nrow(section_rows) > 0L, "manifest_all.tsv references a manifest with no valid rows")
  for (output_relative in section_rows$output_file) {
    output_path <- file.path(combined_root, output_relative)
    expect_true(file.exists(output_path), "A section manifest references a missing output target")
    expect_true(file.info(output_path)$size > 0, "A section manifest references an empty output target")
  }
}

pruned_taxon <- readr::read_tsv(taxon_manifest_path, show_col_types = FALSE, na = "NA")
expect_true(nrow(pruned_taxon) == 0L, "Combined reconciliation did not leave the stale-only manifest empty")

message("PASS: P3 manifest target integrity and aggregate reconciliation")
