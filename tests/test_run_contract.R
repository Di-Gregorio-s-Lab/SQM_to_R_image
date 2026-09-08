source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")

run_id <- script_env$generate_run_id(
  now = as.POSIXct("2026-09-01 14:35:27", tz = "UTC"),
  utc_offset = "+0200",
  collision_suffix = "a7f3"
)
stopifnot(identical(run_id, "20260901T143527_UTCp0200_a7f3"))
stopifnot(identical(
  script_env$add_run_id_to_path("results/barplot.tsv", run_id),
  paste0("results/barplot__", run_id, ".tsv")
))
stopifnot(identical(
  script_env$add_run_id_to_path(paste0("results/barplot__", run_id, ".tsv"), run_id),
  paste0("results/barplot__", run_id, ".tsv")
))

cli_selection <- script_env$resolve_sample_selection("S2,S1", c("S1", "S2", "S3"))
stopifnot(identical(cli_selection$samples, c("S2", "S1")))
stopifnot(identical(cli_selection$basis, "cli"))
sqm_selection <- script_env$resolve_sample_selection(NULL, c("S3", "S1", "S2"))
stopifnot(identical(sqm_selection$samples, c("S3", "S1", "S2")))
stopifnot(identical(sqm_selection$basis, "sqm_column_order"))

fallback_root <- tempfile("fallback_tsv_contract_")
dir.create(fallback_root, recursive = TRUE)
on.exit(unlink(fallback_root, recursive = TRUE, force = TRUE), add = TRUE)
fallback_source <- data.frame(
  field = c("tab\tinside", "line one\nline two", "quote \"inside\""),
  stringsAsFactors = FALSE
)
fallback_path <- script_env$write_tsv_safe(
  fallback_source,
  file.path(fallback_root, "fallback.tsv"),
  readr_available = FALSE
)
fallback_roundtrip <- utils::read.delim(
  fallback_path,
  sep = "\t",
  header = TRUE,
  quote = "\"",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(identical(fallback_roundtrip$field, fallback_source$field))

subset_calls <- list()
fake_subset_fun <- function(
    SQM, fun, columns, ignore_case, fixed, allow_empty) {
  subset_calls[[length(subset_calls) + 1L]] <<- list(
    fun = fun,
    allow_empty = allow_empty
  )
  if (identical(fun, "Empty pathway")) {
    return(list(orfs = list(table = data.frame(), tpm = data.frame())))
  }
  list(
    orfs = list(
      table = data.frame(value = 1, row.names = "orf_1"),
      tpm = data.frame(S1 = 1, row.names = "orf_1")
    )
  )
}
pathway_entries <- list(
  list(pathway_name = "Present pathway", pathway_id = "00001", pathway_selection = "defined"),
  list(pathway_name = "Empty pathway", pathway_id = "00002", pathway_selection = "defined")
)
warnings <- character()
prepared <- withCallingHandlers(
  script_env$prepare_context_pathway_subsets(
    context_sqm = list(),
    pathway_entries = pathway_entries,
    context_label = "global",
    subset_fun = fake_subset_fun
  ),
  warning = function(condition) {
    warnings <<- c(warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
)
stopifnot(all(vapply(subset_calls, `[[`, logical(1), "allow_empty")))
stopifnot(script_env$is_empty_pathway_subset(list()))
stopifnot(length(prepared$pathway_sqms) == 1L)
stopifnot(nrow(prepared$skips) == 1L)
stopifnot(identical(prepared$skips$pathway[[1L]], "Empty pathway"))
stopifnot(any(grepl("Empty pathway", warnings, fixed = TRUE)))

manifest_root <- tempfile("run_manifest_contract_")
dir.create(manifest_root, recursive = TRUE)
on.exit(unlink(manifest_root, recursive = TRUE, force = TRUE), add = TRUE)

script_env$set_run_context(
  run_id = run_id,
  started_at = as.POSIXct("2026-09-01 14:35:27", tz = "UTC"),
  sample_order_basis = "cli",
  samples = c("S2", "S1"),
  cli_args = c("--samples=S2,S1")
)
artifact_a <- script_env$add_run_id_to_path(file.path(manifest_root, "taxonomy.tsv"), run_id)
writeLines("current", artifact_a)
manifest_a <- script_env$write_section_manifest(
  tibble::tibble(
    run_id = run_id,
    output_file = basename(artifact_a),
    marker = "current-a"
  ),
  manifest_root,
  "",
  "manifest_taxon.tsv"
)
stopifnot(identical(basename(manifest_a), "manifest_taxon.tsv"))

run_id_b <- sub("a7f3$", "b8e4", run_id)
script_env$set_run_context(
  run_id = run_id_b,
  started_at = as.POSIXct("2026-09-01 14:35:27", tz = "UTC"),
  sample_order_basis = "sqm_column_order",
  samples = c("S1", "S2"),
  cli_args = character()
)
artifact_b <- script_env$add_run_id_to_path(file.path(manifest_root, "taxonomy.tsv"), run_id_b)
writeLines("current", artifact_b)
manifest_b <- script_env$write_section_manifest(
  tibble::tibble(
    run_id = run_id_b,
    output_file = basename(artifact_b),
    marker = "current-b"
  ),
  manifest_root,
  "",
  "manifest_taxon.tsv"
)
rows_b <- readr::read_tsv(manifest_b, show_col_types = FALSE)
stopifnot(nrow(rows_b) == 1L)
stopifnot(identical(rows_b$run_id, run_id_b))
stopifnot(file.exists(manifest_a), file.exists(manifest_b))

combined_b <- script_env$write_combined_manifest(
  manifest_root,
  section_manifest_paths = c(taxon = manifest_b)
)
combined_rows <- readr::read_tsv(combined_b, show_col_types = FALSE)
stopifnot(identical(combined_rows$run_id, run_id_b))
stopifnot(identical(combined_rows$section, "taxon"))
stopifnot(identical(combined_rows$manifest_file, basename(manifest_b)))

failure_root <- tempfile("failed_run_contract_")
dir.create(failure_root, recursive = TRUE)
on.exit(unlink(failure_root, recursive = TRUE, force = TRUE), add = TRUE)
script_env$set_run_context(
  run_id = run_id,
  started_at = as.POSIXct("2026-09-01 14:35:27", tz = "UTC"),
  sample_order_basis = "cli",
  samples = c("S2", "S1"),
  cli_args = c("--samples=S2,S1")
)
partial_path <- script_env$add_run_id_to_path(file.path(failure_root, "partial.tsv"), run_id)
writeLines("partial", partial_path)
script_env$record_pathway_skips(tibble::tibble(
  context = "global",
  pathway = "Empty pathway",
  reason = "empty_subset"
))
failure_manifests <- script_env$write_failed_run_manifests(
  output_dir = failure_root,
  project_dir = "project",
  mode = "all",
  tax_mode = "prokfilter",
  error_message = "controlled failure"
)
stopifnot(file.exists(failure_manifests$artifacts), file.exists(failure_manifests$run))
failed_artifacts <- readr::read_tsv(failure_manifests$artifacts, show_col_types = FALSE)
run_metadata <- readr::read_tsv(failure_manifests$run, show_col_types = FALSE, na = "NA")
stopifnot(identical(failed_artifacts$output_file, basename(partial_path)))
stopifnot(identical(run_metadata$status, "failed"))
stopifnot(identical(run_metadata$error_message, "controlled failure"))
stopifnot(grepl("Empty pathway", run_metadata$skipped_pathways, fixed = TRUE))

script_env$clear_run_context()
message("PASS: run identity and empty-pathway isolation are enforced")
