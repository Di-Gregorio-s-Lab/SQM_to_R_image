args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript tests/verify_p3_candidate_outputs.R WINDOWS_DIR PIE_DIR ALL_DIR",
    call. = FALSE
  )
}

windows_root <- args[[1L]]
pie_root <- args[[2L]]
all_root <- args[[3L]]
for (candidate_root in args) {
  if (!dir.exists(candidate_root)) {
    stop("Candidate directory does not exist: ", candidate_root, call. = FALSE)
  }
}

is_safe_relative <- function(path) {
  !is.na(path) && nzchar(path) &&
    !grepl("^(?:[A-Za-z]:[/\\\\]|[/\\\\]{1,2})", path, perl = TRUE) &&
    !".." %in% strsplit(gsub("\\\\", "/", path), "/+", perl = TRUE)[[1L]]
}

latest_successful_run_id <- function(root) {
  run_manifests <- list.files(
    root,
    pattern = "^manifest_run__.*\\.tsv$",
    full.names = TRUE
  )
  if (length(run_manifests) == 0L) {
    stop("No run-specific manifest found in ", root, call. = FALSE)
  }
  metadata <- dplyr::bind_rows(lapply(
    run_manifests,
    readr::read_tsv,
    show_col_types = FALSE,
    na = "NA"
  )) |>
    dplyr::filter(.data$status == "success") |>
    dplyr::arrange(.data$finished_at)
  if (nrow(metadata) == 0L) {
    stop("No successful run found in ", root, call. = FALSE)
  }
  metadata$run_id[[nrow(metadata)]]
}

verify_manifest_tree <- function(root, expected_sections = NULL) {
  run_id <- latest_successful_run_id(root)
  manifest_all_path <- file.path(root, paste0("manifest_all__", run_id, ".tsv"))
  if (!file.exists(manifest_all_path) || file.info(manifest_all_path)$size <= 0) {
    stop("Missing or empty manifest_all.tsv in ", root, call. = FALSE)
  }
  manifest_all <- readr::read_tsv(manifest_all_path, show_col_types = FALSE, na = "NA")
  if (any(manifest_all$run_id != run_id)) {
    stop("manifest_all contains rows from another run in ", root, call. = FALSE)
  }
  if (!is.null(expected_sections) && !setequal(manifest_all$section, expected_sections)) {
    stop("manifest_all.tsv has unexpected sections in ", root, call. = FALSE)
  }
  if (nrow(manifest_all) == 0L || !all(vapply(manifest_all$manifest_file, is_safe_relative, logical(1)))) {
    stop("manifest_all.tsv contains unsafe or no manifest paths in ", root, call. = FALSE)
  }

  total_targets <- 0L
  for (manifest_relative in manifest_all$manifest_file) {
    manifest_path <- file.path(root, manifest_relative)
    if (!file.exists(manifest_path) || file.info(manifest_path)$size <= 0) {
      stop("Missing section manifest: ", manifest_relative, call. = FALSE)
    }
    section_manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE, na = "NA")
    if (any(section_manifest$run_id != run_id)) {
      stop("Section manifest mixes run IDs: ", manifest_relative, call. = FALSE)
    }
    if (nrow(section_manifest) == 0L || !all(vapply(section_manifest$output_file, is_safe_relative, logical(1)))) {
      stop("Section manifest is empty or unsafe: ", manifest_relative, call. = FALSE)
    }
    targets <- file.path(root, section_manifest$output_file)
    target_info <- file.info(targets)
    if (any(!file.exists(targets)) || anyNA(target_info$size) || any(target_info$size <= 0) || any(target_info$isdir)) {
      stop("Section manifest contains missing or empty targets: ", manifest_relative, call. = FALSE)
    }
    total_targets <- total_targets + length(targets)
  }
  structure(total_targets, run_id = run_id)
}

windows_targets <- verify_manifest_tree(windows_root)
windows_run_id <- attr(windows_targets, "run_id")
windows_png <- list.files(
  windows_root,
  pattern = paste0("__", windows_run_id, "\\.png$"),
  recursive = TRUE,
  full.names = TRUE
)
if (length(windows_png) == 0L || any(file.info(windows_png)$size <= 0)) {
  stop("Windows candidate has no non-empty PNG files.", call. = FALSE)
}
windows_png_abs <- normalizePath(windows_png, winslash = "/", mustWork = TRUE)
if (any(nchar(windows_png_abs, type = "chars") > 240L)) {
  stop("Windows candidate contains a PNG path over 240 characters.", call. = FALSE)
}
compact_pattern <- paste0(
  "__[[:xdigit:]]{12}_[^/]+__", windows_run_id, "\\.png$"
)
if (!any(grepl(compact_pattern, windows_png_abs, perl = TRUE))) {
  stop("Windows candidate did not exercise deterministic PNG compaction.", call. = FALSE)
}

pie_targets <- verify_manifest_tree(pie_root)
pie_run_id <- attr(pie_targets, "run_id")
pie_tsv <- list.files(
  pie_root,
  pattern = paste0("_data__", pie_run_id, "\\.tsv$"),
  recursive = TRUE,
  full.names = TRUE
)
if (length(pie_tsv) == 0L) {
  stop("PIE candidate has no exported data TSV.", call. = FALSE)
}
required_pie_columns <- c(
  "pathway", "pathway_id", "pathway_selection", "rank", "ko_name",
  "ko_sample_tpm", "pathway_sample_tpm", "ko_pathway_percent", "taxon_order"
)
for (pie_path in pie_tsv) {
  pie_data <- readr::read_tsv(pie_path, show_col_types = FALSE, na = "NA")
  if (!all(required_pie_columns %in% colnames(pie_data))) {
    stop("Incomplete PIE TSV: ", pie_path, call. = FALSE)
  }
  if (!isTRUE(all.equal(sum(pie_data$tpm), unique(pie_data$ko_sample_tpm), tolerance = 1e-10)) ||
      !isTRUE(all.equal(sum(pie_data$pct), 1, tolerance = 1e-10))) {
    stop("PIE denominator invariant failed: ", pie_path, call. = FALSE)
  }
}

all_targets <- verify_manifest_tree(
  all_root,
  expected_sections = c("flow", "funz", "taxon", "pathview", "pie")
)

message(
  "PASS: P3 candidate outputs verified | windows_targets=", windows_targets,
  " | windows_png=", length(windows_png),
  " | pie_targets=", pie_targets,
  " | pie_tsv=", length(pie_tsv),
  " | all_targets=", all_targets
)
