project_dir <- Sys.getenv("SQM_CS8_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) {
  stop("SQM_CS8_PROJECT_DIR is required for the Tranche 4 integration test.", call. = FALSE)
}
if (!dir.exists(project_dir)) {
  stop("SQM_CS8_PROJECT_DIR does not exist: ", project_dir, call. = FALSE)
}

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) {
  rscript <- file.path(R.home("bin"), "Rscript")
}
script_path <- normalizePath("sqm_plots.R", winslash = "/", mustWork = TRUE)
output_dir <- tempfile("t4_cs8_lifecycle_")
dir.create(output_dir, recursive = TRUE)
on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

run_flow <- function() {
  args <- c(
    shQuote(script_path),
    shQuote(paste0("--project_dir=", project_dir)),
    shQuote(paste0("--output_dir=", output_dir)),
    "--mode=flow",
    "--pathways=00361",
    "--pathway_selection_modes=defined",
    "--samples=CS8T0",
    "--taxonomy_ranks=phylum",
    "--flowplot_formats=png",
    "--top_n_ko=3",
    "--top_n_taxa=3",
    "--dimensions=3x2",
    "--plot_dpi=72"
  )
  output <- system2(rscript, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("CS8 FLOW run failed:\n", paste(output, collapse = "\n"), call. = FALSE)
  }
  invisible(output)
}

run_flow()
manifest_path <- file.path(output_dir, "flowplot", "manifest_flow.tsv")
stopifnot(file.exists(manifest_path))
first_manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE)
first_targets <- file.path(dirname(manifest_path), first_manifest$output_file)
stopifnot(length(first_targets) > 0L, all(file.exists(first_targets)))
stopifnot(!any(grepl("__[0-9]{8}T[0-9]{6}_UTC", basename(first_targets))))

sentinel <- file.path(output_dir, "user_sentinel.txt")
writeLines("preserve me", sentinel)
target_to_replace <- first_targets[[1L]]
writeLines("must be overwritten", target_to_replace)

run_flow()
second_manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE)
second_targets <- file.path(dirname(manifest_path), second_manifest$output_file)
stopifnot(file.exists(sentinel), identical(readLines(sentinel), "preserve me"))
stopifnot(file.exists(target_to_replace), !identical(readLines(target_to_replace), "must be overwritten"))
stopifnot(setequal(first_targets, second_targets))
stopifnot(length(list.files(output_dir, pattern = "^[0-9]{8}T.*\\.log$")) == 2L)
stopifnot(!file.exists(file.path(output_dir, "manifest_all.tsv")))
stopifnot(!file.exists(file.path(output_dir, "manifest_taxon.tsv")))
stopifnot(!file.exists(file.path(output_dir, "pathview", "manifest_pathview.tsv")))

message("PASS: CS8 rerun overwrites stable FLOW outputs and preserves unrelated files")
