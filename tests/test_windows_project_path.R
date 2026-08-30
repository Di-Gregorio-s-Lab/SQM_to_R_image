source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")
project_dir_native <- file.path(
  tempdir(), "Synthetic User", "Project Data", "SqueezeMeta", "Au_sip"
)
dir.create(project_dir_native, recursive = TRUE, showWarnings = FALSE)
on.exit(
  unlink(file.path(tempdir(), "Synthetic User"), recursive = TRUE, force = TRUE),
  add = TRUE
)

project_dir_forward <- normalizePath(project_dir_native, winslash = "/", mustWork = TRUE)
project_dir <- if (.Platform$OS.type == "windows") {
  gsub("/", "\\\\", project_dir_forward, fixed = TRUE)
} else {
  project_dir_forward
}

normalized_dir <- script_env$normalize_sqm_project_dir(project_dir)
stopifnot(identical(normalized_dir, project_dir_forward))

message("PASS: Windows SQM project paths are normalized with forward slashes")
