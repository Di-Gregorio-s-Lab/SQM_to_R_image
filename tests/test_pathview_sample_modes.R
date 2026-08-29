source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")

stopifnot(identical(
  script_env$normalize_pathview_sample_modes(NULL),
  c("insieme", "separato")
))
stopifnot(identical(
  script_env$normalize_pathview_sample_modes("separato"),
  "separato"
))

invalid_mode_error <- tryCatch(
  {
    script_env$normalize_pathview_sample_modes("non_valido")
    NULL
  },
  error = function(error) error$message
)
stopifnot(grepl("pathview_sample_modes", invalid_mode_error, fixed = TRUE))

test_root <- tempfile("pathview_sample_modes_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

calls <- list()
fake_export_pathway <- function(
    SQM,
    pathway_id,
    count,
    samples,
    split_samples,
    output_dir,
    output_suffix) {
  calls[[length(calls) + 1L]] <<- list(
    SQM = SQM,
    pathway_id = pathway_id,
    count = count,
    samples = samples,
    split_samples = split_samples,
    output_dir = output_dir,
    output_suffix = output_suffix
  )
  file_name <- if (split_samples) "sample_CS8T0.png" else "combined.png"
  writeLines("fake pathview output", file.path(output_dir, file_name))
}

result <- script_env$run_pathview_mode(
  sqm_object = list(fake = TRUE),
  output_dir = test_root,
  manifest_base_dir = test_root,
  output_manifests = list(pathview = tibble::tibble()),
  script_name = "sqm_plots.R",
  project_dir = "project",
  tax_mode = "prokfilter",
  pathway_name = "Nitrogen metabolism",
  pathway_id = "00643",
  selected_samples = c("CS8T0", "CS8T2"),
  top_n_taxa = 15L,
  top_n_ko = 20L,
  pathview_sample_modes = c("insieme", "separato"),
  export_pathway_fn = fake_export_pathway
)

stopifnot(length(calls) == 2L)
stopifnot(identical(vapply(calls, `[[`, logical(1), "split_samples"), c(FALSE, TRUE)))
stopifnot(identical(vapply(calls, `[[`, character(1), "count"), c("tpm", "tpm")))
stopifnot(identical(
  vapply(calls, `[[`, character(1), "output_dir"),
  file.path(
    normalizePath(test_root, winslash = "/", mustWork = TRUE),
    "pathview",
    "definiti",
    c("insieme", "separato"),
    "Nitrogen_metabolism"
  )
))
stopifnot(identical(
  sort(result$pathview$output_scope),
  c("pathway_defined_insieme", "pathway_defined_separato")
))
stopifnot(all(file.exists(file.path(
  test_root,
  "pathview",
  "definiti",
  c("insieme", "separato"),
  "Nitrogen_metabolism",
  c("combined.png", "sample_CS8T0.png")
))))

message("PASS: Pathview creates combined and split-sample output branches")
