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
    log_scale,
    output_dir,
    output_suffix) {
  calls[[length(calls) + 1L]] <<- list(
    SQM = SQM,
    pathway_id = pathway_id,
    count = count,
    samples = samples,
    split_samples = split_samples,
    log_scale = log_scale,
    output_dir = output_dir,
    output_suffix = output_suffix
  )
  file_name <- if (split_samples) "sample_CS8T0.png" else "combined.png"
  writeLines("fake pathview output", file.path(output_dir, file_name))
}

result <- script_env$run_pathview_mode(
  sqm_object = list(
    fake = TRUE,
    functions = list(KEGG = list(tpm = data.frame(
      CS8T0 = c(1, 2),
      CS8T2 = c(3, 4),
      row.names = c("K00001", "K00002")
    )))
  ),
  output_dir = test_root,
  manifest_base_dir = test_root,
  output_manifests = list(pathview = tibble::tibble()),
  script_name = "sqm_plots.R",
  project_dir = "project",
  tax_mode = "prokfilter",
  pathway_name = "Nitrogen metabolism",
  pathway_id = "00910",
  selected_samples = c("CS8T0", "CS8T2"),
  top_n_taxa = 15L,
  top_n_ko = 20L,
  pathview_sample_modes = c("insieme", "separato"),
  export_pathway_fn = fake_export_pathway
)

stopifnot(length(calls) == 3L)
stopifnot(all(!vapply(calls, `[[`, logical(1), "split_samples")))
stopifnot(all(!vapply(calls, `[[`, logical(1), "log_scale")))
stopifnot(identical(vapply(calls, `[[`, character(1), "count"), rep("tpm", 3L)))
stopifnot(length(unique(vapply(calls, `[[`, character(1), "output_dir"))) == 3L)
stopifnot(identical(
  sort(unique(result$pathview$output_scope)),
  c("pathway_defined_insieme", "pathway_defined_separato")
))
stopifnot(all(file.exists(file.path(test_root, result$pathview$output_file))))
separate_plots <- result$pathview$output_type == "pathview_file" &
  result$pathview$output_scope == "pathway_defined_separato"
stopifnot(setequal(result$pathview$samples[separate_plots], c("CS8T0", "CS8T2")))

input_rows <- result$pathview$output_type == "pathview_input_all_ko_complete_matrix_tsv"
stopifnot(sum(input_rows) == 3L)
stopifnot(all(grepl("pathview_input_all_ko_complete_matrix", result$pathview$output_file[input_rows], fixed = TRUE)))
config_path <- file.path(
  test_root,
  result$pathview$output_file[result$pathview$output_type == "pathview_render_config_tsv"][[1L]]
)
config <- readr::read_tsv(config_path, show_col_types = FALSE, na = "NA")
stopifnot(identical(config$log_scale, FALSE))
stopifnot(is.na(config$pseudocount))
stopifnot(identical(config$color_source, "pathview_native"))
stopifnot(identical(config$input_scope, "complete_all_ko_matrix"))

message("PASS: Pathview creates combined and split-sample output branches")
