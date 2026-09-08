script_env <- if (identical(Sys.getenv("R_COVR"), "true")) environment() else new.env(parent = globalenv())
if (!identical(Sys.getenv("R_COVR"), "true")) sys.source("sqm_plots_lean.R", envir = script_env)

if (!exists("build_config", envir = script_env, mode = "function", inherits = FALSE)) {
  stop("Missing lean CLI helper: build_config", call. = FALSE)
}

config <- script_env$build_config(script_env$parse_args(c(
  "--project_dir", "in/Au_sip", "--output_dir=out/lean", "--mode", "huge"
)))
stopifnot(
  identical(config$mode, "huge"),
  identical(config$workers, script_env$default_workers()),
  identical(config$selection_modes, c("defined", "top20")),
  identical(config$pie_selection_modes, "defined"),
  identical(vapply(config$dimensions, `[[`, character(1L), "name"), c("12x9", "16x9", "12x16")),
  identical(config$formats, c("png", "html")),
  identical(config$pathview_sample_modes, c("insieme", "separato"))
)

normal <- script_env$build_config(script_env$parse_args(c(
  "--project_dir", "in/Au_sip", "--output_dir=out/lean", "--mode", "normal"
)))
stopifnot(
  identical(normal$mode, "normal"),
  identical(normal$selection_modes, "defined")
)

without_pie <- script_env$build_config(script_env$parse_args(c(
  "--project_dir", "in/Au_sip", "--output_dir", "out/lean", "--mode", "funz,flow,taxon,pathview"
)))
stopifnot(identical(without_pie$mode, c("funz", "flow", "taxon", "pathview")))

failed <- FALSE
tryCatch(
  script_env$build_config(script_env$parse_args(c(
    "--project_dir", "in/Au_sip", "--output_dir", "out/lean", "--mode", "enzimi"
  ))),
  error = function(...) failed <<- TRUE
)
stopifnot(failed)

failed <- FALSE
tryCatch(
  script_env$build_config(script_env$parse_args(c(
    "--project_dir", "in/Au_sip", "--output_dir", "out/lean", "--mode", "all"
  ))),
  error = function(...) failed <<- TRUE
)
stopifnot(failed)

failed <- FALSE
tryCatch(
  script_env$build_config(script_env$parse_args(c(
    "--project_dir", "in/Au_sip", "--output_dir", "out/lean", "--mode", "funz",
    "--enzyme_plot_types", "heatmap"
  ))),
  error = function(...) failed <<- TRUE
)
stopifnot(failed)

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
help <- system2(rscript, c("sqm_plots_lean.R", "--help"), stdout = TRUE, stderr = TRUE)
status <- attr(help, "status")
if (is.null(status)) status <- 0L
stopifnot(status == 0L, any(grepl("--plan_only", help, fixed = TRUE)))

ids <- c("orf1", "orf2")
sqm <- list(
  orfs = list(
    table = data.frame(`KEGG ID` = c("K00001", "K00002"), row.names = ids, check.names = FALSE),
    tax = data.frame(phylum = c("Alpha", "Beta"), row.names = ids),
    tpm = data.frame(S1 = c(2, 3), row.names = ids)
  ),
  functions = list(KEGG = list(tpm = data.frame(S1 = c(2, 3), row.names = c("K00001", "K00002"))))
)
stopifnot(identical(script_env$validate_sqm(sqm, "S1", "phylum"), sqm))
contexts <- script_env$make_contexts(sqm, list(taxa = "Alpha", output_dir = tempdir()))
stopifnot(
  length(contexts) == 1L,
  identical(contexts[[1L]]$filtered_taxon_rank, "phylum"),
  identical(rownames(contexts[[1L]]$sqm$orfs$table), "orf1")
)
global <- script_env$make_contexts(sqm, list(taxa = character(), output_dir = tempdir()))
stopifnot(length(global) == 1L, is.na(global[[1L]]$filtered_taxon))
stopifnot(
  identical(script_env$required_packages(list(plan_only = TRUE, mode = "huge", formats = "html")), c("SQMtools", "ggplot2")),
  all(c("ggalluvial", "plotly", "htmlwidgets", "pathview") %in%
        script_env$required_packages(list(plan_only = FALSE, mode = "normal", formats = "html")))
)
help_text <- capture.output(script_env$print_help())
stopifnot(
  any(grepl("--workers", help_text, fixed = TRUE)),
  any(grepl("huge", help_text, fixed = TRUE)),
  any(grepl("normal", help_text, fixed = TRUE)),
  !any(grepl("Modes: all", help_text, fixed = TRUE))
)

message("PASS: lean CLI exposes validated defaults without the legacy enzimi mode")
