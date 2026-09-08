source("tests/helpers/t1_sqmtools_oracles.R")

script_env <- t1_source_sqm_plots_without_main()

flow_fixture <- tibble::tibble(
  pathway = "Synthetic pathway",
  rank = "family",
  sample = "S_flow",
  taxon = factor(
    c("Tax_Z", "Tax_A", "Tax_M", "Tax_Z"),
    levels = c("Tax_Z", "Tax_A", "Tax_M")
  ),
  KO = factor(
    c("K00002", "K00001", "K00002", "Other"),
    levels = c("Other", "K00002", "K00001")
  ),
  KO_name = c("Function two", "Function one", "Function two", "Other KOs"),
  ec_codes = c("2.2.2.2", "1.1.1.1", "2.2.2.2", NA_character_),
  TPM = c(25, 20, 30, 25),
  taxon_percent = c(50, 20, 30, 50),
  KO_percent = c(55, 20, 55, 25),
  flow_percent = c(25, 20, 30, 25)
)

widget <- script_env$make_flow_sankey(
  flow_fixture,
  "Synthetic pathway",
  "family",
  "S_flow"
)
trace <- plotly::plotly_build(widget)$x$data[[1L]]

if (!identical(as.character(trace$arrangement), "snap")) {
  stop(
    "Simple FLOW HTML must use Plotly's automatic snap layout; observed ",
    as.character(trace$arrangement), ".",
    call. = FALSE
  )
}
if (length(trace$node$x) > 0L || length(trace$node$y) > 0L) {
  stop("Simple FLOW HTML must not force node coordinates.", call. = FALSE)
}

observed_labels <- as.character(unlist(trace$node$label, use.names = FALSE))
expected_labels <- c("Tax_Z", "Tax_A", "Tax_M", "Other KOs", "K00002", "K00001")
if (!identical(observed_labels, expected_labels)) {
  stop(
    "Simple FLOW HTML node labels changed; expected ",
    paste(expected_labels, collapse = ","),
    "; observed ", paste(observed_labels, collapse = ","),
    call. = FALSE
  )
}

message("PASS: simple FLOW HTML uses automatic layout and plain labels")

output_root <- tempfile("t5_flow_html_")
dir.create(output_root, recursive = TRUE)
on.exit(unlink(output_root, recursive = TRUE, force = TRUE), add = TRUE)
html_path <- file.path(output_root, "flowplot.html")
dependency_dir <- file.path(output_root, "flowplot_files")

script_env$save_html_widget(widget, html_path)
if (!file.exists(html_path) || file.info(html_path)$size <= 0L) {
  stop("FLOW HTML writer did not create a non-empty HTML file.", call. = FALSE)
}
dependency_files <- list.files(
  dependency_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
)
if (!dir.exists(dependency_dir) || length(dependency_files) == 0L) {
  stop(
    "FLOW HTML must create its support directory: ",
    dependency_dir,
    call. = FALSE
  )
}
html_text <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
if (!grepl("flowplot_files", html_text, fixed = TRUE)) {
  stop("FLOW HTML does not reference its support directory.", call. = FALSE)
}

first_dependency_count <- length(dependency_files)
script_env$save_html_widget(widget, html_path)
second_dependency_count <- length(list.files(
  dependency_dir,
  recursive = TRUE,
  all.files = TRUE,
  no.. = TRUE
))
if (second_dependency_count != first_dependency_count) {
  stop("FLOW HTML support directory is not stable across overwrite.", call. = FALSE)
}

message("PASS: FLOW HTML creates and overwrites its support directory")
