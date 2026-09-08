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
