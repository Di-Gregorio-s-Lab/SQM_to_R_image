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
    c("K00002", "K00001", "K00002", "K00003"),
    levels = c("K00003", "K00002", "K00001")
  ),
  KO_name = c("Function two", "Function one", "Function two", "Function three"),
  ec_codes = NA_character_,
  TPM = c(25, 20, 30, 25),
  taxon_percent = c(50, 20, 30, 50),
  KO_percent = c(55, 20, 55, 25),
  flow_percent = c(25, 20, 30, 25)
)

plot_object <- script_env$make_flow_plot(
  flow_fixture,
  "Synthetic pathway",
  "family",
  "S_flow"
)
built_plot <- ggplot2::ggplot_build(plot_object)

flow_intervals <- built_plot$data[[1L]] |>
  dplyr::filter(.data$x == 1) |>
  dplyr::group_by(.data$stratum) |>
  dplyr::summarise(
    flow_ymin = min(.data$ymin),
    flow_ymax = max(.data$ymax),
    .groups = "drop"
  ) |>
  dplyr::mutate(stratum = as.character(.data$stratum))

stratum_intervals <- built_plot$data[[2L]] |>
  dplyr::transmute(
    stratum = as.character(.data$stratum),
    stratum_ymin = .data$ymin,
    stratum_ymax = .data$ymax
  )

comparison <- dplyr::left_join(
  flow_intervals,
  stratum_intervals,
  by = "stratum",
  relationship = "one-to-one"
) |>
  dplyr::mutate(
    delta = pmax(
      abs(.data$flow_ymin - .data$stratum_ymin),
      abs(.data$flow_ymax - .data$stratum_ymax)
    )
  )

if (nrow(comparison) != length(levels(flow_fixture$taxon)) ||
    anyNA(comparison$delta) ||
    any(comparison$delta > 1e-8)) {
  details <- apply(
    comparison,
    1L,
    function(row) paste0(
      row[["stratum"]],
      ": flow=[", row[["flow_ymin"]], ",", row[["flow_ymax"]], "]",
      ", stratum=[", row[["stratum_ymin"]], ",", row[["stratum_ymax"]], "]",
      ", delta=", row[["delta"]]
    )
  )
  stop(
    "FLOW ribbons are vertically displaced from their taxonomy strata: ",
    paste(details, collapse = "; "),
    call. = FALSE
  )
}

if (!identical(plot_object$theme$legend.position, "none")) {
  stop("Simple FLOW PNG must not render a legend.", call. = FALSE)
}

message("PASS: simple FLOW ribbons align with their taxonomy strata")
