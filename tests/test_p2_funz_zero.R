source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

expect_true <- function(condition, label) {
  if (!isTRUE(condition)) {
    stop(label, call. = FALSE)
  }
}

failures <- character()
run_case <- function(name, code) {
  tryCatch(
    {
      force(code)
      message("PASS: ", name)
    },
    error = function(error) {
      failures <<- c(failures, paste0(name, ": ", conditionMessage(error)))
      message("FAIL: ", name, " -- ", conditionMessage(error))
    }
  )
}

script_env <- source_without_main("sqm_plots.R")

run_case("FUNZ preserves a zero-TPM sample as an auditable sentinel", {
  orf_long <- tibble::tibble(
    orf_id = c("orf_a", "orf_b"),
    sample = c("S_positive", "S_positive"),
    tpm = c(60, 40),
    ko_id = c("K00001", "K00002"),
    kegg_function = c("Function A", "Function B"),
    ec_codes = c("1.1.1.1", NA_character_)
  )
  selected_samples <- c("S_positive", "S_zero")
  warnings <- character()

  plot_tbl <- withCallingHandlers(
    script_env$build_ko_plot_table(
      orf_long = orf_long,
      selected_samples = selected_samples,
      top_n_ko = 2L,
      ko_lookup = c(K00001 = "Lookup A", K00002 = "Lookup B"),
      pathway_name = "Synthetic pathway"
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )

  required_columns <- c(
    "sample", "ko_id", "tpm", "sample_pathway_percent",
    "denominator", "status", "plotted"
  )
  expect_true(
    all(required_columns %in% colnames(plot_tbl)),
    "The FUNZ table is missing denominator/status/plotted provenance columns"
  )

  positive_tbl <- plot_tbl[
    as.character(plot_tbl$sample) == "S_positive",
    required_columns,
    drop = FALSE
  ]
  expect_true(nrow(positive_tbl) == 2L, "The positive sample lost a KO row")
  expect_true(
    isTRUE(all.equal(sum(positive_tbl$tpm), 100, tolerance = 1e-12)),
    "The positive sample TPM mass changed"
  )
  expect_true(
    isTRUE(all.equal(
      sum(positive_tbl$sample_pathway_percent),
      100,
      tolerance = 1e-6
    )),
    "Positive-sample FUNZ percentages do not sum to 100"
  )
  expect_true(all(positive_tbl$denominator == 100), "The FUNZ denominator is not the expanded KO TPM total")
  expect_true(all(positive_tbl$status == "ok"), "Positive FUNZ rows are not marked ok")
  expect_true(all(positive_tbl$plotted), "Positive FUNZ rows are not marked for plotting")

  zero_tbl <- plot_tbl[
    as.character(plot_tbl$sample) == "S_zero",
    required_columns,
    drop = FALSE
  ]
  expect_true(nrow(zero_tbl) == 1L, "A zero-denominator sample must have exactly one sentinel row")
  expect_true(is.na(zero_tbl$ko_id[[1]]), "The zero-denominator sentinel must use ko_id=NA")
  expect_true(zero_tbl$tpm[[1]] == 0, "The zero-denominator sentinel must use tpm=0")
  expect_true(zero_tbl$denominator[[1]] == 0, "The zero-denominator sentinel must use denominator=0")
  expect_true(
    is.na(zero_tbl$sample_pathway_percent[[1]]),
    "The zero-denominator sentinel percent must be NA"
  )
  expect_true(
    identical(as.character(zero_tbl$status[[1]]), "zero_denominator"),
    "The zero-denominator sentinel has the wrong status"
  )
  expect_true(
    identical(as.logical(zero_tbl$plotted[[1]]), FALSE),
    "The zero-denominator sentinel must not be plotted"
  )
  expect_true(
    any(grepl("Synthetic pathway", warnings, fixed = TRUE)) &&
      any(grepl("S_zero", warnings, fixed = TRUE)),
    "The zero-denominator warning must identify both pathway and sample"
  )

  plot_object <- script_env$make_ko_barplot(
    plot_tbl = plot_tbl,
    pathway_name = "Synthetic pathway",
    selected_samples = selected_samples
  )
  expect_true(inherits(plot_object, "ggplot"), "FUNZ did not return a ggplot object")
  expect_true(
    all(plot_object$data$plotted),
    "The FUNZ plot layer must be controlled only by rows marked plotted=TRUE"
  )
  plot_build <- ggplot2::ggplot_build(plot_object)
  x_labels <- as.character(plot_build$layout$panel_params[[1]]$x$get_labels())
  expect_true(
    identical(x_labels, selected_samples),
    "The zero-denominator sample disappeared from the FUNZ x axis"
  )
})

run_case("all-zero FUNZ writes sentinels and skips only the PNG", {
  orf_ids <- c("orf_zero_a", "orf_zero_b")
  selected_samples <- c("S_zero_a", "S_zero_b")
  pathway_sqm <- list(
    orfs = list(
      table = data.frame(
        `KEGG ID` = c("K00001", "K00002"),
        KEGGFUN = c("Function A [EC:1.1.1.1]", "Function B"),
        KEGGPATH = rep("All zero pathway", 2L),
        row.names = orf_ids,
        check.names = FALSE
      ),
      tax = data.frame(
        phylum = c("Alpha", "Beta"),
        row.names = orf_ids,
        check.names = FALSE
      ),
      tpm = data.frame(
        S_zero_a = c(0, 0),
        S_zero_b = c(0, 0),
        row.names = orf_ids,
        check.names = FALSE
      )
    ),
    misc = list(
      KEGG_names = c(K00001 = "Function A", K00002 = "Function B")
    )
  )

  test_root <- tempfile("p2_funz_all_zero_")
  dir.create(test_root, recursive = TRUE)
  on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)
  warnings <- character()

  result <- withCallingHandlers(
    script_env$run_funz_mode(
      output_dir = test_root,
      manifest_base_dir = test_root,
      output_manifests = list(funz = tibble::tibble()),
      script_name = "sqm_plots.R",
      project_dir = "synthetic_project",
      tax_mode = "prokfilter",
      pathway_name = "All zero pathway",
      pathway_sqm = pathway_sqm,
      selected_samples = selected_samples,
      dimensions = list(`2x2` = c(width = 2, height = 2)),
      plot_dpi = 72,
      top_n_taxa = 2L,
      top_n_ko = 2L,
      pathway_id = "00000"
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )

  pathway_dir <- file.path(
    test_root,
    "funz",
    "pathway",
    "definiti",
    script_env$safe_output_component("All zero pathway", max_length = 28L)
  )
  tsv_path <- file.path(pathway_dir, "barplot_ko_data.tsv")
  expect_true(file.exists(tsv_path), "All-zero FUNZ did not write its sentinel TSV")
  expect_true(
    length(list.files(pathway_dir, pattern = "\\.png$", recursive = TRUE)) == 0L,
    "All-zero FUNZ must skip PNG generation"
  )

  exported <- readr::read_tsv(tsv_path, show_col_types = FALSE)
  required_columns <- c(
    "sample", "ko_id", "tpm", "sample_pathway_percent",
    "denominator", "status", "plotted"
  )
  expect_true(
    all(required_columns %in% colnames(exported)),
    "The all-zero FUNZ TSV is missing sentinel provenance columns"
  )
  expect_true(nrow(exported) == 2L, "All-zero FUNZ must write one sentinel per selected sample")
  expect_true(
    identical(as.character(exported$sample), selected_samples),
    "All-zero FUNZ sentinel sample order changed"
  )
  expect_true(all(is.na(exported$ko_id)), "All-zero FUNZ sentinels must use ko_id=NA")
  expect_true(all(exported$tpm == 0), "All-zero FUNZ sentinels must use tpm=0")
  expect_true(all(exported$denominator == 0), "All-zero FUNZ sentinels must use denominator=0")
  expect_true(
    all(is.na(exported$sample_pathway_percent)),
    "All-zero FUNZ sentinel percentages must be NA"
  )
  expect_true(all(exported$status == "zero_denominator"), "All-zero FUNZ sentinels have the wrong status")
  expect_true(!any(exported$plotted), "All-zero FUNZ sentinels must not be plotted")
  expect_true(
    all(vapply(
      selected_samples,
      function(sample_name) any(grepl(sample_name, warnings, fixed = TRUE)),
      logical(1)
    )) && any(grepl("All zero pathway", warnings, fixed = TRUE)),
    "All-zero warnings must identify the pathway and every affected sample"
  )
  expect_true(
    any(result$funz$output_type == "data_tsv") &&
      !any(result$funz$output_type == "plot_png"),
    "All-zero FUNZ manifest must retain the TSV and omit a PNG target"
  )
})

if (length(failures) > 0L) {
  stop(
    paste(
      "P2 FUNZ zero-denominator regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: P2 FUNZ handles zero-denominator samples explicitly")
