source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)
  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")
failures <- character()
run_case <- function(label, code) {
  tryCatch(
    {
      force(code)
      message("PASS: ", label)
    },
    error = function(error) {
      failures <<- c(failures, paste0(label, ": ", conditionMessage(error)))
      message("FAIL: ", label, " -- ", conditionMessage(error))
    }
  )
}

test_root <- tempfile("audit_remediation_outputs_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

run_case("PIE exports every positive KO and marks top_n_ko not applicable", {
  pie_sqm <- list(
    orfs = list(
      table = data.frame(
        "KEGG ID" = c("K00001", "K00002"),
        KEGGFUN = c("Function A [EC:1.1.1.1]", "Function B [EC:1.1.1.2]"),
        KEGGPATH = rep("Synthetic pathway", 2L),
        row.names = c("orf_a", "orf_b"),
        check.names = FALSE
      ),
      tax = data.frame(
        phylum = c("Alpha", "Beta"),
        row.names = c("orf_a", "orf_b"),
        check.names = FALSE
      ),
      tpm = data.frame(
        S0 = c(10, 20),
        row.names = c("orf_a", "orf_b"),
        check.names = FALSE
      )
    ),
    misc = list(KEGG_names = c(K00001 = "Function A", K00002 = "Function B"))
  )
  result <- script_env$run_pie_mode(
    output_dir = test_root,
    manifest_base_dir = test_root,
    output_manifests = list(pie = tibble::tibble()),
    script_name = "sqm_plots.R",
    project_dir = "project",
    tax_mode = "prokfilter",
    pathway_name = "Synthetic pathway",
    pathway_sqm = pie_sqm,
    selected_samples = "S0",
    taxonomy_ranks = "phylum",
    dimensions = list("2x2" = c(width = 2, height = 2)),
    plot_dpi = 72,
    top_n_taxa = 5L,
    top_n_ko = 1L,
    pathway_id = "12345"
  )
  stopifnot(setequal(unique(result$pie$ko_id), c("K00001", "K00002")))
  stopifnot(all(is.na(result$pie$top_n_ko)))
  stopifnot(identical(unique(result$pie$ko_selection_policy), "all_positive_ko"))
})

run_case("Pathview isolates each export and manifests exact source data", {
  pathview_sqm <- list(
    functions = list(
      KEGG = list(
        tpm = data.frame(
          S0 = c(10, 20),
          S1 = c(30, 40),
          row.names = c("K00001", "K00002"),
          check.names = FALSE
        )
      )
    )
  )
  stale_dir <- file.path(
    test_root, "pathview", "definiti", "insieme", "Nitrogen_metabolism"
  )
  dir.create(stale_dir, recursive = TRUE)
  writeLines("stale", file.path(stale_dir, "stale_old.png"))

  calls <- list()
  fake_export <- function(
      SQM, pathway_id, count, samples, split_samples, log_scale, output_dir, output_suffix) {
    calls[[length(calls) + 1L]] <<- list(
      samples = samples,
      split_samples = split_samples,
      log_scale = log_scale,
      output_dir = output_dir
    )
    writeLines("current", file.path(output_dir, "map.png"))
  }
  result <- script_env$run_pathview_mode(
    sqm_object = pathview_sqm,
    output_dir = test_root,
    manifest_base_dir = test_root,
    output_manifests = list(pathview = tibble::tibble()),
    script_name = "sqm_plots.R",
    project_dir = "project",
    tax_mode = "prokfilter",
    pathway_name = "Nitrogen metabolism",
    pathway_id = "00910",
    selected_samples = c("S0", "S1"),
    top_n_taxa = 15L,
    top_n_ko = 5L,
    pathview_sample_modes = c("insieme", "separato"),
    export_pathway_fn = fake_export
  )

  stopifnot(length(calls) == 3L)
  stopifnot(all(!vapply(calls, `[[`, logical(1), "split_samples")))
  stopifnot(all(!vapply(calls, `[[`, logical(1), "log_scale")))
  stopifnot(identical(calls[[1L]]$samples, c("S0", "S1")))
  stopifnot(identical(calls[[2L]]$samples, "S0"))
  stopifnot(identical(calls[[3L]]$samples, "S1"))
  stopifnot(!any(grepl("stale_old", result$pathview$output_file, fixed = TRUE)))
  input_rows <- result$pathview$output_type == "pathview_input_all_ko_complete_matrix_tsv"
  stopifnot(sum(input_rows) == 3L)
  stopifnot(sum(result$pathview$output_type == "pathview_render_config_tsv") == 3L)
  plot_rows <- result$pathview$output_type == "pathview_file"
  stopifnot(all(!is.na(result$pathview$source_data_file[plot_rows])))
  stopifnot(all(is.na(result$pathview$top_n_ko)))
  stopifnot(all(result$pathview$ko_selection_policy[input_rows] == "all_ko_complete_matrix"))
  stopifnot(all(result$pathview$ko_selection_policy[plot_rows] == "pathview_native_mapping"))
  stopifnot(all(file.exists(file.path(test_root, result$pathview$output_file))))

  separate_rows <- result$pathview$output_scope == "pathway_defined_separato" & plot_rows
  stopifnot(setequal(result$pathview$samples[separate_rows], c("S0", "S1")))
  source_path <- file.path(
    test_root,
    result$pathview$output_file[input_rows][[1L]]
  )
  source_data <- readr::read_tsv(source_path, show_col_types = FALSE)
  stopifnot(identical(source_data$ko_id, c("K00001", "K00002")))
  stopifnot(identical(source_data$S0, c(10, 20)))
})

run_case("FLOW HTML is self-contained", {
  html_path <- file.path(test_root, "selfcontained", "widget.html")
  widget <- plotly::plot_ly(x = 1:2, y = c(2, 1), type = "scatter", mode = "lines")
  script_env$save_html_widget(widget, html_path)
  stopifnot(file.exists(html_path), file.info(html_path)$size > 0)
  stopifnot(!dir.exists(file.path(dirname(html_path), "widget_files")))

  protected_html_path <- file.path(test_root, "protected", "widget.html")
  protected_dependency_dir <- file.path(dirname(protected_html_path), "widget_files")
  dir.create(protected_dependency_dir, recursive = TRUE)
  protected_file <- file.path(protected_dependency_dir, "existing-output.txt")
  writeLines("preserve", protected_file)
  protected_error <- tryCatch(
    {
      script_env$save_html_widget(widget, protected_html_path)
      NA_character_
    },
    error = function(error) conditionMessage(error)
  )
  stopifnot(!is.na(protected_error))
  stopifnot(file.exists(protected_file))
})

run_case("zero-signal enzymes have status TSV but no PNG", {
  enzyme_sqm <- list(
    orfs = list(
      table = data.frame(
        "KEGG ID" = c("K00001", "K00002"),
        KEGGFUN = c("Observed [EC:1.1.1.1]", "Other [EC:9.9.9.9]"),
        row.names = c("orf_a", "orf_b"),
        check.names = FALSE
      ),
      tpm = data.frame(
        S0 = c(10, 20),
        S1 = c(30, 40),
        row.names = c("orf_a", "orf_b"),
        check.names = FALSE
      )
    ),
    functions = list(
      KEGG = list(
        tpm = data.frame(
          S0 = c(10, 20),
          S1 = c(30, 40),
          row.names = c("K00001", "K00002"),
          check.names = FALSE
        )
      )
    ),
    misc = list(
      KEGG_names = c(
        K00001 = "Observed [EC:1.1.1.1]",
        K00002 = "Other [EC:9.9.9.9]"
      )
    )
  )
  result <- suppressWarnings(script_env$run_enzyme_mode(
    sqm_object = enzyme_sqm,
    output_dir = test_root,
    manifest_base_dir = test_root,
    output_manifests = list(funz = tibble::tibble()),
    script_name = "sqm_plots.R",
    project_dir = "project",
    tax_mode = "prokfilter",
    selected_samples = c("S0", "S1"),
    dimensions = list("2x2" = c(width = 2, height = 2)),
    plot_dpi = 72,
    top_n_taxa = 15L,
    top_n_ko = 20L,
    enzyme_ecs = c("1.1.1.1", "2.2.2.2"),
    enzyme_plot_types = c("bar", "line"),
    sample_order_basis = "sqm_column_order"
  ))

  zero_dir <- file.path(test_root, "funz", "enzimi", "separato", "2.2.2.2")
  zero_tsv <- file.path(zero_dir, "enzima_data.tsv")
  stopifnot(file.exists(zero_tsv))
  zero_data <- readr::read_tsv(zero_tsv, show_col_types = FALSE)
  stopifnot(all(zero_data$status == "no_positive_tpm"), !any(zero_data$plotted))
  stopifnot(length(list.files(zero_dir, pattern = "\\.png$")) == 0L)

  png_rows <- result$funz$format == "png"
  stopifnot(all(result$funz$width[png_rows] == 2))
  stopifnot(all(result$funz$height[png_rows] == 2))
  stopifnot(all(result$funz$sample_order_basis == "sqm_column_order"))
  combined_data <- readr::read_tsv(
    file.path(test_root, "funz", "enzimi", "insieme", "enzimi_data.tsv"),
    show_col_types = FALSE
  )
  stopifnot(identical(unique(combined_data$sample_order), c(1, 2)))
})

if (length(failures) > 0L) {
  stop(
    paste(
      "Audit remediation output regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: audit remediation output invariants are enforced")
