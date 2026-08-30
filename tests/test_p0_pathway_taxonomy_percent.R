source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")

# This fixture is already restricted to one pathway. Its ORF TPM total is
# therefore the required per-sample denominator, not SQM$total_reads.
orf_ids <- c(
  "orf_alpha", "orf_beta", "orf_gamma", "orf_delta", "orf_unclassified"
)

fake_pathway_sqm <- list(
  orfs = list(
    table = data.frame(
      KEGGPATH = rep("Synthetic pathway", length(orf_ids)),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tax = data.frame(
      phylum = c("Alpha", "Beta", "Gamma", "Delta", NA_character_),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tpm = data.frame(
      S_positive = c(45, 25, 15, 5, 10),
      S_zero = rep(0, length(orf_ids)),
      row.names = orf_ids,
      check.names = FALSE
    )
  )
)

selected_samples <- c("S_positive", "S_zero")
captured_warnings <- character()

percent_tbl <- withCallingHandlers(
  script_env$build_pathway_taxonomy_percent_table(
    sqm_object = fake_pathway_sqm,
    rank = "phylum",
    selected_samples = selected_samples,
    top_n_taxa = 2L,
    pathway_name = "Synthetic pathway"
  ),
  warning = function(condition) {
    captured_warnings <<- c(captured_warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
)

required_columns <- c(
  "sample", "taxon", "taxon_tpm", "pathway_tpm", "value",
  "denominator", "status", "plotted"
)
stopifnot(all(required_columns %in% colnames(percent_tbl)))

# Top N is selected once for the pathway. Classified taxa below that cut are
# combined into Other, while missing taxonomy is retained as Unclassified.
positive_tbl <- percent_tbl[
  as.character(percent_tbl$sample) == "S_positive",
  required_columns,
  drop = FALSE
]
positive_tbl$taxon <- as.character(positive_tbl$taxon)
positive_tbl <- positive_tbl[order(positive_tbl$taxon), , drop = FALSE]

expected_values <- c(Alpha = 45, Beta = 25, Other = 20, Unclassified = 10)
stopifnot(identical(positive_tbl$taxon, sort(names(expected_values))))
stopifnot(isTRUE(all.equal(
  unname(positive_tbl$value),
  unname(expected_values[positive_tbl$taxon]),
  tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  unname(positive_tbl$taxon_tpm),
  unname(expected_values[positive_tbl$taxon]),
  tolerance = 1e-12
)))
stopifnot(all(positive_tbl$pathway_tpm == 100))
stopifnot(all(positive_tbl$denominator == 100))
stopifnot(all(positive_tbl$status == "ok"))
stopifnot(all(positive_tbl$plotted))
stopifnot(abs(sum(positive_tbl$value) - 100) <= 1e-6)
stopifnot("Unclassified" %in% positive_tbl$taxon)
stopifnot(!"Gamma" %in% positive_tbl$taxon)
stopifnot(!"Delta" %in% positive_tbl$taxon)

# A zero denominator remains auditable in the TSV table, emits a warning, and
# is not silently turned into zero percent or dropped from the sample set.
zero_tbl <- percent_tbl[
  as.character(percent_tbl$sample) == "S_zero",
  required_columns,
  drop = FALSE
]
stopifnot(nrow(zero_tbl) == 1L)
stopifnot(is.na(zero_tbl$taxon[[1]]))
stopifnot(zero_tbl$taxon_tpm[[1]] == 0)
stopifnot(zero_tbl$pathway_tpm[[1]] == 0)
stopifnot(zero_tbl$denominator[[1]] == 0)
stopifnot(is.na(zero_tbl$value[[1]]))
stopifnot(identical(as.character(zero_tbl$status[[1]]), "zero_denominator"))
stopifnot(identical(as.logical(zero_tbl$plotted[[1]]), FALSE))
stopifnot(any(grepl("S_zero", captured_warnings, fixed = TRUE)))
stopifnot(any(grepl("denominator", captured_warnings, ignore.case = TRUE)))

plot_object <- script_env$make_pathway_taxonomy_percent_plot(
  plot_tbl = percent_tbl,
  pathway_name = "Synthetic pathway",
  rank = "phylum",
  selected_samples = selected_samples
)
stopifnot(inherits(plot_object, "ggplot"))

plot_build <- ggplot2::ggplot_build(plot_object)
x_labels <- as.character(plot_build$layout$panel_params[[1]]$x$get_labels())
stopifnot(identical(x_labels, selected_samples))

# The zero sample may be represented by a zero-height/blank layer to preserve
# its axis label, but it must never receive a positive bar.
zero_x <- match("S_zero", x_labels)
zero_has_positive_bar <- any(vapply(
  plot_build$data,
  function(layer_data) {
    if (!all(c("x", "y") %in% colnames(layer_data))) {
      return(FALSE)
    }
    any(
      as.numeric(layer_data$x) == zero_x &
        !is.na(layer_data$y) &
        layer_data$y > 0
    )
  },
  logical(1)
))
stopifnot(!zero_has_positive_bar)

# Exercise the same table/plot pair through the existing output helpers.
test_root <- tempfile("p0_pathway_taxonomy_percent_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

tsv_path <- file.path(test_root, "taxonomy_percent_data.tsv")
script_env$write_tsv_safe(percent_tbl, tsv_path)
png_paths <- script_env$save_png_dimensions(
  plot_object = plot_object,
  output_dir = test_root,
  file_stem = "taxonomy_percent",
  dimensions = list(`2x2` = c(width = 2, height = 2)),
  dpi = 72
)

roundtrip_tbl <- readr::read_tsv(tsv_path, show_col_types = FALSE)
stopifnot(all(required_columns %in% colnames(roundtrip_tbl)))
stopifnot(nrow(roundtrip_tbl) == nrow(percent_tbl))
stopifnot(length(png_paths) == 1L)
stopifnot(file.exists(unname(png_paths[[1]])))
stopifnot(file.info(unname(png_paths[[1]]))$size > 0)

message("PASS: pathway taxonomy percentages use pathway/sample TPM denominators")
