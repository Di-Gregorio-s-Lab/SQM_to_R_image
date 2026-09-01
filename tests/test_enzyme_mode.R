source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")
expected_ecs <- c(
  "1.14.12.11", "1.14.12.12", "1.14.12.-",
  "3.8.1.2", "3.8.1.3", "1.13.11.-",
  "1.21.99.5", "1.14.13.243", "1.14.13.236", "1.14.13.25",
  "1.14.18.3", "1.14.99.39", "1.14.13.244", "1.14.13.7",
  "1.14.13.227", "1.14.13.230", "1.14.13.69", "2.5.1.18",
  "4.4.1.34", "5.2.1.2"
)

stopifnot(identical(script_env$normalize_enzyme_ecs(NULL), expected_ecs))
stopifnot(identical(script_env$normalize_enzyme_plot_types(NULL), c("bar", "line")))

invalid_ec_error <- tryCatch(
  {
    script_env$normalize_enzyme_ecs("not-an-ec")
    NULL
  },
  error = function(error) error$message
)
stopifnot(grepl("enzyme_ecs", invalid_ec_error, fixed = TRUE))

invalid_plot_type_error <- tryCatch(
  {
    script_env$normalize_enzyme_plot_types("boxplot")
    NULL
  },
  error = function(error) error$message
)
stopifnot(grepl("enzyme_plot_types", invalid_plot_type_error, fixed = TRUE))

fake_sqm <- list(
  orfs = list(
    table = data.frame(
      "KEGG ID" = c("K00001;K00002", "K00003", "K00004", "K00005"),
      KEGGFUN = c(
        "Function A [EC:1.14.12.11]",
        "Function B [EC:1.14.12.-]",
        "Function C [EC:1.14.12.11 1.14.12.12]",
        "Function D [EC:1.14.12.99]"
      ),
      check.names = FALSE,
      row.names = c("orf_a", "orf_b", "orf_c", "orf_d")
    ),
    tpm = data.frame(
      S0 = c(10, 5, 3, 100),
      S1 = c(20, 0, 4, 100),
      row.names = c("orf_a", "orf_b", "orf_c", "orf_d")
    )
  )
)

enzyme_tbl <- script_env$build_enzyme_plot_table(fake_sqm, c("S0", "S1"), expected_ecs)
value_for <- function(sample_name, ec_code) {
  enzyme_tbl$tpm[
    as.character(enzyme_tbl$sample) == sample_name &
      as.character(enzyme_tbl$ec_code) == ec_code
  ]
}

stopifnot(identical(value_for("S0", "1.14.12.11"), 13))
stopifnot(identical(value_for("S1", "1.14.12.11"), 24))
stopifnot(identical(value_for("S0", "1.14.12.-"), 5))
stopifnot(identical(value_for("S1", "1.14.12.-"), 0))
stopifnot(identical(value_for("S0", "1.14.12.12"), 3))
stopifnot(identical(value_for("S1", "1.14.12.12"), 4))
stopifnot(identical(value_for("S0", "3.8.1.2"), 0))

test_root <- tempfile("enzyme_mode_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

result <- suppressWarnings(script_env$run_enzyme_mode(
  sqm_object = fake_sqm,
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
  enzyme_ecs = expected_ecs,
  enzyme_plot_types = c("bar", "line")
))

combined_dir <- file.path(test_root, "funz", "enzimi", "insieme")
separate_dir <- file.path(test_root, "funz", "enzimi", "separato", "1.14.12.11")
stopifnot(file.exists(file.path(combined_dir, "enzimi_data.tsv")))
stopifnot(file.exists(file.path(combined_dir, "barplot_enzimi_2x2.png")))
stopifnot(file.exists(file.path(combined_dir, "lineplot_enzimi_2x2.png")))
stopifnot(file.exists(file.path(separate_dir, "enzima_data.tsv")))
stopifnot(file.exists(file.path(separate_dir, "barplot_enzima_2x2.png")))
stopifnot(file.exists(file.path(separate_dir, "lineplot_enzima_2x2.png")))
positive_ec_count <- length(unique(as.character(enzyme_tbl$ec_code[enzyme_tbl$plotted])))
stopifnot(nrow(result$funz) == 3L + length(expected_ecs) + positive_ec_count * 2L)
stopifnot(all(c("ec_code", "output_scope") %in% colnames(result$funz)))
stopifnot(identical(
  sort(unique(result$funz$output_scope)),
  c("enzyme_insieme", "enzyme_separato")
))
stopifnot(all(result$funz$ec_code[result$funz$output_scope == "enzyme_separato"] %in% expected_ecs))

fake_pathway_sqm <- fake_sqm
fake_pathway_sqm$orfs$table$KEGGPATH <- "Test pathway"
fake_pathway_sqm$orfs$tax <- data.frame(
  superkingdom = rep("Bacteria", 4),
  phylum = rep("Proteobacteria", 4),
  class = rep("Gammaproteobacteria", 4),
  order = rep("Pseudomonadales", 4),
  family = rep("Pseudomonadaceae", 4),
  genus = rep("Pseudomonas", 4),
  species = rep("Pseudomonas sp.", 4),
  row.names = rownames(fake_sqm$orfs$table),
  check.names = FALSE
)
fake_pathway_sqm$misc <- list(
  KEGG_names = c(
    K00001 = "Function A", K00002 = "Function B", K00003 = "Function C",
    K00004 = "Function D", K00005 = "Function E"
  )
)
pathway_result <- script_env$run_funz_mode(
  output_dir = test_root,
  manifest_base_dir = test_root,
  output_manifests = list(funz = tibble::tibble()),
  script_name = "sqm_plots.R",
  project_dir = "project",
  tax_mode = "prokfilter",
  pathway_name = "Test pathway",
  pathway_sqm = fake_pathway_sqm,
  selected_samples = c("S0", "S1"),
  dimensions = list("2x2" = c(width = 2, height = 2)),
  plot_dpi = 72,
  top_n_taxa = 15L,
  top_n_ko = 20L,
  pathway_id = "00000"
)
stopifnot(file.exists(file.path(
  test_root,
  "funz",
  "pathway",
  "definiti",
  script_env$safe_output_component("Test pathway", max_length = 28L),
  "barplot_ko_data.tsv"
)))
stopifnot(!dir.exists(file.path(test_root, "funz", "Test_pathway")))
stopifnot(nrow(pathway_result$funz) == 2L)

message("PASS: Enzyme mode aggregates exact EC matches and writes both plot layouts")
