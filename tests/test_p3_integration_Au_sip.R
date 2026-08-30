source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)
  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")
project_dir <- "in/Au_sip"
if (!dir.exists(project_dir)) {
  stop("Integration fixture not found: ", project_dir, call. = FALSE)
}

# Load the real SQM exactly once; the compatibility warning must remain visible.
sqm <- SQMtools::loadSQM(
  project_path = script_env$normalize_sqm_project_dir(project_dir),
  tax_mode = "prokfilter",
  trusted_functions_only = FALSE,
  load_sequences = FALSE
)
script_env$validate_sqm_object(sqm)

resolved <- script_env$resolve_pathways(sqm, "00633")[[1L]]
stopifnot(identical(resolved$pathway_id, "00633"))
pathway_sqm <- script_env$subset_pathway(sqm, resolved$canonical_pathway_name)
sample_name <- "S13_1_8"
orf_result <- script_env$build_orf_long_result(pathway_sqm, sample_name)
orf_long <- orf_result$data
stopifnot(nrow(orf_long) > 0L)

sample_orfs <- orf_long |>
  dplyr::filter(.data$sample == sample_name)
pathway_sample_tpm <- sum(sample_orfs$tpm)
ko_reference <- sample_orfs |>
  dplyr::group_by(.data$ko_id) |>
  dplyr::summarise(ko_sample_tpm = sum(.data$tpm), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(.data$ko_sample_tpm), .data$ko_id)
ko_id <- ko_reference$ko_id[[1L]]
expected_ko_tpm <- ko_reference$ko_sample_tpm[[1L]]
expected_contribution <- expected_ko_tpm / pathway_sample_tpm * 100

pie_table <- script_env$build_pie_chart_table(
  orf_long = orf_long,
  sample_name = sample_name,
  ko_id_filter = ko_id,
  rank_name = "phylum",
  top_n_taxa = 3L,
  pathway_name = resolved$canonical_pathway_name,
  pathway_id = resolved$pathway_id,
  pathway_selection = "defined",
  pathway_sample_tpm = pathway_sample_tpm
)

required_columns <- c(
  "pathway", "pathway_id", "pathway_selection", "rank", "ko_name",
  "ko_sample_tpm", "pathway_sample_tpm", "ko_pathway_percent", "taxon_order"
)
stopifnot(all(required_columns %in% colnames(pie_table)))
stopifnot(isTRUE(all.equal(sum(pie_table$tpm), expected_ko_tpm, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(unique(pie_table$ko_sample_tpm), expected_ko_tpm, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(unique(pie_table$pathway_sample_tpm), pathway_sample_tpm, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(unique(pie_table$ko_pathway_percent), expected_contribution, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(sum(pie_table$pct), 1, tolerance = 1e-10)))
stopifnot(identical(as.integer(pie_table$taxon_order), seq_len(nrow(pie_table))))

roundtrip_root <- tempfile("p3_integration_pie_")
dir.create(roundtrip_root, recursive = TRUE)
on.exit(unlink(roundtrip_root, recursive = TRUE, force = TRUE), add = TRUE)
roundtrip_path <- file.path(roundtrip_root, "pie.tsv")
readr::write_tsv(pie_table, roundtrip_path, na = "NA")
roundtrip <- readr::read_tsv(roundtrip_path, show_col_types = FALSE, na = "NA")
pie_plot <- script_env$make_pie_plot(roundtrip)
stopifnot(isTRUE(all.equal(sum(pie_plot$data$tpm), expected_ko_tpm, tolerance = 1e-10)))
stopifnot(grepl(resolved$canonical_pathway_name, pie_plot$labels$title, fixed = TRUE))
stopifnot(grepl(ko_id, pie_plot$labels$title, fixed = TRUE))

message(
  "PASS: P3 integration checks completed | SQMtools=",
  as.character(utils::packageVersion("SQMtools")),
  " | pathway=00633 | sample=", sample_name,
  " | KO=", ko_id,
  " | KO_TPM=", format(expected_ko_tpm, digits = 12),
  " | pathway_TPM=", format(pathway_sample_tpm, digits = 12)
)
