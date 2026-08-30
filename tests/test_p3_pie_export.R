source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

expect_identical <- function(actual, expected, label) {
  if (!identical(actual, expected)) {
    stop(
      label,
      "; expected ", paste(expected, collapse = ", "),
      ", observed ", paste(actual, collapse = ", "),
      call. = FALSE
    )
  }
}

expect_number <- function(actual, expected, label, tolerance = 1e-10) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance))) {
    stop(label, "; expected ", expected, ", observed ", actual, call. = FALSE)
  }
}

script_env <- source_without_main("sqm_plots.R")

# The two target-ORF TPM values are 30 and 10; a third ORF contributes 60 TPM
# to the same pathway but belongs to another KO. These worked-example literals
# independently fix the KO denominator (40) and pathway denominator (100).
orf_fixture <- tibble::tibble(
  orf_id = c("orf_alpha", "orf_beta", "orf_other_ko"),
  sample = rep("S_pie_export", 3L),
  tpm = c(30, 10, 60),
  ko_id = c("K00001", "K00001", "K99999"),
  kegg_function = c(
    "Synthetic KO label",
    "Synthetic KO label",
    "Other synthetic KO"
  ),
  ec_codes = c("1.1.1.-;1.1.1.1", "1.1.1.-;1.1.1.1", NA_character_),
  phylum = c("Alpha", "Beta", "Gamma")
)

pie_table <- script_env$build_pie_chart_table(
  orf_long = orf_fixture,
  sample_name = "S_pie_export",
  ko_id_filter = "K00001",
  rank_name = "phylum",
  top_n_taxa = 2L,
  pathway_name = "Synthetic pathway",
  pathway_id = "12345",
  pathway_selection = "defined",
  pathway_sample_tpm = 100
)

legacy_columns <- c(
  "taxon_rank", "tpm", "sample", "ko_id", "ec_codes",
  "total_tpm", "pct", "label"
)
p3_columns <- c(
  "pathway", "pathway_id", "pathway_selection", "rank", "ko_name",
  "ko_sample_tpm", "pathway_sample_tpm", "ko_pathway_percent", "taxon_order"
)
expect_identical(
  colnames(pie_table),
  c(legacy_columns, p3_columns),
  "PIE TSV columns must remain append-only"
)

expect_identical(as.character(pie_table$taxon_rank), c("Alpha", "Beta"), "PIE taxa changed")
expect_identical(as.numeric(pie_table$tpm), c(30, 10), "PIE taxon TPM changed")
expect_identical(as.numeric(pie_table$pct), c(0.75, 0.25), "PIE KO fractions changed")
expect_identical(as.character(pie_table$pathway), rep("Synthetic pathway", 2L), "PIE pathway missing")
expect_identical(as.character(pie_table$pathway_id), rep("12345", 2L), "PIE pathway ID missing")
expect_identical(as.character(pie_table$pathway_selection), rep("defined", 2L), "PIE selection missing")
expect_identical(as.character(pie_table$rank), rep("phylum", 2L), "PIE rank missing")
expect_identical(as.character(pie_table$ko_name), rep("Synthetic KO label", 2L), "PIE KO name missing")
expect_identical(as.numeric(pie_table$ko_sample_tpm), c(40, 40), "PIE KO denominator is not literal 40")
expect_identical(as.numeric(pie_table$pathway_sample_tpm), c(100, 100), "PIE pathway denominator is not literal 100")
expect_identical(as.numeric(pie_table$ko_pathway_percent), c(40, 40), "PIE KO contribution is not literal 40 percent")
expect_identical(as.integer(pie_table$taxon_order), c(1L, 2L), "PIE taxon order is not explicit")

roundtrip_root <- tempfile("p3_pie_roundtrip_")
dir.create(roundtrip_root, recursive = TRUE)
on.exit(unlink(roundtrip_root, recursive = TRUE, force = TRUE), add = TRUE)
roundtrip_path <- file.path(roundtrip_root, "pie.tsv")
readr::write_tsv(pie_table, roundtrip_path, na = "NA")
roundtrip_table <- readr::read_tsv(roundtrip_path, show_col_types = FALSE, na = "NA")

# This call intentionally supplies no parallel metadata. All plot semantics must
# survive the TSV boundary and be reconstructed from roundtrip_table alone.
pie_plot <- script_env$make_pie_plot(roundtrip_table)

expect_identical(
  names(formals(script_env$make_pie_plot)),
  "plot_data",
  "make_pie_plot must expose the exported table as its only input"
)
expect_identical(
  pie_plot$labels$title,
  "Pathway: Synthetic pathway - KO K00001",
  "PIE title is not derived from the exported table"
)
expect_identical(
  pie_plot$labels$subtitle,
  "Sample: S_pie_export | Rank: phylum | total TPM = 40.00 | pathway contribution: 40.00%",
  "PIE subtitle is not derived from the exported table"
)
expect_identical(
  pie_plot$labels$caption,
  "KO name: Synthetic KO label | EC: 1.1.1.-;1.1.1.1",
  "PIE caption is not derived from the exported table"
)
expect_identical(pie_plot$labels$fill, "phylum", "PIE legend title is not derived from rank")
expect_identical(
  levels(pie_plot$data$taxon_rank),
  c("Alpha", "Beta"),
  "PIE factor order did not survive TSV round-trip"
)
expect_number(sum(pie_plot$data$tpm), 40, "PIE plot changed the exported TPM mass")

message("PASS: P3 PIE TSV is complete and fully controls the plot")
