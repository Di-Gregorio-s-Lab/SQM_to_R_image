source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

script_env <- t1_source_sqm_plots_without_main()
orf_ids <- c("orf_single", "orf_multi", "orf_zero_ko", "orf_outside")
sqm <- list(
  orfs = list(
    table = data.frame(
      "KEGG ID" = c("K00001", "K00001;K99999*", "K00003", "K77777"),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tax = data.frame(
      superkingdom = rep("Bacteria", 4L),
      phylum = c("Alpha", "Beta", "Gamma", "Outside"),
      class = rep("Synthetic class", 4L),
      order = rep("Synthetic order", 4L),
      family = rep("Synthetic family", 4L),
      genus = rep("Synthetic genus", 4L),
      species = rep("Synthetic species", 4L),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tpm = data.frame(
      S1 = c(10, 20, 4, 1000),
      row.names = orf_ids,
      check.names = FALSE
    )
  ),
  functions = list(
    KEGG = list(
      tpm = data.frame(
        S1 = c(100, 0, 40, 1000),
        row.names = c("K00001", "K00003", "K99999", "K77777"),
        check.names = FALSE
      )
    )
  ),
  misc = list(
    KEGG_names = c(
      K00001 = "Shared enzyme [EC:1.1.1.1 2.2.2.2]",
      K00003 = "Zero pathway KO [EC:1.1.1.1]",
      K99999 = "Companion KO [EC:1.1.1.1]",
      K77777 = "Outside KO [EC:7.7.7.7]"
    )
  )
)

node_data <- list(
  kegg.names = list(`1` = c("K00001", "K00003"), `2` = "K00001", `3` = "C00001"),
  type = c(`1` = "ortholog", `2` = "ortholog", `3` = "compound")
)
t1_expect_identical(
  script_env$extract_pathway_ko_ids(node_data),
  c("K00001", "K00003"),
  "Repeated KGML nodes must not duplicate pathway KOs"
)

result <- script_env$build_pathway_ko_result(
  context_sqm = sqm,
  selected_samples = "S1",
  pathway_ko_ids = c("K00003", "K00001", "K00001")
)
t1_expect_identical(
  result$orf_ids,
  c("orf_single", "orf_multi", "orf_zero_ko"),
  "KGML KO membership must select unique ORFs without KEGGPATH"
)
t1_expect_equal(
  result$audit$raw_conservation_max_abs_error,
  0,
  "Multi-KO raw allocation conservation"
)
t1_expect_equal(
  result$audit$official_margin_max_abs_error,
  0,
  "Official KO margin conservation"
)

k1 <- result$data |>
  dplyr::filter(.data$ko_id == "K00001") |>
  dplyr::arrange(.data$phylum)
t1_expect_identical(
  k1$phylum,
  c("Alpha", "Beta"),
  "K00001 taxonomic allocation keys"
)
t1_expect_equal(
  k1$tpm,
  c(50, 50),
  "Multi-KO split must occur before pathway filtering"
)
t1_expect_equal(sum(k1$tpm), 100, "K00001 official margin")
t1_expect_true(
  !any(result$data$ko_id == "K00003"),
  "A zero official KO produced a positive allocation"
)
t1_expect_equal(
  result$totals$tpm[result$totals$ko_id == "K00003"],
  0,
  "A zero official KO disappeared from canonical totals"
)

funz_output <- script_env$build_ko_plot_table(
  result$totals,
  "S1",
  10L,
  script_env$get_ko_name_lookup(sqm),
  "Synthetic pathway"
)
t1_expect_identical(
  colnames(funz_output),
  c(
    "sample", "ko_id", "tpm", "kegg_function",
    "sample_pathway_total_tpm", "sample_pathway_percent", "denominator",
    "status", "plotted", "ec_codes"
  ),
  "FUNZ TSV columns changed"
)
flow_rank <- script_env$build_flow_table_for_rank(
  result$data,
  "phylum",
  "S1",
  10L,
  10L,
  script_env$get_ko_name_lookup(sqm)
)
flow_output <- script_env$build_flow_table_for_sample(
  flow_rank,
  "Synthetic pathway",
  "phylum",
  "S1"
)
t1_expect_identical(
  colnames(flow_output),
  c(
    "pathway", "rank", "sample", "taxon", "KO", "KO_name", "ec_codes",
    "TPM", "taxon_percent", "KO_percent", "flow_percent"
  ),
  "FLOW TSV columns changed"
)
pie_output <- script_env$build_pie_chart_table(
  result$data,
  "S1",
  "K00001",
  "phylum",
  10L
)
t1_expect_identical(
  colnames(pie_output),
  c(
    "taxon_rank", "tpm", "sample", "ko_id", "ec_codes", "total_tpm",
    "pct", "label", "pathway", "pathway_id", "pathway_selection", "rank",
    "ko_name", "ko_sample_tpm", "pathway_sample_tpm",
    "ko_pathway_percent", "taxon_order"
  ),
  "PIE TSV columns changed"
)

positive_without_orf <- sqm
positive_without_orf$functions$KEGG$tpm <- rbind(
  positive_without_orf$functions$KEGG$tpm,
  K00004 = c(S1 = 5)
)
positive_without_orf$misc$KEGG_names <- c(
  positive_without_orf$misc$KEGG_names,
  K00004 = "No supporting ORF"
)
allocation_error <- tryCatch(
  {
    script_env$build_pathway_ko_result(
      positive_without_orf,
      "S1",
      "K00004"
    )
    NA_character_
  },
  error = function(error) conditionMessage(error)
)
t1_expect_true(
  !is.na(allocation_error) && grepl("without ORFs", allocation_error, fixed = TRUE),
  "Positive official KO without allocatable ORFs must fail explicitly"
)

enzyme <- script_env$build_enzyme_plot_table(
  sqm,
  "S1",
  c("1.1.1.1", "2.2.2.2", "9.9.9.9")
)
enzyme_value <- function(ec) {
  enzyme$tpm[as.character(enzyme$ec_code) == ec]
}
t1_expect_equal(enzyme_value("1.1.1.1"), 140, "Multiple KOs per EC")
t1_expect_equal(enzyme_value("2.2.2.2"), 100, "Multiple ECs per KO")
t1_expect_equal(enzyme_value("9.9.9.9"), 0, "Missing EC zero state")
t1_expect_identical(
  as.character(enzyme$status[as.character(enzyme$ec_code) == "9.9.9.9"]),
  "no_positive_tpm",
  "Missing EC status"
)

catalog <- script_env$parse_kegg_pathway_catalog(c(
  "path:ko12345\tExample pathway - Reference pathway",
  "ko54321\tSecond   pathway"
))
t1_expect_identical(
  script_env$resolve_pathway_id_from_catalog(" example  pathway ", catalog),
  "12345",
  "KEGG catalog normalization"
)
t1_expect_identical(
  script_env$resolve_kegg_pathway_id("Nitrogen metabolism"),
  "00910",
  "Curated pathway IDs must take precedence"
)
catalog_loads <- 0L
catalog_loader <- function() {
  catalog_loads <<- catalog_loads + 1L
  catalog
}
t1_expect_identical(
  script_env$resolve_kegg_pathway_id("Example pathway", catalog_loader),
  "12345",
  "Live catalog pathway resolution"
)
t1_expect_identical(
  script_env$resolve_kegg_pathway_id("Second pathway", catalog_loader),
  "54321",
  "Cached live catalog pathway resolution"
)
t1_expect_identical(catalog_loads, 1L, "KEGG catalog must download once per run")

ambiguous_catalog <- dplyr::bind_rows(
  catalog,
  tibble::tibble(
    pathway_id = "99999",
    pathway_name = "Example pathway",
    normalized_name = "example pathway"
  )
)
ambiguous_error <- tryCatch(
  script_env$resolve_pathway_id_from_catalog("Example pathway", ambiguous_catalog),
  error = function(error) conditionMessage(error)
)
t1_expect_true(
  grepl("ambiguously", ambiguous_error, fixed = TRUE),
  "Ambiguous KEGG catalog mapping must fail"
)
missing_error <- tryCatch(
  script_env$resolve_pathway_id_from_catalog("Absent pathway", catalog),
  error = function(error) conditionMessage(error)
)
t1_expect_true(
  grepl("could not be resolved", missing_error, fixed = TRUE),
  "Absent KEGG catalog mapping must fail"
)

unreachable <- paste0("file:///", tempfile("missing_kegg_catalog_"))
download_error <- tryCatch(
  script_env$download_kegg_pathway_catalog(unreachable),
  error = function(error) conditionMessage(error)
)
t1_expect_true(
  grepl("Unable to download", download_error, fixed = TRUE),
  "Unreachable KEGG catalog must fail explicitly"
)

message("PASS: canonical KEGG engine conserves allocations and resolves metadata")
