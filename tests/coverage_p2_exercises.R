# This file is sourced by covr::file_coverage() after sqm_plots.R has been
# instrumented in the same environment. It exercises the P2 public helpers
# without sourcing a second, uninstrumented copy of the production script.

taxa_fixture <- tibble::tibble(
  taxon = c("Alpha", "Unclassified", "Beta", "Gamma", NA_character_),
  tpm = c(50, 40, 30, 20, 0)
)
top_taxa <- select_top_classified_taxa(taxa_fixture, "taxon", "tpm", 2L)
collapse_taxa_preserving_unclassified(taxa_fixture$taxon, top_taxa)
try(select_top_classified_taxa(taxa_fixture, "missing", "tpm", 2L), silent = TRUE)
try(
  select_top_classified_taxa(
    tibble::tibble(taxon = "Other", tpm = 1),
    "taxon",
    "tpm",
    1L
  ),
  silent = TRUE
)
try(collapse_taxa_preserving_unclassified(c("Alpha", "Other"), "Alpha"), silent = TRUE)

parse_positive_integer_arg("0007", "top_n_ko")
parse_positive_integer_arg(as.character(.Machine$integer.max), "top_n_ko")
for (value in list("1.5", "0", "2147483648", NA_character_, 1L, character())) {
  try(parse_positive_integer_arg(value, "top_n_ko"), silent = TRUE)
}

validate_tax_mode("prokfilter")
try(validate_tax_mode("invalid"), silent = TRUE)
validate_tpm_matrix(data.frame(S = c(0, 1)), "valid TPM")
try(validate_tpm_matrix(data.frame(S = "invalid"), "invalid TPM"), silent = TRUE)
try(validate_tpm_matrix(data.frame(S = Inf), "invalid TPM"), silent = TRUE)

orf_ids <- c("orf_none", "orf_mono", "orf_multi")
pathway_sqm <- list(
  orfs = list(
    table = data.frame(
      `KEGG ID` = c(NA_character_, "K00001", "K00002; K00003*"),
      KEGGFUN = c(
        "No annotation",
        "Function one [EC:2.2.2.2 1.1.1.-]",
        NA_character_
      ),
      KEGGPATH = c(NA_character_, "Synthetic pathway", "Synthetic pathway"),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tax = data.frame(
      phylum = c("Unclassified", "Alpha", "Beta"),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tpm = data.frame(
      S_positive = c(4, 10, 20),
      S_zero = c(0, 0, 0),
      row.names = orf_ids,
      check.names = FALSE
    )
  ),
  misc = list(
    KEGG_names = c(
      K00001 = "Lookup one",
      K00002 = "Lookup two",
      K00003 = "Lookup three"
    )
  )
)

validate_sqm_tpm_inputs(pathway_sqm, c("S_positive", "S_zero"))
try(validate_sqm_tpm_inputs(pathway_sqm, "missing"), silent = TRUE)
pathway_sqm$functions <- list(KEGG = list(tpm = data.frame(
  S_positive = c(10, 20),
  S_zero = c(0, 0),
  row.names = c("K00001", "K00002")
)))
validate_sqm_tpm_inputs(
  pathway_sqm,
  c("S_positive", "S_zero"),
  require_kegg_tpm = TRUE
)
try(
  validate_sqm_tpm_inputs(
    within(pathway_sqm, rm(functions)),
    "S_positive",
    require_kegg_tpm = TRUE
  ),
  silent = TRUE
)

enzyme_sqm <- list(
  orfs = list(
    table = data.frame(
      `KEGG ID` = c("K00001", "K00002"),
      KEGGFUN = c("Observed [EC:1.1.1.1]", "Other [EC:9.9.9.9]"),
      row.names = c("orf_a", "orf_b"),
      check.names = FALSE
    ),
    tpm = data.frame(
      S_positive = c(10, 20),
      S_zero = c(5, 0),
      row.names = c("orf_a", "orf_b")
    )
  )
)
enzyme_table <- build_enzyme_plot_table(
  enzyme_sqm,
  c("S_positive", "S_zero"),
  c("1.1.1.1", "2.2.2.2")
)
make_enzyme_barplot(enzyme_table, "Synthetic enzymes")
make_enzyme_lineplot(enzyme_table, "Synthetic enzymes")
duplicate_enzyme_sqm <- enzyme_sqm
attr(duplicate_enzyme_sqm$orfs$table, "row.names") <- c("orf_a", "orf_a")
try(
  build_enzyme_plot_table(duplicate_enzyme_sqm, "S_positive", "1.1.1.1"),
  silent = TRUE
)
mismatched_enzyme_sqm <- enzyme_sqm
rownames(mismatched_enzyme_sqm$orfs$tpm)[[2L]] <- "different"
try(
  build_enzyme_plot_table(mismatched_enzyme_sqm, "S_positive", "1.1.1.1"),
  silent = TRUE
)
missing_function_sqm <- enzyme_sqm
missing_function_sqm$orfs$table$KEGGFUN <- NULL
try(
  build_enzyme_plot_table(missing_function_sqm, "S_positive", "1.1.1.1"),
  silent = TRUE
)
try(
  build_enzyme_plot_table(enzyme_sqm, "missing", "1.1.1.1"),
  silent = TRUE
)

audit <- build_ko_expansion_audit(pathway_sqm$orfs$table)
try(build_ko_expansion_audit(data.frame(value = 1)), silent = TRUE)
normalize_manifest_ko_audit()
normalize_manifest_ko_audit(audit)
try(normalize_manifest_ko_audit(list()), silent = TRUE)
try(normalize_manifest_ko_audit(audit["input_orf_count"]), silent = TRUE)

orf_result <- build_orf_long_result(pathway_sqm, c("S_positive", "S_zero"))
orf_long <- orf_result$data

mismatched_sqm <- pathway_sqm
rownames(mismatched_sqm$orfs$tpm)[[1L]] <- "different_orf"
try(build_orf_long_result(mismatched_sqm, "S_positive"), silent = TRUE)
try(build_orf_long_result(pathway_sqm, "missing_sample"), silent = TRUE)

extract_ko_ec_lookup(orf_long)
extract_ko_ec_lookup(tibble::tibble(ko_id = c("K1", "K1")))
extract_ko_ec_lookup(tibble::tibble(
  ko_id = c("K1", "K1", "K2"),
  ec_codes = c("2.2.2.2; 1.1.1.-", "", NA_character_)
))
try(extract_ko_ec_lookup(tibble::tibble(value = 1)), silent = TRUE)

funz_orfs <- tibble::tibble(
  orf_id = c("orf_a", "orf_b", "orf_c"),
  sample = rep("S_positive", 3L),
  tpm = c(50, 30, 20),
  ko_id = c("K00001", "K00002", "K00003"),
  kegg_function = c("Function one", NA_character_, ""),
  ec_codes = c("1.1.1.1", NA_character_, "3.3.3.-"),
  phylum = c("Alpha", "Unclassified", "Beta")
)
suppressWarnings(build_ko_plot_table(
  orf_long = funz_orfs,
  selected_samples = c("S_positive", "S_zero"),
  top_n_ko = 2L,
  ko_lookup = c(K00002 = "Lookup two"),
  pathway_name = "Synthetic pathway"
))
suppressWarnings(build_ko_plot_table(
  orf_long = funz_orfs[0, ],
  selected_samples = c("S_zero_a", "S_zero_b"),
  top_n_ko = 2L,
  ko_lookup = character(),
  pathway_name = NA_character_
))

build_pie_chart_table(
  orf_long = funz_orfs,
  sample_name = "S_positive",
  ko_id_filter = "K00001",
  rank_name = "phylum",
  top_n_taxa = 1L
)
build_pie_chart_table(
  orf_long = funz_orfs,
  sample_name = "missing_sample",
  ko_id_filter = "K00001",
  rank_name = "phylum",
  top_n_taxa = 1L
)
