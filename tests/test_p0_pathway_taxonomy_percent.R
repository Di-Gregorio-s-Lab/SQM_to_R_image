source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

script_env <- t1_source_sqm_plots_without_main()
taxonomy_abund <- rbind(
  Alpha = c(S_positive = 10, S_zero = 0),
  Beta = c(S_positive = 5, S_zero = 0)
)
# These percentages retain the denominator of the complete sample library.
# The pathway contributes 1.5%, so its own bars must not be rescaled to 100%.
taxonomy_percent <- rbind(
  Alpha = c(S_positive = 1.0, S_zero = 0),
  Beta = c(S_positive = 0.5, S_zero = 0)
)
pathway_sqm <- structure(
  list(
    misc = list(samples = c("S_positive", "S_zero")),
    taxa = list(
      phylum = list(
        abund = taxonomy_abund,
        percent = taxonomy_percent
      )
    )
  ),
  class = "SQM"
)

oracle <- SQMtools::plotTaxonomy(
  SQM = pathway_sqm,
  rank = "phylum",
  count = "percent",
  N = 2L,
  others = TRUE,
  samples = c("S_positive", "S_zero"),
  ignore_unmapped = FALSE,
  ignore_unclassified = FALSE,
  no_partial_classifications = FALSE,
  rescale = FALSE
)
observed <- script_env$make_taxonomy_plot(
  sqm_object = pathway_sqm,
  rank = "phylum",
  count = "percent",
  selected_samples = c("S_positive", "S_zero"),
  top_n_taxa = 2L,
  ignore_unmapped = FALSE,
  ignore_unclassified = FALSE
)
t1_assert_plot_data(
  observed$data,
  oracle$data,
  "Pathway plotTaxonomy percent oracle",
  tolerance = 1e-12
)

positive_sum <- sum(observed$data$abun[observed$data$sample == "S_positive"])
t1_expect_equal(
  positive_sum,
  1.5,
  "Pathway percent must retain the whole-sample denominator",
  tolerance = 1e-12
)
t1_expect_true(
  abs(positive_sum - 100) > 1e-8,
  "Pathway percent was incorrectly forced to 100"
)

tsv_data <- script_env$extract_taxonomy_plot_data(observed, "percent") |>
  dplyr::select(all_of(c("sample", "taxon", "value", "count")))
t1_expect_identical(
  colnames(tsv_data),
  c("sample", "taxon", "value", "count"),
  "Pathway taxonomy TSV schema"
)

message("PASS: pathway taxonomy percent delegates to non-rescaled plotTaxonomy")
