# This file is sourced by covr::file_coverage() after sqm_plots.R has been
# instrumented in the same environment. It intentionally does not source the
# production script itself.

mapping_sqm <- list(
  misc = list(
    KEGG_paths = c(
      "Metabolism; Energy metabolism; Carbon fixation in photosynthetic organisms",
      "Metabolism; Xenobiotics biodegradation and metabolism; Nitrotoluene degradation",
      "Metabolism; Energy metabolism; Nitrogen metabolism"
    )
  )
)
resolve_pathways(mapping_sqm, c("00710", "00633", "00910"))
resolve_pathways(mapping_sqm, "Nitrogen metabolism")
try(resolve_pathways(mapping_sqm, "00643"), silent = TRUE)
try(resolve_pathways(mapping_sqm, "p00910"), silent = TRUE)
try(resolve_pathways(mapping_sqm, "Nitrogen"), silent = TRUE)
try(resolve_pathways(mapping_sqm, "No matching pathway"), silent = TRUE)
resolve_pathways(
  list(misc = list(KEGG_paths = "Metabolism; Synthetic category; Uncurated pathway")),
  "Uncurated pathway"
)
pathview_is_exportable("defined", "00910")
pathview_is_exportable("defined", NA_character_)
pathview_is_exportable("top20", NA_character_)

orf_ids <- c("orf_a", "orf_b", "orf_c", "orf_d", "orf_e")
filter_sqm <- list(
  orfs = list(
    table = data.frame(value = seq_along(orf_ids), row.names = orf_ids),
    tax = data.frame(
      superkingdom = rep("Bacteria", length(orf_ids)),
      phylum = c("Bacillota", "Other phylum", "Bacillota", "Third", "Fourth"),
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
resolved_taxon <- resolve_taxa_filters(filter_sqm, "Bacillota")[[1]]
try(resolve_taxa_filters(filter_sqm, "Missing taxon"), silent = TRUE)

fake_subset_orfs <- function(
    SQM,
    orfs,
    tax_source = "orfs",
    trusted_functions_only = FALSE,
    ignore_unclassified_functions = FALSE,
    rescale_tpm = FALSE,
    rescale_copy_number = FALSE,
    recalculate_bin_stats = TRUE,
    contigs_override = NULL,
    allow_empty = FALSE) {
  result <- SQM
  result$orfs$table <- SQM$orfs$table[orfs, , drop = FALSE]
  result$orfs$tax <- SQM$orfs$tax[orfs, , drop = FALSE]
  result$orfs$tpm <- SQM$orfs$tpm[orfs, , drop = FALSE]
  result
}
subset_sqm_by_taxon(filter_sqm, resolved_taxon$orf_ids, fake_subset_orfs)
try(
  subset_sqm_by_taxon(
    filter_sqm,
    resolved_taxon$orf_ids,
    function(...) fake_subset_orfs(filter_sqm, resolved_taxon$orf_ids[[1]])
  ),
  silent = TRUE
)

taxonomy_sqm <- filter_sqm
taxonomy_sqm$orfs$tax$phylum <- c(
  "Alpha", "Beta", "Gamma", "Delta", NA_character_
)
taxonomy_table <- suppressWarnings(build_pathway_taxonomy_percent_table(
  sqm_object = taxonomy_sqm,
  rank = "phylum",
  selected_samples = c("S_positive", "S_zero"),
  top_n_taxa = 2L,
  pathway_name = "Synthetic pathway"
))
make_pathway_taxonomy_percent_plot(
  plot_tbl = taxonomy_table,
  pathway_name = "Synthetic pathway",
  rank = "phylum",
  selected_samples = c("S_positive", "S_zero")
)
try(
  build_pathway_taxonomy_percent_table(
    taxonomy_sqm,
    "missing_rank",
    "S_positive",
    2L,
    "Synthetic pathway"
  ),
  silent = TRUE
)

global_taxonomy_sqm <- list(
  total_reads = c(S_positive = 100),
  taxa = list(phylum = list(percent = matrix(
    c(70, 10, 20),
    ncol = 1L,
    dimnames = list(c("Unmapped", "Unclassified", "Alpha"), "S_positive")
  )))
)
global_plot_data <- tibble::tibble(
  sample = "S_positive",
  taxon = "Alpha",
  value = 20,
  count = "percent"
)
add_global_taxonomy_percent_metadata(
  global_plot_data,
  global_taxonomy_sqm,
  "phylum",
  "S_positive"
)
no_exclusion_sqm <- global_taxonomy_sqm
no_exclusion_sqm$taxa$phylum$percent <- matrix(
  100,
  ncol = 1L,
  dimnames = list("Alpha", "S_positive")
)
add_global_taxonomy_percent_metadata(
  tibble::tibble(
    sample = "S_positive",
    taxon = "Alpha",
    value = 100,
    count = "percent"
  ),
  no_exclusion_sqm,
  "phylum",
  "S_positive"
)
try(
  add_global_taxonomy_percent_metadata(
    global_plot_data |> dplyr::select(-.data$count),
    global_taxonomy_sqm,
    "phylum",
    "S_positive"
  ),
  silent = TRUE
)
