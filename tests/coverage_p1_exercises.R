# This file is sourced by covr::file_coverage() after sqm_plots.R has been
# instrumented in the same environment. It intentionally does not source the
# production script itself.

valid_hierarchies <- c(
  "Metabolism; Carbohydrate metabolism; Alpha pathway",
  "Environmental Information Processing; Signal transduction; Beta pathway",
  "Brite Hierarchies; Protein families; Transporters",
  "Metabolism; Incomplete hierarchy"
)
parse_kegg_pathway_entries(valid_hierarchies)
parse_kegg_pathway_entries(NA_character_)
parse_kegg_pathway_entries(character())
split_kegg_pathway_field(valid_hierarchies)

selection_sqm <- list(
  orfs = list(
    table = data.frame(
      KEGGPATH = c(valid_hierarchies[[1]], valid_hierarchies[[2]]),
      row.names = c("orf_alpha", "orf_beta"),
      check.names = FALSE
    ),
    tpm = data.frame(
      S0 = c(60, 5),
      S1 = c(40, 5),
      row.names = c("orf_alpha", "orf_beta"),
      check.names = FALSE
    )
  )
)
selected_top <- select_top_pathways(selection_sqm, c("S0", "S1"), 2L)

ambiguous_sqm <- selection_sqm
ambiguous_sqm$orfs$table$KEGGPATH[[1]] <- paste(
  "Metabolism; Category A; Shared pathway",
  "Human Diseases; Category B; Shared pathway",
  sep = " | "
)
try(select_top_pathways(ambiguous_sqm, c("S0", "S1"), 2L), silent = TRUE)

empty_sqm <- selection_sqm
empty_sqm$orfs$table$KEGGPATH <- rep(
  "Brite Hierarchies; Protein families; Transporters",
  2L
)
suppressWarnings(select_top_pathways(empty_sqm, c("S0", "S1"), 2L))

try(
  select_top_pathways(
    list(orfs = list(
      table = data.frame(value = 1, row.names = "orf"),
      tpm = data.frame(S0 = 1, row.names = "orf")
    )),
    "S0",
    1L
  ),
  silent = TRUE
)
mismatched_sqm <- selection_sqm
rownames(mismatched_sqm$orfs$tpm)[[1]] <- "different_orf"
try(select_top_pathways(mismatched_sqm, "S0", 1L), silent = TRUE)

resolved_defined <- list(list(
  input_value = "00361",
  pathway_id = "00361",
  canonical_pathway_name = "Chlorocyclohexane and chlorobenzene degradation"
))
select_context_pathway_groups(
  context_sqm = selection_sqm,
  resolved_defined_pathways = resolved_defined,
  pathway_selection_modes = c("defined", "top20"),
  selected_samples = c("S0", "S1"),
  pathway_top_n = 1L
)
suppressWarnings(select_context_pathway_groups(
  context_sqm = empty_sqm,
  resolved_defined_pathways = list(),
  pathway_selection_modes = "top20",
  selected_samples = c("S0", "S1"),
  pathway_top_n = 1L
))

flow_orfs <- tibble::tibble(
  orf_id = c("orf_b", "orf_a", "orf_c", "orf_fallback", "orf_id"),
  sample = c("S0", "S0", "S0", "S1", "S1"),
  tpm = c(10, 5, 15, 12, 8),
  ko_id = c("K00001", "K00001", "K00001", "K00002", "K00003"),
  kegg_function = c(
    "Description B", "Description A", "Description A", " ", NA_character_
  ),
  ec_codes = c("2.2.2.2", "1.1.1.1", "2.2.2.2", NA_character_, NA_character_),
  phylum = c("Alpha", "Alpha", "Beta", "Alpha", "Beta")
)
ko_lookup <- c(K00002 = "Lookup description")
flow_meta <- build_flow_ko_metadata(flow_orfs, ko_lookup)
try(build_flow_ko_metadata(flow_orfs["ko_id"], ko_lookup), silent = TRUE)

summary_tbl <- tibble::tibble(
  sample = c("S0", "S0"),
  taxon = c("Alpha", "Beta"),
  KO = c("K00001", "K00001"),
  TPM = c(15, 15)
)
join_flow_ko_metadata(summary_tbl, flow_meta)
try(join_flow_ko_metadata(summary_tbl[1:3], flow_meta), silent = TRUE)
try(join_flow_ko_metadata(summary_tbl, flow_meta["ko_id"]), silent = TRUE)
duplicate_meta <- dplyr::bind_rows(flow_meta, flow_meta[1, , drop = FALSE])
try(join_flow_ko_metadata(summary_tbl, duplicate_meta), silent = TRUE)
duplicate_summary <- dplyr::bind_rows(summary_tbl, summary_tbl[1, , drop = FALSE])
try(join_flow_ko_metadata(duplicate_summary, flow_meta), silent = TRUE)

flow_rank <- build_flow_table_for_rank(
  orf_long = flow_orfs,
  rank = "phylum",
  selected_samples = c("S0", "S1"),
  top_n_taxa = 1L,
  top_n_ko = 2L,
  ko_lookup = ko_lookup
)
flow_sample <- build_flow_table_for_sample(
  flow_rank,
  "Synthetic pathway",
  "phylum",
  "S0"
)
build_flow_legend_spec(flow_sample)
make_flow_plot(flow_sample, "Synthetic pathway", "phylum", "S0")
make_flow_sankey(flow_sample, "Synthetic pathway", "phylum", "S0")
