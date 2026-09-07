source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

script_env <- t1_source_sqm_plots_without_main()
selected_samples <- c("S1", "S2")
orf_ids <- c("orf_in_pathway", "orf_without_keggpath")

full_sqm <- list(
  orfs = list(
    table = data.frame(
      "KEGG ID" = rep("K01563", 2L),
      KEGGFUN = rep("haloalkane dehalogenase [EC:3.8.1.5]", 2L),
      KEGGPATH = c("Synthetic degradation pathway", ""),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tax = data.frame(
      superkingdom = rep("Bacteria", 2L),
      phylum = c("Alpha", "Beta"),
      class = rep("Synthetic class", 2L),
      order = rep("Synthetic order", 2L),
      family = rep("Synthetic family", 2L),
      genus = c("Genus A", "Genus B"),
      species = c("Species A", "Species B"),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tpm = data.frame(
      S1 = c(10, 20),
      S2 = c(5, 7),
      row.names = orf_ids,
      check.names = FALSE
    )
  ),
  functions = list(
    KEGG = list(
      tpm = data.frame(
        S1 = 30,
        S2 = 12,
        row.names = "K01563",
        check.names = FALSE
      )
    )
  ),
  misc = list(
    KEGG_names = c(K01563 = "haloalkane dehalogenase [EC:3.8.1.5]")
  )
)

fake_subset_fun <- function(
    SQM,
    fun,
    columns,
    ignore_case,
    fixed,
    allow_empty) {
  keep <- grepl(
    fun,
    as.character(SQM$orfs$table[[columns]]),
    ignore.case = ignore_case,
    fixed = fixed
  )
  subset_sqm <- SQM
  subset_sqm$orfs$table <- SQM$orfs$table[keep, , drop = FALSE]
  subset_sqm$orfs$tax <- SQM$orfs$tax[keep, , drop = FALSE]
  subset_sqm$orfs$tpm <- SQM$orfs$tpm[keep, , drop = FALSE]
  subset_sqm
}

pathway_sqm <- script_env$subset_pathway(
  full_sqm,
  "Synthetic degradation pathway",
  subset_fun = fake_subset_fun
)
t1_expect_identical(
  rownames(pathway_sqm$orfs$table),
  "orf_in_pathway",
  "Synthetic KEGGPATH subset"
)

pathway_analysis <- script_env$build_pathway_analysis(
  list(
    pathway_name = "Synthetic degradation pathway",
    pathway_id = "99999",
    pathway_selection = "defined",
    pathway_sqm = pathway_sqm,
    context_sqm = full_sqm,
    pathway_ko_ids = "K01563"
  ),
  selected_samples
)
flow_input <- pathway_analysis$flow_orf_long_result
flow_table <- script_env$build_flow_table_for_rank(
  orf_long = flow_input$data,
  rank = "phylum",
  selected_samples = selected_samples,
  top_n_taxa = 2L,
  top_n_ko = 1L,
  ko_lookup = script_env$get_ko_name_lookup(full_sqm)
)
flow_values <- t1_flow_ko_totals(flow_table, "K01563", selected_samples)

node_fixture <- data.frame(
  node_id = "synthetic_k01563",
  kegg_names = "K01563",
  label = "K01563",
  type = "ortholog",
  stringsAsFactors = FALSE
)
oracle_values <- unlist(
  t1_map_pathview_tpm(full_sqm, node_fixture, selected_samples)[1L, ],
  use.names = FALSE
)
comparison <- t1_flow_pathview_comparison(
  context = "synthetic_KEGGPATH_gap",
  selected_samples = selected_samples,
  oracle_values = oracle_values,
  flow_values = flow_values
)

# The textual pathway subset still contains only one ORF, but FLOW must use
# pathview KO membership and the complete SQM KEGG functional margin.
t1_assert_flow_pathview_parity(comparison)
