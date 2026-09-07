source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

script_env <- t1_source_sqm_plots_without_main()

node_data <- list(
  kegg.names = list(
    `1` = c("K00001", "K00002"),
    `2` = "K00002",
    `3` = "C00001"
  ),
  type = c(`1` = "ortholog", `2` = "ortholog", `3` = "compound")
)

t1_expect_identical(
  script_env$extract_pathway_ko_ids(node_data),
  c("K00001", "K00002"),
  "FLOW pathway KO membership must follow pathview ortholog nodes"
)

orf_ids <- c("orf_multi", "orf_single", "orf_other")
sqm <- list(
  orfs = list(
    table = data.frame(
      "KEGG ID" = c("K00001;K00002", "K00001", "K99999"),
      KEGGFUN = c("multi", "single", "other"),
      KEGGPATH = c("", "Synthetic pathway", "Synthetic pathway"),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tax = data.frame(
      superkingdom = rep("Bacteria", 3L),
      phylum = c("Alpha", "Beta", "Gamma"),
      class = rep("Synthetic class", 3L),
      order = rep("Synthetic order", 3L),
      family = rep("Synthetic family", 3L),
      genus = c("Genus A", "Genus B", "Genus C"),
      species = c("Species A", "Species B", "Species C"),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tpm = data.frame(
      S1 = c(30, 10, 100),
      row.names = orf_ids,
      check.names = FALSE
    )
  ),
  functions = list(
    KEGG = list(
      tpm = data.frame(
        S1 = c(25, 15, 100),
        row.names = c("K00001", "K00002", "K99999"),
        check.names = FALSE
      )
    )
  ),
  misc = list(
    KEGG_names = c(K00001 = "KO one", K00002 = "KO two", K99999 = "Other KO")
  )
)

prepared <- script_env$prepare_context_pathway_subsets(
  context_sqm = sqm,
  pathway_entries = list(list(
    pathway_name = "Synthetic pathway",
    pathway_id = "99999",
    pathway_selection = "defined"
  )),
  context_label = "synthetic",
  subset_fun = function(SQM, ...) SQM,
  include_kegg_oracle = TRUE,
  pathway_ko_resolver = function(pathway_id) c("K00001", "K00002")
)
pathway_analysis <- script_env$build_pathway_analysis(
  prepared$pathway_sqms[[1L]],
  "S1"
)
flow_input <- pathway_analysis$orf_long_result

flow_rank <- script_env$build_flow_table_for_rank(
  orf_long = flow_input$data,
  rank = "phylum",
  selected_samples = "S1",
  top_n_taxa = 10L,
  top_n_ko = 10L,
  ko_lookup = script_env$get_ko_name_lookup(sqm)
)

observed <- flow_rank |>
  dplyr::arrange(.data$KO, .data$taxon) |>
  dplyr::select(all_of(c("sample", "taxon", "KO", "TPM")))

expected <- tibble::tibble(
  sample = rep("S1", 3L),
  taxon = c("Alpha", "Beta", "Alpha"),
  KO = c("K00001", "K00001", "K00002"),
  TPM = c(15, 10, 15)
)

t1_expect_identical(
  observed[c("sample", "taxon", "KO")],
  expected[c("sample", "taxon", "KO")],
  "FLOW SQMtools multi-KO taxonomic allocation keys"
)
t1_expect_equal(
  observed$TPM,
  expected$TPM,
  "FLOW SQMtools multi-KO taxonomic allocation values",
  tolerance = 1e-12
)

empty_flow <- script_env$build_pathway_ko_result(
  context_sqm = sqm,
  selected_samples = "S1",
  pathway_ko_ids = character()
)
t1_expect_identical(
  nrow(empty_flow$data),
  0L,
  "FLOW pathway without ortholog nodes must be empty"
)

message("PASS: FLOW uses pathview KO membership and SQMtools multi-KO TPM semantics")
