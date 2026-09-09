source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

script_env <- t1_source_sqm_plots_without_main()
selected_samples <- c("S1", "S2")
orf_ids <- c("orf_explicit", "orf_without_fields", "orf_multi_ko")

full_sqm <- list(
  orfs = list(
    table = data.frame(
      "KEGG ID" = c("K01563", "K01563", "K01563;K11991"),
      KEGGFUN = c(
        "haloalkane dehalogenase [EC:3.8.1.5]",
        "",
        ""
      ),
      KEGGPATH = c("Synthetic degradation pathway", "", ""),
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
      S1 = c(10, 20, 12),
      S2 = c(5, 7, 8),
      row.names = orf_ids,
      check.names = FALSE
    )
  ),
  functions = list(
    KEGG = list(
      tpm = data.frame(
        S1 = c(100, 30),
        S2 = c(50, 10),
        row.names = c("K01563", "K11991"),
        check.names = FALSE
      )
    )
  ),
  misc = list(
    KEGG_names = c(
      K01563 = "haloalkane dehalogenase [EC:3.8.1.5]",
      K11991 = "synthetic companion [EC:9.9.9.9]"
    )
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

legacy_pathway_sqm <- script_env$subset_pathway(
  full_sqm,
  "Synthetic degradation pathway",
  subset_fun = fake_subset_fun
)
legacy_pathway_sqm$functions$KEGG$tpm <- data.frame(
  S1 = 10,
  S2 = 5,
  row.names = "K01563",
  check.names = FALSE
)
pathway_analysis <- script_env$build_pathway_analysis(
  list(
    pathway_name = "Synthetic degradation pathway",
    pathway_id = "99999",
    pathway_selection = "defined",
    pathway_sqm = legacy_pathway_sqm,
    context_sqm = full_sqm,
    pathway_ko_ids = "K01563"
  ),
  selected_samples
)

canonical_input <- pathway_analysis$orf_long_result$data
funz <- script_env$build_ko_plot_table(
  orf_long = canonical_input,
  selected_samples = selected_samples,
  top_n_ko = 10L,
  ko_lookup = script_env$get_ko_name_lookup(full_sqm),
  pathway_name = "Synthetic degradation pathway"
)
funz_values <- vapply(selected_samples, function(sample_name) {
  sum(funz$tpm[
    as.character(funz$sample) == sample_name &
      as.character(funz$ko_id) == "K01563"
  ])
}, numeric(1))

pie_values <- vapply(selected_samples, function(sample_name) {
  pie <- script_env$build_pie_chart_table(
    orf_long = canonical_input,
    sample_name = sample_name,
    ko_id_filter = "K01563",
    rank_name = "phylum",
    top_n_taxa = 10L,
    pathway_name = "Synthetic degradation pathway",
    pathway_id = "99999",
    pathway_selection = "defined"
  )
  sum(pie$tpm)
}, numeric(1))

enzyme <- script_env$build_enzyme_plot_table(
  full_sqm,
  selected_samples,
  "3.8.1.5"
)
enzyme_values <- vapply(selected_samples, function(sample_name) {
  sum(enzyme$tpm[
    as.character(enzyme$sample) == sample_name &
      as.character(enzyme$ec_code) == "3.8.1.5"
  ])
}, numeric(1))

oracle_values <- as.numeric(
  legacy_pathway_sqm$functions$KEGG$tpm["K01563", selected_samples]
)
comparisons <- rbind(
  data.frame(
    mode = rep("FUNZ", length(selected_samples)),
    sample = selected_samples,
    oracle = oracle_values,
    observed = funz_values
  ),
  data.frame(
    mode = rep("PIE", length(selected_samples)),
    sample = selected_samples,
    oracle = oracle_values,
    observed = pie_values
  ),
  data.frame(
    mode = "TAXONOMY_MEMBERSHIP",
    sample = "all",
    oracle = nrow(legacy_pathway_sqm$orfs$table),
    observed = nrow(pathway_analysis$pathway_sqm$orfs$table)
  )
)
comparisons$delta <- comparisons$observed - comparisons$oracle

mismatched <- abs(comparisons$delta) > 1e-8
if (any(mismatched)) {
  details <- apply(comparisons[mismatched, , drop = FALSE], 1L, function(row) {
    paste0(
      "mode=", row[["mode"]],
      " sample=", row[["sample"]],
      " oracle=", format(as.numeric(row[["oracle"]]), digits = 15),
      " observed=", format(as.numeric(row[["observed"]]), digits = 15),
      " delta=", format(as.numeric(row[["delta"]]), digits = 15)
    )
  })
  stop(
    "CROSS-MODE SQMTOOLS PARITY FAILURE\n",
    paste(details, collapse = "\n"),
    call. = FALSE
  )
}

stopifnot(isTRUE(all.equal(
  enzyme_values,
  as.numeric(full_sqm$functions$KEGG$tpm["K01563", selected_samples])
)))
message("PASS: FUNZ, PIE and taxonomy share the pathway subset; ENZIMI stays global")
