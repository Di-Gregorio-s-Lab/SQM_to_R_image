source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

parse_oracle_mode <- function(args) {
  matches <- args[startsWith(args, "--oracle=")]
  if (length(args) != 1L || length(matches) != 1L) {
    stop(
      "Usage: Rscript tests/test_t3_integration_CS8.R --oracle=fixture|live",
      call. = FALSE
    )
  }
  mode <- substring(matches, nchar("--oracle=") + 1L)
  if (!mode %in% c("fixture", "live")) {
    stop("oracle must be fixture or live.", call. = FALSE)
  }
  mode
}

oracle_mode <- parse_oracle_mode(commandArgs(trailingOnly = TRUE))
project_dir <- Sys.getenv("SQM_CS8_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) {
  stop("SQM_CS8_PROJECT_DIR is required.", call. = FALSE)
}
if (!dir.exists(project_dir)) {
  stop("SQM_CS8_PROJECT_DIR does not exist: ", project_dir, call. = FALSE)
}

script_env <- t1_source_sqm_plots_without_main()
sqm <- script_env$load_sqm_project(project_dir, "prokfilter")
selected_samples <- c("CS8T0", "CS8T2", "CS8T3", "CS8T4", "CS8T6")
expected_k01563 <- c(
  20.461747379587,
  64.241150618838,
  48.990246839864,
  111.090699101642,
  38.836884026958
)
official_k01563 <- unname(sqm$functions$KEGG$tpm["K01563", selected_samples])
t1_expect_equal(
  official_k01563,
  expected_k01563,
  "CS8 K01563 official SQM TPM drifted",
  tolerance = 1e-8
)

fixture <- t1_read_k01563_fixture()
pathway_ids <- c("00361", "00625")
pathway_names <- c(
  `00361` = "Chlorocyclohexane and chlorobenzene degradation",
  `00625` = "Chloroalkane and chloroalkene degradation"
)

pathway_ko_ids <- stats::setNames(vector("list", length(pathway_ids)), pathway_ids)
if (oracle_mode == "fixture") {
  for (pathway_id in pathway_ids) {
    nodes <- fixture[fixture$pathway_id == pathway_id, , drop = FALSE]
    t1_expect_true(nrow(nodes) > 0L, paste0("No K01563 fixture nodes for ", pathway_id))
    mapped <- t1_map_pathview_tpm(sqm, nodes, selected_samples)
    for (node_id in rownames(mapped)) {
      t1_expect_equal(
        mapped[node_id, selected_samples],
        official_k01563,
        paste0("Pathview fixture node ", pathway_id, "/", node_id),
        tolerance = 1e-8
      )
    }
    pathway_ko_ids[[pathway_id]] <- "K01563"
  }
} else {
  for (pathway_id in pathway_ids) {
    node_data <- script_env$download_pathway_node_data(pathway_id)
    pathway_ko_ids[[pathway_id]] <- script_env$extract_pathway_ko_ids(node_data)
    k01563_nodes <- names(node_data$type)[
      node_data$type == "ortholog" &
        vapply(
          node_data$kegg.names,
          function(values) "K01563" %in% as.character(values),
          logical(1)
        )
    ]
    expected_nodes <- fixture$node_id[fixture$pathway_id == pathway_id]
    t1_expect_identical(
      sort(k01563_nodes),
      sort(expected_nodes),
      paste0("Live KEGG K01563 node drift for ", pathway_id)
    )
  }
}

all_k01563_orfs <- rownames(sqm$orfs$table)[vapply(
  as.character(sqm$orfs$table[["KEGG ID"]]),
  function(value) "K01563" %in% script_env$extract_ko_ids(value),
  logical(1)
)]
t1_expect_identical(length(all_k01563_orfs), 4L, "CS8 K01563 ORF count")

enzyme <- script_env$build_enzyme_plot_table(
  sqm,
  selected_samples,
  "3.8.1.5"
)
enzyme_k01563 <- enzyme$tpm[
  as.character(enzyme$ec_code) == "3.8.1.5" &
    as.character(enzyme$sample) %in% selected_samples
]
t1_expect_equal(
  enzyme_k01563,
  official_k01563,
  "CS8 ENZIMI 3.8.1.5 parity",
  tolerance = 1e-8
)

for (pathway_id in pathway_ids) {
  prepared <- script_env$prepare_context_pathway_subsets(
    context_sqm = sqm,
    pathway_entries = list(list(
      pathway_name = unname(pathway_names[[pathway_id]]),
      pathway_id = pathway_id,
      pathway_selection = "defined"
    )),
    context_label = "CS8",
    include_kegg_oracle = TRUE,
    pathway_ko_resolver = function(resolved_id) pathway_ko_ids[[resolved_id]],
    selected_samples = selected_samples
  )
  t1_expect_identical(
    length(prepared$pathway_sqms),
    1L,
    paste0("Productive pathway preparation ", pathway_id)
  )
  analysis <- script_env$build_pathway_analysis(
    prepared$pathway_sqms[[1L]],
    selected_samples
  )

  t1_expect_true(
    all(all_k01563_orfs %in% analysis$orf_long_result$orf_ids),
    paste0("Canonical membership omitted K01563 ORFs for ", pathway_id)
  )
  t1_expect_identical(
    anyDuplicated(analysis$orf_long_result$orf_ids),
    0L,
    paste0("Canonical membership duplicated ORFs for ", pathway_id)
  )

  funz <- script_env$build_ko_plot_table(
    orf_long = analysis$orf_long_result$totals,
    selected_samples = selected_samples,
    top_n_ko = max(1L, length(pathway_ko_ids[[pathway_id]])),
    ko_lookup = analysis$ko_lookup,
    pathway_name = analysis$pathway_name
  )
  funz_values <- vapply(selected_samples, function(sample_name) {
    sum(funz$tpm[
      as.character(funz$sample) == sample_name &
        as.character(funz$ko_id) == "K01563"
    ])
  }, numeric(1))

  flow <- script_env$build_flow_table_for_rank(
    orf_long = analysis$orf_long_result$data,
    rank = "phylum",
    selected_samples = selected_samples,
    top_n_taxa = max(1L, length(unique(analysis$orf_long_result$data$phylum))),
    top_n_ko = max(1L, length(pathway_ko_ids[[pathway_id]])),
    ko_lookup = analysis$ko_lookup
  )
  flow_values <- t1_flow_ko_totals(flow, "K01563", selected_samples)

  pie_values <- vapply(selected_samples, function(sample_name) {
    pie <- script_env$build_pie_chart_table(
      orf_long = analysis$orf_long_result$data,
      sample_name = sample_name,
      ko_id_filter = "K01563",
      rank_name = "phylum",
      top_n_taxa = 100L,
      pathway_name = analysis$pathway_name,
      pathway_id = pathway_id,
      pathway_selection = "defined"
    )
    sum(pie$tpm)
  }, numeric(1))

  t1_expect_equal(funz_values, official_k01563, paste0("FUNZ parity ", pathway_id))
  t1_expect_equal(flow_values, official_k01563, paste0("FLOW parity ", pathway_id))
  t1_expect_equal(pie_values, official_k01563, paste0("PIE parity ", pathway_id))

  for (count in c("abund", "percent")) {
    oracle <- SQMtools::plotTaxonomy(
      SQM = analysis$pathway_sqm,
      rank = "phylum",
      count = count,
      N = 15L,
      others = TRUE,
      samples = selected_samples,
      ignore_unmapped = FALSE,
      ignore_unclassified = FALSE,
      no_partial_classifications = FALSE,
      rescale = FALSE
    )
    observed <- script_env$make_taxonomy_plot(
      sqm_object = analysis$pathway_sqm,
      rank = "phylum",
      count = count,
      selected_samples = selected_samples,
      top_n_taxa = 15L,
      ignore_unmapped = FALSE,
      ignore_unclassified = FALSE
    )
    t1_assert_plot_data(
      observed$data,
      oracle$data,
      paste0("plotTaxonomy ", count, " parity ", pathway_id),
      tolerance = 1e-8
    )
    exported_taxonomy <- script_env$extract_taxonomy_plot_data(observed, count) |>
      dplyr::select(all_of(c("sample", "taxon", "value", "count")))
    t1_expect_identical(
      colnames(exported_taxonomy),
      c("sample", "taxon", "value", "count"),
      paste0("Pathway taxonomy TSV schema ", count, " ", pathway_id)
    )
    if (count == "percent") {
      percent_sums <- stats::aggregate(
        observed$data$abun,
        by = list(sample = as.character(observed$data$sample)),
        FUN = sum
      )
      t1_expect_true(
        all(percent_sums$x < 100 - 1e-8),
        paste0("Pathway percentages were incorrectly rescaled to 100 for ", pathway_id)
      )
    }
  }
}

message(
  "PASS: SQMtools=", as.character(utils::packageVersion("SQMtools")),
  " CS8 cross-mode oracle=", oracle_mode
)
