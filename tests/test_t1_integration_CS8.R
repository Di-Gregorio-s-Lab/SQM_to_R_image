source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

parse_choice <- function(args, option_name, choices) {
  prefix <- paste0("--", option_name, "=")
  matches <- args[startsWith(args, prefix)]
  if (length(matches) != 1L) {
    stop("Exactly one ", prefix, " option is required.", call. = FALSE)
  }
  value <- substring(matches, nchar(prefix) + 1L)
  if (!value %in% choices) {
    stop(
      option_name, " must be one of: ", paste(choices, collapse = ", "),
      call. = FALSE
    )
  }
  value
}

args <- commandArgs(trailingOnly = TRUE)
known_args <- grepl("^--(?:oracle|check)=", args)
if (length(args) != 2L || any(!known_args)) {
  stop(
    "Usage: Rscript tests/test_t1_integration_CS8.R ",
    "--oracle=fixture|live --check=oracle|parity",
    call. = FALSE
  )
}
oracle_mode <- parse_choice(args, "oracle", c("fixture", "live"))
check_mode <- parse_choice(args, "check", c("oracle", "parity"))

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
script_env$validate_samples(selected_samples, colnames(sqm$functions$KEGG$tpm))

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

map_and_assert_k01563 <- function(nodes, pathway_id, source_label) {
  mapped <- t1_map_pathview_tpm(sqm, nodes, selected_samples)
  for (node_id in rownames(mapped)) {
    t1_expect_equal(
      mapped[node_id, selected_samples],
      official_k01563,
      paste0(
        "CS8 K01563 pathview node ", pathway_id, "/", node_id,
        " from ", source_label
      ),
      tolerance = 1e-8
    )
  }
  invisible(mapped)
}

if (oracle_mode == "fixture") {
  for (pathway_id in pathway_ids) {
    nodes <- fixture[fixture$pathway_id == pathway_id, , drop = FALSE]
    t1_expect_true(nrow(nodes) > 0L, paste0("No fixture nodes for ", pathway_id))
    map_and_assert_k01563(nodes, pathway_id, "fixture")
  }
} else {
  download_dir <- tempfile("t1_kegg_live_")
  dir.create(download_dir, recursive = TRUE)
  on.exit(unlink(download_dir, recursive = TRUE, force = TRUE), add = TRUE)

  for (pathway_id in pathway_ids) {
    status <- pathview::download.kegg(
      pathway.id = pathway_id,
      species = "ko",
      kegg.dir = download_dir
    )
    if (!identical(unname(status), "succeed")) {
      stop("KEGG download failed for pathway ", pathway_id, call. = FALSE)
    }

    node_data <- pathview::node.info(
      file.path(download_dir, paste0("ko", pathway_id, ".xml"))
    )
    contains_k01563 <- vapply(
      node_data$kegg.names,
      function(values) any(values == "K01563"),
      logical(1)
    )
    contains_k01563 <- contains_k01563 & node_data$type == "ortholog"
    live_ids <- names(node_data$kegg.names)[contains_k01563]
    expected_ids <- fixture$node_id[fixture$pathway_id == pathway_id]

    if (!identical(sort(live_ids), sort(expected_ids))) {
      stop(
        "KEGG K01563 node drift for pathway ", pathway_id,
        "; expected=", paste(sort(expected_ids), collapse = ","),
        "; observed=", paste(sort(live_ids), collapse = ","),
        call. = FALSE
      )
    }

    live_nodes <- data.frame(
      node_id = live_ids,
      kegg_names = vapply(
        node_data$kegg.names[live_ids],
        paste,
        character(1),
        collapse = ";"
      ),
      label = unname(node_data$labels[live_ids]),
      type = unname(node_data$type[live_ids]),
      stringsAsFactors = FALSE
    )
    map_and_assert_k01563(live_nodes, pathway_id, "live KEGG")
  }
}

message(
  "PASS: SQMtools=", as.character(utils::packageVersion("SQMtools")),
  " CS8 K01563 oracle=", oracle_mode
)

if (check_mode == "parity") {
  pathway_names <- c(
    "00361" = "Chlorocyclohexane and chlorobenzene degradation",
    "00625" = "Chloroalkane and chloroalkene degradation"
  )
  all_k01563 <- grepl(
    "K01563",
    as.character(sqm$orfs$table[["KEGG ID"]]),
    fixed = TRUE
  )
  t1_expect_identical(
    sum(all_k01563),
    4L,
    "CS8 complete K01563 ORF count"
  )

  comparisons <- vector("list", length(pathway_ids))
  for (pathway_index in seq_along(pathway_ids)) {
    pathway_id <- pathway_ids[[pathway_index]]
    pathway_name <- unname(pathway_names[[pathway_id]])
    pathway_sqm <- script_env$subset_pathway(sqm, pathway_name)
    subset_k01563 <- grepl(
      "K01563",
      as.character(pathway_sqm$orfs$table[["KEGG ID"]]),
      fixed = TRUE
    )
    t1_expect_identical(
      sum(subset_k01563),
      1L,
      paste0("CS8 ", pathway_id, " KEGGPATH-selected K01563 ORF count")
    )

    orf_long <- script_env$build_orf_long_table(pathway_sqm, selected_samples)
    flow_table <- script_env$build_flow_table_for_rank(
      orf_long = orf_long,
      rank = "phylum",
      selected_samples = selected_samples,
      top_n_taxa = length(unique(orf_long$phylum)),
      top_n_ko = length(unique(orf_long$ko_id)),
      ko_lookup = script_env$get_ko_name_lookup(pathway_sqm)
    )
    flow_values <- t1_flow_ko_totals(flow_table, "K01563", selected_samples)
    comparisons[[pathway_index]] <- t1_flow_pathview_comparison(
      context = pathway_id,
      selected_samples = selected_samples,
      oracle_values = official_k01563,
      flow_values = flow_values
    )
  }

  # Intentional RED until tranche 2.
  t1_assert_flow_pathview_parity(do.call(rbind, comparisons))
}
