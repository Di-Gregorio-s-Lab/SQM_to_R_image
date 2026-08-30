source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")

if (!exists(
  "select_context_pathway_groups",
  envir = script_env,
  mode = "function",
  inherits = FALSE
)) {
  stop(
    paste0(
      "BUG-P1-02 RED: select_context_pathway_groups() is missing; ",
      "top20 cannot yet be resolved from each filtered SQM context."
    ),
    call. = FALSE
  )
}

make_context_sqm <- function(alpha_tpm, beta_tpm) {
  stopifnot(length(alpha_tpm) == 2L, length(beta_tpm) == 2L)

  list(
    orfs = list(
      table = data.frame(
        KEGGPATH = c(
          "Metabolism; Carbohydrate metabolism; Alpha pathway",
          "Environmental Information Processing; Signal transduction; Beta pathway"
        ),
        row.names = c("orf_alpha", "orf_beta"),
        check.names = FALSE
      ),
      tpm = data.frame(
        S0 = c(alpha_tpm[[1]], beta_tpm[[1]]),
        S1 = c(alpha_tpm[[2]], beta_tpm[[2]]),
        row.names = c("orf_alpha", "orf_beta"),
        check.names = FALSE
      )
    )
  )
}

alpha_context_sqm <- make_context_sqm(
  alpha_tpm = c(60, 40),
  beta_tpm = c(3, 2)
)
beta_context_sqm <- make_context_sqm(
  alpha_tpm = c(1, 2),
  beta_tpm = c(50, 40)
)

resolved_defined_pathways <- list(list(
  input_value = "00361",
  pathway_id = "00361",
  canonical_pathway_name = "Chlorocyclohexane and chlorobenzene degradation"
))

select_for_context <- function(context_sqm) {
  script_env$select_context_pathway_groups(
    context_sqm = context_sqm,
    resolved_defined_pathways = resolved_defined_pathways,
    pathway_selection_modes = c("defined", "top20"),
    selected_samples = c("S0", "S1"),
    pathway_top_n = 1L
  )
}

alpha_groups <- select_for_context(alpha_context_sqm)
beta_groups <- select_for_context(beta_context_sqm)

stopifnot(identical(names(alpha_groups), c("defined", "top20")))
stopifnot(identical(names(beta_groups), c("defined", "top20")))

# Defined pathways are resolved once on the complete SQM object and reused
# unchanged, independently of the filtered context used for top20 ranking.
stopifnot(identical(alpha_groups$defined, resolved_defined_pathways))
stopifnot(identical(beta_groups$defined, resolved_defined_pathways))

stopifnot(length(alpha_groups$top20) == 1L)
stopifnot(length(beta_groups$top20) == 1L)

alpha_top <- alpha_groups$top20[[1]]
beta_top <- beta_groups$top20[[1]]

stopifnot(identical(alpha_top$canonical_pathway_name, "Alpha pathway"))
stopifnot(identical(alpha_top$pathway_root, "Metabolism"))
stopifnot(identical(alpha_top$pathway_category, "Carbohydrate metabolism"))
stopifnot(isTRUE(all.equal(alpha_top$total_tpm, 100, tolerance = 1e-12)))

stopifnot(identical(beta_top$canonical_pathway_name, "Beta pathway"))
stopifnot(identical(
  beta_top$pathway_root,
  "Environmental Information Processing"
))
stopifnot(identical(beta_top$pathway_category, "Signal transduction"))
stopifnot(isTRUE(all.equal(beta_top$total_tpm, 90, tolerance = 1e-12)))

stopifnot(!identical(
  alpha_top$canonical_pathway_name,
  beta_top$canonical_pathway_name
))

# A top20 artifact emitted under a taxon-filtered context must retain both the
# selected pathway and the context that determined its ranking.
manifest_row <- script_env$new_manifest_row(
  script_name = "sqm_plots.R",
  project_dir = "project",
  tax_mode = "prokfilter",
  pathway = alpha_top$canonical_pathway_name,
  pathway_id = alpha_top$pathway_id,
  samples = c("S0", "S1"),
  metric = "tpm",
  top_n_taxa = 5L,
  top_n_ko = 5L,
  output_type = "data_tsv",
  output_file = file.path(
    "taxon_filter",
    "phylum",
    "Alpha_taxon",
    "flow",
    "top20",
    "Alpha_pathway",
    "flow_data.tsv"
  ),
  mode = "flow",
  rank = "phylum",
  format = "tsv",
  output_scope = "pathway_top20",
  filtered_taxon = "Alpha taxon",
  filtered_taxon_rank = "phylum"
)

stopifnot(nrow(manifest_row) == 1L)
stopifnot(identical(manifest_row$pathway, "Alpha pathway"))
stopifnot(identical(manifest_row$output_scope, "pathway_top20"))
stopifnot(identical(manifest_row$filtered_taxon, "Alpha taxon"))
stopifnot(identical(manifest_row$filtered_taxon_rank, "phylum"))
stopifnot(grepl(
  "taxon_filter/phylum/Alpha_taxon/flow/top20/Alpha_pathway",
  gsub("\\\\", "/", manifest_row$output_file),
  fixed = TRUE
))

message(
  "PASS: Defined pathways stay global while top20 is ranked within each SQM context"
)
