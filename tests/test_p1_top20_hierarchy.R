source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")

pathway_fields <- c(
  paste(
    "Metabolism; Carbohydrate metabolism; Glycolysis / Gluconeogenesis",
    "Metabolism; Carbohydrate metabolism; Glycolysis / Gluconeogenesis",
    sep = " | "
  ),
  "Genetic Information Processing; Translation; Ribosome",
  "Environmental Information Processing; Signal transduction; Two-component system",
  "Cellular Processes; Cell growth and death; Cell cycle - prokaryotes",
  "Organismal Systems; Immune system; Toll and Imd signaling pathway",
  "Human Diseases; Cancer: overview; Pathways in cancer",
  "Metabolism; Amino acid metabolism; Aardvark pathway",
  "Metabolism; Energy metabolism; Zebra pathway",
  "Brite Hierarchies; Protein families: signaling and cellular processes; Transporters",
  "Not Included in Pathway or Brite; Poorly characterized; Function unknown",
  "Unknown KEGG root; Unknown category; Unknown leaf",
  "Metabolism; Incomplete hierarchy"
)

orf_ids <- sprintf("orf_%02d", seq_along(pathway_fields))
fake_sqm <- list(
  orfs = list(
    table = data.frame(
      KEGGPATH = pathway_fields,
      row.names = orf_ids,
      check.names = FALSE
    ),
    tpm = data.frame(
      S0 = c(10, 60, 50, 40, 30, 20, 15, 15, 1000, 900, 800, 700),
      S1 = rep(0, length(pathway_fields)),
      row.names = orf_ids,
      check.names = FALSE
    )
  )
)

top_pathways <- script_env$select_top_pathways(
  fake_sqm,
  selected_samples = c("S0", "S1"),
  pathway_top_n = 20L
)

expected <- data.frame(
  canonical_pathway_name = c(
    "Ribosome",
    "Two-component system",
    "Cell cycle - prokaryotes",
    "Toll and Imd signaling pathway",
    "Pathways in cancer",
    "Aardvark pathway",
    "Zebra pathway",
    "Glycolysis / Gluconeogenesis"
  ),
  pathway_root = c(
    "Genetic Information Processing",
    "Environmental Information Processing",
    "Cellular Processes",
    "Organismal Systems",
    "Human Diseases",
    "Metabolism",
    "Metabolism",
    "Metabolism"
  ),
  pathway_category = c(
    "Translation",
    "Signal transduction",
    "Cell growth and death",
    "Immune system",
    "Cancer: overview",
    "Amino acid metabolism",
    "Energy metabolism",
    "Carbohydrate metabolism"
  ),
  total_tpm = c(60, 50, 40, 30, 20, 15, 15, 10),
  stringsAsFactors = FALSE
)

failures <- character()
expect_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    failures <<- c(failures, message)
  }
}

selected_names <- vapply(
  top_pathways,
  `[[`,
  character(1),
  "canonical_pathway_name"
)

expect_true(
  identical(selected_names, expected$canonical_pathway_name),
  paste0(
    "top20 must contain only valid three-level KEGG PATHWAY leaves in deterministic order; got: ",
    paste(selected_names, collapse = ", ")
  )
)
expect_true(
  !any(c("Transporters", "Function unknown", "Unknown leaf", "Incomplete hierarchy") %in% selected_names),
  "BRITE, Not Included, unknown-root, and malformed entries must be excluded"
)

required_fields <- c(
  "pathway_root",
  "pathway_category",
  "canonical_pathway_name",
  "total_tpm"
)
has_required_fields <- all(vapply(
  top_pathways,
  function(entry) all(required_fields %in% names(entry)),
  logical(1)
))
expect_true(
  has_required_fields,
  paste("each selected pathway must expose", paste(required_fields, collapse = ", "))
)

if (has_required_fields && identical(selected_names, expected$canonical_pathway_name)) {
  expect_true(
    identical(
      vapply(top_pathways, `[[`, character(1), "pathway_root"),
      expected$pathway_root
    ),
    "all six accepted KEGG PATHWAY roots must be retained"
  )
  expect_true(
    identical(
      vapply(top_pathways, `[[`, character(1), "pathway_category"),
      expected$pathway_category
    ),
    "the second hierarchy level must be retained as pathway_category"
  )
  expect_true(
    identical(
      vapply(top_pathways, `[[`, numeric(1), "total_tpm"),
      expected$total_tpm
    ),
    "pathway TPM totals must be exact and duplicate ORF/pathway membership must not be counted twice"
  )
}

ambiguous_sqm <- list(
  orfs = list(
    table = data.frame(
      KEGGPATH = paste(
        "Metabolism; Category A; Shared pathway",
        "Human Diseases; Category B; Shared pathway",
        sep = " | "
      ),
      row.names = "orf_ambiguous",
      check.names = FALSE
    ),
    tpm = data.frame(
      S0 = 5,
      row.names = "orf_ambiguous",
      check.names = FALSE
    )
  )
)

ambiguous_error <- tryCatch(
  {
    script_env$select_top_pathways(
      ambiguous_sqm,
      selected_samples = "S0",
      pathway_top_n = 1L
    )
    NULL
  },
  error = function(error) error$message
)

expect_true(
  !is.null(ambiguous_error) &&
    grepl("Shared pathway", ambiguous_error, fixed = TRUE) &&
    grepl("hierarch|incompat|ambig", ambiguous_error, ignore.case = TRUE),
  "a pathway leaf assigned to incompatible hierarchies must fail with a controlled diagnostic"
)

if (length(failures) > 0L) {
  stop(
    paste(c("P1-01 hierarchy expectations failed:", paste0("- ", failures)), collapse = "\n"),
    call. = FALSE
  )
}

message("PASS: top20 accepts only unambiguous three-level KEGG PATHWAY entries")
