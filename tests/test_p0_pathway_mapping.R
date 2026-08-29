source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")

failures <- character()

record_failure <- function(label, detail = NULL) {
  message <- if (is.null(detail) || !nzchar(detail)) {
    label
  } else {
    paste0(label, ": ", detail)
  }
  failures <<- c(failures, message)
}

expect_true <- function(condition, label, detail = NULL) {
  if (!isTRUE(condition)) {
    record_failure(label, detail)
  }
}

expected_pathways <- c(
  "Carbon fixation in photosynthetic organisms" = "00710",
  "Nitrotoluene degradation" = "00633",
  "Nitrogen metabolism" = "00910"
)
legacy_wrong_ids <- c("00622", "00642", "00643")

fake_sqm <- list(
  misc = list(KEGG_paths = names(expected_pathways))
)

for (canonical_name in names(expected_pathways)) {
  expected_id <- unname(expected_pathways[[canonical_name]])
  curated_rows <- script_env$known_pathways[
    script_env$known_pathways$canonical_pathway_name == canonical_name,
    ,
    drop = FALSE
  ]

  expect_true(
    nrow(curated_rows) == 1L && identical(curated_rows$pathway_id[[1]], expected_id),
    paste0("official curated mapping for '", canonical_name, "'"),
    paste0("expected ", expected_id)
  )

  numeric_resolution <- tryCatch(
    script_env$resolve_pathways(fake_sqm, expected_id),
    error = identity
  )
  if (inherits(numeric_resolution, "error")) {
    record_failure(
      paste0("numeric resolution for ", expected_id),
      conditionMessage(numeric_resolution)
    )
  } else {
    expect_true(
      length(numeric_resolution) == 1L &&
        identical(numeric_resolution[[1]]$pathway_id, expected_id) &&
        identical(numeric_resolution[[1]]$canonical_pathway_name, canonical_name),
      paste0("numeric resolution for ", expected_id),
      "did not return the official ID/name pair"
    )
  }

  name_resolution <- tryCatch(
    script_env$resolve_pathways(fake_sqm, canonical_name),
    error = identity
  )
  if (inherits(name_resolution, "error")) {
    record_failure(
      paste0("canonical-name resolution for '", canonical_name, "'"),
      conditionMessage(name_resolution)
    )
  } else {
    expect_true(
      length(name_resolution) == 1L &&
        identical(name_resolution[[1]]$pathway_id, expected_id) &&
        identical(name_resolution[[1]]$canonical_pathway_name, canonical_name),
      paste0("canonical-name resolution for '", canonical_name, "'"),
      paste0("expected pathway_id ", expected_id)
    )
  }
}

expect_true(
  !any(script_env$known_pathways$pathway_id %in% legacy_wrong_ids),
  "legacy wrong IDs are absent from the curated pathway table",
  paste(legacy_wrong_ids, collapse = ", ")
)

for (legacy_id in legacy_wrong_ids) {
  legacy_resolution <- tryCatch(
    script_env$resolve_pathways(fake_sqm, legacy_id),
    error = identity
  )
  expect_true(
    inherits(legacy_resolution, "error"),
    paste0("legacy wrong ID ", legacy_id, " is rejected")
  )
}

expect_true(
  identical(script_env$pathview_is_exportable("defined", NA_character_), FALSE),
  "Pathview gate rejects defined + NA"
)
expect_true(
  identical(script_env$pathview_is_exportable("top20", NA_character_), FALSE),
  "Pathview gate rejects top20 + NA"
)
expect_true(
  identical(script_env$pathview_is_exportable("defined", "00910"), TRUE),
  "Pathview gate accepts defined + 00910"
)

test_root <- tempfile("p0_pathway_mapping_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

exported_pathway_ids <- character()
fake_export_pathway <- function(
    SQM,
    pathway_id,
    count,
    samples,
    split_samples,
    output_dir,
    output_suffix) {
  exported_pathway_ids <<- c(exported_pathway_ids, pathway_id)
  writeLines("fake pathview output", file.path(output_dir, "map.png"))
}

pathview_result <- tryCatch(
  script_env$run_pathview_mode(
    sqm_object = list(fake = TRUE),
    output_dir = test_root,
    manifest_base_dir = test_root,
    output_manifests = list(pathview = tibble::tibble()),
    script_name = "sqm_plots.R",
    project_dir = "project",
    tax_mode = "prokfilter",
    pathway_name = "Nitrogen metabolism",
    pathway_id = "00910",
    selected_samples = c("S0", "S1"),
    top_n_taxa = 15L,
    top_n_ko = 20L,
    pathway_selection = "defined",
    pathview_sample_modes = "insieme",
    export_pathway_fn = fake_export_pathway
  ),
  error = identity
)

if (inherits(pathview_result, "error")) {
  record_failure("fake Pathview export with 00910", conditionMessage(pathview_result))
} else {
  expect_true(
    identical(exported_pathway_ids, "00910"),
    "fake exportPathway receives 00910"
  )
}

if (length(failures) > 0L) {
  writeLines(
    c("P0 pathway mapping regression failures:", paste0("- ", failures)),
    con = stderr()
  )
  stop(length(failures), " P0 pathway mapping assertions failed.", call. = FALSE)
}

message("PASS: Official KEGG mappings and Pathview export gate are correct")
