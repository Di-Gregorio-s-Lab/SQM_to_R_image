script_env <- if (identical(Sys.getenv("R_COVR"), "true")) environment() else new.env(parent = globalenv())
if (!identical(Sys.getenv("R_COVR"), "true")) sys.source("sqm_plots_lean.R", envir = script_env)

if (!exists("run_pipeline", envir = script_env, mode = "function", inherits = FALSE)) {
  stop("Missing lean public API: run_pipeline", call. = FALSE)
}
run_pipeline <- get("run_pipeline", envir = script_env, inherits = FALSE)

orf_ids <- paste0("orf", 1:4)
sqm <- list(
  orfs = list(
    table = data.frame(
      "KEGG ID" = c("K00001", "K00001;K00002", "K00002", NA),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tax = data.frame(
      phylum = c("Alpha", "Beta", "Alpha", NA),
      row.names = orf_ids,
      check.names = FALSE
    ),
    tpm = data.frame(
      S1 = c(30, 20, 10, 5),
      row.names = orf_ids,
      check.names = FALSE
    )
  ),
  functions = list(
    KEGG = list(
      tpm = data.frame(
        S1 = c(80, 40),
        row.names = c("K00001", "K00002"),
        check.names = FALSE
      )
    )
  ),
  misc = list(
    KEGG_names = c(
      K00001 = "First enzyme [EC:1.1.1.1]",
      K00002 = "Second enzyme [EC:2.2.2.2]"
    ),
    KEGG_paths = c(
      K00001 = "Metabolism; Carbohydrate metabolism; First pathway",
      K00002 = "Metabolism; Energy metabolism; Second pathway"
    )
  )
)

catalog <- data.frame(
  pathway_id = c("00001", "00002"),
  pathway_name = c("First pathway", "Second pathway"),
  pathway_root = rep("Metabolism", 2L),
  pathway_category = c("Carbohydrate metabolism", "Energy metabolism"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
resolved_by_name <- script_env$resolve_defined_pathways("first PATHWAY", catalog)
stopifnot(
  identical(resolved_by_name$pathway_id, "00001"),
  identical(resolved_by_name$pathway_name, "First pathway")
)
alias_catalog <- data.frame(pathway_id = "00720", pathway_name = "Other carbon fixation pathways")
stopifnot(identical(
  script_env$catalog_id("Carbon fixation pathways in prokaryotes", alias_catalog),
  "00720"
))
pathway_kos <- list(`00001` = "K00001", `00002` = "K00002")
kgml_loader <- function(pathway_id) pathway_kos[[as.character(pathway_id)]]

fake_plot_taxonomy <- function(SQM, rank, samples, rescale, ...) {
  stopifnot(identical(rank, "phylum"), identical(samples, "S1"), identical(rescale, FALSE))
  data <- data.frame(phylum = c("Alpha", "Beta"), TPM = c(80, 40))
  ggplot2::ggplot(data, ggplot2::aes(phylum, TPM)) + ggplot2::geom_col()
}

fake_export_pathway <- function(SQM, pathway_id, samples, output_dir, ...) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0("ko", pathway_id, ".png"))
  grDevices::png(path, width = 320, height = 320)
  graphics::plot.new()
  grDevices::dev.off()
  path
}

make_config <- function(output_dir, mode = "huge", plan_only = FALSE) list(
  output_dir = output_dir,
  mode = mode,
  samples = "S1",
  ranks = "phylum",
  counts = "abund",
  formats = "png",
  dimensions = list(list(name = "tiny", width = 2, height = 2, dpi = 72)),
  top_n_pathways = 20L,
  top_n_ko = 2L,
  top_n_taxa = 2L,
  workers = 1L,
  plan_only = plan_only,
  refresh_kegg = FALSE,
  selection_modes = c("defined", "top20"),
  defined_pathways = c("00001", "00002")
)

invoke <- function(config) run_pipeline(
  sqm = sqm,
  config = config,
  catalog = catalog,
  kgml_loader = kgml_loader,
  plot_taxonomy_fn = fake_plot_taxonomy,
  export_pathway_fn = fake_export_pathway
)

top20_columns <- c(
  "rank", "pathway_id", "pathway_name", "pathway_root",
  "pathway_category", "total_tpm", "samples", "filtered_taxon",
  "filtered_taxon_rank"
)
graphic_roots <- c(
  "funz", "flowplot", "taxonomy_global", "taxonomy_by_pathway", "pie", "pathview"
)

sqm$orfs$table$KEGGPATH <- rep("Human Diseases; Synthetic; Wrong pathway", nrow(sqm$orfs$table))
plan_dir <- tempfile("lean_pipeline_plan_")
dir.create(plan_dir)
on.exit(unlink(plan_dir, recursive = TRUE, force = TRUE), add = TRUE)
invoke(make_config(plan_dir, plan_only = TRUE))
stopifnot(
  file.exists(file.path(plan_dir, "errors.tsv")),
  identical(readLines(file.path(plan_dir, "errors.tsv"), warn = FALSE), "task_id\tmode\tpathway_id\terror")
)

planned <- utils::read.delim(
  file.path(plan_dir, "top20.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  colClasses = c(pathway_id = "character")
)
stopifnot(
  identical(names(planned), top20_columns),
  identical(planned$rank, 1:2),
  identical(as.character(planned$pathway_id), c("00001", "00002")),
  isTRUE(all.equal(planned$total_tpm, c(80, 40), tolerance = 1e-12)),
  !any(dir.exists(file.path(plan_dir, graphic_roots)))
)

empty_dir <- tempfile("lean_pipeline_empty_")
dir.create(empty_dir)
on.exit(unlink(empty_dir, recursive = TRUE, force = TRUE), add = TRUE)
empty_config <- make_config(empty_dir)
empty_config$selection_modes <- "defined"
empty_config$defined_pathways <- "00003"
empty_result <- run_pipeline(
  sqm, empty_config,
  rbind(catalog, data.frame(
    pathway_id = "00003", pathway_name = "Empty pathway", pathway_root = "Metabolism",
    pathway_category = "Synthetic", stringsAsFactors = FALSE
  )),
  function(...) "K00003", fake_plot_taxonomy, fake_export_pathway
)
stopifnot(nrow(empty_result$errors) == 0L)

output_dir <- tempfile("lean_pipeline_all_")
dir.create(output_dir)
on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
pipeline_result <- invoke(make_config(output_dir))

missing_roots <- graphic_roots[!dir.exists(file.path(output_dir, graphic_roots))]
if (length(missing_roots)) stop(
  "Missing graphic roots: ", paste(missing_roots, collapse = ", "),
  "; errors: ", paste(paste(pipeline_result$errors$task_id, pipeline_result$errors$error, sep = "="), collapse = " | "), call. = FALSE
)
stopifnot(length(list.files(file.path(output_dir, "pathview"), pattern = "\\.tsv$", recursive = TRUE)) > 0L)
stopifnot(length(list.files(file.path(output_dir, "funz"), pattern = "enzim", recursive = TRUE, ignore.case = TRUE)) > 0L)
stopifnot(length(list.files(file.path(output_dir, "funz", "enzimi", "separato"), pattern = "\\.png$", recursive = TRUE)) > 0L)

manifests <- sort(basename(list.files(
  output_dir,
  pattern = "^manifest_.*\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)))
stopifnot(identical(manifests, sort(c(
  "manifest_flow.tsv", "manifest_funz.tsv", "manifest_pie.tsv"
))))

normal_dir <- tempfile("lean_pipeline_normal_")
dir.create(normal_dir)
on.exit(unlink(normal_dir, recursive = TRUE, force = TRUE), add = TRUE)
normal_result <- invoke(make_config(normal_dir, mode = "normal"))
normal_roots <- setdiff(graphic_roots, "pie")
stopifnot(
  nrow(normal_result$errors) == 0L,
  all(dir.exists(file.path(normal_dir, normal_roots))),
  !dir.exists(file.path(normal_dir, "pie")),
  file.exists(file.path(normal_dir, "top20.tsv")),
  !any(basename(list.dirs(normal_dir, recursive = TRUE, full.names = TRUE)) == "top20")
)

errors_path <- file.path(output_dir, "errors.tsv")
errors <- utils::read.delim(errors_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(
  file.exists(errors_path),
  length(readLines(errors_path, warn = FALSE)) == 1L,
  nrow(errors) == 0L,
  ncol(errors) > 0L
)

sentinel <- file.path(output_dir, "keep-me.txt")
writeLines("owned by user", sentinel)
writeLines("stale", file.path(output_dir, "top20.tsv"))
invoke(make_config(output_dir))

rewritten <- utils::read.delim(
  file.path(output_dir, "top20.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  colClasses = c(pathway_id = "character")
)
stopifnot(
  identical(readLines(sentinel, warn = FALSE), "owned by user"),
  identical(names(rewritten), top20_columns),
  identical(as.character(rewritten$pathway_id), c("00001", "00002"))
)

message("PASS: lean pipeline public contract is deterministic and network-free")
