script_env <- if (identical(Sys.getenv("R_COVR"), "true")) environment() else new.env(parent = globalenv())
if (!identical(Sys.getenv("R_COVR"), "true")) sys.source("sqm_plots_lean.R", envir = script_env)

need <- function(name) {
  if (!exists(name, envir = script_env, mode = "function", inherits = FALSE)) {
    stop("Missing lean helper: ", name, call. = FALSE)
  }
  get(name, envir = script_env, inherits = FALSE)
}

parse_args <- need("parse_args")
default_workers <- need("default_workers")
equals <- parse_args(c(
  "--project_dir=project", "--output_dir=output", "--mode=funz",
  "--workers=2", "--plan_only", "--refresh_kegg"
))
separate <- parse_args(c(
  "--project_dir", "project", "--output_dir", "output", "--mode", "funz",
  "--workers", "2", "--plan_only", "--refresh_kegg"
))
fields <- c(
  "project_dir", "output_dir", "mode", "workers", "plan_only", "refresh_kegg"
)
stopifnot(
  identical(equals[fields], separate[fields]),
  identical(equals$workers, 2L),
  isTRUE(equals$plan_only),
  isTRUE(equals$refresh_kegg)
)
stopifnot(identical(
  parse_args(c("--project_dir", "project", "--output_dir", "output", "--mode", "funz"))$workers,
  default_workers()
))

stopifnot(identical(
  need("extract_ko_ids")("K00002; K00001*; K00002; invalid"),
  c("K00002", "K00001")
))
stopifnot(identical(
  need("extract_ec_codes")(c(
    "enzyme [EC:2.2.2.2 1.1.1.-]",
    "no EC"
  )),
  c("2.2.2.2;1.1.1.-", NA_character_)
))

top20_columns <- c(
  "rank", "pathway_id", "pathway_name", "pathway_root", "pathway_category",
  "total_tpm", "samples", "filtered_taxon", "filtered_taxon_rank"
)
ranked <- need("rank_top20")(
  data.frame(
    pathway_id = c("00003", "00002", "00001"),
    pathway_name = c("Low", "Tie B", "Tie A"),
    pathway_root = rep("Metabolism", 3L),
    pathway_category = rep("Synthetic", 3L),
    total_tpm = c(5, 10, 10),
    stringsAsFactors = FALSE
  ),
  n = 2L,
  samples = c("S2", "S1")
)
stopifnot(
  identical(names(ranked), top20_columns),
  identical(ranked$rank, 1:2),
  identical(ranked$pathway_id, c("00001", "00002")),
  identical(ranked$total_tpm, c(10, 10)),
  identical(ranked$samples, rep("S2,S1", 2L))
)

orf_table <- data.frame(
  "KEGG ID" = c("K00001", "K00001;K99999*"),
  row.names = c("orf_single", "orf_multi"),
  check.names = FALSE
)
orf_tpm <- data.frame(
  S1 = c(20, 30),
  row.names = row.names(orf_table),
  check.names = FALSE
)
official_ko_tpm <- data.frame(
  S1 = c(35, 15),
  row.names = c("K00001", "K99999"),
  check.names = FALSE
)
allocated <- need("allocate_ko_tpm")(
  orf_table,
  orf_tpm,
  official_ko_tpm,
  pathway_ko_ids = "K00001",
  samples = "S1",
  context = "defined:00001 Test pathway"
)
allocated <- allocated[order(allocated$orf_id), , drop = FALSE]
stopifnot(
  identical(as.character(allocated$orf_id), c("orf_multi", "orf_single")),
  identical(as.character(allocated$ko_id), c("K00001", "K00001")),
  isTRUE(all.equal(as.numeric(allocated$tpm), c(15, 20), tolerance = 1e-12)),
  isTRUE(all.equal(sum(allocated$tpm), 35, tolerance = 1e-12))
)

mismatched_ko_tpm <- official_ko_tpm
mismatched_ko_tpm["K00001", "S1"] <- 70
mismatch_error <- tryCatch(
  need("allocate_ko_tpm")(
    orf_table, orf_tpm, mismatched_ko_tpm,
    pathway_ko_ids = "K00001", samples = "S1", context = "defined:00001 Test pathway"
  ),
  error = conditionMessage
)
stopifnot(
  grepl("defined:00001 Test pathway", mismatch_error, fixed = TRUE),
  grepl("S1/K00001", mismatch_error, fixed = TRUE),
  grepl("observed=35", mismatch_error, fixed = TRUE),
  grepl("target=70", mismatch_error, fixed = TRUE),
  grepl("factor=2", mismatch_error, fixed = TRUE)
)

rounded_ko_tpm <- official_ko_tpm
rounded_ko_tpm["K00001", "S1"] <- 35 + 1e-8
rounded <- need("allocate_ko_tpm")(
  orf_table, orf_tpm, rounded_ko_tpm,
  pathway_ko_ids = "K00001", samples = "S1", context = "defined:00001 Test pathway"
)
stopifnot(isTRUE(all.equal(sum(rounded$tpm), 35, tolerance = 1e-12)))

local({
  output_dir <- tempfile("sqm_plots_lean_")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  global_path <- need("top20_path")(output_dir)
  taxon_path <- need("top20_path")(output_dir, "phylum", "Alpha beta")
  stopifnot(
    identical(global_path, file.path(output_dir, "top20.tsv")),
    identical(
      taxon_path,
      file.path(output_dir, "taxon_filter", "phylum", "Alpha_beta", "top20.tsv")
    )
  )
  need("write_top20")(taxon_path, ranked)
  written <- read.delim(taxon_path, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(
    identical(names(written), top20_columns),
    identical(sprintf("%05d", written$pathway_id), ranked$pathway_id),
    isTRUE(all.equal(written$total_tpm, ranked$total_tpm, tolerance = 1e-12))
  )
})

tasks <- list(
  list(id = "b", value = "B"),
  list(id = "a", value = "A"),
  list(id = "c", value = "C")
)
task_result <- need("run_tasks")(
  tasks,
  function(task) {
    if (task$id == "b") stop("synthetic failure", call. = FALSE)
    task$value
  },
  workers = 1L
)
stopifnot(
  identical(names(task_result), c("results", "errors")),
  identical(names(task_result$results), c("a", "c")),
  identical(unname(unlist(task_result$results)), c("A", "C")),
  is.data.frame(task_result$errors),
  identical(task_result$errors$task_id, "b"),
  grepl("synthetic failure", task_result$errors$error, fixed = TRUE)
)

parallel_result <- need("run_tasks")(
  tasks,
  function(task) {
    if (task$id == "b") stop("synthetic failure", call. = FALSE)
    task$value
  },
  workers = 2L
)
stopifnot(
  identical(parallel_result$results, task_result$results),
  identical(parallel_result$errors, task_result$errors)
)

local({
  render_value <- function(task) task$value
  render_task <- function(task) render_value(task)
  exported_result <- need("run_tasks")(tasks, render_task, workers = 2L)
  stopifnot(identical(exported_result$results, list(a = "A", b = "B", c = "C")))
})

evalq({
  cluster_render_value <- function(task) task$value
  cluster_render_task <- function(task) cluster_render_value(task)
}, envir = script_env)
script_env_result <- need("run_tasks")(
  tasks, get("cluster_render_task", envir = script_env), workers = 2L
)
stopifnot(identical(script_env_result$results, list(a = "A", b = "B", c = "C")))

message("PASS: lean core contracts are deterministic and network-free")
