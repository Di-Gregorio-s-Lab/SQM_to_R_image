script_env <- if (identical(Sys.getenv("R_COVR"), "true")) environment() else new.env(parent = globalenv())
if (!identical(Sys.getenv("R_COVR"), "true")) sys.source("sqm_plots_lean.R", envir = script_env)

need <- function(name) {
  if (!exists(name, envir = script_env, mode = "function", inherits = FALSE)) {
    stop("Missing lean cache helper: ", name, call. = FALSE)
  }
  get(name, envir = script_env, inherits = FALSE)
}

expect_error <- function(expr) {
  failed <- FALSE
  tryCatch(force(expr), error = function(...) failed <<- TRUE)
  stopifnot(failed)
}

output_dir <- tempfile("lean_kegg_cache_")
dir.create(output_dir)
on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

get_kegg_catalog <- need("get_kegg_catalog")
catalog_calls <- 0L
catalog_downloader <- function() {
  catalog_calls <<- catalog_calls + 1L
  data.frame(
    pathway_id = "00010",
    pathway_name = paste("Catalog", catalog_calls),
    stringsAsFactors = FALSE
  )
}

catalog_path <- file.path(output_dir, "_cache", "kegg", "pathway_catalog.tsv")
catalog <- get_kegg_catalog(output_dir, refresh = FALSE, catalog_downloader)
stopifnot(
  catalog_calls == 1L,
  file.exists(catalog_path),
  identical(as.character(catalog$pathway_id), "00010"),
  identical(as.character(catalog$pathway_name), "Catalog 1")
)

catalog_hit <- get_kegg_catalog(output_dir, refresh = FALSE, catalog_downloader)
stopifnot(
  catalog_calls == 1L,
  identical(as.character(catalog_hit$pathway_id), "00010"),
  identical(as.character(catalog_hit$pathway_name), "Catalog 1")
)

writeLines("corrupt", catalog_path)
catalog_repaired <- get_kegg_catalog(output_dir, refresh = FALSE, catalog_downloader)
stopifnot(
  catalog_calls == 2L,
  identical(as.character(catalog_repaired$pathway_name), "Catalog 2"),
  identical(as.character(utils::read.delim(
    catalog_path,
    colClasses = "character",
    check.names = FALSE
  )$pathway_id), "00010")
)

catalog_refreshed <- get_kegg_catalog(output_dir, refresh = TRUE, catalog_downloader)
stopifnot(
  catalog_calls == 3L,
  identical(as.character(catalog_refreshed$pathway_name), "Catalog 3"),
  grepl("Catalog 3", paste(readLines(catalog_path, warn = FALSE), collapse = "\n"), fixed = TRUE)
)

get_pathway_kos <- need("get_pathway_kos")
kgml_calls <- 0L
kgml_downloader <- function(pathway_id) {
  kgml_calls <<- kgml_calls + 1L
  sprintf(
    '<pathway name="path:ko%s"><entry type="ortholog" name="ko:K%05d"/></pathway>',
    pathway_id,
    kgml_calls
  )
}
kgml_parser <- function(path) {
  xml <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!grepl("^<pathway .*?</pathway>$", xml)) stop("Malformed KGML", call. = FALSE)
  matches <- regmatches(xml, gregexpr("K[0-9]{5}", xml))[[1L]]
  if (identical(matches, character(0))) stop("KGML contains no KO", call. = FALSE)
  unique(matches)
}

kgml_path <- file.path(output_dir, "_cache", "kegg", "ko00010.xml")
kos <- get_pathway_kos(output_dir, "00010", refresh = FALSE, kgml_downloader, kgml_parser)
stopifnot(kgml_calls == 1L, file.exists(kgml_path), identical(kos, "K00001"))

kos_hit <- get_pathway_kos(output_dir, "00010", refresh = FALSE, kgml_downloader, kgml_parser)
stopifnot(kgml_calls == 1L, identical(kos_hit, "K00001"))

writeLines("<broken>", kgml_path)
kos_repaired <- get_pathway_kos(output_dir, "00010", refresh = FALSE, kgml_downloader, kgml_parser)
stopifnot(kgml_calls == 2L, identical(kos_repaired, "K00002"))

kos_refreshed <- get_pathway_kos(output_dir, "00010", refresh = TRUE, kgml_downloader, kgml_parser)
stopifnot(
  kgml_calls == 3L,
  identical(kos_refreshed, "K00003"),
  grepl("K00003", paste(readLines(kgml_path, warn = FALSE), collapse = "\n"), fixed = TRUE)
)

calls_before_invalid_ids <- kgml_calls
for (invalid_id in list("10", "ko00010", "abcde", "000000", NA_character_)) {
  expect_error(get_pathway_kos(
    output_dir,
    invalid_id,
    refresh = FALSE,
    kgml_downloader,
    kgml_parser
  ))
}
stopifnot(kgml_calls == calls_before_invalid_ids)

message("PASS: lean KEGG cache contracts are deterministic and network-free")
