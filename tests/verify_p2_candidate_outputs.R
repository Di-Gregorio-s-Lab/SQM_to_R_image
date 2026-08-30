args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript tests/verify_p2_candidate_outputs.R FLOW_DIR FUNZ_DIR PIE_DIR",
    call. = FALSE
  )
}

flow_root <- args[[1L]]
funz_root <- args[[2L]]
pie_root <- args[[3L]]
for (candidate_root in args) {
  if (!dir.exists(candidate_root)) {
    stop("Candidate directory does not exist: ", candidate_root, call. = FALSE)
  }
}

audit_columns <- c(
  "input_orf_count",
  "excluded_orfs_without_ko",
  "multi_ko_orf_count",
  "orf_ko_association_count",
  "multi_ko_policy",
  "ko_denominator_basis"
)

assert_nonempty_files <- function(paths, label) {
  if (length(paths) == 0L || anyNA(file.info(paths)$size) || any(file.info(paths)$size <= 0)) {
    stop(label, " files are missing or empty.", call. = FALSE)
  }
}

assert_complete_audit <- function(manifest, label) {
  missing_columns <- setdiff(audit_columns, colnames(manifest))
  if (length(missing_columns) > 0L || anyNA(manifest[audit_columns])) {
    stop(label, " manifest has incomplete KO provenance.", call. = FALSE)
  }
  if (any(manifest$multi_ko_policy != "full_tpm_per_ko") ||
      any(manifest$ko_denominator_basis != "expanded_orf_sample_ko_tpm")) {
    stop(label, " manifest reports an unexpected KO expansion policy.", call. = FALSE)
  }
}

# FLOW: three selected samples, unique Sankey keys, complete mass, and both
# reserved aggregate categories represented without conflation.
flow_files <- list.files(
  flow_root,
  pattern = "_data\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(flow_files) != 3L) {
  stop("FLOW candidate must contain three sample TSV files.", call. = FALSE)
}
flow <- dplyr::bind_rows(lapply(
  flow_files,
  readr::read_tsv,
  show_col_types = FALSE
))
if (anyDuplicated(flow[c("pathway", "rank", "sample", "taxon", "KO")]) > 0L) {
  stop("FLOW candidate contains duplicate edge keys.", call. = FALSE)
}
flow_sums <- flow |>
  dplyr::group_by(.data$sample) |>
  dplyr::summarise(total = sum(.data$flow_percent), .groups = "drop")
if (any(abs(flow_sums$total - 100) > 1e-6)) {
  stop("FLOW candidate percentages do not sum to 100.", call. = FALSE)
}
if (!all(c("Unclassified", "Other") %in% unique(flow$taxon))) {
  stop("FLOW candidate does not keep Unclassified distinct from Other.", call. = FALSE)
}
flow_manifest <- readr::read_tsv(
  file.path(flow_root, "flowplot", "manifest_flow.tsv"),
  show_col_types = FALSE
)
if (any(flow_manifest$pathway_id != "00361") || any(flow_manifest$top_n_taxa != 3L)) {
  stop("FLOW candidate manifest does not describe pathway 00361 Top 3.", call. = FALSE)
}
assert_complete_audit(flow_manifest, "FLOW")
flow_png <- list.files(flow_root, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
if (length(flow_png) != 3L) {
  stop("FLOW candidate must contain three PNG files.", call. = FALSE)
}
assert_nonempty_files(flow_png, "FLOW PNG")

# FUNZ: one row per sample/KO, exact denominators, complete multi-EC lookup,
# and no artificial Other row because top_n_ko=999 includes every KO.
funz_files <- list.files(
  funz_root,
  pattern = "barplot_ko_data\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(funz_files) != 1L) {
  stop("FUNZ candidate must contain one pathway KO TSV.", call. = FALSE)
}
funz <- readr::read_tsv(funz_files[[1L]], show_col_types = FALSE)
if (anyDuplicated(funz[c("sample", "ko_id")]) > 0L) {
  stop("FUNZ candidate contains duplicate sample/KO keys.", call. = FALSE)
}
if (any(as.character(funz$ko_id) == "Other", na.rm = TRUE)) {
  stop("FUNZ candidate unexpectedly collapsed KOs with top_n_ko=999.", call. = FALSE)
}
funz_sums <- funz |>
  dplyr::filter(.data$status == "ok") |>
  dplyr::group_by(.data$sample) |>
  dplyr::summarise(
    percent = sum(.data$sample_pathway_percent),
    tpm = sum(.data$tpm),
    denominator = dplyr::first(.data$denominator),
    .groups = "drop"
  )
if (any(abs(funz_sums$percent - 100) > 1e-6) ||
    any(abs(funz_sums$tpm - funz_sums$denominator) > 1e-10)) {
  stop("FUNZ candidate percentages or denominators are inconsistent.", call. = FALSE)
}
multi_ec_funz <- funz |>
  dplyr::filter(!is.na(.data$ec_codes), grepl(";", .data$ec_codes, fixed = TRUE)) |>
  dplyr::distinct(.data$ko_id, .data$ec_codes)
if (nrow(multi_ec_funz) != 3L) {
  stop("FUNZ candidate did not expose the three observed multi-EC KOs.", call. = FALSE)
}
funz_manifest <- readr::read_tsv(
  file.path(funz_root, "funz", "manifest_funz.tsv"),
  show_col_types = FALSE
) |>
  dplyr::filter(.data$mode == "funz")
if (any(funz_manifest$pathway_id != "00710") || any(funz_manifest$top_n_ko != 999L)) {
  stop("FUNZ candidate manifest does not describe pathway 00710/top_n_ko=999.", call. = FALSE)
}
assert_complete_audit(funz_manifest, "FUNZ")
funz_png <- list.files(
  file.path(funz_root, "funz", "pathway"),
  pattern = "barplot_ko_.*\\.png$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(funz_png) != 1L) {
  stop("FUNZ candidate must contain one pathway PNG.", call. = FALSE)
}
assert_nonempty_files(funz_png, "FUNZ PNG")

# PIE: each KO table preserves full mass, and the observed K10679 multi-EC
# association is identical in its TSV and manifest.
pie_files <- list.files(
  pie_root,
  pattern = "_data\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(pie_files) != 11L) {
  stop("PIE candidate must contain the 11 observed KO TSV files.", call. = FALSE)
}
pie <- dplyr::bind_rows(lapply(
  pie_files,
  readr::read_tsv,
  show_col_types = FALSE
))
pie_sums <- pie |>
  dplyr::group_by(.data$sample, .data$ko_id) |>
  dplyr::summarise(total = sum(.data$pct), .groups = "drop")
if (any(abs(pie_sums$total - 1) > 1e-10)) {
  stop("PIE candidate fractions do not sum to 1.", call. = FALSE)
}
expected_multi_ec <- "1.-.-.-;1.5.1.34"
if (!any(pie$ko_id == "K10679" & pie$ec_codes == expected_multi_ec)) {
  stop("PIE candidate TSV truncated the K10679 EC list.", call. = FALSE)
}
pie_manifest <- readr::read_tsv(
  file.path(pie_root, "pie", "manifest_pie.tsv"),
  show_col_types = FALSE
)
if (any(pie_manifest$pathway_id != "00633") ||
    any(pie_manifest$samples != "S13_1_8") ||
    !any(pie_manifest$ko_id == "K10679" & pie_manifest$ec_code == expected_multi_ec)) {
  stop("PIE candidate manifest does not retain the requested scope or full EC list.", call. = FALSE)
}
assert_complete_audit(pie_manifest, "PIE")
pie_png <- list.files(pie_root, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
if (length(pie_png) != 11L) {
  stop("PIE candidate must contain 11 PNG files.", call. = FALSE)
}
assert_nonempty_files(pie_png, "PIE PNG")

message(
  "PASS: P2 candidate outputs verified | FLOW rows=", nrow(flow),
  " | FUNZ rows=", nrow(funz),
  " | FUNZ multi-EC KOs=", nrow(multi_ec_funz),
  " | PIE KOs=", nrow(pie_sums)
)
