check_required_packages(
  c("available"),
  "synthetic",
  availability_fn = function(package_name, quietly = TRUE) TRUE
)
try(
  check_required_packages(
    c("missing"),
    "synthetic",
    availability_fn = function(package_name, quietly = TRUE) FALSE
  ),
  silent = TRUE
)
required_packages_for_mode("flow", "png")
required_packages_for_mode("flow", "html")
required_packages_for_mode("pie", "png")
required_packages_for_mode("pathview", "png")
required_packages_for_mode("all", c("png", "html"))
check_flow_html_preflight("flow", "html", pandoc_available_fn = function() TRUE)
check_flow_html_preflight("flow", "png", pandoc_available_fn = function() FALSE)
try(
  check_flow_html_preflight("flow", "html", pandoc_available_fn = function() FALSE),
  silent = TRUE
)

coverage_root <- tempfile("p3_cov_")
dir.create(coverage_root, recursive = TRUE)
on.exit(unlink(coverage_root, recursive = TRUE, force = TRUE), add = TRUE)

plot_object <- ggplot2::ggplot(
  data.frame(x = 1, y = 1),
  ggplot2::aes(x = x, y = y)
) + ggplot2::geom_point()
dimensions <- list(`2x2` = c(width = 2, height = 2))
short_dir <- file.path(coverage_root, "short")
short_files <- save_png_dimensions(
  plot_object,
  short_dir,
  "short_plot",
  dimensions,
  72,
  max_path_length = 1000L
)
assert_output_artifact(short_files[[1L]])
try(assert_output_artifact(file.path(coverage_root, "missing.png")), silent = TRUE)
stable_path_token("logical_plot_2x2.png")
safe_output_component("S1")
safe_output_component("sample/one")
safe_output_component(
  "Chlorocyclohexane and chlorobenzene degradation",
  max_length = 28L
)
try(safe_output_component(""), silent = TRUE)
try(safe_output_component("sample", max_length = 10L), silent = TRUE)
taxonomy_file_stem("taxonomy_global", "percent", "phylum")
taxonomy_file_stem(
  "taxonomy_by_pathway",
  "abund",
  "species",
  "Synthetic pathway"
)
try(taxonomy_file_stem("invalid", "percent", "phylum"), silent = TRUE)

widget_file <- file.path(coverage_root, "widget.html")
save_html_widget(
  plotly::plot_ly(x = 1:2, y = c(1, 2), type = "scatter", mode = "lines"),
  widget_file
)
protected_widget_file <- file.path(coverage_root, "protected_widget.html")
dir.create(file.path(coverage_root, "protected_widget_files"))
try(
  save_html_widget(
    plotly::plot_ly(x = 1:2, y = c(1, 2), type = "scatter", mode = "lines"),
    protected_widget_file
  ),
  silent = TRUE
)

pathview_sqm <- list(functions = list(KEGG = list(tpm = data.frame(
  S1 = c(1, 2),
  S2 = c(3, 4),
  row.names = c("K00001", "K00002")
))))
build_pathview_input_table(pathview_sqm, c("S1", "S2"))
try(build_pathview_input_table(list(), "S1"), silent = TRUE)
fake_pathview_export <- function(
    SQM, pathway_id, count, samples, split_samples, output_dir, output_suffix) {
  writeLines("pathview", file.path(output_dir, "map.png"))
}
export_pathview_isolated(
  fake_pathview_export,
  pathview_sqm,
  "00910",
  c("S1", "S2"),
  file.path(coverage_root, "pathview")
)
try(
  export_pathview_isolated(
    function(...) invisible(NULL),
    pathview_sqm,
    "00910",
    "S1",
    file.path(coverage_root, "empty_pathview")
  ),
  silent = TRUE
)

long_dir <- file.path(coverage_root, "long")
dir.create(long_dir, recursive = TRUE)
long_limit <- nchar(normalizePath(long_dir, winslash = "/")) + 50L
save_png_dimensions(
  plot_object,
  long_dir,
  paste(rep("long", 30L), collapse = "_"),
  dimensions,
  72,
  max_path_length = long_limit
)
try(
  portable_png_output_path(
    long_dir,
    paste(rep("long", 30L), collapse = "_"),
    "2x2",
    max_path_length = nchar(normalizePath(long_dir, winslash = "/")) + 20L
  ),
  silent = TRUE
)

valid_relative <- "taxon/valid.tsv"
valid_target <- file.path(coverage_root, valid_relative)
dir.create(dirname(valid_target), recursive = TRUE)
writeLines("valid", valid_target)
manifest_target_status(
  c(valid_relative, "taxon/missing.tsv", "../outside.tsv", NA_character_, ""),
  coverage_root
)
prune_stale_manifest_targets(tibble::tibble(), coverage_root, "empty")
suppressWarnings(
  prune_stale_manifest_targets(tibble::tibble(marker = "no-path"), coverage_root, "missing-column")
)
valid_row <- tibble::tibble(output_file = valid_relative, marker = "new")
validate_current_manifest_targets(valid_row, coverage_root)
try(
  validate_current_manifest_targets(
    tibble::tibble(output_file = "taxon/missing.tsv"),
    coverage_root
  ),
  silent = TRUE
)

manifest_path <- file.path(coverage_root, "manifest_taxon.tsv")
write_tsv_safe(
  tibble::tibble(
    output_file = c(valid_relative, "taxon/stale.tsv"),
    marker = c("old", "stale")
  ),
  manifest_path
)
suppressWarnings(write_section_manifest(valid_row, coverage_root, "", "manifest_taxon.tsv"))
write_combined_manifest(coverage_root)
merge_section_manifest(
  tibble::tibble(output_file = "new.tsv", marker = "new"),
  tibble::tibble(output_file = "old.tsv", marker = NA)
)
merge_section_manifest(
  tibble::tibble(output_file = "new.tsv", marker = NA),
  tibble::tibble(output_file = "old.tsv", marker = "old")
)
merge_section_manifest(
  tibble::tibble(marker = "new"),
  tibble::tibble(marker = "old")
)

orf_fixture <- tibble::tibble(
  orf_id = c("orf_a", "orf_b", "orf_c"),
  sample = rep("S1", 3L),
  tpm = c(30, 10, 60),
  ko_id = c("K00001", "K00001", "K99999"),
  kegg_function = c("KO one", "KO one", "KO other"),
  ec_codes = c("1.1.1.1", "1.1.1.1", NA_character_),
  phylum = c("Alpha", "Beta", "Gamma")
)
pie_table <- build_pie_chart_table(
  orf_fixture,
  "S1",
  "K00001",
  "phylum",
  2L,
  pathway_name = "Synthetic pathway",
  pathway_id = "12345",
  pathway_selection = "defined",
  pathway_sample_tpm = 100
)
make_pie_plot(pie_table)
build_pie_chart_table(
  orf_fixture,
  "S1",
  "K11111",
  "phylum",
  1L,
  pathway_name = "Synthetic pathway",
  pathway_id = "12345",
  pathway_selection = "defined",
  pathway_sample_tpm = 100
)
