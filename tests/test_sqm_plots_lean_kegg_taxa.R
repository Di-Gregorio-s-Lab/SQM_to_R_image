#!/usr/bin/env Rscript

# Regression checks for the per-pathway, per-sample KEGG taxon barplots.
run_tests <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  test_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else
    "SQM_to_R_image_dev/tests/test_sqm_plots_lean_kegg_taxa.R"
  project_dir <- dirname(dirname(normalizePath(test_file, mustWork = TRUE)))
  source(file.path(project_dir, "sqm_plots_lean.R"), local = TRUE)

  assert_equal <- function(actual, expected, label) {
    if (!isTRUE(all.equal(actual, expected, check.attributes = FALSE))) {
      stop(label, ": expected ", paste(expected, collapse = ", "),
           "; got ", paste(actual, collapse = ", "), call. = FALSE)
    }
  }
  value <- function(table, sample, taxon) {
    hit <- table$sample == sample & table$taxon == taxon
    if (!any(hit)) return(0)
    sum(table$tpm[hit])
  }
  make_allocated <- function(sample, taxon, tpm) {
    data.frame(orf_id = paste0("orf", seq_along(tpm)), ko_id = "K00001",
               sample = sample, genus = taxon, tpm = tpm,
               stringsAsFactors = FALSE)
  }

  # Selection is made independently for each sample. A taxon selected in one
  # sample remains visible with its actual TPM in the other sample.
  allocated <- make_allocated(
    c("large", "large", "large", "small", "small", "small", "small"),
    c("Alpha", "Beta", "Unclassified", "Alpha", "Gamma", "Delta", "Minor"),
    c(900, 99, 1, 0.5, 60, 39, 0.5)
  )
  table <- build_kegg_taxa_pathway_table(
    allocated, "00361", "Chlorocyclohexane degradation",
    c("large", "small", "empty"), "genus", top_n = 2L, min_percent = 1
  )
  stopifnot(all(c("pathway_id", "pathway_name", "sample", "taxon", "tpm",
                  "pathway_tpm", "percent", "rank") %in% names(table)))
  stopifnot(identical(unique(table$sample), c("large", "small", "empty")))
  stopifnot(all(table$pathway_id == "00361"), all(table$rank == "genus"))
  assert_equal(value(table, "large", "Alpha"), 900, "large sample's principal taxon")
  assert_equal(value(table, "small", "Gamma"), 60, "small sample's principal taxon")
  assert_equal(value(table, "small", "Delta"), 39, "second principal taxon")
  assert_equal(value(table, "small", "Alpha"), 0.5, "union keeps a subthreshold taxon")
  assert_equal(value(table, "large", "Unclassified"), 1, "unclassified separate")
  assert_equal(value(table, "small", "Other"), 0.5, "other collects unselected")
  assert_equal(value(table, "empty", "Other"), 0, "zero TPM bar")
  for (sample in c("large", "small", "empty")) {
    rows <- table[table$sample == sample, , drop = FALSE]
    expected <- switch(sample, large = 1000, small = 100, empty = 0)
    assert_equal(sum(rows$tpm), expected, paste("TPM conserved for", sample))
    stopifnot(all(rows$pathway_tpm == expected))
    if (expected > 0) assert_equal(sum(rows$percent), 100, paste("percent total", sample))
  }

  # The 1% cutoff is inclusive; a classified taxon below it joins Other.
  threshold <- make_allocated("S", c("A", "AtCutoff", "Below", "Unclassified"),
                              c(97.51, 1, 0.49, 1))
  threshold_table <- build_kegg_taxa_pathway_table(
    threshold, "00625", "Chloroalkane degradation", "S", "genus",
    top_n = 10L, min_percent = 1
  )
  assert_equal(value(threshold_table, "S", "AtCutoff"), 1, "inclusive 1% threshold")
  assert_equal(value(threshold_table, "S", "Other"), 0.49, "below 1% goes to Other")
  assert_equal(value(threshold_table, "S", "Unclassified"), 1, "Unclassified retained")

  # The cap is per bar, so union across samples may exceed ten distinct taxa.
  taxa_a <- sprintf("A%02d", seq_len(11L))
  taxa_b <- sprintf("B%02d", seq_len(11L))
  many <- make_allocated(c(rep("S1", 11L), rep("S2", 11L)),
                         c(taxa_a, taxa_b), rep(10, 22L))
  many_table <- build_kegg_taxa_pathway_table(
    many, "00361", "Chlorocyclohexane degradation", c("S1", "S2"),
    "genus", top_n = 10L, min_percent = 1
  )
  visible <- setdiff(unique(many_table$taxon), c("Other", "Unclassified"))
  stopifnot(length(visible) == 20L)
  assert_equal(value(many_table, "S1", "Other"), 10, "cap leaves one taxon in Other S1")
  assert_equal(value(many_table, "S2", "Other"), 10, "cap leaves one taxon in Other S2")
  assert_equal(sum(many_table$tpm[many_table$sample == "S1"]), 110, "many S1 total")
  assert_equal(sum(many_table$tpm[many_table$sample == "S2"]), 110, "many S2 total")

  # Colors are stable for the same taxon across pathway pages and input order.
  colored <- rbind(table, transform(table, pathway_id = "00625",
                                     pathway_name = "Chloroalkane degradation"))
  palette <- build_kegg_taxa_color_map(colored)
  reversed_palette <- build_kegg_taxa_color_map(colored[nrow(colored):1L, ])
  stopifnot(is.character(palette), !is.null(names(palette)),
            all(unique(colored$taxon) %in% names(palette)),
            identical(palette[sort(names(palette))], reversed_palette[sort(names(reversed_palette))]))
  stopifnot(all(grepl("^#[[:xdigit:]]{6}$", unname(palette))))

  # Pagination preserves pathway order and gives each selected pathway one page.
  pathway_ids <- sprintf("%05d", seq_len(30L))
  pages_input <- do.call(rbind, lapply(seq_along(pathway_ids), function(i) {
    transform(table, pathway_id = pathway_ids[[i]],
              pathway_name = paste("Pathway", pathway_ids[[i]]))
  }))
  pages <- paginate_kegg_taxa(pages_input, c("large", "small", "empty"),
                              width = 12, height = 9)
  stopifnot(is.list(pages), length(pages) > 1L,
            identical(unlist(pages, use.names = FALSE), pathway_ids),
            all(lengths(pages) >= 1L))

  # SQMtools may omit misc$samples even though every TPM matrix names its samples.
  sqm_stub <- list(
    misc = list(samples = character()),
    orfs = list(tpm = matrix(1, nrow = 1L, dimnames = list("orf1", "Sample1"))),
    functions = list(KEGG = list(tpm = matrix(1, nrow = 1L,
                                            dimnames = list("K00001", "Sample1"))))
  )
  repaired <- ensure_sqm_sample_metadata(sqm_stub)
  stopifnot(identical(repaired$misc$samples, "Sample1"),
            length(sqm_stub$misc$samples) == 0L)
  sqm_stub$misc$samples <- "Wrong"
  stopifnot(inherits(try(ensure_sqm_sample_metadata(sqm_stub), silent = TRUE),
                   "try-error"))

  # CLI defaults and overrides are independent of FLOW and PIE options.
  base <- list(project_dir = "project", output_dir = "output", mode = "kegg_taxa")
  defaults <- build_config(base)
  stopifnot("kegg_taxa" %in% defaults$mode)
  assert_equal(defaults$top_n_kegg_taxa, 10L, "default KEGG taxa cap")
  assert_equal(defaults$min_kegg_taxon_percent, 1, "default KEGG threshold")
  overridden <- build_config(c(base, list(top_n_kegg_taxa = "3",
                                          min_kegg_taxon_percent = "2.5")))
  assert_equal(overridden$top_n_kegg_taxa, 3L, "KEGG taxa cap override")
  assert_equal(overridden$min_kegg_taxon_percent, 2.5, "KEGG threshold override")
  stopifnot(inherits(try(build_config(c(base, list(top_n_kegg_taxa = "0"))), silent = TRUE),
                   "try-error"))
  stopifnot(inherits(try(build_config(c(base, list(min_kegg_taxon_percent = "101"))),
                       silent = TRUE), "try-error"))


  # Rendering writes the same values represented by the plot and one PNG per page.
  output_dir <- tempfile("kegg_taxa_render_")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  context <- list(output_dir = output_dir, selection = "defined")
  dimension <- list(list(name = "12x9", width = 12, height = 9, dpi = 72))
  rendered <- render_kegg_taxa(pages_input, context, "genus",
                               c("large", "small", "empty"), dimension,
                               build_kegg_taxa_color_map(pages_input))
  png_files <- rendered$output_file[rendered$output_type == "plot_png"]
  tsv_file <- rendered$output_file[rendered$output_type == "data_tsv"]
  stopifnot(length(tsv_file) == 1L, length(png_files) == length(pages),
            all(file.exists(rendered$output_file)),
            all(file.info(png_files)$size > 0))
  saved <- utils::read.delim(tsv_file, check.names = FALSE, stringsAsFactors = FALSE,
                             colClasses = c(pathway_id = "character"))
  assert_equal(saved$tpm, pages_input$tpm, "rendered TSV TPM")
  stopifnot(identical(saved$pathway_id, pages_input$pathway_id),
            identical(saved$sample, pages_input$sample),
            identical(saved$taxon, pages_input$taxon))
  png_header <- readBin(png_files[[1L]], what = "raw", n = 8L)
  stopifnot(identical(as.integer(png_header), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))

  # Many samples are divided into repeated legend blocks on a single-pathway page.
  sample_groups <- kegg_taxa_sample_groups(paste0("S", seq_len(15L)), width = 12)
  stopifnot(length(sample_groups) > 1L,
            identical(unlist(sample_groups, use.names = FALSE), paste0("S", seq_len(15L))))

  cat("KEGG taxon barplot tests passed\n")
}

run_tests()
