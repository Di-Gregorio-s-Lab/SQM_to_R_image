source(file.path("tests", "helpers", "t1_sqmtools_oracles.R"))

suppressPackageStartupMessages(library(SQMtools))
suppressPackageStartupMessages(library(pathview))

sqmtools_version <- as.character(utils::packageVersion("SQMtools"))

# exportPathway() starts from the complete KEGG TPM matrix, maps orthologs with
# node.sum="sum", and replaces values that did not map with zero.
pathview_sqm <- list(
  functions = list(
    KEGG = list(
      tpm = data.frame(
        S1 = c(10, 3, 2),
        S2 = c(5, 7, 1),
        row.names = c("K00001", "K00002", "K00003"),
        check.names = FALSE
      )
    )
  )
)
pathview_nodes <- data.frame(
  node_id = c("single", "multi", "missing", "repeated"),
  kegg_names = c("K00001", "K00002;K00003", "K99999", "K00001"),
  label = c("single", "multi", "missing", "repeated"),
  type = rep("ortholog", 4L),
  stringsAsFactors = FALSE
)
mapped <- t1_map_pathview_tpm(pathview_sqm, pathview_nodes, c("S2", "S1"))

t1_expect_equal(mapped["single", ], c(5, 10), "Single-KO node mapping")
t1_expect_equal(mapped["multi", ], c(8, 5), "Multi-KO node sum mapping")
t1_expect_equal(mapped["missing", ], c(0, 0), "Missing-KO node mapping")
t1_expect_equal(mapped["repeated", ], c(5, 10), "Repeated node mapping")

# plotTaxonomy() starts from SQM$taxa[[rank]][[count]]. These cases freeze the
# exact Top-N/Other ordering and the treatment of special categories used by
# the script (no_partial_classifications=FALSE, rescale=FALSE).
taxonomy_abund <- rbind(
  Alpha = c(S1 = 10, S2 = 0),
  Beta = c(S1 = 5, S2 = 20),
  Gamma = c(S1 = 1, S2 = 1),
  Unmapped = c(S1 = 2, S2 = 3),
  Unclassified = c(S1 = 4, S2 = 5)
)
taxonomy_percent <- 100 * t(t(taxonomy_abund) / colSums(taxonomy_abund))
taxonomy_sqm <- structure(
  list(
    misc = list(samples = c("S1", "S2")),
    taxa = list(
      phylum = list(
        abund = taxonomy_abund,
        percent = taxonomy_percent
      )
    )
  ),
  class = "SQM"
)

top_n_plot <- SQMtools::plotTaxonomy(
  SQM = taxonomy_sqm,
  rank = "phylum",
  count = "abund",
  N = 2L,
  others = TRUE,
  samples = c("S2", "S1"),
  ignore_unmapped = FALSE,
  ignore_unclassified = FALSE,
  no_partial_classifications = FALSE,
  rescale = FALSE
)
top_n_expected <- data.frame(
  sample = rep(c("S2", "S1"), 3L),
  item = rep(c("Other", "Beta", "Alpha"), each = 2L),
  abun = c(9, 7, 20, 5, 0, 10),
  stringsAsFactors = FALSE
)
t1_assert_plot_data(top_n_plot$data, top_n_expected, "plotTaxonomy Top-N")

excluded_plot <- suppressWarnings(SQMtools::plotTaxonomy(
  SQM = taxonomy_sqm,
  rank = "phylum",
  count = "percent",
  N = 5L,
  others = TRUE,
  samples = c("S2", "S1"),
  ignore_unmapped = TRUE,
  ignore_unclassified = TRUE,
  no_partial_classifications = FALSE,
  rescale = FALSE
))
excluded_expected <- data.frame(
  sample = rep(c("S2", "S1"), 4L),
  item = rep(c("Other", "Beta", "Alpha", "Gamma"), each = 2L),
  abun = c(
    0, 0,
    68.9655172413793, 22.7272727272727,
    0, 45.4545454545455,
    3.44827586206897, 4.54545454545455
  ),
  stringsAsFactors = FALSE
)
t1_assert_plot_data(
  excluded_plot$data,
  excluded_expected,
  "plotTaxonomy exclusions"
)

fixture <- t1_read_k01563_fixture()
t1_expect_identical(
  fixture$node_id[fixture$pathway_id == "00361"],
  c("75", "77", "131", "132"),
  "Pathway 00361 K01563 fixture"
)
t1_expect_identical(
  fixture$node_id[fixture$pathway_id == "00625"],
  c("58", "61", "103"),
  "Pathway 00625 K01563 fixture"
)
t1_expect_true(
  all(fixture$kegg_names == "K01563" & fixture$type == "ortholog"),
  "K01563 fixture contains a non-ortholog or unexpected KO."
)

message(
  "PASS: SQMtools=", sqmtools_version,
  " pathview and plotTaxonomy oracle behavior characterized"
)
