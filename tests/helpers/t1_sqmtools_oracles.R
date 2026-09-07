t1_expect_true <- function(condition, label) {
  if (!isTRUE(condition)) {
    stop(label, call. = FALSE)
  }
}

t1_expect_identical <- function(actual, expected, label) {
  if (!identical(actual, expected)) {
    stop(
      label,
      "; expected=", paste(expected, collapse = ","),
      "; observed=", paste(actual, collapse = ","),
      call. = FALSE
    )
  }
}

t1_expect_equal <- function(actual, expected, label, tolerance = 1e-8) {
  if (!isTRUE(all.equal(
    as.numeric(actual),
    as.numeric(expected),
    tolerance = tolerance,
    check.attributes = FALSE
  ))) {
    stop(
      label,
      "; expected=", paste(format(expected, digits = 15), collapse = ","),
      "; observed=", paste(format(actual, digits = 15), collapse = ","),
      call. = FALSE
    )
  }
}

t1_make_pathview_node_data <- function(node_table) {
  required_columns <- c("node_id", "kegg_names", "label", "type")
  missing_columns <- setdiff(required_columns, colnames(node_table))
  if (length(missing_columns) > 0L) {
    stop(
      "Pathview node fixture is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  node_ids <- as.character(node_table$node_id)
  if (length(node_ids) == 0L || anyNA(node_ids) || anyDuplicated(node_ids)) {
    stop("Pathview node IDs must be non-empty and unique.", call. = FALSE)
  }

  named_values <- function(values) stats::setNames(values, node_ids)
  kegg_names <- strsplit(
    as.character(node_table$kegg_names),
    ";",
    fixed = TRUE
  )
  kegg_names <- lapply(kegg_names, trimws)
  names(kegg_names) <- node_ids

  list(
    kegg.names = kegg_names,
    type = named_values(as.character(node_table$type)),
    component = named_values(node_ids),
    size = named_values(rep.int(1L, length(node_ids))),
    labels = named_values(as.character(node_table$label)),
    shape = named_values(rep.int("rectangle", length(node_ids))),
    x = named_values(seq_along(node_ids)),
    y = named_values(rep.int(1, length(node_ids))),
    width = named_values(rep.int(46, length(node_ids))),
    height = named_values(rep.int(17, length(node_ids)))
  )
}

t1_map_pathview_tpm <- function(sqm_object, node_table, selected_samples) {
  kegg_tpm <- sqm_object$functions$KEGG$tpm
  if (is.null(kegg_tpm)) {
    stop("SQM pathview oracle requires functions$KEGG$tpm.", call. = FALSE)
  }

  selected_samples <- as.character(selected_samples)
  missing_samples <- setdiff(selected_samples, colnames(kegg_tpm))
  if (length(missing_samples) > 0L) {
    stop(
      "SQM pathview oracle is missing samples: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  mapped <- pathview::node.map(
    mol.data = as.matrix(kegg_tpm[, selected_samples, drop = FALSE]),
    node.data = t1_make_pathview_node_data(node_table),
    node.types = "ortholog",
    node.sum = "sum",
    entrez.gnodes = FALSE
  )
  values <- as.data.frame(mapped[, selected_samples, drop = FALSE])
  values[] <- lapply(values, as.numeric)
  values[is.na(values)] <- 0
  values
}

t1_read_k01563_fixture <- function(
    path = file.path("tests", "fixtures", "pathview_k01563_nodes.tsv")) {
  if (!file.exists(path)) {
    stop("Missing K01563 pathview fixture: ", path, call. = FALSE)
  }

  fixture <- utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    colClasses = "character",
    check.names = FALSE
  )
  fixture$pathway_id <- sprintf("%05d", as.integer(fixture$pathway_id))
  fixture
}

t1_assert_plot_data <- function(actual, expected, label, tolerance = 1e-8) {
  required_columns <- c("sample", "item", "abun")
  t1_expect_identical(
    colnames(actual),
    required_columns,
    paste0(label, " column contract changed")
  )
  t1_expect_identical(
    as.character(actual$sample),
    as.character(expected$sample),
    paste0(label, " sample order changed")
  )
  t1_expect_identical(
    as.character(actual$item),
    as.character(expected$item),
    paste0(label, " taxon order changed")
  )
  t1_expect_equal(
    actual$abun,
    expected$abun,
    paste0(label, " values changed"),
    tolerance = tolerance
  )
  invisible(actual)
}
