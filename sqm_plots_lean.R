# Lean SqueezeMeta plotting pipeline.
# Source this file to reuse its pure helpers; main() runs only via Rscript.

`%||%` <- function(x, fallback) if (is.null(x)) fallback else x

colors_hex <- c(
  "#5d8aa8", "#e32636", "#efdecd", "#ffbf00", "#9966cc", "#a4c639",
  "#cd9575", "#915c83", "#008000", "#fbceb1", "#00ffff", "#4b5320",
  "#b2beb5", "#87a96b", "#ff9966", "#a52a2a", "#6e7f80", "#ff2052",
  "#007fff", "#f0ffff", "#89cff0", "#f4c2c2", "#21abcd", "#fae7b5",
  "#ffe135", "#848482", "#98777b", "#f5f5dc", "#3d2b1f", "#fe6f5e",
  "#000000", "#ffebcd", "#318ce7", "#ace5ee", "#faf0be", "#0000ff",
  "#a2a2d0", "#6699cc", "#0d98ba", "#8a2be2", "#de5d83", "#79443b",
  "#0095b6", "#e3dac9", "#cc0000", "#006a4e", "#873260", "#0070ff",
  "#b5a642", "#cb4154", "#1dacd6", "#66ff00", "#bf94e4", "#c32148",
  "#ff007f", "#08e8de", "#d19fe8", "#f4bbff", "#ff55a3", "#fb607f",
  "#004225", "#cd7f32", "#ffc1cc", "#e7feff", "#f0dc82"
)

default_workers <- function() {
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(cores) || cores < 1L) return(1L)
  as.integer(min(4L, cores))
}

safe_name <- function(x) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", trimws(as.character(x)))
  value <- gsub("^_+|_+$", "", value)
  ifelse(nzchar(value), value, "unnamed")
}

parse_args <- function(args) {
  flags <- c("plan_only", "refresh_kegg", "help")
  values <- list(plan_only = FALSE, refresh_kegg = FALSE, help = FALSE)
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (token == "-h") token <- "--help"
    if (!startsWith(token, "--")) {
      stop("Unexpected argument: ", token, call. = FALSE)
    }
    option <- substring(token, 3L)
    if (grepl("=", option, fixed = TRUE)) {
      pair <- strsplit(option, "=", fixed = TRUE)[[1L]]
      name <- pair[[1L]]
      value <- paste(pair[-1L], collapse = "=")
      if (!nzchar(value)) stop("Missing value for --", name, call. = FALSE)
    } else {
      name <- option
      if (name %in% flags) {
        values[[name]] <- TRUE
        i <- i + 1L
        next
      }
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
        stop("Missing value for --", name, call. = FALSE)
      }
      value <- args[[i + 1L]]
      i <- i + 1L
    }
    values[[name]] <- value
    i <- i + 1L
  }
  if (!is.null(values$workers)) {
    workers <- suppressWarnings(as.integer(values$workers))
    if (is.na(workers) || workers < 1L) stop("workers must be a positive integer.", call. = FALSE)
    values$workers <- workers
  } else {
    values$workers <- default_workers()
  }
  values
}

extract_ko_ids <- function(x) {
  if (length(x) != 1L || is.na(x)) return(character())
  hit <- regmatches(as.character(x), gregexpr("K[0-9]{5}", as.character(x), perl = TRUE))[[1L]]
  if (identical(hit, "")) character() else unique(hit)
}

extract_ec_codes <- function(x) {
  vapply(as.character(x), function(value) {
    hit <- regexec("\\[EC:([^]]+)\\]", value, perl = TRUE)
    parts <- regmatches(value, hit)[[1L]]
    if (length(parts) < 2L) return(NA_character_)
    paste(strsplit(trimws(parts[[2L]]), "[[:space:]]+")[[1L]], collapse = ";")
  }, character(1L), USE.NAMES = FALSE)
}

rank_top20 <- function(pathways, n = 20L, samples = character(),
                       filtered_taxon = NA_character_, filtered_taxon_rank = NA_character_) {
  required <- c("pathway_id", "pathway_name", "pathway_root", "pathway_category", "total_tpm")
  missing <- setdiff(required, names(pathways))
  if (length(missing)) stop("Top20 input missing: ", paste(missing, collapse = ", "), call. = FALSE)
  n <- suppressWarnings(as.integer(n))
  if (is.na(n) || n < 1L) stop("n must be a positive integer.", call. = FALSE)
  pathways <- pathways[order(-pathways$total_tpm, pathways$pathway_name,
                             pathways$pathway_root, pathways$pathway_category), , drop = FALSE]
  pathways <- head(pathways, n)
  data.frame(
    rank = seq_len(nrow(pathways)),
    pathway_id = as.character(pathways$pathway_id),
    pathway_name = as.character(pathways$pathway_name),
    pathway_root = as.character(pathways$pathway_root),
    pathway_category = as.character(pathways$pathway_category),
    total_tpm = as.numeric(pathways$total_tpm),
    samples = rep(paste(samples, collapse = ","), nrow(pathways)),
    filtered_taxon = rep(filtered_taxon, nrow(pathways)),
    filtered_taxon_rank = rep(filtered_taxon_rank, nrow(pathways)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

top20_path <- function(output_dir, taxon_rank = NULL, taxon = NULL) {
  if (is.null(taxon_rank) || is.null(taxon)) return(file.path(output_dir, "top20.tsv"))
  file.path(output_dir, "taxon_filter", safe_name(taxon_rank), safe_name(taxon), "top20.tsv")
}

write_table <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  path
}

write_top20 <- function(path, table) write_table(table, path)

validate_matrix <- function(x, samples, label) {
  frame <- as.data.frame(x, check.names = FALSE)
  missing <- setdiff(samples, names(frame))
  if (length(missing)) stop(label, " missing samples: ", paste(missing, collapse = ", "), call. = FALSE)
  values <- unlist(frame[samples], use.names = FALSE)
  if (!all(vapply(frame[samples], is.numeric, logical(1L))) || anyNA(values) ||
      any(!is.finite(values)) || any(values < 0)) {
    stop(label, " must contain finite non-negative numbers.", call. = FALSE)
  }
  invisible(frame)
}

allocate_ko_tpm <- function(orf_table, orf_tpm, official_ko_tpm,
                            pathway_ko_ids, samples, tolerance = 1e-8, context = "analysis") {
  orf_table <- as.data.frame(orf_table, check.names = FALSE)
  orf_tpm <- validate_matrix(orf_tpm, samples, "ORF TPM")
  official <- validate_matrix(official_ko_tpm, samples, "Official KO TPM")
  if (!"KEGG ID" %in% names(orf_table)) stop("ORF table missing KEGG ID.", call. = FALSE)
  ids <- list(metadata = rownames(orf_table), tpm = rownames(orf_tpm), official = rownames(official))
  if (any(vapply(ids, is.null, logical(1L))) || anyDuplicated(ids$metadata) || anyDuplicated(ids$tpm) ||
      anyDuplicated(ids$official) || !setequal(ids$metadata, ids$tpm)) {
    stop("ORF and KO matrices require unique matching row names.", call. = FALSE)
  }
  pathway_ko_ids <- sort(unique(as.character(pathway_ko_ids)))
  if (!length(pathway_ko_ids)) return(data.frame(
    orf_id = character(), sample = character(), ko_id = character(), tpm = numeric()
  ))
  all_kos <- lapply(as.character(orf_table[ids$metadata, "KEGG ID"]), extract_ko_ids)
  selected <- lapply(all_kos, intersect, pathway_ko_ids)
  source <- rep.int(seq_along(selected), lengths(selected))
  ko_id <- unlist(selected, use.names = FALSE)
  allocated <- do.call(rbind, lapply(samples, function(sample) {
    share <- as.numeric(orf_tpm[ids$metadata[source], sample]) / lengths(all_kos)[source]
    keep <- share > 0
    data.frame(orf_id = ids$metadata[source][keep], sample = rep(sample, sum(keep)), ko_id = ko_id[keep],
               raw_share = share[keep], stringsAsFactors = FALSE)
  }))
  targets <- expand.grid(sample = samples, ko_id = pathway_ko_ids,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  targets$target <- mapply(function(sample, ko_id) {
    if (ko_id %in% rownames(official)) as.numeric(official[ko_id, sample]) else 0
  }, targets$sample, targets$ko_id)
  if (nrow(allocated)) {
    observed <- aggregate(allocated$raw_share, allocated[c("sample", "ko_id")], sum)
    names(observed)[[3L]] <- "observed"
    targets <- merge(targets, observed, by = c("sample", "ko_id"), all.x = TRUE, sort = FALSE)
    targets$observed[is.na(targets$observed)] <- 0
  } else targets$observed <- 0
  difference <- abs(targets$target - targets$observed)
  bad <- which(difference > tolerance * pmax(1, targets$target, targets$observed))
  if (length(bad)) {
    row <- bad[[1L]]
    factor <- if (targets$observed[[row]] > 0) targets$target[[row]] / targets$observed[[row]] else NA_real_
    stop("KO TPM mismatch [", context, "] at ", targets$sample[[row]], "/", targets$ko_id[[row]],
         ": observed=", format(targets$observed[[row]], digits = 17),
         ", target=", format(targets$target[[row]], digits = 17),
         ", difference=", format(difference[[row]], digits = 17),
         ", factor=", format(factor, digits = 17), ".", call. = FALSE)
  }
  result <- allocated[c("orf_id", "sample", "ko_id", "raw_share")]
  names(result)[[4L]] <- "tpm"
  rownames(result) <- NULL
  result
}

run_tasks <- function(tasks, worker, workers = 1L) {
  execute <- function(task) tryCatch(
    list(id = as.character(task$id), value = worker(task), error = NULL),
    error = function(error) list(id = as.character(task$id), value = NULL,
                                 error = conditionMessage(error))
  )
  workers <- as.integer(workers)
  if (workers > 1L && length(tasks) > 1L) {
    cluster <- parallel::makeCluster(min(workers, length(tasks)))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    worker_environment <- environment(worker)
    parallel::clusterExport(
      cluster,
      ls(envir = worker_environment, all.names = TRUE),
      envir = worker_environment
    )
    completed <- parallel::parLapply(cluster, tasks, execute)
  } else {
    completed <- lapply(tasks, execute)
  }
  completed <- completed[order(vapply(completed, `[[`, character(1L), "id"))]
  failed <- vapply(completed, function(x) !is.null(x$error), logical(1L))
  results <- lapply(completed[!failed], `[[`, "value")
  names(results) <- vapply(completed[!failed], `[[`, character(1L), "id")
  errors <- data.frame(
    task_id = vapply(completed[failed], `[[`, character(1L), "id"),
    error = vapply(completed[failed], `[[`, character(1L), "error"),
    stringsAsFactors = FALSE
  )
  list(results = results, errors = errors)
}

top_values <- function(data, label, value, n, reserved = c("Unclassified", "Unmapped")) {
  keep <- !data[[label]] %in% reserved & data[[value]] > 0
  if (!any(keep)) return(character())
  totals <- aggregate(data[[value]][keep], list(data[[label]][keep]), sum)
  names(totals) <- c("label", "value")
  head(totals$label[order(-totals$value, totals$label)], n)
}

build_funz_table <- function(allocated, official_totals, samples, top_n = 20L) {
  official <- as.data.frame(official_totals, check.names = FALSE)
  required <- c("sample", "ko_id", "tpm")
  if (length(setdiff(required, names(official)))) stop("Official totals require sample, ko_id, tpm.", call. = FALSE)
  official <- official[official$sample %in% samples, , drop = FALSE]
  ko_totals <- aggregate(official$tpm, list(official$ko_id), sum)
  top_kos <- head(ko_totals[[1L]][order(-ko_totals[[2L]], ko_totals[[1L]])], as.integer(top_n))
  metadata <- unique(as.data.frame(allocated)[, intersect(c("ko_id", "kegg_function", "ec_codes"), names(allocated)), drop = FALSE])
  table <- merge(official, metadata, by = "ko_id", all.x = TRUE, sort = FALSE)
  table$ko_id <- ifelse(table$ko_id %in% top_kos, as.character(table$ko_id), "Other")
  if ("kegg_function" %in% names(table)) {
    table$kegg_function[is.na(table$kegg_function)] <- table$ko_id[is.na(table$kegg_function)]
    table$kegg_function[table$ko_id == "Other"] <- "Other KOs"
  }
  if ("ec_codes" %in% names(table)) {
    table$ec_codes[is.na(table$ec_codes)] <- ""
    table$ec_codes[table$ko_id == "Other"] <- ""
  }
  group <- c("sample", "ko_id", intersect(c("kegg_function", "ec_codes"), names(table)))
  table <- aggregate(table$tpm, table[group], sum)
  names(table)[ncol(table)] <- "tpm"
  if ("ec_codes" %in% names(table)) table$ec_codes[!nzchar(table$ec_codes)] <- NA_character_
  denominators <- aggregate(table$tpm, list(sample = table$sample), sum)
  names(denominators)[[2L]] <- "pathway_tpm"
  table <- merge(table, denominators, by = "sample", all.x = TRUE, sort = FALSE)
  table$percent <- ifelse(table$pathway_tpm > 0, 100 * table$tpm / table$pathway_tpm, NA_real_)
  table$sample <- factor(table$sample, levels = samples)
  table <- table[order(table$sample, table$ko_id == "Other", -table$tpm, table$ko_id), , drop = FALSE]
  table$sample <- as.character(table$sample)
  rownames(table) <- NULL
  table
}

build_flow_table <- function(allocated, samples, rank, top_n_ko = 20L, top_n_taxa = 15L) {
  data <- as.data.frame(allocated, check.names = FALSE)
  if (!rank %in% names(data)) stop("Missing taxonomy rank: ", rank, call. = FALSE)
  data <- data[data$sample %in% samples & data$tpm > 0, , drop = FALSE]
  if (!"kegg_function" %in% names(data)) data$kegg_function <- data$ko_id
  if (!"ec_codes" %in% names(data)) data$ec_codes <- NA_character_
  data$taxon <- as.character(data[[rank]])
  data$taxon[is.na(data$taxon) | !nzchar(trimws(data$taxon))] <- "Unclassified"
  if (any(data$taxon == "Other")) stop("Other is reserved for collapsed taxa.", call. = FALSE)
  kos <- top_values(data, "ko_id", "tpm", top_n_ko, reserved = character())
  taxa <- top_values(data, "taxon", "tpm", top_n_taxa)
  metadata <- data[match(kos, data$ko_id), c("ko_id", "kegg_function", "ec_codes"), drop = FALSE]
  metadata$kegg_function[is.na(metadata$kegg_function) | !nzchar(metadata$kegg_function)] <-
    metadata$ko_id[is.na(metadata$kegg_function) | !nzchar(metadata$kegg_function)]
  metadata <- rbind(metadata, data.frame(
    ko_id = "Other", kegg_function = "Other KOs", ec_codes = NA_character_
  ))
  data$ko_id <- ifelse(data$ko_id %in% kos, data$ko_id, "Other")
  data$taxon <- ifelse(data$taxon %in% c(taxa, "Unclassified", "Unmapped"), data$taxon, "Other")
  table <- aggregate(data$tpm, data[c("sample", "taxon", "ko_id")], sum)
  names(table)[[4L]] <- "tpm"
  denominators <- aggregate(table$tpm, list(sample = table$sample), sum)
  names(denominators)[[2L]] <- "pathway_tpm"
  table <- merge(table, denominators, by = "sample", sort = FALSE)
  taxon_totals <- aggregate(table$tpm, table[c("sample", "taxon")], sum)
  names(taxon_totals)[[3L]] <- "taxon_tpm"
  ko_totals <- aggregate(table$tpm, table[c("sample", "ko_id")], sum)
  names(ko_totals)[[3L]] <- "ko_tpm"
  table <- merge(table, taxon_totals, by = c("sample", "taxon"), sort = FALSE)
  table <- merge(table, ko_totals, by = c("sample", "ko_id"), sort = FALSE)
  table <- merge(table, metadata, by = "ko_id", all.x = TRUE, sort = FALSE)
  table$flow_percent <- 100 * table$tpm / table$pathway_tpm
  table$percent <- table$flow_percent
  table$taxon_percent <- 100 * table$taxon_tpm / table$pathway_tpm
  table$ko_percent <- 100 * table$ko_tpm / table$pathway_tpm
  table$taxon_tpm <- NULL
  table$ko_tpm <- NULL
  table$rank <- rank
  table$sample <- factor(table$sample, levels = samples)
  table <- table[order(table$sample, -table$tpm, table$taxon, table$ko_id), , drop = FALSE]
  table$sample <- as.character(table$sample)
  rownames(table) <- NULL
  table
}

build_pie_table <- function(allocated, official_totals, samples, ko_id, rank, top_n_taxa = 15L) {
  if (length(samples) != 1L) stop("PIE requires exactly one sample.", call. = FALSE)
  sample <- samples[[1L]]
  data <- as.data.frame(allocated, check.names = FALSE)
  data <- data[data$sample == sample & data$ko_id == ko_id & data$tpm > 0, , drop = FALSE]
  if (!nrow(data)) return(data.frame())
  data$taxon <- as.character(data[[rank]])
  data$taxon[is.na(data$taxon) | !nzchar(trimws(data$taxon))] <- "Unclassified"
  if (any(data$taxon == "Other")) stop("Other is reserved for collapsed taxa.", call. = FALSE)
  target <- official_totals$tpm[official_totals$sample == sample & official_totals$ko_id == ko_id]
  if (length(target) != 1L || target < 0) stop("Missing official KO total.", call. = FALSE)
  data$tpm <- data$tpm * target / sum(data$tpm)
  taxa <- top_values(data, "taxon", "tpm", top_n_taxa)
  data$taxon <- ifelse(data$taxon %in% c(taxa, "Unclassified", "Unmapped"), data$taxon, "Other")
  table <- aggregate(data$tpm, list(taxon = data$taxon), sum)
  names(table)[[2L]] <- "tpm"
  table <- table[order(-table$tpm, table$taxon), , drop = FALSE]
  table$pct <- table$tpm / sum(table$tpm)
  table$sample <- sample
  table$ko_id <- ko_id
  table$rank <- rank
  table$ko_sample_tpm <- target
  table$pathway_sample_tpm <- sum(official_totals$tpm[official_totals$sample == sample])
  table$ko_pathway_percent <- 100 * target / table$pathway_sample_tpm
  rownames(table) <- NULL
  table
}


build_kegg_taxa_pathway_table <- function(allocated, pathway_id, pathway_name, samples,
                                          rank, top_n = 10L, min_percent = 1) {
  data <- as.data.frame(allocated, check.names = FALSE)
  if (!rank %in% names(data)) {
    if (nrow(data)) stop("Missing taxonomy rank: ", rank, call. = FALSE)
    data[[rank]] <- character()
  }
  if (length(samples) < 1L || anyDuplicated(samples)) stop("Samples must be unique and nonempty.", call. = FALSE)
  if (length(top_n) != 1L || is.na(top_n) || top_n < 1L || top_n != as.integer(top_n))
    stop("top_n must be a positive integer.", call. = FALSE)
  if (length(min_percent) != 1L || is.na(min_percent) || !is.finite(min_percent) ||
      min_percent < 0 || min_percent > 100) stop("min_percent must be between 0 and 100.", call. = FALSE)
  data <- data[data$sample %in% samples & data$tpm > 0, c("sample", rank, "tpm"), drop = FALSE]
  names(data)[[2L]] <- "taxon"
  data$taxon <- trimws(as.character(data$taxon))
  data$taxon[is.na(data$taxon) | !nzchar(data$taxon)] <- "Unclassified"
  if (any(data$taxon == "Other")) stop("Other is reserved for collapsed taxa.", call. = FALSE)
  if (nrow(data)) {
    data <- aggregate(data$tpm, data[c("sample", "taxon")], sum)
    names(data)[[3L]] <- "tpm"
  } else data <- data.frame(sample = character(), taxon = character(), tpm = numeric())
  selected <- character()
  for (sample in samples) {
    current <- data[data$sample == sample, , drop = FALSE]
    total <- sum(current$tpm)
    candidates <- current[!current$taxon %in% c("Unclassified", "Unmapped") &
                            current$tpm >= total * min_percent / 100, , drop = FALSE]
    if (nrow(candidates)) {
      candidates <- candidates[order(-candidates$tpm, candidates$taxon), , drop = FALSE]
      selected <- union(selected, head(candidates$taxon, top_n))
    }
  }
  if (nrow(data)) {
    data$taxon[!data$taxon %in% c(selected, "Unclassified", "Unmapped")] <- "Other"
    data <- aggregate(data$tpm, data[c("sample", "taxon")], sum)
    names(data)[[3L]] <- "tpm"
  }
  taxa <- unique(c(sort(selected), intersect(c("Unclassified", "Unmapped", "Other"), data$taxon)))
  if (!length(taxa)) taxa <- "Other"
  complete <- expand.grid(sample = samples, taxon = taxa, stringsAsFactors = FALSE)
  data <- merge(complete, data, by = c("sample", "taxon"), all.x = TRUE, sort = FALSE)
  data$tpm[is.na(data$tpm)] <- 0
  totals <- aggregate(data$tpm, list(sample = data$sample), sum)
  names(totals)[[2L]] <- "pathway_tpm"
  data <- merge(data, totals, by = "sample", sort = FALSE)
  data$percent <- ifelse(data$pathway_tpm > 0, 100 * data$tpm / data$pathway_tpm, 0)
  data$pathway_id <- as.character(pathway_id)
  data$pathway_name <- as.character(pathway_name)
  data$rank <- rank
  data <- data[c("pathway_id", "pathway_name", "sample", "taxon", "tpm", "pathway_tpm", "percent", "rank")]
  data <- data[order(match(data$sample, samples), match(data$taxon, taxa)), , drop = FALSE]
  rownames(data) <- NULL
  data
}

selection_dir <- function(selection) if (identical(selection, "defined")) "definiti" else "top20"

build_flow_color_map <- function(taxon_levels, ko_levels = character()) {
  categories <- unique(c(as.character(taxon_levels), as.character(ko_levels)))
  categories <- categories[!is.na(categories) & nzchar(categories)]
  non_other <- setdiff(categories, "Other")
  base_colors <- unique(colors_hex)
  category_colors <- if (length(non_other) <= length(base_colors)) {
    base_colors[seq_along(non_other)]
  } else {
    grDevices::hcl.colors(length(non_other), palette = "Dynamic")
  }
  stats::setNames(
    c(category_colors, if ("Other" %in% categories) "grey70" else character()),
    c(non_other, if ("Other" %in% categories) "Other" else character())
  )
}

flow_category_levels <- function(table, column) {
  totals <- aggregate(table$tpm, list(label = as.character(table[[column]])), sum)
  totals$label[order(totals$label == "Other", totals$x, totals$label)]
}

format_flow_percent <- function(x) {
  labels <- paste0(formatC(x, format = "f", digits = 1), "%")
  labels[is.na(x)] <- "NA%"
  labels[!is.na(x) & x > 0 & x < 0.1] <- "<0.1%"
  labels
}

build_flow_legend_labels <- function(table, taxon_levels, ko_levels) {
  taxon_rows <- match(taxon_levels, as.character(table$taxon))
  ko_rows <- match(ko_levels, as.character(table$ko_id))
  ec_codes <- trimws(as.character(table$ec_codes[ko_rows]))
  ec_codes[is.na(ec_codes) | !nzchar(ec_codes)] <- "NA"
  ko_names <- ifelse(ko_levels == "Other", "Other KOs", paste0(ko_levels, " / EC ", ec_codes))
  list(
    taxonomy = stats::setNames(
      paste0(taxon_levels, " | ", format_flow_percent(table$taxon_percent[taxon_rows])),
      taxon_levels
    ),
    functional = stats::setNames(
      paste0(ko_names, " | ", format_flow_percent(table$ko_percent[ko_rows])),
      ko_levels
    )
  )
}

make_flow_plot <- function(table, pathway_name, rank, sample_name) {
  taxon_levels <- flow_category_levels(table, "taxon")
  ko_levels <- flow_category_levels(table, "ko_id")
  palette <- build_flow_color_map(taxon_levels, ko_levels)
  legend_labels <- build_flow_legend_labels(table, taxon_levels, ko_levels)
  plot_table <- transform(
    table,
    taxon = factor(taxon, levels = taxon_levels),
    ko_id = factor(ko_id, levels = ko_levels)
  )
  ko_legend <- data.frame(ko_id = factor(ko_levels, levels = ko_levels))
  ggplot2::ggplot(plot_table, ggplot2::aes(axis1 = taxon, axis2 = ko_id, y = flow_percent)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = taxon), alpha = 0.72, width = 1 / 12) +
    ggalluvial::geom_stratum(
      ggplot2::aes(fill = ggplot2::after_stat(stratum)),
      width = 1 / 5,
      color = "grey35",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      ggplot2::aes(label = ggplot2::after_stat(stratum)),
      size = 2.8
    ) +
    ggplot2::scale_fill_manual(
      name = "Taxonomy", values = palette, breaks = taxon_levels,
      labels = unname(legend_labels$taxonomy[taxon_levels]), drop = FALSE
    ) +
    ggplot2::geom_point(
      data = ko_legend,
      ggplot2::aes(x = 1, y = 0, colour = ko_id),
      inherit.aes = FALSE,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_colour_manual(
      name = "Function (KO / EC)", values = palette, breaks = ko_levels,
      labels = unname(legend_labels$functional[ko_levels]), drop = FALSE
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1),
      colour = ggplot2::guide_legend(
        order = 2, override.aes = list(alpha = 1, shape = 15, size = 4)
      )
    ) +
    ggplot2::scale_x_discrete(limits = c("Taxon", "KO"), expand = c(0.08, 0.08)) +
    ggplot2::labs(
      title = paste0("Flowplot - ", pathway_name, " - ", rank, " - ", sample_name),
      subtitle = "Taxonomy-to-KO flow from the ORF x sample x KO table",
      x = NULL,
      y = "Relative flow (%)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(size = 8),
      legend.key.size = grid::unit(0.55, "lines"),
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

make_flow_sankey <- function(table, pathway_name, rank, sample_name) {
  taxon_levels <- flow_category_levels(table, "taxon")
  ko_levels <- flow_category_levels(table, "ko_id")
  palette <- build_flow_color_map(taxon_levels, ko_levels)
  legend_labels <- build_flow_legend_labels(table, taxon_levels, ko_levels)
  taxon_percents <- stats::setNames(
    table$taxon_percent[match(taxon_levels, as.character(table$taxon))], taxon_levels
  )
  ko_percents <- stats::setNames(
    table$ko_percent[match(ko_levels, as.character(table$ko_id))], ko_levels
  )
  plotly::plot_ly(
    type = "sankey",
    arrangement = "snap",
    valueformat = ".2f",
    valuesuffix = "%",
    node = list(
      pad = 18,
      thickness = 18,
      line = list(color = "rgba(70,70,70,0.35)", width = 0.5),
      label = c(
        unname(legend_labels$taxonomy[taxon_levels]),
        unname(legend_labels$functional[ko_levels])
      ),
      color = unname(palette[c(taxon_levels, ko_levels)])
    ),
    link = list(
      source = match(as.character(table$taxon), taxon_levels) - 1L,
      target = length(taxon_levels) + match(as.character(table$ko_id), ko_levels) - 1L,
      value = table$flow_percent,
      color = unname(grDevices::adjustcolor(
        palette[as.character(table$taxon)], alpha.f = 0.65
      )),
      customdata = paste0(
        "Taxon: ", table$taxon,
        "<br>Taxon share: ",
        format_flow_percent(unname(taxon_percents[as.character(table$taxon)])),
        "<br>KO: ", table$ko_id,
        "<br>EC: ", ifelse(is.na(table$ec_codes), "NA", table$ec_codes),
        "<br>Function: ", table$kegg_function,
        "<br>Function share: ",
        format_flow_percent(unname(ko_percents[as.character(table$ko_id)])),
        "<br>TPM: ", sprintf("%.3f", table$tpm),
        "<br>Flow: ", sprintf("%.2f", table$flow_percent), "%"
      ),
      hovertemplate = "%{customdata}<extra></extra>"
    )
  ) |>
    plotly::layout(
      title = list(text = paste0(
        "Flowplot - ", pathway_name, " - ", rank, " - ", sample_name,
        "<br><sup>Taxonomy-to-KO flow based on TPM</sup>"
      )),
      font = list(size = 11),
      margin = list(l = 20, r = 20, t = 60, b = 20)
    )
}

save_plot <- function(plot, path, width, height, dpi) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi, units = "in")
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Plot was not written: ", path, call. = FALSE)
  path
}

artifact_manifest <- function(mode, type, path, pathway_id, pathway_name, context,
                              sample = NA_character_, rank = NA_character_) {
  data.frame(
    mode, output_type = type, output_file = normalizePath(path, winslash = "/", mustWork = TRUE),
    pathway_id, pathway = pathway_name, samples = sample, rank,
    pathway_selection = context$selection %||% "defined",
    filtered_taxon = context$filtered_taxon %||% NA_character_,
    filtered_taxon_rank = context$filtered_taxon_rank %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

make_funz_plot <- function(table, pathway_name) {
  totals <- aggregate(table$tpm, list(ko_id = as.character(table$ko_id)), sum)
  ko_levels <- totals$ko_id[order(totals$ko_id == "Other", -totals$x, totals$ko_id)]
  samples <- unique(as.character(table$sample))
  labels <- vapply(ko_levels, function(ko) {
    rows <- table[as.character(table$ko_id) == ko, , drop = FALSE]
    rows <- rows[match(samples, as.character(rows$sample)), , drop = FALSE]
    ecs <- rows$ec_codes[!is.na(rows$ec_codes) & nzchar(rows$ec_codes)]
    ec <- if (ko == "Other" || !length(ecs)) "NA" else ecs[[1L]]
    name <- if (ko == "Other") "Other KOs" else ko
    paste0(name, " / EC ", ec, " | ", paste(format_flow_percent(rows$percent), collapse = " | "))
  }, character(1L))
  palette <- c(
    stats::setNames(rep(colors_hex, length.out = length(setdiff(ko_levels, "Other"))),
                    setdiff(ko_levels, "Other")),
    if ("Other" %in% ko_levels) c(Other = "grey70") else NULL
  )
  plot_table <- transform(table, ko_id = factor(ko_id, levels = ko_levels))
  ggplot2::ggplot(plot_table, ggplot2::aes(x = sample, y = tpm, fill = ko_id)) +
    ggplot2::geom_col(color = "grey25", linewidth = 0.15, width = 0.78) +
    ggplot2::scale_fill_manual(values = palette, labels = stats::setNames(labels, ko_levels), drop = FALSE) +
    ggplot2::scale_x_discrete(limits = samples, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    ggplot2::labs(title = paste0("KO barplot - ", pathway_name), x = "Sample", y = "TPM",
                  fill = "KO / EC | % per sample") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"), legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"), legend.text = ggplot2::element_text(size = 8),
      legend.key.size = grid::unit(0.55, "lines")
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, byrow = TRUE))
}


build_kegg_taxa_color_map <- function(table) {
  taxa <- sort(setdiff(unique(as.character(table$taxon)), c("Other", "Unclassified", "Unmapped")))
  base <- unique(colors_hex)
  extra <- if (length(taxa) > length(base)) {
    grDevices::hcl.colors(length(taxa) - length(base), palette = "Dynamic")
  } else character()
  palette <- stats::setNames(c(base, extra)[seq_along(taxa)], taxa)
  c(palette,
    if ("Unmapped" %in% table$taxon) c(Unmapped = "#8a8a8a"),
    if ("Unclassified" %in% table$taxon) c(Unclassified = "#555555"),
    if ("Other" %in% table$taxon) c(Other = "#b3b3b3"))
}

kegg_taxa_sample_groups <- function(samples, width) {
  columns <- max(1L, floor((width - 3.5) / 1.1))
  split(samples, ceiling(seq_along(samples) / columns))
}

kegg_taxa_legend_rows <- function(table, pathway_id, samples, width) {
  current <- table[table$pathway_id == pathway_id, , drop = FALSE]
  groups <- length(kegg_taxa_sample_groups(samples, width))
  1L + groups * (length(unique(current$taxon)) + 1L)
}

paginate_kegg_taxa <- function(table, samples, width, height) {
  ids <- unique(as.character(table$pathway_id))
  if (!length(ids)) return(list())
  if (!is.finite(width) || width <= 0 || !is.finite(height) || height <= 0)
    stop("Image dimensions must be positive.", call. = FALSE)
  max_rows <- max(1L, floor((height - 3.1) / 0.19))
  max_pathways <- max(1L, floor(width / (1.3 + 0.65 * length(samples))))
  pages <- list()
  current <- character()
  used <- 0L
  for (id in ids) {
    needed <- kegg_taxa_legend_rows(table, id, samples, width)
    if (length(current) && (used + needed > max_rows || length(current) >= max_pathways)) {
      pages[[length(pages) + 1L]] <- current
      current <- character()
      used <- 0L
    }
    current <- c(current, id)
    used <- used + needed
  }
  if (length(current)) pages[[length(pages) + 1L]] <- current
  pages
}

make_kegg_taxa_plot <- function(table, samples, palette) {
  ids <- unique(as.character(table$pathway_id))
  plot_table <- table
  plot_table$sample <- factor(plot_table$sample, levels = samples)
  plot_table$pathway_id <- factor(plot_table$pathway_id, levels = ids)
  plot_table$taxon <- factor(plot_table$taxon, levels = names(palette))
  ggplot2::ggplot(plot_table, ggplot2::aes(x = sample, y = tpm, fill = taxon)) +
    ggplot2::geom_col(colour = "grey25", linewidth = 0.15, width = 0.78) +
    ggplot2::facet_grid(. ~ pathway_id, scales = "free_x", space = "free_x", drop = FALSE) +
    ggplot2::scale_x_discrete(limits = samples, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::labs(title = paste0("KEGG pathway taxonomy - ", unique(table$rank)[[1L]]),
                  x = "Sample within pathway", y = "TPM") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "none")
}

draw_kegg_taxa_legend <- function(table, ids, samples, palette, width, height) {
  row_height <- 0.19
  y <- height - 0.15
  label_width <- min(width * 0.42, width - 1.2)
  values_width <- width - label_width - 0.2
  sample_groups <- kegg_taxa_sample_groups(samples, width)
  fmt <- function(value) formatC(value, digits = 3, format = "fg", flag = "#")
  for (id in ids) {
    current <- table[table$pathway_id == id, , drop = FALSE]
    pathway_name <- unique(current$pathway_name)[[1L]]
    grid::grid.text(paste0(id, "  ", pathway_name), x = grid::unit(0.12, "in"),
                    y = grid::unit(y, "in"), just = c("left", "top"),
                    gp = grid::gpar(fontsize = 9, fontface = "bold"))
    y <- y - row_height
    taxa <- unique(as.character(current$taxon))
    taxa <- taxa[order(match(taxa, names(palette)))]
    for (group in sample_groups) {
      grid::grid.text("Taxon", x = grid::unit(0.37, "in"), y = grid::unit(y, "in"),
                      just = c("left", "top"), gp = grid::gpar(fontsize = 8, fontface = "bold"))
      for (i in seq_along(group)) {
        x <- label_width + values_width * (i - 0.5) / length(group)
        grid::grid.text(group[[i]], x = grid::unit(x, "in"), y = grid::unit(y, "in"),
                        just = c("centre", "top"), gp = grid::gpar(fontsize = 8, fontface = "bold"))
      }
      y <- y - row_height
      for (taxon in taxa) {
        grid::grid.rect(x = grid::unit(0.18, "in"), y = grid::unit(y - 0.055, "in"),
                        width = grid::unit(0.12, "in"), height = grid::unit(0.12, "in"),
                        just = c("centre", "centre"),
                        gp = grid::gpar(fill = palette[[taxon]], col = "grey30"))
        grid::grid.text(taxon, x = grid::unit(0.37, "in"), y = grid::unit(y, "in"),
                        just = c("left", "top"), gp = grid::gpar(fontsize = 8))
        for (i in seq_along(group)) {
          value <- current$tpm[current$sample == group[[i]] & current$taxon == taxon]
          x <- label_width + values_width * (i - 0.5) / length(group)
          grid::grid.text(fmt(value), x = grid::unit(x, "in"), y = grid::unit(y, "in"),
                          just = c("centre", "top"), gp = grid::gpar(fontsize = 8))
        }
        y <- y - row_height
      }
    }
  }
}

save_kegg_taxa_page <- function(table, ids, samples, palette, path, width, height, dpi) {
  page <- table[table$pathway_id %in% ids, , drop = FALSE]
  rows <- sum(vapply(ids, function(id) kegg_taxa_legend_rows(page, id, samples, width), integer(1L)))
  legend_height <- 0.3 + rows * 0.19
  actual_height <- max(height, legend_height + 3.1)
  plot_height <- actual_height - legend_height
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = width, height = actual_height, units = "in", res = dpi)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    2L, 1L, heights = grid::unit(c(plot_height, legend_height), "in")
  )))
  grid::pushViewport(grid::viewport(layout.pos.row = 1L))
  grid::grid.draw(ggplot2::ggplotGrob(make_kegg_taxa_plot(page, samples, palette)))
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2L))
  draw_kegg_taxa_legend(page, ids, samples, palette, width, legend_height)
  grid::popViewport(2L)
  path
}

render_kegg_taxa <- function(table, context, rank, samples, dimensions, palette) {
  directory <- file.path(context$output_dir, "kegg_taxa", selection_dir(context$selection), safe_name(rank))
  data_path <- write_table(table, file.path(directory, "kegg_taxa_data.tsv"))
  rows <- list(artifact_manifest("kegg_taxa", "data_tsv", data_path, NA_character_,
                                 "combined", context, paste(samples, collapse = ","), rank))
  for (dimension in dimensions) {
    pages <- paginate_kegg_taxa(table, samples, dimension$width, dimension$height)
    for (index in seq_along(pages)) {
      path <- file.path(directory, paste0("kegg_taxa_", safe_name(dimension$name),
                                           "_page_", sprintf("%02d", index), ".png"))
      save_kegg_taxa_page(table, pages[[index]], samples, palette, path,
                          dimension$width, dimension$height, dimension$dpi)
      rows[[length(rows) + 1L]] <- artifact_manifest(
        "kegg_taxa", "plot_png", path, NA_character_,
        paste(pages[[index]], collapse = ","), context, paste(samples, collapse = ","), rank
      )
    }
  }
  do.call(rbind, rows)
}

render_funz <- function(table, context, pathway_id, pathway_name, width, height, dpi,
                        dimension_name = paste0(width, "x", height)) {
  directory <- file.path(context$output_dir, "funz", "pathway", selection_dir(context$selection), safe_name(pathway_name))
  data_path <- write_table(table, file.path(directory, "barplot_ko_data.tsv"))
  plot <- make_funz_plot(table, pathway_name)
  plot_path <- save_plot(plot, file.path(directory, paste0("barplot_ko_", safe_name(dimension_name), ".png")), width, height, dpi)
  rbind(
    artifact_manifest("funz", "data_tsv", data_path, pathway_id, pathway_name, context),
    artifact_manifest("funz", "plot_png", plot_path, pathway_id, pathway_name, context)
  )
}

render_flow <- function(table, context, pathway_id, pathway_name, width, height, dpi,
                        dimension_name = paste0(width, "x", height)) {
  rank <- unique(table$rank)[[1L]]
  sample <- paste(unique(as.character(table$sample)), collapse = "_")
  directory <- file.path(context$output_dir, "flowplot", selection_dir(context$selection), safe_name(pathway_name), safe_name(rank))
  stem <- paste0("flow_", safe_name(sample))
  data_path <- write_table(table, file.path(directory, paste0(stem, "_data.tsv")))
  plot <- make_flow_plot(table, pathway_name, rank, sample)
  plot_path <- save_plot(plot, file.path(directory, paste0(stem, "_", safe_name(dimension_name), ".png")), width, height, dpi)
  rbind(
    artifact_manifest("flow", "data_tsv", data_path, pathway_id, pathway_name, context, rank = rank),
    artifact_manifest("flow", "plot_png", plot_path, pathway_id, pathway_name, context, rank = rank)
  )
}

render_flow_html <- function(table, context, pathway_id, pathway_name) {
  rank <- unique(table$rank)[[1L]]
  sample <- paste(unique(as.character(table$sample)), collapse = "_")
  directory <- file.path(context$output_dir, "flowplot", selection_dir(context$selection), safe_name(pathway_name), safe_name(rank))
  widget <- make_flow_sankey(table, pathway_name, rank, sample)
  path <- file.path(directory, paste0("flow_", safe_name(sample), ".html"))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  htmlwidgets::saveWidget(widget, path, selfcontained = TRUE)
  if (!file.exists(path) || file.info(path)$size <= 0) stop("FLOW HTML was not written.", call. = FALSE)
  artifact_manifest("flow", "plot_html", path, pathway_id, pathway_name, context, sample, rank)
}

render_pie <- function(table, context, pathway_id, pathway_name, width, height, dpi,
                       dimension_name = paste0(width, "x", height)) {
  sample <- unique(table$sample)[[1L]]
  ko_id <- unique(table$ko_id)[[1L]]
  rank <- unique(table$rank)[[1L]]
  directory <- file.path(context$output_dir, "pie", selection_dir(context$selection), safe_name(pathway_name),
                         safe_name(sample), safe_name(ko_id), safe_name(rank))
  data_path <- write_table(table, file.path(directory, "pie_data.tsv"))
  plot <- ggplot2::ggplot(table, ggplot2::aes(x = "", y = tpm, fill = taxon)) +
    ggplot2::geom_col(width = 1, colour = "white") + ggplot2::coord_polar(theta = "y") +
    ggplot2::theme_void()
  plot_path <- save_plot(plot, file.path(directory, paste0("pie_", safe_name(dimension_name), ".png")), width, height, dpi)
  rbind(
    artifact_manifest("pie", "data_tsv", data_path, pathway_id, pathway_name, context, sample, rank),
    artifact_manifest("pie", "plot_png", plot_path, pathway_id, pathway_name, context, sample, rank)
  )
}

render_taxonomy <- function(sqm_object, samples, rank, context, pathway_id, pathway_name,
                            width, height, dpi, plot_fun = SQMtools::plotTaxonomy, count = "abund",
                            dimension_name = paste0(width, "x", height)) {
  plot <- plot_fun(SQM = sqm_object, rank = rank, count = count, N = 15L, samples = samples,
                   ignore_unmapped = TRUE, ignore_unclassified = TRUE,
                   no_partial_classifications = FALSE, rescale = FALSE)
  if (is.null(plot$data)) stop("plotTaxonomy did not expose plot data.", call. = FALSE)
  root <- if (is.na(pathway_id)) "taxonomy_global" else "taxonomy_by_pathway"
  directory <- if (root == "taxonomy_global") {
    file.path(context$output_dir, root, count, safe_name(rank))
  } else {
    file.path(context$output_dir, root, selection_dir(context$selection), safe_name(pathway_name), count, safe_name(rank))
  }
  unlink(file.path(directory, "taxonomy_data.tsv"), force = TRUE)
  plot_path <- save_plot(plot, file.path(directory, paste0("taxonomy_", safe_name(dimension_name), ".png")), width, height, dpi)
  artifact_manifest("taxon", "plot_png", plot_path, pathway_id, pathway_name, context,
                    paste(samples, collapse = ","), rank)
}

export_pathview_isolated <- function(export_pathway_fn, sqm_object, pathway_id,
                                     selected_samples, final_dir) {
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(final_dir, c(
    "pathview_input_all_ko_complete_matrix.tsv", "pathview_render_config.tsv"
  )), force = TRUE)
  export_pathway_fn(
    SQM = sqm_object, pathway_id = pathway_id, count = "tpm", samples = selected_samples,
    split_samples = FALSE, log_scale = FALSE, output_dir = normalizePath(final_dir, winslash = "/"),
    output_suffix = paste0("pathview_", pathway_id)
  )
  files <- list.files(final_dir, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE & file.info(files)$size > 0]
  if (!length(files)) stop("Pathview produced no files.", call. = FALSE)
  data.frame(output_file = normalizePath(files, winslash = "/", mustWork = TRUE), stringsAsFactors = FALSE)
}

write_atomic <- function(contents, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("write_", tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  if (is.data.frame(contents)) {
    utils::write.table(contents, temporary, sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "NA")
  } else {
    writeLines(as.character(contents), temporary, useBytes = TRUE)
  }
  if (!file.copy(temporary, path, overwrite = TRUE)) stop("Cannot write ", path, call. = FALSE)
  path
}

validate_catalog <- function(catalog) {
  catalog <- as.data.frame(catalog, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("pathway_id", "pathway_name")
  if (length(setdiff(required, names(catalog))) || !nrow(catalog)) stop("Invalid KEGG catalog.", call. = FALSE)
  catalog$pathway_id <- as.character(catalog$pathway_id)
  if (anyNA(catalog$pathway_id) || any(!grepl("^[0-9]{5}$", catalog$pathway_id)) ||
      anyNA(catalog$pathway_name) || any(!nzchar(catalog$pathway_name))) {
    stop("Invalid KEGG catalog rows.", call. = FALSE)
  }
  catalog
}

get_kegg_catalog <- function(output_dir, refresh = FALSE, downloader) {
  path <- file.path(output_dir, "_cache", "kegg", "pathway_catalog.tsv")
  read_cache <- function() validate_catalog(utils::read.delim(
    path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = c(pathway_id = "character")
  ))
  if (!refresh && file.exists(path) && file.info(path)$size > 0) {
    cached <- tryCatch(read_cache(), error = function(...) NULL)
    if (!is.null(cached)) return(cached)
  }
  catalog <- validate_catalog(downloader())
  write_atomic(catalog, path)
  read_cache()
}

get_pathway_kos <- function(output_dir, pathway_id, refresh = FALSE, downloader, parser) {
  if (length(pathway_id) != 1L || is.na(pathway_id) || !grepl("^[0-9]{5}$", pathway_id)) {
    stop("pathway_id must contain exactly five digits.", call. = FALSE)
  }
  path <- file.path(output_dir, "_cache", "kegg", paste0("ko", pathway_id, ".xml"))
  parse_cache <- function() {
    kos <- sort(unique(as.character(parser(path))))
    if (!length(kos) || any(!grepl("^K[0-9]{5}$", kos))) stop("Invalid KGML KO membership.", call. = FALSE)
    kos
  }
  if (!refresh && file.exists(path) && file.info(path)$size > 0) {
    cached <- tryCatch(parse_cache(), error = function(...) NULL)
    if (!is.null(cached)) return(cached)
  }
  xml <- downloader(pathway_id)
  if (!length(xml) || !any(nzchar(xml))) stop("Empty KGML download.", call. = FALSE)
  write_atomic(xml, path)
  parse_cache()
}

kegg_roots <- c(
  "Metabolism", "Genetic Information Processing", "Environmental Information Processing",
  "Cellular Processes", "Organismal Systems", "Human Diseases"
)

split_pathways <- function(values) {
  entries <- strsplit(as.character(values), "\\s*\\|\\s*", perl = TRUE)
  source_index <- rep.int(seq_along(entries), lengths(entries))
  parts <- strsplit(unlist(entries, use.names = FALSE), "\\s*;\\s*", perl = TRUE)
  valid <- lengths(parts) == 3L
  if (!any(valid)) return(data.frame(
    source_index = integer(), pathway_root = character(), pathway_category = character(),
    pathway_name = character(), stringsAsFactors = FALSE
  ))
  fields <- do.call(rbind, parts[valid])
  valid_rows <- fields[, 1L] %in% kegg_roots & !grepl("not included", fields[, 2L], ignore.case = TRUE)
  unique(data.frame(
    source_index = source_index[valid][valid_rows], pathway_root = fields[valid_rows, 1L],
    pathway_category = fields[valid_rows, 2L],
    pathway_name = sub("\\s*-\\s*Reference pathway\\s*$", "", fields[valid_rows, 3L], ignore.case = TRUE),
    stringsAsFactors = FALSE
  ))
}

catalog_id <- function(names, catalog) {
  normalized <- function(x) tolower(trimws(sub("\\s*-\\s*Reference pathway\\s*$", "", x, ignore.case = TRUE)))
  catalog_names <- normalized(catalog$pathway_name)
  vapply(names, function(name) {
    key <- normalized(name)
    if (key == "carbon fixation pathways in prokaryotes") key <- "other carbon fixation pathways"
    hit <- unique(as.character(catalog$pathway_id[catalog_names == key]))
    if (length(hit) > 1L) stop("Ambiguous KEGG pathway name: ", name, call. = FALSE)
    if (length(hit)) hit else NA_character_
  }, character(1L), USE.NAMES = FALSE)
}

resolve_defined_pathways <- function(requested, catalog) {
  normalized <- function(x) tolower(trimws(sub("\\s*-\\s*Reference pathway\\s*$", "", x, ignore.case = TRUE)))
  rows <- lapply(requested, function(value) {
    hit <- if (grepl("^[0-9]{5}$", value)) {
      which(as.character(catalog$pathway_id) == value)
    } else {
      which(normalized(catalog$pathway_name) == normalized(value))
    }
    if (length(hit) != 1L) stop("Pathway must resolve exactly once: ", value, call. = FALSE)
    data.frame(
      pathway_id = as.character(catalog$pathway_id[[hit]]),
      pathway_name = as.character(catalog$pathway_name[[hit]]),
      stringsAsFactors = FALSE
    )
  })
  unique(do.call(rbind, rows))
}

compute_top20 <- function(sqm, samples, catalog, n = 20L,
                          filtered_taxon = NA_character_, filtered_taxon_rank = NA_character_) {
  paths <- sqm$misc$KEGG_paths %||% character()
  membership <- split_pathways(paths)
  if (!nrow(membership)) {
    empty <- data.frame(pathway_id = character(), pathway_name = character(), pathway_root = character(),
                        pathway_category = character(), total_tpm = numeric())
    return(rank_top20(empty, n, samples, filtered_taxon, filtered_taxon_rank))
  }
  official <- validate_matrix(sqm$functions$KEGG$tpm, samples, "Official KO TPM")
  rows <- match(names(paths)[membership$source_index], rownames(official))
  membership$total_tpm <- 0
  membership$total_tpm[!is.na(rows)] <- rowSums(official[rows[!is.na(rows)], samples, drop = FALSE])
  hierarchy <- unique(membership[c("pathway_name", "pathway_root", "pathway_category")])
  duplicates <- duplicated(hierarchy$pathway_name) | duplicated(hierarchy$pathway_name, fromLast = TRUE)
  if (any(duplicates)) stop("A pathway name belongs to multiple KEGG hierarchies.", call. = FALSE)
  totals <- aggregate(membership$total_tpm, membership[c("pathway_name", "pathway_root", "pathway_category")], sum)
  names(totals)[[4L]] <- "total_tpm"
  totals$pathway_id <- catalog_id(totals$pathway_name, catalog)
  rank_top20(totals[c("pathway_id", "pathway_name", "pathway_root", "pathway_category", "total_tpm")],
             n, samples, filtered_taxon, filtered_taxon_rank)
}

ko_metadata <- function(sqm, ko_ids) {
  labels <- sqm$misc$KEGG_names %||% character()
  values <- unname(as.character(labels[ko_ids]))
  values[is.na(values) | !nzchar(values)] <- ko_ids[is.na(values) | !nzchar(values)]
  data.frame(
    ko_id = ko_ids,
    kegg_function = values,
    ec_codes = extract_ec_codes(values),
    stringsAsFactors = FALSE
  )
}

subset_sqm_ids <- function(sqm, orf_ids) {
  orf_ids <- unique(as.character(orf_ids))
  if (!length(orf_ids)) return(NULL)
  if (inherits(sqm, "SQM")) {
    return(SQMtools::subsetORFs(
      SQM = sqm, orfs = orf_ids, tax_source = "orfs", trusted_functions_only = FALSE,
      ignore_unclassified_functions = FALSE, rescale_tpm = FALSE, rescale_copy_number = FALSE,
      recalculate_bin_stats = FALSE, contigs_override = NULL, allow_empty = FALSE
    ))
  }
  subset <- sqm
  for (name in names(subset$orfs)) {
    value <- subset$orfs[[name]]
    if (!is.null(rownames(value))) subset$orfs[[name]] <- value[orf_ids, , drop = FALSE]
  }
  subset
}

prepare_analysis <- function(sqm, pathway_ko_ids, samples, context = "analysis") {
  allocated <- allocate_ko_tpm(
    sqm$orfs$table, sqm$orfs$tpm, sqm$functions$KEGG$tpm,
    pathway_ko_ids, samples, context = context
  )
  if (nrow(allocated)) {
    tax <- as.data.frame(sqm$orfs$tax, check.names = FALSE)
    tax <- tax[match(allocated$orf_id, rownames(tax)), , drop = FALSE]
    for (rank in names(tax)) {
      values <- as.character(tax[[rank]])
      values[is.na(values) | !nzchar(trimws(values))] <- "Unclassified"
      allocated[[rank]] <- values
    }
    metadata <- ko_metadata(sqm, pathway_ko_ids)
    allocated <- merge(allocated, metadata, by = "ko_id", all.x = TRUE, sort = FALSE)
  } else {
    metadata <- ko_metadata(sqm, pathway_ko_ids)
  }
  official_matrix <- as.data.frame(sqm$functions$KEGG$tpm, check.names = FALSE)
  official <- expand.grid(ko_id = pathway_ko_ids, sample = samples, stringsAsFactors = FALSE)
  official$tpm <- mapply(function(ko, sample) {
    if (ko %in% rownames(official_matrix)) as.numeric(official_matrix[ko, sample]) else 0
  }, official$ko_id, official$sample)
  list(data = allocated, totals = official, metadata = metadata,
       pathway_sqm = sqm)
}

build_enzyme_table <- function(sqm, samples, requested_ecs = NULL) {
  matrix <- as.data.frame(sqm$functions$KEGG$tpm, check.names = FALSE)
  metadata <- ko_metadata(sqm, rownames(matrix))
  relations <- list()
  index <- 0L
  for (row in seq_len(nrow(metadata))) {
    ecs <- if (is.na(metadata$ec_codes[[row]])) character() else strsplit(metadata$ec_codes[[row]], ";", fixed = TRUE)[[1L]]
    for (ec in ecs) {
      if (!is.null(requested_ecs) && !ec %in% requested_ecs) next
      index <- index + 1L
      relations[[index]] <- data.frame(ko_id = metadata$ko_id[[row]], ec_code = ec, stringsAsFactors = FALSE)
    }
  }
  if (!length(relations)) return(data.frame())
  relation <- unique(do.call(rbind, relations))
  totals <- expand.grid(ko_id = rownames(matrix), sample = samples, stringsAsFactors = FALSE)
  totals$tpm <- mapply(function(ko, sample) as.numeric(matrix[ko, sample]), totals$ko_id, totals$sample)
  table <- merge(totals, relation, by = "ko_id")
  table <- aggregate(table$tpm, table[c("sample", "ec_code")], sum)
  names(table)[[3L]] <- "tpm"
  table[order(match(table$sample, samples), table$ec_code), , drop = FALSE]
}

enzyme_palette <- function(enzyme_ecs) {
  stats::setNames(rep(colors_hex, length.out = length(enzyme_ecs)), enzyme_ecs)
}

make_enzyme_barplot <- function(table, title) {
  ec_levels <- unique(as.character(table$ec_code))
  plot_table <- transform(table, ec_code = factor(ec_code, levels = ec_levels))
  ggplot2::ggplot(plot_table, ggplot2::aes(x = sample, y = tpm, fill = ec_code)) +
    ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single"),
                      color = "grey25", linewidth = 0.15, width = 0.78) +
    ggplot2::scale_fill_manual(values = enzyme_palette(ec_levels), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    ggplot2::labs(title = title, x = "Sample", y = "TPM", fill = "EC") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"), legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold")
    )
}

make_enzyme_lineplot <- function(table, title) {
  ec_levels <- unique(as.character(table$ec_code))
  plot_table <- transform(table, ec_code = factor(ec_code, levels = ec_levels))
  ggplot2::ggplot(plot_table, ggplot2::aes(x = sample, y = tpm, colour = ec_code, group = ec_code)) +
    ggplot2::geom_line(linewidth = 0.75) + ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = enzyme_palette(ec_levels), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    ggplot2::labs(title = title, x = "Sample", y = "TPM", colour = "EC") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"), legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold")
    )
}

render_enzymes <- function(table, context, dimensions, plot_types = c("bar", "line")) {
  if (!nrow(table)) return(data.frame())
  groups <- c(list(list(name = "enzimi", data = table, directory = file.path(context$output_dir, "funz", "enzimi", "insieme"))),
              lapply(unique(table$ec_code), function(ec) list(
                name = ec, data = table[table$ec_code == ec, , drop = FALSE],
                directory = file.path(context$output_dir, "funz", "enzimi", "separato", safe_name(ec))
              )))
  manifests <- list()
  for (group in groups) {
    data_path <- write_table(group$data, file.path(group$directory, "enzimi_data.tsv"))
    manifests[[length(manifests) + 1L]] <- artifact_manifest(
      "funz", "enzyme_data_tsv", data_path, NA_character_, group$name, context
    )
    for (dimension in dimensions) for (type in intersect(c("bar", "line"), plot_types)) {
      scope <- if (group$name == "enzimi") "combined" else paste("EC", group$name)
      title <- paste(if (type == "bar") "Enzyme barplot" else "Enzyme line chart", scope, sep = " - ")
      plot <- if (type == "bar") make_enzyme_barplot(group$data, title) else make_enzyme_lineplot(group$data, title)
      path <- save_plot(plot, file.path(group$directory, paste0("enzimi_", type, "_", safe_name(dimension$name), ".png")),
                        dimension$width, dimension$height, dimension$dpi)
      manifests[[length(manifests) + 1L]] <- artifact_manifest(
        "funz", "enzyme_plot_png", path, NA_character_, group$name, context
      )
    }
  }
  do.call(rbind, manifests)
}

deduplicate_manifest <- function(rows) {
  if (!length(rows)) return(data.frame())
  table <- do.call(rbind, rows)
  table[!duplicated(table$output_file), , drop = FALSE]
}

render_pathway_mode <- function(task) {
  analysis <- task$analysis
  context <- task$context
  entry <- task$entry
  config <- task$config
  manifests <- list()
  if (task$mode %in% c("funz", "flow", "pie") && !nrow(analysis$data)) return(data.frame())
  if (task$mode == "funz") {
    table <- build_funz_table(analysis$data, analysis$totals, config$samples, config$top_n_ko)
    if (nrow(table)) for (dimension in config$dimensions) {
      manifests[[length(manifests) + 1L]] <- render_funz(
        table, context, entry$pathway_id, entry$pathway_name,
        dimension$width, dimension$height, dimension$dpi, dimension$name
      )
    }
  }
  if (task$mode == "flow") for (rank in config$ranks) {
    table <- build_flow_table(analysis$data, config$samples, rank, config$top_n_ko, config$top_n_taxa)
    for (sample in config$samples) {
      current <- table[table$sample == sample, , drop = FALSE]
      if (!nrow(current)) next
      for (dimension in config$dimensions) manifests[[length(manifests) + 1L]] <- render_flow(
        current, context, entry$pathway_id, entry$pathway_name,
        dimension$width, dimension$height, dimension$dpi, dimension$name
      )
      if ("html" %in% config$formats) {
        manifests[[length(manifests) + 1L]] <- render_flow_html(
          current, context, entry$pathway_id, entry$pathway_name
        )
      }
    }
  }
  if (task$mode == "taxon" && !is.null(analysis$pathway_sqm)) for (count in config$counts) for (rank in config$ranks) {
    for (dimension in config$dimensions) render_taxonomy(
      analysis$pathway_sqm, config$samples, rank, context, entry$pathway_id, entry$pathway_name,
      dimension$width, dimension$height, dimension$dpi, config$plot_taxonomy_fn,
      count, dimension$name
    )
  }
  if (task$mode == "pie") for (sample in config$samples) for (ko_id in entry$ko_ids) for (rank in config$ranks) {
    table <- build_pie_table(analysis$data, analysis$totals, sample, ko_id, rank, config$top_n_taxa)
    if (!nrow(table)) next
    metadata <- analysis$metadata[analysis$metadata$ko_id == ko_id, , drop = FALSE]
    table$kegg_function <- metadata$kegg_function[[1L]]
    table$ec_codes <- metadata$ec_codes[[1L]]
    for (dimension in config$dimensions) manifests[[length(manifests) + 1L]] <- render_pie(
      table, context, entry$pathway_id, entry$pathway_name,
      dimension$width, dimension$height, dimension$dpi, dimension$name
    )
  }
  deduplicate_manifest(manifests)
}

write_mode_manifest <- function(rows, output_dir, mode) {
  if (!length(rows)) return(NULL)
  table <- deduplicate_manifest(rows)
  root <- c(funz = "funz", flow = "flowplot", pie = "pie", kegg_taxa = "kegg_taxa")[[mode]]
  path <- file.path(output_dir, root, paste0("manifest_", mode, ".tsv"))
  table$output_file <- vapply(table$output_file, function(file) {
    root_path <- normalizePath(file.path(output_dir, root), winslash = "/", mustWork = TRUE)
    substring(normalizePath(file, winslash = "/", mustWork = TRUE), nchar(root_path) + 2L)
  }, character(1L))
  write_table(table, path)
}

run_pipeline <- function(sqm, config, catalog, kgml_loader,
                         plot_taxonomy_fn = SQMtools::plotTaxonomy,
                         export_pathway_fn = SQMtools::exportPathway,
                         subset_fun = SQMtools::subsetFun) {
  required <- c("output_dir", "mode", "samples", "ranks", "counts", "dimensions",
                "top_n_pathways", "top_n_ko", "top_n_taxa", "workers", "plan_only")
  missing <- setdiff(required, names(config))
  if (length(missing)) stop("Pipeline config missing: ", paste(missing, collapse = ", "), call. = FALSE)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  context <- list(output_dir = config$output_dir, selection = NA_character_,
                  filtered_taxon = config$filtered_taxon %||% NA_character_,
                  filtered_taxon_rank = config$filtered_taxon_rank %||% NA_character_)
  top <- compute_top20(
    sqm, config$samples, catalog, config$top_n_pathways,
    context$filtered_taxon, context$filtered_taxon_rank
  )
  write_top20(top20_path(config$output_dir), top)
  empty_errors <- data.frame(task_id = character(), mode = character(), pathway_id = character(),
                             error = character(), stringsAsFactors = FALSE)
  if (isTRUE(config$plan_only)) {
    if (isTRUE(config$write_errors %||% TRUE)) write_table(empty_errors, file.path(config$output_dir, "errors.tsv"))
    return(invisible(list(errors = empty_errors)))
  }

  selection_modes <- if ("normal" %in% config$mode) "defined" else config$selection_modes %||% c("defined", "top20")
  entries <- list()
  if ("defined" %in% selection_modes && length(config$defined_pathways %||% character())) {
    defined <- resolve_defined_pathways(config$defined_pathways, catalog)
    for (row_index in seq_len(nrow(defined))) entries[[length(entries) + 1L]] <- list(
      pathway_id = defined$pathway_id[[row_index]], pathway_name = defined$pathway_name[[row_index]],
      selection = "defined"
    )
  }
  if ("top20" %in% selection_modes) for (row_index in seq_len(nrow(top))) {
    if (is.na(top$pathway_id[[row_index]])) next
    entries[[length(entries) + 1L]] <- list(pathway_id = as.character(top$pathway_id[[row_index]]),
                                            pathway_name = as.character(top$pathway_name[[row_index]]), selection = "top20")
  }

  errors <- list()
  prepared <- list()
  for (entry in entries) {
    key <- paste(entry$selection, entry$pathway_id, sep = ":")
    tryCatch({
      pathway_sqm <- subset_fun(
        SQM = sqm, fun = entry$pathway_name, columns = "KEGGPATH",
        ignore_case = FALSE, fixed = TRUE, allow_empty = FALSE
      )
      entry$ko_ids <- rownames(pathway_sqm$functions$KEGG$tpm)
      entry$analysis <- prepare_analysis(
        pathway_sqm, entry$ko_ids, config$samples,
        context = paste(key, entry$pathway_name)
      )
      prepared[[key]] <- entry
    }, error = function(error) {
      errors[[length(errors) + 1L]] <<- data.frame(
        task_id = paste0(key, ":prepare"), mode = "prepare", pathway_id = entry$pathway_id,
        error = conditionMessage(error), stringsAsFactors = FALSE
      )
    })
  }

  mode_list <- if ("huge" %in% config$mode) {
    c("funz", "flow", "taxon", "pie", "kegg_taxa")
  } else if ("normal" %in% config$mode) {
    c("funz", "flow", "taxon", "kegg_taxa")
  } else {
    intersect(config$mode, c("funz", "flow", "taxon", "pie", "kegg_taxa"))
  }
  tasks <- list()
  for (key in names(prepared)) for (mode in setdiff(mode_list, "kegg_taxa")) {
    entry <- prepared[[key]]
    if (mode == "pie" && !entry$selection %in% (config$pie_selection_modes %||% selection_modes)) next
    task_context <- context
    task_context$selection <- entry$selection
    task_config <- config
    task_config$plot_taxonomy_fn <- plot_taxonomy_fn
    tasks[[length(tasks) + 1L]] <- list(
      id = paste(key, mode, sep = ":"), mode = mode, entry = entry,
      analysis = entry$analysis, context = task_context, config = task_config
    )
  }
  completed <- run_tasks(tasks, render_pathway_mode, config$workers)
  if (nrow(completed$errors)) for (row in seq_len(nrow(completed$errors))) errors[[length(errors) + 1L]] <- data.frame(
    task_id = completed$errors$task_id[[row]], mode = sub("^.*:", "", completed$errors$task_id[[row]]),
    pathway_id = sub("^[^:]+:([^:]+):.*$", "\\1", completed$errors$task_id[[row]]),
    error = completed$errors$error[[row]], stringsAsFactors = FALSE
  )

  manifest_rows <- list(funz = list(), flow = list(), pie = list(), kegg_taxa = list())
  for (id in names(completed$results)) {
    table <- completed$results[[id]]
    if (!is.data.frame(table) || !nrow(table)) next
    mode <- as.character(table$mode[[1L]])
    if (mode %in% names(manifest_rows)) manifest_rows[[mode]][[length(manifest_rows[[mode]]) + 1L]] <- table
  }


  if ("kegg_taxa" %in% mode_list && length(prepared)) {
    for (rank in config$ranks) {
      selected_entries <- lapply(c("defined", "top20"), function(selection) {
        Filter(function(entry) identical(entry$selection, selection), prepared)
      })
      names(selected_entries) <- c("defined", "top20")
      tables <- lapply(selected_entries, function(group) {
        if (!length(group)) return(NULL)
        do.call(rbind, lapply(group, function(entry) {
          build_kegg_taxa_pathway_table(
            entry$analysis$data, entry$pathway_id, entry$pathway_name,
            config$samples, rank, config$top_n_kegg_taxa %||% 10L,
            config$min_kegg_taxon_percent %||% 1
          )
        }))
      })
      present <- Filter(Negate(is.null), tables)
      if (!length(present)) next
      palette <- build_kegg_taxa_color_map(do.call(rbind, present))
      for (selection in names(tables)) {
        table <- tables[[selection]]
        if (is.null(table)) next
        task_context <- context
        task_context$selection <- selection
        tryCatch({
          rows <- render_kegg_taxa(table, task_context, rank, config$samples,
                                   config$dimensions, palette)
          manifest_rows$kegg_taxa[[length(manifest_rows$kegg_taxa) + 1L]] <- rows
        }, error = function(error) errors[[length(errors) + 1L]] <<- data.frame(
          task_id = paste(selection, "kegg_taxa", rank, sep = ":"),
          mode = "kegg_taxa", pathway_id = NA_character_,
          error = conditionMessage(error), stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (any(config$mode %in% c("huge", "normal", "taxon"))) for (count in config$counts) for (rank in config$ranks) for (dimension in config$dimensions) {
    tryCatch(render_taxonomy(
      sqm, config$samples, rank, context, NA_character_, NA_character_,
      dimension$width, dimension$height, dimension$dpi, plot_taxonomy_fn, count, dimension$name
    ), error = function(error) errors[[length(errors) + 1L]] <<- data.frame(
      task_id = paste("global", "taxon", count, rank, sep = ":"), mode = "taxon",
      pathway_id = NA_character_, error = conditionMessage(error), stringsAsFactors = FALSE
    ))
  }

  if (any(config$mode %in% c("huge", "normal", "funz"))) tryCatch({
    enzyme <- build_enzyme_table(sqm, config$samples, config$enzyme_ecs %||% NULL)
    rows <- render_enzymes(enzyme, context, config$dimensions, config$enzyme_plot_types %||% c("bar", "line"))
    if (nrow(rows)) manifest_rows$funz[[length(manifest_rows$funz) + 1L]] <- rows
  }, error = function(error) errors[[length(errors) + 1L]] <<- data.frame(
    task_id = "global:funz:enzimi", mode = "funz", pathway_id = NA_character_,
    error = conditionMessage(error), stringsAsFactors = FALSE
  ))

  if (any(config$mode %in% c("huge", "normal", "pathview"))) for (key in names(prepared)) {
    entry <- prepared[[key]]
    if (is.null(entry$analysis$pathway_sqm)) next
    task_context <- context
    task_context$selection <- entry$selection
    modes <- config$pathview_sample_modes %||% "insieme"
    for (sample_mode in modes) {
      groups <- if (sample_mode == "separato") as.list(config$samples) else list(config$samples)
      for (samples in groups) tryCatch({
        directory <- file.path(config$output_dir, "pathview", selection_dir(entry$selection), sample_mode,
                               safe_name(entry$pathway_name))
        if (sample_mode == "separato") directory <- file.path(directory, safe_name(samples))
        export_pathview_isolated(export_pathway_fn, sqm, entry$pathway_id, samples, directory)
      }, error = function(error) errors[[length(errors) + 1L]] <<- data.frame(
        task_id = paste(key, "pathview", sample_mode, paste(samples, collapse = ","), sep = ":"),
        mode = "pathview", pathway_id = entry$pathway_id, error = conditionMessage(error),
        stringsAsFactors = FALSE
      ))
    }
  }

  for (mode in names(manifest_rows)) write_mode_manifest(manifest_rows[[mode]], config$output_dir, mode)
  error_table <- if (length(errors)) do.call(rbind, errors) else empty_errors
  if (isTRUE(config$write_errors %||% TRUE)) {
    write_table(error_table, file.path(config$output_dir, "errors.tsv"))
  }
  invisible(list(results = completed$results, errors = error_table))
}

known_pathways <- c(
  `00361` = "Chlorocyclohexane and chlorobenzene degradation",
  `00622` = "Xylene degradation",
  `00623` = "Toluene degradation",
  `00621` = "Dioxin degradation",
  `00625` = "Chloroalkane and chloroalkene degradation",
  `00630` = "Glyoxylate and dicarboxylate metabolism",
  `00633` = "Nitrotoluene degradation",
  `00910` = "Nitrogen metabolism",
  `00980` = "Metabolism of xenobiotics by cytochrome P450",
  `00710` = "Carbon fixation in photosynthetic organisms"
)
default_pathway_ids <- setdiff(names(known_pathways), "00710")

default_enzyme_ecs <- c(
  "1.14.12.11", "1.14.12.12", "1.14.12.-", "3.8.1.2", "3.8.1.3", "1.13.11.-",
  "1.21.99.5", "1.14.13.243", "1.14.13.236", "1.14.13.25", "1.14.18.3",
  "1.14.99.39", "1.14.13.244", "1.14.13.7", "1.14.13.227", "1.14.13.230",
  "1.14.13.69", "2.5.1.18", "4.4.1.34", "5.2.1.2", "3.8.1.5"
)

csv <- function(value, fallback = character()) {
  if (is.null(value)) return(fallback)
  output <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
  unique(output[nzchar(output)])
}

positive_integer <- function(value, name, fallback) {
  if (is.null(value)) return(as.integer(fallback))
  if (!grepl("^[0-9]+$", value)) stop(name, " must be a positive integer.", call. = FALSE)
  result <- suppressWarnings(as.integer(value))
  if (is.na(result) || result < 1L) stop(name, " must be a positive integer.", call. = FALSE)
  result
}

bounded_percent <- function(value, name, fallback) {
  if (is.null(value)) return(fallback)
  result <- suppressWarnings(as.numeric(value))
  if (length(result) != 1L || is.na(result) || !is.finite(result) ||
      result < 0 || result > 100) {
    stop(name, " must be between 0 and 100.", call. = FALSE)
  }
  result
}

parse_dimensions_lean <- function(value, dpi) {
  labels <- csv(value, c("12x9", "16x9", "12x16"))
  lapply(labels, function(label) {
    parts <- strsplit(tolower(label), "x", fixed = TRUE)[[1L]]
    numbers <- suppressWarnings(as.numeric(parts))
    if (length(numbers) != 2L || anyNA(numbers) || any(numbers <= 0)) {
      stop("Invalid dimension: ", label, call. = FALSE)
    }
    list(name = paste0(parts[[1L]], "x", parts[[2L]]), width = numbers[[1L]],
         height = numbers[[2L]], dpi = dpi)
  })
}

build_config <- function(args) {
  allowed <- c(
    "project_dir", "output_dir", "mode", "pathways", "pathway_selection_modes",
    "pathway_top_n", "samples", "tax_mode", "top_n_ko", "top_n_taxa", "taxa",
    "taxonomy_ranks", "taxonomy_counts", "flowplot_formats", "pathview_sample_modes",
    "enzyme_ecs", "enzyme_plot_types", "dimensions", "plot_dpi", "workers",
    "top_n_kegg_taxa", "min_kegg_taxon_percent",
    "plan_only", "refresh_kegg", "help"
  )
  unknown <- setdiff(names(args), allowed)
  if (length(unknown)) stop("Unknown options: ", paste(unknown, collapse = ", "), call. = FALSE)
  for (name in c("project_dir", "output_dir", "mode")) {
    if (is.null(args[[name]]) || !nzchar(args[[name]])) stop("Missing --", name, call. = FALSE)
  }
  modes <- c("huge", "normal", "funz", "flow", "taxon", "pie", "pathview", "kegg_taxa")
  selected_modes <- csv(args$mode)
  if (!length(selected_modes) || length(setdiff(selected_modes, modes)) ||
      (any(selected_modes %in% c("huge", "normal")) && length(selected_modes) > 1L)) {
    stop("mode accepts a comma-separated subset of: ", paste(modes, collapse = ", "), call. = FALSE)
  }
  tax_mode <- args$tax_mode %||% "prokfilter"
  if (!tax_mode %in% c("prokfilter", "allfilter", "nofilter")) stop("Invalid tax_mode.", call. = FALSE)
  dpi <- suppressWarnings(as.numeric(args$plot_dpi %||% 600))
  if (length(dpi) != 1L || is.na(dpi) || dpi <= 0) stop("plot_dpi must be positive.", call. = FALSE)
  selection_modes <- csv(args$pathway_selection_modes, c("defined", "top20"))
  if (!length(selection_modes) || length(setdiff(selection_modes, c("defined", "top20")))) {
    stop("pathway_selection_modes accepts defined,top20.", call. = FALSE)
  }
  if ("normal" %in% selected_modes) selection_modes <- "defined"
  ranks <- csv(args$taxonomy_ranks, c("phylum", "class", "order", "family", "genus", "species"))
  counts <- csv(args$taxonomy_counts, c("abund", "percent"))
  formats <- csv(args$flowplot_formats, c("png", "html"))
  pathview_modes <- csv(args$pathview_sample_modes, c("insieme", "separato"))
  enzyme_plot_types <- csv(args$enzyme_plot_types, c("bar", "line"))
  if (length(setdiff(counts, c("abund", "percent")))) stop("Invalid taxonomy_counts.", call. = FALSE)
  if (length(setdiff(formats, c("png", "html")))) stop("Invalid flowplot_formats.", call. = FALSE)
  if (length(setdiff(pathview_modes, c("insieme", "separato")))) stop("Invalid pathview_sample_modes.", call. = FALSE)
  if (length(setdiff(enzyme_plot_types, c("bar", "line")))) stop("Invalid enzyme_plot_types.", call. = FALSE)
  list(
    project_dir = args$project_dir, output_dir = args$output_dir, mode = selected_modes,
    samples = if (is.null(args$samples)) NULL else csv(args$samples),
    taxa = csv(args$taxa), tax_mode = tax_mode, ranks = ranks, counts = counts,
    formats = formats, pathview_sample_modes = pathview_modes,
    dimensions = parse_dimensions_lean(args$dimensions, dpi),
    top_n_pathways = positive_integer(args$pathway_top_n, "pathway_top_n", 20L),
    top_n_ko = positive_integer(args$top_n_ko, "top_n_ko", 20L),
    top_n_taxa = positive_integer(args$top_n_taxa, "top_n_taxa", 15L),
    top_n_kegg_taxa = positive_integer(args$top_n_kegg_taxa, "top_n_kegg_taxa", 10L),
    min_kegg_taxon_percent = bounded_percent(args$min_kegg_taxon_percent, "min_kegg_taxon_percent", 1),
    workers = args$workers %||% default_workers(), plan_only = isTRUE(args$plan_only),
    refresh_kegg = isTRUE(args$refresh_kegg), selection_modes = selection_modes,
    pie_selection_modes = if (is.null(args$pathway_selection_modes)) "defined" else selection_modes,
    defined_pathways = csv(args$pathways, default_pathway_ids),
    enzyme_ecs = csv(args$enzyme_ecs, default_enzyme_ecs),
    enzyme_plot_types = enzyme_plot_types
  )
}

print_help <- function() cat(
  "Usage: Rscript sqm_plots_lean.R --project_dir PATH --output_dir PATH --mode MODE [options]\n",
  "Modes: huge (all outputs), normal (without PIE or Top 20 plots), or: funz,flow,taxon,pie,pathview,kegg_taxa\n",
  "Options:\n",
  "  --plan_only       Write contextual top20.tsv files and stop.\n",
  "  --workers=N       Parallel rendering workers; default min(4, physical cores).\n",
  "  --refresh_kegg    Refresh the KEGG pathway catalog cache.\n",
  "  --top_n_kegg_taxa=N  Maximum selected taxa per pathway/sample; default 10.\n",
  "  --min_kegg_taxon_percent=P  Minimum relative abundance; default 1%.\n",
  "  --help, -h        Show this help.\n",
  sep = ""
)


ensure_sqm_sample_metadata <- function(sqm) {
  samples <- colnames(sqm$orfs$tpm)
  if (is.null(samples) || !length(samples) || anyNA(samples) ||
      any(!nzchar(samples)) || anyDuplicated(samples)) {
    stop("ORF TPM requires unique named sample columns.", call. = FALSE)
  }
  if (!setequal(samples, colnames(sqm$functions$KEGG$tpm))) {
    stop("ORF and KEGG TPM sample columns disagree.", call. = FALSE)
  }
  metadata <- sqm$misc$samples
  if (!length(metadata)) {
    misc <- sqm$misc %||% list()
    misc$samples <- samples
    sqm$misc <- misc
  } else if (!setequal(metadata, samples)) {
    stop("SQM sample metadata disagrees with TPM columns.", call. = FALSE)
  }
  sqm
}

validate_sqm <- function(sqm, samples, ranks) {
  missing <- setdiff(c("table", "tax", "tpm"), names(sqm$orfs))
  if (length(missing)) stop("sqm$orfs missing: ", paste(missing, collapse = ", "), call. = FALSE)
  ids <- lapply(sqm$orfs[c("table", "tax", "tpm")], rownames)
  if (any(vapply(ids, is.null, logical(1L))) || any(vapply(ids, anyDuplicated, integer(1L))) ||
      !setequal(ids$table, ids$tax) || !setequal(ids$table, ids$tpm)) {
    stop("SQM ORF tables require unique matching row names.", call. = FALSE)
  }
  if (length(setdiff(samples, colnames(sqm$orfs$tpm)))) stop("Unknown samples.", call. = FALSE)
  if (length(setdiff(ranks, colnames(sqm$orfs$tax)))) stop("Unknown taxonomy ranks.", call. = FALSE)
  validate_matrix(sqm$orfs$tpm, samples, "ORF TPM")
  validate_matrix(sqm$functions$KEGG$tpm, samples, "Official KO TPM")
  invisible(sqm)
}

make_contexts <- function(sqm, config) {
  if (!length(config$taxa)) return(list(list(
    sqm = sqm, output_dir = config$output_dir,
    filtered_taxon = NA_character_, filtered_taxon_rank = NA_character_
  )))
  taxonomy <- as.data.frame(sqm$orfs$tax, check.names = FALSE)
  lapply(config$taxa, function(taxon) {
    hits <- names(taxonomy)[vapply(taxonomy, function(values) {
      any(tolower(trimws(as.character(values))) == tolower(trimws(taxon)), na.rm = TRUE)
    }, logical(1L))]
    if (length(hits) != 1L) stop("Taxon must match exactly one rank: ", taxon, call. = FALSE)
    values <- tolower(trimws(as.character(taxonomy[[hits]])))
    ids <- rownames(taxonomy)[!is.na(values) & values == tolower(trimws(taxon))]
    list(
      sqm = subset_sqm_ids(sqm, ids),
      output_dir = file.path(config$output_dir, "taxon_filter", safe_name(hits), safe_name(taxon)),
      filtered_taxon = taxon, filtered_taxon_rank = hits
    )
  })
}

download_catalog <- function() {
  lines <- readLines(url("https://rest.kegg.jp/list/pathway/ko", open = "rb"), warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  valid <- lengths(fields) >= 2L
  data.frame(
    pathway_id = sub("^(path:)?ko", "", vapply(fields[valid], `[[`, character(1L), 1L)),
    pathway_name = sub("\\s*-\\s*Reference pathway\\s*$", "", vapply(fields[valid], `[[`, character(1L), 2L), ignore.case = TRUE),
    stringsAsFactors = FALSE
  )
}

download_kgml <- function(pathway_id) {
  readLines(url(paste0("https://rest.kegg.jp/get/ko", pathway_id, "/kgml"), open = "rb"), warn = FALSE)
}

parse_kgml <- function(path) {
  nodes <- pathview::node.info(path)
  ids <- names(nodes$type)[as.character(nodes$type) == "ortholog"]
  names <- nodes$kegg.names[ids]
  sort(unique(unlist(lapply(names, function(value) extract_ko_ids(paste(value, collapse = ";"))), use.names = FALSE)))
}

required_packages <- function(config) {
  packages <- c("SQMtools", "ggplot2")
  if (!config$plan_only && any(config$mode %in% c("huge", "normal", "flow"))) packages <- c(packages, "ggalluvial")
  if (!config$plan_only && "html" %in% config$formats && any(config$mode %in% c("huge", "normal", "flow"))) {
    packages <- c(packages, "plotly", "htmlwidgets")
  }
  if (!config$plan_only) packages <- c(packages, "pathview")
  unique(packages)
}

main_impl <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || any(args %in% c("--help", "-h"))) {
    print_help()
    return(0L)
  }
  config <- build_config(parse_args(args))
  if (!dir.exists(config$project_dir)) stop("Project directory does not exist: ", config$project_dir, call. = FALSE)
  missing <- required_packages(config)[!vapply(required_packages(config), requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  # SQMtools 1.7.2 resolves its MGOGs dataset through the attached package environment.
  suppressPackageStartupMessages(library("SQMtools", character.only = TRUE))
  sqm <- SQMtools::loadSQM(
    project_path = normalizePath(config$project_dir, winslash = "/", mustWork = TRUE),
    tax_mode = config$tax_mode, trusted_functions_only = FALSE, load_sequences = FALSE
  )
  sqm <- ensure_sqm_sample_metadata(sqm)
  if (is.null(config$samples)) config$samples <- colnames(sqm$orfs$tpm)
  validate_sqm(sqm, config$samples, config$ranks)
  catalog <- get_kegg_catalog(config$output_dir, config$refresh_kegg, download_catalog)
  kgml_loader <- function(id) get_pathway_kos(
    config$output_dir, id, config$refresh_kegg, download_kgml, parse_kgml
  )
  contexts <- make_contexts(sqm, config)
  errors <- list()
  for (context in contexts) {
    current <- config
    current$output_dir <- context$output_dir
    current$filtered_taxon <- context$filtered_taxon
    current$filtered_taxon_rank <- context$filtered_taxon_rank
    current$write_errors <- FALSE
    result <- run_pipeline(context$sqm, current, catalog, kgml_loader)
    if (nrow(result$errors)) errors[[length(errors) + 1L]] <- result$errors
  }
  error_table <- if (length(errors)) do.call(rbind, errors) else data.frame(
    task_id = character(), mode = character(), pathway_id = character(), error = character()
  )
  write_table(error_table, file.path(config$output_dir, "errors.tsv"))
  if (nrow(error_table)) 1L else 0L
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  output_dir <- tryCatch(parse_args(args)$output_dir, error = function(...) NULL)
  run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", sprintf("%04x", sample.int(65535L, 1L)))
  log_path <- if (!is.null(output_dir)) file.path(output_dir, paste0(run_id, ".log")) else NULL
  log_line <- function(level, text) if (!is.null(log_path)) {
    dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    cat(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), level, text, "\n", file = log_path, append = TRUE)
  }
  tryCatch(withCallingHandlers({
    log_line("INFO", paste("args", paste(args, collapse = " ")))
    status <- main_impl(args)
    log_line("INFO", paste("status", status))
    status
  }, warning = function(warning) log_line("WARNING", conditionMessage(warning))), error = function(error) {
    log_line("ERROR", conditionMessage(error))
    message("ERROR: ", conditionMessage(error))
    1L
  })
}

if (identical(environment(), globalenv()) && sys.nframe() == 0L) {
  quit(status = main())
}
