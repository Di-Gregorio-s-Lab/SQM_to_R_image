# SqueezeMeta plotting command-line program.
#
# The script loads a validated SQM project once, derives ORF-level data from
# that object, and writes reproducible plots plus manifest files for each mode.
# CLI mode names and output directory names are kept stable for compatibility.

common_required_packages <- c(
  "SQMtools", "readr", "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "glue", "purrr", "scales"
)

required_packages_for_mode <- function(mode, flowplot_formats = c("png", "html")) {
  required <- common_required_packages
  if (mode %in% c("all", "flow")) {
    required <- c(required, "ggalluvial")
    if ("html" %in% tolower(as.character(flowplot_formats))) {
      required <- c(required, "plotly", "htmlwidgets")
    }
  }
  if (mode %in% c("all", "pie")) {
    required <- c(required, "forcats", "rlang")
  }
  if (mode %in% c("all", "flow", "funz", "pie", "taxon", "pathview")) {
    required <- c(required, "pathview")
  }
  unique(required)
}

check_flow_html_preflight <- function(
    mode,
    flowplot_formats,
    pandoc_available_fn = function() TRUE) {
  invisible(TRUE)
}

check_required_packages <- function(required, mode, availability_fn = requireNamespace) {
  available <- vapply(
    required,
    function(package_name) availability_fn(package_name, quietly = TRUE),
    logical(1)
  )
  missing <- required[!available]
  if (length(missing) > 0L) {
    stop(
      "Missing required R packages for mode ", mode, ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(required)
}

bootstrap_cli_value <- function(args, option_name) {
  equals_prefix <- paste0("--", option_name, "=")
  equals_match <- startsWith(args, equals_prefix)
  if (any(equals_match)) {
    return(sub(equals_prefix, "", args[which(equals_match)[[1L]]], fixed = TRUE))
  }
  separate_match <- which(args == paste0("--", option_name))
  if (length(separate_match) > 0L && separate_match[[1L]] < length(args)) {
    candidate <- args[[separate_match[[1L]] + 1L]]
    if (!startsWith(candidate, "--")) {
      return(candidate)
    }
  }
  NULL
}

direct_cli_execution <- sys.nframe() == 0L
bootstrap_args <- commandArgs(trailingOnly = TRUE)
bootstrap_help <- direct_cli_execution &&
  (length(bootstrap_args) == 0L || any(bootstrap_args %in% c("--help", "-h")))

if (!bootstrap_help && !direct_cli_execution) {
  bootstrap_mode <- if (direct_cli_execution) bootstrap_cli_value(bootstrap_args, "mode") else NULL
  if (is.null(bootstrap_mode) || !nzchar(bootstrap_mode)) {
    bootstrap_mode <- if (direct_cli_execution) "unspecified" else "source"
  }
  bootstrap_format_value <- if (direct_cli_execution) {
    bootstrap_cli_value(bootstrap_args, "flowplot_formats")
  } else {
    NULL
  }
  bootstrap_formats <- if (is.null(bootstrap_format_value)) {
    c("png", "html")
  } else {
    trimws(strsplit(bootstrap_format_value, ",", fixed = TRUE)[[1L]])
  }
  bootstrap_required <- if (direct_cli_execution) {
    required_packages_for_mode(bootstrap_mode, bootstrap_formats)
  } else {
    common_required_packages
  }
  check_required_packages(bootstrap_required, bootstrap_mode)
  suppressPackageStartupMessages(
    invisible(lapply(common_required_packages, library, character.only = TRUE))
  )
}

# Shared palette used by every plot to keep categories visually consistent.
colors_hex <- c(
  "#5d8aa8", "#e32636", "#efdecd", "#ffbf00",
  "#9966cc", "#a4c639", "#cd9575", "#915c83",
  "#008000", "#fbceb1", "#00ffff",
  "#4b5320", "#b2beb5", "#87a96b", "#ff9966", "#a52a2a",
  "#6e7f80", "#ff2052", "#007fff", "#f0ffff", "#89cff0",
  "#f4c2c2", "#21abcd", "#fae7b5", "#ffe135", "#848482",
  "#98777b", "#f5f5dc", "#3d2b1f",
  "#fe6f5e", "#000000", "#ffebcd", "#318ce7", "#ace5ee", "#faf0be",
  "#0000ff", "#a2a2d0", "#6699cc", "#0d98ba", "#8a2be2", "#8a2be2",
  "#de5d83", "#79443b", "#0095b6", "#e3dac9", "#cc0000", "#006a4e",
  "#873260", "#0070ff", "#b5a642", "#cb4154", "#1dacd6", "#66ff00",
  "#bf94e4", "#c32148", "#ff007f", "#08e8de", "#d19fe8", "#f4bbff",
  "#ff55a3", "#fb607f", "#004225", "#cd7f32", "#a52a2a", "#ffc1cc",
  "#e7feff", "#f0dc82"
)

default_dimensions <- list(
  "12x9" = c(width = 12, height = 9),
  "16x9" = c(width = 16, height = 9),
  "12x16" = c(width = 12, height = 16)
)

default_taxonomy_ranks <- c("phylum", "class", "order", "family", "genus", "species")
default_enzyme_ecs <- c(
  "1.14.12.11", "1.14.12.12", "1.14.12.-",
  "3.8.1.2", "3.8.1.3", "1.13.11.-",
  "1.21.99.5", "1.14.13.243", "1.14.13.236", "1.14.13.25",
  "1.14.18.3", "1.14.99.39", "1.14.13.244", "1.14.13.7",
  "1.14.13.227", "1.14.13.230", "1.14.13.69", "2.5.1.18",
  "4.4.1.34", "5.2.1.2", "3.8.1.5"
)
default_enzyme_plot_types <- c("bar", "line")
all_taxonomy_columns <- c(
  "superkingdom", "phylum", "class", "order", "family", "genus", "species"
)

# Curated pathway IDs supported by the explicit `defined` selection mode.
known_pathways <- data.frame(
  pathway_id = c(
    "00361", "00710", "00623", "00621", "00625",
    "00630", "00633", "00910", "00980"
  ),
  canonical_pathway_name = c(
    "Chlorocyclohexane and chlorobenzene degradation",
    "Carbon fixation in photosynthetic organisms",
    "Toluene degradation",
    "Dioxin degradation",
    "Chloroalkane and chloroalkene degradation",
    "Glyoxylate and dicarboxylate metabolism",
    "Nitrotoluene degradation",
    "Nitrogen metabolism",
    "Metabolism of xenobiotics by cytochrome P450"
  ),
  stringsAsFactors = FALSE
)
pathway_name_aliases <- data.frame(
  pathway_id = "00720",
  alias_name = "Carbon fixation pathways in prokaryotes",
  current_kegg_name = "Other carbon fixation pathways",
  stringsAsFactors = FALSE
)
default_pathway_ids <- known_pathways$pathway_id
default_pathway_selection_modes <- c("defined", "top20")
default_pathway_top_n <- 20L
kegg_pathway_roots <- c(
  "Metabolism",
  "Genetic Information Processing",
  "Environmental Information Processing",
  "Cellular Processes",
  "Organismal Systems",
  "Human Diseases"
)

# ---- Command-line parsing and shared configuration -----------------------

parse_named_args <- function(args) {
  named <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("--help", "-h")) {
      named[["help"]] <- "TRUE"
      i <- i + 1L
      next
    }

    if (grepl("^--[^=]+=.*", arg)) {
      parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
      if (length(parts) != 2L || !nzchar(parts[[1]])) {
        stop("Invalid argument: ", arg, call. = FALSE)
      }
      named[[parts[[1]]]] <- parts[[2]]
      i <- i + 1L
      next
    }

    if (grepl("^--[^=]+$", arg)) {
      if (i == length(args)) {
        stop("Missing value for argument: ", arg, call. = FALSE)
      }
      if (grepl("^--", args[[i + 1L]])) {
        stop("Missing value for argument: ", arg, call. = FALSE)
      }
      key <- sub("^--", "", arg)
      named[[key]] <- args[[i + 1L]]
      i <- i + 2L
      next
    }

    stop("Invalid argument: ", arg, call. = FALSE)
  }

  named
}

print_help <- function() {
  cat(
    "Usage: Rscript sqm_plots.R --project_dir PATH --output_dir PATH --mode MODE [options]\n",
    "\n",
    "Required arguments:\n",
    "  --project_dir=PATH          SqueezeMeta project directory.\n",
    "  --output_dir=PATH           Output directory.\n",
    "  --mode=MODE                 Mode: all, flow, funz, enzimi, taxon, pathview, pie.\n",
    "\n",
    "Optional arguments:\n",
    "  --pathways=LIST            Comma-separated KEGG pathways.\n",
    "                             Accepts numeric codes (for example, 00361) or full names.\n",
    "  --pathway_selection_modes=LIST Default: defined,top20.\n",
    "                             top20 ranks only three-level KEGG PATHWAY entries.\n",
    "                             With --taxa, top20 is recalculated inside each taxon context.\n",
    "  --pathway_top_n=NUM        Number of KEGG pathways ranked in top20; default: 20.\n",
    "  --samples=LIST             Comma-separated samples.\n",
    "  --tax_mode=MODE            Default: prokfilter.\n",
    "  --top_n_ko=NUM             Default: 20; applies to FUNZ and FLOW, not PIE.\n",
    "  --top_n_taxa=NUM           Default: 15.\n",
    "  --taxa=LIST                Taxa to use as a filter.\n",
    "  --taxonomy_ranks=LIST      Default: phylum,class,order,family,genus,species.\n",
    "  --taxonomy_counts=LIST     Default: abund,percent.\n",
    "  --flowplot_formats=LIST    Default: png,html.\n",
    "  --pathview_sample_modes=LIST Default: insieme,separato.\n",
    "                             Allowed values: insieme,separato.\n",
    "  --enzyme_ecs=LIST          Default: ", paste(default_enzyme_ecs, collapse = ","), ".\n",
    "  --enzyme_plot_types=LIST   Default: bar,line.\n",
    "  --dimensions=LIST          Default: 12x9,16x9,12x16.\n",
    "  --plot_dpi=NUM             Default: 600.\n",
    "  --help, -h                 Show this help message.\n",
    "\n",
    "Pie chart notes:\n",
    "  Pie mode creates one plot for every positive sample x KO x rank combination.\n",
    "  Pie mode does not apply --top_n_ko; its manifest records all_positive_ko.\n",
    "  By default pie charts use only defined pathways; top20 requires --pathway_selection_modes=top20.\n",
    "\n",
    "Known pathway codes:\n",
    paste0(
      "  ",
      known_pathways$pathway_id,
      " = ",
      known_pathways$canonical_pathway_name,
      "\n",
      collapse = ""
    ),
    sep = ""
  )
}

split_csv_arg <- function(x) {
  values <- str_split(x, ",\\s*")[[1]]
  values <- trimws(values)
  values[nzchar(values)]
}

normalize_nonempty_cli_list <- function(value, default, option_name) {
  values <- if (is.null(value)) {
    default
  } else if (length(value) == 1L) {
    split_csv_arg(value)
  } else {
    trimws(as.character(value))
  }
  values <- unique(values[nzchar(values)])

  if (length(values) == 0L) {
    stop(option_name, " cannot be empty.", call. = FALSE)
  }

  values
}

normalize_taxonomy_counts <- function(value = NULL) {
  allowed_counts <- c("abund", "percent")
  counts <- normalize_nonempty_cli_list(
    value,
    default = allowed_counts,
    option_name = "taxonomy_counts"
  )
  if (!all(counts %in% allowed_counts)) {
    stop("taxonomy_counts must contain only abund and/or percent.", call. = FALSE)
  }
  counts
}

normalize_flowplot_formats <- function(value = NULL) {
  allowed_formats <- c("png", "html")
  formats <- normalize_nonempty_cli_list(
    value,
    default = allowed_formats,
    option_name = "flowplot_formats"
  )
  if (!all(formats %in% allowed_formats)) {
    stop("flowplot_formats must contain only png and/or html.", call. = FALSE)
  }
  formats
}

normalize_pathview_sample_modes <- function(value = NULL) {
  allowed_modes <- c("insieme", "separato")
  modes <- if (is.null(value)) {
    allowed_modes
  } else if (length(value) == 1L) {
    split_csv_arg(value)
  } else {
    trimws(as.character(value))
  }
  modes <- unique(modes[nzchar(modes)])

  if (length(modes) == 0L || !all(modes %in% allowed_modes)) {
    stop(
      "pathview_sample_modes must contain only insieme and/or separato.",
      call. = FALSE
    )
  }

  modes
}

normalize_enzyme_ecs <- function(value = NULL) {
  ecs <- if (is.null(value)) {
    default_enzyme_ecs
  } else if (length(value) == 1L) {
    split_csv_arg(value)
  } else {
    trimws(as.character(value))
  }
  ecs <- unique(ecs[nzchar(ecs)])
  valid_ec_pattern <- "^[0-9]+\\.[0-9]+\\.[0-9]+\\.(?:[0-9]+|-)$"

  if (length(ecs) == 0L || !all(grepl(valid_ec_pattern, ecs))) {
    stop(
      "enzyme_ecs must contain EC codes in the 1.2.3.4 or 1.2.3.- format.",
      call. = FALSE
    )
  }

  ecs
}

normalize_enzyme_plot_types <- function(value = NULL) {
  allowed_types <- c("bar", "line")
  plot_types <- if (is.null(value)) {
    default_enzyme_plot_types
  } else if (length(value) == 1L) {
    split_csv_arg(value)
  } else {
    trimws(as.character(value))
  }
  plot_types <- unique(plot_types[nzchar(plot_types)])

  if (length(plot_types) == 0L || !all(plot_types %in% allowed_types)) {
    stop("enzyme_plot_types must contain only bar and/or line.", call. = FALSE)
  }

  plot_types
}

normalize_pathway_selection_modes <- function(value = NULL) {
  allowed_modes <- c("defined", "top20")
  modes <- if (is.null(value)) {
    default_pathway_selection_modes
  } else if (length(value) == 1L) {
    split_csv_arg(value)
  } else {
    trimws(as.character(value))
  }
  modes <- unique(modes[nzchar(modes)])

  if (length(modes) == 0L || !all(modes %in% allowed_modes)) {
    stop("pathway_selection_modes must contain only defined and/or top20.", call. = FALSE)
  }

  modes
}

pie_pathway_selection_modes <- function(pathway_selection_modes, explicitly_requested = FALSE) {
  if (isTRUE(explicitly_requested)) {
    return(pathway_selection_modes)
  }
  "defined"
}

pathway_selection_directory <- function(pathway_selection) {
  selection_dirs <- c(defined = "definiti", top20 = "top20")
  if (!pathway_selection %in% names(selection_dirs)) {
    stop("Invalid pathway selection: ", pathway_selection, call. = FALSE)
  }
  unname(selection_dirs[[pathway_selection]])
}

pathview_is_exportable <- function(pathway_selection, pathway_id) {
  length(pathway_id) == 1L &&
    !is.na(pathway_id) &&
    isTRUE(grepl("^[0-9]{5}$", pathway_id))
}

# Normalize user-derived labels before using them as directory or file names.
sanitize_name <- function(x) {
  x |>
    str_replace_all("[^A-Za-z0-9_.-]", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("(^_+|_+$)", "")
}

# Preserve already-portable labels and add a stable suffix whenever sanitizing
# could otherwise collapse distinct user values to the same path component.
safe_output_component <- function(value, max_length = 80L) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("Output path components must be one non-empty value.", call. = FALSE)
  }
  if (
    length(max_length) != 1L || is.na(max_length) || !is.finite(max_length) ||
      max_length < 15 || max_length != as.integer(max_length)
  ) {
    stop("max_length must be one integer of at least 15.", call. = FALSE)
  }
  max_length <- as.integer(max_length)
  portable_unchanged <- grepl("^[A-Za-z0-9_.-]+$", value)
  if (portable_unchanged && nchar(value, type = "chars") <= max_length) {
    return(value)
  }
  sanitized <- sanitize_name(value)
  if (!nzchar(sanitized)) {
    sanitized <- "item"
  }
  suffix <- paste0("__", stable_path_token(value))
  prefix_length <- max_length - nchar(suffix, type = "chars")
  paste0(substr(sanitized, 1L, prefix_length), suffix)
}

format_dimension_label <- function(x) {
  str_replace(format(x, trim = TRUE, scientific = FALSE), "\\.0+$", "")
}

# Parse named width-by-height presets once so all exports share the same sizes.
parse_dimensions <- function(named_args) {
  if (is.null(named_args$dimensions)) {
    return(default_dimensions)
  }

  labels <- split_csv_arg(named_args$dimensions)
  if (length(labels) == 0L) {
    stop("dimensions cannot be empty.", call. = FALSE)
  }
  parsed <- map(labels, function(label) {
    parts <- str_split(label, "x", simplify = TRUE)
    if (ncol(parts) != 2L) {
      stop("Invalid dimension: ", label, call. = FALSE)
    }
    width <- as.numeric(parts[[1]])
    height <- as.numeric(parts[[2]])
    if (
      is.na(width) || is.na(height) ||
        !is.finite(width) || !is.finite(height) ||
        width <= 0 || height <= 0
    ) {
      stop("Invalid dimension: ", label, call. = FALSE)
    }
    c(width = width, height = height)
  })

  setNames(parsed, labels)
}

relative_to_output <- function(path, output_dir) {
  path_abs <- normalizePath(path, winslash = "/", mustWork = FALSE)
  out_abs <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(out_abs, "/")
  if (startsWith(path_abs, prefix)) {
    substr(path_abs, nchar(prefix) + 1L, nchar(path_abs))
  } else {
    path_abs
  }
}

# ---- Output and manifest helpers -----------------------------------------

.sqm_run_context <- new.env(parent = emptyenv())

clear_run_context <- function() {
  context_fields <- ls(envir = .sqm_run_context, all.names = TRUE)
  if (length(context_fields) > 0L) {
    rm(list = context_fields, envir = .sqm_run_context)
  }
  invisible(TRUE)
}

validate_run_id <- function(run_id) {
  run_id <- as.character(run_id)
  pattern <- "^[0-9]{8}T[0-9]{6}_UTC[pm][0-9]{4}_[0-9a-f]{4}$"
  if (length(run_id) != 1L || is.na(run_id) || !grepl(pattern, run_id)) {
    stop("Invalid run_id: ", paste(run_id, collapse = ","), call. = FALSE)
  }
  run_id
}

generate_collision_suffix <- function() {
  paste0(sample(c(0:9, letters[1:6]), 4L, replace = TRUE), collapse = "")
}

generate_run_id <- function(
    now = Sys.time(),
    utc_offset = format(now, "%z"),
    collision_suffix = generate_collision_suffix()) {
  utc_offset <- as.character(utc_offset)
  collision_suffix <- tolower(as.character(collision_suffix))
  if (length(utc_offset) != 1L || !grepl("^[+-][0-9]{4}$", utc_offset)) {
    stop("UTC offset must use +HHMM or -HHMM format.", call. = FALSE)
  }
  if (length(collision_suffix) != 1L || !grepl("^[0-9a-f]{4}$", collision_suffix)) {
    stop("Run collision suffix must contain four hexadecimal characters.", call. = FALSE)
  }
  offset_label <- paste0(
    "UTC",
    if (startsWith(utc_offset, "+")) "p" else "m",
    substring(utc_offset, 2L)
  )
  validate_run_id(paste0(
    format(now, "%Y%m%dT%H%M%S"),
    "_", offset_label, "_", collision_suffix
  ))
}

run_id_exists <- function(output_dir, run_id) {
  if (!dir.exists(output_dir)) {
    return(FALSE)
  }
  existing_files <- list.files(
    output_dir,
    recursive = TRUE,
    full.names = FALSE,
    all.files = TRUE,
    include.dirs = FALSE
  )
  run_id <- validate_run_id(run_id)
  any(
    basename(existing_files) == paste0(run_id, ".log") |
      grepl(paste0("__", run_id), basename(existing_files), fixed = TRUE)
  )
}

allocate_run_id <- function(
    output_dir,
    now = Sys.time(),
    suffix_fn = generate_collision_suffix,
    max_attempts = 1024L) {
  for (attempt in seq_len(max_attempts)) {
    candidate <- generate_run_id(
      now = now,
      utc_offset = format(now, "%z"),
      collision_suffix = suffix_fn()
    )
    if (!run_id_exists(output_dir, candidate)) {
      return(candidate)
    }
  }
  stop("Unable to allocate a unique run_id in output_dir.", call. = FALSE)
}

set_run_context <- function(
    run_id,
    started_at = Sys.time(),
    sample_order_basis = NA_character_,
    samples = character(),
    cli_args = character()) {
  clear_run_context()
  .sqm_run_context$run_id <- validate_run_id(run_id)
  .sqm_run_context$started_at <- started_at
  .sqm_run_context$sample_order_basis <- as.character(sample_order_basis)
  .sqm_run_context$samples <- as.character(samples)
  .sqm_run_context$cli_args <- as.character(cli_args)
  .sqm_run_context$warnings <- character()
  .sqm_run_context$completed_artifacts <- character()
  .sqm_run_context$completed_manifests <- character()
  .sqm_run_context$active_manifest_registry <- NULL
  .sqm_run_context$log_path <- NA_character_
  .sqm_run_context$pathway_skips <- data.frame(
    context = character(),
    pathway = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  invisible(.sqm_run_context$run_id)
}

current_run_id <- function(required = FALSE) {
  run_id <- .sqm_run_context$run_id
  if (is.null(run_id)) {
    if (isTRUE(required)) {
      stop("No analytical run context is active.", call. = FALSE)
    }
    return(NA_character_)
  }
  run_id
}

current_sample_order_basis <- function() {
  basis <- .sqm_run_context$sample_order_basis
  if (is.null(basis) || length(basis) != 1L) NA_character_ else basis
}

update_run_sample_order <- function(samples, sample_order_basis) {
  if (!is.na(current_run_id())) {
    .sqm_run_context$samples <- as.character(samples)
    .sqm_run_context$sample_order_basis <- as.character(sample_order_basis)
  }
  invisible(samples)
}

record_run_warning <- function(message) {
  if (!is.na(current_run_id())) {
    .sqm_run_context$warnings <- unique(c(
      .sqm_run_context$warnings,
      as.character(message)
    ))
    append_run_log("WARN", message)
  }
  invisible(message)
}

record_pathway_skips <- function(skips) {
  if (nrow(skips) == 0L || is.na(current_run_id())) {
    return(invisible(skips))
  }
  required_columns <- c("context", "pathway", "reason")
  if (!all(required_columns %in% colnames(skips))) {
    stop("Pathway skip audit is missing required columns.", call. = FALSE)
  }
  .sqm_run_context$pathway_skips <- unique(rbind(
    .sqm_run_context$pathway_skips,
    as.data.frame(skips[required_columns], stringsAsFactors = FALSE)
  ))
  apply(skips[required_columns], 1L, function(skip) {
    append_run_log(
      "SKIP",
      paste0(
        "context=", skip[["context"]],
        " pathway=", skip[["pathway"]],
        " reason=", skip[["reason"]]
      )
    )
  })
  invisible(skips)
}

add_run_id_to_path <- function(path, run_id) {
  run_id <- validate_run_id(run_id)
  extension <- tools::file_ext(path)
  stem <- if (nzchar(extension)) tools::file_path_sans_ext(path) else path
  suffix <- paste0("__", run_id)
  if (endsWith(stem, suffix)) {
    return(path)
  }
  if (nzchar(extension)) {
    paste0(stem, suffix, ".", extension)
  } else {
    paste0(stem, suffix)
  }
}

run_artifact_path <- function(path) {
  path
}

append_run_log <- function(level, message) {
  log_path <- .sqm_run_context$log_path
  if (is.null(log_path) || length(log_path) != 1L || is.na(log_path) || !nzchar(log_path)) {
    return(invisible(FALSE))
  }
  line <- paste0(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    " [", as.character(level), "] ",
    paste0(as.character(message), collapse = "")
  )
  cat(line, "\n", file = log_path, append = TRUE, sep = "")
  invisible(TRUE)
}

record_completed_artifact <- function(path) {
  if (!is.na(current_run_id())) {
    normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
    .sqm_run_context$completed_artifacts <- unique(c(
      .sqm_run_context$completed_artifacts,
      normalized
    ))
  }
  invisible(path)
}

record_completed_manifest <- function(path) {
  if (!is.na(current_run_id())) {
    normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
    .sqm_run_context$completed_manifests <- unique(c(
      .sqm_run_context$completed_manifests,
      normalized
    ))
  }
  invisible(path)
}

write_tsv_contents <- function(
    data,
    path,
    readr_available = requireNamespace("readr", quietly = TRUE)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(readr_available)) {
    readr::write_tsv(data, path, na = "NA")
  } else {
    utils::write.table(
      as.data.frame(data, stringsAsFactors = FALSE),
      file = path,
      sep = "\t",
      quote = TRUE,
      qmethod = "double",
      row.names = FALSE,
      col.names = TRUE,
      na = "NA",
      fileEncoding = "UTF-8"
    )
  }
  path
}

write_tsv_safe <- function(
    data,
    path,
    readr_available = requireNamespace("readr", quietly = TRUE)) {
  write_tsv_contents(data, path, readr_available = readr_available)
  record_completed_artifact(path)
  path
}

progress_message <- function(..., .prefix = "[sqm_plots]") {
  text <- paste0(..., collapse = "")
  message(.prefix, " ", text)
  append_run_log("INFO", paste0(.prefix, " ", text))
}

installed_package_version <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    return("not-installed")
  }
  as.character(utils::packageVersion(package))
}

initialize_run_log <- function(output_dir, project_dir, mode, tax_mode) {
  run_id <- current_run_id(required = TRUE)
  if (!dir.exists(output_dir)) {
    stop("Unable to initialize run log; output_dir does not exist: ", output_dir, call. = FALSE)
  }
  log_path <- file.path(output_dir, paste0(run_id, ".log"))
  if (file.exists(log_path)) {
    stop("Refusing to overwrite an existing run log: ", log_path, call. = FALSE)
  }

  header <- c(
    paste0("RUN_ID=", run_id),
    paste0("STARTED_AT=", format(.sqm_run_context$started_at, "%Y-%m-%dT%H:%M:%S%z")),
    paste0("PROJECT_DIR=", as.character(project_dir)),
    paste0("OUTPUT_DIR=", normalizePath(output_dir, winslash = "/", mustWork = TRUE)),
    paste0("MODE=", as.character(mode)),
    paste0("TAX_MODE=", as.character(tax_mode)),
    paste0("CLI_ARGUMENTS=", paste(.sqm_run_context$cli_args, collapse = " ")),
    paste0("R_VERSION=", as.character(getRversion())),
    paste0("SQMTOOLS_VERSION=", installed_package_version("SQMtools")),
    paste0("PATHVIEW_VERSION=", installed_package_version("pathview"))
  )
  writeLines(header, log_path, useBytes = TRUE)
  if (!file.exists(log_path) || is.na(file.info(log_path)$size)) {
    stop("Unable to create run log: ", log_path, call. = FALSE)
  }
  .sqm_run_context$log_path <- log_path
  append_run_log("INFO", "Run log initialized")
  log_path
}

record_legacy_outputs <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    return(invisible(character()))
  }
  files <- list.files(
    output_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE
  )
  relative <- vapply(files, relative_to_output, character(1), output_dir = output_dir)
  normalized <- gsub("\\\\", "/", relative)
  legacy_manifest <- grepl(
    "(^|/)(manifest_(all|run|failed_artifacts|taxon|pathview)(__[0-9]{8}T[0-9]{6}_UTC[pm][0-9]{4}_[0-9a-f]{4})?\\.tsv)$",
    normalized,
    perl = TRUE
  )
  legacy_run_suffix <- grepl(
    "__[0-9]{8}T[0-9]{6}_UTC[pm][0-9]{4}_[0-9a-f]{4}(\\.[^/]+)?$",
    basename(normalized),
    perl = TRUE
  )
  legacy <- sort(unique(normalized[legacy_manifest | legacy_run_suffix]))
  if (length(legacy) > 0L) {
    preview <- paste(utils::head(legacy, 10L), collapse = " | ")
    suffix <- if (length(legacy) > 10L) paste0(" | ... ", length(legacy) - 10L, " more") else ""
    append_run_log(
      "LEGACY",
      paste0(length(legacy), " legacy manifest/artifact(s) left unchanged: ", preview, suffix)
    )
  }
  invisible(legacy)
}

finalize_run_log <- function(status, error_message = NA_character_) {
  status <- toupper(as.character(status))
  if (!status %in% c("SUCCESS", "FAILED")) {
    stop("Run log status must be SUCCESS or FAILED.", call. = FALSE)
  }
  if (identical(status, "FAILED") && !is.na(error_message) && nzchar(error_message)) {
    append_run_log("ERROR", error_message)
  }
  append_run_log(
    "SUMMARY",
    paste0(
      "artifacts=", length(.sqm_run_context$completed_artifacts),
      " manifests=", length(.sqm_run_context$completed_manifests),
      " warnings=", length(.sqm_run_context$warnings),
      " skipped_pathways=", nrow(.sqm_run_context$pathway_skips)
    )
  )
  finished_at <- Sys.time()
  append_run_log("STATUS", paste0("STATUS=", status))
  append_run_log("STATUS", paste0("FINISHED_AT=", format(finished_at, "%Y-%m-%dT%H:%M:%S%z")))
  append_run_log(
    "STATUS",
    paste0("DURATION_SECONDS=", round(as.numeric(difftime(
      finished_at,
      .sqm_run_context$started_at,
      units = "secs"
    )), 3L))
  )
  invisible(.sqm_run_context$log_path)
}

stable_path_token <- function(value) {
  bytes <- as.integer(charToRaw(enc2utf8(as.character(value))))
  rolling_hash <- function(seed, multiplier, modulus) {
    hash_value <- seed
    for (byte in bytes) {
      hash_value <- (hash_value * multiplier + byte + 1) %% modulus
    }
    sprintf("%08x", as.integer(hash_value))
  }

  paste0(
    substr(rolling_hash(17, 131, 2147483629), 1L, 6L),
    substr(rolling_hash(23, 137, 2147483587), 1L, 6L)
  )
}

portable_png_output_path <- function(
    output_dir,
    file_stem,
    dimension_name,
    max_path_length = 240L) {
  if (length(max_path_length) != 1L || is.na(max_path_length) || max_path_length <= 0) {
    stop("max_path_length must be one positive number.", call. = FALSE)
  }

  output_dir_abs <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  logical_filename <- paste0(file_stem, "_", dimension_name, ".png")
  logical_path <- file.path(output_dir, logical_filename)
  logical_path_abs <- normalizePath(logical_path, winslash = "/", mustWork = FALSE)
  if (nchar(logical_path_abs, type = "chars") <= max_path_length) {
    return(logical_path)
  }

  readable_stem <- sanitize_name(file_stem)
  if (!nzchar(readable_stem)) {
    readable_stem <- "plot"
  }
  token <- stable_path_token(logical_filename)
  semantic_suffix <- paste0("_", dimension_name, ".png")
  fixed_length <- 2L + nchar(token, type = "chars") + nchar(semantic_suffix, type = "chars")
  filename_budget <- as.integer(max_path_length) - nchar(output_dir_abs, type = "chars") - 1L
  prefix_budget <- filename_budget - fixed_length
  if (prefix_budget < 1L) {
    stop(
      "output_dir is too long for a portable PNG filename; choose a shorter output_dir: ",
      output_dir_abs,
      call. = FALSE
    )
  }

  compact_filename <- paste0(
    substr(readable_stem, 1L, prefix_budget),
    "__", token, semantic_suffix
  )
  compact_path <- file.path(output_dir, compact_filename)
  compact_path_abs <- normalizePath(compact_path, winslash = "/", mustWork = FALSE)
  if (nchar(compact_path_abs, type = "chars") > max_path_length) {
    stop(
      "output_dir is too long for the compact PNG filename; choose a shorter output_dir: ",
      output_dir_abs,
      call. = FALSE
    )
  }

  progress_message("Compacted PNG filename: ", logical_filename, " -> ", compact_filename)
  compact_path
}

assert_output_artifact <- function(path, label = "Output artifact") {
  if (!file.exists(path)) {
    stop(label, " was not created at the requested path: ", path, call. = FALSE)
  }
  artifact_info <- suppressWarnings(file.info(path))
  if (isTRUE(artifact_info$isdir[[1L]]) || is.na(artifact_info$size[[1L]]) || artifact_info$size[[1L]] <= 0) {
    stop(label, " is not a non-empty regular file: ", path, call. = FALSE)
  }
  record_completed_artifact(path)
  invisible(path)
}

save_png_dimensions <- function(
    plot_object,
    output_dir,
    file_stem,
    dimensions,
    dpi,
    max_path_length = 240L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_files <- character()
  for (dim_name in names(dimensions)) {
    dims <- dimensions[[dim_name]]
    output_file <- portable_png_output_path(
      output_dir,
      file_stem,
      dim_name,
      max_path_length = max_path_length
    )
    ggplot2::ggsave(
      filename = output_file,
      plot = plot_object,
      width = dims[["width"]],
      height = dims[["height"]],
      units = "in",
      dpi = dpi,
      bg = "white"
    )
    assert_output_artifact(output_file, "PNG output")
    output_files[[dim_name]] <- output_file
  }
  output_files
}

save_html_widget <- function(widget, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  dependency_dir <- paste0(tools::file_path_sans_ext(path), "_files")
  htmlwidgets::saveWidget(
    widget,
    file = path,
    selfcontained = FALSE,
    libdir = basename(dependency_dir)
  )
  assert_output_artifact(path, "HTML output")
  dependency_files <- list.files(
    dependency_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  if (!dir.exists(dependency_dir) || length(dependency_files) == 0L) {
    stop(
      "HTML support directory was not created or is empty: ",
      dependency_dir,
      call. = FALSE
    )
  }
  path
}

normalize_manifest_ko_audit <- function(ko_audit = NULL) {
  required_fields <- c(
    "input_orf_count",
    "excluded_orfs_without_ko",
    "multi_ko_orf_count",
    "orf_ko_association_count",
    "multi_ko_policy",
    "ko_denominator_basis"
  )
  if (is.null(ko_audit)) {
    return(list(
      input_orf_count = NA_integer_,
      excluded_orfs_without_ko = NA_integer_,
      multi_ko_orf_count = NA_integer_,
      orf_ko_association_count = NA_integer_,
      multi_ko_policy = NA_character_,
      ko_denominator_basis = NA_character_
    ))
  }
  if (!is.data.frame(ko_audit) || nrow(ko_audit) != 1L) {
    stop("KO manifest audit must be a one-row table.", call. = FALSE)
  }
  missing_fields <- setdiff(required_fields, colnames(ko_audit))
  if (length(missing_fields) > 0L) {
    stop(
      "KO manifest audit is missing fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    input_orf_count = as.integer(ko_audit$input_orf_count[[1L]]),
    excluded_orfs_without_ko = as.integer(ko_audit$excluded_orfs_without_ko[[1L]]),
    multi_ko_orf_count = as.integer(ko_audit$multi_ko_orf_count[[1L]]),
    orf_ko_association_count = as.integer(ko_audit$orf_ko_association_count[[1L]]),
    multi_ko_policy = as.character(ko_audit$multi_ko_policy[[1L]]),
    ko_denominator_basis = as.character(ko_audit$ko_denominator_basis[[1L]])
  )
}

new_manifest_row <- function(
    script_name,
    project_dir,
    tax_mode,
    pathway,
    ec_code = NA_character_,
    ko_id = NA_character_,
    samples,
    metric,
    top_n_taxa,
    top_n_ko,
    output_type,
    output_file,
    mode,
    rank = NA_character_,
    count = NA_character_,
    format = NA_character_,
    width = NA_real_,
    height = NA_real_,
    dpi = NA_real_,
    output_scope = NA_character_,
    filtered_taxon = NA_character_,
    filtered_taxon_rank = NA_character_,
    pathway_id = NA_character_,
    ko_selection_policy = NA_character_,
    taxonomy_display_policy = NA_character_,
    denominator_type = NA_character_,
    source_data_file = NA_character_,
    sample_order_basis = NA_character_,
    ko_audit = NULL) {
  ko_audit_values <- normalize_manifest_ko_audit(ko_audit)
  if (length(sample_order_basis) != 1L || is.na(sample_order_basis)) {
    sample_order_basis <- current_sample_order_basis()
  }
  tibble::tibble(
    run_id = current_run_id(),
    script = script_name,
    project_dir = project_dir,
    tax_mode = tax_mode,
    pathway = pathway,
    ec_code = ec_code,
    ko_id = ko_id,
    pathway_id = pathway_id,
    samples = paste(samples, collapse = ","),
    metric = metric,
    top_n_taxa = top_n_taxa,
    top_n_ko = top_n_ko,
    output_type = output_type,
    output_file = output_file,
    mode = mode,
    rank = rank,
    count = count,
    format = format,
    width = width,
    height = height,
    dpi = dpi,
    output_scope = output_scope,
    filtered_taxon = filtered_taxon,
    filtered_taxon_rank = filtered_taxon_rank,
    ko_selection_policy = ko_selection_policy,
    taxonomy_display_policy = taxonomy_display_policy,
    denominator_type = denominator_type,
    source_data_file = source_data_file,
    sample_order_basis = sample_order_basis,
    input_orf_count = ko_audit_values$input_orf_count,
    excluded_orfs_without_ko = ko_audit_values$excluded_orfs_without_ko,
    multi_ko_orf_count = ko_audit_values$multi_ko_orf_count,
    orf_ko_association_count = ko_audit_values$orf_ko_association_count,
    multi_ko_policy = ko_audit_values$multi_ko_policy,
    ko_denominator_basis = ko_audit_values$ko_denominator_basis
  )
}

# ---- SQM data normalization and validation --------------------------------

# Clean KEGG IDs before expanding multi-KO ORFs into individual observations.
clean_ko_field <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\*", "")
  str_trim(x)
}

extract_ko_ids <- function(x) {
  if (is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(character())
  }
  unique(str_extract_all(clean_ko_field(x), "K[0-9]{5}")[[1]])
}

extract_ec_codes <- function(x) {
  match <- str_match(as.character(x), "\\[EC:([^]]+)\\]")[, 2]
  cleaned <- stringr::str_replace_all(stringr::str_squish(match), "\\s+", ";")
  cleaned[is.na(match) | !nzchar(trimws(match))] <- NA_character_
  cleaned
}

extract_ec_ids <- function(x) {
  ec_groups <- stringr::str_extract_all(as.character(x), "\\[EC:[^]]+\\]")
  lapply(ec_groups, function(groups) {
    ec_ids <- stringr::str_extract_all(
      paste(groups, collapse = " "),
      "[0-9]+\\.[0-9]+\\.[0-9]+\\.(?:[0-9]+|-)"
    )[[1L]]
    sort(unique(as.character(ec_ids)))
  })
}

build_ko_metadata <- function(sqm, ko_ids = NULL) {
  kegg_names <- sqm$misc$KEGG_names
  if (is.null(kegg_names) || is.null(names(kegg_names))) {
    stop("SQM KEGG metadata requires a named misc$KEGG_names vector.", call. = FALSE)
  }
  if (anyNA(names(kegg_names)) || any(!nzchar(names(kegg_names))) ||
      anyDuplicated(names(kegg_names)) > 0L) {
    stop("SQM misc$KEGG_names requires unique, non-empty KO names.", call. = FALSE)
  }

  if (is.null(ko_ids)) {
    ko_ids <- names(kegg_names)
  } else {
    ko_ids <- unique(as.character(ko_ids))
  }
  labels <- unname(as.character(kegg_names[ko_ids]))
  labels[is.na(labels) | !nzchar(trimws(labels))] <- ko_ids[
    is.na(labels) | !nzchar(trimws(labels))
  ]
  ec_lists <- extract_ec_ids(labels)

  tibble::tibble(
    ko_id = ko_ids,
    kegg_function = labels,
    ec_codes = vapply(
      ec_lists,
      function(ec_ids) {
        if (length(ec_ids) == 0L) NA_character_ else paste(ec_ids, collapse = ";")
      },
      character(1)
    )
  )
}

build_ko_ec_map <- function(sqm, ko_ids = NULL) {
  metadata <- build_ko_metadata(sqm, ko_ids)
  metadata |>
    mutate(ec_code = extract_ec_ids(.data$kegg_function)) |>
    tidyr::unnest_longer("ec_code") |>
    filter(!is.na(.data$ec_code) & nzchar(.data$ec_code)) |>
    select("ko_id", "ec_code") |>
    distinct()
}

normalize_taxon_value <- function(x) {
  x <- as.character(x)
  empty <- is.na(x) | !nzchar(trimws(x))
  x[empty] <- "Unclassified"
  x
}

select_top_classified_taxa <- function(data, taxon_col, value_col, top_n) {
  missing_columns <- setdiff(c(taxon_col, value_col), colnames(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Taxonomy ranking requires columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  validate_positive_integer(top_n, "top_n_taxa")

  taxon_values <- normalize_taxon_value(data[[taxon_col]])
  if (any(taxon_values == "Other")) {
    stop(
      "The source taxonomy label 'Other' is reserved for classified taxa outside Top N.",
      call. = FALSE
    )
  }

  tibble::tibble(
    taxon = taxon_values,
    value = as.numeric(data[[value_col]])
  ) |>
    filter(.data$taxon != "Unclassified", !is.na(.data$value), .data$value > 0) |>
    group_by(.data$taxon) |>
    summarise(total_value = sum(.data$value), .groups = "drop") |>
    arrange(desc(.data$total_value), .data$taxon) |>
    slice_head(n = top_n) |>
    pull(.data$taxon)
}

collapse_taxa_preserving_unclassified <- function(taxon_values, top_taxa) {
  taxon_values <- normalize_taxon_value(taxon_values)
  if (any(taxon_values == "Other")) {
    stop(
      "The source taxonomy label 'Other' is reserved for classified taxa outside Top N.",
      call. = FALSE
    )
  }

  dplyr::case_when(
    taxon_values == "Unclassified" ~ "Unclassified",
    taxon_values %in% top_taxa ~ taxon_values,
    TRUE ~ "Other"
  )
}

format_display_number <- function(x, suffix = "") {
  dplyr::case_when(
    is.na(x) ~ paste0("NA", suffix),
    x > 0 & x < 0.01 ~ paste0("<0.01", suffix),
    TRUE ~ paste0(formatC(x, format = "f", digits = 2, big.mark = ","), suffix)
  )
}

format_display_percent <- function(x) {
  format_display_number(x * 100, suffix = "%")
}

# Validate CLI numeric parameters before they affect data selection or output size.
validate_positive_integer <- function(x, arg_name) {
  if (
    length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x <= 0 ||
      x > .Machine$integer.max ||
      x != floor(x)
  ) {
    stop(arg_name, " must be a positive integer.", call. = FALSE)
  }
}

parse_positive_integer_arg <- function(value, arg_name) {
  if (
    length(value) != 1L ||
      is.na(value) ||
      !is.character(value) ||
      !grepl("^[0-9]+$", value)
  ) {
    stop(arg_name, " must be a positive integer string.", call. = FALSE)
  }

  numeric_value <- suppressWarnings(as.numeric(value))
  if (
    is.na(numeric_value) ||
      !is.finite(numeric_value) ||
      numeric_value <= 0 ||
      numeric_value > .Machine$integer.max
  ) {
    stop(arg_name, " must be between 1 and ", .Machine$integer.max, ".", call. = FALSE)
  }

  as.integer(numeric_value)
}

# Resolve user-supplied pathway codes or names to canonical SQM pathway names.
resolve_pathways <- function(sqm, requested_pathways) {
  raw_candidates <- unique(as.character(sqm$misc$KEGG_paths))
  raw_candidates <- raw_candidates[!is.na(raw_candidates) & nzchar(raw_candidates)]
  parsed_candidates <- map_dfr(raw_candidates, parse_kegg_pathway_entries)
  valid_candidates <- parsed_candidates |>
    filter(.data$is_pathway_map) |>
    distinct(
      .data$pathway_root,
      .data$pathway_category,
      .data$canonical_pathway_name
    )
  canonical_candidates <- unique(valid_candidates$canonical_pathway_name)

  resolved <- map(requested_pathways, function(requested) {
    if (grepl("^p[0-9]{5}$", requested, ignore.case = TRUE)) {
      stop(
        "Use ", substring(requested, 2), " instead of ", requested, ".",
        call. = FALSE
      )
    }

    if (grepl("^[0-9]{5}$", requested)) {
      row <- known_pathways |>
        filter(.data$pathway_id == requested)
      if (nrow(row) == 0L) {
        stop(
          "Unsupported pathway code: ", requested, ".",
          call. = FALSE
        )
      }
      return(list(
        input_value = requested,
        pathway_id = row$pathway_id[[1]],
        canonical_pathway_name = row$canonical_pathway_name[[1]]
      ))
    }

    matching_rows <- valid_candidates |>
      filter(tolower(.data$canonical_pathway_name) == tolower(requested))
    if (nrow(matching_rows) == 0L) {
      invalid_match <- parsed_candidates |>
        filter(
          tolower(.data$canonical_pathway_name) == tolower(requested),
          !.data$is_pathway_map
        )
      if (nrow(invalid_match) > 0L) {
        stop(
          "'", requested,
          "' is not a KEGG PATHWAY; it belongs to a BRITE/non-pathway hierarchy.",
          call. = FALSE
        )
      }
      candidates <- canonical_candidates[grepl(tolower(requested), tolower(canonical_candidates), fixed = TRUE)]
      candidate_text <- if (length(candidates) > 0L) {
        paste(candidates, collapse = "; ")
      } else {
        "no useful candidates"
      }
      stop(
        "Pathway not found: ", requested, ". Candidates: ", candidate_text,
        call. = FALSE
      )
    }
    if (nrow(matching_rows) > 1L) {
      hierarchy_labels <- paste(
        matching_rows$pathway_root,
        matching_rows$pathway_category,
        matching_rows$canonical_pathway_name,
        sep = "; "
      )
      stop(
        "Ambiguous pathway: ", requested, ". Matches: ",
        paste(hierarchy_labels, collapse = " | "),
        call. = FALSE
      )
    }

    canonical_name <- matching_rows$canonical_pathway_name[[1L]]
    row <- known_pathways |>
      filter(tolower(.data$canonical_pathway_name) == tolower(canonical_name))
    pathway_id <- if (nrow(row) > 0L) row$pathway_id[[1]] else NA_character_

    list(
      input_value = requested,
      pathway_id = pathway_id,
      canonical_pathway_name = canonical_name
    )
  })

  unique_by_name <- !duplicated(map_chr(resolved, "canonical_pathway_name"))
  resolved[unique_by_name]
}

parse_kegg_pathway_entries <- function(pathway_field) {
  empty_entries <- tibble::tibble(
    pathway_root = character(),
    pathway_category = character(),
    canonical_pathway_name = character(),
    is_pathway_map = logical()
  )
  if (length(pathway_field) == 0L || all(is.na(pathway_field))) {
    return(empty_entries)
  }

  pathway_values <- unlist(
    str_split(as.character(pathway_field), "\\s*\\|\\s*"),
    use.names = FALSE
  )
  pathway_values <- pathway_values[
    !is.na(pathway_values) & nzchar(trimws(pathway_values))
  ]
  if (length(pathway_values) == 0L) {
    return(empty_entries)
  }

  hierarchy_parts <- str_split(pathway_values, "\\s*;\\s*")
  hierarchy_parts <- hierarchy_parts[lengths(hierarchy_parts) == 3L]
  if (length(hierarchy_parts) == 0L) {
    return(empty_entries)
  }

  entries <- tibble::tibble(
    pathway_root = map_chr(hierarchy_parts, ~ trimws(.x[[1L]])),
    pathway_category = map_chr(hierarchy_parts, ~ trimws(.x[[2L]])),
    canonical_pathway_name = map_chr(hierarchy_parts, ~ trimws(.x[[3L]]))
  ) |>
    filter(
      nzchar(.data$pathway_root),
      nzchar(.data$pathway_category),
      nzchar(.data$canonical_pathway_name)
    ) |>
    mutate(is_pathway_map = .data$pathway_root %in% kegg_pathway_roots) |>
    distinct()

  entries
}

split_kegg_pathway_field <- function(pathway_field) {
  parse_kegg_pathway_entries(pathway_field) |>
    filter(.data$is_pathway_map) |>
    pull(.data$canonical_pathway_name) |>
    unique()
}

parse_kegg_pathway_membership <- function(orf_ids, pathway_fields) {
  if (length(orf_ids) != length(pathway_fields)) {
    stop("orf_ids and pathway_fields must have the same length.", call. = FALSE)
  }

  empty_membership <- tibble::tibble(
    orf_id = character(),
    pathway_root = character(),
    pathway_category = character(),
    canonical_pathway_name = character(),
    is_pathway_map = logical()
  )
  if (length(orf_ids) == 0L) {
    return(empty_membership)
  }

  entries_by_orf <- str_split(as.character(pathway_fields), "\\s*\\|\\s*")
  membership_index <- tibble::tibble(
    orf_id = rep(as.character(orf_ids), lengths(entries_by_orf)),
    hierarchy_entry = trimws(unlist(entries_by_orf, use.names = FALSE))
  ) |>
    filter(!is.na(.data$hierarchy_entry) & nzchar(.data$hierarchy_entry))
  if (nrow(membership_index) == 0L) {
    return(empty_membership)
  }

  unique_entries <- unique(membership_index$hierarchy_entry)
  parsed_entries <- map_dfr(unique_entries, function(hierarchy_entry) {
    parse_kegg_pathway_entries(hierarchy_entry) |>
      mutate(hierarchy_entry = hierarchy_entry, .before = 1L)
  })
  if (nrow(parsed_entries) == 0L) {
    return(empty_membership)
  }

  membership_index |>
    inner_join(
      parsed_entries,
      by = "hierarchy_entry",
      relationship = "many-to-one"
    ) |>
    select(-"hierarchy_entry") |>
    distinct()
}

pathway_id_for_name <- function(pathway_name) {
  matched_pathway <- known_pathways |>
    filter(tolower(.data$canonical_pathway_name) == tolower(pathway_name))
  if (nrow(matched_pathway) > 0L) {
    return(matched_pathway$pathway_id[[1L]])
  }
  matched_alias <- pathway_name_aliases |>
    filter(tolower(.data$alias_name) == tolower(pathway_name))
  if (nrow(matched_alias) > 0L) {
    return(matched_alias$pathway_id[[1L]])
  }
  NA_character_
}

normalize_kegg_pathway_name <- function(pathway_name) {
  pathway_name <- as.character(pathway_name)
  pathway_name <- sub(
    "\\s*-\\s*Reference pathway\\s*$",
    "",
    pathway_name,
    ignore.case = TRUE
  )
  tolower(trimws(gsub("\\s+", " ", pathway_name)))
}

parse_kegg_pathway_catalog <- function(lines) {
  lines <- as.character(lines)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  valid <- lengths(fields) >= 2L
  if (!any(valid)) {
    stop("KEGG pathway catalog is empty or malformed.", call. = FALSE)
  }
  fields <- fields[valid]
  catalog <- tibble::tibble(
    pathway_id = sub("^(?:path:)?ko", "", vapply(fields, `[[`, character(1), 1L)),
    pathway_name = vapply(fields, `[[`, character(1), 2L)
  ) |>
    filter(grepl("^[0-9]{5}$", .data$pathway_id)) |>
    mutate(normalized_name = normalize_kegg_pathway_name(.data$pathway_name)) |>
    distinct(.data$pathway_id, .data$normalized_name, .keep_all = TRUE)
  if (nrow(catalog) == 0L) {
    stop("KEGG pathway catalog contains no valid ko pathway entries.", call. = FALSE)
  }
  catalog
}

download_kegg_pathway_catalog <- function(
    endpoint = "https://rest.kegg.jp/list/pathway/ko") {
  connection <- NULL
  lines <- tryCatch({
    connection <- url(endpoint, open = "rb")
    on.exit(close(connection), add = TRUE)
    readLines(connection, warn = FALSE, encoding = "UTF-8")
  }, error = function(error) {
    stop(
      "Unable to download the KEGG pathway catalog from ", endpoint,
      ": ", conditionMessage(error),
      call. = FALSE
    )
  })
  parse_kegg_pathway_catalog(lines)
}

resolve_pathway_id_from_catalog <- function(pathway_name, catalog) {
  required_columns <- c("pathway_id", "normalized_name")
  missing_columns <- setdiff(required_columns, colnames(catalog))
  if (length(missing_columns) > 0L) {
    stop(
      "KEGG pathway catalog is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  normalized_name <- normalize_kegg_pathway_name(pathway_name)
  matches <- catalog |>
    filter(.data$normalized_name == .env$normalized_name) |>
    distinct(.data$pathway_id)
  if (nrow(matches) == 0L) {
    stop(
      "KEGG pathway name could not be resolved to an ID: ", pathway_name,
      call. = FALSE
    )
  }
  if (nrow(matches) > 1L) {
    stop(
      "KEGG pathway name maps ambiguously to IDs ",
      paste(matches$pathway_id, collapse = ", "),
      ": ", pathway_name,
      call. = FALSE
    )
  }
  matches$pathway_id[[1L]]
}

resolve_kegg_pathway_id <- local({
  catalog <- NULL

  function(pathway_name, catalog_loader = download_kegg_pathway_catalog) {
    curated_id <- pathway_id_for_name(pathway_name)
    if (!is.na(curated_id)) {
      return(curated_id)
    }
    if (is.null(catalog)) {
      catalog <<- catalog_loader()
    }
    resolve_pathway_id_from_catalog(pathway_name, catalog)
  }
})

# Rank pathways by total TPM across the selected samples for `top20` mode.
select_top_pathways <- function(sqm, selected_samples, pathway_top_n = default_pathway_top_n) {
  validate_positive_integer(pathway_top_n, "pathway_top_n")
  orf_table <- as.data.frame(sqm$orfs$table, check.names = FALSE) |>
    tibble::rownames_to_column("orf_id")
  tpm_table <- as.data.frame(sqm$orfs$tpm, check.names = FALSE) |>
    tibble::rownames_to_column("orf_id")

  if (!"KEGGPATH" %in% colnames(orf_table)) {
    stop("The ORF table does not contain the KEGGPATH column.", call. = FALSE)
  }
  if (anyDuplicated(orf_table$orf_id) > 0L || anyDuplicated(tpm_table$orf_id) > 0L) {
    stop("orf_id keys must be unique in the ORF and TPM tables.", call. = FALSE)
  }
  if (!setequal(orf_table$orf_id, tpm_table$orf_id)) {
    stop("orf_id values do not match between the ORF and TPM tables.", call. = FALSE)
  }
  validate_samples(selected_samples, colnames(tpm_table))

  pathway_membership <- parse_kegg_pathway_membership(
    orf_table$orf_id,
    orf_table$KEGGPATH
  ) |>
    filter(.data$is_pathway_map) |>
    select(-"is_pathway_map") |>
    distinct(
      .data$orf_id,
      .data$pathway_root,
      .data$pathway_category,
      .data$canonical_pathway_name
    )

  hierarchy_by_name <- pathway_membership |>
    distinct(
      .data$canonical_pathway_name,
      .data$pathway_root,
      .data$pathway_category
    ) |>
    count(.data$canonical_pathway_name, name = "hierarchy_count") |>
    filter(.data$hierarchy_count > 1L)
  if (nrow(hierarchy_by_name) > 0L) {
    stop(
      "Ambiguous KEGG pathway hierarchy for: ",
      paste(hierarchy_by_name$canonical_pathway_name, collapse = "; "),
      call. = FALSE
    )
  }

  if (nrow(pathway_membership) == 0L) {
    warning(
      "No valid three-level KEGG PATHWAY entries were found for top20 selection.",
      call. = FALSE
    )
    return(list())
  }

  pathway_totals <- tpm_table |>
    pivot_longer(
      cols = -all_of("orf_id"),
      names_to = "sample",
      values_to = "tpm"
    ) |>
    filter(.data$sample %in% selected_samples) |>
    inner_join(pathway_membership, by = "orf_id", relationship = "many-to-many") |>
    group_by(
      .data$pathway_root,
      .data$pathway_category,
      .data$canonical_pathway_name
    ) |>
    summarise(total_tpm = sum(as.numeric(.data$tpm), na.rm = TRUE), .groups = "drop") |>
    arrange(
      desc(.data$total_tpm),
      .data$canonical_pathway_name,
      .data$pathway_root,
      .data$pathway_category
    ) |>
    slice_head(n = pathway_top_n)

  map(seq_len(nrow(pathway_totals)), function(index) {
    pathway_name <- pathway_totals$canonical_pathway_name[[index]]
    list(
      input_value = pathway_name,
      pathway_id = pathway_id_for_name(pathway_name),
      pathway_root = pathway_totals$pathway_root[[index]],
      pathway_category = pathway_totals$pathway_category[[index]],
      canonical_pathway_name = pathway_name,
      total_tpm = pathway_totals$total_tpm[[index]]
    )
  })
}

select_context_pathway_groups <- function(
    context_sqm,
    resolved_defined_pathways,
    pathway_selection_modes,
    selected_samples,
    pathway_top_n = default_pathway_top_n) {
  pathway_selection_modes <- normalize_pathway_selection_modes(pathway_selection_modes)
  groups <- list()

  if ("defined" %in% pathway_selection_modes) {
    groups[["defined"]] <- resolved_defined_pathways
  }
  if ("top20" %in% pathway_selection_modes) {
    context_top_pathways <- select_top_pathways(
      context_sqm,
      selected_samples,
      pathway_top_n
    )
    if (length(context_top_pathways) > 0L) {
      groups[["top20"]] <- context_top_pathways
    }
  }

  groups
}

select_pathway_groups <- function(
    sqm,
    requested_pathways,
    pathway_selection_modes,
    selected_samples,
    pathway_top_n = default_pathway_top_n) {
  pathway_selection_modes <- normalize_pathway_selection_modes(pathway_selection_modes)
  resolved_defined_pathways <- list()
  if ("defined" %in% pathway_selection_modes) {
    defined_pathways <- if (length(requested_pathways) == 0L) {
      default_pathway_ids
    } else {
      requested_pathways
    }
    resolved_defined_pathways <- resolve_pathways(sqm, defined_pathways)
  }

  select_context_pathway_groups(
    context_sqm = sqm,
    resolved_defined_pathways = resolved_defined_pathways,
    pathway_selection_modes = pathway_selection_modes,
    selected_samples = selected_samples,
    pathway_top_n = pathway_top_n
  )
}

normalize_sqm_project_dir <- function(project_dir) {
  normalizePath(project_dir, winslash = "/", mustWork = TRUE)
}

load_sqm_project <- function(project_dir, tax_mode) {
  SQMtools::loadSQM(
    project_path = normalize_sqm_project_dir(project_dir),
    tax_mode = tax_mode,
    trusted_functions_only = FALSE,
    load_sequences = FALSE
  )
}

# Ensure the SQM object exposes the ORF components required by every analysis mode.
validate_sqm_object <- function(sqm) {
  required_orf_parts <- c("table", "tax", "tpm")
  missing_parts <- setdiff(required_orf_parts, names(sqm$orfs))
  if (length(missing_parts) > 0L) {
    stop(
      "Missing components in sqm$orfs: ",
      paste(missing_parts, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(rownames(sqm$orfs$table)) ||
      is.null(rownames(sqm$orfs$tax)) ||
      is.null(rownames(sqm$orfs$tpm))) {
    stop("ORF tables must have row names that can be used as orf_id values.", call. = FALSE)
  }
}

validate_tax_mode <- function(tax_mode) {
  allowed_tax_modes <- c("prokfilter", "allfilter", "nofilter")
  if (
    length(tax_mode) != 1L || is.na(tax_mode) ||
      !tax_mode %in% allowed_tax_modes
  ) {
    stop(
      "tax_mode must be one of: ",
      paste(allowed_tax_modes, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  tax_mode
}

validate_tpm_matrix <- function(tpm_matrix, label) {
  tpm_frame <- as.data.frame(tpm_matrix, check.names = FALSE)
  numeric_columns <- vapply(tpm_frame, is.numeric, logical(1))
  if (ncol(tpm_frame) == 0L || !all(numeric_columns)) {
    stop(label, " must contain only numeric sample columns.", call. = FALSE)
  }
  values <- unlist(tpm_frame, use.names = FALSE)
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
    stop(label, " values must be finite and non-negative.", call. = FALSE)
  }
  TRUE
}

validate_sqm_tpm_inputs <- function(
    sqm,
    selected_samples,
    require_kegg_tpm = FALSE) {
  validate_samples(selected_samples, colnames(sqm$orfs$tpm))
  validate_tpm_matrix(
    as.data.frame(sqm$orfs$tpm, check.names = FALSE)[selected_samples],
    "sqm$orfs$tpm"
  )

  if (isTRUE(require_kegg_tpm)) {
    kegg_tpm <- sqm$functions$KEGG$tpm
    if (is.null(kegg_tpm)) {
      stop("sqm$functions$KEGG$tpm is required for KEGG plot data.", call. = FALSE)
    }
    validate_samples(selected_samples, colnames(kegg_tpm))
    validate_tpm_matrix(
      as.data.frame(kegg_tpm, check.names = FALSE)[selected_samples],
      "sqm$functions$KEGG$tpm"
    )
  }
  TRUE
}

validate_samples <- function(requested_samples, available_samples) {
  if (anyDuplicated(requested_samples) > 0L) {
    stop("Samples must be unique and retain one explicit display position each.", call. = FALSE)
  }
  missing_samples <- setdiff(requested_samples, available_samples)
  if (length(missing_samples) > 0L) {
    stop(
      "Samples not found: ", paste(missing_samples, collapse = ", "),
      ". Available: ", paste(available_samples, collapse = ", "),
      call. = FALSE
    )
  }
}

resolve_sample_selection <- function(samples_argument, available_samples) {
  selected_samples <- if (is.null(samples_argument)) {
    as.character(available_samples)
  } else {
    split_csv_arg(samples_argument)
  }
  if (length(selected_samples) == 0L) {
    stop("samples cannot be empty.", call. = FALSE)
  }
  validate_samples(selected_samples, available_samples)
  list(
    samples = selected_samples,
    basis = if (is.null(samples_argument)) "sqm_column_order" else "cli"
  )
}

validate_taxonomy_ranks <- function(requested_ranks, available_ranks) {
  missing_ranks <- setdiff(requested_ranks, available_ranks)
  if (length(missing_ranks) > 0L) {
    stop(
      "Taxonomic ranks not found: ", paste(missing_ranks, collapse = ", "),
      call. = FALSE
    )
  }
}

subset_pathway <- function(
    sqm,
    canonical_pathway,
    subset_fun = SQMtools::subsetFun) {
  subset_fun(
    SQM = sqm,
    fun = canonical_pathway,
    columns = "KEGGPATH",
    ignore_case = FALSE,
    fixed = TRUE,
    allow_empty = TRUE
  )
}

is_empty_pathway_subset <- function(pathway_sqm) {
  if (is.null(pathway_sqm)) {
    return(TRUE)
  }
  if (is.null(pathway_sqm$orfs)) {
    return(TRUE)
  }
  orf_table <- pathway_sqm$orfs$table
  orf_tpm <- pathway_sqm$orfs$tpm
  if (is.null(orf_table) && is.null(orf_tpm)) {
    return(TRUE)
  }
  if (is.null(orf_table) || is.null(orf_tpm)) {
    stop("Pathway subset is missing ORF table or TPM data.", call. = FALSE)
  }
  nrow(as.data.frame(orf_table, check.names = FALSE)) == 0L ||
    nrow(as.data.frame(orf_tpm, check.names = FALSE)) == 0L
}

prepare_context_pathway_subsets <- function(
    context_sqm,
    pathway_entries,
    context_label,
    subset_fun = SQMtools::subsetFun,
    include_kegg_oracle = FALSE,
    pathway_ko_resolver = resolve_pathway_ko_ids,
    pathway_id_resolver = resolve_kegg_pathway_id,
    selected_samples = colnames(context_sqm$orfs$tpm)) {
  pathway_sqms <- list()
  skips <- tibble::tibble(
    context = character(),
    pathway = character(),
    reason = character()
  )

  for (pathway_info in pathway_entries) {
    if (is.na(pathway_info$pathway_id) ||
        !grepl("^[0-9]{5}$", pathway_info$pathway_id)) {
      pathway_info$pathway_id <- pathway_id_resolver(pathway_info$pathway_name)
    }
    pathway_ko_ids <- NULL
    if (isTRUE(include_kegg_oracle)) {
      if (!pathview_is_exportable(
          pathway_info$pathway_selection,
          pathway_info$pathway_id
      )) {
        stop(
          "KEGG pathway analysis requires a resolvable pathway ID for '",
          pathway_info$pathway_name,
          "'.",
          call. = FALSE
        )
      }
      pathway_ko_ids <- pathway_ko_resolver(pathway_info$pathway_id)
      if (length(pathway_ko_ids) == 0L) {
        warning(
          "No ortholog nodes in KEGG pathway ",
          pathway_info$pathway_id,
          " (",
          pathway_info$pathway_name,
          ").",
          call. = FALSE
        )
        skips <- bind_rows(
          skips,
          tibble::tibble(
            context = as.character(context_label),
            pathway = pathway_info$pathway_name,
            reason = "no_ortholog_nodes"
          )
        )
        next
      }
    }

    pathway_sqm <- if (isTRUE(include_kegg_oracle)) {
      subset_orfs_by_ko_membership(context_sqm, pathway_ko_ids)
    } else {
      subset_pathway(
        context_sqm,
        pathway_info$pathway_name,
        subset_fun = subset_fun
      )
    }
    if (is_empty_pathway_subset(pathway_sqm)) {
      if (isTRUE(include_kegg_oracle)) {
        kegg_tpm <- as.data.frame(
          context_sqm$functions$KEGG$tpm,
          check.names = FALSE
        )
        available_ko_ids <- intersect(pathway_ko_ids, rownames(kegg_tpm))
        if (length(available_ko_ids) > 0L &&
            any(as.matrix(kegg_tpm[
              available_ko_ids,
              selected_samples,
              drop = FALSE
            ]) > 0)) {
          stop(
            "Cannot allocate positive official SQM KEGG TPM without ORFs for pathway ",
            pathway_info$pathway_id, " (", pathway_info$pathway_name, ").",
            call. = FALSE
          )
        }
      }
      warning(
        "Skipping empty context x pathway combination: context=", context_label,
        " | pathway=", pathway_info$pathway_name,
        call. = FALSE
      )
      skips <- bind_rows(
        skips,
        tibble::tibble(
          context = as.character(context_label),
          pathway = pathway_info$pathway_name,
          reason = if (isTRUE(include_kegg_oracle)) "empty_ko_subset" else "empty_subset"
        )
      )
      next
    }
    pathway_key <- paste(
      pathway_info$pathway_selection,
      pathway_info$pathway_name,
      sep = "::"
    )
    pathway_sqms[[pathway_key]] <- list(
      pathway_name = pathway_info$pathway_name,
      pathway_id = pathway_info$pathway_id,
      pathway_selection = pathway_info$pathway_selection,
      pathway_sqm = pathway_sqm,
      context_sqm = context_sqm,
      pathway_ko_ids = pathway_ko_ids
    )
  }

  list(pathway_sqms = pathway_sqms, skips = skips)
}

resolve_taxa_filters <- function(sqm, requested_taxa) {
  tax_table <- as.data.frame(sqm$orfs$tax, check.names = FALSE)
  tax_columns <- intersect(all_taxonomy_columns, colnames(tax_table))

  if (is.null(rownames(tax_table)) || anyDuplicated(rownames(tax_table)) > 0L) {
    stop("ORF taxonomy must have unique row names for taxon filtering.", call. = FALSE)
  }

  map(requested_taxa, function(requested_taxon) {
    taxon_lower <- tolower(trimws(requested_taxon))
    matched_ranks <- tax_columns[vapply(
      tax_columns,
      function(rank) {
        tax_values <- tolower(trimws(as.character(tax_table[[rank]])))
        any(tax_values == taxon_lower, na.rm = TRUE)
      },
      logical(1)
    )]

    if (length(matched_ranks) == 0L) {
      stop(
        "Taxon not found in the taxonomy table: ", requested_taxon,
        call. = FALSE
      )
    }

    if (length(matched_ranks) > 1L) {
      stop(
        "Ambiguous taxon: ", requested_taxon, ". Found in multiple columns: ",
        paste(matched_ranks, collapse = ", "),
        call. = FALSE
      )
    }

    matched_rank <- matched_ranks[[1]]
    rank_values <- tolower(trimws(as.character(tax_table[[matched_rank]])))
    orf_ids <- rownames(tax_table)[!is.na(rank_values) & rank_values == taxon_lower]

    list(
      taxon = requested_taxon,
      rank = matched_rank,
      orf_ids = orf_ids
    )
  })
}

subset_sqm_by_taxon <- function(
  sqm,
  orf_ids,
  subset_orfs_fn = SQMtools::subsetORFs
) {
  if (length(orf_ids) == 0L || anyNA(orf_ids) || anyDuplicated(orf_ids) > 0L) {
    stop("Taxon filtering requires a non-empty set of unique ORF IDs.", call. = FALSE)
  }

  subset_sqm <- subset_orfs_fn(
    SQM = sqm,
    orfs = orf_ids,
    tax_source = "orfs",
    trusted_functions_only = FALSE,
    ignore_unclassified_functions = FALSE,
    rescale_tpm = FALSE,
    rescale_copy_number = FALSE,
    recalculate_bin_stats = FALSE,
    contigs_override = NULL,
    allow_empty = FALSE
  )

  actual_orf_ids <- rownames(subset_sqm$orfs$table)
  if (!identical(actual_orf_ids, orf_ids)) {
    stop(
      "ORF subset postcondition failed: requested and returned ORF IDs differ.",
      call. = FALSE
    )
  }

  subset_sqm
}

select_orf_ids_by_ko <- function(sqm, ko_ids) {
  ko_ids <- sort(unique(as.character(ko_ids)))
  if (length(ko_ids) == 0L) {
    return(character())
  }
  if (anyNA(ko_ids) || any(!grepl("^K[0-9]{5}$", ko_ids))) {
    stop("Pathway KO IDs must be valid KEGG ortholog IDs.", call. = FALSE)
  }

  orf_table <- as.data.frame(sqm$orfs$table, check.names = FALSE)
  if (!"KEGG ID" %in% colnames(orf_table)) {
    stop("Pathway membership requires the ORF annotation column: KEGG ID", call. = FALSE)
  }
  if (is.null(rownames(orf_table)) || anyDuplicated(rownames(orf_table)) > 0L) {
    stop("Pathway membership requires unique ORF row names.", call. = FALSE)
  }

  orf_ko_ids <- lapply(as.character(orf_table[["KEGG ID"]]), extract_ko_ids)
  keep <- vapply(
    orf_ko_ids,
    function(orf_kos) any(orf_kos %in% ko_ids),
    logical(1)
  )
  rownames(orf_table)[keep]
}

subset_sqm_by_orf_ids <- function(
    sqm,
    orf_ids,
    subset_orfs_fn = SQMtools::subsetORFs) {
  orf_ids <- unique(as.character(orf_ids))
  if (length(orf_ids) == 0L || anyNA(orf_ids)) {
    stop("ORF subsetting requires a non-empty set of ORF IDs.", call. = FALSE)
  }

  available_ids <- rownames(sqm$orfs$table)
  missing_ids <- setdiff(orf_ids, available_ids)
  if (length(missing_ids) > 0L) {
    stop(
      "ORF subset contains IDs absent from sqm$orfs$table: ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }

  if (inherits(sqm, "SQM")) {
    subset_sqm <- subset_orfs_fn(
      SQM = sqm,
      orfs = orf_ids,
      tax_source = "orfs",
      trusted_functions_only = FALSE,
      ignore_unclassified_functions = FALSE,
      rescale_tpm = FALSE,
      rescale_copy_number = FALSE,
      recalculate_bin_stats = FALSE,
      contigs_override = NULL,
      allow_empty = FALSE
    )
  } else {
    subset_sqm <- sqm
    for (component_name in names(subset_sqm$orfs)) {
      component <- subset_sqm$orfs[[component_name]]
      component_ids <- rownames(component)
      if (!is.null(component_ids) && all(orf_ids %in% component_ids)) {
        subset_sqm$orfs[[component_name]] <- component[orf_ids, , drop = FALSE]
      }
    }
  }

  if (!identical(rownames(subset_sqm$orfs$table), orf_ids)) {
    stop("ORF subset postcondition failed: requested and returned IDs differ.", call. = FALSE)
  }
  subset_sqm
}

build_ko_expansion_audit <- function(orf_table) {
  if (!"KEGG ID" %in% colnames(orf_table)) {
    stop("KO expansion audit requires column: KEGG ID", call. = FALSE)
  }
  ko_ids <- map(as.character(orf_table[["KEGG ID"]]), extract_ko_ids)
  ko_counts <- lengths(ko_ids)
  tibble::tibble(
    input_orf_count = as.integer(nrow(orf_table)),
    excluded_orfs_without_ko = as.integer(sum(ko_counts == 0L)),
    multi_ko_orf_count = as.integer(sum(ko_counts > 1L)),
    orf_ko_association_count = as.integer(sum(ko_counts)),
    multi_ko_policy = "full_tpm_per_ko",
    ko_denominator_basis = "expanded_orf_sample_ko_tpm"
  )
}

# Expand the raw ORF × sample × KO associations before canonical allocation.
# build_pathway_ko_result() divides each ORF TPM by its complete KO count,
# filters to KGML membership, then aligns each KO to the official SQM margin.
build_orf_long_result <- function(pathway_sqm, selected_samples) {
  orf_table <- as.data.frame(pathway_sqm$orfs$table, check.names = FALSE) |>
    tibble::rownames_to_column("orf_id")
  tax_table <- as.data.frame(pathway_sqm$orfs$tax, check.names = FALSE) |>
    tibble::rownames_to_column("orf_id")
  tpm_table <- as.data.frame(pathway_sqm$orfs$tpm, check.names = FALSE) |>
    tibble::rownames_to_column("orf_id")

  if (anyDuplicated(orf_table$orf_id) > 0L ||
      anyDuplicated(tax_table$orf_id) > 0L ||
      anyDuplicated(tpm_table$orf_id) > 0L) {
    stop("orf_id keys must be unique in all source tables.", call. = FALSE)
  }

  ids_table <- orf_table$orf_id
  ids_tax <- tax_table$orf_id
  ids_tpm <- tpm_table$orf_id
  if (!setequal(ids_table, ids_tax) || !setequal(ids_table, ids_tpm)) {
    stop("orf_id values do not match between the table, tax, and tpm components.", call. = FALSE)
  }

  missing_samples <- setdiff(selected_samples, colnames(tpm_table))
  if (length(missing_samples) > 0L) {
    stop(
      "Samples missing from sqm$orfs$tpm: ", paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }
  validate_tpm_matrix(tpm_table[selected_samples], "sqm$orfs$tpm")

  fun_lookup <- sqm_misc_names <- pathway_sqm$misc$KEGG_names

  if (!"KEGG ID" %in% colnames(orf_table)) {
    stop("ORF KO allocation requires column: KEGG ID", call. = FALSE)
  }
  keggfun_values <- if ("KEGGFUN" %in% colnames(orf_table)) {
    as.character(orf_table[["KEGGFUN"]])
  } else {
    rep(NA_character_, nrow(orf_table))
  }
  keggpath_values <- if ("KEGGPATH" %in% colnames(orf_table)) {
    as.character(orf_table[["KEGGPATH"]])
  } else {
    rep(NA_character_, nrow(orf_table))
  }
  annotations <- tibble::tibble(
    orf_id = orf_table$orf_id,
    `KEGG ID` = as.character(orf_table[["KEGG ID"]]),
    kegg_function_raw = keggfun_values,
    ec_codes = extract_ec_codes(keggfun_values),
    KEGGPATH = keggpath_values,
    ko_ids = map(as.character(orf_table[["KEGG ID"]]), extract_ko_ids)
  )

  ko_audit <- build_ko_expansion_audit(orf_table)

  tpm_long <- tpm_table |>
    select(all_of(c("orf_id", selected_samples))) |>
    pivot_longer(
      cols = -all_of("orf_id"),
      names_to = "sample",
      values_to = "tpm"
    ) |>
    mutate(tpm = as.numeric(.data$tpm)) |>
    filter(!is.na(.data$tpm) & .data$tpm > 0)

  tax_subset <- tax_table |>
    select(any_of(c("orf_id", all_taxonomy_columns)))

  for (rank in setdiff(all_taxonomy_columns, names(tax_subset))) {
    tax_subset[[rank]] <- "Unclassified"
  }

  tax_subset <- tax_subset |>
    select(all_of(c("orf_id", all_taxonomy_columns))) |>
    mutate(across(all_of(all_taxonomy_columns), normalize_taxon_value))

  joined <- tpm_long |>
    left_join(annotations, by = "orf_id") |>
    left_join(tax_subset, by = "orf_id")

  if (any(is.na(joined$KEGGPATH))) {
    joined$KEGGPATH[is.na(joined$KEGGPATH)] <- NA_character_
  }

  expanded <- joined |>
    filter(lengths(.data$ko_ids) > 0L) |>
    tidyr::unnest_longer(ko_ids, values_to = "ko_id") |>
    mutate(
      ko_id = as.character(.data$ko_id),
      kegg_function = if_else(
        is.na(.data$kegg_function_raw) | !nzchar(trimws(.data$kegg_function_raw)),
        unname(as.character(fun_lookup[.data$ko_id])),
        .data$kegg_function_raw
      ),
      kegg_function = if_else(
        is.na(.data$kegg_function) | !nzchar(trimws(.data$kegg_function)),
        .data$ko_id,
        .data$kegg_function
      )
    ) |>
    select(
      "orf_id", "sample", "tpm", "ko_id", "kegg_function",
      "ec_codes", "KEGGPATH", all_of(all_taxonomy_columns)
    )

  list(data = expanded, audit = ko_audit)
}

build_orf_long_table <- function(pathway_sqm, selected_samples) {
  build_orf_long_result(pathway_sqm, selected_samples)$data
}

# ---- Functional, flow, and taxonomy plot data ----------------------------

# Use SQM KEGG names only as a fallback when the ORF annotation is absent.
# Match the ortholog membership used by pathview::node.map(). Compound nodes
# and other KGML node types are deliberately excluded from FLOW.
extract_pathway_ko_ids <- function(node_data) {
  if (!is.list(node_data) ||
      is.null(node_data$kegg.names) ||
      is.null(node_data$type)) {
    stop("Pathview node data must contain kegg.names and type.", call. = FALSE)
  }

  node_ids <- names(node_data$type)
  node_types <- as.character(node_data$type)
  if (is.null(node_ids) || anyNA(node_ids) || anyDuplicated(node_ids) > 0L) {
    stop("Pathview node types must have unique node IDs.", call. = FALSE)
  }

  ortholog_ids <- node_ids[node_types == "ortholog"]
  ortholog_names <- node_data$kegg.names[ortholog_ids]
  ko_ids <- unlist(
    lapply(ortholog_names, function(names_for_node) {
      stringr::str_extract_all(
        paste(as.character(names_for_node), collapse = ";"),
        "K[0-9]{5}"
      )[[1L]]
    }),
    use.names = FALSE
  )

  sort(unique(as.character(ko_ids)))
}

download_pathway_node_data <- function(
    pathway_id,
    download_fun = pathview::download.kegg,
    node_info_fun = pathview::node.info) {
  pathway_id <- as.character(pathway_id)
  if (length(pathway_id) != 1L || is.na(pathway_id) ||
      !grepl("^[0-9]{5}$", pathway_id)) {
    stop("KEGG pathway mapping requires a five-digit pathway ID.", call. = FALSE)
  }

  download_dir <- tempfile(paste0("sqm_flow_ko", pathway_id, "_"))
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(download_dir, recursive = TRUE, force = TRUE), add = TRUE)

  download_fun(
    pathway.id = pathway_id,
    species = "ko",
    kegg.dir = download_dir
  )
  xml_file <- file.path(download_dir, paste0("ko", pathway_id, ".xml"))
  if (!file.exists(xml_file)) {
    stop(
      "KEGG KGML download did not create expected file for pathway ",
      pathway_id,
      ".",
      call. = FALSE
    )
  }

  node_info_fun(xml_file)
}

resolve_pathway_ko_ids <- local({
  cache <- new.env(parent = emptyenv())

  function(pathway_id, node_data_loader = download_pathway_node_data) {
    pathway_id <- as.character(pathway_id)
    if (length(pathway_id) != 1L || is.na(pathway_id) ||
        !grepl("^[0-9]{5}$", pathway_id)) {
      stop("KEGG pathway mapping requires a five-digit pathway ID.", call. = FALSE)
    }

    if (!exists(pathway_id, envir = cache, inherits = FALSE)) {
      ko_ids <- extract_pathway_ko_ids(node_data_loader(pathway_id))
      assign(pathway_id, ko_ids, envir = cache)
    }
    get(pathway_id, envir = cache, inherits = FALSE)
  }
})

subset_orfs_by_ko_membership <- function(context_sqm, pathway_ko_ids) {
  orf_table <- as.data.frame(context_sqm$orfs$table, check.names = FALSE)
  if (!"KEGG ID" %in% colnames(orf_table)) {
    stop("KO membership requires the ORF annotation column: KEGG ID", call. = FALSE)
  }

  orf_ko_ids <- lapply(as.character(orf_table[["KEGG ID"]]), extract_ko_ids)
  keep <- vapply(
    orf_ko_ids,
    function(ko_ids) any(ko_ids %in% pathway_ko_ids),
    logical(1)
  )
  selected_orf_ids <- rownames(orf_table)[keep]

  pathway_sqm <- context_sqm
  pathway_sqm$orfs$table <- context_sqm$orfs$table[
    selected_orf_ids,
    ,
    drop = FALSE
  ]
  pathway_sqm$orfs$tax <- context_sqm$orfs$tax[
    selected_orf_ids,
    ,
    drop = FALSE
  ]
  pathway_sqm$orfs$tpm <- context_sqm$orfs$tpm[
    selected_orf_ids,
    ,
    drop = FALSE
  ]
  pathway_sqm
}

# Reproduce SQMtools KEGG aggregation for taxonomic allocation while keeping
# SQM$functions$KEGG$tpm as the authoritative functional margin. SQMtools
# divides an ORF equally among all of its KO annotations before aggregation.
build_pathway_ko_result <- function(
    context_sqm,
    selected_samples,
    pathway_ko_ids,
    tolerance = 1e-8) {
  pathway_ko_ids <- sort(unique(as.character(pathway_ko_ids)))
  if (anyNA(pathway_ko_ids) ||
      any(!grepl("^K[0-9]{5}$", pathway_ko_ids))) {
    stop("Pathway KO IDs must be valid KEGG ortholog IDs.", call. = FALSE)
  }

  kegg_tpm <- context_sqm$functions$KEGG$tpm
  if (is.null(kegg_tpm)) {
    stop("SQMtools oracle requires SQM$functions$KEGG$tpm.", call. = FALSE)
  }
  kegg_tpm <- as.data.frame(kegg_tpm, check.names = FALSE)
  if (is.null(rownames(kegg_tpm)) || anyDuplicated(rownames(kegg_tpm)) > 0L ||
      any(!nzchar(rownames(kegg_tpm)))) {
    stop("SQM$functions$KEGG$tpm requires unique, non-empty KO row names.", call. = FALSE)
  }
  validate_samples(selected_samples, colnames(kegg_tpm))
  validate_tpm_matrix(kegg_tpm[selected_samples], "SQM$functions$KEGG$tpm")

  available_ko_ids <- intersect(pathway_ko_ids, rownames(kegg_tpm))
  pathway_orf_ids <- select_orf_ids_by_ko(context_sqm, pathway_ko_ids)
  metadata <- build_ko_metadata(context_sqm, pathway_ko_ids)
  official <- tidyr::expand_grid(
    ko_id = pathway_ko_ids,
    sample = as.character(selected_samples)
  ) |>
    dplyr::left_join(
      kegg_tpm[available_ko_ids, selected_samples, drop = FALSE] |>
        tibble::rownames_to_column("ko_id") |>
        tidyr::pivot_longer(
          cols = -all_of("ko_id"),
          names_to = "sample",
          values_to = "official_tpm"
        ) |>
        dplyr::mutate(
          ko_id = as.character(.data$ko_id),
          sample = as.character(.data$sample),
          official_tpm = as.numeric(.data$official_tpm)
        ),
      by = c("ko_id", "sample"),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(official_tpm = tidyr::replace_na(.data$official_tpm, 0))

  if (length(pathway_orf_ids) == 0L) {
    if (any(official$official_tpm > tolerance)) {
      missing_keys <- paste0(
        official$sample[official$official_tpm > tolerance],
        "/",
        official$ko_id[official$official_tpm > tolerance],
        collapse = ", "
      )
      stop(
        "Cannot allocate positive official SQM KEGG TPM without ORFs: ",
        missing_keys,
        call. = FALSE
      )
    }
    empty_sqm <- subset_orfs_by_ko_membership(context_sqm, character())
    full_result <- build_orf_long_result(empty_sqm, selected_samples)
    return(list(
      data = full_result$data[0, , drop = FALSE],
      totals = official |>
        dplyr::rename(tpm = "official_tpm") |>
        dplyr::left_join(metadata, by = "ko_id", relationship = "many-to-one"),
      metadata = metadata,
      orf_ids = character(),
      pathway_sqm = empty_sqm,
      audit = dplyr::mutate(
          full_result$audit,
          multi_ko_policy = "split_tpm_equally_per_ko",
          ko_denominator_basis = "sqm_functions_kegg_tpm",
          membership_basis = "kgml_ortholog_nodes",
          pathway_ko_count = as.integer(length(pathway_ko_ids)),
          pathway_orf_count = 0L,
          raw_conservation_max_abs_error = 0,
          official_margin_max_abs_error = 0
      )
    ))
  }

  pathway_sqm <- subset_sqm_by_orf_ids(context_sqm, pathway_orf_ids)
  full_result <- build_orf_long_result(pathway_sqm, selected_samples)
  raw_allocated <- full_result$data |>
    dplyr::group_by(.data$orf_id, .data$sample) |>
    dplyr::mutate(
      ko_count = dplyr::n(),
      allocated_tpm = .data$tpm / .data$ko_count
    ) |>
    dplyr::ungroup()

  raw_conservation <- raw_allocated |>
    dplyr::group_by(.data$orf_id, .data$sample) |>
    dplyr::summarise(
      original_tpm = dplyr::first(.data$tpm),
      allocated_tpm = sum(.data$allocated_tpm),
      .groups = "drop"
    )
  if (nrow(raw_conservation) > 0L && any(
    abs(raw_conservation$allocated_tpm - raw_conservation$original_tpm) > tolerance
  )) {
    stop("Raw multi-KO allocation failed to conserve ORF TPM.", call. = FALSE)
  }

  allocated <- raw_allocated |>
    dplyr::filter(.data$ko_id %in% pathway_ko_ids)

  allocated_totals <- allocated |>
    dplyr::group_by(.data$sample, .data$ko_id) |>
    dplyr::summarise(
      allocated_total = sum(.data$allocated_tpm),
      .groups = "drop"
    )
  missing_allocation <- official |>
    dplyr::filter(.data$official_tpm > tolerance) |>
    dplyr::anti_join(allocated_totals, by = c("sample", "ko_id"))
  if (nrow(missing_allocation) > 0L) {
    missing_keys <- paste0(
      missing_allocation$sample,
      "/",
      missing_allocation$ko_id,
      collapse = ", "
    )
    stop(
      "Cannot allocate positive official SQM KEGG TPM without ORFs: ",
      missing_keys,
      call. = FALSE
    )
  }

  allocated <- allocated |>
    dplyr::left_join(
      official,
      by = c("sample", "ko_id"),
      relationship = "many-to-one"
    ) |>
    dplyr::group_by(.data$sample, .data$ko_id) |>
    dplyr::mutate(
      allocated_total = sum(.data$allocated_tpm),
      tpm = dplyr::if_else(
        .data$allocated_total > 0,
        .data$allocated_tpm * .data$official_tpm / .data$allocated_total,
        0
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$tpm > 0) |>
    dplyr::select(-all_of(c("kegg_function", "ec_codes"))) |>
    dplyr::left_join(metadata, by = "ko_id", relationship = "many-to-one") |>
    dplyr::select(
      -all_of(c(
        "ko_count", "allocated_tpm", "official_tpm", "allocated_total"
      ))
    )

  observed <- allocated |>
    dplyr::group_by(.data$sample, .data$ko_id) |>
    dplyr::summarise(flow_tpm = sum(.data$tpm), .groups = "drop") |>
    dplyr::left_join(
      official,
      by = c("sample", "ko_id"),
      relationship = "one-to-one"
  )
  if (nrow(observed) > 0L &&
      any(abs(observed$flow_tpm - observed$official_tpm) > tolerance)) {
    stop("KO allocation failed to conserve the official SQM KEGG TPM margin.", call. = FALSE)
  }

  list(
    data = allocated,
    totals = official |>
      dplyr::rename(tpm = "official_tpm") |>
      dplyr::left_join(metadata, by = "ko_id", relationship = "many-to-one"),
    metadata = metadata,
    orf_ids = pathway_orf_ids,
    pathway_sqm = pathway_sqm,
    audit = dplyr::mutate(
        full_result$audit,
        multi_ko_policy = "split_tpm_equally_per_ko",
        ko_denominator_basis = "sqm_functions_kegg_tpm",
        membership_basis = "kgml_ortholog_nodes",
        pathway_ko_count = as.integer(length(pathway_ko_ids)),
        pathway_orf_count = as.integer(length(pathway_orf_ids)),
        raw_conservation_max_abs_error = if (nrow(raw_conservation) == 0L) {
          0
        } else {
          max(abs(raw_conservation$allocated_tpm - raw_conservation$original_tpm))
        },
        official_margin_max_abs_error = if (nrow(observed) == 0L) {
          0
        } else {
          max(abs(observed$flow_tpm - observed$official_tpm))
        }
      )
  )
}

get_ko_name_lookup <- function(pathway_sqm) {
  kegg_names <- pathway_sqm$misc$KEGG_names
  if (is.null(kegg_names)) {
    return(setNames(character(), character()))
  }
  stats::setNames(as.character(kegg_names), names(kegg_names))
}

validate_pathway_analysis <- function(pathway_analysis) {
  required_fields <- c(
    "pathway_name",
    "pathway_id",
    "pathway_selection",
    "pathway_sqm",
    "selected_samples",
    "orf_long_result",
    "ko_lookup"
  )
  if (!inherits(pathway_analysis, "sqm_pathway_analysis")) {
    stop("pathway_analysis must be created by build_pathway_analysis().", call. = FALSE)
  }
  missing_fields <- setdiff(required_fields, names(pathway_analysis))
  if (length(missing_fields) > 0L) {
    stop(
      "pathway_analysis is missing fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  if (
    !is.list(pathway_analysis$orf_long_result) ||
      !all(c("data", "audit") %in% names(pathway_analysis$orf_long_result))
  ) {
    stop("pathway_analysis has an invalid ORF-long result.", call. = FALSE)
  }
  invisible(pathway_analysis)
}

build_pathway_analysis <- function(pathway_info, selected_samples) {
  required_fields <- c(
    "pathway_name", "pathway_id", "pathway_selection", "pathway_sqm"
  )
  if (!is.list(pathway_info)) {
    stop("pathway_info must be a list.", call. = FALSE)
  }
  missing_fields <- setdiff(required_fields, names(pathway_info))
  if (length(missing_fields) > 0L) {
    stop(
      "pathway_info is missing fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  selected_samples <- as.character(selected_samples)
  if (length(selected_samples) == 0L || anyNA(selected_samples)) {
    stop("pathway_analysis requires at least one selected sample.", call. = FALSE)
  }

  context_sqm <- pathway_info$context_sqm
  if (is.null(context_sqm)) {
    context_sqm <- pathway_info$pathway_sqm
  }
  pathway_ko_ids <- pathway_info$pathway_ko_ids
  if (!is.null(pathway_ko_ids)) {
    orf_long_result <- build_pathway_ko_result(
      context_sqm = context_sqm,
      selected_samples = selected_samples,
      pathway_ko_ids = pathway_ko_ids
    )
    pathway_sqm <- orf_long_result$pathway_sqm
  } else {
    # Compatibility seam for callers that construct analyses without KGML.
    # Productive pathway modes always populate pathway_ko_ids in preflight.
    pathway_sqm <- pathway_info$pathway_sqm
    orf_long_result <- build_orf_long_result(pathway_sqm, selected_samples)
    orf_long_result$totals <- orf_long_result$data |>
      group_by(.data$sample, .data$ko_id) |>
      summarise(
        tpm = sum(.data$tpm),
        kegg_function = dplyr::first(.data$kegg_function),
        ec_codes = dplyr::first(.data$ec_codes),
        .groups = "drop"
      )
    orf_long_result$metadata <- orf_long_result$totals |>
      select("ko_id", "kegg_function", "ec_codes") |>
      distinct()
    orf_long_result$orf_ids <- unique(as.character(orf_long_result$data$orf_id))
  }

  pathway_analysis <- structure(
    list(
      pathway_name = as.character(pathway_info$pathway_name),
      pathway_id = as.character(pathway_info$pathway_id),
      pathway_selection = as.character(pathway_info$pathway_selection),
      pathway_sqm = pathway_sqm,
      context_sqm = context_sqm,
      pathway_ko_ids = pathway_ko_ids,
      selected_samples = selected_samples,
      orf_long_result = orf_long_result,
      ko_lookup = get_ko_name_lookup(context_sqm)
    ),
    class = c("sqm_pathway_analysis", "list")
  )
  validate_pathway_analysis(pathway_analysis)
  pathway_analysis
}

resolve_pathway_analysis <- function(
    pathway_analysis,
    pathway_name,
    pathway_sqm,
    selected_samples,
    pathway_id,
    pathway_selection) {
  if (is.null(pathway_analysis)) {
    pathway_analysis <- build_pathway_analysis(
      list(
        pathway_name = pathway_name,
        pathway_id = pathway_id,
        pathway_selection = pathway_selection,
        pathway_sqm = pathway_sqm
      ),
      selected_samples
    )
  }
  validate_pathway_analysis(pathway_analysis)
  pathway_analysis
}

extract_ko_ec_lookup <- function(orf_long) {
  if (!"ko_id" %in% colnames(orf_long)) {
    stop("KO EC lookup requires column: ko_id", call. = FALSE)
  }
  if (!"ec_codes" %in% colnames(orf_long)) {
    return(
      orf_long |>
        transmute(ko_id = as.character(.data$ko_id)) |>
        distinct(.data$ko_id) |>
        mutate(ec_codes = NA_character_) |>
        arrange(.data$ko_id)
    )
  }

  orf_long |>
    transmute(
      ko_id = as.character(.data$ko_id),
      ec_code = str_split(as.character(.data$ec_codes), ";")
    ) |>
    tidyr::unnest_longer("ec_code", keep_empty = TRUE) |>
    mutate(
      ec_code = trimws(as.character(.data$ec_code)),
      ec_code = if_else(
        is.na(.data$ec_code) | !nzchar(.data$ec_code),
        NA_character_,
        .data$ec_code
      )
    ) |>
    group_by(.data$ko_id) |>
    summarise(
      ec_codes = {
        distinct_codes <- sort(unique(stats::na.omit(.data$ec_code)))
        if (length(distinct_codes) == 0L) {
          NA_character_
        } else {
          paste(distinct_codes, collapse = ";")
        }
      },
      .groups = "drop"
    ) |>
    arrange(.data$ko_id)
}

get_ko_dir_name <- function(ko_id, ko_ec) {
  if (is.null(ko_ec) || length(ko_ec) == 0L || is.na(ko_ec) || !nzchar(ko_ec)) {
    ko_id
  } else {
    paste0(ko_id, "_EC", sanitize_name(ko_ec))
  }
}

validate_percent_sum <- function(data, group_col, value_col, tolerance = 1e-6, expected = 100) {
  summary_tbl <- data |>
    group_by(.data[[group_col]]) |>
    summarise(total = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  bad <- summary_tbl |>
    filter(abs(.data$total - expected) > tolerance)
  if (nrow(bad) > 0L) {
    stop(
      "Controllo percentuali fallito per ", group_col, ".",
      call. = FALSE
    )
  }
}

build_ko_plot_table <- function(
    orf_long,
    selected_samples,
    top_n_ko,
    ko_lookup,
    pathway_name = NA_character_) {
  validate_positive_integer(top_n_ko, "top_n_ko")

  first_text <- function(values, fallback = NA_character_) {
    values <- trimws(as.character(values))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values) == 0L) fallback else values[[1L]]
  }

  selected_orfs <- orf_long |>
    filter(.data$sample %in% selected_samples)
  ko_ec_lookup <- if (nrow(selected_orfs) == 0L) {
    tibble::tibble(ko_id = character(), ec_codes = character())
  } else {
    extract_ko_ec_lookup(selected_orfs)
  }
  summary_tbl <- selected_orfs |>
    group_by(.data$sample, .data$ko_id) |>
    summarise(
      tpm = sum(.data$tpm),
      kegg_function = first_text(.data$kegg_function),
      .groups = "drop"
    )

  top_ko_ids <- summary_tbl |>
    group_by(.data$ko_id) |>
    summarise(total_tpm = sum(.data$tpm), .groups = "drop") |>
    arrange(desc(.data$total_tpm), .data$ko_id) |>
    slice_head(n = top_n_ko) |>
    pull(.data$ko_id)

  collapsed <- summary_tbl |>
    mutate(
      plot_ko_id = if_else(.data$ko_id %in% top_ko_ids, .data$ko_id, "Other"),
      plot_function = if_else(
        .data$ko_id %in% top_ko_ids,
        coalesce(
          .data$kegg_function,
          unname(as.character(ko_lookup[.data$ko_id])),
          .data$ko_id
        ),
        "Other KO outside top N"
      )
    ) |>
    group_by(.data$sample, .data$plot_ko_id) |>
    summarise(
      tpm = sum(.data$tpm),
      kegg_function = first_text(.data$plot_function),
      .groups = "drop"
    ) |>
    rename(ko_id = plot_ko_id)

  sample_totals <- collapsed |>
    group_by(.data$sample) |>
    summarise(sample_pathway_total_tpm = sum(.data$tpm), .groups = "drop")
  positive_samples <- selected_samples[selected_samples %in% sample_totals$sample]
  zero_samples <- selected_samples[!selected_samples %in% positive_samples]
  plot_ko_ids <- c(
    top_ko_ids,
    if (any(collapsed$ko_id == "Other")) "Other" else character()
  )

  ko_info <- collapsed |>
    filter(.data$ko_id %in% plot_ko_ids) |>
    group_by(.data$ko_id) |>
    summarise(
      kegg_function = first_text(.data$kegg_function),
      .groups = "drop"
    )

  positive_rows <- if (length(positive_samples) == 0L || length(plot_ko_ids) == 0L) {
    tibble::tibble(
      sample = character(),
      ko_id = character(),
      tpm = double(),
      kegg_function = character(),
      sample_pathway_total_tpm = double(),
      sample_pathway_percent = double(),
      denominator = double(),
      status = character(),
      plotted = logical()
    )
  } else {
    tidyr::expand_grid(sample = positive_samples, ko_id = plot_ko_ids) |>
      left_join(
        collapsed |>
          filter(.data$ko_id %in% plot_ko_ids),
        by = c("sample", "ko_id")
      ) |>
      left_join(ko_info, by = "ko_id", suffix = c("", "_info")) |>
      mutate(
        tpm = replace_na(.data$tpm, 0),
        kegg_function = coalesce(
          .data$kegg_function,
          .data$kegg_function_info,
          .data$ko_id
        )
      ) |>
      select("sample", "ko_id", "tpm", "kegg_function") |>
      left_join(sample_totals, by = "sample") |>
      mutate(
        sample_pathway_percent = 100 * .data$tpm / .data$sample_pathway_total_tpm,
        denominator = .data$sample_pathway_total_tpm,
        status = "ok",
        plotted = TRUE
      )
  }

  if (length(zero_samples) > 0L) {
    pathway_label <- if (
      length(pathway_name) == 1L && !is.na(pathway_name) && nzchar(pathway_name)
    ) {
      pathway_name
    } else {
      "<unknown>"
    }
    for (sample_name in zero_samples) {
      warning(
        "FUNZ TPM denominator is zero for pathway '", pathway_label,
        "', sample '", sample_name, "'.",
        call. = FALSE
      )
    }
  }
  zero_rows <- tibble::tibble(
    sample = zero_samples,
    ko_id = NA_character_,
    tpm = 0,
    kegg_function = NA_character_,
    sample_pathway_total_tpm = 0,
    sample_pathway_percent = NA_real_,
    denominator = 0,
    status = "zero_denominator",
    plotted = FALSE
  )
  plot_tbl <- bind_rows(positive_rows, zero_rows)

  if (nrow(positive_rows) > 0L) {
    validate_percent_sum(
      positive_rows,
      "sample",
      "sample_pathway_percent"
    )
  }

  metadata_keys <- plot_tbl |>
    transmute(sample = as.character(.data$sample), ko_id = as.character(.data$ko_id))
  metadata_mass <- plot_tbl |>
    group_by(.data$sample) |>
    summarise(tpm = sum(.data$tpm), .groups = "drop") |>
    arrange(.data$sample)
  plot_tbl <- plot_tbl |>
    left_join(
      ko_ec_lookup,
      by = "ko_id",
      relationship = "many-to-one",
      na_matches = "never"
    ) |>
    mutate(
      ec_codes = if_else(
        is.na(.data$ko_id) | .data$ko_id == "Other",
        NA_character_,
        .data$ec_codes
      )
    )
  joined_keys <- plot_tbl |>
    transmute(sample = as.character(.data$sample), ko_id = as.character(.data$ko_id))
  joined_mass <- plot_tbl |>
    group_by(.data$sample) |>
    summarise(tpm = sum(.data$tpm), .groups = "drop") |>
    arrange(.data$sample)
  if (
    nrow(plot_tbl) != nrow(metadata_keys) ||
      !identical(joined_keys, metadata_keys) ||
      !isTRUE(all.equal(
        joined_mass,
        metadata_mass,
        tolerance = 1e-10,
        check.attributes = FALSE
      ))
  ) {
    stop("FUNZ EC metadata join changed row keys or TPM mass.", call. = FALSE)
  }

  ko_levels <- plot_tbl |>
    filter(!is.na(.data$ko_id)) |>
    mutate(is_other = .data$ko_id == "Other") |>
    group_by(.data$ko_id, .data$is_other) |>
    summarise(total_tpm = sum(.data$tpm), .groups = "drop") |>
    arrange(.data$is_other, desc(.data$total_tpm), .data$ko_id) |>
    pull(.data$ko_id)

  plot_tbl |>
    mutate(
      sample = factor(.data$sample, levels = selected_samples),
      ko_id = factor(.data$ko_id, levels = ko_levels)
    ) |>
    arrange(.data$sample, desc(.data$plotted), desc(.data$tpm), .data$ko_id)
}

format_ko_sample_percent <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "NA%",
    x > 0 & x < 0.1 ~ "<0.1%",
    TRUE ~ paste0(formatC(x, format = "f", digits = 1), "%")
  )
}

build_ko_legend_labels <- function(plot_tbl, selected_samples) {
  plot_tbl |>
    filter(.data$plotted, !is.na(.data$ko_id)) |>
    mutate(
      sample = as.character(.data$sample),
      ko_id = as.character(.data$ko_id),
      ko_ec = replace_na(.data$ec_codes, "NA")
    ) |>
    arrange(.data$ko_id, match(.data$sample, selected_samples)) |>
    group_by(.data$ko_id, .data$ko_ec) |>
    summarise(
      percents = paste(format_ko_sample_percent(.data$sample_pathway_percent), collapse = " | "),
      .groups = "drop"
    ) |>
    mutate(legend_label = glue("{ko_id} / EC {ko_ec} | {percents}")) |>
    select("ko_id", "legend_label") |>
    tibble::deframe()
}

# Build visualizations only from the accompanying TSV tables written to disk.
make_ko_barplot <- function(plot_tbl, pathway_name, selected_samples) {
  plotted_rows <- plot_tbl |>
    filter(.data$plotted, !is.na(.data$ko_id))
  legend_labels <- build_ko_legend_labels(plotted_rows, selected_samples)
  ko_levels <- levels(plotted_rows$ko_id)
  non_other <- setdiff(ko_levels, "Other")
  palette <- c(
    setNames(rep(colors_hex, length.out = length(non_other)), non_other),
    if ("Other" %in% ko_levels) c(Other = "grey70") else NULL
  )

  ggplot(plotted_rows, aes(x = .data$sample, y = .data$tpm, fill = .data$ko_id)) +
    geom_col(color = "grey25", linewidth = 0.15, width = 0.78) +
    scale_fill_manual(values = palette, labels = legend_labels, drop = FALSE) +
    scale_x_discrete(limits = selected_samples, drop = FALSE) +
    scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    labs(
      title = paste0("KO barplot - ", pathway_name),
      x = "Sample",
      y = "TPM",
      fill = "KO / EC | % per sample"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = grid::unit(0.55, "lines")
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))
}

split_ec_code_field <- function(ec_codes) {
  codes <- str_split(ec_codes, ";", simplify = TRUE)
  codes <- trimws(as.character(codes))
  unique(codes[nzchar(codes)])
}

build_enzyme_plot_table <- function(sqm_object, selected_samples, enzyme_ecs) {
  enzyme_ecs <- normalize_enzyme_ecs(enzyme_ecs)
  kegg_tpm <- sqm_object$functions$KEGG$tpm
  if (is.null(kegg_tpm)) {
    stop("ENZIMI requires SQM$functions$KEGG$tpm.", call. = FALSE)
  }
  kegg_tpm <- as.data.frame(kegg_tpm, check.names = FALSE)
  validate_samples(selected_samples, colnames(kegg_tpm))
  validate_tpm_matrix(kegg_tpm[selected_samples], "SQM$functions$KEGG$tpm")
  if (is.null(rownames(kegg_tpm)) || anyDuplicated(rownames(kegg_tpm)) > 0L) {
    stop("SQM$functions$KEGG$tpm requires unique KO row names.", call. = FALSE)
  }

  matched_ecs <- build_ko_ec_map(sqm_object) |>
    filter(.data$ec_code %in% enzyme_ecs)

  tpm_by_enzyme <- kegg_tpm |>
    tibble::rownames_to_column("ko_id") |>
    pivot_longer(
      cols = all_of(selected_samples),
      names_to = "sample",
      values_to = "tpm"
    ) |>
    inner_join(matched_ecs, by = "ko_id", relationship = "many-to-many") |>
    group_by(.data$sample, .data$ec_code) |>
    summarise(tpm = sum(as.numeric(.data$tpm), na.rm = TRUE), .groups = "drop")

  tidyr::expand_grid(sample = selected_samples, ec_code = enzyme_ecs) |>
    left_join(tpm_by_enzyme, by = c("sample", "ec_code")) |>
    mutate(
      tpm = replace_na(.data$tpm, 0)
    ) |>
    group_by(.data$ec_code) |>
    mutate(
      enzyme_total_tpm = sum(.data$tpm),
      status = if_else(
        .data$enzyme_total_tpm > 0,
        "positive_tpm",
        "no_positive_tpm"
      ),
      plotted = .data$enzyme_total_tpm > 0
    ) |>
    ungroup() |>
    mutate(
      sample_order = match(.data$sample, selected_samples),
      sample = factor(.data$sample, levels = selected_samples),
      ec_code = factor(.data$ec_code, levels = enzyme_ecs)
    )
}

enzyme_palette <- function(enzyme_ecs) {
  stats::setNames(rep(colors_hex, length.out = length(enzyme_ecs)), enzyme_ecs)
}

make_enzyme_barplot <- function(enzyme_tbl, title) {
  plotted_tbl <- enzyme_tbl |>
    filter(.data$plotted) |>
    mutate(ec_code = droplevels(.data$ec_code))
  if (nrow(plotted_tbl) == 0L) {
    stop("Cannot build an enzyme barplot without positive TPM.", call. = FALSE)
  }
  enzyme_ecs <- levels(plotted_tbl$ec_code)
  ggplot(plotted_tbl, aes(x = .data$sample, y = .data$tpm, fill = .data$ec_code)) +
    geom_col(
      position = position_dodge2(preserve = "single"),
      color = "grey25",
      linewidth = 0.15,
      width = 0.78
    ) +
    scale_fill_manual(values = enzyme_palette(enzyme_ecs), drop = FALSE) +
    scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    labs(title = title, x = "Sample", y = "TPM", fill = "EC") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

make_enzyme_lineplot <- function(enzyme_tbl, title) {
  plotted_tbl <- enzyme_tbl |>
    filter(.data$plotted) |>
    mutate(ec_code = droplevels(.data$ec_code))
  if (nrow(plotted_tbl) == 0L) {
    stop("Cannot build an enzyme line plot without positive TPM.", call. = FALSE)
  }
  enzyme_ecs <- levels(plotted_tbl$ec_code)
  ggplot(plotted_tbl, aes(x = .data$sample, y = .data$tpm, color = .data$ec_code, group = .data$ec_code)) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2) +
    scale_color_manual(values = enzyme_palette(enzyme_ecs), drop = FALSE) +
    scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    labs(title = title, x = "Sample", y = "TPM", color = "EC") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

build_flow_ko_metadata <- function(orf_long, ko_lookup) {
  required_columns <- c("ko_id", "kegg_function")
  missing_columns <- setdiff(required_columns, colnames(orf_long))
  if (length(missing_columns) > 0L) {
    stop(
      "FLOW KO metadata requires columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  ko_metadata <- orf_long |>
    transmute(
      ko_id = as.character(.data$ko_id),
      description = trimws(as.character(.data$kegg_function))
    ) |>
    mutate(
      description = if_else(
        is.na(.data$description) | !nzchar(.data$description),
        NA_character_,
        .data$description
      )
    ) |>
    group_by(.data$ko_id) |>
    summarise(
      KO_name = {
        descriptions <- sort(unique(stats::na.omit(.data$description)))
        if (length(descriptions) == 0L) {
          NA_character_
        } else {
          paste(descriptions, collapse = "; ")
        }
      },
      .groups = "drop"
    ) |>
    arrange(.data$ko_id)

  lookup_names <- unname(as.character(ko_lookup[ko_metadata$ko_id]))
  lookup_names <- trimws(lookup_names)
  lookup_names[is.na(lookup_names) | !nzchar(lookup_names)] <- NA_character_
  ko_metadata$KO_name <- coalesce(
    ko_metadata$KO_name,
    lookup_names,
    ko_metadata$ko_id
  )
  ko_metadata <- ko_metadata |>
    left_join(
      extract_ko_ec_lookup(orf_long),
      by = "ko_id",
      relationship = "one-to-one"
    )

  if (anyDuplicated(ko_metadata$ko_id) > 0L) {
    stop("FLOW KO metadata keys must be unique.", call. = FALSE)
  }

  ko_metadata
}

join_flow_ko_metadata <- function(summary_tbl, ko_meta, tolerance = 1e-10) {
  required_summary_columns <- c("sample", "taxon", "KO", "TPM")
  missing_summary_columns <- setdiff(
    required_summary_columns,
    colnames(summary_tbl)
  )
  if (length(missing_summary_columns) > 0L) {
    stop(
      "FLOW summary is missing columns: ",
      paste(missing_summary_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!all(c("ko_id", "KO_name", "ec_codes") %in% colnames(ko_meta))) {
    stop(
      "FLOW KO metadata must contain ko_id, KO_name, and ec_codes.",
      call. = FALSE
    )
  }
  if (anyDuplicated(ko_meta$ko_id) > 0L) {
    stop(
      "FLOW KO metadata keys must be unique for a many-to-one join.",
      call. = FALSE
    )
  }

  key_columns <- c("sample", "taxon", "KO")
  if (anyDuplicated(summary_tbl[key_columns]) > 0L) {
    stop("FLOW summary keys must be unique before metadata join.", call. = FALSE)
  }
  before_keys <- summary_tbl[key_columns]
  before_mass <- summary_tbl |>
    group_by(.data$sample) |>
    summarise(TPM = sum(.data$TPM), .groups = "drop") |>
    arrange(.data$sample)

  joined <- summary_tbl |>
    left_join(
      ko_meta,
      by = c("KO" = "ko_id"),
      relationship = "many-to-one"
    )

  after_mass <- joined |>
    group_by(.data$sample) |>
    summarise(TPM = sum(.data$TPM), .groups = "drop") |>
    arrange(.data$sample)
  keys_unchanged <- identical(joined[key_columns], before_keys)
  mass_unchanged <- identical(after_mass$sample, before_mass$sample) &&
    isTRUE(all.equal(
      after_mass$TPM,
      before_mass$TPM,
      tolerance = tolerance,
      check.attributes = FALSE
    ))

  if (nrow(joined) != nrow(summary_tbl) || !keys_unchanged || !mass_unchanged) {
    stop(
      "FLOW metadata join postcondition failed: keys, rows, or TPM mass changed.",
      call. = FALSE
    )
  }

  joined
}

build_flow_table_for_rank <- function(orf_long, rank, selected_samples, top_n_taxa, top_n_ko, ko_lookup) {
  selected_orfs <- orf_long |>
    filter(.data$sample %in% selected_samples)
  top_taxa <- select_top_classified_taxa(
    data = selected_orfs,
    taxon_col = rank,
    value_col = "tpm",
    top_n = top_n_taxa
  )

  top_kos <- orf_long |>
    filter(.data$sample %in% selected_samples) |>
    group_by(.data$ko_id) |>
    summarise(total_tpm = sum(.data$tpm), .groups = "drop") |>
    arrange(desc(.data$total_tpm), .data$ko_id) |>
    slice_head(n = top_n_ko) |>
    pull(.data$ko_id)

  summary_tbl <- selected_orfs |>
    mutate(
      taxon = collapse_taxa_preserving_unclassified(.data[[rank]], top_taxa),
      KO = if_else(.data$ko_id %in% top_kos, .data$ko_id, "Other")
    ) |>
    group_by(.data$sample, .data$taxon, .data$KO) |>
    summarise(TPM = sum(.data$tpm), .groups = "drop") |>
    filter(.data$TPM > 0)

  ko_meta <- build_flow_ko_metadata(orf_long, ko_lookup)

  summary_tbl |>
    join_flow_ko_metadata(ko_meta) |>
    mutate(
      KO_name = if_else(
        .data$KO == "Other",
        "Other KOs",
        coalesce(.data$KO_name, .data$KO)
      )
    )
}

build_flow_table_for_sample <- function(flow_rank_table, pathway_name, rank, sample_name) {
  sample_tbl <- flow_rank_table |>
    filter(.data$sample == sample_name) |>
    group_by(.data$taxon, .data$KO, .data$KO_name, .data$ec_codes) |>
    summarise(TPM = sum(.data$TPM), .groups = "drop") |>
    filter(.data$TPM > 0)

  if (nrow(sample_tbl) == 0L) {
    return(sample_tbl)
  }

  sample_total <- sum(sample_tbl$TPM)
  tax_totals <- sample_tbl |>
    group_by(.data$taxon) |>
    summarise(taxon_tpm = sum(.data$TPM), .groups = "drop")
  ko_totals <- sample_tbl |>
    group_by(.data$KO) |>
    summarise(ko_tpm = sum(.data$TPM), .groups = "drop")

  taxon_levels <- tax_totals |>
    mutate(is_other = .data$taxon == "Other") |>
    arrange(.data$is_other, .data$taxon_tpm, .data$taxon) |>
    pull(.data$taxon)

  ko_levels <- ko_totals |>
    mutate(is_other = .data$KO == "Other") |>
    arrange(.data$is_other, .data$ko_tpm, .data$KO) |>
    pull(.data$KO)

  flow_tbl <- sample_tbl |>
    left_join(tax_totals, by = "taxon") |>
    left_join(ko_totals, by = "KO") |>
    mutate(
      pathway = pathway_name,
      rank = rank,
      sample = sample_name,
      flow_percent = 100 * .data$TPM / sample_total,
      taxon_percent = 100 * .data$taxon_tpm / sample_total,
      KO_percent = 100 * .data$ko_tpm / sample_total,
      taxon = factor(.data$taxon, levels = taxon_levels),
      KO = factor(.data$KO, levels = ko_levels)
    ) |>
    select(
      "pathway", "rank", "sample", "taxon", "KO",
      "KO_name", "ec_codes", "TPM", "taxon_percent", "KO_percent",
      "flow_percent"
    )

  validate_percent_sum(flow_tbl, "sample", "flow_percent")
  flow_tbl
}

build_flow_legend_spec <- function(flow_tbl, tolerance = 1e-6) {
  required_columns <- c(
    "taxon", "KO", "ec_codes", "taxon_percent", "KO_percent"
  )
  missing_columns <- setdiff(required_columns, colnames(flow_tbl))
  if (length(missing_columns) > 0L) {
    stop(
      "FLOW legend data requires columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(flow_tbl) == 0L) {
    stop("FLOW legend data cannot be empty.", call. = FALSE)
  }

  taxon_levels <- if (is.factor(flow_tbl$taxon)) {
    levels(flow_tbl$taxon)
  } else {
    unique(as.character(flow_tbl$taxon))
  }
  ko_levels <- if (is.factor(flow_tbl$KO)) {
    levels(flow_tbl$KO)
  } else {
    unique(as.character(flow_tbl$KO))
  }

  taxonomy <- flow_tbl |>
    transmute(
      node = as.character(.data$taxon),
      percent = as.numeric(.data$taxon_percent)
    ) |>
    distinct(.data$node, .data$percent) |>
    arrange(match(.data$node, taxon_levels))
  functional <- flow_tbl |>
    transmute(
      node = as.character(.data$KO),
      ec_codes = as.character(.data$ec_codes),
      percent = as.numeric(.data$KO_percent)
    ) |>
    distinct(.data$node, .data$ec_codes, .data$percent) |>
    arrange(match(.data$node, ko_levels))

  if (anyDuplicated(taxonomy$node) > 0L || anyDuplicated(functional$node) > 0L) {
    stop(
      "FLOW legend nodes must have one stable percentage and EC association.",
      call. = FALSE
    )
  }
  if (
    anyNA(taxonomy$percent) || anyNA(functional$percent) ||
      abs(sum(taxonomy$percent) - 100) > tolerance ||
      abs(sum(functional$percent) - 100) > tolerance
  ) {
    stop("FLOW legend percentages must sum to 100 per section.", call. = FALSE)
  }

  taxonomy <- taxonomy |>
    mutate(
      display_name = .data$node,
      label = paste0(
        .data$display_name,
        " | ",
        format_ko_sample_percent(.data$percent)
      )
    )
  functional <- functional |>
    mutate(
      display_name = case_when(
        .data$node == "Other" ~ "Other KOs",
        is.na(.data$ec_codes) | !nzchar(trimws(.data$ec_codes)) ~
          paste0(.data$node, " / EC NA"),
        TRUE ~ paste0(.data$node, " / EC ", .data$ec_codes)
      ),
      label = paste0(
        .data$display_name,
        " | ",
        format_ko_sample_percent(.data$percent)
      )
    )

  list(taxonomy = taxonomy, functional = functional)
}

build_flow_color_map <- function(taxon_levels, ko_levels) {
  categories <- unique(c(as.character(taxon_levels), as.character(ko_levels)))
  categories <- categories[!is.na(categories) & nzchar(categories)]
  non_other <- setdiff(categories, "Other")
  base_colors <- unique(as.character(colors_hex))
  base_colors <- base_colors[!is.na(base_colors) & nzchar(base_colors)]

  if (length(non_other) <= length(base_colors)) {
    category_colors <- base_colors[seq_along(non_other)]
  } else {
    category_colors <- grDevices::hcl.colors(length(non_other), palette = "Dynamic")
    if (anyNA(category_colors) || anyDuplicated(toupper(category_colors))) {
      stop(
        "Unable to generate distinct FLOW colors for ",
        length(non_other), " displayed categories.",
        call. = FALSE
      )
    }
  }

  stats::setNames(
    c(category_colors, if ("Other" %in% categories) "grey70" else character()),
    c(non_other, if ("Other" %in% categories) "Other" else character())
  )
}

make_flow_plot_with_legends <- function(flow_tbl, pathway_name, rank, sample_name) {
  taxon_levels <- levels(flow_tbl$taxon)
  ko_levels <- levels(flow_tbl$KO)
  legend_spec <- build_flow_legend_spec(flow_tbl)
  lodes_tbl <- ggalluvial::to_lodes_form(flow_tbl, axes = c("taxon", "KO"), discern = FALSE) |>
    mutate(
      x = factor(.data$x, levels = c("taxon", "KO"), labels = c("Taxon", "KO")),
      stratum = as.character(.data$stratum)
    )
  flow_palette <- build_flow_color_map(taxon_levels, ko_levels)

  ggplot(
    flow_tbl,
    aes(axis1 = .data$taxon, axis2 = .data$KO, y = .data$flow_percent)
  ) +
    ggalluvial::geom_alluvium(aes(fill = .data$taxon), alpha = 0.78, width = 1 / 12) +
    ggalluvial::geom_stratum(
      data = dplyr::filter(lodes_tbl, .data$x == "Taxon"),
      aes(x = .data$x, stratum = .data$stratum, alluvium = .data$alluvium, y = .data$flow_percent, fill = after_stat(stratum)),
      inherit.aes = FALSE,
      width = 1 / 5,
      color = "grey35",
      linewidth = 0.25
    ) +
    scale_fill_manual(
      name = "Taxonomy | % of sample",
      values = flow_palette,
      breaks = legend_spec$taxonomy$node,
      labels = legend_spec$taxonomy$label,
      drop = FALSE
    ) +
    guides(fill = guide_legend(order = 1, ncol = 1, byrow = TRUE)) +
    ggnewscale::new_scale_fill() +
    ggalluvial::geom_stratum(
      data = dplyr::filter(lodes_tbl, .data$x == "KO"),
      aes(x = .data$x, stratum = .data$stratum, alluvium = .data$alluvium, y = .data$flow_percent, fill = after_stat(stratum)),
      inherit.aes = FALSE,
      width = 1 / 5,
      color = "grey35",
      linewidth = 0.25
    ) +
    scale_fill_manual(
      name = "Function (KO / EC) | % of sample",
      values = flow_palette,
      breaks = legend_spec$functional$node,
      labels = legend_spec$functional$label,
      drop = FALSE
    ) +
    guides(fill = guide_legend(order = 2, ncol = 1, byrow = TRUE)) +
    geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(stratum)),
      size = 2.8
    ) +
    scale_x_discrete(limits = c("Taxon", "KO"), expand = c(0.08, 0.08)) +
    labs(
      title = paste0("Flowplot - ", pathway_name, " - ", rank, " - ", sample_name),
      subtitle = "Taxonomy-to-KO/EC flow from the ORF x sample x KO table",
      x = NULL,
      y = "Relative flow (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = grid::unit(0.55, "lines"),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

make_flow_plot <- function(flow_tbl, pathway_name, rank, sample_name) {
  taxon_levels <- levels(flow_tbl$taxon)
  ko_levels <- levels(flow_tbl$KO)
  ko_display_levels <- ifelse(ko_levels == "Other", "Other KOs", ko_levels)
  plot_tbl <- flow_tbl |>
    mutate(
      taxon = factor(as.character(.data$taxon), levels = taxon_levels),
      KO_display = if_else(
        as.character(.data$KO) == "Other",
        "Other KOs",
        as.character(.data$KO)
      ),
      KO_display = factor(.data$KO_display, levels = ko_display_levels)
    )
  taxon_palette <- build_flow_color_map(taxon_levels, character())

  ggplot(
    plot_tbl,
    aes(axis1 = .data$taxon, axis2 = .data$KO_display, y = .data$flow_percent)
  ) +
    ggalluvial::geom_alluvium(
      aes(fill = .data$taxon),
      alpha = 0.72,
      width = 1 / 12
    ) +
    ggalluvial::geom_stratum(
      width = 1 / 5,
      fill = "grey95",
      color = "grey35",
      linewidth = 0.25
    ) +
    geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(stratum)),
      size = 2.8
    ) +
    scale_fill_manual(values = taxon_palette, drop = FALSE, guide = "none") +
    scale_x_discrete(limits = c("Taxon", "KO"), expand = c(0.08, 0.08)) +
    labs(
      title = paste0("Flowplot - ", pathway_name, " - ", rank, " - ", sample_name),
      subtitle = "Taxonomy-to-KO flow from the ORF x sample x KO table",
      x = NULL,
      y = "Relative flow (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

make_flow_sankey_with_detailed_labels <- function(flow_tbl, pathway_name, rank, sample_name) {
  tax_labels <- levels(flow_tbl$taxon)
  ko_labels <- levels(flow_tbl$KO)
  legend_spec <- build_flow_legend_spec(flow_tbl)
  tax_node_labels <- stats::setNames(
    legend_spec$taxonomy$label,
    legend_spec$taxonomy$node
  )
  functional_node_labels <- stats::setNames(
    legend_spec$functional$label,
    legend_spec$functional$node
  )
  taxon_percents <- stats::setNames(
    legend_spec$taxonomy$percent,
    legend_spec$taxonomy$node
  )
  functional_percents <- stats::setNames(
    legend_spec$functional$percent,
    legend_spec$functional$node
  )
  functional_names <- stats::setNames(
    legend_spec$functional$display_name,
    legend_spec$functional$node
  )
  flow_palette <- build_flow_color_map(tax_labels, ko_labels)
  tax_y <- if (length(tax_labels) == 1L) 0.5 else seq(0.98, 0.02, length.out = length(tax_labels))
  ko_y <- if (length(ko_labels) == 1L) 0.5 else seq(0.98, 0.02, length.out = length(ko_labels))

  plotly::plot_ly(
    type = "sankey",
    arrangement = "fixed",
    valueformat = ".2f",
    valuesuffix = "%",
    node = list(
      pad = 18,
      thickness = 18,
      line = list(color = "rgba(70,70,70,0.35)", width = 0.5),
      label = c(
        unname(tax_node_labels[tax_labels]),
        unname(functional_node_labels[ko_labels])
      ),
      color = unname(flow_palette[c(tax_labels, ko_labels)]),
      x = c(rep(0.02, length(tax_labels)), rep(0.98, length(ko_labels))),
      y = c(tax_y, ko_y)
    ),
    link = list(
      source = match(as.character(flow_tbl$taxon), tax_labels) - 1L,
      target = length(tax_labels) + match(as.character(flow_tbl$KO), ko_labels) - 1L,
      value = flow_tbl$flow_percent,
      color = unname(grDevices::adjustcolor(flow_palette[as.character(flow_tbl$taxon)], alpha.f = 0.65)),
      customdata = paste0(
        "Taxon: ", flow_tbl$taxon,
        "<br>Taxon share: ",
        format_ko_sample_percent(unname(taxon_percents[as.character(flow_tbl$taxon)])),
        "<br>Function: ", unname(functional_names[as.character(flow_tbl$KO)]),
        "<br>Function share: ",
        format_ko_sample_percent(unname(functional_percents[as.character(flow_tbl$KO)])),
        "<br>Function name: ", flow_tbl$KO_name,
        "<br>TPM: ", sprintf("%.3f", flow_tbl$TPM),
        "<br>Flow: ", sprintf("%.2f", flow_tbl$flow_percent), "%"
      ),
      hovertemplate = "%{customdata}<extra></extra>"
    )
  ) |>
    plotly::layout(
      title = list(
        text = paste0(
          "Flowplot - ", pathway_name, " - ", rank, " - ", sample_name,
          "<br><sup>ORF-linked taxon -> KO / EC flow based on TPM</sup>"
        )
      ),
      font = list(size = 11),
      margin = list(l = 20, r = 20, t = 60, b = 20)
    )
}

make_flow_sankey <- function(flow_tbl, pathway_name, rank, sample_name) {
  tax_labels <- levels(flow_tbl$taxon)
  ko_labels <- levels(flow_tbl$KO)
  ko_display_labels <- ifelse(ko_labels == "Other", "Other KOs", ko_labels)
  taxon_palette <- build_flow_color_map(tax_labels, character())
  taxon_percents <- stats::setNames(
    flow_tbl$taxon_percent[match(tax_labels, as.character(flow_tbl$taxon))],
    tax_labels
  )
  functional_percents <- stats::setNames(
    flow_tbl$KO_percent[match(ko_labels, as.character(flow_tbl$KO))],
    ko_labels
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
      label = c(tax_labels, ko_display_labels),
      color = rep("rgba(245,245,245,1)", length(tax_labels) + length(ko_labels))
    ),
    link = list(
      source = match(as.character(flow_tbl$taxon), tax_labels) - 1L,
      target = length(tax_labels) + match(as.character(flow_tbl$KO), ko_labels) - 1L,
      value = flow_tbl$flow_percent,
      color = unname(grDevices::adjustcolor(
        taxon_palette[as.character(flow_tbl$taxon)],
        alpha.f = 0.65
      )),
      customdata = paste0(
        "Taxon: ", flow_tbl$taxon,
        "<br>Taxon share: ",
        format_ko_sample_percent(unname(taxon_percents[as.character(flow_tbl$taxon)])),
        "<br>KO: ", flow_tbl$KO,
        "<br>EC: ", ifelse(is.na(flow_tbl$ec_codes), "NA", flow_tbl$ec_codes),
        "<br>Function: ", flow_tbl$KO_name,
        "<br>Function share: ",
        format_ko_sample_percent(unname(functional_percents[as.character(flow_tbl$KO)])),
        "<br>TPM: ", sprintf("%.3f", flow_tbl$TPM),
        "<br>Flow: ", sprintf("%.2f", flow_tbl$flow_percent), "%"
      ),
      hovertemplate = "%{customdata}<extra></extra>"
    )
  ) |>
    plotly::layout(
      title = list(
        text = paste0(
          "Flowplot - ", pathway_name, " - ", rank, " - ", sample_name,
          "<br><sup>Taxonomy-to-KO flow based on TPM</sup>"
        )
      ),
      font = list(size = 11),
      margin = list(l = 20, r = 20, t = 60, b = 20)
    )
}

make_taxonomy_plot <- function(
    sqm_object,
    rank,
    count,
    selected_samples,
    top_n_taxa,
    ignore_unmapped,
    ignore_unclassified) {
  plot_object <- SQMtools::plotTaxonomy(
    SQM = sqm_object,
    rank = rank,
    count = count,
    N = top_n_taxa,
    samples = selected_samples,
    ignore_unmapped = ignore_unmapped,
    ignore_unclassified = ignore_unclassified,
    no_partial_classifications = FALSE,
    rescale = FALSE
  )

  plot_object +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

extract_taxonomy_plot_data <- function(plot_object, count) {
  if (is.null(plot_object$data)) {
    stop("plotTaxonomy() did not expose the plot data.", call. = FALSE)
  }

  plot_data <- as_tibble(plot_object$data)
  names(plot_data)[names(plot_data) == "item"] <- "taxon"
  names(plot_data)[names(plot_data) == "abun"] <- "value"
  plot_data <- plot_data |>
    mutate(
      count = count,
      sample = as.character(.data$sample),
      taxon = as.character(.data$taxon),
      value = as.numeric(.data$value)
    )

  plot_data
}

taxonomy_category_subsets <- function(categories) {
  categories <- unique(as.character(categories))
  subsets <- list(character())
  if (length(categories) == 0L) {
    return(subsets)
  }
  for (subset_size in seq_along(categories)) {
    combinations <- utils::combn(categories, subset_size, simplify = FALSE)
    subsets <- c(subsets, combinations)
  }
  subsets
}

resolve_effective_taxonomy_exclusions <- function(
    percent_frame,
    displayed_percent_sum,
    selected_samples,
    requested_excluded_categories = c("Unmapped", "Unclassified"),
    tolerance = 1e-6) {
  selected_samples <- as.character(selected_samples)
  requested_excluded_categories <- unique(as.character(requested_excluded_categories))
  displayed_percent_sum <- as.numeric(displayed_percent_sum[selected_samples])
  names(displayed_percent_sum) <- selected_samples

  if (anyNA(displayed_percent_sum) || any(!is.finite(displayed_percent_sum))) {
    stop(
      "Displayed global taxonomy percentages are missing or non-finite for selected samples.",
      call. = FALSE
    )
  }

  raw_percent_sum <- colSums(percent_frame[, selected_samples, drop = FALSE])
  unexplained_percent <- raw_percent_sum - displayed_percent_sum
  unexplained_percent[abs(unexplained_percent) <= tolerance] <- 0
  if (any(unexplained_percent < -tolerance)) {
    stop(
      "Global taxonomy percentages cannot be reconciled: displayed values exceed raw percentages.",
      call. = FALSE
    )
  }

  present_requested <- intersect(
    requested_excluded_categories,
    rownames(percent_frame)
  )
  positive_requested <- present_requested[vapply(
    present_requested,
    function(category) {
      any(percent_frame[category, selected_samples, drop = TRUE] > tolerance)
    },
    logical(1)
  )]
  candidate_subsets <- taxonomy_category_subsets(positive_requested)
  candidate_matches <- vapply(
    candidate_subsets,
    function(categories) {
      candidate_percent <- if (length(categories) == 0L) {
        stats::setNames(rep(0, length(selected_samples)), selected_samples)
      } else {
        colSums(percent_frame[categories, selected_samples, drop = FALSE])
      }
      all(abs(candidate_percent[selected_samples] - unexplained_percent[selected_samples]) <= tolerance)
    },
    logical(1)
  )

  matching_subsets <- candidate_subsets[candidate_matches]
  if (length(matching_subsets) == 0L) {
    stop(
      "Global taxonomy percentages cannot be reconciled with requested exclusions.",
      call. = FALSE
    )
  }
  if (length(matching_subsets) > 1L) {
    stop(
      "Global taxonomy effective exclusions are ambiguous across selected samples.",
      call. = FALSE
    )
  }

  effective_categories <- matching_subsets[[1L]]
  retained_categories <- setdiff(positive_requested, effective_categories)
  excluded_percent <- if (length(effective_categories) == 0L) {
    stats::setNames(rep(0, length(selected_samples)), selected_samples)
  } else {
    colSums(percent_frame[effective_categories, selected_samples, drop = FALSE])
  }
  accounted_percent_sum <- displayed_percent_sum + excluded_percent
  if (any(abs(accounted_percent_sum - raw_percent_sum) > tolerance)) {
    stop(
      "Global taxonomy accounted percentages do not match raw percentages.",
      call. = FALSE
    )
  }

  resolution_status <- if (length(positive_requested) == 0L) {
    "no_positive_requested_exclusions"
  } else if (length(retained_categories) > 0L) {
    "requested_retained"
  } else {
    "matched_requested"
  }

  list(
    raw_percent_sum = raw_percent_sum,
    displayed_percent_sum = displayed_percent_sum,
    excluded_percent = excluded_percent,
    accounted_percent_sum = accounted_percent_sum,
    requested_excluded_categories = requested_excluded_categories,
    effective_excluded_categories = effective_categories,
    retained_requested_categories = retained_categories,
    exclusion_resolution_status = resolution_status
  )
}

add_global_taxonomy_percent_metadata <- function(
    plot_data,
    sqm_object,
    rank,
    selected_samples,
    excluded_categories = c("Unmapped", "Unclassified"),
    tolerance = 1e-6,
    context_label = "global") {
  required_columns <- c("sample", "value", "count")
  missing_columns <- setdiff(required_columns, colnames(plot_data))
  if (length(missing_columns) > 0L) {
    stop(
      "Global taxonomy plot data is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!all(as.character(plot_data$count) == "percent")) {
    stop("Global taxonomy metadata requires count='percent'.", call. = FALSE)
  }

  percent_matrix <- sqm_object$taxa[[rank]]$percent
  if (is.null(percent_matrix)) {
    stop("Missing SQM taxonomy percent matrix for rank: ", rank, call. = FALSE)
  }
  percent_frame <- as.data.frame(percent_matrix, check.names = FALSE)
  validate_samples(selected_samples, colnames(percent_frame))
  validate_tpm_matrix(percent_frame[selected_samples], paste0("sqm$taxa$", rank, "$percent"))

  total_reads <- as.numeric(sqm_object$total_reads)
  total_read_names <- names(sqm_object$total_reads)
  if (is.null(total_read_names) && !is.null(sqm_object$misc$samples) &&
      length(total_reads) == length(sqm_object$misc$samples)) {
    total_read_names <- as.character(sqm_object$misc$samples)
  }
  if (is.null(total_read_names)) {
    stop("sqm$total_reads must be named by sample.", call. = FALSE)
  }
  names(total_reads) <- total_read_names
  denominator_values <- total_reads[selected_samples]
  if (
    anyNA(denominator_values) || any(!is.finite(denominator_values)) ||
      any(denominator_values <= 0)
  ) {
    stop("sqm$total_reads must be finite and positive for selected samples.", call. = FALSE)
  }

  displayed_by_sample <- plot_data |>
    mutate(sample = as.character(.data$sample)) |>
    group_by(.data$sample) |>
    summarise(displayed_percent_sum = sum(.data$value), .groups = "drop")
  unknown_samples <- setdiff(displayed_by_sample$sample, selected_samples)
  if (length(unknown_samples) > 0L) {
    stop(
      "Global taxonomy plot data contains unexpected samples: ",
      paste(unknown_samples, collapse = ", "),
      call. = FALSE
    )
  }
  displayed_percent_sum <- stats::setNames(
    displayed_by_sample$displayed_percent_sum,
    displayed_by_sample$sample
  )
  exclusion_audit <- resolve_effective_taxonomy_exclusions(
    percent_frame = percent_frame,
    displayed_percent_sum = displayed_percent_sum,
    selected_samples = selected_samples,
    requested_excluded_categories = excluded_categories,
    tolerance = tolerance
  )
  serialize_categories <- function(categories) {
    if (length(categories) == 0L) NA_character_ else paste(categories, collapse = ";")
  }
  if (length(exclusion_audit$retained_requested_categories) > 0L) {
    warning(
      "SQMtools::plotTaxonomy retained requested exclusion categories in displayed data",
      " | context=", as.character(context_label),
      " | rank=", as.character(rank),
      ": ", paste(exclusion_audit$retained_requested_categories, collapse = ", "),
      "; values may be included in Other.",
      call. = FALSE
    )
  }
  audit_by_sample <- tibble::tibble(
    sample = selected_samples,
    denominator_value = as.numeric(denominator_values),
    raw_percent_sum = as.numeric(exclusion_audit$raw_percent_sum[selected_samples]),
    displayed_percent_sum = as.numeric(
      exclusion_audit$displayed_percent_sum[selected_samples]
    ),
    excluded_percent = as.numeric(exclusion_audit$excluded_percent[selected_samples]),
    accounted_percent_sum = as.numeric(
      exclusion_audit$accounted_percent_sum[selected_samples]
    ),
    requested_excluded_categories = serialize_categories(
      exclusion_audit$requested_excluded_categories
    ),
    effective_excluded_categories = serialize_categories(
      exclusion_audit$effective_excluded_categories
    ),
    retained_requested_categories = serialize_categories(
      exclusion_audit$retained_requested_categories
    ),
    exclusion_resolution_status = exclusion_audit$exclusion_resolution_status,
    excluded_categories = serialize_categories(
      exclusion_audit$effective_excluded_categories
    )
  )

  plot_data |>
    mutate(sample = as.character(.data$sample)) |>
    left_join(audit_by_sample, by = "sample", relationship = "many-to-one") |>
    mutate(
      percent_of_total_library = .data$value,
      denominator_type = "total_reads"
    )
}

build_pathway_taxonomy_percent_table <- function(
    sqm_object,
    rank,
    selected_samples,
    top_n_taxa,
    pathway_name) {
  tax_table <- as.data.frame(sqm_object$orfs$tax, check.names = FALSE)
  tpm_table <- as.data.frame(sqm_object$orfs$tpm, check.names = FALSE)

  if (!rank %in% colnames(tax_table)) {
    stop("Taxonomic rank not found in ORF taxonomy: ", rank, call. = FALSE)
  }
  validate_samples(selected_samples, colnames(tpm_table))
  validate_positive_integer(top_n_taxa, "top_n_taxa")

  tax_ids <- rownames(tax_table)
  tpm_ids <- rownames(tpm_table)
  if (is.null(tax_ids) || is.null(tpm_ids) ||
      anyDuplicated(tax_ids) > 0L || anyDuplicated(tpm_ids) > 0L) {
    stop("ORF taxonomy and TPM tables must have unique row names.", call. = FALSE)
  }
  if (!setequal(tax_ids, tpm_ids)) {
    stop("ORF taxonomy and TPM tables contain different orf_id keys.", call. = FALSE)
  }

  tax_long <- tibble(
    orf_id = tax_ids,
    taxon = normalize_taxon_value(tax_table[[rank]])
  )
  if ("Other" %in% tax_long$taxon) {
    stop("ORF taxonomy contains the reserved label 'Other'.", call. = FALSE)
  }

  tpm_long <- tpm_table |>
    rownames_to_column("orf_id") |>
    select("orf_id", all_of(selected_samples)) |>
    pivot_longer(
      cols = all_of(selected_samples),
      names_to = "sample",
      values_to = "tpm"
    ) |>
    mutate(tpm = as.numeric(.data$tpm))

  invalid_tpm <- !is.na(tpm_long$tpm) &
    (!is.finite(tpm_long$tpm) | tpm_long$tpm < 0)
  if (any(invalid_tpm)) {
    stop("ORF TPM values must be finite and non-negative.", call. = FALSE)
  }

  base_table <- tpm_long |>
    mutate(tpm = replace_na(.data$tpm, 0)) |>
    left_join(tax_long, by = "orf_id")

  if (nrow(base_table) != nrow(tpm_long) || anyNA(base_table$taxon)) {
    stop("ORF taxonomy join did not preserve every TPM row.", call. = FALSE)
  }

  taxon_totals <- base_table |>
    group_by(.data$sample, .data$taxon) |>
    summarise(taxon_tpm = sum(.data$tpm), .groups = "drop")
  sample_denominators <- base_table |>
    group_by(.data$sample) |>
    summarise(pathway_tpm = sum(.data$tpm), .groups = "drop")

  reserved_taxa <- c("Unclassified", "Unmapped")
  top_taxa <- taxon_totals |>
    filter(!.data$taxon %in% reserved_taxa, .data$taxon_tpm > 0) |>
    group_by(.data$taxon) |>
    summarise(total_tpm = sum(.data$taxon_tpm), .groups = "drop") |>
    arrange(desc(.data$total_tpm), .data$taxon) |>
    slice_head(n = top_n_taxa) |>
    pull(.data$taxon)

  positive_rows <- taxon_totals |>
    left_join(sample_denominators, by = "sample") |>
    filter(.data$pathway_tpm > 0, .data$taxon_tpm > 0) |>
    mutate(
      taxon = case_when(
        .data$taxon %in% reserved_taxa ~ .data$taxon,
        .data$taxon %in% top_taxa ~ .data$taxon,
        TRUE ~ "Other"
      )
    ) |>
    group_by(.data$sample, .data$taxon, .data$pathway_tpm) |>
    summarise(taxon_tpm = sum(.data$taxon_tpm), .groups = "drop") |>
    mutate(
      value = 100 * .data$taxon_tpm / .data$pathway_tpm,
      denominator = .data$pathway_tpm,
      status = "ok",
      plotted = TRUE
    )

  if (nrow(positive_rows) > 0L) {
    percent_sums <- positive_rows |>
      group_by(.data$sample) |>
      summarise(percent_sum = sum(.data$value), .groups = "drop")
    invalid_sums <- abs(percent_sums$percent_sum - 100) > 1e-6
    if (any(invalid_sums)) {
      stop(
        "Pathway taxonomy percentages do not sum to 100 for samples: ",
        paste(percent_sums$sample[invalid_sums], collapse = ", "),
        call. = FALSE
      )
    }
  }

  zero_samples <- sample_denominators |>
    filter(.data$pathway_tpm <= 0) |>
    pull(.data$sample)
  if (length(zero_samples) > 0L) {
    for (sample_name in zero_samples) {
      warning(
        "Pathway TPM denominator is zero for pathway '", pathway_name,
        "', sample '", sample_name, "', rank '", rank, "'.",
        call. = FALSE
      )
    }
  }

  zero_rows <- tibble(
    sample = zero_samples,
    taxon = NA_character_,
    pathway_tpm = 0,
    taxon_tpm = 0,
    value = NA_real_,
    denominator = 0,
    status = "zero_denominator",
    plotted = FALSE
  )

  bind_rows(positive_rows, zero_rows) |>
    mutate(
      sample = factor(.data$sample, levels = selected_samples),
      count = "percent",
      rank = rank
    ) |>
    arrange(.data$sample, desc(.data$taxon_tpm), .data$taxon)
}

make_pathway_taxonomy_percent_plot <- function(
    plot_tbl,
    pathway_name,
    rank,
    selected_samples) {
  plotted_rows <- plot_tbl |>
    filter(.data$plotted, !is.na(.data$taxon), !is.na(.data$value)) |>
    mutate(sample = factor(as.character(.data$sample), levels = selected_samples))

  taxon_order <- plotted_rows |>
    group_by(.data$taxon) |>
    summarise(total_tpm = sum(.data$taxon_tpm), .groups = "drop") |>
    arrange(desc(.data$total_tpm), .data$taxon) |>
    pull(.data$taxon)
  taxon_palette <- stats::setNames(
    rep(colors_hex, length.out = length(taxon_order)),
    taxon_order
  )

  plot_object <- ggplot(
    plotted_rows,
    aes(x = .data$sample, y = .data$value, fill = .data$taxon)
  ) +
    geom_col() +
    scale_x_discrete(limits = selected_samples, drop = FALSE) +
    labs(
      title = paste0("Taxonomy - ", pathway_name),
      subtitle = paste0("Pathway TPM composition at ", rank, " rank"),
      x = "Sample",
      y = "Percent of pathway TPM",
      fill = rank
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  if (length(taxon_palette) > 0L) {
    plot_object <- plot_object + scale_fill_manual(values = taxon_palette)
  }

  plot_object
}

build_pie_chart_table <- function(
    orf_long,
    sample_name,
    ko_id_filter,
    rank_name,
    top_n_taxa,
    pathway_name = NA_character_,
    pathway_id = NA_character_,
    pathway_selection = NA_character_,
    pathway_sample_tpm = NA_real_) {
  rank_sym <- rlang::sym(rank_name)
  ko_ec_lookup <- extract_ko_ec_lookup(orf_long)
  ko_ec <- ko_ec_lookup |>
    filter(.data$ko_id == ko_id_filter) |>
    pull(.data$ec_codes)
  if (length(ko_ec) == 0L) {
    ko_ec <- NA_character_
  } else {
    ko_ec <- ko_ec[[1L]]
  }
  ko_descriptions <- orf_long |>
    filter(.data$ko_id == ko_id_filter) |>
    pull(.data$kegg_function) |>
    as.character() |>
    trimws()
  ko_descriptions <- sort(unique(ko_descriptions[!is.na(ko_descriptions) & nzchar(ko_descriptions)]))
  ko_name <- if (length(ko_descriptions) == 0L) {
    ko_id_filter
  } else {
    paste(ko_descriptions, collapse = "; ")
  }

  base_tbl <- orf_long |>
    filter(.data$sample == sample_name, .data$ko_id == ko_id_filter) |>
    group_by(taxon_rank = !!rank_sym) |>
    summarise(tpm = sum(.data$tpm), .groups = "drop") |>
    mutate(taxon_rank = replace_na(as.character(.data$taxon_rank), "Unclassified")) |>
    group_by(.data$taxon_rank) |>
    summarise(tpm = sum(.data$tpm), .groups = "drop") |>
    arrange(desc(.data$tpm), .data$taxon_rank)

  if (nrow(base_tbl) == 0L || sum(base_tbl$tpm) <= 0) {
    return(base_tbl)
  }

  top_taxa <- select_top_classified_taxa(
    data = base_tbl,
    taxon_col = "taxon_rank",
    value_col = "tpm",
    top_n = top_n_taxa
  )
  plot_tbl <- base_tbl |>
    mutate(
      taxon_rank = collapse_taxa_preserving_unclassified(
        .data$taxon_rank,
        top_taxa
      )
    ) |>
    group_by(.data$taxon_rank) |>
    summarise(tpm = sum(.data$tpm), .groups = "drop")

  if (!isTRUE(all.equal(sum(plot_tbl$tpm), sum(base_tbl$tpm), tolerance = 1e-10))) {
    stop("Taxonomy Top N collapse changed TPM mass.", call. = FALSE)
  }

  plot_tbl |>
    arrange(desc(.data$tpm), .data$taxon_rank) |>
    mutate(
      sample = sample_name,
      ko_id = ko_id_filter,
      ec_codes = ko_ec,
      total_tpm = sum(.data$tpm),
      pct = if_else(.data$total_tpm > 0, .data$tpm / .data$total_tpm, 0),
      label = if_else(.data$pct >= 0.03, as.character(.data$taxon_rank), ""),
      pathway = pathway_name,
      pathway_id = pathway_id,
      pathway_selection = pathway_selection,
      rank = rank_name,
      ko_name = ko_name,
      ko_sample_tpm = .data$total_tpm,
      pathway_sample_tpm = as.numeric(pathway_sample_tpm),
      ko_pathway_percent = if_else(
        !is.na(.data$pathway_sample_tpm) & .data$pathway_sample_tpm > 0,
        .data$ko_sample_tpm / .data$pathway_sample_tpm * 100,
        NA_real_
      ),
      taxon_order = dplyr::row_number(),
      taxon_rank = factor(
        as.character(.data$taxon_rank),
        levels = unique(as.character(.data$taxon_rank))
      )
    )
}

make_pie_plot <- function(plot_data) {
  required_columns <- c(
    "taxon_rank", "tpm", "sample", "ko_id", "ec_codes", "pct", "label",
    "pathway", "rank", "ko_name", "ko_sample_tpm", "pathway_sample_tpm",
    "ko_pathway_percent", "taxon_order"
  )
  missing_columns <- setdiff(required_columns, colnames(plot_data))
  if (length(missing_columns) > 0L) {
    stop("PIE plot data missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  if (nrow(plot_data) == 0L) {
    stop("PIE plot data cannot be empty.", call. = FALSE)
  }

  scalar_value <- function(column_name, allow_na = FALSE) {
    values <- unique(plot_data[[column_name]])
    if (allow_na && length(values) == 1L && is.na(values[[1L]])) {
      return(values[[1L]])
    }
    if (length(values) != 1L || (!allow_na && is.na(values[[1L]]))) {
      stop("PIE plot metadata must have one value for ", column_name, ".", call. = FALSE)
    }
    values[[1L]]
  }

  pathway_name <- as.character(scalar_value("pathway"))
  sample_name <- as.character(scalar_value("sample"))
  ko_id <- as.character(scalar_value("ko_id"))
  ko_name <- as.character(scalar_value("ko_name"))
  rank_name <- as.character(scalar_value("rank"))
  ko_ec <- scalar_value("ec_codes", allow_na = TRUE)
  total_tpm <- as.numeric(scalar_value("ko_sample_tpm"))
  pathway_contribution <- as.numeric(scalar_value("ko_pathway_percent", allow_na = TRUE))

  if (!isTRUE(all.equal(sum(as.numeric(plot_data$tpm)), total_tpm, tolerance = 1e-10))) {
    stop("PIE exported TPM does not match ko_sample_tpm.", call. = FALSE)
  }
  taxon_order <- as.integer(plot_data$taxon_order)
  if (anyNA(taxon_order) || anyDuplicated(taxon_order) > 0L) {
    stop("PIE taxon_order must contain unique integers.", call. = FALSE)
  }
  ordered_levels <- as.character(plot_data$taxon_rank[order(taxon_order)])
  plot_data <- plot_data |>
    mutate(
      taxon_rank = factor(
        as.character(.data$taxon_rank),
        levels = unique(ordered_levels)
      )
    ) |>
    arrange(.data$taxon_order)

  plot_data_with_legend <- plot_data |>
    mutate(
      legend_label = glue::glue(
        "{taxon_rank} | {format_display_number(tpm)} TPM | {format_display_percent(pct)}"
      )
    )

  taxon_levels <- levels(plot_data_with_legend$taxon_rank)
  non_other <- setdiff(taxon_levels, "Other")
  palette <- c(
    stats::setNames(rep(colors_hex, length.out = length(non_other)), non_other),
    if ("Other" %in% taxon_levels) c(Other = "grey70") else NULL
  )

  legend_map <- plot_data_with_legend |>
    mutate(taxon_rank_chr = as.character(.data$taxon_rank)) |>
    distinct(.data$taxon_rank_chr, .data$legend_label)
  legend_labels <- stats::setNames(legend_map$legend_label, legend_map$taxon_rank_chr)

  ggplot(plot_data_with_legend, aes(x = "", y = .data$tpm, fill = .data$taxon_rank)) +
    geom_col(width = 1, color = "white", linewidth = 0.2) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = .data$label),
      position = position_stack(vjust = 0.5),
      size = 3.2,
      color = "black"
    ) +
    scale_fill_manual(values = palette, labels = legend_labels, drop = FALSE) +
    labs(
      title = as.character(glue::glue("Pathway: {pathway_name} - KO {ko_id}")),
      subtitle = as.character(glue::glue(
        "Sample: {sample_name} | Rank: {rank_name} | total TPM = {format_display_number(total_tpm)} | pathway contribution: {format_display_number(pathway_contribution, suffix = '%')}"
      )),
      fill = rank_name,
      x = NULL,
      y = NULL,
      caption = as.character(glue::glue("KO name: {ko_name}{if (!is.na(ko_ec) && nzchar(ko_ec)) paste0(' | EC: ', ko_ec) else ''}"))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 9, face = "italic"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

# ---- Analysis modes --------------------------------------------------------

# Each mode appends rows to its manifest instead of deleting prior output.
run_funz_mode <- function(
    output_dir,
    manifest_base_dir,
    output_manifests,
    script_name,
    project_dir,
    tax_mode,
    pathway_name,
    pathway_sqm,
    selected_samples,
    dimensions,
    plot_dpi,
    top_n_taxa,
    top_n_ko,
    pathway_id = NA_character_,
    pathway_selection = "defined",
    filtered_taxon = NA_character_,
    filtered_taxon_rank = NA_character_,
    pathway_analysis = NULL) {
  pathway_analysis <- resolve_pathway_analysis(
    pathway_analysis = pathway_analysis,
    pathway_name = pathway_name,
    pathway_sqm = pathway_sqm,
    selected_samples = selected_samples,
    pathway_id = pathway_id,
    pathway_selection = pathway_selection
  )
  pathway_name <- pathway_analysis$pathway_name
  pathway_sqm <- pathway_analysis$pathway_sqm
  selected_samples <- pathway_analysis$selected_samples
  pathway_id <- pathway_analysis$pathway_id
  pathway_selection <- pathway_analysis$pathway_selection
  progress_message(
    "FUNZ | pathway=", pathway_name,
    " | samples=", paste(selected_samples, collapse = ",")
  )
  pathway_dir <- file.path(
    output_dir,
    "funz",
    "pathway",
    pathway_selection_directory(pathway_selection),
    safe_output_component(pathway_name, max_length = 28L)
  )
  dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

  orf_long_result <- pathway_analysis$orf_long_result
  ko_totals <- orf_long_result$totals
  ko_lookup <- pathway_analysis$ko_lookup
  plot_tbl <- build_ko_plot_table(
    orf_long = ko_totals,
    selected_samples = selected_samples,
    top_n_ko = top_n_ko,
    ko_lookup = ko_lookup,
    pathway_name = pathway_name
  )
  plot_path <- file.path(pathway_dir, "barplot_ko_data.tsv")
  progress_message("FUNZ | writing data: ", plot_path)
  plot_path <- write_tsv_safe(plot_tbl, plot_path)

  output_manifests$funz <- bind_rows(
    output_manifests$funz,
    new_manifest_row(
      script_name = script_name,
      project_dir = project_dir,
      tax_mode = tax_mode,
      pathway = pathway_name,
      samples = selected_samples,
      metric = "tpm",
      top_n_taxa = top_n_taxa,
      top_n_ko = top_n_ko,
      output_type = "data_tsv",
      output_file = relative_to_output(plot_path, manifest_base_dir),
      mode = "funz",
      format = "tsv",
      dpi = plot_dpi,
      output_scope = paste0("pathway_", pathway_selection),
      filtered_taxon = filtered_taxon,
      filtered_taxon_rank = filtered_taxon_rank,
      pathway_id = pathway_id,
      ko_audit = orf_long_result$audit
    )
  )

  if (!any(plot_tbl$plotted)) {
    progress_message(
      "FUNZ | skipping PNG because every selected sample has zero denominator | pathway=",
      pathway_name
    )
    return(output_manifests)
  }

  plot_object <- make_ko_barplot(plot_tbl, pathway_name, selected_samples)
  progress_message("FUNZ | saving PNG | pathway=", pathway_name)
  png_files <- save_png_dimensions(plot_object, pathway_dir, "barplot_ko", dimensions, plot_dpi)

  for (dim_name in names(png_files)) {
    dims <- dimensions[[dim_name]]
    output_manifests$funz <- bind_rows(
      output_manifests$funz,
      new_manifest_row(
        script_name = script_name,
        project_dir = project_dir,
        tax_mode = tax_mode,
        pathway = pathway_name,
        samples = selected_samples,
        metric = "tpm",
        top_n_taxa = top_n_taxa,
        top_n_ko = top_n_ko,
        output_type = "plot_png",
        output_file = relative_to_output(png_files[[dim_name]], manifest_base_dir),
        mode = "funz",
        format = "png",
        width = dims[["width"]],
        height = dims[["height"]],
        dpi = plot_dpi,
        output_scope = paste0("pathway_", pathway_selection),
        filtered_taxon = filtered_taxon,
        filtered_taxon_rank = filtered_taxon_rank,
        pathway_id = pathway_id,
        ko_audit = orf_long_result$audit
      )
    )
  }

  output_manifests
}

run_enzyme_mode <- function(
    sqm_object,
    output_dir,
    manifest_base_dir,
    output_manifests,
    script_name,
    project_dir,
    tax_mode,
    selected_samples,
    dimensions,
    plot_dpi,
    top_n_taxa,
    top_n_ko,
    enzyme_ecs = default_enzyme_ecs,
    enzyme_plot_types = default_enzyme_plot_types,
    sample_order_basis = current_sample_order_basis(),
    filtered_taxon = NA_character_,
    filtered_taxon_rank = NA_character_) {
  enzyme_ecs <- normalize_enzyme_ecs(enzyme_ecs)
  enzyme_plot_types <- normalize_enzyme_plot_types(enzyme_plot_types)
  progress_message(
    "ENZIMI | EC=", paste(enzyme_ecs, collapse = ","),
    " | samples=", paste(selected_samples, collapse = ","),
    " | plot_types=", paste(enzyme_plot_types, collapse = ",")
  )

  enzyme_tbl <- build_enzyme_plot_table(sqm_object, selected_samples, enzyme_ecs)
  enzyme_root <- file.path(output_dir, "funz", "enzimi")

  append_manifest_entry <- function(
      manifest_tbl,
      output_file,
      output_type,
      output_scope,
      ec_code = NA_character_,
      width = NA_real_,
      height = NA_real_,
      source_data_file = NA_character_) {
    bind_rows(
      manifest_tbl,
      new_manifest_row(
        script_name = script_name,
        project_dir = project_dir,
        tax_mode = tax_mode,
        pathway = NA_character_,
        ec_code = ec_code,
        samples = selected_samples,
        metric = "tpm",
        top_n_taxa = top_n_taxa,
        top_n_ko = NA_integer_,
        output_type = output_type,
        output_file = relative_to_output(output_file, manifest_base_dir),
        mode = "enzimi",
        format = infer_format_from_path(output_file),
        width = width,
        height = height,
        dpi = plot_dpi,
        output_scope = output_scope,
        filtered_taxon = filtered_taxon,
        filtered_taxon_rank = filtered_taxon_rank,
        ko_selection_policy = "not_applicable",
        source_data_file = source_data_file,
        sample_order_basis = sample_order_basis
      )
    )
  }

  combined_dir <- file.path(enzyme_root, "insieme")
  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  combined_data_file <- file.path(combined_dir, "enzimi_data.tsv")
  combined_data_file <- write_tsv_safe(enzyme_tbl, combined_data_file)
  output_manifests$funz <- append_manifest_entry(
    output_manifests$funz,
    combined_data_file,
    "enzyme_data_tsv",
    "enzyme_insieme"
  )

  combined_data_relative <- relative_to_output(combined_data_file, manifest_base_dir)
  if (!any(enzyme_tbl$plotted)) {
    warning(
      "No requested enzyme has positive TPM; skipping combined enzyme PNG files.",
      call. = FALSE
    )
  }

  if (any(enzyme_tbl$plotted) && "bar" %in% enzyme_plot_types) {
    bar_files <- save_png_dimensions(
      make_enzyme_barplot(enzyme_tbl, "Enzyme barplot - combined"),
      combined_dir,
      "barplot_enzimi",
      dimensions,
      plot_dpi
    )
    for (dim_name in names(bar_files)) {
      dims <- dimensions[[dim_name]]
      output_manifests$funz <- append_manifest_entry(
        output_manifests$funz,
        bar_files[[dim_name]],
        "enzyme_barplot_png",
        "enzyme_insieme",
        width = dims[["width"]],
        height = dims[["height"]],
        source_data_file = combined_data_relative
      )
    }
  }
  if (any(enzyme_tbl$plotted) && "line" %in% enzyme_plot_types) {
    line_files <- save_png_dimensions(
      make_enzyme_lineplot(enzyme_tbl, "Enzyme line chart - combined"),
      combined_dir,
      "lineplot_enzimi",
      dimensions,
      plot_dpi
    )
    for (dim_name in names(line_files)) {
      dims <- dimensions[[dim_name]]
      output_manifests$funz <- append_manifest_entry(
        output_manifests$funz,
        line_files[[dim_name]],
        "enzyme_lineplot_png",
        "enzyme_insieme",
        width = dims[["width"]],
        height = dims[["height"]],
        source_data_file = combined_data_relative
      )
    }
  }

  for (current_ec_code in enzyme_ecs) {
    ec_dir <- file.path(enzyme_root, "separato", sanitize_name(current_ec_code))
    dir.create(ec_dir, recursive = TRUE, showWarnings = FALSE)
    ec_tbl <- enzyme_tbl |>
      filter(as.character(.data$ec_code) == current_ec_code) |>
      mutate(ec_code = factor(as.character(.data$ec_code), levels = current_ec_code))
    ec_data_file <- file.path(ec_dir, "enzima_data.tsv")
    ec_data_file <- write_tsv_safe(ec_tbl, ec_data_file)
    output_manifests$funz <- append_manifest_entry(
      output_manifests$funz,
      ec_data_file,
      "enzyme_data_tsv",
      "enzyme_separato",
      current_ec_code
    )
    ec_data_relative <- relative_to_output(ec_data_file, manifest_base_dir)
    if (!any(ec_tbl$plotted)) {
      warning(
        "Skipping enzyme PNG files for EC ", current_ec_code,
        " because it has no positive TPM.",
        call. = FALSE
      )
      next
    }

    if ("bar" %in% enzyme_plot_types) {
      bar_files <- save_png_dimensions(
        make_enzyme_barplot(ec_tbl, paste0("Enzyme barplot - EC ", current_ec_code)),
        ec_dir,
        "barplot_enzima",
        dimensions,
        plot_dpi
      )
      for (dim_name in names(bar_files)) {
        dims <- dimensions[[dim_name]]
        output_manifests$funz <- append_manifest_entry(
          output_manifests$funz,
          bar_files[[dim_name]],
          "enzyme_barplot_png",
          "enzyme_separato",
          current_ec_code,
          width = dims[["width"]],
          height = dims[["height"]],
          source_data_file = ec_data_relative
        )
      }
    }
    if ("line" %in% enzyme_plot_types) {
      line_files <- save_png_dimensions(
        make_enzyme_lineplot(ec_tbl, paste0("Enzyme line chart - EC ", current_ec_code)),
        ec_dir,
        "lineplot_enzima",
        dimensions,
        plot_dpi
      )
      for (dim_name in names(line_files)) {
        dims <- dimensions[[dim_name]]
        output_manifests$funz <- append_manifest_entry(
          output_manifests$funz,
          line_files[[dim_name]],
          "enzyme_lineplot_png",
          "enzyme_separato",
          current_ec_code,
          width = dims[["width"]],
          height = dims[["height"]],
          source_data_file = ec_data_relative
        )
      }
    }
  }

  output_manifests
}

run_flow_mode <- function(
    output_dir,
    manifest_base_dir,
    output_manifests,
    script_name,
    project_dir,
    tax_mode,
    pathway_name,
    pathway_sqm,
    selected_samples,
    dimensions,
    plot_dpi,
    top_n_taxa,
    top_n_ko,
    taxonomy_ranks,
    flowplot_formats,
    pathway_id = NA_character_,
    pathway_selection = "defined",
    filtered_taxon = NA_character_,
    filtered_taxon_rank = NA_character_,
    pathway_analysis = NULL) {
  pathway_analysis <- resolve_pathway_analysis(
    pathway_analysis = pathway_analysis,
    pathway_name = pathway_name,
    pathway_sqm = pathway_sqm,
    selected_samples = selected_samples,
    pathway_id = pathway_id,
    pathway_selection = pathway_selection
  )
  pathway_name <- pathway_analysis$pathway_name
  pathway_sqm <- pathway_analysis$pathway_sqm
  selected_samples <- pathway_analysis$selected_samples
  pathway_id <- pathway_analysis$pathway_id
  pathway_selection <- pathway_analysis$pathway_selection
  progress_message(
    "FLOW | pathway=", pathway_name,
    " | ranks=", paste(taxonomy_ranks, collapse = ","),
    " | samples=", paste(selected_samples, collapse = ",")
  )
  orf_long_result <- pathway_analysis$orf_long_result
  orf_long <- orf_long_result$data
  if (nrow(orf_long) == 0L) {
    warning("No positive ORF data for pathway ", pathway_name, ".", call. = FALSE)
    return(output_manifests)
  }

  ko_lookup <- pathway_analysis$ko_lookup

  for (rank in taxonomy_ranks) {
    progress_message("FLOW | pathway=", pathway_name, " | rank=", rank)
    rank_dir <- file.path(
      output_dir,
      "flowplot",
      pathway_selection_directory(pathway_selection),
      safe_output_component(pathway_name, max_length = 28L),
      rank
    )
    dir.create(rank_dir, recursive = TRUE, showWarnings = FALSE)
    flow_rank_tbl <- build_flow_table_for_rank(
      orf_long = orf_long,
      rank = rank,
      selected_samples = selected_samples,
      top_n_taxa = top_n_taxa,
      top_n_ko = top_n_ko,
      ko_lookup = ko_lookup
    )

    for (sample_name in selected_samples) {
      progress_message("FLOW | pathway=", pathway_name, " | rank=", rank, " | sample=", sample_name)
      flow_tbl <- build_flow_table_for_sample(flow_rank_tbl, pathway_name, rank, sample_name)
      if (nrow(flow_tbl) == 0L) {
        warning(
          "No positive flow for pathway ", pathway_name,
          ", rank ", rank, ", sample ", sample_name, ".",
          call. = FALSE
        )
        next
      }

      sample_component <- safe_output_component(sample_name)
      data_file <- file.path(
        rank_dir,
        paste0("flowplot_", rank, "_", sample_component, "_data.tsv")
      )
      progress_message("FLOW | writing data: ", data_file)
      data_file <- write_tsv_safe(flow_tbl, data_file)
      output_manifests$flow <- bind_rows(
        output_manifests$flow,
        new_manifest_row(
          script_name = script_name,
          project_dir = project_dir,
          tax_mode = tax_mode,
          pathway = pathway_name,
          samples = sample_name,
          metric = "flow_percent",
          top_n_taxa = top_n_taxa,
          top_n_ko = top_n_ko,
          output_type = "data_tsv",
          output_file = relative_to_output(data_file, manifest_base_dir),
          mode = "flow",
          rank = rank,
          format = "tsv",
          dpi = plot_dpi,
          output_scope = paste0("pathway_", pathway_selection),
          filtered_taxon = filtered_taxon,
          filtered_taxon_rank = filtered_taxon_rank,
          pathway_id = pathway_id,
          ko_audit = orf_long_result$audit
        )
      )

      if ("png" %in% flowplot_formats) {
        progress_message("FLOW | saving PNG | pathway=", pathway_name, " | rank=", rank, " | sample=", sample_name)
        plot_object <- make_flow_plot(flow_tbl, pathway_name, rank, sample_name)
        png_files <- save_png_dimensions(
          plot_object,
          rank_dir,
          paste0("flowplot_", rank, "_", sample_component),
          dimensions,
          plot_dpi
        )

        for (dim_name in names(png_files)) {
          dims <- dimensions[[dim_name]]
          output_manifests$flow <- bind_rows(
            output_manifests$flow,
            new_manifest_row(
              script_name = script_name,
              project_dir = project_dir,
              tax_mode = tax_mode,
              pathway = pathway_name,
              samples = sample_name,
              metric = "flow_percent",
              top_n_taxa = top_n_taxa,
              top_n_ko = top_n_ko,
              output_type = "plot_png",
              output_file = relative_to_output(png_files[[dim_name]], manifest_base_dir),
              mode = "flow",
              rank = rank,
              format = "png",
              width = dims[["width"]],
              height = dims[["height"]],
              dpi = plot_dpi,
              output_scope = paste0("pathway_", pathway_selection),
              filtered_taxon = filtered_taxon,
              filtered_taxon_rank = filtered_taxon_rank,
              pathway_id = pathway_id,
              ko_audit = orf_long_result$audit
            )
          )
        }
      }

      if ("html" %in% flowplot_formats) {
        html_file <- file.path(
          rank_dir,
          paste0("flowplot_", rank, "_", sample_component, ".html")
        )
        progress_message("FLOW | saving HTML | pathway=", pathway_name, " | rank=", rank, " | sample=", sample_name)
        widget <- make_flow_sankey(flow_tbl, pathway_name, rank, sample_name)
        html_file <- save_html_widget(widget, html_file)
        output_manifests$flow <- bind_rows(
          output_manifests$flow,
          new_manifest_row(
            script_name = script_name,
            project_dir = project_dir,
            tax_mode = tax_mode,
            pathway = pathway_name,
            samples = sample_name,
            metric = "flow_percent",
            top_n_taxa = top_n_taxa,
            top_n_ko = top_n_ko,
            output_type = "plot_html",
            output_file = relative_to_output(html_file, manifest_base_dir),
            mode = "flow",
            rank = rank,
            format = "html",
            dpi = plot_dpi,
            output_scope = paste0("pathway_", pathway_selection),
            filtered_taxon = filtered_taxon,
            filtered_taxon_rank = filtered_taxon_rank,
            pathway_id = pathway_id,
            ko_audit = orf_long_result$audit
          )
        )
      }
    }
  }

  output_manifests
}

taxonomy_file_stem <- function(scope_name, count, rank, pathway_name = NA_character_) {
  if (identical(scope_name, "taxonomy_global")) {
    return(paste0("taxonomy_global_", count, "_", rank))
  }
  if (identical(scope_name, "taxonomy_by_pathway")) {
    return(paste0("taxonomy_", count, "_", rank))
  }
  stop("Unsupported taxonomy scope: ", scope_name, call. = FALSE)
}

run_taxonomy_scope <- function(
    sqm_object,
    output_dir,
    manifest_base_dir,
    output_manifests,
    script_name,
    project_dir,
    tax_mode,
    pathway_name,
    selected_samples,
    dimensions,
    plot_dpi,
    top_n_taxa,
    top_n_ko,
    taxonomy_ranks,
    taxonomy_counts,
    scope_name,
    ignore_unmapped,
    ignore_unclassified,
    pathway_id = NA_character_,
    pathway_selection = "defined",
    filtered_taxon = NA_character_,
    filtered_taxon_rank = NA_character_) {
  scope_label <- if (is.na(pathway_name)) "global" else "pathway"
  progress_message(
    "TAXON | scope=", scope_name,
    if (!is.na(pathway_name)) paste0(" | pathway=", pathway_name) else "",
    " | counts=", paste(taxonomy_counts, collapse = ","),
    " | ranks=", paste(taxonomy_ranks, collapse = ",")
  )

  for (count in taxonomy_counts) {
    for (rank in taxonomy_ranks) {
      progress_message(
        "TAXON | scope=", scope_name,
        if (!is.na(pathway_name)) paste0(" | pathway=", pathway_name) else "",
        " | count=", count,
        " | rank=", rank
      )
      rank_dir <- if (scope_name == "taxonomy_global") {
        file.path(output_dir, scope_name, count, rank)
      } else {
      file.path(
        output_dir,
        scope_name,
        pathway_selection_directory(pathway_selection),
        safe_output_component(pathway_name, max_length = 28L),
        count,
        rank
      )
      }
      dir.create(rank_dir, recursive = TRUE, showWarnings = FALSE)

      plot_object <- make_taxonomy_plot(
        sqm_object = sqm_object,
        rank = rank,
        count = count,
        selected_samples = selected_samples,
        top_n_taxa = top_n_taxa,
        ignore_unmapped = ignore_unmapped,
        ignore_unclassified = ignore_unclassified
      )
      plot_data <- extract_taxonomy_plot_data(plot_object, count)
      if (scope_name == "taxonomy_by_pathway") {
        plot_data <- plot_data |>
          select(all_of(c("sample", "taxon", "value", "count")))
      }
      if (scope_name == "taxonomy_global" && identical(count, "percent")) {
        taxonomy_context_label <- if (is.na(filtered_taxon)) {
          "global"
        } else {
          paste0(filtered_taxon, "@", filtered_taxon_rank)
        }
        plot_data <- add_global_taxonomy_percent_metadata(
          plot_data = plot_data,
          sqm_object = sqm_object,
          rank = rank,
          selected_samples = selected_samples,
          context_label = taxonomy_context_label
        )
      }
      file_stem <- taxonomy_file_stem(scope_name, count, rank, pathway_name)
      data_file <- file.path(rank_dir, paste0(file_stem, "_data.tsv"))
      progress_message("TAXON | writing data: ", data_file)
      data_file <- write_tsv_safe(plot_data, data_file)

      output_manifests$taxon <- bind_rows(
        output_manifests$taxon,
        new_manifest_row(
          script_name = script_name,
          project_dir = project_dir,
          tax_mode = tax_mode,
          pathway = ifelse(is.na(pathway_name), NA_character_, pathway_name),
          samples = selected_samples,
          metric = count,
          top_n_taxa = top_n_taxa,
          top_n_ko = top_n_ko,
          output_type = "data_tsv",
          output_file = relative_to_output(data_file, manifest_base_dir),
          mode = "taxon",
          rank = rank,
          count = count,
          format = "tsv",
          dpi = plot_dpi,
        output_scope = if (
          is.na(pathway_name)
        ) scope_label else paste0(scope_label, "_", pathway_selection),
          filtered_taxon = filtered_taxon,
          filtered_taxon_rank = filtered_taxon_rank,
          pathway_id = pathway_id,
          ko_selection_policy = "not_applicable",
          taxonomy_display_policy = if (
            scope_name == "taxonomy_global"
          ) "SQMtools_non_rescaled_excluding_unmapped_unclassified" else NA_character_,
          denominator_type = if (
            scope_name == "taxonomy_global" && identical(count, "percent")
          ) "total_reads" else NA_character_
        )
      )

      progress_message(
        "TAXON | saving PNG | scope=", scope_name,
        if (!is.na(pathway_name)) paste0(" | pathway=", pathway_name) else "",
        " | count=", count,
        " | rank=", rank
      )
      png_files <- save_png_dimensions(plot_object, rank_dir, file_stem, dimensions, plot_dpi)

      for (dim_name in names(png_files)) {
        dims <- dimensions[[dim_name]]
        output_manifests$taxon <- bind_rows(
          output_manifests$taxon,
          new_manifest_row(
            script_name = script_name,
            project_dir = project_dir,
            tax_mode = tax_mode,
            pathway = ifelse(is.na(pathway_name), NA_character_, pathway_name),
            samples = selected_samples,
            metric = count,
            top_n_taxa = top_n_taxa,
            top_n_ko = top_n_ko,
            output_type = "plot_png",
            output_file = relative_to_output(png_files[[dim_name]], manifest_base_dir),
            mode = "taxon",
            rank = rank,
            count = count,
            format = "png",
            width = dims[["width"]],
            height = dims[["height"]],
            dpi = plot_dpi,
          output_scope = if (
            is.na(pathway_name)
          ) scope_label else paste0(scope_label, "_", pathway_selection),
            filtered_taxon = filtered_taxon,
            filtered_taxon_rank = filtered_taxon_rank,
            pathway_id = pathway_id,
            ko_selection_policy = "not_applicable",
            taxonomy_display_policy = if (
              scope_name == "taxonomy_global"
            ) "SQMtools_non_rescaled_excluding_unmapped_unclassified" else NA_character_,
            denominator_type = if (
              scope_name == "taxonomy_global" && identical(count, "percent")
            ) "total_reads" else NA_character_,
            source_data_file = relative_to_output(data_file, manifest_base_dir)
          )
        )
      }
    }
  }

  output_manifests
}

# ---- Manifest compatibility ------------------------------------------------

# Keep legacy pie manifests readable while normalizing new entries.
infer_format_from_path <- function(path) {
  ext <- tools::file_ext(path)
  if (!nzchar(ext)) {
    return(NA_character_)
  }
  tolower(ext)
}

build_pathview_input_table <- function(sqm_object, selected_samples) {
  kegg_tpm <- sqm_object$functions$KEGG$tpm
  if (is.null(kegg_tpm)) {
    stop("sqm$functions$KEGG$tpm is required for Pathview.", call. = FALSE)
  }
  kegg_frame <- as.data.frame(kegg_tpm, check.names = FALSE)
  validate_samples(selected_samples, colnames(kegg_frame))
  validate_tpm_matrix(kegg_frame[selected_samples], "sqm$functions$KEGG$tpm")
  ko_ids <- rownames(kegg_frame)
  if (is.null(ko_ids) || anyDuplicated(ko_ids) > 0L || any(!nzchar(ko_ids))) {
    stop("Pathview KEGG TPM rows require unique, non-empty KO IDs.", call. = FALSE)
  }
  kegg_frame |>
    tibble::rownames_to_column("ko_id") |>
    select("ko_id", all_of(selected_samples))
}

export_pathview_isolated <- function(
    export_pathway_fn,
    sqm_object,
    pathway_id,
    selected_samples,
    final_dir) {
  temporary_dir <- tempfile("sqm_pathview_export_")
  dir.create(temporary_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temporary_dir, recursive = TRUE, force = TRUE), add = TRUE)
  temporary_dir_abs <- normalizePath(temporary_dir, winslash = "/", mustWork = TRUE)

  export_pathway_fn(
    SQM = sqm_object,
    pathway_id = pathway_id,
    count = "tpm",
    samples = selected_samples,
    split_samples = FALSE,
    log_scale = FALSE,
    output_dir = temporary_dir_abs,
    output_suffix = paste0("pathview_", pathway_id)
  )

  produced_files <- list.files(
    temporary_dir_abs,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE,
    include.dirs = FALSE
  )
  if (length(produced_files) == 0L) {
    stop("Pathview did not produce any output files.", call. = FALSE)
  }

  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  relative_files <- substring(produced_files, nchar(temporary_dir_abs) + 2L)
  final_files <- file.path(final_dir, relative_files)
  for (file_index in seq_along(produced_files)) {
    dir.create(dirname(final_files[[file_index]]), recursive = TRUE, showWarnings = FALSE)
    copied <- file.copy(
      produced_files[[file_index]],
      final_files[[file_index]],
      overwrite = TRUE,
      copy.mode = TRUE,
      copy.date = TRUE
    )
    if (!isTRUE(copied)) {
      stop("Failed to transfer Pathview artifact: ", relative_files[[file_index]], call. = FALSE)
    }
    assert_output_artifact(final_files[[file_index]], "Pathview output")
  }
  sort(final_files)
}

run_pathview_mode <- function(
    sqm_object,
    output_dir,
    manifest_base_dir,
    output_manifests,
    script_name,
    project_dir,
    tax_mode,
    pathway_name,
    pathway_id,
    selected_samples,
    top_n_taxa,
    top_n_ko,
    pathway_selection = "defined",
    pathview_sample_modes = c("insieme", "separato"),
    filtered_taxon = NA_character_,
    filtered_taxon_rank = NA_character_,
    export_pathway_fn = SQMtools::exportPathway) {
  pathview_sample_modes <- normalize_pathview_sample_modes(pathview_sample_modes)
  progress_message(
    "PATHVIEW | pathway=", pathway_name,
    " | pathway_id=", pathway_id,
    " | samples=", paste(selected_samples, collapse = ","),
    " | sample_modes=", paste(pathview_sample_modes, collapse = ",")
  )
  if (is.na(pathway_id) || !nzchar(pathway_id)) {
    stop(
      "Missing pathway_id for pathway '", pathway_name,
      "'. Pathview requires a known numeric code.",
      call. = FALSE
    )
  }

  append_pathview_row <- function(
      manifest,
      invocation_samples,
      output_type,
      output_file,
      output_scope,
      source_data_file = NA_character_,
      ko_selection_policy = "pathview_native_mapping") {
    bind_rows(
      manifest,
      new_manifest_row(
        script_name = script_name,
        project_dir = project_dir,
        tax_mode = tax_mode,
        pathway = pathway_name,
        pathway_id = pathway_id,
        samples = invocation_samples,
        metric = "tpm",
        top_n_taxa = top_n_taxa,
        top_n_ko = NA_integer_,
        output_type = output_type,
        output_file = relative_to_output(output_file, manifest_base_dir),
        mode = "pathview",
        format = infer_format_from_path(output_file),
        output_scope = output_scope,
        filtered_taxon = filtered_taxon,
        filtered_taxon_rank = filtered_taxon_rank,
        ko_selection_policy = ko_selection_policy,
        source_data_file = source_data_file
      )
    )
  }

  for (pathview_sample_mode in pathview_sample_modes) {
    invocation_samples <- if (identical(pathview_sample_mode, "insieme")) {
      list(selected_samples)
    } else {
      lapply(selected_samples, function(sample_name) sample_name)
    }

    for (current_samples in invocation_samples) {
      pathview_dir <- file.path(
        output_dir,
        "pathview",
        pathway_selection_directory(pathway_selection),
        pathview_sample_mode,
        safe_output_component(pathway_name, max_length = 28L)
      )
      if (identical(pathview_sample_mode, "separato")) {
        pathview_dir <- file.path(pathview_dir, safe_output_component(current_samples[[1L]]))
      }
      dir.create(pathview_dir, recursive = TRUE, showWarnings = FALSE)
      output_scope <- paste0(
        "pathway_", pathway_selection, "_", pathview_sample_mode
      )

      input_table <- build_pathview_input_table(sqm_object, current_samples)
      input_file <- file.path(pathview_dir, "pathview_input_all_ko_complete_matrix.tsv")
      config_file <- file.path(pathview_dir, "pathview_render_config.tsv")
      input_file <- write_tsv_safe(input_table, input_file)
      config_file <- write_tsv_safe(
        tibble::tibble(
          pathway_id = pathway_id,
          metric = "tpm",
          samples = paste(current_samples, collapse = ","),
          pathview_sample_mode = pathview_sample_mode,
          split_samples = FALSE,
          log_scale = FALSE,
          pseudocount = NA_real_,
          color_bins = 10L,
          max_scale_value = "automatic",
          color_source = "pathview_native",
          input_scope = "complete_all_ko_matrix"
        ),
        config_file
      )
      input_relative <- relative_to_output(input_file, manifest_base_dir)
      output_manifests$pathview <- append_pathview_row(
        output_manifests$pathview,
        current_samples,
        "pathview_input_all_ko_complete_matrix_tsv",
        input_file,
        output_scope,
        ko_selection_policy = "all_ko_complete_matrix"
      )
      output_manifests$pathview <- append_pathview_row(
        output_manifests$pathview,
        current_samples,
        "pathview_render_config_tsv",
        config_file,
        output_scope,
        source_data_file = input_relative,
        ko_selection_policy = "all_ko_complete_matrix"
      )

      progress_message(
        "PATHVIEW | mode=", pathview_sample_mode,
        " | samples=", paste(current_samples, collapse = ","),
        " | isolated export"
      )
      produced_files <- tryCatch(
        export_pathview_isolated(
          export_pathway_fn = export_pathway_fn,
          sqm_object = sqm_object,
          pathway_id = pathway_id,
          selected_samples = current_samples,
          final_dir = pathview_dir
        ),
        error = function(error) {
          warning(
            "Skipping PATHVIEW export for pathway ", pathway_id,
            " (mode=", pathview_sample_mode,
            ", samples=", paste(current_samples, collapse = ","),
            "): ", conditionMessage(error),
            call. = FALSE
          )
          character()
        }
      )
      for (output_file in produced_files) {
        output_manifests$pathview <- append_pathview_row(
          output_manifests$pathview,
          current_samples,
          "pathview_file",
          output_file,
          output_scope,
          source_data_file = input_relative
        )
      }
    }
  }

  output_manifests
}

run_pie_mode <- function(
    output_dir,
    manifest_base_dir,
    output_manifests,
    script_name,
    project_dir,
    tax_mode,
    pathway_name,
    pathway_sqm,
    selected_samples,
    taxonomy_ranks,
    dimensions,
    plot_dpi,
    top_n_taxa,
    top_n_ko,
    pathway_id = NA_character_,
    pathway_selection = "defined",
    filtered_taxon = NA_character_,
    filtered_taxon_rank = NA_character_,
    pathway_analysis = NULL) {
  pathway_analysis <- resolve_pathway_analysis(
    pathway_analysis = pathway_analysis,
    pathway_name = pathway_name,
    pathway_sqm = pathway_sqm,
    selected_samples = selected_samples,
    pathway_id = pathway_id,
    pathway_selection = pathway_selection
  )
  pathway_name <- pathway_analysis$pathway_name
  pathway_sqm <- pathway_analysis$pathway_sqm
  selected_samples <- pathway_analysis$selected_samples
  pathway_id <- pathway_analysis$pathway_id
  pathway_selection <- pathway_analysis$pathway_selection
  progress_message(
    "PIE | pathway=", pathway_name,
    " | ranks=", paste(taxonomy_ranks, collapse = ","),
    " | samples=", paste(selected_samples, collapse = ",")
  )
  pie_root <- file.path(
    output_dir,
    "pie",
    pathway_selection_directory(pathway_selection),
    safe_output_component(pathway_name, max_length = 28L)
  )
  dir.create(pie_root, recursive = TRUE, showWarnings = FALSE)

  orf_long_result <- pathway_analysis$orf_long_result
  orf_long <- orf_long_result$data
  if (nrow(orf_long) == 0L) {
    warning("No positive ORF data for pathway ", pathway_name, ".", call. = FALSE)
    return(output_manifests)
  }

  ko_ec_lookup <- orf_long_result$metadata |>
    transmute(ko_id = .data$ko_id, ec_codes = .data$ec_codes)
  pathway_sample_totals <- orf_long_result$totals |>
    group_by(.data$sample) |>
    summarise(pathway_sample_tpm = sum(.data$tpm), .groups = "drop")

  positive_ko_totals <- orf_long_result$totals |>
    filter(.data$tpm > 0)
  positive_sample_names <- positive_ko_totals |>
    distinct(.data$sample) |>
    pull(.data$sample)
  sample_names <- selected_samples[selected_samples %in% positive_sample_names]
  ko_ids <- positive_ko_totals |>
    distinct(.data$ko_id) |>
    pull(.data$ko_id) |>
    sort()

  for (sample_name in sample_names) {
    progress_message("PIE | pathway=", pathway_name, " | sample=", sample_name)
    sample_dir <- file.path(pie_root, safe_output_component(sample_name))
    for (ko_id_value in ko_ids) {
      ko_ec_row <- ko_ec_lookup |>
        filter(.data$ko_id == ko_id_value)
      ko_ec <- if (nrow(ko_ec_row) > 0L) ko_ec_row$ec_codes[[1]] else NA_character_
      ko_dir <- file.path(
        sample_dir,
        safe_output_component(get_ko_dir_name(ko_id_value, ko_ec), max_length = 20L)
      )
      dir.create(ko_dir, recursive = TRUE, showWarnings = FALSE)

      pathway_sample_tpm <- pathway_sample_totals |>
        filter(.data$sample == sample_name) |>
        pull(.data$pathway_sample_tpm)
      if (length(pathway_sample_tpm) == 0L) {
        pathway_sample_tpm <- 0
      }

      for (rank_name in taxonomy_ranks) {
        progress_message("PIE | pathway=", pathway_name, " | sample=", sample_name, " | KO=", ko_id_value, " | rank=", rank_name)
        plot_data <- build_pie_chart_table(
          orf_long = orf_long,
          sample_name = sample_name,
          ko_id_filter = ko_id_value,
          rank_name = rank_name,
          top_n_taxa = top_n_taxa,
          pathway_name = pathway_name,
          pathway_id = pathway_id,
          pathway_selection = pathway_selection,
          pathway_sample_tpm = pathway_sample_tpm
        )

        if (nrow(plot_data) == 0L || sum(plot_data$tpm) <= 0) {
          next
        }

        data_file_stem <- paste0(
          rank_name, "_", sanitize_name(ko_id_value), "_", sanitize_name(coalesce(ko_ec, "NA")), "_", sanitize_name(sample_name)
        )
        data_file <- file.path(ko_dir, paste0(data_file_stem, "_data.tsv"))
        progress_message("PIE | writing data: ", data_file)
        data_file <- write_tsv_safe(plot_data, data_file)

        output_manifests$pie <- bind_rows(
          output_manifests$pie,
          new_manifest_row(
            script_name = script_name,
            project_dir = project_dir,
          tax_mode = tax_mode,
          pathway = pathway_name,
          ec_code = ko_ec,
          ko_id = ko_id_value,
          pathway_id = pathway_id,
            samples = sample_name,
            metric = "tpm",
            top_n_taxa = top_n_taxa,
            top_n_ko = NA_integer_,
            output_type = "data_tsv",
            output_file = relative_to_output(data_file, manifest_base_dir),
            mode = "pie",
            rank = rank_name,
            format = "tsv",
            dpi = plot_dpi,
            output_scope = paste0("pathway_", pathway_selection),
            filtered_taxon = filtered_taxon,
            filtered_taxon_rank = filtered_taxon_rank,
            ko_selection_policy = "all_positive_ko",
            ko_audit = orf_long_result$audit
          )
        )

        plot_object <- make_pie_plot(plot_data)
        progress_message("PIE | saving PNG | pathway=", pathway_name, " | sample=", sample_name, " | KO=", ko_id_value, " | rank=", rank_name)
        png_files <- save_png_dimensions(plot_object, ko_dir, data_file_stem, dimensions, plot_dpi)

        for (dim_name in names(png_files)) {
          dims <- dimensions[[dim_name]]
          output_manifests$pie <- bind_rows(
            output_manifests$pie,
            new_manifest_row(
              script_name = script_name,
              project_dir = project_dir,
            tax_mode = tax_mode,
            pathway = pathway_name,
            ec_code = ko_ec,
            ko_id = ko_id_value,
            pathway_id = pathway_id,
              samples = sample_name,
              metric = "tpm",
              top_n_taxa = top_n_taxa,
              top_n_ko = NA_integer_,
              output_type = "plot_png",
              output_file = relative_to_output(png_files[[dim_name]], manifest_base_dir),
              mode = "pie",
              rank = rank_name,
              format = "png",
              width = dims[["width"]],
              height = dims[["height"]],
              dpi = plot_dpi,
              output_scope = paste0("pathway_", pathway_selection),
              filtered_taxon = filtered_taxon,
              filtered_taxon_rank = filtered_taxon_rank,
              ko_selection_policy = "all_positive_ko",
              ko_audit = orf_long_result$audit
            )
          )
        }
      }
    }
  }

  output_manifests
}

normalize_legacy_pie_manifest <- function(manifest_tbl) {
  if (nrow(manifest_tbl) == 0L || !"output_file" %in% colnames(manifest_tbl)) {
    return(manifest_tbl)
  }

  output_file <- as.character(manifest_tbl$output_file)
  legacy_pie <- !is.na(output_file) &
    startsWith(output_file, "pie/") &
    !startsWith(output_file, "pie/definiti/") &
    !startsWith(output_file, "pie/top20/")

  if (!any(legacy_pie)) {
    return(manifest_tbl)
  }

  manifest_tbl |>
    mutate(
      output_file = if_else(
        legacy_pie,
        paste0("pie/definiti/", sub("^pie/", "", output_file)),
        output_file
      ),
      output_scope = if_else(legacy_pie, "pathway_defined", as.character(.data$output_scope))
    )
}

new_context_manifest_registry <- function(context_output_dir) {
  registry <- new.env(parent = emptyenv())
  registry$context_output_dir <- context_output_dir
  registry$flow <- tibble::tibble()
  registry$funz <- tibble::tibble()
  registry$taxon <- tibble::tibble()
  registry$pathview <- tibble::tibble()
  registry$pie <- tibble::tibble()
  registry
}

activate_context_manifest_registry <- function(registry) {
  if (!is.environment(registry) || is.null(registry$context_output_dir)) {
    stop("Invalid context manifest registry.", call. = FALSE)
  }
  .sqm_run_context$active_manifest_registry <- registry
  invisible(registry)
}

clear_active_manifest_registry <- function() {
  .sqm_run_context$active_manifest_registry <- NULL
  invisible(TRUE)
}

read_section_manifest <- function(manifest_path, relative_dir) {
  if (!file.exists(manifest_path)) {
    return(tibble::tibble())
  }

  manifest_tbl <- readr::read_tsv(manifest_path, show_col_types = FALSE, na = "NA")
  if (identical(relative_dir, "pie")) {
    normalize_legacy_pie_manifest(manifest_tbl)
  } else {
    manifest_tbl
  }
}

manifest_target_status <- function(output_files, output_dir) {
  output_files <- as.character(output_files)
  output_root <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  output_root_prefix <- paste0(output_root, "/")
  case_normalize <- if (.Platform$OS.type == "windows") tolower else identity

  inspected <- lapply(output_files, function(output_file) {
    normalized_relative <- if (is.na(output_file)) NA_character_ else gsub("\\\\", "/", output_file)
    reason <- NA_character_
    target_path <- NA_character_

    if (is.na(normalized_relative) || !nzchar(normalized_relative)) {
      reason <- "empty output_file"
    } else if (grepl("^(?:[A-Za-z]:[/\\\\]|[/\\\\]{1,2})", normalized_relative, perl = TRUE)) {
      reason <- "output_file must be relative"
    } else if (".." %in% strsplit(normalized_relative, "/+", perl = TRUE)[[1L]]) {
      reason <- "output_file must not contain '..' path traversal"
    } else {
      target_path <- normalizePath(
        file.path(output_root, normalized_relative),
        winslash = "/",
        mustWork = FALSE
      )
      inside_root <- startsWith(
        case_normalize(target_path),
        case_normalize(output_root_prefix)
      )
      if (!inside_root) {
        reason <- "output_file resolves outside output_dir"
      } else if (!file.exists(target_path)) {
        reason <- "target does not exist"
      } else {
        target_info <- suppressWarnings(file.info(target_path))
        if (isTRUE(target_info$isdir[[1L]])) {
          reason <- "target is not a regular file"
        } else if (is.na(target_info$size[[1L]]) || target_info$size[[1L]] <= 0) {
          reason <- "target is empty"
        }
      }
    }

    tibble::tibble(
      output_file = normalized_relative,
      target_path = target_path,
      valid = is.na(reason),
      reason = reason
    )
  })

  dplyr::bind_rows(inspected)
}

validate_current_manifest_targets <- function(manifest_tbl, output_dir) {
  if (nrow(manifest_tbl) == 0L) {
    return(invisible(manifest_tbl))
  }
  if (!"output_file" %in% colnames(manifest_tbl)) {
    stop("Current-run manifest rows require an output_file column.", call. = FALSE)
  }

  status <- manifest_target_status(manifest_tbl$output_file, output_dir)
  invalid <- status[!status$valid, , drop = FALSE]
  if (nrow(invalid) > 0L) {
    details <- paste0(invalid$output_file, " (", invalid$reason, ")")
    stop(
      "Current-run manifest target validation failed: ",
      paste(details, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(manifest_tbl)
}

prune_stale_manifest_targets <- function(manifest_tbl, output_dir, section_label) {
  if (nrow(manifest_tbl) == 0L) {
    return(manifest_tbl)
  }
  if (!"output_file" %in% colnames(manifest_tbl)) {
    warning(
      "Pruned ", nrow(manifest_tbl), " stale ", section_label,
      " manifest row(s): missing output_file column.",
      call. = FALSE
    )
    return(manifest_tbl[0, , drop = FALSE])
  }

  status <- manifest_target_status(manifest_tbl$output_file, output_dir)
  if (any(!status$valid)) {
    warning(
      "Pruned ", sum(!status$valid), " stale ", section_label,
      " manifest row(s): ",
      paste0(status$output_file[!status$valid], " (", status$reason[!status$valid], ")", collapse = ", "),
      call. = FALSE
    )
  }
  kept <- manifest_tbl[status$valid, , drop = FALSE]
  if (nrow(kept) > 0L) {
    kept$output_file <- status$output_file[status$valid]
  }
  kept
}

merge_section_manifest <- function(new_manifest_tbl, existing_manifest_tbl) {
  if (nrow(new_manifest_tbl) > 0L && "output_file" %in% colnames(new_manifest_tbl)) {
    new_manifest_tbl$output_file <- gsub("\\\\", "/", as.character(new_manifest_tbl$output_file))
  }
  common_columns <- intersect(colnames(new_manifest_tbl), colnames(existing_manifest_tbl))
  for (column_name in common_columns) {
    new_column <- new_manifest_tbl[[column_name]]
    existing_column <- existing_manifest_tbl[[column_name]]
    if (is.logical(existing_column) && all(is.na(existing_column)) && !is.logical(new_column)) {
      existing_manifest_tbl[[column_name]] <- rep(new_column[NA_integer_], length(existing_column))
    } else if (is.logical(new_column) && all(is.na(new_column)) && !is.logical(existing_column)) {
      new_manifest_tbl[[column_name]] <- rep(existing_column[NA_integer_], length(new_column))
    }
  }
  merged_manifest_tbl <- bind_rows(new_manifest_tbl, existing_manifest_tbl)
  if (!"output_file" %in% colnames(merged_manifest_tbl)) {
    return(merged_manifest_tbl[0, , drop = FALSE])
  }
  merged_manifest_tbl |>
    distinct(.data$output_file, .keep_all = TRUE)
}

write_tsv_atomic <- function(
    data,
    path,
    readr_available = requireNamespace("readr", quietly = TRUE)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  backup_path <- tempfile(
    pattern = paste0(".", basename(path), ".backup."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(c(temporary_path, backup_path), force = TRUE), add = TRUE)
  write_tsv_contents(data, temporary_path, readr_available = readr_available)
  temporary_info <- suppressWarnings(file.info(temporary_path))
  if (!file.exists(temporary_path) || is.na(temporary_info$size[[1L]]) || temporary_info$size[[1L]] <= 0) {
    stop("Atomic TSV staging did not create a non-empty file: ", path, call. = FALSE)
  }

  had_existing <- file.exists(path)
  if (had_existing && !file.rename(path, backup_path)) {
    stop("Unable to stage the existing TSV for replacement: ", path, call. = FALSE)
  }
  replaced <- file.rename(temporary_path, path)
  if (!isTRUE(replaced)) {
    if (had_existing) {
      file.rename(backup_path, path)
    }
    stop("Unable to replace TSV atomically: ", path, call. = FALSE)
  }
  if (had_existing) {
    unlink(backup_path, force = TRUE)
  }
  path
}

write_section_manifest <- function(manifest_tbl, output_dir, relative_dir, filename) {
  dir_path <- file.path(output_dir, relative_dir)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  validate_current_manifest_targets(manifest_tbl, dir_path)
  manifest_path <- file.path(dir_path, filename)
  write_tsv_atomic(manifest_tbl, manifest_path)
  record_completed_manifest(manifest_path)
  append_run_log(
    "MANIFEST",
    paste0("Wrote ", nrow(manifest_tbl), " row(s): ", normalizePath(manifest_path, winslash = "/", mustWork = TRUE))
  )
  manifest_path
}

flush_context_manifests <- function(registry) {
  if (!is.environment(registry) || is.null(registry$context_output_dir)) {
    stop("Invalid context manifest registry.", call. = FALSE)
  }
  section_contract <- list(
    flow = c(directory = "flowplot", filename = "manifest_flow.tsv"),
    funz = c(directory = "funz", filename = "manifest_funz.tsv"),
    pie = c(directory = "pie", filename = "manifest_pie.tsv")
  )
  written <- character()
  for (section in names(section_contract)) {
    manifest_tbl <- registry[[section]]
    if (is.null(manifest_tbl) || nrow(manifest_tbl) == 0L) {
      next
    }
    contract <- section_contract[[section]]
    written[[section]] <- write_section_manifest(
      manifest_tbl = manifest_tbl,
      output_dir = registry$context_output_dir,
      relative_dir = unname(contract[["directory"]]),
      filename = unname(contract[["filename"]])
    )
  }
  written
}

flush_active_context_manifests <- function() {
  registry <- .sqm_run_context$active_manifest_registry
  if (is.null(registry)) {
    return(character())
  }
  flush_context_manifests(registry)
}

write_combined_manifest <- function(
    output_dir,
    section_manifest_paths = character()) {
  section_names <- names(section_manifest_paths)
  if (is.null(section_names)) {
    section_names <- character()
  }
  if (length(section_manifest_paths) > 0L) {
    if (any(!nzchar(section_names))) {
      stop("section_manifest_paths must be named by section.", call. = FALSE)
    }
    valid_manifest <- vapply(
      section_manifest_paths,
      function(path) file.exists(path) && !dir.exists(path) && file.info(path)$size > 0,
      logical(1)
    )
    if (!all(valid_manifest)) {
      stop(
        "Combined manifest received missing or empty section manifests: ",
        paste(names(section_manifest_paths)[!valid_manifest], collapse = ", "),
        call. = FALSE
      )
    }
  }

  manifest_all <- tibble::tibble(
    run_id = rep(current_run_id(), length(section_manifest_paths)),
    section = section_names,
    manifest_file = vapply(
      section_manifest_paths,
      relative_to_output,
      character(1),
      output_dir = output_dir
    )
  )
  manifest_all_path <- file.path(output_dir, "manifest_all.tsv")
  write_tsv_safe(manifest_all, manifest_all_path)
}

format_run_time <- function(value) {
  if (length(value) != 1L || is.na(value)) {
    return(NA_character_)
  }
  format(value, "%Y-%m-%dT%H:%M:%S%z")
}

serialize_pathway_skips <- function(skips) {
  if (is.null(skips) || nrow(skips) == 0L) {
    return(NA_character_)
  }
  paste(
    paste0(
      "context=", skips$context,
      ";pathway=", skips$pathway,
      ";reason=", skips$reason
    ),
    collapse = " | "
  )
}

list_run_artifacts <- function(output_dir, run_id = current_run_id(required = TRUE)) {
  if (!dir.exists(output_dir)) {
    return(data.frame(
      run_id = character(),
      output_file = character(),
      size_bytes = numeric(),
      modified_at = character(),
      stringsAsFactors = FALSE
    ))
  }
  paths <- list.files(
    output_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE
  )
  paths <- paths[grepl(paste0("__", run_id), basename(paths), fixed = TRUE)]
  if (length(paths) == 0L) {
    return(data.frame(
      run_id = character(),
      output_file = character(),
      size_bytes = numeric(),
      modified_at = character(),
      stringsAsFactors = FALSE
    ))
  }
  info <- file.info(paths)
  keep <- !info$isdir & !is.na(info$size) & info$size > 0
  paths <- paths[keep]
  info <- info[keep, , drop = FALSE]
  artifacts <- data.frame(
    run_id = rep(run_id, length(paths)),
    output_file = vapply(paths, relative_to_output, character(1), output_dir = output_dir),
    size_bytes = as.numeric(info$size),
    modified_at = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  artifacts[order(artifacts$output_file), , drop = FALSE]
}

write_run_manifest <- function(
    output_dir,
    project_dir,
    mode,
    tax_mode,
    status,
    error_message = NA_character_,
    artifact_manifest = NA_character_,
    combined_manifest = NA_character_) {
  run_id <- current_run_id(required = TRUE)
  warnings <- .sqm_run_context$warnings
  if (is.null(warnings)) warnings <- character()
  skips <- .sqm_run_context$pathway_skips
  if (is.null(skips)) {
    skips <- data.frame(
      context = character(),
      pathway = character(),
      reason = character(),
      stringsAsFactors = FALSE
    )
  }
  relative_manifest <- function(path) {
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
      return(NA_character_)
    }
    relative_to_output(path, output_dir)
  }
  run_manifest <- data.frame(
    run_id = run_id,
    status = status,
    started_at = format_run_time(.sqm_run_context$started_at),
    finished_at = format_run_time(Sys.time()),
    project_dir = as.character(project_dir),
    output_dir = as.character(output_dir),
    tax_mode = as.character(tax_mode),
    mode = as.character(mode),
    samples = paste(.sqm_run_context$samples, collapse = ","),
    sample_order_basis = current_sample_order_basis(),
    cli_arguments = paste(.sqm_run_context$cli_args, collapse = " "),
    warning_count = length(warnings),
    warnings = if (length(warnings) == 0L) NA_character_ else paste(warnings, collapse = " | "),
    skipped_pathway_count = nrow(skips),
    skipped_pathways = serialize_pathway_skips(skips),
    error_message = error_message,
    artifact_manifest = relative_manifest(artifact_manifest),
    combined_manifest = relative_manifest(combined_manifest),
    stringsAsFactors = FALSE
  )
  write_tsv_safe(run_manifest, file.path(output_dir, "manifest_run.tsv"))
}

write_failed_run_manifests <- function(
    output_dir,
    project_dir,
    mode,
    tax_mode,
    error_message) {
  partial_artifacts <- list_run_artifacts(output_dir)
  artifact_manifest <- write_tsv_safe(
    partial_artifacts,
    file.path(output_dir, "manifest_failed_artifacts.tsv")
  )
  run_manifest <- write_run_manifest(
    output_dir = output_dir,
    project_dir = project_dir,
    mode = mode,
    tax_mode = tax_mode,
    status = "failed",
    error_message = as.character(error_message),
    artifact_manifest = artifact_manifest
  )
  list(artifacts = artifact_manifest, run = run_manifest)
}

write_success_run_manifest <- function(
    output_dir,
    project_dir,
    mode,
    tax_mode,
    combined_manifest) {
  write_run_manifest(
    output_dir = output_dir,
    project_dir = project_dir,
    mode = mode,
    tax_mode = tax_mode,
    status = "success",
    combined_manifest = combined_manifest
  )
}

# ---- Command-line entry point ---------------------------------------------

# Validate inputs before loading SQM, then run only the requested analysis modes.
main_impl <- function() {
  progress_message("Starting script")

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L || any(args %in% c("--help", "-h"))) {
    print_help()
    return(invisible(0))
  }

  named_args <- parse_named_args(args)
  allowed_args <- c(
    "help", "project_dir", "output_dir", "mode", "pathways", "pathway_selection_modes", "pathway_top_n", "samples", "tax_mode",
    "top_n_ko", "top_n_taxa", "taxa", "taxonomy_ranks", "taxonomy_counts",
    "flowplot_formats", "pathview_sample_modes", "enzyme_ecs", "enzyme_plot_types",
    "dimensions", "plot_dpi"
  )
  unknown_args <- setdiff(names(named_args), allowed_args)
  if (length(unknown_args) > 0L) {
    stop("Unknown options: ", paste(unknown_args, collapse = ", "), call. = FALSE)
  }

  if (is.null(named_args$project_dir) || !nzchar(named_args$project_dir)) {
    stop("Missing required argument: --project_dir", call. = FALSE)
  }
  if (is.null(named_args$output_dir) || !nzchar(named_args$output_dir)) {
    stop("Missing required argument: --output_dir", call. = FALSE)
  }
  if (is.null(named_args$mode) || !nzchar(named_args$mode)) {
    stop("Missing required argument: --mode", call. = FALSE)
  }

  project_dir <- named_args$project_dir
  output_dir <- named_args$output_dir
  mode <- named_args$mode
  allowed_modes <- c("all", "flow", "funz", "enzimi", "taxon", "pathview", "pie")
  if (!mode %in% allowed_modes) {
    stop("Invalid mode: ", mode, ". Allowed values: all, flow, funz, enzimi, taxon, pathview, pie.", call. = FALSE)
  }
  if (!dir.exists(project_dir)) {
    stop("Project directory does not exist: ", project_dir, call. = FALSE)
  }

  runtime_flow_formats <- if (is.null(named_args$flowplot_formats)) {
    c("png", "html")
  } else {
    trimws(strsplit(named_args$flowplot_formats, ",", fixed = TRUE)[[1L]])
  }
  check_required_packages(
    required_packages_for_mode(mode, runtime_flow_formats),
    mode
  )
  suppressPackageStartupMessages(
    invisible(lapply(common_required_packages, library, character.only = TRUE))
  )

  tax_mode <- validate_tax_mode(
    if (is.null(named_args$tax_mode)) "prokfilter" else named_args$tax_mode
  )
  top_n_ko <- if (is.null(named_args$top_n_ko)) {
    20L
  } else {
    parse_positive_integer_arg(named_args$top_n_ko, "top_n_ko")
  }
  top_n_taxa <- if (is.null(named_args$top_n_taxa)) {
    15L
  } else {
    parse_positive_integer_arg(named_args$top_n_taxa, "top_n_taxa")
  }
  pathway_top_n <- if (is.null(named_args$pathway_top_n)) {
    default_pathway_top_n
  } else {
    parse_positive_integer_arg(named_args$pathway_top_n, "pathway_top_n")
  }
  plot_dpi <- if (is.null(named_args$plot_dpi)) {
    600
  } else {
    suppressWarnings(as.numeric(named_args$plot_dpi))
  }
  validate_positive_integer(top_n_ko, "top_n_ko")
  validate_positive_integer(top_n_taxa, "top_n_taxa")
  validate_positive_integer(pathway_top_n, "pathway_top_n")
  if (
    length(plot_dpi) != 1L || is.na(plot_dpi) ||
      !is.finite(plot_dpi) || plot_dpi <= 0
  ) {
    stop("plot_dpi must be a positive number.", call. = FALSE)
  }

  taxonomy_ranks <- if (is.null(named_args$taxonomy_ranks)) {
    default_taxonomy_ranks
  } else {
    split_csv_arg(named_args$taxonomy_ranks)
  }
  if (length(taxonomy_ranks) == 0L) {
    stop("taxonomy_ranks cannot be empty.", call. = FALSE)
  }

  taxonomy_counts <- normalize_taxonomy_counts(named_args$taxonomy_counts)
  flowplot_formats <- normalize_flowplot_formats(named_args$flowplot_formats)
  check_flow_html_preflight(mode, flowplot_formats)
  pathview_sample_modes <- normalize_pathview_sample_modes(named_args$pathview_sample_modes)
  pathway_selection_modes <- normalize_pathway_selection_modes(named_args$pathway_selection_modes)
  pie_selection_modes <- pie_pathway_selection_modes(
    pathway_selection_modes,
    explicitly_requested = !is.null(named_args$pathway_selection_modes)
  )
  enzyme_ecs <- normalize_enzyme_ecs(named_args$enzyme_ecs)
  enzyme_plot_types <- normalize_enzyme_plot_types(named_args$enzyme_plot_types)

  requested_pathways <- if (is.null(named_args$pathways)) character() else split_csv_arg(named_args$pathways)
  requested_taxa <- if (is.null(named_args$taxa)) character() else split_csv_arg(named_args$taxa)
  
  dimensions <- parse_dimensions(named_args)

  progress_message(
    "Configurazione | mode=", mode,
    " | tax_mode=", tax_mode,
    " | top_n_ko=", top_n_ko,
    " | top_n_taxa=", top_n_taxa,
    " | pathway_selection_modes=", paste(pathway_selection_modes, collapse = ","),
    " | pathway_top_n=", pathway_top_n,
    " | pathview_sample_modes=", paste(pathview_sample_modes, collapse = ","),
    " | enzyme_ecs=", paste(enzyme_ecs, collapse = ","),
    " | enzyme_plot_types=", paste(enzyme_plot_types, collapse = ","),
    " | dpi=", plot_dpi
  )

  progress_message("Loading SQM project: ", project_dir)
  sqm <- load_sqm_project(project_dir, tax_mode)
  validate_sqm_object(sqm)
  validate_taxonomy_ranks(taxonomy_ranks, colnames(sqm$orfs$tax))

  available_samples <- colnames(sqm$orfs$tpm)
  sample_selection <- resolve_sample_selection(named_args$samples, available_samples)
  selected_samples <- sample_selection$samples
  sample_order_basis <- sample_selection$basis
  update_run_sample_order(selected_samples, sample_order_basis)
  invisible(vapply(selected_samples, safe_output_component, character(1)))
  validate_sqm_tpm_inputs(
    sqm,
    selected_samples,
    require_kegg_tpm = mode %in% c(
      "all", "flow", "funz", "pie", "taxon", "pathview", "enzimi"
    )
  )
  progress_message("Selected samples: ", paste(selected_samples, collapse = ", "))

  resolved_defined_pathways <- list()
  if (mode != "enzimi" && "defined" %in% pathway_selection_modes) {
    defined_pathways <- if (length(requested_pathways) == 0L) {
      default_pathway_ids
    } else {
      requested_pathways
    }
    resolved_defined_pathways <- resolve_pathways(sqm, defined_pathways)
  }

  filter_contexts <- if (length(requested_taxa) == 0L) {
    list(list(
      sqm = sqm,
      output_dir = output_dir,
      filtered_taxon = NA_character_,
      filtered_taxon_rank = NA_character_
    ))
  } else {
    progress_message("Resolving taxonomic filters: ", paste(requested_taxa, collapse = ", "))
    resolved_taxa <- resolve_taxa_filters(sqm, requested_taxa)
    map(resolved_taxa, function(taxon_info) {
      progress_message("Applying taxon filter: ", taxon_info$taxon, " @ ", taxon_info$rank)
      filtered_sqm <- subset_sqm_by_taxon(sqm, taxon_info$orf_ids)
      validate_sqm_object(filtered_sqm)
      validate_sqm_tpm_inputs(
        filtered_sqm,
        selected_samples,
        require_kegg_tpm = mode %in% c(
          "all", "flow", "funz", "pie", "taxon", "pathview", "enzimi"
        )
      )
      list(
        sqm = filtered_sqm,
        output_dir = file.path(
          output_dir,
          "taxon_filter",
          taxon_info$rank,
          safe_output_component(taxon_info$taxon, max_length = 24L)
        ),
        filtered_taxon = taxon_info$taxon,
        filtered_taxon_rank = taxon_info$rank
      )
    })
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  script_name <- "sqm_plots.R"

  for (context in filter_contexts) {
    context_output_dir <- context$output_dir
    context_sqm <- context$sqm
    manifests <- new_context_manifest_registry(context_output_dir)
    activate_context_manifest_registry(manifests)
    progress_message(
      "Output context: ", context_output_dir,
      if (!is.na(context$filtered_taxon)) paste0(" | taxon filter=", context$filtered_taxon, " @ ", context$filtered_taxon_rank) else ""
    )

    pathway_groups <- if (mode == "enzimi") {
      list()
    } else {
      select_context_pathway_groups(
        context_sqm = context_sqm,
        resolved_defined_pathways = resolved_defined_pathways,
        pathway_selection_modes = pathway_selection_modes,
        selected_samples = selected_samples,
        pathway_top_n = pathway_top_n
      )
    }
    if (length(pathway_groups) > 0L) {
      progress_message(
        "Resolved pathways for context ", context_output_dir, ": ",
        paste(
          unlist(map(
            pathway_groups,
            ~ vapply(.x, `[[`, character(1), "canonical_pathway_name")
          )),
          collapse = " | "
        )
      )
    }
    pathway_entries <- imap(pathway_groups, function(pathways, pathway_selection) {
      map(pathways, function(pathway_info) {
        list(
          pathway_name = pathway_info$canonical_pathway_name,
          pathway_id = pathway_info$pathway_id,
          pathway_selection = pathway_selection
        )
      })
    }) |>
      purrr::flatten()

    context_label <- if (is.na(context$filtered_taxon)) {
      "global"
    } else {
      paste0(context$filtered_taxon, "@", context$filtered_taxon_rank)
    }
    prepared_pathways <- prepare_context_pathway_subsets(
      context_sqm = context_sqm,
      pathway_entries = pathway_entries,
      context_label = context_label,
      include_kegg_oracle = mode %in% c("all", "flow", "funz", "pie", "taxon"),
      selected_samples = selected_samples
    )
    pathway_sqms <- prepared_pathways$pathway_sqms
    record_pathway_skips(prepared_pathways$skips)

    if (mode %in% c("all", "funz")) {
      progress_message("Starting FUNZ section")
    }
    if (mode %in% c("all", "flow")) {
      progress_message("Starting FLOW section")
    }
    if (mode %in% c("all", "taxon")) {
      progress_message("Starting TAXON section")
      manifests <- run_taxonomy_scope(
        sqm_object = context_sqm,
        output_dir = context_output_dir,
        manifest_base_dir = context_output_dir,
        output_manifests = manifests,
        script_name = script_name,
        project_dir = project_dir,
        tax_mode = tax_mode,
        pathway_name = NA_character_,
        selected_samples = selected_samples,
        dimensions = dimensions,
        plot_dpi = plot_dpi,
        top_n_taxa = top_n_taxa,
        top_n_ko = top_n_ko,
        taxonomy_ranks = taxonomy_ranks,
        taxonomy_counts = taxonomy_counts,
        scope_name = "taxonomy_global",
        ignore_unmapped = TRUE,
        ignore_unclassified = TRUE,
        pathway_id = NA_character_,
        filtered_taxon = context$filtered_taxon,
        filtered_taxon_rank = context$filtered_taxon_rank
      )
    }
    if (mode %in% c("all", "pathview")) {
      progress_message("Starting PATHVIEW section")
    }
    if (mode %in% c("all", "pie")) {
      progress_message("Starting PIE section")
    }

    for (pathway_key in names(pathway_sqms)) {
      pathway_info <- pathway_sqms[[pathway_key]]
      pathway_name <- pathway_info$pathway_name
      pie_enabled <- mode %in% c("all", "pie") &&
        pathway_info$pathway_selection %in% pie_selection_modes
      needs_pathway_analysis <- mode %in% c("all", "funz", "flow", "taxon") ||
        pie_enabled
      pathway_analysis <- if (needs_pathway_analysis) {
        build_pathway_analysis(pathway_info, selected_samples)
      } else {
        NULL
      }

      if (mode %in% c("all", "funz")) {
        manifests <- run_funz_mode(
          output_dir = context_output_dir,
          manifest_base_dir = file.path(context_output_dir, "funz"),
          output_manifests = manifests,
          script_name = script_name,
          project_dir = project_dir,
          tax_mode = tax_mode,
          pathway_name = pathway_name,
          pathway_sqm = pathway_info$pathway_sqm,
          selected_samples = selected_samples,
          dimensions = dimensions,
          plot_dpi = plot_dpi,
          top_n_taxa = top_n_taxa,
          top_n_ko = top_n_ko,
          pathway_id = pathway_info$pathway_id,
          pathway_selection = pathway_info$pathway_selection,
          filtered_taxon = context$filtered_taxon,
          filtered_taxon_rank = context$filtered_taxon_rank,
          pathway_analysis = pathway_analysis
        )
      }

      if (mode %in% c("all", "flow")) {
        manifests <- run_flow_mode(
          output_dir = context_output_dir,
          manifest_base_dir = file.path(context_output_dir, "flowplot"),
          output_manifests = manifests,
          script_name = script_name,
          project_dir = project_dir,
          tax_mode = tax_mode,
          pathway_name = pathway_name,
          pathway_sqm = pathway_info$pathway_sqm,
          selected_samples = selected_samples,
          dimensions = dimensions,
          plot_dpi = plot_dpi,
          top_n_taxa = top_n_taxa,
          top_n_ko = top_n_ko,
          taxonomy_ranks = taxonomy_ranks,
          flowplot_formats = flowplot_formats,
          pathway_id = pathway_info$pathway_id,
          pathway_selection = pathway_info$pathway_selection,
          filtered_taxon = context$filtered_taxon,
          filtered_taxon_rank = context$filtered_taxon_rank,
          pathway_analysis = pathway_analysis
        )
      }

      if (mode %in% c("all", "taxon")) {
          manifests <- run_taxonomy_scope(
            sqm_object = pathway_analysis$pathway_sqm,
            output_dir = context_output_dir,
            manifest_base_dir = context_output_dir,
          output_manifests = manifests,
          script_name = script_name,
          project_dir = project_dir,
          tax_mode = tax_mode,
          pathway_name = pathway_name,
          selected_samples = selected_samples,
          dimensions = dimensions,
          plot_dpi = plot_dpi,
          top_n_taxa = top_n_taxa,
          top_n_ko = top_n_ko,
          taxonomy_ranks = taxonomy_ranks,
          taxonomy_counts = taxonomy_counts,
          scope_name = "taxonomy_by_pathway",
          ignore_unmapped = FALSE,
          ignore_unclassified = FALSE,
          pathway_id = pathway_info$pathway_id,
          pathway_selection = pathway_info$pathway_selection,
          filtered_taxon = context$filtered_taxon,
          filtered_taxon_rank = context$filtered_taxon_rank
        )
      }

      if (mode %in% c("all", "pathview")) {
        if (!pathview_is_exportable(
          pathway_info$pathway_selection,
          pathway_info$pathway_id
        )) {
          warning(
            "Skipping PATHVIEW: pathway has no resolvable KEGG ID: ",
            pathway_name,
            call. = FALSE
          )
        } else {
          manifests <- run_pathview_mode(
            sqm_object = context_sqm,
            output_dir = context_output_dir,
            manifest_base_dir = context_output_dir,
            output_manifests = manifests,
            script_name = script_name,
            project_dir = project_dir,
            tax_mode = tax_mode,
            pathway_name = pathway_name,
            pathway_id = pathway_info$pathway_id,
            pathway_selection = pathway_info$pathway_selection,
            selected_samples = selected_samples,
            top_n_taxa = top_n_taxa,
            top_n_ko = top_n_ko,
            pathview_sample_modes = pathview_sample_modes,
            filtered_taxon = context$filtered_taxon,
            filtered_taxon_rank = context$filtered_taxon_rank
          )
        }
      }

      if (pie_enabled) {
        manifests <- run_pie_mode(
          output_dir = context_output_dir,
          manifest_base_dir = file.path(context_output_dir, "pie"),
          output_manifests = manifests,
          script_name = script_name,
          project_dir = project_dir,
          tax_mode = tax_mode,
          pathway_name = pathway_name,
          pathway_sqm = pathway_info$pathway_sqm,
          selected_samples = selected_samples,
          taxonomy_ranks = taxonomy_ranks,
          dimensions = dimensions,
          plot_dpi = plot_dpi,
          top_n_taxa = top_n_taxa,
          top_n_ko = top_n_ko,
          pathway_id = pathway_info$pathway_id,
          pathway_selection = pathway_info$pathway_selection,
          filtered_taxon = context$filtered_taxon,
          filtered_taxon_rank = context$filtered_taxon_rank,
          pathway_analysis = pathway_analysis
        )
      }
    }

    if (mode %in% c("all", "funz", "enzimi")) {
      progress_message("Starting ENZIMI section")
      manifests <- run_enzyme_mode(
        sqm_object = context_sqm,
        output_dir = context_output_dir,
        manifest_base_dir = file.path(context_output_dir, "funz"),
        output_manifests = manifests,
        script_name = script_name,
        project_dir = project_dir,
        tax_mode = tax_mode,
        selected_samples = selected_samples,
        dimensions = dimensions,
        plot_dpi = plot_dpi,
        top_n_taxa = top_n_taxa,
        top_n_ko = top_n_ko,
        enzyme_ecs = enzyme_ecs,
        enzyme_plot_types = enzyme_plot_types,
        sample_order_basis = sample_order_basis,
        filtered_taxon = context$filtered_taxon,
        filtered_taxon_rank = context$filtered_taxon_rank
      )
    }
    flush_context_manifests(manifests)
    clear_active_manifest_registry()
  }

  message("Output directory: ", output_dir)
  invisible(list(status = 0L))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L || any(args %in% c("--help", "-h"))) {
    return(main_impl())
  }

  output_dir <- bootstrap_cli_value(args, "output_dir")
  if (is.null(output_dir) || !nzchar(output_dir)) {
    return(main_impl())
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("Unable to create output directory: ", output_dir, call. = FALSE)
  }
  run_id <- allocate_run_id(output_dir)
  preliminary_sample_arg <- bootstrap_cli_value(args, "samples")
  sample_basis <- if (is.null(preliminary_sample_arg)) "sqm_column_order" else "cli"
  preliminary_samples <- if (is.null(preliminary_sample_arg)) {
    character()
  } else {
    values <- trimws(strsplit(preliminary_sample_arg, ",", fixed = TRUE)[[1L]])
    values[nzchar(values)]
  }
  set_run_context(
    run_id = run_id,
    started_at = Sys.time(),
    sample_order_basis = sample_basis,
    samples = preliminary_samples,
    cli_args = args
  )
  on.exit(clear_run_context(), add = TRUE)

  project_dir <- bootstrap_cli_value(args, "project_dir")
  if (is.null(project_dir)) project_dir <- NA_character_
  mode <- bootstrap_cli_value(args, "mode")
  if (is.null(mode)) mode <- NA_character_
  tax_mode <- bootstrap_cli_value(args, "tax_mode")
  if (is.null(tax_mode)) tax_mode <- "prokfilter"

  log_path <- initialize_run_log(
    output_dir = output_dir,
    project_dir = project_dir,
    mode = mode,
    tax_mode = tax_mode
  )
  record_legacy_outputs(output_dir)

  run_error <- NULL
  result <- tryCatch(
    withCallingHandlers(
      main_impl(),
      warning = function(condition) {
        record_run_warning(conditionMessage(condition))
      }
    ),
    error = function(error) {
      run_error <<- error
      NULL
    }
  )

  if (!is.null(run_error)) {
    partial_manifest_error <- tryCatch(
      {
        flush_active_context_manifests()
        NULL
      },
      error = function(error) error
    )
    if (!is.null(partial_manifest_error)) {
      append_run_log(
        "ERROR",
        paste0("Unable to flush partial section manifests: ", conditionMessage(partial_manifest_error))
      )
    }
    finalize_run_log("FAILED", conditionMessage(run_error))
    stop(run_error)
  }

  finalize_run_log("SUCCESS")
  message("Run ID: ", run_id)
  message("Run log: ", log_path)
  invisible(result$status)
}

main()
