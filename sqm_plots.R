# SqueezeMeta plotting command-line program.
#
# The script loads a validated SQM project once, derives ORF-level data from
# that object, and writes reproducible plots plus manifest files for each mode.
# CLI mode names and output directory names are kept stable for compatibility.
library(SQMtools)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggalluvial)
library(glue)
library(purrr)
library(plotly)
library(htmlwidgets)
library(scales)

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
  "4.4.1.34", "5.2.1.2"
)
default_enzyme_plot_types <- c("bar", "line")
all_taxonomy_columns <- c(
  "superkingdom", "phylum", "class", "order", "family", "genus", "species"
)

# Curated pathway IDs supported by the explicit `defined` selection mode.
known_pathways <- tibble::tribble(
  ~pathway_id, ~canonical_pathway_name,
  "00361", "Chlorocyclohexane and chlorobenzene degradation",
  "00710", "Carbon fixation in photosynthetic organisms",
  "00623", "Toluene degradation",
  "00621", "Dioxin degradation",
  "00625", "Chloroalkane and chloroalkene degradation",
  "00630", "Glyoxylate and dicarboxylate metabolism",
  "00633", "Nitrotoluene degradation",
  "00910", "Nitrogen metabolism",
  "00980", "Metabolism of xenobiotics by cytochrome P450"
)
default_pathway_ids <- known_pathways$pathway_id
default_pathway_selection_modes <- c("defined", "top20")
default_pathway_top_n <- 20L

# ---- Command-line parsing and shared configuration -----------------------

check_required_packages <- function() {
  required <- c(
    "SQMtools", "readr", "dplyr", "tidyr", "tibble", "stringr", "ggplot2",
    "ggalluvial", "glue", "purrr", "plotly", "htmlwidgets", "scales"
  )
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

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
    "  --pathway_top_n=NUM        Default: 20.\n",
    "  --samples=LIST             Comma-separated samples.\n",
    "  --tax_mode=MODE            Default: prokfilter.\n",
    "  --top_n_ko=NUM             Default: 20.\n",
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
    "  Pie mode creates one plot for each sample x KO x rank combination.\n",
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

format_dimension_label <- function(x) {
  str_replace(format(x, trim = TRUE, scientific = FALSE), "\\.0+$", "")
}

# Parse named width-by-height presets once so all exports share the same sizes.
parse_dimensions <- function(named_args) {
  if (is.null(named_args$dimensions)) {
    return(default_dimensions)
  }

  labels <- split_csv_arg(named_args$dimensions)
  parsed <- map(labels, function(label) {
    parts <- str_split(label, "x", simplify = TRUE)
    if (ncol(parts) != 2L) {
      stop("Invalid dimension: ", label, call. = FALSE)
    }
    width <- as.numeric(parts[[1]])
    height <- as.numeric(parts[[2]])
    if (is.na(width) || is.na(height) || width <= 0 || height <= 0) {
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

write_tsv_safe <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(data, path, na = "NA")
}

progress_message <- function(..., .prefix = "[sqm_plots]") {
  message(.prefix, " ", paste0(..., collapse = ""))
}

save_png_dimensions <- function(plot_object, output_dir, file_stem, dimensions, dpi) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_files <- character()
  for (dim_name in names(dimensions)) {
    dims <- dimensions[[dim_name]]
    output_file <- file.path(output_dir, paste0(file_stem, "_", dim_name, ".png"))
    ggplot2::ggsave(
      filename = output_file,
      plot = plot_object,
      width = dims[["width"]],
      height = dims[["height"]],
      units = "in",
      dpi = dpi,
      bg = "white"
    )
    output_files[[dim_name]] <- output_file
  }
  output_files
}

save_html_widget <- function(widget, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  htmlwidgets::saveWidget(widget, file = path, selfcontained = FALSE)
  path
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
    pathway_id = NA_character_) {
  tibble::tibble(
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
    filtered_taxon_rank = filtered_taxon_rank
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

normalize_taxon_value <- function(x) {
  x <- as.character(x)
  empty <- is.na(x) | !nzchar(trimws(x))
  x[empty] <- "Unclassified"
  x
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
  if (is.na(x) || x <= 0 || x != as.integer(x)) {
    stop(arg_name, " must be a positive integer.", call. = FALSE)
  }
}

# Resolve user-supplied pathway codes or names to canonical SQM pathway names.
resolve_pathways <- function(sqm, requested_pathways) {
  raw_candidates <- unique(as.character(sqm$misc$KEGG_paths))
  raw_candidates <- raw_candidates[!is.na(raw_candidates) & nzchar(raw_candidates)]
  canonical_candidates <- raw_candidates |>
    str_split("\\s*\\|\\s*") |>
    unlist(use.names = FALSE) |>
    str_split("\\s*;\\s*") |>
    map_chr(~ .x[[length(.x)]]) |>
    unique()

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

    matches <- canonical_candidates[tolower(canonical_candidates) == tolower(requested)]
    matches <- unique(matches)
    if (length(matches) == 0L) {
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
    if (length(matches) > 1L) {
      stop(
        "Ambiguous pathway: ", requested, ". Matches: ",
        paste(matches, collapse = "; "),
        call. = FALSE
      )
    }

    canonical_name <- matches[[1]]
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

split_kegg_pathway_field <- function(pathway_field) {
  pathway_values <- unlist(str_split(as.character(pathway_field), "\\s*\\|\\s*"), use.names = FALSE)
  pathway_values <- pathway_values[!is.na(pathway_values) & nzchar(trimws(pathway_values))]
  canonical_values <- map_chr(
    str_split(pathway_values, "\\s*;\\s*"),
    ~ trimws(.x[[length(.x)]])
  )
  unique(canonical_values[nzchar(canonical_values)])
}

pathway_id_for_name <- function(pathway_name) {
  match <- known_pathways |>
    filter(tolower(.data$canonical_pathway_name) == tolower(pathway_name))
  if (nrow(match) == 0L) {
    NA_character_
  } else {
    match$pathway_id[[1]]
  }
}

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

  pathway_membership <- orf_table |>
    transmute(orf_id = .data$orf_id, KEGGPATH = as.character(.data$KEGGPATH)) |>
    filter(!is.na(.data$KEGGPATH)) |>
    mutate(canonical_pathway_name = map(.data$KEGGPATH, split_kegg_pathway_field)) |>
    tidyr::unnest_longer("canonical_pathway_name") |>
    distinct(.data$orf_id, .data$canonical_pathway_name)

  pathway_totals <- tpm_table |>
    pivot_longer(
      cols = -all_of("orf_id"),
      names_to = "sample",
      values_to = "tpm"
    ) |>
    filter(.data$sample %in% selected_samples) |>
    inner_join(pathway_membership, by = "orf_id", relationship = "many-to-many") |>
    group_by(.data$canonical_pathway_name) |>
    summarise(total_tpm = sum(as.numeric(.data$tpm), na.rm = TRUE), .groups = "drop") |>
    arrange(desc(.data$total_tpm), .data$canonical_pathway_name) |>
    slice_head(n = pathway_top_n)

  map(seq_len(nrow(pathway_totals)), function(index) {
    pathway_name <- pathway_totals$canonical_pathway_name[[index]]
    list(
      input_value = pathway_name,
      pathway_id = pathway_id_for_name(pathway_name),
      canonical_pathway_name = pathway_name,
      total_tpm = pathway_totals$total_tpm[[index]]
    )
  })
}

select_pathway_groups <- function(
    sqm,
    requested_pathways,
    pathway_selection_modes,
    selected_samples,
    pathway_top_n = default_pathway_top_n) {
  pathway_selection_modes <- normalize_pathway_selection_modes(pathway_selection_modes)
  groups <- list()

  if ("defined" %in% pathway_selection_modes) {
    defined_pathways <- if (length(requested_pathways) == 0L) default_pathway_ids else requested_pathways
    groups[["defined"]] <- resolve_pathways(sqm, defined_pathways)
  }
  if ("top20" %in% pathway_selection_modes) {
    groups[["top20"]] <- select_top_pathways(sqm, selected_samples, pathway_top_n)
  }

  groups
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

validate_samples <- function(requested_samples, available_samples) {
  missing_samples <- setdiff(requested_samples, available_samples)
  if (length(missing_samples) > 0L) {
    stop(
      "Samples not found: ", paste(missing_samples, collapse = ", "),
      ". Available: ", paste(available_samples, collapse = ", "),
      call. = FALSE
    )
  }
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

subset_pathway <- function(sqm, canonical_pathway) {
  SQMtools::subsetFun(
    SQM = sqm,
    fun = canonical_pathway,
    columns = "KEGGPATH",
    ignore_case = FALSE,
    fixed = TRUE
  )
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
    ignore_unclassified_functions = TRUE,
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

# Build the canonical ORF × sample × KO table used by functional and flow plots.
build_orf_long_table <- function(pathway_sqm, selected_samples) {
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

  fun_lookup <- sqm_misc_names <- pathway_sqm$misc$KEGG_names

  annotations <- orf_table |>
    transmute(
      orf_id = .data$orf_id,
      `KEGG ID` = as.character(.data[["KEGG ID"]]),
      kegg_function_raw = as.character(.data[["KEGGFUN"]]),
      ec_codes = extract_ec_codes(.data[["KEGGFUN"]]),
      KEGGPATH = as.character(.data[["KEGGPATH"]])
    )

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

  no_ko_orfs <- joined |>
    distinct(.data$orf_id, .data[["KEGG ID"]]) |>
    mutate(has_ko = lengths(map(.data[["KEGG ID"]], extract_ko_ids)) > 0L) |>
    filter(!.data$has_ko) |>
    nrow()

  expanded <- joined |>
    mutate(ko_ids = map(.data[["KEGG ID"]], extract_ko_ids)) |>
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

  attr(expanded, "excluded_orfs_without_ko") <- no_ko_orfs
  expanded
}

# ---- Functional, flow, and taxonomy plot data ----------------------------

# Use SQM KEGG names only as a fallback when the ORF annotation is absent.
get_ko_name_lookup <- function(pathway_sqm) {
  kegg_names <- pathway_sqm$misc$KEGG_names
  if (is.null(kegg_names)) {
    return(setNames(character(), character()))
  }
  stats::setNames(as.character(kegg_names), names(kegg_names))
}

extract_ko_ec_lookup <- function(orf_long) {
  orf_long |>
    distinct(.data$ko_id, .data$ec_codes) |>
    group_by(.data$ko_id) |>
    summarise(ec_codes = first(stats::na.omit(.data$ec_codes)), .groups = "drop")
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

build_ko_plot_table <- function(orf_long, selected_samples, top_n_ko, ko_lookup) {
  summary_tbl <- orf_long |>
    filter(.data$sample %in% selected_samples) |>
    group_by(.data$sample, .data$ko_id) |>
    summarise(
      tpm = sum(.data$tpm),
      kegg_function = first(na.omit(.data$kegg_function)),
      ec_codes = first(na.omit(.data$ec_codes)),
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
        coalesce(.data$kegg_function, unname(as.character(ko_lookup[.data$ko_id])), .data$ko_id),
        "Other KO outside top N"
      ),
      plot_ec = if_else(.data$ko_id %in% top_ko_ids, coalesce(.data$ec_codes, NA_character_), NA_character_)
    ) |>
    group_by(.data$sample, .data$plot_ko_id) |>
    summarise(
      tpm = sum(.data$tpm),
      kegg_function = first(na.omit(.data$plot_function)),
      ec_codes = first(na.omit(.data$plot_ec)),
      .groups = "drop"
    ) |>
    rename(ko_id = plot_ko_id)

  sample_totals <- collapsed |>
    group_by(.data$sample) |>
    summarise(sample_pathway_total_tpm = sum(.data$tpm), .groups = "drop")

  plot_ko_ids <- c(
    top_ko_ids,
    if (any(collapsed$ko_id == "Other")) "Other" else character()
  )

  ko_info <- collapsed |>
    filter(.data$ko_id %in% plot_ko_ids) |>
    group_by(.data$ko_id) |>
    summarise(
      kegg_function = first(na.omit(.data$kegg_function)),
      ec_codes = first(na.omit(.data$ec_codes)),
      .groups = "drop"
    )

  plot_tbl <- tidyr::expand_grid(sample = selected_samples, ko_id = plot_ko_ids) |>
    left_join(
      collapsed |>
        filter(.data$ko_id %in% plot_ko_ids),
      by = c("sample", "ko_id")
    ) |>
    left_join(ko_info, by = "ko_id", suffix = c("", "_info")) |>
    mutate(
      tpm = replace_na(.data$tpm, 0),
      kegg_function = coalesce(.data$kegg_function, .data$kegg_function_info, .data$ko_id),
      ec_codes = coalesce(.data$ec_codes, .data$ec_codes_info, "NA")
    ) |>
    select(.data$sample, .data$ko_id, .data$tpm, .data$kegg_function, .data$ec_codes) |>
    left_join(sample_totals, by = "sample") |>
    mutate(
      sample_pathway_percent = if_else(
        .data$sample_pathway_total_tpm > 0,
        100 * .data$tpm / .data$sample_pathway_total_tpm,
        0
      )
    )

  validate_percent_sum(plot_tbl, "sample", "sample_pathway_percent")

  ko_levels <- plot_tbl |>
    mutate(is_other = .data$ko_id == "Other") |>
    group_by(.data$ko_id, .data$is_other) |>
    summarise(total_tpm = sum(.data$tpm), .groups = "drop") |>
    arrange(.data$is_other, desc(.data$total_tpm), .data$ko_id) |>
    pull(.data$ko_id)

  plot_tbl |>
    mutate(
      sample = factor(.data$sample, levels = selected_samples),
      ko_id = factor(.data$ko_id, levels = ko_levels)
    )
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
    select(.data$ko_id, .data$legend_label) |>
    tibble::deframe()
}

# Build visualizations only from the accompanying TSV tables written to disk.
make_ko_barplot <- function(plot_tbl, pathway_name, selected_samples) {
  legend_labels <- build_ko_legend_labels(plot_tbl, selected_samples)
  ko_levels <- levels(plot_tbl$ko_id)
  non_other <- setdiff(ko_levels, "Other")
  palette <- c(
    setNames(rep(colors_hex, length.out = length(non_other)), non_other),
    if ("Other" %in% ko_levels) c(Other = "grey70") else NULL
  )

  ggplot(plot_tbl, aes(x = .data$sample, y = .data$tpm, fill = .data$ko_id)) +
    geom_col(color = "grey25", linewidth = 0.15, width = 0.78) +
    scale_fill_manual(values = palette, labels = legend_labels, drop = FALSE) +
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
  orf_table <- as.data.frame(sqm_object$orfs$table, check.names = FALSE) |>
    tibble::rownames_to_column("orf_id")
  tpm_table <- as.data.frame(sqm_object$orfs$tpm, check.names = FALSE) |>
    tibble::rownames_to_column("orf_id")

  if (anyDuplicated(orf_table$orf_id) > 0L || anyDuplicated(tpm_table$orf_id) > 0L) {
    stop("orf_id keys must be unique in the ORF and TPM tables.", call. = FALSE)
  }
  if (!setequal(orf_table$orf_id, tpm_table$orf_id)) {
    stop("orf_id values do not match between the ORF and TPM tables.", call. = FALSE)
  }
  if (!"KEGGFUN" %in% colnames(orf_table)) {
    stop("The ORF table does not contain the KEGGFUN column.", call. = FALSE)
  }

  missing_samples <- setdiff(selected_samples, colnames(tpm_table))
  if (length(missing_samples) > 0L) {
    stop(
      "Samples missing from sqm$orfs$tpm: ", paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  matched_ecs <- orf_table |>
    transmute(
      orf_id = .data$orf_id,
      ec_codes = extract_ec_codes(.data$KEGGFUN)
    ) |>
    filter(!is.na(.data$ec_codes)) |>
    mutate(ec_code = map(.data$ec_codes, split_ec_code_field)) |>
    tidyr::unnest_longer("ec_code") |>
    filter(.data$ec_code %in% enzyme_ecs) |>
    distinct(.data$orf_id, .data$ec_code)

  tpm_by_enzyme <- tpm_table |>
    pivot_longer(
      cols = -all_of("orf_id"),
      names_to = "sample",
      values_to = "tpm"
    ) |>
    filter(.data$sample %in% selected_samples) |>
    inner_join(matched_ecs, by = "orf_id", relationship = "many-to-many") |>
    group_by(.data$sample, .data$ec_code) |>
    summarise(tpm = sum(as.numeric(.data$tpm), na.rm = TRUE), .groups = "drop")

  tidyr::expand_grid(sample = selected_samples, ec_code = enzyme_ecs) |>
    left_join(tpm_by_enzyme, by = c("sample", "ec_code")) |>
    mutate(
      tpm = replace_na(.data$tpm, 0),
      sample = factor(.data$sample, levels = selected_samples),
      ec_code = factor(.data$ec_code, levels = enzyme_ecs)
    )
}

enzyme_palette <- function(enzyme_ecs) {
  stats::setNames(rep(colors_hex, length.out = length(enzyme_ecs)), enzyme_ecs)
}

make_enzyme_barplot <- function(enzyme_tbl, title) {
  enzyme_ecs <- levels(enzyme_tbl$ec_code)
  ggplot(enzyme_tbl, aes(x = .data$sample, y = .data$tpm, fill = .data$ec_code)) +
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
  enzyme_ecs <- levels(enzyme_tbl$ec_code)
  ggplot(enzyme_tbl, aes(x = .data$sample, y = .data$tpm, color = .data$ec_code, group = .data$ec_code)) +
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

build_flow_table_for_rank <- function(orf_long, rank, selected_samples, top_n_taxa, top_n_ko, ko_lookup) {
  top_taxa <- orf_long |>
    filter(.data$sample %in% selected_samples) |>
    group_by(.data[[rank]]) |>
    summarise(total_tpm = sum(.data$tpm), .groups = "drop") |>
    arrange(desc(.data$total_tpm), .data[[rank]]) |>
    slice_head(n = top_n_taxa) |>
    pull(.data[[rank]])

  top_kos <- orf_long |>
    filter(.data$sample %in% selected_samples) |>
    group_by(.data$ko_id) |>
    summarise(total_tpm = sum(.data$tpm), .groups = "drop") |>
    arrange(desc(.data$total_tpm), .data$ko_id) |>
    slice_head(n = top_n_ko) |>
    pull(.data$ko_id)

  summary_tbl <- orf_long |>
    filter(.data$sample %in% selected_samples) |>
    mutate(
      taxon = if_else(.data[[rank]] %in% top_taxa, .data[[rank]], "Other"),
      KO = if_else(.data$ko_id %in% top_kos, .data$ko_id, "Other")
    ) |>
    group_by(.data$sample, .data$taxon, .data$KO) |>
    summarise(TPM = sum(.data$tpm), .groups = "drop") |>
    filter(.data$TPM > 0)

  ko_meta <- orf_long |>
    distinct(.data$ko_id, .data$kegg_function)

  summary_tbl |>
    left_join(ko_meta, by = c("KO" = "ko_id")) |>
    mutate(
      KO_name = if_else(
        .data$KO == "Other",
        "Other KOs",
        coalesce(.data$kegg_function, unname(as.character(ko_lookup[.data$KO])), .data$KO)
      )
    ) |>
    select(-"kegg_function")
}

build_flow_table_for_sample <- function(flow_rank_table, pathway_name, rank, sample_name) {
  sample_tbl <- flow_rank_table |>
    filter(.data$sample == sample_name) |>
    group_by(.data$taxon, .data$KO, .data$KO_name) |>
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
      "KO_name", "TPM", "taxon_percent", "KO_percent",
      "flow_percent"
    )

  validate_percent_sum(flow_tbl, "sample", "flow_percent")
  flow_tbl
}

make_flow_plot <- function(flow_tbl, pathway_name, rank, sample_name) {
  taxon_levels <- levels(flow_tbl$taxon)
  ko_levels <- levels(flow_tbl$KO)
  lodes_tbl <- ggalluvial::to_lodes_form(flow_tbl, axes = c("taxon", "KO"), discern = FALSE) |>
    mutate(x = factor(.data$x, levels = c("taxon", "KO"), labels = c("Taxon", "KO")))
  tax_palette <- c(
    stats::setNames(rep(colors_hex, length.out = length(setdiff(taxon_levels, "Other"))), setdiff(taxon_levels, "Other")),
    if ("Other" %in% taxon_levels) c(Other = "grey70") else NULL
  )
  ko_palette <- c(
    stats::setNames(rep(colors_hex, length.out = length(setdiff(ko_levels, "Other"))), setdiff(ko_levels, "Other")),
    if ("Other" %in% ko_levels) c(Other = "grey70") else NULL
  )

  ggplot(
    flow_tbl,
    aes(axis1 = .data$taxon, axis2 = .data$KO, y = .data$flow_percent)
  ) +
    ggalluvial::geom_alluvium(aes(fill = .data$taxon), alpha = 0.78, width = 1 / 12) +
    ggalluvial::geom_stratum(
      data = dplyr::filter(lodes_tbl, .data$x == "Taxon"),
      aes(x = .data$x, stratum = .data$stratum, alluvium = .data$alluvium, y = .data$flow_percent, fill = .data$stratum),
      inherit.aes = FALSE,
      width = 1 / 5,
      color = "grey35",
      linewidth = 0.25
    ) +
    ggalluvial::geom_stratum(
      data = dplyr::filter(lodes_tbl, .data$x == "KO"),
      aes(x = .data$x, stratum = .data$stratum, alluvium = .data$alluvium, y = .data$flow_percent, fill = .data$stratum),
      inherit.aes = FALSE,
      width = 1 / 5,
      color = "grey35",
      linewidth = 0.25
    ) +
    geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(stratum)),
      size = 2.8
    ) +
    scale_x_discrete(limits = c("Taxon", "KO"), expand = c(0.08, 0.08)) +
    scale_fill_manual(values = c(tax_palette, ko_palette[setdiff(names(ko_palette), names(tax_palette))]), drop = FALSE) +
    labs(
      title = paste0("Flowplot - ", pathway_name, " - ", rank, " - ", sample_name),
      subtitle = "Taxonomy-to-KO flow from the ORF x sample x KO table",
      x = NULL,
      y = "Relative flow (%)",
      fill = rank
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

make_flow_sankey <- function(flow_tbl, pathway_name, rank, sample_name) {
  tax_labels <- levels(flow_tbl$taxon)
  ko_labels <- levels(flow_tbl$KO)
  tax_palette <- c(
    stats::setNames(rep(colors_hex, length.out = length(setdiff(tax_labels, "Other"))), setdiff(tax_labels, "Other")),
    if ("Other" %in% tax_labels) c(Other = "grey70") else NULL
  )
  ko_palette <- c(
    stats::setNames(rep(colors_hex, length.out = length(setdiff(ko_labels, "Other"))), setdiff(ko_labels, "Other")),
    if ("Other" %in% ko_labels) c(Other = "grey70") else NULL
  )
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
      label = c(tax_labels, ko_labels),
      color = c(unname(tax_palette[tax_labels]), unname(ko_palette[ko_labels])),
      x = c(rep(0.02, length(tax_labels)), rep(0.98, length(ko_labels))),
      y = c(tax_y, ko_y)
    ),
    link = list(
      source = match(as.character(flow_tbl$taxon), tax_labels) - 1L,
      target = length(tax_labels) + match(as.character(flow_tbl$KO), ko_labels) - 1L,
      value = flow_tbl$flow_percent,
      color = unname(grDevices::adjustcolor(tax_palette[as.character(flow_tbl$taxon)], alpha.f = 0.65)),
      customdata = paste0(
        "Taxon: ", flow_tbl$taxon,
        "<br>KO: ", flow_tbl$KO,
        "<br>Name: ", flow_tbl$KO_name,
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
          "<br><sup>ORF-linked taxon -> KO flow based on TPM</sup>"
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
  sqm_subset <- SQMtools::subsetSamples(
    SQM = sqm_object,
    samples = selected_samples
  )

  plot_object <- SQMtools::plotTaxonomy(
    SQM = sqm_subset,
    rank = rank,
    count = count,
    N = top_n_taxa,
    ignore_unmapped = ignore_unmapped,
    ignore_unclassified = ignore_unclassified,
    no_partial_classifications = FALSE,
    color = colors_hex
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

build_pie_chart_table <- function(orf_long, sample_name, ko_id_filter, rank_name, top_n_taxa) {
  rank_sym <- rlang::sym(rank_name)

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

  top_taxa <- head(base_tbl$taxon_rank, top_n_taxa)
  kept_rows <- base_tbl |>
    filter(.data$taxon_rank %in% top_taxa)
  other_rows <- base_tbl |>
    filter(!.data$taxon_rank %in% top_taxa)

  plot_tbl <- if (nrow(other_rows) > 0L) {
    bind_rows(
      kept_rows,
      tibble::tibble(
        taxon_rank = "Other",
        tpm = sum(other_rows$tpm)
      )
    )
  } else {
    kept_rows
  }

  plot_tbl |>
    mutate(
      total_tpm = sum(.data$tpm),
      pct = if_else(.data$total_tpm > 0, .data$tpm / .data$total_tpm, 0),
      label = if_else(.data$pct >= 0.03, as.character(.data$taxon_rank), ""),
      taxon_rank = forcats::fct_reorder(.data$taxon_rank, .data$tpm, .desc = TRUE)
    )
}

make_pie_plot <- function(
    plot_data,
    pathway_name,
    sample_name,
    ko_id,
    ko_name,
    rank_name,
    ko_ec,
    pathway_sample_tpm) {
  total_tpm <- sum(plot_data$tpm)
  pathway_contribution <- if_else(pathway_sample_tpm > 0, total_tpm / pathway_sample_tpm * 100, 0)

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
      title = glue::glue("Pathway: {pathway_name} - KO {ko_id}"),
      subtitle = glue::glue(
        "Sample: {sample_name} | Rank: {rank_name} | total TPM = {format_display_number(total_tpm)} | pathway contribution: {format_display_number(pathway_contribution, suffix = '%')}"
      ),
      fill = rank_name,
      x = NULL,
      y = NULL,
      caption = glue::glue("KO name: {ko_name}{if (!is.na(ko_ec) && nzchar(ko_ec)) paste0(' | EC: ', ko_ec) else ''}")
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
    filtered_taxon_rank = NA_character_) {
  progress_message(
    "FUNZ | pathway=", pathway_name,
    " | samples=", paste(selected_samples, collapse = ",")
  )
  pathway_dir <- file.path(
    output_dir,
    "funz",
    "pathway",
    pathway_selection_directory(pathway_selection),
    sanitize_name(pathway_name)
  )
  dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

  orf_long <- build_orf_long_table(pathway_sqm, selected_samples)
  if (nrow(orf_long) == 0L) {
    warning("No positive ORF data for pathway ", pathway_name, ".", call. = FALSE)
    return(output_manifests)
  }

  ko_lookup <- get_ko_name_lookup(pathway_sqm)
  plot_tbl <- build_ko_plot_table(orf_long, selected_samples, top_n_ko, ko_lookup)
  plot_path <- file.path(pathway_dir, "barplot_ko_data.tsv")
  progress_message("FUNZ | writing data: ", plot_path)
  write_tsv_safe(plot_tbl, plot_path)

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
      pathway_id = pathway_id
    )
  )

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
        pathway_id = pathway_id
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

  append_manifest_entry <- function(manifest_tbl, output_file, output_type, output_scope, ec_code = NA_character_) {
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
        top_n_ko = top_n_ko,
        output_type = output_type,
        output_file = relative_to_output(output_file, manifest_base_dir),
        mode = "enzimi",
        format = infer_format_from_path(output_file),
        dpi = plot_dpi,
        output_scope = output_scope,
        filtered_taxon = filtered_taxon,
        filtered_taxon_rank = filtered_taxon_rank
      )
    )
  }

  combined_dir <- file.path(enzyme_root, "insieme")
  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  combined_data_file <- file.path(combined_dir, "enzimi_data.tsv")
  write_tsv_safe(enzyme_tbl, combined_data_file)
  output_manifests$funz <- append_manifest_entry(
    output_manifests$funz,
    combined_data_file,
    "enzyme_data_tsv",
    "enzyme_insieme"
  )

  if ("bar" %in% enzyme_plot_types) {
    bar_files <- save_png_dimensions(
      make_enzyme_barplot(enzyme_tbl, "Enzyme barplot - combined"),
      combined_dir,
      "barplot_enzimi",
      dimensions,
      plot_dpi
    )
    for (output_file in unname(bar_files)) {
      output_manifests$funz <- append_manifest_entry(
        output_manifests$funz,
        output_file,
        "enzyme_barplot_png",
        "enzyme_insieme"
      )
    }
  }
  if ("line" %in% enzyme_plot_types) {
    line_files <- save_png_dimensions(
      make_enzyme_lineplot(enzyme_tbl, "Enzyme line chart - combined"),
      combined_dir,
      "lineplot_enzimi",
      dimensions,
      plot_dpi
    )
    for (output_file in unname(line_files)) {
      output_manifests$funz <- append_manifest_entry(
        output_manifests$funz,
        output_file,
        "enzyme_lineplot_png",
        "enzyme_insieme"
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
    write_tsv_safe(ec_tbl, ec_data_file)
    output_manifests$funz <- append_manifest_entry(
      output_manifests$funz,
      ec_data_file,
      "enzyme_data_tsv",
      "enzyme_separato",
      current_ec_code
    )

    if ("bar" %in% enzyme_plot_types) {
      bar_files <- save_png_dimensions(
        make_enzyme_barplot(ec_tbl, paste0("Enzyme barplot - EC ", current_ec_code)),
        ec_dir,
        "barplot_enzima",
        dimensions,
        plot_dpi
      )
      for (output_file in unname(bar_files)) {
        output_manifests$funz <- append_manifest_entry(
          output_manifests$funz,
          output_file,
          "enzyme_barplot_png",
          "enzyme_separato",
          current_ec_code
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
      for (output_file in unname(line_files)) {
        output_manifests$funz <- append_manifest_entry(
          output_manifests$funz,
          output_file,
          "enzyme_lineplot_png",
          "enzyme_separato",
          current_ec_code
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
    filtered_taxon_rank = NA_character_) {
  progress_message(
    "FLOW | pathway=", pathway_name,
    " | ranks=", paste(taxonomy_ranks, collapse = ","),
    " | samples=", paste(selected_samples, collapse = ",")
  )
  orf_long <- build_orf_long_table(pathway_sqm, selected_samples)
  if (nrow(orf_long) == 0L) {
    warning("No positive ORF data for pathway ", pathway_name, ".", call. = FALSE)
    return(output_manifests)
  }

  ko_lookup <- get_ko_name_lookup(pathway_sqm)

  for (rank in taxonomy_ranks) {
    progress_message("FLOW | pathway=", pathway_name, " | rank=", rank)
    rank_dir <- file.path(
      output_dir,
      "flowplot",
      pathway_selection_directory(pathway_selection),
      sanitize_name(pathway_name),
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

      data_file <- file.path(rank_dir, paste0("flowplot_", rank, "_", sample_name, "_data.tsv"))
      progress_message("FLOW | writing data: ", data_file)
      write_tsv_safe(flow_tbl, data_file)
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
          pathway_id = pathway_id
        )
      )

      if ("png" %in% flowplot_formats) {
        progress_message("FLOW | saving PNG | pathway=", pathway_name, " | rank=", rank, " | sample=", sample_name)
        plot_object <- make_flow_plot(flow_tbl, pathway_name, rank, sample_name)
        png_files <- save_png_dimensions(
          plot_object,
          rank_dir,
          paste0("flowplot_", rank, "_", sample_name),
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
              pathway_id = pathway_id
            )
          )
        }
      }

      if ("html" %in% flowplot_formats) {
        html_file <- file.path(rank_dir, paste0("flowplot_", rank, "_", sample_name, ".html"))
        progress_message("FLOW | saving HTML | pathway=", pathway_name, " | rank=", rank, " | sample=", sample_name)
        widget <- make_flow_sankey(flow_tbl, pathway_name, rank, sample_name)
        save_html_widget(widget, html_file)
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
            pathway_id = pathway_id
          )
        )
      }
    }
  }

  output_manifests
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
        sanitize_name(pathway_name),
        count,
        rank
      )
      }
      dir.create(rank_dir, recursive = TRUE, showWarnings = FALSE)

      use_pathway_tpm_percent <- scope_name != "taxonomy_global" && identical(count, "percent")
      if (use_pathway_tpm_percent) {
        plot_data <- build_pathway_taxonomy_percent_table(
          sqm_object = sqm_object,
          rank = rank,
          selected_samples = selected_samples,
          top_n_taxa = top_n_taxa,
          pathway_name = pathway_name
        )
        plot_object <- make_pathway_taxonomy_percent_plot(
          plot_tbl = plot_data,
          pathway_name = pathway_name,
          rank = rank,
          selected_samples = selected_samples
        )
      } else {
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
      }
      data_file <- if (scope_name == "taxonomy_global") {
        file.path(rank_dir, paste0("taxonomy_global_", count, "_", rank, "_data.tsv"))
      } else {
        file.path(rank_dir, paste0("taxonomy_", sanitize_name(pathway_name), "_", count, "_", rank, "_data.tsv"))
      }
      progress_message("TAXON | writing data: ", data_file)
      write_tsv_safe(plot_data, data_file)

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
          pathway_id = pathway_id
        )
      )

      file_stem <- if (scope_name == "taxonomy_global") {
        paste0("taxonomy_global_", count, "_", rank)
      } else {
        paste0("taxonomy_", sanitize_name(pathway_name), "_", count, "_", rank)
      }
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
            pathway_id = pathway_id
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

  for (pathview_sample_mode in pathview_sample_modes) {
    split_samples <- identical(pathview_sample_mode, "separato")
    pathview_dir <- file.path(
      output_dir,
      "pathview",
      pathway_selection_directory(pathway_selection),
      pathview_sample_mode,
      sanitize_name(pathway_name)
    )
    dir.create(pathview_dir, recursive = TRUE, showWarnings = FALSE)
    pathview_dir_abs <- normalizePath(pathview_dir, winslash = "/", mustWork = FALSE)
    before_files <- list.files(pathview_dir_abs, recursive = TRUE, full.names = TRUE, all.files = FALSE)

    progress_message(
      "PATHVIEW | mode=", pathview_sample_mode,
      " | exportPathway -> ", pathview_dir_abs
    )
    export_pathway_fn(
      SQM = sqm_object,
      pathway_id = pathway_id,
      count = "tpm",
      samples = selected_samples,
      split_samples = split_samples,
      output_dir = pathview_dir_abs,
      output_suffix = paste0("pathview_", pathway_id)
    )

    after_files <- list.files(pathview_dir_abs, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    new_files <- setdiff(after_files, before_files)
    if (length(new_files) == 0L) {
      new_files <- after_files
    }

    for (output_file in sort(new_files)) {
      output_manifests$pathview <- bind_rows(
        output_manifests$pathview,
        new_manifest_row(
          script_name = script_name,
          project_dir = project_dir,
          tax_mode = tax_mode,
          pathway = pathway_name,
          pathway_id = pathway_id,
          samples = selected_samples,
          metric = "tpm",
          top_n_taxa = top_n_taxa,
          top_n_ko = top_n_ko,
          output_type = "pathview_file",
          output_file = relative_to_output(output_file, manifest_base_dir),
          mode = "pathview",
          format = infer_format_from_path(output_file),
          output_scope = paste0("pathway_", pathway_selection, "_", pathview_sample_mode),
          filtered_taxon = filtered_taxon,
          filtered_taxon_rank = filtered_taxon_rank
        )
      )
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
    filtered_taxon_rank = NA_character_) {
  progress_message(
    "PIE | pathway=", pathway_name,
    " | ranks=", paste(taxonomy_ranks, collapse = ","),
    " | samples=", paste(selected_samples, collapse = ",")
  )
  pie_root <- file.path(
    output_dir,
    "pie",
    pathway_selection_directory(pathway_selection),
    sanitize_name(pathway_name)
  )
  dir.create(pie_root, recursive = TRUE, showWarnings = FALSE)

  orf_long <- build_orf_long_table(pathway_sqm, selected_samples)
  if (nrow(orf_long) == 0L) {
    warning("No positive ORF data for pathway ", pathway_name, ".", call. = FALSE)
    return(output_manifests)
  }

  ko_names <- orf_long |>
    distinct(.data$ko_id, .data$kegg_function)
  ko_ec_lookup <- extract_ko_ec_lookup(orf_long)
  pathway_sample_totals <- orf_long |>
    group_by(.data$sample) |>
    summarise(pathway_sample_tpm = sum(.data$tpm), .groups = "drop")

  sample_names <- orf_long |>
    distinct(.data$sample) |>
    pull(.data$sample) |>
    sort()
  ko_ids <- orf_long |>
    distinct(.data$ko_id) |>
    pull(.data$ko_id) |>
    sort()

  for (sample_name in sample_names) {
    progress_message("PIE | pathway=", pathway_name, " | sample=", sample_name)
    sample_dir <- file.path(pie_root, sanitize_name(sample_name))
    for (ko_id_value in ko_ids) {
      ko_name_row <- ko_names |>
        filter(.data$ko_id == ko_id_value)
      ko_name <- if (nrow(ko_name_row) > 0L && !is.na(ko_name_row$kegg_function[[1]])) {
        ko_name_row$kegg_function[[1]]
      } else {
        ko_id_value
      }
      ko_ec_row <- ko_ec_lookup |>
        filter(.data$ko_id == ko_id_value)
      ko_ec <- if (nrow(ko_ec_row) > 0L) ko_ec_row$ec_codes[[1]] else NA_character_
      ko_dir <- file.path(sample_dir, sanitize_name(get_ko_dir_name(ko_id_value, ko_ec)))
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
          top_n_taxa = top_n_taxa
        )

        if (nrow(plot_data) == 0L || sum(plot_data$tpm) <= 0) {
          next
        }

        data_file_stem <- paste0(
          rank_name, "_", sanitize_name(ko_id_value), "_", sanitize_name(coalesce(ko_ec, "NA")), "_", sanitize_name(sample_name)
        )
        data_file <- file.path(ko_dir, paste0(data_file_stem, "_data.tsv"))
        progress_message("PIE | writing data: ", data_file)
        write_tsv_safe(plot_data, data_file)

        output_manifests$pie <- bind_rows(
          output_manifests$pie,
          new_manifest_row(
            script_name = script_name,
            project_dir = project_dir,
            tax_mode = tax_mode,
            pathway = pathway_name,
            ko_id = ko_id_value,
            pathway_id = pathway_id,
            samples = sample_name,
            metric = "tpm",
            top_n_taxa = top_n_taxa,
            top_n_ko = top_n_ko,
            output_type = "data_tsv",
            output_file = relative_to_output(data_file, manifest_base_dir),
            mode = "pie",
            rank = rank_name,
            format = "tsv",
            dpi = plot_dpi,
            output_scope = paste0("pathway_", pathway_selection),
            filtered_taxon = filtered_taxon,
            filtered_taxon_rank = filtered_taxon_rank
          )
        )

        plot_object <- make_pie_plot(
          plot_data = plot_data,
          pathway_name = pathway_name,
          sample_name = sample_name,
          ko_id = ko_id_value,
          ko_name = ko_name,
          rank_name = rank_name,
          ko_ec = ko_ec,
          pathway_sample_tpm = pathway_sample_tpm
        )
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
              ko_id = ko_id_value,
              pathway_id = pathway_id,
              samples = sample_name,
              metric = "tpm",
              top_n_taxa = top_n_taxa,
              top_n_ko = top_n_ko,
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
              filtered_taxon_rank = filtered_taxon_rank
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

merge_section_manifest <- function(new_manifest_tbl, existing_manifest_tbl) {
  bind_rows(new_manifest_tbl, existing_manifest_tbl) |>
    distinct(.data$output_file, .keep_all = TRUE)
}

section_manifest_paths <- function(output_dir) {
  paths <- c(
    flow = file.path(output_dir, "flowplot", "manifest_flow.tsv"),
    funz = file.path(output_dir, "funz", "manifest_funz.tsv"),
    taxon = file.path(output_dir, "manifest_taxon.tsv"),
    pathview = file.path(output_dir, "pathview", "manifest_pathview.tsv"),
    pie = file.path(output_dir, "pie", "manifest_pie.tsv")
  )
  paths[file.exists(paths)]
}

write_section_manifest <- function(manifest_tbl, output_dir, relative_dir, filename) {
  dir_path <- file.path(output_dir, relative_dir)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(dir_path, filename)
  existing_manifest_tbl <- read_section_manifest(manifest_path, relative_dir)
  merged_manifest_tbl <- merge_section_manifest(manifest_tbl, existing_manifest_tbl)
  write_tsv_safe(merged_manifest_tbl, manifest_path)
  manifest_path
}

# ---- Command-line entry point ---------------------------------------------

# Validate inputs before loading SQM, then run only the requested analysis modes.
main <- function() {
  check_required_packages()
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

  tax_mode <- if (is.null(named_args$tax_mode)) "prokfilter" else named_args$tax_mode
  top_n_ko <- if (is.null(named_args$top_n_ko)) 20L else as.integer(named_args$top_n_ko)
  top_n_taxa <- if (is.null(named_args$top_n_taxa)) 15L else as.integer(named_args$top_n_taxa)
  pathway_top_n <- if (is.null(named_args$pathway_top_n)) default_pathway_top_n else as.integer(named_args$pathway_top_n)
  plot_dpi <- if (is.null(named_args$plot_dpi)) 600 else as.numeric(named_args$plot_dpi)
  validate_positive_integer(top_n_ko, "top_n_ko")
  validate_positive_integer(top_n_taxa, "top_n_taxa")
  validate_positive_integer(pathway_top_n, "pathway_top_n")
  if (is.na(plot_dpi) || plot_dpi <= 0) {
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

  taxonomy_counts <- if (is.null(named_args$taxonomy_counts)) {
    c("abund", "percent")
  } else {
    split_csv_arg(named_args$taxonomy_counts)
  }
  if (!all(taxonomy_counts %in% c("abund", "percent"))) {
    stop("taxonomy_counts must contain only abund and/or percent.", call. = FALSE)
  }

  flowplot_formats <- if (is.null(named_args$flowplot_formats)) {
    c("png", "html")
  } else {
    split_csv_arg(named_args$flowplot_formats)
  }
  if (!all(flowplot_formats %in% c("png", "html"))) {
    stop("flowplot_formats must contain only png and/or html.", call. = FALSE)
  }
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

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
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
  selected_samples <- if (is.null(named_args$samples)) available_samples else split_csv_arg(named_args$samples)
  if (length(selected_samples) == 0L) {
    stop("samples cannot be empty.", call. = FALSE)
  }
  validate_samples(selected_samples, available_samples)
  progress_message("Selected samples: ", paste(selected_samples, collapse = ", "))

  pathway_groups <- if (mode == "enzimi") {
    list()
  } else {
    select_pathway_groups(
      sqm = sqm,
      requested_pathways = requested_pathways,
      pathway_selection_modes = pathway_selection_modes,
      selected_samples = selected_samples,
      pathway_top_n = pathway_top_n
    )
  }
  if (length(pathway_groups) > 0L) {
    progress_message(
      "Resolved pathways: ",
      paste(
        unlist(map(pathway_groups, ~ vapply(.x, `[[`, character(1), "canonical_pathway_name"))),
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
  }) |> purrr::flatten()

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
      list(
        sqm = filtered_sqm,
        output_dir = file.path(
          output_dir,
          "taxon_filter",
          taxon_info$rank,
          sanitize_name(taxon_info$taxon)
        ),
        filtered_taxon = taxon_info$taxon,
        filtered_taxon_rank = taxon_info$rank
      )
    })
  }

  manifests <- list(
    flow = tibble(),
    funz = tibble(),
    taxon = tibble(),
    pathview = tibble(),
    pie = tibble()
  )
  script_name <- "sqm_plots.R"

  for (context in filter_contexts) {
    context_output_dir <- context$output_dir
    context_sqm <- context$sqm
    progress_message(
      "Output context: ", context_output_dir,
      if (!is.na(context$filtered_taxon)) paste0(" | taxon filter=", context$filtered_taxon, " @ ", context$filtered_taxon_rank) else ""
    )
    pathway_sqms <- map(
      pathway_entries,
      function(pathway_info) {
        list(
          pathway_name = pathway_info$pathway_name,
          pathway_id = pathway_info$pathway_id,
          pathway_selection = pathway_info$pathway_selection,
          pathway_sqm = subset_pathway(context_sqm, pathway_info$pathway_name)
        )
      }
    )
    names(pathway_sqms) <- vapply(
      pathway_entries,
      function(pathway_info) paste(pathway_info$pathway_selection, pathway_info$pathway_name, sep = "::"),
      character(1)
    )

    if (mode %in% c("all", "funz")) {
      progress_message("Starting FUNZ section")
    for (pathway_key in names(pathway_sqms)) {
      pathway_info <- pathway_sqms[[pathway_key]]
      pathway_name <- pathway_info$pathway_name
      manifests <- run_funz_mode(
          output_dir = context_output_dir,
          manifest_base_dir = output_dir,
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
          filtered_taxon_rank = context$filtered_taxon_rank
        )
      }
    }

    if (mode %in% c("all", "funz", "enzimi")) {
      progress_message("Starting ENZIMI section")
      manifests <- run_enzyme_mode(
        sqm_object = context_sqm,
        output_dir = context_output_dir,
        manifest_base_dir = output_dir,
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
        filtered_taxon = context$filtered_taxon,
        filtered_taxon_rank = context$filtered_taxon_rank
      )
    }

    if (mode %in% c("all", "flow")) {
      progress_message("Starting FLOW section")
    for (pathway_key in names(pathway_sqms)) {
      pathway_info <- pathway_sqms[[pathway_key]]
      pathway_name <- pathway_info$pathway_name
      manifests <- run_flow_mode(
          output_dir = context_output_dir,
          manifest_base_dir = output_dir,
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
          filtered_taxon_rank = context$filtered_taxon_rank
        )
      }
    }

    if (mode %in% c("all", "taxon")) {
      progress_message("Starting TAXON section")
      manifests <- run_taxonomy_scope(
        sqm_object = context_sqm,
        output_dir = context_output_dir,
        manifest_base_dir = output_dir,
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

      if (length(pathway_sqms) > 0L) {
      for (pathway_key in names(pathway_sqms)) {
        pathway_info <- pathway_sqms[[pathway_key]]
        pathway_name <- pathway_info$pathway_name
        manifests <- run_taxonomy_scope(
          sqm_object = pathway_info$pathway_sqm,
            output_dir = context_output_dir,
            manifest_base_dir = output_dir,
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
      }
    }

    if (mode %in% c("all", "pathview")) {
      progress_message("Starting PATHVIEW section")
    for (pathway_key in names(pathway_sqms)) {
      pathway_info <- pathway_sqms[[pathway_key]]
      pathway_name <- pathway_info$pathway_name
      if (!pathview_is_exportable(pathway_info$pathway_selection, pathway_info$pathway_id)) {
        warning(
          "Skipping PATHVIEW: pathway has no resolvable KEGG ID: ",
          pathway_name,
          call. = FALSE
        )
        next
      }
      manifests <- run_pathview_mode(
          sqm_object = context_sqm,
          output_dir = context_output_dir,
          manifest_base_dir = output_dir,
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

    if (mode %in% c("all", "pie")) {
      progress_message("Starting PIE section")
    for (pathway_key in names(pathway_sqms)) {
      pathway_info <- pathway_sqms[[pathway_key]]
      if (!pathway_info$pathway_selection %in% pie_selection_modes) {
        next
      }
      pathway_name <- pathway_info$pathway_name
      manifests <- run_pie_mode(
          output_dir = context_output_dir,
          manifest_base_dir = output_dir,
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
          filtered_taxon_rank = context$filtered_taxon_rank
        )
      }
    }
  }

  manifest_paths <- character()
  if (nrow(manifests$flow) > 0L) {
    manifest_paths[["flow"]] <- write_section_manifest(manifests$flow, output_dir, "flowplot", "manifest_flow.tsv")
  }
  if (nrow(manifests$funz) > 0L) {
    manifest_paths[["funz"]] <- write_section_manifest(manifests$funz, output_dir, "funz", "manifest_funz.tsv")
  }
  if (nrow(manifests$taxon) > 0L) {
    manifest_paths[["taxon"]] <- write_section_manifest(manifests$taxon, output_dir, "", "manifest_taxon.tsv")
  }
  if (nrow(manifests$pathview) > 0L) {
    manifest_paths[["pathview"]] <- write_section_manifest(manifests$pathview, output_dir, "pathview", "manifest_pathview.tsv")
  }
  if (nrow(manifests$pie) > 0L) {
    manifest_paths[["pie"]] <- write_section_manifest(manifests$pie, output_dir, "pie", "manifest_pie.tsv")
  }

  manifest_paths <- section_manifest_paths(output_dir)
  manifest_all <- tibble::tibble(
    section = names(manifest_paths),
    manifest_file = vapply(manifest_paths, relative_to_output, character(1), output_dir = output_dir)
  )
  manifest_all_path <- file.path(output_dir, "manifest_all.tsv")
  write_tsv_safe(manifest_all, manifest_all_path)

  message("Output directory: ", output_dir)
  message("Combined manifest: ", manifest_all_path)
  invisible(0)
}

main()
