# Diagnostic replica of the FLOW KO rescaling calculation in sqm_plots_lean.R.

TOLERANCE <- 1e-8
PROJECT_DIR <- file.path("in", "Au_sip")
OUTPUT_FILE <- file.path("out", "flow_scale_factors.tsv")

extract_ko_ids <- function(x) {
  if (length(x) != 1L || is.na(x)) return(character())
  hits <- regmatches(as.character(x), gregexpr("K[0-9]{5}", as.character(x), perl = TRUE))[[1L]]
  if (identical(hits, "")) character() else unique(hits)
}

pathway_names <- function(values) {
  roots <- c(
    "Metabolism", "Genetic Information Processing", "Environmental Information Processing",
    "Cellular Processes", "Organismal Systems", "Human Diseases"
  )
  entries <- strsplit(as.character(values), "\\s*\\|\\s*", perl = TRUE)
  parts <- strsplit(unlist(entries, use.names = FALSE), "\\s*;\\s*", perl = TRUE)
  valid <- lengths(parts) == 3L
  fields <- do.call(rbind, parts[valid])
  keep <- fields[, 1L] %in% roots & !grepl("not included", fields[, 2L], ignore.case = TRUE)
  sort(unique(sub("\\s*-\\s*Reference pathway\\s*$", "", fields[keep, 3L], ignore.case = TRUE)))
}

make_factor_table <- function(pathway_name, targets, observed, tolerance = TOLERANCE) {
  target_key <- paste(targets$sample, targets$ko_id, sep = "\r")
  observed_key <- paste(observed$sample, observed$ko_id, sep = "\r")
  raw_tpm <- observed$raw_tpm[match(target_key, observed_key)]
  raw_tpm[is.na(raw_tpm)] <- 0

  bad <- which(targets$sqm_tpm > tolerance & raw_tpm <= 0)
  if (length(bad)) {
    row <- bad[[1L]]
    stop(
      "Cannot allocate positive SQM KO TPM for ", pathway_name, "/",
      targets$sample[[row]], "/", targets$ko_id[[row]], ".",
      call. = FALSE
    )
  }

  factor <- rep(NA_real_, nrow(targets))
  defined <- targets$sqm_tpm > 0 & raw_tpm > 0
  factor[defined] <- targets$sqm_tpm[defined] / raw_tpm[defined]
  status <- ifelse(
    defined, "scaled",
    ifelse(targets$sqm_tpm <= 0 & raw_tpm <= 0, "both_zero",
           ifelse(targets$sqm_tpm <= 0, "target_zero", "raw_zero"))
  )

  data.frame(
    pathway_name = pathway_name,
    sample = targets$sample,
    ko_id = targets$ko_id,
    raw_tpm = raw_tpm,
    sqm_tpm = targets$sqm_tpm,
    factor = factor,
    difference = targets$sqm_tpm - raw_tpm,
    status = status,
    stringsAsFactors = FALSE
  )
}

check_pathway <- function(sqm, pathway_name, samples) {
  pathway_sqm <- SQMtools::subsetFun(
    SQM = sqm, fun = pathway_name, columns = "KEGGPATH",
    ignore_case = FALSE, fixed = TRUE, allow_empty = FALSE
  )
  orf_table <- as.data.frame(pathway_sqm$orfs$table, check.names = FALSE)
  orf_tpm <- as.data.frame(pathway_sqm$orfs$tpm, check.names = FALSE)
  official <- as.data.frame(pathway_sqm$functions$KEGG$tpm, check.names = FALSE)
  ko_ids <- rownames(official)

  all_kos <- lapply(as.character(orf_table[, "KEGG ID"]), extract_ko_ids)
  selected <- lapply(all_kos, intersect, ko_ids)
  source <- rep.int(seq_along(selected), lengths(selected))
  selected_ko <- unlist(selected, use.names = FALSE)

  allocated <- do.call(rbind, lapply(samples, function(sample) {
    share <- as.numeric(orf_tpm[rownames(orf_table)[source], sample]) / lengths(all_kos)[source]
    keep <- share > 0
    data.frame(
      sample = rep(sample, sum(keep)), ko_id = selected_ko[keep], raw_tpm = share[keep],
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(allocated)) {
    observed <- aggregate(allocated$raw_tpm, allocated[c("sample", "ko_id")], sum)
    names(observed)[[3L]] <- "raw_tpm"
  } else {
    observed <- data.frame(sample = character(), ko_id = character(), raw_tpm = numeric())
  }

  targets <- expand.grid(sample = samples, ko_id = ko_ids,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  targets$sqm_tpm <- mapply(
    function(sample, ko_id) as.numeric(official[ko_id, sample]),
    targets$sample, targets$ko_id
  )
  make_factor_table(pathway_name, targets, observed)
}

main <- function() {
  if (!dir.exists(PROJECT_DIR)) stop("Project directory does not exist: ", PROJECT_DIR, call. = FALSE)
  if (!requireNamespace("SQMtools", quietly = TRUE)) stop("Missing package: SQMtools", call. = FALSE)
  suppressPackageStartupMessages(library("SQMtools", character.only = TRUE))
  sqm <- SQMtools::loadSQM(
    project_path = normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE),
    tax_mode = "prokfilter", trusted_functions_only = FALSE, load_sequences = FALSE
  )
  samples <- colnames(sqm$orfs$tpm)
  pathways <- pathway_names(sqm$orfs$table[, "KEGGPATH"])
  cat("Pathway da controllare:", length(pathways), "\n")

  results <- vector("list", length(pathways))
  for (index in seq_along(pathways)) {
    results[[index]] <- check_pathway(sqm, pathways[[index]], samples)
    if (index %% 25L == 0L || index == length(pathways)) {
      cat("Completati:", index, "/", length(pathways), "\n")
    }
  }
  result <- do.call(rbind, results)
  rownames(result) <- NULL
  dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(result, OUTPUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  factors <- result$factor[!is.na(result$factor)]
  cat("Fattori definiti:", length(factors), "\n")
  cat("Minimo:", format(min(factors), digits = 17), "\n")
  cat("Massimo:", format(max(factors), digits = 17), "\n")
  cat("Diversi da 1 oltre", TOLERANCE, ":", sum(abs(factors - 1) > TOLERANCE), "\n")
  cat("Output:", OUTPUT_FILE, "\n")
  invisible(result)
}

if (identical(environment(), globalenv()) && sys.nframe() == 0L) main()
