test_files <- list.files(
  "tests",
  pattern = "^test_.*\\.R$",
  full.names = TRUE
)

excluded_tests <- c(
  "test_p0_integration_Au_sip.R",
  "test_p1_integration_Au_sip.R",
  "test_p2_integration_Au_sip.R",
  "test_p3_integration_Au_sip.R",
  "test_t1_integration_CS8.R",
  "test_run_contract_integration_Au_sip.R"
)
test_files <- test_files[!basename(test_files) %in% excluded_tests]
test_files <- sort(test_files)

if (length(test_files) == 0L) {
  stop("No fast R tests were found.", call. = FALSE)
}

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) {
  rscript <- file.path(R.home("bin"), "Rscript")
}

failed_tests <- character()
for (test_file in test_files) {
  message("[fast-tests] Running ", test_file)
  output <- system2(
    command = rscript,
    args = shQuote(normalizePath(test_file, winslash = "/", mustWork = TRUE)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }

  if (!identical(as.integer(status), 0L)) {
    failed_tests <- c(failed_tests, test_file)
    writeLines(output)
  }
}

if (length(failed_tests) > 0L) {
  stop(
    "Fast tests failed: ", paste(failed_tests, collapse = ", "),
    call. = FALSE
  )
}

message("PASS: ", length(test_files), " fast R test files completed")
