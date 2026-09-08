source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

expect_true <- function(condition, label) {
  if (!isTRUE(condition)) {
    stop(label, call. = FALSE)
  }
}

failures <- character()
run_case <- function(name, code) {
  tryCatch(
    {
      force(code)
      message("PASS: ", name)
    },
    error = function(error) {
      failures <<- c(failures, paste0(name, ": ", conditionMessage(error)))
      message("FAIL: ", name, " -- ", conditionMessage(error))
    }
  )
}

script_env <- source_without_main("sqm_plots.R")

run_case("raw positive-integer parser accepts only canonical decimal digits", {
  expect_true(
    exists(
      "parse_positive_integer_arg",
      envir = script_env,
      mode = "function",
      inherits = FALSE
    ),
    "Required CLI boundary helper is absent: parse_positive_integer_arg"
  )
  parse_integer <- get(
    "parse_positive_integer_arg",
    envir = script_env,
    mode = "function",
    inherits = FALSE
  )

  accepted <- c("1", "0007", as.character(.Machine$integer.max))
  expected <- c(1L, 7L, .Machine$integer.max)
  observed <- vapply(
    accepted,
    function(value) parse_integer(value, "top_n_ko"),
    integer(1)
  )
  expect_true(
    identical(unname(observed), expected),
    "The raw CLI integer parser rejected a valid boundary or changed its value"
  )

  rejected <- c(
    "1.9", "2.5", "3.1", "1e2", "+1", "-1", " 1", "1 ",
    "0", "000", "", "2147483648", "999999999999999999999999999"
  )
  for (value in rejected) {
    parse_error <- tryCatch(
      {
        parse_integer(value, "top_n_ko")
        NULL
      },
      error = function(error) error
    )
    display_value <- if (nzchar(value)) value else "<empty>"
    expect_true(
      inherits(parse_error, "error"),
      paste0("Invalid raw integer was accepted: ", display_value)
    )
    expect_true(
      grepl("top_n_ko", conditionMessage(parse_error), fixed = TRUE),
      paste0("Invalid integer error omitted its argument name: ", display_value)
    )
  }
})

run_case("invalid integer CLI values fail before SQM loading and write a failed-run manifest", {
  test_root <- tempfile("p2_cli_integer_")
  project_dir <- file.path(test_root, "synthetic_project")
  dir.create(project_dir, recursive = TRUE)
  on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

  script_path <- normalizePath("sqm_plots.R", winslash = "/", mustWork = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  if (!file.exists(rscript)) {
    rscript <- file.path(R.home("bin"), "Rscript")
  }

  invalid_options <- c(
    top_n_ko = "1.9",
    top_n_taxa = "2.5",
    pathway_top_n = "3.1"
  )
  for (option_name in names(invalid_options)) {
    output_dir <- file.path(test_root, paste0("failed_run_", option_name))
    cli_output <- suppressWarnings(system2(
      rscript,
      args = c(
        shQuote(script_path),
        shQuote(paste0("--project_dir=", project_dir)),
        shQuote(paste0("--output_dir=", output_dir)),
        "--mode=funz",
        paste0("--", option_name, "=", invalid_options[[option_name]])
      ),
      stdout = TRUE,
      stderr = TRUE
    ))
    exit_status <- attr(cli_output, "status")
    if (is.null(exit_status)) {
      exit_status <- 0L
    }

    expect_true(
      !identical(as.integer(exit_status), 0L),
      paste0("Invalid --", option_name, " unexpectedly exited successfully")
    )
    expect_true(
      dir.exists(output_dir),
      paste0("Invalid --", option_name, " did not initialize its failed-run output")
    )
    expect_true(
      !any(grepl("Loading SQM project", cli_output, fixed = TRUE)),
      paste0("Invalid --", option_name, " reached loadSQM before failing")
    )
    expect_true(
      any(grepl(option_name, cli_output, fixed = TRUE)),
      paste0("Invalid --", option_name, " error did not identify the option")
    )
    run_logs <- list.files(output_dir, pattern = "^[0-9]{8}T.*\\.log$", full.names = TRUE)
    expect_true(
      length(run_logs) == 1L,
      paste0("Invalid --", option_name, " did not write one isolated failure log")
    )
    run_log <- paste(readLines(run_logs[[1L]], warn = FALSE), collapse = "\n")
    expect_true(
      grepl("STATUS=FAILED", run_log, fixed = TRUE) &&
        grepl(option_name, run_log, fixed = TRUE),
      paste0("Invalid --", option_name, " wrote incorrect failure log metadata")
    )
  }
})

if (length(failures) > 0L) {
  stop(
    paste(
      "P2 CLI integer regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: P2 CLI integer parsing validates raw input before side effects")
