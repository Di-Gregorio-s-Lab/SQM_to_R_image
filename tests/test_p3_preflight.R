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

expect_package_set <- function(observed, expected, label) {
  expect_true(
    length(observed) == length(expected) && setequal(observed, expected),
    paste0(
      label,
      "\nExpected: ", paste(expected, collapse = ", "),
      "\nObserved: ", paste(observed, collapse = ", ")
    )
  )
}

run_isolated_rscript <- function(rscript, args, user_library, site_library) {
  old_user <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  old_site <- Sys.getenv("R_LIBS_SITE", unset = NA_character_)
  on.exit({
    if (is.na(old_user)) Sys.unsetenv("R_LIBS_USER") else Sys.setenv(R_LIBS_USER = old_user)
    if (is.na(old_site)) Sys.unsetenv("R_LIBS_SITE") else Sys.setenv(R_LIBS_SITE = old_site)
  }, add = TRUE)

  Sys.setenv(R_LIBS_USER = user_library, R_LIBS_SITE = site_library)
  suppressWarnings(system2(
    command = rscript,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  ))
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

run_case("required packages are selected from mode and FLOW format", {
  expect_true(
    exists(
      "required_packages_for_mode",
      envir = script_env,
      mode = "function",
      inherits = FALSE
    ),
    "Required P3 helper is absent: required_packages_for_mode"
  )
  required_packages <- get(
    "required_packages_for_mode",
    envir = script_env,
    mode = "function",
    inherits = FALSE
  )

  common <- c(
    "SQMtools", "readr", "dplyr", "tidyr", "tibble", "stringr",
    "ggplot2", "glue", "purrr", "scales"
  )

  expect_package_set(
    required_packages("flow", "png"),
    c(common, "ggalluvial"),
    "FLOW PNG package selection is not minimal"
  )
  expect_package_set(
    required_packages("flow", "html"),
    c(common, "ggalluvial", "plotly", "htmlwidgets", "rmarkdown"),
    "FLOW HTML package selection is incomplete"
  )
  expect_package_set(
    required_packages("pie", "png"),
    c(common, "forcats", "rlang"),
    "PIE package selection is incomplete"
  )
  expect_package_set(
    required_packages("pathview", "png"),
    c(common, "pathview"),
    "Pathview package selection is incomplete"
  )
  expect_package_set(
    required_packages("all", c("png", "html")),
    c(
      common, "ggalluvial", "plotly", "htmlwidgets", "rmarkdown", "forcats",
      "rlang", "pathview"
    ),
    "All-mode package selection is incomplete"
  )

  all_selected <- required_packages("all", c("png", "html"))
  expect_true(
    !any(c("ggpattern", "magick") %in% all_selected),
    "Optional SQMtools packages became mandatory"
  )
})

run_case("FLOW HTML preflight requires Pandoc", {
  pandoc_error <- tryCatch(
    {
      script_env$check_flow_html_preflight(
        "flow",
        "html",
        pandoc_available_fn = function() FALSE
      )
      NA_character_
    },
    error = function(error) conditionMessage(error)
  )
  expect_true(
    !is.na(pandoc_error) && grepl("Pandoc", pandoc_error, fixed = TRUE),
    "FLOW HTML accepted an environment without Pandoc"
  )
  expect_true(
    isTRUE(script_env$check_flow_html_preflight(
      "flow",
      "png",
      pandoc_available_fn = function() FALSE
    )),
    "PNG-only FLOW unnecessarily required Pandoc"
  )
})

run_case("help and no-argument CLI work with no optional R libraries", {
  test_root <- tempfile("p3_preflight_help_")
  empty_user_library <- file.path(test_root, "empty-user-library")
  empty_site_library <- file.path(test_root, "empty-site-library")
  dir.create(empty_user_library, recursive = TRUE)
  dir.create(empty_site_library, recursive = TRUE)
  on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

  script_path <- normalizePath("sqm_plots.R", winslash = "/", mustWork = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  if (!file.exists(rscript)) {
    rscript <- file.path(R.home("bin"), "Rscript")
  }
  if (.Platform$OS.type == "windows") {
    rscript <- utils::shortPathName(rscript)
  }
  isolated_user_library <- normalizePath(empty_user_library, winslash = "/", mustWork = TRUE)
  isolated_site_library <- normalizePath(empty_site_library, winslash = "/", mustWork = TRUE)

  invocations <- list(
    no_arguments = character(),
    help = "--help"
  )
  for (invocation_name in names(invocations)) {
    output <- run_isolated_rscript(
      rscript = rscript,
      args = c("--vanilla", shQuote(script_path), invocations[[invocation_name]]),
      user_library = isolated_user_library,
      site_library = isolated_site_library
    )
    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0L
    }

    expect_true(
      identical(as.integer(status), 0L),
      paste0(
        "CLI ", invocation_name,
        " failed before showing help:\n", paste(output, collapse = "\n")
      )
    )
    expect_true(
      any(grepl("Usage:", output, fixed = TRUE)),
      paste0("CLI ", invocation_name, " did not show Usage")
    )
    expect_true(
      !any(grepl("there is no package called", output, fixed = TRUE)),
      paste0("CLI ", invocation_name, " tried to load an optional package")
    )
    expect_true(
      !any(grepl("Missing required R packages", output, fixed = TRUE)),
      paste0("CLI ", invocation_name, " ran dependency preflight before help")
    )
  }
})

run_case("analytical CLI reports all missing packages before side effects", {
  test_root <- tempfile("p3_preflight_run_")
  empty_user_library <- file.path(test_root, "empty-user-library")
  empty_site_library <- file.path(test_root, "empty-site-library")
  project_dir <- file.path(test_root, "synthetic-project")
  output_dir <- file.path(test_root, "must-not-exist")
  dir.create(empty_user_library, recursive = TRUE)
  dir.create(empty_site_library, recursive = TRUE)
  dir.create(project_dir, recursive = TRUE)
  on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

  script_path <- normalizePath("sqm_plots.R", winslash = "/", mustWork = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  if (!file.exists(rscript)) {
    rscript <- file.path(R.home("bin"), "Rscript")
  }
  if (.Platform$OS.type == "windows") {
    rscript <- utils::shortPathName(rscript)
  }
  isolated_user_library <- normalizePath(empty_user_library, winslash = "/", mustWork = TRUE)
  isolated_site_library <- normalizePath(empty_site_library, winslash = "/", mustWork = TRUE)

  output <- run_isolated_rscript(
    rscript = rscript,
    args = c(
      "--vanilla",
      shQuote(script_path),
      shQuote(paste0("--project_dir=", project_dir)),
      shQuote(paste0("--output_dir=", output_dir)),
      "--mode=flow",
      "--flowplot_formats=png"
    ),
    user_library = isolated_user_library,
    site_library = isolated_site_library
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  output_text <- paste(output, collapse = "\n")

  expect_true(
    !identical(as.integer(status), 0L),
    "Analytical CLI unexpectedly succeeded without required packages"
  )
  expect_true(
    sum(grepl("Missing required R packages", output, fixed = TRUE)) == 1L,
    paste0(
      "Missing-package failure was not reported as one controlled error:\n",
      output_text
    )
  )
  expect_true(
    grepl("mode.*flow|flow.*mode", output_text, ignore.case = TRUE),
    paste0("Missing-package error did not identify mode flow:\n", output_text)
  )
  expect_true(
    !dir.exists(output_dir),
    "Dependency preflight created output_dir before failing"
  )
  expect_true(
    !grepl("Loading SQM project", output_text, fixed = TRUE),
    "Dependency preflight reached loadSQM before failing"
  )
  expect_true(
    !grepl("there is no package called", output_text, fixed = TRUE),
    paste0("A library() call failed before controlled preflight:\n", output_text)
  )
})

if (length(failures) > 0L) {
  stop(
    paste(
      "P3 dependency preflight regressions failed:",
      paste0("- ", failures, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: P3 dependency preflight is mode-aware and side-effect free")
