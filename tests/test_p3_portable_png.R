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

expect_error_matching <- function(code, patterns, label) {
  captured <- NULL
  tryCatch(
    force(code),
    error = function(error) {
      captured <<- conditionMessage(error)
    }
  )

  if (is.null(captured)) {
    stop(label, ": expected an error", call. = FALSE)
  }
  for (pattern in patterns) {
    if (!grepl(pattern, captured, ignore.case = TRUE, perl = TRUE)) {
      stop(
        label, ": error did not match '", pattern, "': ", captured,
        call. = FALSE
      )
    }
  }
  invisible(captured)
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

test_root <- tempfile("p3 portable png ")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

plot_object <- ggplot2::ggplot(
  data.frame(x = 1, y = 1),
  ggplot2::aes(x = x, y = y)
) + ggplot2::geom_point()
dimensions <- list(`2x2` = c(width = 2, height = 2))

run_case("a short PNG path keeps its historical physical filename", {
  short_dir <- file.path(test_root, "short")
  expected_path <- file.path(short_dir, "short_plot_2x2.png")

  output_files <- script_env$save_png_dimensions(
    plot_object = plot_object,
    output_dir = short_dir,
    file_stem = "short_plot",
    dimensions = dimensions,
    dpi = 72,
    max_path_length = 1000L
  )

  expect_true(length(output_files) == 1L, "Expected exactly one short-path PNG")
  expect_true(
    identical(
      normalizePath(unname(output_files[[1]]), winslash = "/", mustWork = TRUE),
      normalizePath(expected_path, winslash = "/", mustWork = TRUE)
    ),
    "A path below the limit must keep the historical filename"
  )
  expect_true(file.info(expected_path)$size > 0, "The short-path PNG is empty")
})

run_case("long PNG names compact deterministically within an injected limit", {
  compact_dir <- file.path(test_root, "compact")
  dir.create(compact_dir, recursive = TRUE)
  compact_dir_abs <- normalizePath(compact_dir, winslash = "/", mustWork = TRUE)
  max_path_length <- nchar(compact_dir_abs, type = "chars") + 50L
  common_prefix <- paste(rep("readable", 20L), collapse = "_")
  first_stem <- paste0(common_prefix, "_alpha")
  second_stem <- paste0(common_prefix, "_beta")

  first_files <- script_env$save_png_dimensions(
    plot_object = plot_object,
    output_dir = compact_dir,
    file_stem = first_stem,
    dimensions = dimensions,
    dpi = 72,
    max_path_length = max_path_length
  )
  repeated_files <- script_env$save_png_dimensions(
    plot_object = plot_object,
    output_dir = compact_dir,
    file_stem = first_stem,
    dimensions = dimensions,
    dpi = 72,
    max_path_length = max_path_length
  )
  second_files <- script_env$save_png_dimensions(
    plot_object = plot_object,
    output_dir = compact_dir,
    file_stem = second_stem,
    dimensions = dimensions,
    dpi = 72,
    max_path_length = max_path_length
  )

  first_path <- normalizePath(
    unname(first_files[[1]]), winslash = "/", mustWork = TRUE
  )
  repeated_path <- normalizePath(
    unname(repeated_files[[1]]), winslash = "/", mustWork = TRUE
  )
  second_path <- normalizePath(
    unname(second_files[[1]]), winslash = "/", mustWork = TRUE
  )
  first_name <- basename(first_path)
  second_name <- basename(second_path)
  compact_pattern <- "^.+__[[:xdigit:]]{12}_2x2\\.png$"

  expect_true(
    identical(first_path, repeated_path),
    "The same logical PNG name must produce the same physical path"
  )
  expect_true(
    grepl(compact_pattern, first_name),
    paste0("The compact filename lacks a readable prefix or 12-hex token: ", first_name)
  )
  expect_true(
    grepl(compact_pattern, second_name),
    paste0("The second compact filename has the wrong shape: ", second_name)
  )
  expect_true(
    !identical(first_name, second_name),
    "Different logical stems must not collide after compaction"
  )
  expect_true(
    nchar(first_path, type = "chars") <= max_path_length &&
      nchar(second_path, type = "chars") <= max_path_length,
    "A compact physical path exceeds the configured character limit"
  )
  expect_true(
    all(file.info(c(first_path, second_path))$size > 0),
    "A returned compact PNG is missing or empty"
  )

  physical_files <- normalizePath(
    list.files(compact_dir, full.names = TRUE),
    winslash = "/",
    mustWork = TRUE
  )
  expect_true(
    identical(sort(physical_files), sort(c(first_path, second_path))),
    "The filesystem does not contain exactly the returned physical PNG paths"
  )
  expect_true(
    all(grepl("\\.png$", basename(physical_files), ignore.case = TRUE)),
    "A truncated output without the .png extension was created"
  )
})

run_case("an overlong directory fails before attempting to save a PNG", {
  impossible_dir <- file.path(test_root, "directory_without_filename_budget")
  dir.create(impossible_dir, recursive = TRUE)
  impossible_dir_abs <- normalizePath(
    impossible_dir, winslash = "/", mustWork = TRUE
  )
  # One readable character, '__', twelve hexadecimal characters and the
  # '_2x2.png' suffix are the shortest valid compact filename.
  minimum_compact_filename <- 1L + 2L + 12L + nchar("_2x2.png")
  impossible_limit <- nchar(impossible_dir_abs, type = "chars") +
    minimum_compact_filename

  expect_error_matching(
    script_env$save_png_dimensions(
      plot_object = plot_object,
      output_dir = impossible_dir,
      file_stem = paste(rep("too_long", 20L), collapse = "_"),
      dimensions = dimensions,
      dpi = 72,
      max_path_length = impossible_limit
    ),
    patterns = c("output_dir", "short"),
    label = "An output directory without filename budget must be rejected"
  )

  expect_true(
    length(list.files(impossible_dir, all.files = FALSE)) == 0L,
    "The impossible-path case wrote a partial or truncated file"
  )
})

if (length(failures) > 0L) {
  stop(
    "P3 portable PNG RED failures:\n- ",
    paste(failures, collapse = "\n- "),
    call. = FALSE
  )
}

message("PASS: portable PNG filenames are deterministic and filesystem-safe")
