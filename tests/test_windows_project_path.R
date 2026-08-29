source_without_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  main_call <- which(trimws(lines) == "main()")
  stopifnot(length(main_call) == 1L)

  env <- new.env(parent = globalenv())
  eval(parse(text = paste(lines[-main_call], collapse = "\n")), envir = env)
  env
}

script_env <- source_without_main("sqm_plots.R")
project_dir <- "C:\\Users\\unico\\OneDrive - University of Pisa\\Documenti\\UNIPI\\Grani\\Progetti\\TCE\\Shotgun\\minion\\montescudaio\\CS8_All\\CS8_All"

normalized_dir <- script_env$normalize_sqm_project_dir(project_dir)
stopifnot(identical(normalized_dir, gsub("\\\\", "/", project_dir)))

message("PASS: Windows SQM project paths are normalized with forward slashes")
