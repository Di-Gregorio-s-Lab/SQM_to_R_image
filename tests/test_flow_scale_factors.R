source("check_flow_scale_factors.R")

targets <- data.frame(
  sample = c("S1", "S1", "S1"),
  ko_id = c("K00001", "K00002", "K00003"),
  sqm_tpm = c(12, 0, 0)
)
observed <- data.frame(
  sample = c("S1", "S1", "S1"),
  ko_id = c("K00001", "K00002", "K00003"),
  raw_tpm = c(12, 0, 4)
)

result <- make_factor_table("Pathway", targets, observed)
stopifnot(
  identical(result$status, c("scaled", "both_zero", "target_zero")),
  identical(result$factor[[1L]], 1),
  is.na(result$factor[[2L]]),
  is.na(result$factor[[3L]])
)

bad <- targets[1L, ]
bad_observed <- transform(observed[1L, ], raw_tpm = 0)
error <- tryCatch(make_factor_table("Pathway", bad, bad_observed), error = conditionMessage)
stopifnot(grepl("Cannot allocate positive SQM KO TPM", error, fixed = TRUE))

cat("test_flow_scale_factors: OK\n")
