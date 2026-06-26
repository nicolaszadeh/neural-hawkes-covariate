#!/usr/bin/env Rscript

# ============================================================
# Collect SPBDBS2025 results into one project summary
# ============================================================

setwd(normalizePath("~/GitHub/biomedical-hawkes-covariates"))

library(readr)
library(dplyr)

safe_read <- function(file) {
  if (file.exists(file)) read_csv(file, show_col_types = FALSE) else NULL
}

out <- tibble(
  script = character(),
  message = character(),
  metric = character(),
  value = character()
)

add_row_safe <- function(script, message, metric, value) {
  tibble(
    script = script,
    message = message,
    metric = metric,
    value = as.character(value)
  )
}

# 26
x <- safe_read("results/SPBDBS2025_algorithm1_fit_summary.csv")
if (!is.null(x)) {
  out <- bind_rows(out, add_row_safe(
    "26", "Algorithm 1 Wald test",
    "branching_hat",
    round(x$branching_hat[1], 3)
  ))
}

# 27
x <- safe_read("results/SPBDBS2025_algorithm2_test_summary.csv")
if (!is.null(x)) {
  out <- bind_rows(out, add_row_safe(
    "27", "Algorithm 2 corrected GOF",
    "corrected_p",
    round(x$corrected_p[1], 4)
  ))
}

# 28
x <- safe_read("results/SPBDBS2025_algorithm3_equality_wald_results.csv")
if (!is.null(x)) {
  out <- bind_rows(out, add_row_safe(
    "28", "Algorithm 3 equality Wald test",
    "min_p_value",
    round(min(x$p_value, na.rm = TRUE), 4)
  ))
}

# 29 algorithm 1
x <- safe_read("results/SPBDBS2025_wald_calibration_algorithm1_summary.csv")
if (!is.null(x)) {
  out <- bind_rows(out, add_row_safe(
    "29a", "Wald calibration Algorithm 1",
    "mean_rejection_rate",
    round(mean(x$rejection_rate, na.rm = TRUE), 3)
  ))
}

# 29 algorithm 3
x <- safe_read("results/SPBDBS2025_wald_calibration_algorithm3_summary.csv")
if (!is.null(x)) {
  out <- bind_rows(out, add_row_safe(
    "29b", "Wald calibration Algorithm 3",
    "rejection_rate",
    round(x$rejection_rate[1], 3)
  ))
}

# 30
x <- safe_read("results/SPBDBS2025_algorithm2_gof_calibration_summary.csv")
if (!is.null(x)) {
  out <- bind_rows(out, add_row_safe(
    "30", "Algorithm 2 GOF calibration",
    "subsampled_rejection_rate",
    round(x$rejection_rate[x$test == "subsampled_corrected"][1], 3)
  ))
}

# 31
x <- safe_read("results/SPBDBS2025_covariate_model_fit_summary.csv")
if (!is.null(x)) {
  out <- bind_rows(out, add_row_safe(
    "31", "Covariate-Hawkes model example",
    "branching_hat",
    round(x$branching_hat[1], 3)
  ))
}

# 32
x <- safe_read("results/SPBDBS2025_covariate_wald_calibration_summary.csv")
if (!is.null(x)) {
  out <- bind_rows(out, add_row_safe(
    "32", "Covariate Wald calibration",
    "mean_rejection_rate",
    round(mean(x$rejection_rate, na.rm = TRUE), 3)
  ))
}

# 33
x <- safe_read("results/SPBDBS2025_covariate_algorithm2_gof_summary.csv")
if (!is.null(x)) {
  out <- bind_rows(
    out,
    add_row_safe(
      "33", "Covariate Algorithm 2 GOF",
      "usual_p",
      round(x$usual_p[1], 4)
    ),
    add_row_safe(
      "33", "Covariate Algorithm 2 GOF",
      "corrected_p",
      round(x$corrected_p[1], 4)
    ),
    add_row_safe(
      "33", "Covariate Algorithm 2 GOF",
      "subsampled_p",
      round(x$subsampled_p[1], 4)
    ),
    add_row_safe(
      "33", "Covariate Algorithm 2 GOF",
      "mean_rescaled_gap",
      round(x$mean_z[1], 3)
    ),
    add_row_safe(
      "33", "Covariate Algorithm 2 GOF",
      "var_rescaled_gap",
      round(x$var_z[1], 3)
    )
  )
}

print(out)

write_csv(out, "results/SPBDBS2025_project_summary.csv")

cat("\nSaved:\n")
cat("- results/SPBDBS2025_project_summary.csv\n")
