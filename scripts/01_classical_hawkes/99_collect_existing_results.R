# 18_collect_existing_results.R
#
# Collector only.
# Does not rerun simulations.
# Reads existing CSV result files and writes one project summary.

dir.create("results", showWarnings = FALSE)

safe_read <- function(path) {
  if (!file.exists(path)) {
    warning(paste("Missing:", path))
    return(NULL)
  }
  read.csv(path)
}

add_row <- function(rows, script, message, metric, value) {
  rows[[length(rows) + 1]] <- data.frame(
    script = script,
    message = message,
    metric = metric,
    value = as.character(value)
  )
  rows
}

rows <- list()

# 06
x <- safe_read("results/covariate_baseline_model_comparison.csv")
if (!is.null(x)) {
  best <- x$model[which.min(x$AIC)]
  rows <- add_row(rows, "06", "Covariate baseline control",
                  "best_AIC_model", best)
}

# 07
x <- safe_read("results/covariate_hawkes_positive_estimates.csv")
if (!is.null(x)) {
  br <- x$value[x$parameter == "branching_hat"]
  rows <- add_row(rows, "07", "Covariate-Hawkes positive control",
                  "branching_hat", round(br, 3))
}

# 08
x <- safe_read("results/covariate_hawkes_bootstrap_summary.csv")
if (!is.null(x)) {
  rows <- add_row(rows, "08", "Bootstrap LR test",
                  "p_value", x$p_value[1])
}

# 09
x <- safe_read("results/omitted_covariate_bias_summary.csv")
if (!is.null(x)) {
  br <- x$mean_branching[x$gamma1_true == max(x$gamma1_true)]
  rows <- add_row(rows, "09", "Omitted deterministic covariate",
                  "max_fake_branching", round(br, 3))
}

# 10
x <- safe_read("results/covariate_hawkes_recovery_summary.csv")
if (!is.null(x)) {
  br <- x$mean_estimate[x$parameter == "branching"]
  rows <- add_row(rows, "10", "Monte Carlo recovery",
                  "mean_branching_hat", round(br, 3))
}

# 11
x <- safe_read("results/sample_size_effect_summary.csv")
if (!is.null(x)) {
  br <- x$sd_branching[x$T_end == max(x$T_end)]
  rows <- add_row(rows, "11", "Sample-size effect",
                  "sd_branching_at_max_T", round(br, 3))
}

# 12
x <- safe_read("results/power_vs_branching_summary.csv")
if (!is.null(x)) {
  pow <- x$rejection_rate[
    x$branching_true == max(x$branching_true)
  ]
  rows <- add_row(rows, "12", "Power versus branching",
                  "power_at_max_branching", round(pow, 3))
}

# 13
x <- safe_read("results/false_positive_calibration_summary.csv")
if (!is.null(x)) {
  fpr <- mean(x$false_positive_rate)
  rows <- add_row(rows, "13", "False-positive calibration",
                  "mean_false_positive_rate", round(fpr, 3))
}

# 14
x <- safe_read("results/ou_covariate_hawkes_estimates.csv")
if (!is.null(x)) {
  br <- x$value[x$parameter == "branching_hat"]
  rows <- add_row(rows, "14", "OU covariate-Hawkes positive control",
                  "branching_hat", round(br, 3))
}

# 15
x <- safe_read("results/ou_omitted_covariate_bias_summary.csv")
if (!is.null(x)) {
  br <- x$mean_branching[x$gamma1_true == max(x$gamma1_true)]
  rows <- add_row(rows, "15", "Omitted OU covariate",
                  "max_fake_branching", round(br, 3))
}

# 16
x <- safe_read("results/ou_covariate_hawkes_bootstrap_summary.csv")
if (!is.null(x)) {
  rows <- add_row(rows, "16", "OU bootstrap LR test",
                  "corrected_p_value", round(x$p_corrected[1], 4))
}

# 17
x <- safe_read("results/ou_false_positive_calibration_summary.csv")
if (!is.null(x)) {
  fpr <- mean(x$false_positive_rate)
  rows <- add_row(rows, "17", "OU false-positive calibration",
                  "mean_false_positive_rate", round(fpr, 3))
}

project_summary <- do.call(rbind, rows)

write.csv(
  project_summary,
  "results/project_summary.csv",
  row.names = FALSE
)

print(project_summary)

cat("\nSaved:\n")
cat("- results/project_summary.csv\n")