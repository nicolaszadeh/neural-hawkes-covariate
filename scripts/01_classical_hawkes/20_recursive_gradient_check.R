# 20_recursive_gradient_check.R
#
# Refactored version.
# Compares recursive analytical score with finite differences.

set.seed(20)

source("R/load_all.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

T_end <- 500

mu_true <- 0.8
alpha_true <- 0.35
beta_true <- 1.2

theta_true <- log(c(mu_true, alpha_true, beta_true))
names(theta_true) <- c("log_mu", "log_alpha", "log_beta")

events <- simulate_hawkes(
  T_end = T_end,
  mu = mu_true,
  alpha = alpha_true,
  beta = beta_true
)

cat("Number of events:", length(events), "\n")

f <- function(theta) {
  hawkes_loglik_theta(
    theta = theta,
    events = events,
    T_end = T_end
  )
}

score_analytic <- hawkes_score_theta(
  theta = theta_true,
  events = events,
  T_end = T_end
)

score_numeric <- finite_difference_gradient(
  f = f,
  theta = theta_true,
  eps = 1e-5
)

names(score_numeric) <- names(score_analytic)

diff <- score_analytic - score_numeric

comparison <- data.frame(
  parameter = names(score_analytic),
  analytic_score = as.numeric(score_analytic),
  numeric_score = as.numeric(score_numeric),
  difference = as.numeric(diff),
  abs_difference = abs(as.numeric(diff))
)

summary_df <- data.frame(
  T_end = T_end,
  n_events = length(events),
  mu_true = mu_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  max_abs_difference = max(comparison$abs_difference)
)

print(comparison)
print(summary_df)

write.csv(
  comparison,
  "results/recursive_gradient_check_comparison.csv",
  row.names = FALSE
)

write.csv(
  summary_df,
  "results/recursive_gradient_check_summary.csv",
  row.names = FALSE
)

cat("\nSaved:\n")
cat("- results/recursive_gradient_check_comparison.csv\n")
cat("- results/recursive_gradient_check_summary.csv\n")
