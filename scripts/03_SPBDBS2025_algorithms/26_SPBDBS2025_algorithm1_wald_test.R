# 26_SPBDBS2025_algorithm1_wald_test.R
#
# Refactored version using the shared R/ library.
# SPBDBS2025 Algorithm 1: one-parameter Wald tests.

set.seed(26)

source("R/load_all.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

T_end <- 1000

mu_true <- 0.8
alpha_true <- 0.35
beta_true <- 1.2

branching_true <- alpha_true / beta_true

theta_true <- log(c(mu_true, alpha_true, beta_true))
names(theta_true) <- c("log_mu", "log_alpha", "log_beta")

events <- simulate_hawkes(
  T_end = T_end,
  mu = mu_true,
  alpha = alpha_true,
  beta = beta_true
)

cat("Number of events:", length(events), "\n")

fit_out <- fit_hawkes(
  events = events,
  T_end = T_end,
  method = "BFGS",
  use_gradient = TRUE
)

info <- hawkes_observed_information(
  events = events,
  T_end = T_end,
  theta_hat = fit_out$theta_hat
)

wald_results <- wald_all_true_parameters(
  theta_hat = fit_out$theta_hat,
  I_inv_hat = info$I_inv_hat,
  theta_true = theta_true,
  T_end = T_end,
  alpha_level = 0.05
)

fit_summary <- data.frame(
  n_events = length(events),
  T_end = T_end,
  mu_true = mu_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = branching_true,
  mu_hat = fit_out$mu,
  alpha_hat = fit_out$alpha,
  beta_hat = fit_out$beta,
  branching_hat = fit_out$branching,
  neg_loglik = fit_out$fit$value,
  convergence = fit_out$fit$convergence
)

information_df <- information_summary(
  theta_hat = fit_out$theta_hat,
  H_total = info$H_total,
  I_inv_hat = info$I_inv_hat
)

print(fit_summary)
print(wald_results)
print(information_df)

write.csv(
  fit_summary,
  "results/SPBDBS2025_algorithm1_fit_summary.csv",
  row.names = FALSE
)

write.csv(
  wald_results,
  "results/SPBDBS2025_algorithm1_wald_results.csv",
  row.names = FALSE
)

write.csv(
  information_df,
  "results/SPBDBS2025_algorithm1_information_summary.csv",
  row.names = FALSE
)

write.csv(
  info$H_total,
  "results/SPBDBS2025_algorithm1_observed_information_total.csv",
  row.names = TRUE
)

write.csv(
  info$I_hat,
  "results/SPBDBS2025_algorithm1_observed_information_scaled.csv",
  row.names = TRUE
)

cat("\nSaved:\n")
cat("- results/SPBDBS2025_algorithm1_fit_summary.csv\n")
cat("- results/SPBDBS2025_algorithm1_wald_results.csv\n")
cat("- results/SPBDBS2025_algorithm1_information_summary.csv\n")
cat("- results/SPBDBS2025_algorithm1_observed_information_total.csv\n")
cat("- results/SPBDBS2025_algorithm1_observed_information_scaled.csv\n")
