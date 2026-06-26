# 31_SPBDBS2025_covariate_model_example.R
#
# Refactored version using the shared R/ library.
# Covariate-Hawkes model example with Wald tests.

set.seed(3101)

source("R/load_all.R")
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

T_end <- 1000
dt <- 0.01
grid <- seq(0, T_end, by = dt)

gamma0_true <- -0.5
gamma1_true <- 0.9

alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

kappa_true <- 0.4
sigma_true <- 0.8
X0 <- 0

theta_true <- c(
  gamma0 = gamma0_true,
  gamma1 = gamma1_true,
  log_alpha = log(alpha_true),
  log_beta = log(beta_true)
)

X <- simulate_ou(
  grid = grid,
  kappa = kappa_true,
  sigma = sigma_true,
  X0 = X0
)

events <- simulate_covariate_hawkes(
  grid = grid,
  X = X,
  gamma0 = gamma0_true,
  gamma1 = gamma1_true,
  alpha = alpha_true,
  beta = beta_true
)

cat("Number of events:", length(events), "\n")

fit_out <- fit_covariate_hawkes(
  events = events,
  grid = grid,
  X = X,
  method = "Nelder-Mead"
)

cat("MLE:\n")
cat("gamma0_hat =", fit_out$gamma0, "\n")
cat("gamma1_hat =", fit_out$gamma1, "\n")
cat("alpha_hat =", fit_out$alpha, "\n")
cat("beta_hat =", fit_out$beta, "\n")
cat("branching_hat =", fit_out$branching, "\n")

info <- covariate_hawkes_observed_information(
  events = events,
  grid = grid,
  X = X,
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
  gamma0_true = gamma0_true,
  gamma1_true = gamma1_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = branching_true,
  gamma0_hat = fit_out$gamma0,
  gamma1_hat = fit_out$gamma1,
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
  "results/SPBDBS2025_covariate_model_fit_summary.csv",
  row.names = FALSE
)

write.csv(
  wald_results,
  "results/SPBDBS2025_covariate_model_wald_results.csv",
  row.names = FALSE
)

write.csv(
  information_df,
  "results/SPBDBS2025_covariate_model_information_summary.csv",
  row.names = FALSE
)

write.csv(
  info$H_total,
  "results/SPBDBS2025_covariate_model_observed_information_total.csv",
  row.names = TRUE
)

write.csv(
  info$I_hat,
  "results/SPBDBS2025_covariate_model_observed_information_scaled.csv",
  row.names = TRUE
)

baseline_true <- exp(gamma0_true + gamma1_true * X)
baseline_hat <- exp(fit_out$gamma0 + fit_out$gamma1 * X)

plot_df <- data.frame(
  t = grid,
  X = X,
  baseline_true = baseline_true,
  baseline_hat = baseline_hat
)

p_baseline <- ggplot(plot_df, aes(x = t)) +
  geom_line(aes(y = baseline_true), linewidth = 0.6) +
  geom_line(aes(y = baseline_hat), linetype = "dashed", linewidth = 0.6) +
  labs(
    title = "SPBDBS2025 covariate-Hawkes baseline",
    subtitle = "Solid: true baseline; dashed: fitted baseline",
    x = "time",
    y = "baseline intensity"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_covariate_model_baseline.pdf",
  p_baseline,
  width = 8,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_covariate_model_baseline.png",
  p_baseline,
  width = 8,
  height = 4,
  dpi = 300
)

p_covariate <- ggplot(plot_df, aes(x = t, y = X)) +
  geom_line(linewidth = 0.6) +
  labs(
    title = "OU covariate path",
    x = "time",
    y = "X(t)"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_covariate_model_X_path.pdf",
  p_covariate,
  width = 8,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_covariate_model_X_path.png",
  p_covariate,
  width = 8,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/SPBDBS2025_covariate_model_fit_summary.csv\n")
cat("- results/SPBDBS2025_covariate_model_wald_results.csv\n")
cat("- results/SPBDBS2025_covariate_model_information_summary.csv\n")
cat("- results/SPBDBS2025_covariate_model_observed_information_total.csv\n")
cat("- results/SPBDBS2025_covariate_model_observed_information_scaled.csv\n")
cat("- figures/SPBDBS2025_covariate_model_baseline.pdf/png\n")
cat("- figures/SPBDBS2025_covariate_model_X_path.pdf/png\n")
