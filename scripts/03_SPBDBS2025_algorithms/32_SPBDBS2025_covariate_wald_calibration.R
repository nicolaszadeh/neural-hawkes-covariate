# 32_SPBDBS2025_covariate_wald_calibration.R
#
# Refactored version using the shared R/ library.
# Monte Carlo calibration of covariate-Hawkes Wald tests.

set.seed(32)

source("R/load_all.R")
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

T_end <- 500
dt <- 0.01
grid <- seq(0, T_end, by = dt)

n_rep <- 30
alpha_level <- 0.05

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

results_list <- list()
counter <- 1

for (rep in seq_len(n_rep)) {
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

  fit_out <- fit_covariate_hawkes(
    events = events,
    grid = grid,
    X = X,
    method = "Nelder-Mead"
  )

  info <- covariate_hawkes_observed_information(
    events = events,
    grid = grid,
    X = X,
    theta_hat = fit_out$theta_hat
  )

  tests <- wald_all_true_parameters(
    theta_hat = fit_out$theta_hat,
    I_inv_hat = info$I_inv_hat,
    theta_true = theta_true,
    T_end = T_end,
    alpha_level = alpha_level
  )

  tests$reject <- tests$reject_5_percent
  tests$replicate <- rep
  tests$n_events <- length(events)
  tests$convergence <- fit_out$fit$convergence
  tests$gamma0_hat <- fit_out$gamma0
  tests$gamma1_hat <- fit_out$gamma1
  tests$alpha_hat <- fit_out$alpha
  tests$beta_hat <- fit_out$beta
  tests$branching_hat <- fit_out$branching

  for (i in seq_len(nrow(tests))) {
    results_list[[counter]] <- tests[i, ]
    counter <- counter + 1
  }

  cat(
    "rep", rep, "of", n_rep,
    "- n =", length(events),
    "- gamma0 =", round(fit_out$gamma0, 3),
    "- gamma1 =", round(fit_out$gamma1, 3),
    "- branching =", round(fit_out$branching, 3),
    "\n"
  )
}

results_df <- do.call(rbind, results_list)

summary_reject <- aggregate(reject ~ parameter, data = results_df, FUN = mean)
names(summary_reject)[2] <- "rejection_rate"

summary_mean_p <- aggregate(p_value ~ parameter, data = results_df, FUN = mean)
names(summary_mean_p)[2] <- "mean_p_value"

summary_median_p <- aggregate(p_value ~ parameter, data = results_df, FUN = median)
names(summary_median_p)[2] <- "median_p_value"

summary_mean_est <- aggregate(theta_hat ~ parameter, data = results_df, FUN = mean)
names(summary_mean_est)[2] <- "mean_theta_hat"

summary_sd_est <- aggregate(theta_hat ~ parameter, data = results_df, FUN = sd)
names(summary_sd_est)[2] <- "sd_theta_hat"

summary_df <- merge(summary_reject, summary_mean_p, by = "parameter")
summary_df <- merge(summary_df, summary_median_p, by = "parameter")
summary_df <- merge(summary_df, summary_mean_est, by = "parameter")
summary_df <- merge(summary_df, summary_sd_est, by = "parameter")
summary_df$target_level <- alpha_level

truth_df <- data.frame(
  parameter = names(theta_true),
  theta_true = as.numeric(theta_true)
)

summary_df <- merge(summary_df, truth_df, by = "parameter")
summary_df$bias_theta <- summary_df$mean_theta_hat - summary_df$theta_true

summary_df$parameter <- factor(summary_df$parameter, levels = names(theta_true))
summary_df <- summary_df[order(summary_df$parameter), ]
summary_df$parameter <- as.character(summary_df$parameter)
rownames(summary_df) <- NULL

overall_summary <- data.frame(
  n_rep = n_rep,
  T_end = T_end,
  mean_n_events = mean(results_df$n_events[results_df$parameter == "gamma0"]),
  gamma0_true = gamma0_true,
  gamma1_true = gamma1_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = branching_true,
  kappa_true = kappa_true,
  sigma_true = sigma_true
)

print(overall_summary)
print(summary_df)

write.csv(
  results_df,
  "results/SPBDBS2025_covariate_wald_calibration_replicates.csv",
  row.names = FALSE
)

write.csv(
  summary_df,
  "results/SPBDBS2025_covariate_wald_calibration_summary.csv",
  row.names = FALSE
)

write.csv(
  overall_summary,
  "results/SPBDBS2025_covariate_wald_calibration_overall.csv",
  row.names = FALSE
)

p_reject <- ggplot(summary_df, aes(x = parameter, y = rejection_rate)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = alpha_level, linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "SPBDBS2025 covariate Wald calibration",
    subtitle = "Dashed line: nominal 5% level",
    x = "parameter",
    y = "rejection rate"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_covariate_wald_calibration_rejection_rates.pdf",
  p_reject,
  width = 7,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_covariate_wald_calibration_rejection_rates.png",
  p_reject,
  width = 7,
  height = 4,
  dpi = 300
)

p_values <- ggplot(results_df, aes(x = p_value)) +
  geom_histogram(bins = 20, boundary = 0) +
  facet_wrap(~ parameter) +
  labs(
    title = "SPBDBS2025 covariate Wald p-values under H0",
    subtitle = "Under correct calibration, p-values should be roughly uniform",
    x = "p-value",
    y = "count"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_covariate_wald_calibration_pvalues.pdf",
  p_values,
  width = 8,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_covariate_wald_calibration_pvalues.png",
  p_values,
  width = 8,
  height = 4,
  dpi = 300
)

p_estimates <- ggplot(results_df, aes(x = parameter, y = theta_hat)) +
  geom_boxplot() +
  geom_point(
    data = truth_df,
    aes(x = parameter, y = theta_true),
    inherit.aes = FALSE,
    shape = 4,
    size = 3
  ) +
  labs(
    title = "SPBDBS2025 covariate Wald parameter estimates",
    subtitle = "Crosses indicate true parameter values",
    x = "parameter",
    y = "estimate"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_covariate_wald_calibration_estimates.pdf",
  p_estimates,
  width = 8,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_covariate_wald_calibration_estimates.png",
  p_estimates,
  width = 8,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/SPBDBS2025_covariate_wald_calibration_replicates.csv\n")
cat("- results/SPBDBS2025_covariate_wald_calibration_summary.csv\n")
cat("- results/SPBDBS2025_covariate_wald_calibration_overall.csv\n")
cat("- figures/SPBDBS2025_covariate_wald_calibration_rejection_rates.pdf/png\n")
cat("- figures/SPBDBS2025_covariate_wald_calibration_pvalues.pdf/png\n")
cat("- figures/SPBDBS2025_covariate_wald_calibration_estimates.pdf/png\n")
