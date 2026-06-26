# 29_SPBDBS2025_wald_calibration.R
#
# Refactored version using the shared R/ library.
# Monte Carlo calibration of SPBDBS2025 Algorithms 1 and 3.

set.seed(29)

source("R/load_all.R")
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

T_end <- 1000
n_rep <- 50
alpha_level <- 0.05

mu_true_A1 <- 0.8
alpha_true_A1 <- 0.35
beta_true_A1 <- 1.2

mu_true_A3 <- 1.2
alpha_true_A3 <- 0.35
beta_true_A3 <- 1.2

run_wald_fit <- function(events, T_end) {
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

  list(fit_out = fit_out, info = info)
}

theta_true_A1 <- log(c(mu_true_A1, alpha_true_A1, beta_true_A1))
names(theta_true_A1) <- c("log_mu", "log_alpha", "log_beta")

results_A1 <- list()
counter <- 1

for (rep in seq_len(n_rep)) {
  events <- simulate_hawkes(
    T_end = T_end,
    mu = mu_true_A1,
    alpha = alpha_true_A1,
    beta = beta_true_A1
  )

  out <- run_wald_fit(events, T_end)

  tests <- wald_all_true_parameters(
    theta_hat = out$fit_out$theta_hat,
    I_inv_hat = out$info$I_inv_hat,
    theta_true = theta_true_A1,
    T_end = T_end,
    alpha_level = alpha_level
  )

  tests$replicate <- rep
  tests$n_events <- length(events)
  tests$convergence <- out$fit_out$fit$convergence

  for (i in seq_len(nrow(tests))) {
    results_A1[[counter]] <- tests[i, ]
    counter <- counter + 1
  }

  cat(
    "Algorithm 1 rep", rep, "of", n_rep,
    "- n =", length(events),
    "\n"
  )
}

A1_df <- do.call(rbind, results_A1)
A1_df$reject <- A1_df$reject_5_percent

theta_true_A3 <- log(c(mu_true_A3, alpha_true_A3, beta_true_A3))
names(theta_true_A3) <- c("log_mu", "log_alpha", "log_beta")

results_A3 <- list()

for (rep in seq_len(n_rep)) {
  events <- simulate_hawkes(
    T_end = T_end,
    mu = mu_true_A3,
    alpha = alpha_true_A3,
    beta = beta_true_A3
  )

  out <- run_wald_fit(events, T_end)

  test <- wald_equality(
    theta_hat = out$fit_out$theta_hat,
    I_inv_hat = out$info$I_inv_hat,
    parameter_i = "log_mu",
    parameter_j = "log_beta",
    T_end = T_end,
    alpha_level = alpha_level
  )

  test$replicate <- rep
  test$n_events <- length(events)
  test$convergence <- out$fit_out$fit$convergence
  test$reject <- test$reject_5_percent

  results_A3[[rep]] <- test

  cat(
    "Algorithm 3 rep", rep, "of", n_rep,
    "- n =", length(events),
    "- p =", round(test$p_value, 3),
    "- reject =", test$reject,
    "\n"
  )
}

A3_df <- do.call(rbind, results_A3)

A1_summary <- aggregate(
  reject ~ parameter,
  data = A1_df,
  FUN = mean
)
names(A1_summary)[2] <- "rejection_rate"
A1_summary$target_level <- alpha_level

A1_mean_p <- aggregate(p_value ~ parameter, data = A1_df, FUN = mean)
names(A1_mean_p)[2] <- "mean_p_value"

A1_median_p <- aggregate(p_value ~ parameter, data = A1_df, FUN = median)
names(A1_median_p)[2] <- "median_p_value"

A1_summary <- merge(A1_summary, A1_mean_p, by = "parameter")
A1_summary <- merge(A1_summary, A1_median_p, by = "parameter")

A3_summary <- data.frame(
  test = "log_mu_equals_log_beta",
  rejection_rate = mean(A3_df$reject),
  target_level = alpha_level,
  mean_p_value = mean(A3_df$p_value),
  median_p_value = median(A3_df$p_value),
  mean_n_events = mean(A3_df$n_events)
)

print(A1_summary)
print(A3_summary)

write.csv(
  A1_df,
  "results/SPBDBS2025_wald_calibration_algorithm1_replicates.csv",
  row.names = FALSE
)

write.csv(
  A1_summary,
  "results/SPBDBS2025_wald_calibration_algorithm1_summary.csv",
  row.names = FALSE
)

write.csv(
  A3_df,
  "results/SPBDBS2025_wald_calibration_algorithm3_replicates.csv",
  row.names = FALSE
)

write.csv(
  A3_summary,
  "results/SPBDBS2025_wald_calibration_algorithm3_summary.csv",
  row.names = FALSE
)

p_A1 <- ggplot(A1_summary, aes(x = parameter, y = rejection_rate)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = alpha_level, linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "SPBDBS2025 Algorithm 1 calibration",
    subtitle = "Dashed line: nominal 5% level",
    x = "parameter",
    y = "rejection rate"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_wald_calibration_algorithm1.pdf",
  p_A1,
  width = 7,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_wald_calibration_algorithm1.png",
  p_A1,
  width = 7,
  height = 4,
  dpi = 300
)

p_A3 <- ggplot(A3_summary, aes(x = test, y = rejection_rate)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = alpha_level, linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "SPBDBS2025 Algorithm 3 calibration",
    subtitle = "Null hypothesis: log_mu = log_beta",
    x = "test",
    y = "rejection rate"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_wald_calibration_algorithm3.pdf",
  p_A3,
  width = 7,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_wald_calibration_algorithm3.png",
  p_A3,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/SPBDBS2025_wald_calibration_algorithm1_replicates.csv\n")
cat("- results/SPBDBS2025_wald_calibration_algorithm1_summary.csv\n")
cat("- results/SPBDBS2025_wald_calibration_algorithm3_replicates.csv\n")
cat("- results/SPBDBS2025_wald_calibration_algorithm3_summary.csv\n")
cat("- figures/SPBDBS2025_wald_calibration_algorithm1.pdf/png\n")
cat("- figures/SPBDBS2025_wald_calibration_algorithm3.pdf/png\n")
