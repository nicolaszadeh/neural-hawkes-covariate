# 30_SPBDBS2025_algorithm2_gof_calibration.R
#
# Refactored version using the shared R/ library.
# Monte Carlo calibration of SPBDBS2025 Algorithm 2.

set.seed(30)

source("R/load_all.R")
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

T_end <- 1000
n_rep <- 50
alpha_level <- 0.05

mu_true <- 0.8
alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

algorithm2_test_once <- function(events, T_end) {
  fit_out <- fit_hawkes(
    events = events,
    T_end = T_end,
    method = "BFGS",
    use_gradient = TRUE
  )

  test <- hawkes_algorithm2_test(
    events = events,
    T_end = T_end,
    mu = fit_out$mu,
    alpha = fit_out$alpha,
    beta = fit_out$beta,
    alpha_level = alpha_level
  )

  cbind(
    test,
    data.frame(
      mu_hat = fit_out$mu,
      alpha_hat = fit_out$alpha,
      beta_hat = fit_out$beta,
      branching_hat = fit_out$branching,
      convergence = fit_out$fit$convergence
    )
  )
}

results_list <- list()

for (rep in seq_len(n_rep)) {
  events <- simulate_hawkes(
    T_end = T_end,
    mu = mu_true,
    alpha = alpha_true,
    beta = beta_true
  )

  out <- algorithm2_test_once(events, T_end)
  out$replicate <- rep

  results_list[[rep]] <- out

  cat(
    "rep", rep, "of", n_rep,
    "- n =", out$n_events,
    "- usual p =", round(out$usual_p, 3),
    "- corrected p =", round(out$corrected_p, 3),
    "- subsampled p =", round(out$subsampled_p, 3),
    "\n"
  )
}

results_df <- do.call(rbind, results_list)

summary_df <- data.frame(
  test = c("usual", "corrected", "subsampled_corrected"),
  rejection_rate = c(
    mean(results_df$usual_reject),
    mean(results_df$corrected_reject),
    mean(results_df$subsampled_reject)
  ),
  target_level = alpha_level,
  mean_p_value = c(
    mean(results_df$usual_p),
    mean(results_df$corrected_p),
    mean(results_df$subsampled_p)
  ),
  median_p_value = c(
    median(results_df$usual_p),
    median(results_df$corrected_p),
    median(results_df$subsampled_p)
  ),
  mean_D = c(
    mean(results_df$usual_D),
    mean(results_df$corrected_D),
    mean(results_df$subsampled_D)
  ),
  mean_n_events = mean(results_df$n_events),
  mean_m_subsample = mean(results_df$m_subsample),
  mean_nonpositive_corrected =
    mean(results_df$n_nonpositive_corrected)
)

print(summary_df)

write.csv(
  results_df,
  "results/SPBDBS2025_algorithm2_gof_calibration_replicates.csv",
  row.names = FALSE
)

write.csv(
  summary_df,
  "results/SPBDBS2025_algorithm2_gof_calibration_summary.csv",
  row.names = FALSE
)

p_rejection <- ggplot(summary_df, aes(x = test, y = rejection_rate)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = alpha_level, linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "SPBDBS2025 Algorithm 2 GOF calibration",
    subtitle = "Dashed line: nominal 5% level",
    x = "test",
    y = "rejection rate"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_algorithm2_gof_calibration_rejection_rates.pdf",
  p_rejection,
  width = 7,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_algorithm2_gof_calibration_rejection_rates.png",
  p_rejection,
  width = 7,
  height = 4,
  dpi = 300
)

p_values_long <- rbind(
  data.frame(test = "usual", p_value = results_df$usual_p),
  data.frame(test = "corrected", p_value = results_df$corrected_p),
  data.frame(test = "subsampled_corrected", p_value = results_df$subsampled_p)
)

p_hist <- ggplot(p_values_long, aes(x = p_value)) +
  geom_histogram(bins = 20, boundary = 0) +
  facet_wrap(~ test) +
  labs(
    title = "SPBDBS2025 Algorithm 2 GOF p-values under H0",
    subtitle = "Under calibration, p-values should be roughly uniform",
    x = "p-value",
    y = "count"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_algorithm2_gof_calibration_pvalues.pdf",
  p_hist,
  width = 8,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_algorithm2_gof_calibration_pvalues.png",
  p_hist,
  width = 8,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/SPBDBS2025_algorithm2_gof_calibration_replicates.csv\n")
cat("- results/SPBDBS2025_algorithm2_gof_calibration_summary.csv\n")
cat("- figures/SPBDBS2025_algorithm2_gof_calibration_rejection_rates.pdf/png\n")
cat("- figures/SPBDBS2025_algorithm2_gof_calibration_pvalues.pdf/png\n")
