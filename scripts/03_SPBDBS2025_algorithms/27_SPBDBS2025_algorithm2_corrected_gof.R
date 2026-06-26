# 27_SPBDBS2025_algorithm2_corrected_gof.R
#
# Refactored version using the shared R/ library.
# SPBDBS2025 Algorithm 2: corrected goodness-of-fit increments.

set.seed(27)

source("R/load_all.R")
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

T_end <- 1000

mu_true <- 0.8
alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

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

cat("MLE:\n")
cat("mu_hat =", fit_out$mu, "\n")
cat("alpha_hat =", fit_out$alpha, "\n")
cat("beta_hat =", fit_out$beta, "\n")
cat("branching_hat =", fit_out$branching, "\n")

inc <- hawkes_algorithm2_corrected_increments(
  events = events,
  T_end = T_end,
  mu = fit_out$mu,
  alpha = fit_out$alpha,
  beta = fit_out$beta
)

cat("Non-positive corrected increments:", inc$n_nonpositive, "\n")

test_summary <- hawkes_algorithm2_test(
  events = events,
  T_end = T_end,
  mu = fit_out$mu,
  alpha = fit_out$alpha,
  beta = fit_out$beta,
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

rho_df <- data.frame(
  parameter = names(inc$rho_hat),
  rho_hat = as.numeric(inc$rho_hat)
)

increments_df <- data.frame(
  index = seq_along(inc$usual),
  event_time = events,
  delta_t = diff(c(0, events)),
  usual_increment = inc$usual,
  corrected_increment = inc$corrected,
  corrected_positive = inc$corrected > 0
)

print(fit_summary)
print(rho_df)
print(test_summary)

write.csv(
  fit_summary,
  "results/SPBDBS2025_algorithm2_fit_summary.csv",
  row.names = FALSE
)

write.csv(
  rho_df,
  "results/SPBDBS2025_algorithm2_rho_summary.csv",
  row.names = FALSE
)

write.csv(
  test_summary,
  "results/SPBDBS2025_algorithm2_test_summary.csv",
  row.names = FALSE
)

write.csv(
  inc$I_hat,
  "results/SPBDBS2025_algorithm2_I_hat.csv",
  row.names = TRUE
)

write.csv(
  inc$I_inv_sqrt,
  "results/SPBDBS2025_algorithm2_I_inv_sqrt.csv",
  row.names = TRUE
)

write.csv(
  increments_df,
  "results/SPBDBS2025_algorithm2_increments.csv",
  row.names = FALSE
)

N_events <- length(events)
m_subsample <- floor(N_events^(2 / 3))
subsample_indices <- sample(
  seq_along(inc$positive_corrected),
  size = min(m_subsample, length(inc$positive_corrected)),
  replace = FALSE
)
subsample_corrected <- inc$positive_corrected[subsample_indices]

plot_df <- rbind(
  data.frame(type = "usual", increment = inc$usual),
  data.frame(type = "corrected_positive", increment = inc$positive_corrected),
  data.frame(type = "subsampled_corrected", increment = subsample_corrected)
)

p_hist <- ggplot(plot_df, aes(x = increment)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, boundary = 0) +
  stat_function(fun = dexp, args = list(rate = 1), linewidth = 0.8) +
  facet_wrap(~ type) +
  labs(
    title = "SPBDBS2025 Algorithm 2 corrected increments",
    subtitle = "Histograms compared with Exp(1) density",
    x = "increment",
    y = "density"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_algorithm2_increment_histograms.pdf",
  p_hist,
  width = 9,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_algorithm2_increment_histograms.png",
  p_hist,
  width = 9,
  height = 4,
  dpi = 300
)

make_qq_df <- function(z, label) {
  data.frame(
    type = label,
    theoretical_exp = qexp(ppoints(length(z)), rate = 1),
    empirical_sorted = sort(z)
  )
}

qq_df <- rbind(
  make_qq_df(inc$usual, "usual"),
  make_qq_df(inc$positive_corrected, "corrected_positive"),
  make_qq_df(subsample_corrected, "subsampled_corrected")
)

p_qq <- ggplot(qq_df, aes(x = theoretical_exp, y = empirical_sorted)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ type) +
  labs(
    title = "SPBDBS2025 Algorithm 2 QQ plots",
    subtitle = "Corrected increments should be close to Exp(1)",
    x = "theoretical Exp(1) quantiles",
    y = "empirical quantiles"
  ) +
  theme_minimal()

ggsave(
  "figures/SPBDBS2025_algorithm2_qqplots.pdf",
  p_qq,
  width = 9,
  height = 4
)

ggsave(
  "figures/SPBDBS2025_algorithm2_qqplots.png",
  p_qq,
  width = 9,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/SPBDBS2025_algorithm2_fit_summary.csv\n")
cat("- results/SPBDBS2025_algorithm2_rho_summary.csv\n")
cat("- results/SPBDBS2025_algorithm2_test_summary.csv\n")
cat("- results/SPBDBS2025_algorithm2_I_hat.csv\n")
cat("- results/SPBDBS2025_algorithm2_I_inv_sqrt.csv\n")
cat("- results/SPBDBS2025_algorithm2_increments.csv\n")
cat("- figures/SPBDBS2025_algorithm2_increment_histograms.pdf/png\n")
cat("- figures/SPBDBS2025_algorithm2_qqplots.pdf/png\n")
