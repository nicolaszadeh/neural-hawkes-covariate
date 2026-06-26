# 22_time_rescaling_diagnostics.R
#
# Refactored version.
# Time-rescaling diagnostics for a fitted Hawkes process.

set.seed(22)

source("R/load_all.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

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

fit_out <- fit_hawkes_continuous(
  events = events,
  T_end = T_end,
  method = "BFGS"
)

fit <- fit_out$fit
theta_hat <- fit_out$theta_hat

mu_hat <- fit_out$mu
alpha_hat <- fit_out$alpha
beta_hat <- fit_out$beta
branching_hat <- fit_out$branching

cat("mu_hat:", mu_hat, "\n")
cat("alpha_hat:", alpha_hat, "\n")
cat("beta_hat:", beta_hat, "\n")
cat("branching_hat:", branching_hat, "\n")

z <- hawkes_rescaled_gaps(
  events = events,
  mu = mu_hat,
  alpha = alpha_hat,
  beta = beta_hat
)

residual_summary <- ks_exp1_summary(z)

fit_summary <- data.frame(
  mu_true = mu_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = branching_true,
  mu_hat = mu_hat,
  alpha_hat = alpha_hat,
  beta_hat = beta_hat,
  branching_hat = branching_hat,
  neg_loglik = fit$value,
  convergence = fit$convergence
)

print(fit_summary)
print(residual_summary)

Lambda_events <- hawkes_compensator_at_events(
  events = events,
  mu = mu_hat,
  alpha = alpha_hat,
  beta = beta_hat
)

write.csv(
  fit_summary,
  "results/time_rescaling_fit_summary.csv",
  row.names = FALSE
)

write.csv(
  residual_summary,
  "results/time_rescaling_residual_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    event_index = seq_along(events),
    event_time = events,
    compensator = Lambda_events,
    rescaled_gap = z
  ),
  "results/time_rescaling_residuals.csv",
  row.names = FALSE
)

residual_df <- data.frame(
  z = z,
  theoretical_exp = qexp(
    ppoints(length(z)),
    rate = 1
  ),
  empirical_sorted = sort(z)
)

p_hist <- ggplot(residual_df, aes(x = z)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    boundary = 0
  ) +
  stat_function(
    fun = dexp,
    args = list(rate = 1),
    linewidth = 0.8
  ) +
  labs(
    title = "Time-rescaled residual gaps",
    subtitle = "Histogram should match Exp(1) density",
    x = "rescaled gap",
    y = "density"
  ) +
  theme_minimal()

ggsave(
  "figures/time_rescaling_histogram.pdf",
  p_hist,
  width = 7,
  height = 4
)

ggsave(
  "figures/time_rescaling_histogram.png",
  p_hist,
  width = 7,
  height = 4,
  dpi = 300
)

p_qq <- ggplot(
  residual_df,
  aes(
    x = theoretical_exp,
    y = empirical_sorted
  )
) +
  geom_point(alpha = 0.6) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(
    title = "QQ plot of rescaled gaps",
    subtitle = "Points should follow the diagonal under Exp(1)",
    x = "theoretical Exp(1) quantiles",
    y = "empirical quantiles"
  ) +
  theme_minimal()

ggsave(
  "figures/time_rescaling_qqplot.pdf",
  p_qq,
  width = 5,
  height = 5
)

ggsave(
  "figures/time_rescaling_qqplot.png",
  p_qq,
  width = 5,
  height = 5,
  dpi = 300
)

p_cdf <- ggplot(residual_df, aes(x = z)) +
  stat_ecdf(linewidth = 0.8) +
  stat_function(
    fun = pexp,
    args = list(rate = 1),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(
    title = "Empirical CDF of rescaled gaps",
    subtitle = "Dashed curve: Exp(1) CDF",
    x = "rescaled gap",
    y = "CDF"
  ) +
  theme_minimal()

ggsave(
  "figures/time_rescaling_ecdf.pdf",
  p_cdf,
  width = 7,
  height = 4
)

ggsave(
  "figures/time_rescaling_ecdf.png",
  p_cdf,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/time_rescaling_fit_summary.csv\n")
cat("- results/time_rescaling_residual_summary.csv\n")
cat("- results/time_rescaling_residuals.csv\n")
cat("- figures/time_rescaling_histogram.pdf/png\n")
cat("- figures/time_rescaling_qqplot.pdf/png\n")
cat("- figures/time_rescaling_ecdf.pdf/png\n")
