# 23_time_rescaling_misspecification.R
#
# Refactored version.
# Time-rescaling under correct Hawkes and wrong Poisson models.

set.seed(23)

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

fit_hawkes <- fit_hawkes_continuous(
  events = events,
  T_end = T_end,
  method = "BFGS"
)

mu_hat_hawkes <- fit_hawkes$mu
alpha_hat <- fit_hawkes$alpha
beta_hat <- fit_hawkes$beta
branching_hat <- fit_hawkes$branching

mu_hat_poisson <- length(events) / T_end

cat("Hawkes branching estimate:", branching_hat, "\n")
cat("Poisson mu estimate:", mu_hat_poisson, "\n")

z_hawkes <- hawkes_rescaled_gaps(
  events = events,
  mu = mu_hat_hawkes,
  alpha = alpha_hat,
  beta = beta_hat
)

z_poisson <- rescaled_gaps_from_compensator(
  poisson_compensator_at_events(events, mu_hat_poisson)
)

summary_hawkes <- ks_exp1_summary(z_hawkes)
summary_poisson <- ks_exp1_summary(z_poisson)

summary_df <- rbind(
  data.frame(model = "Hawkes_correct", summary_hawkes),
  data.frame(model = "Poisson_wrong", summary_poisson)
)

fit_summary <- data.frame(
  n_events = length(events),
  mu_true = mu_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = branching_true,
  mu_hat_hawkes = mu_hat_hawkes,
  alpha_hat = alpha_hat,
  beta_hat = beta_hat,
  branching_hat = branching_hat,
  mu_hat_poisson = mu_hat_poisson
)

print(fit_summary)
print(summary_df)

write.csv(
  fit_summary,
  "results/time_rescaling_misspecification_fit_summary.csv",
  row.names = FALSE
)

write.csv(
  summary_df,
  "results/time_rescaling_misspecification_summary.csv",
  row.names = FALSE
)

residual_df <- rbind(
  data.frame(
    model = "Hawkes_correct",
    z = z_hawkes
  ),
  data.frame(
    model = "Poisson_wrong",
    z = z_poisson
  )
)

write.csv(
  residual_df,
  "results/time_rescaling_misspecification_residuals.csv",
  row.names = FALSE
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
  facet_wrap(~ model) +
  labs(
    title = "Time-rescaled residuals: correct versus wrong model",
    subtitle = "Histogram should match Exp(1) density if model is correct",
    x = "rescaled gap",
    y = "density"
  ) +
  theme_minimal()

ggsave(
  "figures/time_rescaling_misspecification_histogram.pdf",
  p_hist,
  width = 8,
  height = 4
)

ggsave(
  "figures/time_rescaling_misspecification_histogram.png",
  p_hist,
  width = 8,
  height = 4,
  dpi = 300
)

make_qq_df <- function(z, model_name) {
  data.frame(
    model = model_name,
    theoretical_exp = qexp(ppoints(length(z)), rate = 1),
    empirical_sorted = sort(z)
  )
}

qq_df <- rbind(
  make_qq_df(z_hawkes, "Hawkes_correct"),
  make_qq_df(z_poisson, "Poisson_wrong")
)

p_qq <- ggplot(
  qq_df,
  aes(
    x = theoretical_exp,
    y = empirical_sorted
  )
) +
  geom_point(alpha = 0.5) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_wrap(~ model) +
  labs(
    title = "QQ plot of rescaled gaps",
    subtitle = "Correct model should follow the diagonal",
    x = "theoretical Exp(1) quantiles",
    y = "empirical quantiles"
  ) +
  theme_minimal()

ggsave(
  "figures/time_rescaling_misspecification_qqplot.pdf",
  p_qq,
  width = 8,
  height = 4
)

ggsave(
  "figures/time_rescaling_misspecification_qqplot.png",
  p_qq,
  width = 8,
  height = 4,
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
  facet_wrap(~ model) +
  labs(
    title = "Empirical CDF of rescaled gaps",
    subtitle = "Dashed curve: Exp(1) CDF",
    x = "rescaled gap",
    y = "CDF"
  ) +
  theme_minimal()

ggsave(
  "figures/time_rescaling_misspecification_ecdf.pdf",
  p_cdf,
  width = 8,
  height = 4
)

ggsave(
  "figures/time_rescaling_misspecification_ecdf.png",
  p_cdf,
  width = 8,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/time_rescaling_misspecification_fit_summary.csv\n")
cat("- results/time_rescaling_misspecification_summary.csv\n")
cat("- results/time_rescaling_misspecification_residuals.csv\n")
cat("- figures/time_rescaling_misspecification_histogram.pdf/png\n")
cat("- figures/time_rescaling_misspecification_qqplot.pdf/png\n")
cat("- figures/time_rescaling_misspecification_ecdf.pdf/png\n")
