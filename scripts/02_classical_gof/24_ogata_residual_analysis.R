# 24_ogata_residual_analysis.R
#
# Refactored version.
# Ogata-style residual analysis for correct Hawkes and wrong Poisson models.

set.seed(24)

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

z_hawkes <- hawkes_rescaled_gaps(
  events = events,
  mu = mu_hat_hawkes,
  alpha = alpha_hat,
  beta = beta_hat
)

z_poisson <- rescaled_gaps_from_compensator(
  poisson_compensator_at_events(events, mu_hat_poisson)
)

residual_df <- rbind(
  data.frame(
    model = "Hawkes_correct",
    index = seq_along(z_hawkes),
    z = z_hawkes
  ),
  data.frame(
    model = "Poisson_wrong",
    index = seq_along(z_poisson),
    z = z_poisson
  )
)

make_summary <- function(z, model_name) {
  ks <- ks.test(z, "pexp", rate = 1)
  
  lb5 <- Box.test(
    z,
    lag = 5,
    type = "Ljung-Box"
  )
  
  lb10 <- Box.test(
    z,
    lag = 10,
    type = "Ljung-Box"
  )
  
  acf_values <- as.numeric(
    acf(
      z,
      lag.max = 10,
      plot = FALSE
    )$acf
  )
  
  data.frame(
    model = model_name,
    n = length(z),
    mean_z = mean(z),
    var_z = var(z),
    median_z = median(z),
    ks_statistic = as.numeric(ks$statistic),
    ks_p_value = as.numeric(ks$p.value),
    acf_lag1 = acf_values[2],
    acf_lag5 = acf_values[6],
    ljung_box_lag5_statistic = as.numeric(lb5$statistic),
    ljung_box_lag5_p_value = as.numeric(lb5$p.value),
    ljung_box_lag10_statistic = as.numeric(lb10$statistic),
    ljung_box_lag10_p_value = as.numeric(lb10$p.value)
  )
}

summary_df <- rbind(
  make_summary(z_hawkes, "Hawkes_correct"),
  make_summary(z_poisson, "Poisson_wrong")
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
  "results/ogata_residual_fit_summary.csv",
  row.names = FALSE
)

write.csv(
  summary_df,
  "results/ogata_residual_independence_summary.csv",
  row.names = FALSE
)

write.csv(
  residual_df,
  "results/ogata_residuals.csv",
  row.names = FALSE
)

make_acf_df <- function(z, model_name, lag_max = 30) {
  a <- acf(z, lag.max = lag_max, plot = FALSE)

  data.frame(
    model = model_name,
    lag = as.numeric(a$lag),
    acf = as.numeric(a$acf)
  )
}

acf_df <- rbind(
  make_acf_df(z_hawkes, "Hawkes_correct"),
  make_acf_df(z_poisson, "Poisson_wrong")
)

write.csv(
  acf_df,
  "results/ogata_residual_acf.csv",
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
    title = "Ogata residual analysis: histogram",
    subtitle = "Correct model should match Exp(1)",
    x = "rescaled gap",
    y = "density"
  ) +
  theme_minimal()

ggsave(
  "figures/ogata_residual_histogram.pdf",
  p_hist,
  width = 8,
  height = 4
)

ggsave(
  "figures/ogata_residual_histogram.png",
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
    title = "Ogata residual analysis: QQ plot",
    subtitle = "Correct model should follow the diagonal",
    x = "theoretical Exp(1) quantiles",
    y = "empirical quantiles"
  ) +
  theme_minimal()

ggsave(
  "figures/ogata_residual_qqplot.pdf",
  p_qq,
  width = 8,
  height = 4
)

ggsave(
  "figures/ogata_residual_qqplot.png",
  p_qq,
  width = 8,
  height = 4,
  dpi = 300
)

lag_df <- rbind(
  data.frame(
    model = "Hawkes_correct",
    z_i = z_hawkes[-length(z_hawkes)],
    z_next = z_hawkes[-1]
  ),
  data.frame(
    model = "Poisson_wrong",
    z_i = z_poisson[-length(z_poisson)],
    z_next = z_poisson[-1]
  )
)

p_lag <- ggplot(
  lag_df,
  aes(
    x = z_i,
    y = z_next
  )
) +
  geom_point(alpha = 0.35) +
  facet_wrap(~ model) +
  labs(
    title = "Ogata residual analysis: lag-1 scatter",
    subtitle = "Independent residuals should show no visible structure",
    x = expression(z[i]),
    y = expression(z[i + 1])
  ) +
  theme_minimal()

ggsave(
  "figures/ogata_residual_lag_scatter.pdf",
  p_lag,
  width = 8,
  height = 4
)

ggsave(
  "figures/ogata_residual_lag_scatter.png",
  p_lag,
  width = 8,
  height = 4,
  dpi = 300
)

acf_plot_df <- acf_df[acf_df$lag > 0, ]

p_acf <- ggplot(
  acf_plot_df,
  aes(
    x = lag,
    y = acf
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_segment(
    aes(
      xend = lag,
      y = 0,
      yend = acf
    )
  ) +
  geom_point(size = 2) +
  facet_wrap(~ model) +
  labs(
    title = "Ogata residual analysis: autocorrelation",
    subtitle = "Correct model should have autocorrelations near zero",
    x = "lag",
    y = "ACF"
  ) +
  theme_minimal()

ggsave(
  "figures/ogata_residual_acf.pdf",
  p_acf,
  width = 8,
  height = 4
)

ggsave(
  "figures/ogata_residual_acf.png",
  p_acf,
  width = 8,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/ogata_residual_fit_summary.csv\n")
cat("- results/ogata_residual_independence_summary.csv\n")
cat("- results/ogata_residuals.csv\n")
cat("- results/ogata_residual_acf.csv\n")
cat("- figures/ogata_residual_histogram.pdf/png\n")
cat("- figures/ogata_residual_qqplot.pdf/png\n")
cat("- figures/ogata_residual_lag_scatter.pdf/png\n")
cat("- figures/ogata_residual_acf.pdf/png\n")
