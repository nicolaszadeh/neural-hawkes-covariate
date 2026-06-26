# 25_bootstrap_gof_time_rescaling.R
#
# Refactored version.
# Parametric bootstrap goodness-of-fit test based on
# time-rescaled residuals.

set.seed(25)

source("R/load_all.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

T_end <- 1000
B <- 100

mu_true <- 0.8
alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

ks_statistic_residuals <- function(events, fit_out) {
  z <- hawkes_rescaled_gaps(
    events = events,
    mu = fit_out$mu,
    alpha = fit_out$alpha,
    beta = fit_out$beta
  )

  ks <- ks.test(z, "pexp", rate = 1)

  list(
    D = as.numeric(ks$statistic),
    raw_ks_p_value = ks$p.value,
    mean_z = mean(z),
    var_z = var(z),
    z = z
  )
}

# ------------------------------------------------------------
# Observed dataset
# ------------------------------------------------------------

events_obs <- simulate_hawkes(
  T_end = T_end,
  mu = mu_true,
  alpha = alpha_true,
  beta = beta_true
)

fit_obs <- fit_hawkes_continuous(
  events_obs,
  T_end,
  method = "BFGS"
)

gof_obs <- ks_statistic_residuals(events_obs, fit_obs)

cat("Observed n:", length(events_obs), "\n")
cat("Observed KS statistic:", gof_obs$D, "\n")
cat("Observed raw KS p-value:", gof_obs$raw_ks_p_value, "\n")
cat("Observed branching estimate:", fit_obs$branching, "\n")

# ------------------------------------------------------------
# Parametric bootstrap
# ------------------------------------------------------------

boot_rows <- list()

for (b in seq_len(B)) {
  boot_events <- simulate_hawkes(
    T_end = T_end,
    mu = fit_obs$mu,
    alpha = fit_obs$alpha,
    beta = fit_obs$beta
  )

  boot_fit <- fit_hawkes_continuous(
    boot_events,
    T_end,
    method = "BFGS"
  )

  boot_gof <- ks_statistic_residuals(
    boot_events,
    boot_fit
  )

  boot_rows[[b]] <- data.frame(
    bootstrap = b,
    n_events = length(boot_events),
    D = boot_gof$D,
    raw_ks_p_value = boot_gof$raw_ks_p_value,
    mean_z = boot_gof$mean_z,
    var_z = boot_gof$var_z,
    mu_hat = boot_fit$mu,
    alpha_hat = boot_fit$alpha,
    beta_hat = boot_fit$beta,
    branching_hat = boot_fit$branching,
    convergence = boot_fit$fit$convergence
  )

  cat(
    "Bootstrap", b, "of", B,
    "- n =", length(boot_events),
    "- D =", round(boot_gof$D, 4),
    "- raw p =", round(boot_gof$raw_ks_p_value, 4),
    "\n"
  )
}

boot_df <- do.call(rbind, boot_rows)

n_exceed <- sum(boot_df$D >= gof_obs$D)

p_empirical <- n_exceed / B
p_corrected <- (1 + n_exceed) / (B + 1)

summary_df <- data.frame(
  observed_n_events = length(events_obs),
  observed_D = gof_obs$D,
  observed_raw_ks_p_value = gof_obs$raw_ks_p_value,
  bootstrap_p_empirical = p_empirical,
  bootstrap_p_corrected = p_corrected,
  observed_mean_z = gof_obs$mean_z,
  observed_var_z = gof_obs$var_z,
  mu_true = mu_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = branching_true,
  mu_hat = fit_obs$mu,
  alpha_hat = fit_obs$alpha,
  beta_hat = fit_obs$beta,
  branching_hat = fit_obs$branching,
  B = B
)

print(summary_df)

write.csv(
  summary_df,
  "results/bootstrap_gof_time_rescaling_summary.csv",
  row.names = FALSE
)

write.csv(
  boot_df,
  "results/bootstrap_gof_time_rescaling_replicates.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    index = seq_along(gof_obs$z),
    z = gof_obs$z
  ),
  "results/bootstrap_gof_time_rescaling_observed_residuals.csv",
  row.names = FALSE
)

p_D <- ggplot(boot_df, aes(x = D)) +
  geom_histogram(bins = 30, boundary = 0) +
  geom_vline(
    xintercept = gof_obs$D,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(
    title = "Parametric bootstrap GOF test",
    subtitle = "Dashed line: observed KS statistic",
    x = "KS statistic after refitting",
    y = "count"
  ) +
  theme_minimal()

ggsave(
  "figures/bootstrap_gof_time_rescaling_KS_histogram.pdf",
  p_D,
  width = 7,
  height = 4
)

ggsave(
  "figures/bootstrap_gof_time_rescaling_KS_histogram.png",
  p_D,
  width = 7,
  height = 4,
  dpi = 300
)

obs_residual_df <- data.frame(z = gof_obs$z)

p_resid <- ggplot(obs_residual_df, aes(x = z)) +
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
    title = "Observed fitted residuals",
    subtitle = "Histogram compared with Exp(1) density",
    x = "rescaled gap",
    y = "density"
  ) +
  theme_minimal()

ggsave(
  "figures/bootstrap_gof_time_rescaling_observed_residuals.pdf",
  p_resid,
  width = 7,
  height = 4
)

ggsave(
  "figures/bootstrap_gof_time_rescaling_observed_residuals.png",
  p_resid,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/bootstrap_gof_time_rescaling_summary.csv\n")
cat("- results/bootstrap_gof_time_rescaling_replicates.csv\n")
cat("- results/bootstrap_gof_time_rescaling_observed_residuals.csv\n")
cat("- figures/bootstrap_gof_time_rescaling_KS_histogram.pdf/png\n")
cat("- figures/bootstrap_gof_time_rescaling_observed_residuals.pdf/png\n")
