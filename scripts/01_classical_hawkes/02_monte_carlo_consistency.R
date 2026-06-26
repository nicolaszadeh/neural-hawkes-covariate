set.seed(2)

source("R/load_all.R")

# ============================================================
# 1. Parameters
# ============================================================

dt <- 0.01
T_values <- c(100, 250, 500)
n_rep <- 20

gamma_true <- 0.8
sigma_true <- 1.0

mu0_true <- -0.5
mu1_true <- 0.7
alpha_true <- 0.7
beta_true <- 2.5

cat("alpha / beta =", alpha_true / beta_true, "\n")
cat("Number of repetitions per T:", n_rep, "\n")
cat("T values:", T_values, "\n")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# ============================================================
# 2. Repeated simulation study
# ============================================================

all_results <- list()
counter <- 1

for (T_end in T_values) {
  for (rep_id in 1:n_rep) {
    cat("Running T =", T_end, "| repetition", rep_id, "\n")

    time <- seq(0, T_end, by = dt)

    X <- simulate_ou(
      time = time,
      gamma = gamma_true,
      sigma = sigma_true,
      stationary_start = TRUE
    )

    sim <- simulate_hawkes_discrete(
      X = X,
      dt = dt,
      mu0 = mu0_true,
      mu1 = mu1_true,
      alpha = alpha_true,
      beta = beta_true
    )

    fit_result <- tryCatch(
      {
        fit_all_models_and_tests(
          X = X,
          dN = sim$dN,
          dt = dt,
          hessian = FALSE
        )
      },
      error = function(e) {
        cat("Error in fit:", conditionMessage(e), "\n")

        data.frame(
          n_events = NA,
          max_dN = NA,
          mu0_hat = NA,
          mu1_hat = NA,
          alpha_hat = NA,
          beta_hat = NA,
          loglik_full = NA,
          loglik_no_cov = NA,
          loglik_no_hawkes = NA,
          LR_cov = NA,
          p_cov = NA,
          LR_hawkes = NA,
          p_hawkes = NA,
          conv_full = NA,
          conv_no_cov = NA,
          conv_no_hawkes = NA
        )
      }
    )

    fit_result$T_end <- T_end
    fit_result$rep_id <- rep_id

    all_results[[counter]] <- fit_result
    counter <- counter + 1
  }
}

results <- do.call(rbind, all_results)

results <- results[
  ,
  c(
    "T_end",
    "rep_id",
    "n_events",
    "max_dN",
    "mu0_hat",
    "mu1_hat",
    "alpha_hat",
    "beta_hat",
    "loglik_full",
    "loglik_no_cov",
    "loglik_no_hawkes",
    "LR_cov",
    "p_cov",
    "LR_hawkes",
    "p_hawkes",
    "conv_full",
    "conv_no_cov",
    "conv_no_hawkes"
  )
]

write.csv(
  results,
  file = "results/monte_carlo_estimates.csv",
  row.names = FALSE
)

cat("\nSaved raw results to:\n")
cat("results/monte_carlo_estimates.csv\n")


# ============================================================
# 3. Summary table
# ============================================================

ok <- complete.cases(
  results[, c("mu0_hat", "mu1_hat", "alpha_hat", "beta_hat")]
)

results_ok <- results[ok, ]

make_summary_one_T <- function(df) {
  data.frame(
    T_end = unique(df$T_end),
    n_runs = nrow(df),

    mean_events = mean(df$n_events),
    sd_events = sd(df$n_events),

    mean_mu0 = mean(df$mu0_hat),
    sd_mu0 = sd(df$mu0_hat),

    mean_mu1 = mean(df$mu1_hat),
    sd_mu1 = sd(df$mu1_hat),

    mean_alpha = mean(df$alpha_hat),
    sd_alpha = sd(df$alpha_hat),

    mean_beta = mean(df$beta_hat),
    sd_beta = sd(df$beta_hat),

    rejection_cov_5pct = mean(df$p_cov < 0.05),
    rejection_hawkes_5pct = mean(df$p_hawkes < 0.05),

    conv_full_rate = mean(df$conv_full == 0),
    conv_no_cov_rate = mean(df$conv_no_cov == 0),
    conv_no_hawkes_rate = mean(df$conv_no_hawkes == 0)
  )
}

summary_table <- do.call(
  rbind,
  lapply(
    split(results_ok, results_ok$T_end),
    make_summary_one_T
  )
)

write.csv(
  summary_table,
  file = "results/monte_carlo_summary.csv",
  row.names = FALSE
)

cat("\nSummary table:\n")
print(summary_table)

cat("\nSaved summary to:\n")
cat("results/monte_carlo_summary.csv\n")


# ============================================================
# 4. Boxplots of parameter estimates
# ============================================================

plot_parameter_boxplots <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  boxplot(
    mu0_hat ~ T_end,
    data = results_ok,
    xlab = "T",
    ylab = "estimate",
    main = "Estimates of mu0"
  )
  abline(h = mu0_true, lty = 2)

  boxplot(
    mu1_hat ~ T_end,
    data = results_ok,
    xlab = "T",
    ylab = "estimate",
    main = "Estimates of mu1"
  )
  abline(h = mu1_true, lty = 2)

  boxplot(
    alpha_hat ~ T_end,
    data = results_ok,
    xlab = "T",
    ylab = "estimate",
    main = "Estimates of alpha"
  )
  abline(h = alpha_true, lty = 2)

  boxplot(
    beta_hat ~ T_end,
    data = results_ok,
    xlab = "T",
    ylab = "estimate",
    main = "Estimates of beta"
  )
  abline(h = beta_true, lty = 2)
}

plot_parameter_boxplots()

save_plot_both(
  file_stem = "figures/monte_carlo_parameter_boxplots",
  plot_function = plot_parameter_boxplots
)


# ============================================================
# 5. Mean estimates as a function of T
# ============================================================

plot_mean_estimates <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  plot(
    summary_table$T_end,
    summary_table$mean_mu0,
    type = "b",
    xlab = "T",
    ylab = "mean estimate",
    main = "Mean estimate of mu0"
  )
  abline(h = mu0_true, lty = 2)

  plot(
    summary_table$T_end,
    summary_table$mean_mu1,
    type = "b",
    xlab = "T",
    ylab = "mean estimate",
    main = "Mean estimate of mu1"
  )
  abline(h = mu1_true, lty = 2)

  plot(
    summary_table$T_end,
    summary_table$mean_alpha,
    type = "b",
    xlab = "T",
    ylab = "mean estimate",
    main = "Mean estimate of alpha"
  )
  abline(h = alpha_true, lty = 2)

  plot(
    summary_table$T_end,
    summary_table$mean_beta,
    type = "b",
    xlab = "T",
    ylab = "mean estimate",
    main = "Mean estimate of beta"
  )
  abline(h = beta_true, lty = 2)
}

plot_mean_estimates()

save_plot_both(
  file_stem = "figures/monte_carlo_mean_estimates",
  plot_function = plot_mean_estimates
)


# ============================================================
# 6. Rejection rates
# ============================================================

plot_rejection_rates <- function() {
  ylim_max <- max(
    summary_table$rejection_cov_5pct,
    summary_table$rejection_hawkes_5pct,
    1
  )

  plot(
    summary_table$T_end,
    summary_table$rejection_cov_5pct,
    type = "b",
    ylim = c(0, ylim_max),
    xlab = "T",
    ylab = "rejection rate",
    main = "Rejection rate at 5 percent"
  )

  lines(
    summary_table$T_end,
    summary_table$rejection_hawkes_5pct,
    type = "b",
    lty = 2
  )

  abline(h = 0.05, lty = 3)

  legend(
    "bottomright",
    legend = c(
      "covariate test",
      "Hawkes test",
      "5 percent reference"
    ),
    lty = c(1, 2, 3),
    bty = "n"
  )
}

plot_rejection_rates()

save_plot_both(
  file_stem = "figures/monte_carlo_rejection_rates",
  plot_function = plot_rejection_rates
)


# ============================================================
# 7. Final message
# ============================================================

cat("\nDone.\n")
cat("Generated files:\n")
cat("results/monte_carlo_estimates.csv\n")
cat("results/monte_carlo_summary.csv\n")
cat("figures/monte_carlo_parameter_boxplots.png\n")
cat("figures/monte_carlo_parameter_boxplots.pdf\n")
cat("figures/monte_carlo_mean_estimates.png\n")
cat("figures/monte_carlo_mean_estimates.pdf\n")
cat("figures/monte_carlo_rejection_rates.png\n")
cat("figures/monte_carlo_rejection_rates.pdf\n")
