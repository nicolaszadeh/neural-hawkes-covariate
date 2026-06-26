set.seed(3)

source("R/load_all.R")

# ============================================================
# 1. Parameters
# ============================================================

T_end <- 250
dt <- 0.01
n_rep <- 100

time <- seq(0, T_end, by = dt)

gamma_true <- 0.8
sigma_true <- 1.0

mu0_true <- -0.5
beta_true <- 2.5

scenarios <- data.frame(
  scenario = c("no_covariate", "no_hawkes"),
  mu0_true = c(mu0_true, mu0_true),
  mu1_true = c(0.0, 0.7),
  alpha_true = c(0.7, 0.0),
  beta_true = c(beta_true, beta_true)
)

cat("Null calibration study\n")
cat("T_end =", T_end, "\n")
cat("dt =", dt, "\n")
cat("n_rep =", n_rep, "\n\n")

print(scenarios)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# ============================================================
# 2. Run null calibration
# ============================================================

all_results <- list()
counter <- 1

for (s in seq_len(nrow(scenarios))) {
  scenario_name <- scenarios$scenario[s]

  mu0_s <- scenarios$mu0_true[s]
  mu1_s <- scenarios$mu1_true[s]
  alpha_s <- scenarios$alpha_true[s]
  beta_s <- scenarios$beta_true[s]

  cat("\nScenario:", scenario_name, "\n")
  cat("mu0_true =", mu0_s, "\n")
  cat("mu1_true =", mu1_s, "\n")
  cat("alpha_true =", alpha_s, "\n")
  cat("beta_true =", beta_s, "\n\n")

  for (rep_id in 1:n_rep) {
    cat(
      "Running",
      scenario_name,
      "| repetition",
      rep_id,
      "\n"
    )

    X <- simulate_ou(
      time = time,
      gamma = gamma_true,
      sigma = sigma_true,
      stationary_start = TRUE
    )

    sim <- simulate_hawkes_discrete(
      X = X,
      dt = dt,
      mu0 = mu0_s,
      mu1 = mu1_s,
      alpha = alpha_s,
      beta = beta_s
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

    fit_result$scenario <- scenario_name
    fit_result$rep_id <- rep_id

    fit_result$mu0_true <- mu0_s
    fit_result$mu1_true <- mu1_s
    fit_result$alpha_true <- alpha_s
    fit_result$beta_true <- beta_s

    all_results[[counter]] <- fit_result
    counter <- counter + 1
  }
}

results <- do.call(rbind, all_results)

results <- results[
  ,
  c(
    "scenario",
    "rep_id",
    "mu0_true",
    "mu1_true",
    "alpha_true",
    "beta_true",
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
  file = "results/null_calibration_estimates.csv",
  row.names = FALSE
)

cat("\nSaved raw results to:\n")
cat("results/null_calibration_estimates.csv\n")


# ============================================================
# 3. Summary table
# ============================================================

ok <- complete.cases(
  results[, c("mu0_hat", "mu1_hat", "alpha_hat", "beta_hat")]
)

results_ok <- results[ok, ]

make_summary_one_scenario <- function(df) {
  data.frame(
    scenario = unique(df$scenario),
    n_runs = nrow(df),

    mu1_true = unique(df$mu1_true),
    alpha_true = unique(df$alpha_true),

    mean_events = mean(df$n_events),
    sd_events = sd(df$n_events),

    mean_mu0 = mean(df$mu0_hat),
    mean_mu1 = mean(df$mu1_hat),
    mean_alpha = mean(df$alpha_hat),
    mean_beta = mean(df$beta_hat),

    sd_mu0 = sd(df$mu0_hat),
    sd_mu1 = sd(df$mu1_hat),
    sd_alpha = sd(df$alpha_hat),
    sd_beta = sd(df$beta_hat),

    rejection_cov_5pct = mean(df$p_cov < 0.05),
    rejection_hawkes_5pct = mean(df$p_hawkes < 0.05),

    median_p_cov = median(df$p_cov),
    median_p_hawkes = median(df$p_hawkes),

    conv_full_rate = mean(df$conv_full == 0),
    conv_no_cov_rate = mean(df$conv_no_cov == 0),
    conv_no_hawkes_rate = mean(df$conv_no_hawkes == 0)
  )
}

summary_table <- do.call(
  rbind,
  lapply(
    split(results_ok, results_ok$scenario),
    make_summary_one_scenario
  )
)

write.csv(
  summary_table,
  file = "results/null_calibration_summary.csv",
  row.names = FALSE
)

cat("\nSummary table:\n")
print(summary_table)

cat("\nSaved summary to:\n")
cat("results/null_calibration_summary.csv\n")


# ============================================================
# 4. Rejection rate plot
# ============================================================

plot_rejection_rates <- function() {
  x <- seq_len(nrow(summary_table))

  ylim_max <- max(
    summary_table$rejection_cov_5pct,
    summary_table$rejection_hawkes_5pct,
    1
  )

  plot(
    x,
    summary_table$rejection_cov_5pct,
    type = "b",
    xaxt = "n",
    ylim = c(0, ylim_max),
    xlab = "scenario",
    ylab = "rejection rate",
    main = "Null calibration: rejection rates at 5 percent"
  )

  lines(
    x,
    summary_table$rejection_hawkes_5pct,
    type = "b",
    lty = 2
  )

  abline(h = 0.05, lty = 3)

  axis(
    side = 1,
    at = x,
    labels = summary_table$scenario
  )

  legend(
    "topright",
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
  file_stem = "figures/null_calibration_rejection_rates",
  plot_function = plot_rejection_rates
)


# ============================================================
# 5. P-value histograms
# ============================================================

plot_pvalue_histograms <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  for (scenario_name in scenarios$scenario) {
    df <- results_ok[results_ok$scenario == scenario_name, ]

    hist(
      df$p_cov,
      breaks = 15,
      xlim = c(0, 1),
      xlab = "p-value",
      main = paste("Covariate test:", scenario_name)
    )

    abline(v = 0.05, lty = 2)

    hist(
      df$p_hawkes,
      breaks = 15,
      xlim = c(0, 1),
      xlab = "p-value",
      main = paste("Hawkes test:", scenario_name)
    )

    abline(v = 0.05, lty = 2)
  }
}

plot_pvalue_histograms()

save_plot_both(
  file_stem = "figures/null_calibration_pvalue_histograms",
  plot_function = plot_pvalue_histograms,
  png_width = 1200,
  png_height = 900,
  pdf_width = 10,
  pdf_height = 7
)


# ============================================================
# 6. Estimate boxplots
# ============================================================

plot_estimate_boxplots <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  boxplot(
    mu1_hat ~ scenario,
    data = results_ok,
    xlab = "scenario",
    ylab = "estimate",
    main = "Estimates of mu1"
  )

  points(
    x = 1:nrow(scenarios),
    y = scenarios$mu1_true,
    pch = 4,
    cex = 1.5,
    lwd = 2
  )

  boxplot(
    alpha_hat ~ scenario,
    data = results_ok,
    xlab = "scenario",
    ylab = "estimate",
    main = "Estimates of alpha"
  )

  points(
    x = 1:nrow(scenarios),
    y = scenarios$alpha_true,
    pch = 4,
    cex = 1.5,
    lwd = 2
  )

  boxplot(
    beta_hat ~ scenario,
    data = results_ok,
    xlab = "scenario",
    ylab = "estimate",
    main = "Estimates of beta"
  )

  points(
    x = 1:nrow(scenarios),
    y = scenarios$beta_true,
    pch = 4,
    cex = 1.5,
    lwd = 2
  )

  plot.new()

  legend(
    "center",
    legend = c("cross = true value"),
    pch = 4,
    pt.cex = 1.5,
    bty = "n"
  )
}

plot_estimate_boxplots()

save_plot_both(
  file_stem = "figures/null_calibration_estimate_boxplots",
  plot_function = plot_estimate_boxplots,
  png_width = 1200,
  png_height = 900,
  pdf_width = 10,
  pdf_height = 7
)


# ============================================================
# 7. Final message
# ============================================================

cat("\nDone.\n")
cat("Generated files:\n")
cat("results/null_calibration_estimates.csv\n")
cat("results/null_calibration_summary.csv\n")
cat("figures/null_calibration_rejection_rates.png\n")
cat("figures/null_calibration_rejection_rates.pdf\n")
cat("figures/null_calibration_pvalue_histograms.png\n")
cat("figures/null_calibration_pvalue_histograms.pdf\n")
cat("figures/null_calibration_estimate_boxplots.png\n")
cat("figures/null_calibration_estimate_boxplots.pdf\n")
