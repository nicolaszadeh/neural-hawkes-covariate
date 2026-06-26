set.seed(1)

source("R/load_all.R")

# ============================================================
# 1. Parameters
# ============================================================

T_end <- 500
dt <- 0.01
time <- seq(0, T_end, by = dt)

gamma_true <- 0.8
sigma_true <- 1.0

mu0_true <- -0.5
mu1_true <- 0.7
alpha_true <- 0.7
beta_true <- 2.5

cat("alpha / beta =", alpha_true / beta_true, "\n")


# ============================================================
# 2. Simulate data
# ============================================================

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

dN <- sim$dN
lambda_true <- sim$lambda

cat("Number of events:", sum(dN), "\n")
cat("Max number of events in one bin:", max(dN), "\n")


# ============================================================
# 3. Fit models and compute LR tests
# ============================================================

fit_result <- fit_all_models_and_tests(
  X = X,
  dN = dN,
  dt = dt,
  hessian = TRUE
)

est_full <- c(
  mu0 = fit_result$mu0_hat,
  mu1 = fit_result$mu1_hat,
  alpha = fit_result$alpha_hat,
  beta = fit_result$beta_hat
)

cat("\nFull model estimates:\n")
print(est_full)

cat("\nTrue values:\n")
print(c(
  mu0 = mu0_true,
  mu1 = mu1_true,
  alpha = alpha_true,
  beta = beta_true
))

cat("\nTest of covariate relevance:\n")
cat("H0: mu1 = 0\n")
cat("LR statistic =", fit_result$LR_cov, "\n")
cat("p-value =", fit_result$p_cov, "\n")

cat("\nTest of Hawkes self-excitation:\n")
cat("H0: alpha = 0\n")
cat("LR statistic =", fit_result$LR_hawkes, "\n")
cat("p-value =", fit_result$p_hawkes, "\n")
cat("Warning: this p-value is approximate because alpha = 0\n")
cat("is on the boundary of the parameter space.\n")

cat("\nOptimization convergence codes:\n")
cat("Full model:", fit_result$conv_full, "\n")
cat("No covariate:", fit_result$conv_no_cov, "\n")
cat("No Hawkes:", fit_result$conv_no_hawkes, "\n")


# ============================================================
# 4. Reconstruct fitted intensity
# ============================================================

lambda_fit <- compute_intensity(
  X = X,
  dN = dN,
  dt = dt,
  mu0 = est_full["mu0"],
  mu1 = est_full["mu1"],
  alpha = est_full["alpha"],
  beta = est_full["beta"]
)


# ============================================================
# 5. Intensity error summaries
# ============================================================

lambda_error <- lambda_fit - lambda_true
relative_error <- lambda_error / pmax(lambda_true, 1e-12)

cat("\nIntensity error summaries:\n")
cat("Mean absolute error:",
    mean(abs(lambda_error)), "\n")
cat("Root mean squared error:",
    sqrt(mean(lambda_error^2)), "\n")
cat("Mean absolute relative error:",
    mean(abs(relative_error)), "\n")

cat("\nRelative error quantiles:\n")
print(quantile(
  relative_error,
  probs = c(0.01, 0.05, 0.5, 0.95, 0.99)
))


# ============================================================
# 6. Main diagnostic plot
# ============================================================

plot_main_diagnostics <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

  plot(
    time,
    X,
    type = "l",
    xlab = "time",
    ylab = "X(t)",
    main = "Exact OU covariate at grid points"
  )

  plot(
    time,
    lambda_true,
    type = "l",
    lwd = 2,
    xlab = "time",
    ylab = "lambda(t)",
    main = "True and fitted intensity"
  )

  lines(time, lambda_fit, lty = 2, lwd = 1)

  legend(
    "topright",
    legend = c("true", "fitted"),
    lty = c(1, 2),
    lwd = c(2, 1),
    bty = "n"
  )

  plot(
    time,
    dN,
    type = "h",
    xlab = "time",
    ylab = "dN",
    main = "Observed events"
  )
}

plot_main_diagnostics()

save_plot_both(
  file_stem = "figures/main_diagnostics",
  plot_function = plot_main_diagnostics
)


# ============================================================
# 7. Zoomed diagnostic plot
# ============================================================

plot_zoom_diagnostics <- function(t_max = 50) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  idx <- time <= t_max

  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

  plot(
    time[idx],
    X[idx],
    type = "l",
    xlab = "time",
    ylab = "X(t)",
    main = "Zoom: OU covariate"
  )

  plot(
    time[idx],
    lambda_true[idx],
    type = "l",
    lwd = 2,
    xlab = "time",
    ylab = "lambda(t)",
    main = "Zoom: true and fitted intensity"
  )

  lines(time[idx], lambda_fit[idx], lty = 2, lwd = 1)

  legend(
    "topright",
    legend = c("true", "fitted"),
    lty = c(1, 2),
    lwd = c(2, 1),
    bty = "n"
  )

  plot(
    time[idx],
    dN[idx],
    type = "h",
    xlab = "time",
    ylab = "dN",
    main = "Zoom: observed events"
  )
}

plot_zoom_diagnostics(t_max = 50)

save_plot_both(
  file_stem = "figures/zoom_diagnostics",
  plot_function = function() plot_zoom_diagnostics(t_max = 50)
)


# ============================================================
# 8. Intensity error plots
# ============================================================

plot_intensity_errors <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

  plot(
    time,
    lambda_error,
    type = "l",
    xlab = "time",
    ylab = "fit - true",
    main = "Difference between fitted and true intensity"
  )

  abline(h = 0, lty = 2)

  plot(
    time,
    relative_error,
    type = "l",
    xlab = "time",
    ylab = "relative error",
    main = "Relative error of fitted intensity"
  )

  abline(h = 0, lty = 2)
}

plot_intensity_errors()

save_plot_both(
  file_stem = "figures/intensity_errors",
  plot_function = plot_intensity_errors
)


# ============================================================
# 9. True versus fitted intensity
# ============================================================

plot_intensity_scatter <- function() {
  plot(
    lambda_true,
    lambda_fit,
    pch = 16,
    cex = 0.3,
    xlab = "true intensity",
    ylab = "fitted intensity",
    main = "Fitted intensity versus true intensity"
  )

  abline(0, 1, lty = 2)
}

plot_intensity_scatter()

save_plot_both(
  file_stem = "figures/intensity_scatter",
  plot_function = plot_intensity_scatter
)


# ============================================================
# 10. Final message
# ============================================================

cat("\nDone.\n")
cat("Generated files:\n")
cat("figures/main_diagnostics.png\n")
cat("figures/main_diagnostics.pdf\n")
cat("figures/zoom_diagnostics.png\n")
cat("figures/zoom_diagnostics.pdf\n")
cat("figures/intensity_errors.png\n")
cat("figures/intensity_errors.pdf\n")
cat("figures/intensity_scatter.png\n")
cat("figures/intensity_scatter.pdf\n")
