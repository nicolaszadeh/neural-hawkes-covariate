set.seed(4)

source("R/hawkes_simulation.R")
source("R/hawkes_likelihood.R")
source("R/hawkes_fitting.R")
source("R/hawkes_bootstrap.R")
source("R/hawkes_plotting.R")

# ============================================================
# 1. Parameters
# ============================================================

T_end <- 250
dt <- 0.01
B_boot <- 100

time <- seq(0, T_end, by = dt)

gamma_true <- 0.8
sigma_true <- 1.0

mu0_true <- -0.5
mu1_true <- 0.7
alpha_true <- 0.0
beta_true <- 2.5

cat("Parametric bootstrap for Hawkes test\n")
cat("T_end =", T_end, "\n")
cat("dt =", dt, "\n")
cat("B_boot =", B_boot, "\n\n")

cat("True parameters:\n")
print(c(
  mu0 = mu0_true,
  mu1 = mu1_true,
  alpha = alpha_true,
  beta = beta_true
))

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# ============================================================
# 2. Simulate observed no-Hawkes dataset
# ============================================================

X <- simulate_ou(
  time = time,
  gamma = gamma_true,
  sigma = sigma_true,
  stationary_start = TRUE
)

sim_obs <- simulate_hawkes_discrete(
  X = X,
  dt = dt,
  mu0 = mu0_true,
  mu1 = mu1_true,
  alpha = alpha_true,
  beta = beta_true
)

dN_obs <- sim_obs$dN

cat("\nObserved dataset:\n")
cat("Number of events:", sum(dN_obs), "\n")
cat("Max number of events in one bin:", max(dN_obs), "\n")


# ============================================================
# 3. Observed Hawkes LR test
# ============================================================

observed_test <- fit_hawkes_test(
  X = X,
  dN = dN_obs,
  dt = dt,
  hessian = TRUE
)

LR_obs <- observed_test$LR

cat("\nObserved test result:\n")
print(observed_test)

cat("\nNaive chi-square p-value:\n")
print(observed_test$p_chisq)


# ============================================================
# 4. Fitted null parameters for bootstrap
# ============================================================

mu0_null_hat <- observed_test$mu0_null
mu1_null_hat <- observed_test$mu1_null

cat("\nFitted null parameters for bootstrap:\n")
print(c(
  mu0_null_hat = mu0_null_hat,
  mu1_null_hat = mu1_null_hat
))


# ============================================================
# 5. Parametric bootstrap under H0: alpha = 0
# ============================================================

boot <- run_hawkes_parametric_bootstrap(
  X = X,
  dt = dt,
  mu0_null_hat = mu0_null_hat,
  mu1_null_hat = mu1_null_hat,
  LR_obs = LR_obs,
  B_boot = B_boot,
  verbose = TRUE
)

cat("\nBootstrap p-value:\n")
print(boot$p_boot)

cat("\nNumber of successful bootstrap fits:\n")
cat(boot$n_success, "out of", B_boot, "\n")

cat("\n95 percent critical values:\n")
cat("bootstrap:", boot$crit_boot_95, "\n")
cat("chi-square:", boot$crit_chisq_95, "\n")


# ============================================================
# 6. Save results
# ============================================================

summary_table <- data.frame(
  observed_test,
  p_boot = boot$p_boot,
  B_boot = B_boot,
  n_success = boot$n_success,
  crit_boot_95 = boot$crit_boot_95,
  crit_chisq_95 = boot$crit_chisq_95,
  mu0_true = mu0_true,
  mu1_true = mu1_true,
  alpha_true = alpha_true,
  beta_true = beta_true
)

write.csv(
  summary_table,
  file = "results/hawkes_bootstrap_summary.csv",
  row.names = FALSE
)

write.csv(
  boot$boot_results,
  file = "results/hawkes_bootstrap_replicates.csv",
  row.names = FALSE
)

cat("\nSaved results:\n")
cat("results/hawkes_bootstrap_summary.csv\n")
cat("results/hawkes_bootstrap_replicates.csv\n")


# ============================================================
# 7. Bootstrap LR histogram
# ============================================================

plot_bootstrap_histogram <- function() {
  hist(
    boot$boot_LR,
    breaks = 20,
    xlab = "bootstrap LR statistic",
    main = "Bootstrap null distribution of Hawkes LR"
  )

  abline(v = LR_obs, lwd = 2)
  abline(v = boot$crit_boot_95, lty = 2, lwd = 2)
  abline(v = boot$crit_chisq_95, lty = 3, lwd = 2)

  legend(
    "topright",
    legend = c(
      "observed LR",
      "bootstrap 95 percent",
      "chi-square 95 percent"
    ),
    lty = c(1, 2, 3),
    lwd = c(2, 2, 2),
    bty = "n"
  )
}

plot_bootstrap_histogram()

save_plot_both(
  file_stem = "figures/hawkes_bootstrap_LR_histogram",
  plot_function = plot_bootstrap_histogram,
  png_width = 1000,
  png_height = 800,
  pdf_width = 8,
  pdf_height = 6
)


# ============================================================
# 8. Final message
# ============================================================

cat("\nDone.\n")
cat("Generated files:\n")
cat("results/hawkes_bootstrap_summary.csv\n")
cat("results/hawkes_bootstrap_replicates.csv\n")
cat("figures/hawkes_bootstrap_LR_histogram.png\n")
cat("figures/hawkes_bootstrap_LR_histogram.pdf\n")