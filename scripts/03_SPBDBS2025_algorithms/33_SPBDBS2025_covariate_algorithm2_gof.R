# ============================================================
# SPBDBS2025 - Covariate Algorithm 2 GOF
# ============================================================

set.seed(3301)

source("R/load_all.R")
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

T_end <- 500
dt <- 0.01
grid <- seq(0, T_end, by = dt)

gamma0_true <- -0.5
gamma1_true <- 0.9
alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

kappa_true <- 0.4
sigma_true <- 0.8

X <- simulate_ou(
  time = grid,
  kappa = kappa_true,
  sigma = sigma_true,
  X0 = 0,
  stationary_start = TRUE
)

events <- simulate_covariate_hawkes(
  grid = grid,
  X = X,
  gamma0 = gamma0_true,
  gamma1 = gamma1_true,
  alpha = alpha_true,
  beta = beta_true
)

cat("Number of events:", length(events), "\n")

fit <- fit_covariate_hawkes(
  events = events,
  grid = grid,
  X = X,
  method = "Nelder-Mead"
)

theta_hat <- fit$theta_hat

gamma0_hat <- theta_hat["gamma0"]
gamma1_hat <- theta_hat["gamma1"]
alpha_hat <- exp(theta_hat["log_alpha"])
beta_hat <- exp(theta_hat["log_beta"])
branching_hat <- alpha_hat / beta_hat

cat("MLE:\n")
cat("gamma0_hat =", gamma0_hat, "\n")
cat("gamma1_hat =", gamma1_hat, "\n")
cat("alpha_hat =", alpha_hat, "\n")
cat("beta_hat =", beta_hat, "\n")
cat("branching_hat =", branching_hat, "\n")

z <- time_rescaling_residuals_covariate_hawkes(
  events = events,
  grid = grid,
  X = X,
  gamma0 = gamma0_hat,
  gamma1 = gamma1_hat,
  alpha = alpha_hat,
  beta = beta_hat
)

usual <- exponential_residual_summary(z)

corrected <- corrected_time_rescaling_gof(z)
subsampled <- corrected_time_rescaling_gof(
  z,
  m_subsample = floor(sqrt(length(z)))
)

fit_summary <- data.frame(
  n_events = length(events),
  T_end = T_end,
  gamma0_true = gamma0_true,
  gamma1_true = gamma1_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = branching_true,
  kappa_true = kappa_true,
  sigma_true = sigma_true,
  gamma0_hat = gamma0_hat,
  gamma1_hat = gamma1_hat,
  alpha_hat = alpha_hat,
  beta_hat = beta_hat,
  branching_hat = branching_hat,
  neg_loglik = fit$fit$value,
  convergence = fit$fit$convergence
)

gof_summary <- data.frame(
  n_events = length(events),
  usual_D = usual$ks_statistic,
  usual_p = usual$ks_p_value,
  corrected_D = corrected$corrected_D,
  corrected_p = corrected$corrected_p,
  subsampled_D = subsampled$subsampled_D,
  subsampled_p = subsampled$subsampled_p,
  mean_z = mean(z),
  var_z = var(z),
  median_z = median(z),
  n_nonpositive_corrected = corrected$n_nonpositive_corrected
)

print(fit_summary)
print(gof_summary)

write.csv(
  fit_summary,
  "results/SPBDBS2025_covariate_algorithm2_gof_fit_summary.csv",
  row.names = FALSE
)

write.csv(
  gof_summary,
  "results/SPBDBS2025_covariate_algorithm2_gof_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(z = z),
  "results/SPBDBS2025_covariate_algorithm2_gof_residuals.csv",
  row.names = FALSE
)

png("figures/SPBDBS2025_covariate_algorithm2_gof_residuals.png")
hist(z, breaks = 40, freq = FALSE, main = "Covariate-Hawkes residuals", xlab = "z")
curve(dexp(x), add = TRUE, lwd = 2)
dev.off()

pdf("figures/SPBDBS2025_covariate_algorithm2_gof_residuals.pdf")
hist(z, breaks = 40, freq = FALSE, main = "Covariate-Hawkes residuals", xlab = "z")
curve(dexp(x), add = TRUE, lwd = 2)
dev.off()

cat("\nSaved:\n")
cat("- results/SPBDBS2025_covariate_algorithm2_gof_fit_summary.csv\n")
cat("- results/SPBDBS2025_covariate_algorithm2_gof_summary.csv\n")
cat("- results/SPBDBS2025_covariate_algorithm2_gof_residuals.csv\n")
cat("- figures/SPBDBS2025_covariate_algorithm2_gof_residuals.pdf/png\n")