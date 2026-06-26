# Smoke test for the refactored helper library.
# Run from the repository root:
#   source("TEST_refactor_smoke_test.R")

set.seed(123)

source("R/load_all.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

T_end <- 100
mu <- 0.8
alpha <- 0.35
beta <- 1.2

events <- simulate_hawkes(T_end, mu, alpha, beta)
cat("Simulated events:", length(events), "\n")

fit <- fit_hawkes(events, T_end)
cat("Fitted branching:", fit$branching, "\n")

info <- hawkes_observed_information(
  events = events,
  T_end = T_end,
  theta_hat = fit$theta_hat
)

wald <- wald_all_true_parameters(
  theta_hat = fit$theta_hat,
  I_inv_hat = info$I_inv_hat,
  theta_true = log(c(
    log_mu = mu,
    log_alpha = alpha,
    log_beta = beta
  )),
  T_end = T_end
)
print(wald)

z <- hawkes_rescaled_gaps(
  events = events,
  mu = fit$mu,
  alpha = fit$alpha,
  beta = fit$beta
)
print(ks_exp1_summary(z))

a2 <- hawkes_algorithm2_test(
  events = events,
  T_end = T_end,
  mu = fit$mu,
  alpha = fit$alpha,
  beta = fit$beta
)
print(a2)

# Covariate-Hawkes smoke test.
grid <- seq(0, 100, by = 0.01)
X <- simulate_ou(grid, kappa = 0.4, sigma = 0.8, X0 = 0)
cevents <- simulate_covariate_hawkes(
  grid = grid,
  X = X,
  gamma0 = -0.5,
  gamma1 = 0.9,
  alpha = 0.35,
  beta = 1.2
)
cat("Covariate-Hawkes events:", length(cevents), "\n")

cfit <- fit_covariate_hawkes(
  events = cevents,
  grid = grid,
  X = X
)
print(cfit$theta_hat)

cat("\nSmoke test completed.\n")
