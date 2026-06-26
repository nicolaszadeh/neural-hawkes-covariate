# Simulate events from:
#   lambda(t) = exp(gamma0 + gamma1 X(t))
# with no Hawkes excitation.
#
# Then compare:
#   1. constant Poisson
#   2. covariate Poisson
#   3. naive Hawkes without covariate
#
# Scientific point:
# clustering can come from X(t), not only from self-excitation.
# 06_covariate_baseline_control.R
#
# Negative control:
# events are generated from a covariate-driven Poisson process
#
#   lambda(t) = exp(gamma0 + gamma1 X(t)),
#
# with no Hawkes excitation.
#
# We compare:
#   1. constant Poisson
#   2. covariate Poisson
#   3. naive Hawkes without covariate
#
# The point is to show that apparent clustering may come from
# the covariate X(t), not from self-excitation.

set.seed(6)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

source("R/load_all.R")

# ------------------------------------------------------------
# 1. Time grid and covariate
# ------------------------------------------------------------

T_end <- 200
dt <- 0.01
grid <- seq(0, T_end, by = dt)

X <- sin(2 * pi * grid / 40) + 0.5 * sin(2 * pi * grid / 13)

gamma0_true <- -0.3
gamma1_true <- 1.1

lambda_true <- exp(gamma0_true + gamma1_true * X)

# ------------------------------------------------------------
# 2. Simulate inhomogeneous Poisson process by thinning
# ------------------------------------------------------------

simulate_inhom_poisson <- function(grid, lambda_values) {
  T_end <- max(grid)
  lambda_max <- max(lambda_values)
  
  t <- 0
  events <- numeric(0)
  
  while (t < T_end) {
    t <- t + rexp(1, rate = lambda_max)
    if (t > T_end) break
    
    lambda_t <- approx(grid, lambda_values, xout = t,
                       rule = 2)$y
    
    if (runif(1) <= lambda_t / lambda_max) {
      events <- c(events, t)
    }
  }
  
  events
}

events <- simulate_inhom_poisson(grid, lambda_true)
n_events <- length(events)

cat("Number of events:", n_events, "\n")

# ------------------------------------------------------------
# 3. Log-likelihoods
# ------------------------------------------------------------

poisson_const_negloglik <- function(theta, events, T_end) {
  mu <- exp(theta[1])
  -(length(events) * log(mu) - mu * T_end)
}

poisson_cov_negloglik <- function(theta, events, grid, X, T_end) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]
  
  lambda_grid <- exp(gamma0 + gamma1 * X)
  integral <- sum(lambda_grid) * (grid[2] - grid[1])
  
  X_events <- approx(grid, X, xout = events, rule = 2)$y
  log_intensity_sum <- sum(gamma0 + gamma1 * X_events)
  
  -(log_intensity_sum - integral)
}

hawkes_naive_negloglik <- function(theta, events, T_end) {
  mu <- exp(theta[1])
  alpha <- exp(theta[2])
  beta <- exp(theta[3])
  
  if (alpha >= beta) {
    return(1e12)
  }
  
  n <- length(events)
  if (n == 0) {
    return(mu * T_end)
  }
  
  R <- numeric(n)
  lambda <- numeric(n)
  
  for (i in seq_len(n)) {
    if (i == 1) {
      R[i] <- 0
    } else {
      dt_i <- events[i] - events[i - 1]
      R[i] <- exp(-beta * dt_i) * (1 + R[i - 1])
    }
    
    lambda[i] <- mu + alpha * R[i]
    
    if (!is.finite(lambda[i]) || lambda[i] <= 0) {
      return(1e12)
    }
  }
  
  integral <- mu * T_end +
    sum(alpha / beta * (1 - exp(-beta * (T_end - events))))
  
  -(sum(log(lambda)) - integral)
}

# ------------------------------------------------------------
# 4. Fit models
# ------------------------------------------------------------

fit_const <- optim(
  par = log(n_events / T_end),
  fn = poisson_const_negloglik,
  events = events,
  T_end = T_end,
  method = "BFGS"
)

fit_cov <- optim(
  par = c(log(n_events / T_end), 0),
  fn = poisson_cov_negloglik,
  events = events,
  grid = grid,
  X = X,
  T_end = T_end,
  method = "BFGS"
)

fit_hawkes <- optim(
  par = c(log(n_events / T_end), log(0.2), log(2.0)),
  fn = hawkes_naive_negloglik,
  events = events,
  T_end = T_end,
  method = "Nelder-Mead",
  control = list(maxit = 5000)
)

# ------------------------------------------------------------
# 5. Summaries
# ------------------------------------------------------------

mu_const_hat <- exp(fit_const$par[1])

gamma0_hat <- fit_cov$par[1]
gamma1_hat <- fit_cov$par[2]

mu_hawkes_hat <- exp(fit_hawkes$par[1])
alpha_hawkes_hat <- exp(fit_hawkes$par[2])
beta_hawkes_hat <- exp(fit_hawkes$par[3])
branching_hat <- alpha_hawkes_hat / beta_hawkes_hat

summary <- data.frame(
  model = c("constant_poisson",
            "covariate_poisson",
            "naive_hawkes"),
  neg_loglik = c(fit_const$value,
                 fit_cov$value,
                 fit_hawkes$value),
  AIC = c(
    2 * 1 + 2 * fit_const$value,
    2 * 2 + 2 * fit_cov$value,
    2 * 3 + 2 * fit_hawkes$value
  )
)

estimates <- data.frame(
  parameter = c(
    "gamma0_true",
    "gamma1_true",
    "gamma0_hat",
    "gamma1_hat",
    "mu_const_hat",
    "mu_hawkes_hat",
    "alpha_hawkes_hat",
    "beta_hawkes_hat",
    "branching_hat"
  ),
  value = c(
    gamma0_true,
    gamma1_true,
    gamma0_hat,
    gamma1_hat,
    mu_const_hat,
    mu_hawkes_hat,
    alpha_hawkes_hat,
    beta_hawkes_hat,
    branching_hat
  )
)

print(summary)
print(estimates)

write.csv(
  summary,
  "results/covariate_baseline_model_comparison.csv",
  row.names = FALSE
)

write.csv(
  estimates,
  "results/covariate_baseline_estimates.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Diagnostic plot
# ------------------------------------------------------------

event_df <- data.frame(t = events)

plot_df <- data.frame(
  t = grid,
  X = X,
  lambda_true = lambda_true
)

p1 <- ggplot(plot_df, aes(x = t)) +
  geom_line(aes(y = lambda_true)) +
  geom_rug(data = event_df, aes(x = t), inherit.aes = FALSE,
           sides = "b", alpha = 0.35) +
  labs(
    title = "Covariate-driven Poisson process",
    subtitle = "Events cluster because lambda(t) depends on X(t), not because of Hawkes excitation",
    x = "time",
    y = "true intensity"
  ) +
  theme_minimal()

ggsave(
  "figures/covariate_baseline_intensity.pdf",
  p1,
  width = 8,
  height = 4
)

ggsave(
  "figures/covariate_baseline_intensity.png",
  p1,
  width = 8,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/covariate_baseline_model_comparison.csv\n")
cat("- results/covariate_baseline_estimates.csv\n")
cat("- figures/covariate_baseline_intensity.pdf/png\n")