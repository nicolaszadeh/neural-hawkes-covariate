# Simulate events from:
#   lambda(t) = exp(gamma0 + gamma1 X(t)) + alpha R(t)
#
# where R(t) is the Hawkes memory.
#
# Then test whether adding Hawkes excitation improves over
# the covariate-only model.
#
# Scientific point:
# excitation should remain detectable after accounting for X(t).
# 07_covariate_hawkes_positive_control.R
#
# Positive control:
# events are generated from
#
#   lambda(t) = exp(gamma0 + gamma1 X(t)) + alpha R(t),
#
# where R(t) is the Hawkes memory:
#
#   R(t) = sum_{t_k < t} exp(-beta * (t - t_k)).
#
# We compare:
#   1. covariate Poisson
#   2. covariate Hawkes
#
# The point is to check whether genuine excitation remains
# detectable after accounting for the covariate X(t).

set.seed(7)

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

gamma0_true <- -0.5
gamma1_true <- 0.9

alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

baseline_grid <- exp(gamma0_true + gamma1_true * X)

# ------------------------------------------------------------
# 2. Simulate covariate Hawkes process by thinning
# ------------------------------------------------------------

simulate_covariate_hawkes <- function(
    grid, X, gamma0, gamma1, alpha, beta, T_end
) {
  t <- 0
  events <- numeric(0)
  R <- 0
  
  baseline_grid <- exp(gamma0 + gamma1 * X)
  baseline_max <- max(baseline_grid)
  
  while (t < T_end) {
    lambda_upper <- baseline_max + alpha * R
    
    if (!is.finite(lambda_upper) || lambda_upper <= 0) {
      stop("Invalid upper intensity.")
    }
    
    waiting_time <- rexp(1, rate = lambda_upper)
    t_candidate <- t + waiting_time
    
    if (t_candidate > T_end) break
    
    R_candidate <- R * exp(-beta * (t_candidate - t))
    
    X_candidate <- approx(grid, X, xout = t_candidate,
                          rule = 2)$y
    
    baseline_candidate <- exp(gamma0 + gamma1 * X_candidate)
    lambda_candidate <- baseline_candidate + alpha * R_candidate
    
    if (runif(1) <= lambda_candidate / lambda_upper) {
      events <- c(events, t_candidate)
      R <- R_candidate + 1
    } else {
      R <- R_candidate
    }
    
    t <- t_candidate
  }
  
  events
}

events <- simulate_covariate_hawkes(
  grid = grid,
  X = X,
  gamma0 = gamma0_true,
  gamma1 = gamma1_true,
  alpha = alpha_true,
  beta = beta_true,
  T_end = T_end
)

n_events <- length(events)
cat("Number of events:", n_events, "\n")

# ------------------------------------------------------------
# 3. Likelihoods
# ------------------------------------------------------------

poisson_cov_negloglik <- function(theta, events, grid, X) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]
  
  lambda_grid <- exp(gamma0 + gamma1 * X)
  integral <- sum(lambda_grid) * (grid[2] - grid[1])
  
  X_events <- approx(grid, X, xout = events, rule = 2)$y
  log_intensity_sum <- sum(gamma0 + gamma1 * X_events)
  
  -(log_intensity_sum - integral)
}

covariate_hawkes_negloglik <- function(theta, events, grid, X, T_end) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]
  alpha <- exp(theta[3])
  beta <- exp(theta[4])
  
  if (alpha >= beta) {
    return(1e12)
  }
  
  n <- length(events)
  if (n == 0) {
    lambda_grid <- exp(gamma0 + gamma1 * X)
    return(sum(lambda_grid) * (grid[2] - grid[1]))
  }
  
  R <- numeric(n)
  lambda_events <- numeric(n)
  
  X_events <- approx(grid, X, xout = events, rule = 2)$y
  baseline_events <- exp(gamma0 + gamma1 * X_events)
  
  for (i in seq_len(n)) {
    if (i == 1) {
      R[i] <- 0
    } else {
      dt_i <- events[i] - events[i - 1]
      R[i] <- exp(-beta * dt_i) * (1 + R[i - 1])
    }
    
    lambda_events[i] <- baseline_events[i] + alpha * R[i]
    
    if (!is.finite(lambda_events[i]) || lambda_events[i] <= 0) {
      return(1e12)
    }
  }
  
  baseline_grid <- exp(gamma0 + gamma1 * X)
  baseline_integral <- sum(baseline_grid) * (grid[2] - grid[1])
  
  hawkes_integral <- sum(
    alpha / beta * (1 - exp(-beta * (T_end - events)))
  )
  
  integral <- baseline_integral + hawkes_integral
  
  -(sum(log(lambda_events)) - integral)
}

# ------------------------------------------------------------
# 4. Fit models
# ------------------------------------------------------------

fit_cov <- optim(
  par = c(log(n_events / T_end), 0),
  fn = poisson_cov_negloglik,
  events = events,
  grid = grid,
  X = X,
  method = "BFGS"
)

fit_cov_hawkes <- optim(
  par = c(fit_cov$par[1], fit_cov$par[2], log(0.2), log(1.0)),
  fn = covariate_hawkes_negloglik,
  events = events,
  grid = grid,
  X = X,
  T_end = T_end,
  method = "Nelder-Mead",
  control = list(maxit = 8000)
)

# ------------------------------------------------------------
# 5. Summaries
# ------------------------------------------------------------

gamma0_cov_hat <- fit_cov$par[1]
gamma1_cov_hat <- fit_cov$par[2]

gamma0_hawkes_hat <- fit_cov_hawkes$par[1]
gamma1_hawkes_hat <- fit_cov_hawkes$par[2]
alpha_hat <- exp(fit_cov_hawkes$par[3])
beta_hat <- exp(fit_cov_hawkes$par[4])
branching_hat <- alpha_hat / beta_hat

summary <- data.frame(
  model = c("covariate_poisson", "covariate_hawkes"),
  neg_loglik = c(fit_cov$value, fit_cov_hawkes$value),
  AIC = c(
    2 * 2 + 2 * fit_cov$value,
    2 * 4 + 2 * fit_cov_hawkes$value
  )
)

estimates <- data.frame(
  parameter = c(
    "gamma0_true",
    "gamma1_true",
    "alpha_true",
    "beta_true",
    "branching_true",
    "gamma0_cov_hat",
    "gamma1_cov_hat",
    "gamma0_hawkes_hat",
    "gamma1_hawkes_hat",
    "alpha_hat",
    "beta_hat",
    "branching_hat"
  ),
  value = c(
    gamma0_true,
    gamma1_true,
    alpha_true,
    beta_true,
    branching_true,
    gamma0_cov_hat,
    gamma1_cov_hat,
    gamma0_hawkes_hat,
    gamma1_hawkes_hat,
    alpha_hat,
    beta_hat,
    branching_hat
  )
)

print(summary)
print(estimates)

write.csv(
  summary,
  "results/covariate_hawkes_positive_model_comparison.csv",
  row.names = FALSE
)

write.csv(
  estimates,
  "results/covariate_hawkes_positive_estimates.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Reconstruct fitted intensity on grid
# ------------------------------------------------------------

compute_hawkes_memory_grid <- function(grid, events, beta) {
  R_grid <- numeric(length(grid))
  event_index <- 1
  R <- 0
  previous_t <- grid[1]
  
  for (i in seq_along(grid)) {
    t <- grid[i]
    R <- R * exp(-beta * (t - previous_t))
    
    while (event_index <= length(events) &&
           events[event_index] <= t) {
      R <- R + exp(-beta * (t - events[event_index]))
      event_index <- event_index + 1
    }
    
    R_grid[i] <- R
    previous_t <- t
  }
  
  R_grid
}

R_true_grid <- compute_hawkes_memory_grid(grid, events, beta_true)
lambda_true_grid <- baseline_grid + alpha_true * R_true_grid

R_hat_grid <- compute_hawkes_memory_grid(grid, events, beta_hat)
lambda_hat_grid <- exp(gamma0_hawkes_hat + gamma1_hawkes_hat * X) +
  alpha_hat * R_hat_grid

plot_df <- data.frame(
  t = grid,
  lambda_true = lambda_true_grid,
  lambda_hat = lambda_hat_grid,
  baseline_true = baseline_grid
)

event_df <- data.frame(t = events)

p1 <- ggplot(plot_df, aes(x = t)) +
  geom_line(aes(y = lambda_true), linewidth = 0.6) +
  geom_line(aes(y = lambda_hat), linetype = "dashed",
            linewidth = 0.6) +
  geom_line(aes(y = baseline_true), linetype = "dotted",
            linewidth = 0.5) +
  geom_rug(data = event_df, aes(x = t), inherit.aes = FALSE,
           sides = "b", alpha = 0.25) +
  labs(
    title = "Covariate Hawkes positive control",
    subtitle = "Solid: true intensity; dashed: fitted covariate Hawkes; dotted: true baseline",
    x = "time",
    y = "intensity"
  ) +
  theme_minimal()

ggsave(
  "figures/covariate_hawkes_positive_intensity.pdf",
  p1,
  width = 8,
  height = 4
)

ggsave(
  "figures/covariate_hawkes_positive_intensity.png",
  p1,
  width = 8,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/covariate_hawkes_positive_model_comparison.csv\n")
cat("- results/covariate_hawkes_positive_estimates.csv\n")
cat("- figures/covariate_hawkes_positive_intensity.pdf/png\n")