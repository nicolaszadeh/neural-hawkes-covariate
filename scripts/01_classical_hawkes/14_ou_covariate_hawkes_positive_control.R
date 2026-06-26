# 14_ou_covariate_hawkes_positive_control.R
#
# OU covariate positive control.
#
# Covariate:
#   dX_t = -kappa X_t dt + sigma dW_t
#
# Event intensity:
#   lambda(t) = exp(gamma0 + gamma1 X_t) + alpha R(t)
#
# Compare:
#   1. covariate Poisson
#   2. covariate Hawkes

set.seed(14)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

source("R/load_all.R")

# ------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------

T_end <- 200
dt <- 0.01
grid <- seq(0, T_end, by = dt)

gamma0_true <- -0.5
gamma1_true <- 0.9

alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

kappa_true <- 0.4
sigma_true <- 0.8
X0 <- 0

# ------------------------------------------------------------
# 2. Simulate OU covariate
# ------------------------------------------------------------

simulate_ou <- function(grid, kappa, sigma, X0) {
  dt <- grid[2] - grid[1]
  X <- numeric(length(grid))
  X[1] <- X0
  
  for (i in 2:length(grid)) {
    X[i] <- X[i - 1] -
      kappa * X[i - 1] * dt +
      sigma * sqrt(dt) * rnorm(1)
  }
  
  X
}

X <- simulate_ou(
  grid = grid,
  kappa = kappa_true,
  sigma = sigma_true,
  X0 = X0
)

# ------------------------------------------------------------
# 3. Simulate covariate-Hawkes process
# ------------------------------------------------------------

simulate_covariate_hawkes <- function(
    grid, X, gamma0, gamma1, alpha, beta
) {
  T_end <- max(grid)
  
  baseline_grid <- exp(gamma0 + gamma1 * X)
  baseline_max <- max(baseline_grid)
  
  t <- 0
  events <- numeric(0)
  R <- 0
  
  while (t < T_end) {
    lambda_upper <- baseline_max + alpha * R
    
    if (!is.finite(lambda_upper) || lambda_upper <= 0) {
      stop("Invalid upper intensity.")
    }
    
    t_candidate <- t + rexp(1, rate = lambda_upper)
    
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
  beta = beta_true
)

cat("Number of events:", length(events), "\n")

# ------------------------------------------------------------
# 4. Likelihoods
# ------------------------------------------------------------

covariate_poisson_negloglik <- function(theta, events, grid, X) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]
  
  lambda_grid <- exp(gamma0 + gamma1 * X)
  integral <- sum(lambda_grid) * (grid[2] - grid[1])
  
  X_events <- approx(grid, X, xout = events,
                     rule = 2)$y
  
  log_sum <- sum(gamma0 + gamma1 * X_events)
  
  -(log_sum - integral)
}

covariate_hawkes_negloglik <- function(theta, events, grid, X) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]
  alpha <- exp(theta[3])
  beta <- exp(theta[4])
  
  if (alpha >= beta) return(1e12)
  
  T_end <- max(grid)
  n <- length(events)
  
  if (n == 0) {
    lambda_grid <- exp(gamma0 + gamma1 * X)
    return(sum(lambda_grid) * (grid[2] - grid[1]))
  }
  
  X_events <- approx(grid, X, xout = events,
                     rule = 2)$y
  baseline_events <- exp(gamma0 + gamma1 * X_events)
  
  R <- numeric(n)
  lambda_events <- numeric(n)
  
  for (i in seq_len(n)) {
    if (i == 1) {
      R[i] <- 0
    } else {
      delta <- events[i] - events[i - 1]
      R[i] <- exp(-beta * delta) * (1 + R[i - 1])
    }
    
    lambda_events[i] <- baseline_events[i] + alpha * R[i]
    
    if (!is.finite(lambda_events[i]) ||
        lambda_events[i] <= 0) {
      return(1e12)
    }
  }
  
  baseline_grid <- exp(gamma0 + gamma1 * X)
  baseline_integral <- sum(baseline_grid) * (grid[2] - grid[1])
  
  hawkes_integral <- sum(
    alpha / beta * (1 - exp(-beta * (T_end - events)))
  )
  
  -(sum(log(lambda_events)) -
      baseline_integral -
      hawkes_integral)
}

# ------------------------------------------------------------
# 5. Fit models
# ------------------------------------------------------------

fit_cov <- optim(
  par = c(log(max(length(events), 1) / T_end), 0),
  fn = covariate_poisson_negloglik,
  events = events,
  grid = grid,
  X = X,
  method = "BFGS"
)

fit_hawkes <- optim(
  par = c(fit_cov$par[1], fit_cov$par[2],
          log(0.2), log(1.0)),
  fn = covariate_hawkes_negloglik,
  events = events,
  grid = grid,
  X = X,
  method = "Nelder-Mead",
  control = list(maxit = 10000)
)

# ------------------------------------------------------------
# 6. Summaries
# ------------------------------------------------------------

gamma0_cov_hat <- fit_cov$par[1]
gamma1_cov_hat <- fit_cov$par[2]

gamma0_hawkes_hat <- fit_hawkes$par[1]
gamma1_hawkes_hat <- fit_hawkes$par[2]
alpha_hat <- exp(fit_hawkes$par[3])
beta_hat <- exp(fit_hawkes$par[4])
branching_hat <- alpha_hat / beta_hat

model_comparison <- data.frame(
  model = c("ou_covariate_poisson", "ou_covariate_hawkes"),
  neg_loglik = c(fit_cov$value, fit_hawkes$value),
  AIC = c(
    2 * 2 + 2 * fit_cov$value,
    2 * 4 + 2 * fit_hawkes$value
  )
)

estimates <- data.frame(
  parameter = c(
    "gamma0_true",
    "gamma1_true",
    "alpha_true",
    "beta_true",
    "branching_true",
    "kappa_true",
    "sigma_true",
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
    kappa_true,
    sigma_true,
    gamma0_cov_hat,
    gamma1_cov_hat,
    gamma0_hawkes_hat,
    gamma1_hawkes_hat,
    alpha_hat,
    beta_hat,
    branching_hat
  )
)

print(model_comparison)
print(estimates)

write.csv(
  model_comparison,
  "results/ou_covariate_hawkes_model_comparison.csv",
  row.names = FALSE
)

write.csv(
  estimates,
  "results/ou_covariate_hawkes_estimates.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Reconstruct fitted intensity
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

baseline_true <- exp(gamma0_true + gamma1_true * X)
R_true_grid <- compute_hawkes_memory_grid(grid, events, beta_true)
lambda_true <- baseline_true + alpha_true * R_true_grid

R_hat_grid <- compute_hawkes_memory_grid(grid, events, beta_hat)
lambda_hat <- exp(gamma0_hawkes_hat + gamma1_hawkes_hat * X) +
  alpha_hat * R_hat_grid

plot_df <- data.frame(
  t = grid,
  X = X,
  baseline_true = baseline_true,
  lambda_true = lambda_true,
  lambda_hat = lambda_hat
)

event_df <- data.frame(t = events)

p_intensity <- ggplot(plot_df, aes(x = t)) +
  geom_line(aes(y = lambda_true), linewidth = 0.6) +
  geom_line(aes(y = lambda_hat), linetype = "dashed",
            linewidth = 0.6) +
  geom_line(aes(y = baseline_true), linetype = "dotted",
            linewidth = 0.5) +
  geom_rug(
    data = event_df,
    aes(x = t),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.25
  ) +
  labs(
    title = "OU covariate Hawkes positive control",
    subtitle = "Solid: true intensity; dashed: fitted Hawkes; dotted: true baseline",
    x = "time",
    y = "intensity"
  ) +
  theme_minimal()

ggsave(
  "figures/ou_covariate_hawkes_intensity.pdf",
  p_intensity,
  width = 8,
  height = 4
)

ggsave(
  "figures/ou_covariate_hawkes_intensity.png",
  p_intensity,
  width = 8,
  height = 4,
  dpi = 300
)

p_covariate <- ggplot(plot_df, aes(x = t, y = X)) +
  geom_line(linewidth = 0.6) +
  labs(
    title = "Simulated Ornstein-Uhlenbeck covariate",
    x = "time",
    y = "X(t)"
  ) +
  theme_minimal()

ggsave(
  "figures/ou_covariate_path.pdf",
  p_covariate,
  width = 8,
  height = 4
)

ggsave(
  "figures/ou_covariate_path.png",
  p_covariate,
  width = 8,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/ou_covariate_hawkes_model_comparison.csv\n")
cat("- results/ou_covariate_hawkes_estimates.csv\n")
cat("- figures/ou_covariate_hawkes_intensity.pdf/png\n")
cat("- figures/ou_covariate_path.pdf/png\n")