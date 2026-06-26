# 08_covariate_hawkes_bootstrap_test.R
#
# Bootstrap likelihood-ratio test:
#
# H0: lambda(t) = exp(gamma0 + gamma1 X(t))
# H1: lambda(t) = exp(gamma0 + gamma1 X(t)) + alpha R(t)
#
# Question:
# after accounting for X(t), is there evidence of Hawkes
# self-excitation?

set.seed(8)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

source("R/load_all.R")

# ------------------------------------------------------------
# 1. Grid, covariate, true positive-control model
# ------------------------------------------------------------

T_end <- 200
dt <- 0.01
grid <- seq(0, T_end, by = dt)

X <- sin(2 * pi * grid / 40) + 0.5 * sin(2 * pi * grid / 13)

gamma0_true <- -0.5
gamma1_true <- 0.9
alpha_true <- 0.35
beta_true <- 1.2

# ------------------------------------------------------------
# 2. Simulation functions
# ------------------------------------------------------------

simulate_covariate_poisson <- function(grid, X, gamma0, gamma1) {
  T_end <- max(grid)
  lambda_grid <- exp(gamma0 + gamma1 * X)
  lambda_max <- max(lambda_grid)
  
  t <- 0
  events <- numeric(0)
  
  while (t < T_end) {
    t <- t + rexp(1, rate = lambda_max)
    if (t > T_end) break
    
    X_t <- approx(grid, X, xout = t, rule = 2)$y
    lambda_t <- exp(gamma0 + gamma1 * X_t)
    
    if (runif(1) <= lambda_t / lambda_max) {
      events <- c(events, t)
    }
  }
  
  events
}

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

# ------------------------------------------------------------
# 3. Likelihoods
# ------------------------------------------------------------

covariate_poisson_negloglik <- function(theta, events, grid, X) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]
  
  lambda_grid <- exp(gamma0 + gamma1 * X)
  integral <- sum(lambda_grid) * (grid[2] - grid[1])
  
  X_events <- approx(grid, X, xout = events, rule = 2)$y
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
  
  X_events <- approx(grid, X, xout = events, rule = 2)$y
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
    
    if (!is.finite(lambda_events[i]) || lambda_events[i] <= 0) {
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
# 4. Fitting helpers
# ------------------------------------------------------------

fit_covariate_poisson <- function(events, grid, X) {
  init <- c(log(length(events) / max(grid)), 0)
  
  optim(
    par = init,
    fn = covariate_poisson_negloglik,
    events = events,
    grid = grid,
    X = X,
    method = "BFGS"
  )
}

fit_covariate_hawkes <- function(events, grid, X, fit_null) {
  init <- c(fit_null$par[1], fit_null$par[2], log(0.2), log(1.0))
  
  optim(
    par = init,
    fn = covariate_hawkes_negloglik,
    events = events,
    grid = grid,
    X = X,
    method = "Nelder-Mead",
    control = list(maxit = 8000)
  )
}

compute_LR <- function(events, grid, X) {
  fit_null <- fit_covariate_poisson(events, grid, X)
  fit_alt <- fit_covariate_hawkes(events, grid, X, fit_null)
  
  LR <- 2 * ( -fit_alt$value - (-fit_null$value) )
  
  alpha_hat <- exp(fit_alt$par[3])
  beta_hat <- exp(fit_alt$par[4])
  
  list(
    LR = LR,
    fit_null = fit_null,
    fit_alt = fit_alt,
    branching_hat = alpha_hat / beta_hat
  )
}

# ------------------------------------------------------------
# 5. Observed positive-control dataset
# ------------------------------------------------------------

events_obs <- simulate_covariate_hawkes(
  grid = grid,
  X = X,
  gamma0 = gamma0_true,
  gamma1 = gamma1_true,
  alpha = alpha_true,
  beta = beta_true
)

cat("Observed number of events:", length(events_obs), "\n")

obs <- compute_LR(events_obs, grid, X)

cat("Observed LR:", obs$LR, "\n")
cat("Observed branching estimate:", obs$branching_hat, "\n")

# ------------------------------------------------------------
# 6. Parametric bootstrap under H0
# ------------------------------------------------------------

B <- 100

gamma0_null_hat <- obs$fit_null$par[1]
gamma1_null_hat <- obs$fit_null$par[2]

bootstrap_LR <- numeric(B)
bootstrap_branching <- numeric(B)
bootstrap_n_events <- numeric(B)

for (b in seq_len(B)) {
  boot_events <- simulate_covariate_poisson(
    grid = grid,
    X = X,
    gamma0 = gamma0_null_hat,
    gamma1 = gamma1_null_hat
  )
  
  boot_fit <- compute_LR(boot_events, grid, X)
  
  bootstrap_LR[b] <- boot_fit$LR
  bootstrap_branching[b] <- boot_fit$branching_hat
  bootstrap_n_events[b] <- length(boot_events)
  
  cat("Bootstrap", b, "of", B,
      "- LR:", round(bootstrap_LR[b], 3),
      "- n:", bootstrap_n_events[b], "\n")
}

p_value <- mean(bootstrap_LR >= obs$LR)

cat("\nBootstrap p-value:", p_value, "\n")

# ------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------

summary <- data.frame(
  observed_LR = obs$LR,
  p_value = p_value,
  observed_branching_hat = obs$branching_hat,
  gamma0_true = gamma0_true,
  gamma1_true = gamma1_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = alpha_true / beta_true,
  gamma0_null_hat = gamma0_null_hat,
  gamma1_null_hat = gamma1_null_hat,
  B = B
)

replicates <- data.frame(
  replicate = seq_len(B),
  LR = bootstrap_LR,
  branching_hat = bootstrap_branching,
  n_events = bootstrap_n_events
)

print(summary)

write.csv(
  summary,
  "results/covariate_hawkes_bootstrap_summary.csv",
  row.names = FALSE
)

write.csv(
  replicates,
  "results/covariate_hawkes_bootstrap_replicates.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Plot null LR distribution
# ------------------------------------------------------------

plot_df <- data.frame(LR = bootstrap_LR)

p <- ggplot(plot_df, aes(x = LR)) +
  geom_histogram(bins = 30, boundary = 0) +
  geom_vline(xintercept = obs$LR, linetype = "dashed",
             linewidth = 0.8) +
  labs(
    title = "Bootstrap LR test: covariate-only vs covariate-Hawkes",
    subtitle = paste0("Dashed line: observed LR; p-value = ",
                      signif(p_value, 3)),
    x = "bootstrap LR under H0",
    y = "count"
  ) +
  theme_minimal()

ggsave(
  "figures/covariate_hawkes_bootstrap_LR_histogram.pdf",
  p,
  width = 7,
  height = 4
)

ggsave(
  "figures/covariate_hawkes_bootstrap_LR_histogram.png",
  p,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/covariate_hawkes_bootstrap_summary.csv\n")
cat("- results/covariate_hawkes_bootstrap_replicates.csv\n")
cat("- figures/covariate_hawkes_bootstrap_LR_histogram.pdf/png\n")