# 16_ou_covariate_hawkes_bootstrap_test.R
#
# OU covariate-adjusted bootstrap LR test.
#
# True model:
#   lambda(t) = exp(gamma0 + gamma1 X_t) + alpha R(t)
#
# Null model:
#   lambda(t) = exp(gamma0 + gamma1 X_t)
#
# Alternative:
#   lambda(t) = exp(gamma0 + gamma1 X_t) + alpha R(t)

set.seed(16)

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

B <- 100

# ------------------------------------------------------------
# 2. OU simulation
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
# 3. Simulation functions
# ------------------------------------------------------------

simulate_covariate_poisson <- function(grid, X, gamma0, gamma1) {
  T_end <- max(grid)
  
  lambda_grid <- exp(gamma0 + gamma1 * X)
  lambda_max <- max(lambda_grid)
  
  t <- 0
  events <- numeric(0)
  
  while (t < T_end) {
    t <- t + rexp(1, rate = lambda_max)
    
    if (t > T_end) {
      break
    }
    
    X_t <- approx(grid, X, xout = t, rule = 2)$y
    lambda_t <- exp(gamma0 + gamma1 * X_t)
    
    if (runif(1) <= lambda_t / lambda_max) {
      events <- c(events, t)
    }
  }
  
  events
}

simulate_covariate_hawkes <- function(
    grid,
    X,
    gamma0,
    gamma1,
    alpha,
    beta
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
    
    if (t_candidate > T_end) {
      break
    }
    
    R_candidate <- R * exp(-beta * (t_candidate - t))
    
    X_candidate <- approx(
      grid,
      X,
      xout = t_candidate,
      rule = 2
    )$y
    
    baseline_candidate <- exp(
      gamma0 + gamma1 * X_candidate
    )
    
    lambda_candidate <- baseline_candidate +
      alpha * R_candidate
    
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
# 4. Likelihoods
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
  
  if (alpha >= beta) {
    return(1e12)
  }
  
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
    
    if (
      !is.finite(lambda_events[i]) ||
      lambda_events[i] <= 0
    ) {
      return(1e12)
    }
  }
  
  baseline_grid <- exp(gamma0 + gamma1 * X)
  
  baseline_integral <- sum(baseline_grid) *
    (grid[2] - grid[1])
  
  hawkes_integral <- sum(
    alpha / beta *
      (1 - exp(-beta * (T_end - events)))
  )
  
  -(
    sum(log(lambda_events)) -
      baseline_integral -
      hawkes_integral
  )
}

# ------------------------------------------------------------
# 5. Fitting helpers
# ------------------------------------------------------------

fit_covariate_poisson <- function(events, grid, X) {
  init <- c(
    log(max(length(events), 1) / max(grid)),
    0
  )
  
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
  init <- c(
    fit_null$par[1],
    fit_null$par[2],
    log(0.2),
    log(1.0)
  )
  
  optim(
    par = init,
    fn = covariate_hawkes_negloglik,
    events = events,
    grid = grid,
    X = X,
    method = "Nelder-Mead",
    control = list(maxit = 10000)
  )
}

compute_LR <- function(events, grid, X) {
  fit_null <- fit_covariate_poisson(events, grid, X)
  fit_alt <- fit_covariate_hawkes(events, grid, X, fit_null)
  
  logL_null <- -fit_null$value
  logL_alt <- -fit_alt$value
  
  LR <- 2 * (logL_alt - logL_null)
  
  alpha_hat <- exp(fit_alt$par[3])
  beta_hat <- exp(fit_alt$par[4])
  
  list(
    LR = max(LR, 0),
    fit_null = fit_null,
    fit_alt = fit_alt,
    branching_hat = alpha_hat / beta_hat,
    alpha_hat = alpha_hat,
    beta_hat = beta_hat
  )
}

# ------------------------------------------------------------
# 6. Observed dataset
# ------------------------------------------------------------

events_obs <- simulate_covariate_hawkes(
  grid = grid,
  X = X,
  gamma0 = gamma0_true,
  gamma1 = gamma1_true,
  alpha = alpha_true,
  beta = beta_true
)

obs <- compute_LR(events_obs, grid, X)

cat("Observed number of events:", length(events_obs), "\n")
cat("Observed LR:", obs$LR, "\n")
cat("Observed branching estimate:", obs$branching_hat, "\n")

# ------------------------------------------------------------
# 7. Bootstrap under null
# ------------------------------------------------------------

gamma0_null_hat <- obs$fit_null$par[1]
gamma1_null_hat <- obs$fit_null$par[2]

boot_results <- vector("list", B)

for (b in seq_len(B)) {
  boot_events <- simulate_covariate_poisson(
    grid = grid,
    X = X,
    gamma0 = gamma0_null_hat,
    gamma1 = gamma1_null_hat
  )
  
  boot <- compute_LR(boot_events, grid, X)
  
  boot_results[[b]] <- data.frame(
    bootstrap = b,
    n_events = length(boot_events),
    LR = boot$LR,
    branching_hat = boot$branching_hat
  )
  
  cat(
    "Bootstrap", b, "of", B,
    "- LR:", round(boot$LR, 3),
    "- n:", length(boot_events),
    "\n"
  )
}

boot_df <- do.call(rbind, boot_results)

n_exceed <- sum(boot_df$LR >= obs$LR)

p_empirical <- n_exceed / B
p_corrected <- (1 + n_exceed) / (B + 1)

cat("\nEmpirical bootstrap p-value:", p_empirical, "\n")
cat("Corrected bootstrap p-value:", p_corrected, "\n")

# ------------------------------------------------------------
# 8. Save summaries
# ------------------------------------------------------------

summary_df <- data.frame(
  observed_n_events = length(events_obs),
  observed_LR = obs$LR,
  p_empirical = p_empirical,
  p_corrected = p_corrected,
  observed_branching_hat = obs$branching_hat,
  observed_alpha_hat = obs$alpha_hat,
  observed_beta_hat = obs$beta_hat,
  gamma0_true = gamma0_true,
  gamma1_true = gamma1_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = branching_true,
  kappa_true = kappa_true,
  sigma_true = sigma_true,
  gamma0_null_hat = gamma0_null_hat,
  gamma1_null_hat = gamma1_null_hat,
  B = B
)

print(summary_df)

write.csv(
  summary_df,
  "results/ou_covariate_hawkes_bootstrap_summary.csv",
  row.names = FALSE
)

write.csv(
  boot_df,
  "results/ou_covariate_hawkes_bootstrap_replicates.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Plots
# ------------------------------------------------------------

p_lr <- ggplot(boot_df, aes(x = LR)) +
  geom_histogram(bins = 30, boundary = 0) +
  geom_vline(
    xintercept = obs$LR,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(
    title = "OU covariate-adjusted Hawkes bootstrap test",
    subtitle =
      "Dashed line: observed LR; histogram: null bootstrap LR values",
    x = "likelihood-ratio statistic",
    y = "count"
  ) +
  theme_minimal()

ggsave(
  "figures/ou_covariate_hawkes_bootstrap_LR_histogram.pdf",
  p_lr,
  width = 7,
  height = 4
)

ggsave(
  "figures/ou_covariate_hawkes_bootstrap_LR_histogram.png",
  p_lr,
  width = 7,
  height = 4,
  dpi = 300
)

p_branching <- ggplot(boot_df, aes(x = branching_hat)) +
  geom_histogram(bins = 30, boundary = 0) +
  geom_vline(
    xintercept = obs$branching_hat,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(
    title = "Bootstrap branching estimates under OU covariate-only null",
    subtitle =
      "Dashed line: observed branching estimate",
    x = "estimated branching ratio",
    y = "count"
  ) +
  theme_minimal()

ggsave(
  "figures/ou_covariate_hawkes_bootstrap_branching_histogram.pdf",
  p_branching,
  width = 7,
  height = 4
)

ggsave(
  "figures/ou_covariate_hawkes_bootstrap_branching_histogram.png",
  p_branching,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/ou_covariate_hawkes_bootstrap_summary.csv\n")
cat("- results/ou_covariate_hawkes_bootstrap_replicates.csv\n")
cat("- figures/ou_covariate_hawkes_bootstrap_LR_histogram.pdf/png\n")
cat("- figures/ou_covariate_hawkes_bootstrap_branching_histogram.pdf/png\n")