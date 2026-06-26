# 10_covariate_hawkes_monte_carlo_recovery.R
#
# Monte Carlo recovery study for the correctly specified model:
#
#   lambda(t) = exp(gamma0 + gamma1 X(t)) + alpha R(t)
#
# Question:
# can we recover the covariate effect and the Hawkes effect
# simultaneously?

set.seed(10)

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

X <- sin(2 * pi * grid / 40) + 0.5 * sin(2 * pi * grid / 13)

gamma0_true <- -0.5
gamma1_true <- 0.9
alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

n_rep <- 50

# ------------------------------------------------------------
# 2. Simulation
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
# 3. Likelihood
# ------------------------------------------------------------

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

fit_covariate_hawkes <- function(events, grid, X) {
  n_events <- length(events)
  
  init <- c(
    log(n_events / max(grid)),
    0,
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

# ------------------------------------------------------------
# 4. Monte Carlo loop
# ------------------------------------------------------------

results_list <- list()

for (rep in seq_len(n_rep)) {
  events <- simulate_covariate_hawkes(
    grid = grid,
    X = X,
    gamma0 = gamma0_true,
    gamma1 = gamma1_true,
    alpha = alpha_true,
    beta = beta_true
  )
  
  fit <- fit_covariate_hawkes(events, grid, X)
  
  gamma0_hat <- fit$par[1]
  gamma1_hat <- fit$par[2]
  alpha_hat <- exp(fit$par[3])
  beta_hat <- exp(fit$par[4])
  branching_hat <- alpha_hat / beta_hat
  
  results_list[[rep]] <- data.frame(
    replicate = rep,
    n_events = length(events),
    gamma0_hat = gamma0_hat,
    gamma1_hat = gamma1_hat,
    alpha_hat = alpha_hat,
    beta_hat = beta_hat,
    branching_hat = branching_hat,
    neg_loglik = fit$value,
    convergence = fit$convergence
  )
  
  cat(
    "rep", rep, "of", n_rep,
    "- n =", length(events),
    "- gamma0 =", round(gamma0_hat, 3),
    "- gamma1 =", round(gamma1_hat, 3),
    "- branching =", round(branching_hat, 3),
    "\n"
  )
}

results <- do.call(rbind, results_list)

# ------------------------------------------------------------
# 5. Long-format estimates
# ------------------------------------------------------------

estimates_long <- rbind(
  data.frame(
    replicate = results$replicate,
    parameter = "gamma0",
    estimate = results$gamma0_hat,
    truth = gamma0_true
  ),
  data.frame(
    replicate = results$replicate,
    parameter = "gamma1",
    estimate = results$gamma1_hat,
    truth = gamma1_true
  ),
  data.frame(
    replicate = results$replicate,
    parameter = "alpha",
    estimate = results$alpha_hat,
    truth = alpha_true
  ),
  data.frame(
    replicate = results$replicate,
    parameter = "beta",
    estimate = results$beta_hat,
    truth = beta_true
  ),
  data.frame(
    replicate = results$replicate,
    parameter = "branching",
    estimate = results$branching_hat,
    truth = branching_true
  )
)

summary_list <- list()

for (param in unique(estimates_long$parameter)) {
  x <- estimates_long$estimate[estimates_long$parameter == param]
  truth <- unique(estimates_long$truth[estimates_long$parameter == param])
  
  summary_list[[param]] <- data.frame(
    parameter = param,
    truth = truth,
    mean_estimate = mean(x),
    bias = mean(x) - truth,
    sd_estimate = sd(x),
    q05 = as.numeric(quantile(x, 0.05)),
    median = as.numeric(quantile(x, 0.50)),
    q95 = as.numeric(quantile(x, 0.95))
  )
}

summary <- do.call(rbind, summary_list)
rownames(summary) <- NULL

print(summary)

write.csv(
  results,
  "results/covariate_hawkes_recovery_replicates.csv",
  row.names = FALSE
)

write.csv(
  summary,
  "results/covariate_hawkes_recovery_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Plots
# ------------------------------------------------------------

p <- ggplot(estimates_long, aes(x = parameter, y = estimate)) +
  geom_boxplot(outlier.alpha = 0.35) +
  geom_point(
    aes(y = truth),
    size = 2,
    shape = 4
  ) +
  labs(
    title = "Monte Carlo recovery: covariate-Hawkes model",
    subtitle = "Crosses indicate true parameter values",
    x = "parameter",
    y = "estimate"
  ) +
  theme_minimal()

ggsave(
  "figures/covariate_hawkes_recovery_boxplots.pdf",
  p,
  width = 7,
  height = 4
)

ggsave(
  "figures/covariate_hawkes_recovery_boxplots.png",
  p,
  width = 7,
  height = 4,
  dpi = 300
)

p_branching <- ggplot(results, aes(x = branching_hat)) +
  geom_histogram(bins = 25, boundary = 0) +
  geom_vline(
    xintercept = branching_true,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(
    title = "Recovery of the branching ratio",
    subtitle = paste0(
      "Dashed line: true branching = ",
      round(branching_true, 3)
    ),
    x = "estimated branching ratio",
    y = "count"
  ) +
  theme_minimal()

ggsave(
  "figures/covariate_hawkes_recovery_branching_histogram.pdf",
  p_branching,
  width = 7,
  height = 4
)

ggsave(
  "figures/covariate_hawkes_recovery_branching_histogram.png",
  p_branching,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/covariate_hawkes_recovery_replicates.csv\n")
cat("- results/covariate_hawkes_recovery_summary.csv\n")
cat("- figures/covariate_hawkes_recovery_boxplots.pdf/png\n")
cat("- figures/covariate_hawkes_recovery_branching_histogram.pdf/png\n")