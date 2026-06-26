# 11_sample_size_effect.R
#
# Sample-size / observation-window effect.
#
# True model:
#   lambda(t) = exp(gamma0 + gamma1 X(t)) + alpha R(t)
#
# Question:
# does recovery of the branching ratio improve when T increases?

set.seed(11)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

source("R/load_all.R")

# ------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------

dt <- 0.01
T_values <- c(50, 100, 200, 400)

gamma0_true <- -0.5
gamma1_true <- 0.9
alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

n_rep <- 30

# ------------------------------------------------------------
# 2. Simulation
# ------------------------------------------------------------

make_covariate <- function(grid) {
  sin(2 * pi * grid / 40) + 0.5 * sin(2 * pi * grid / 13)
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
# 3. Likelihood and fitting
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
    log(max(n_events, 1) / max(grid)),
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
counter <- 1

for (T_end in T_values) {
  grid <- seq(0, T_end, by = dt)
  X <- make_covariate(grid)
  
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
    
    results_list[[counter]] <- data.frame(
      T_end = T_end,
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
      "T =", T_end,
      "- rep", rep, "of", n_rep,
      "- n =", length(events),
      "- branching =", round(branching_hat, 3),
      "\n"
    )
    
    counter <- counter + 1
  }
}

results <- do.call(rbind, results_list)

# ------------------------------------------------------------
# 5. Summary
# ------------------------------------------------------------

summary_list <- list()

for (T_end in sort(unique(results$T_end))) {
  subset_T <- results[results$T_end == T_end, ]
  
  x <- subset_T$branching_hat
  
  summary_list[[as.character(T_end)]] <- data.frame(
    T_end = T_end,
    mean_n_events = mean(subset_T$n_events),
    mean_branching = mean(x),
    bias_branching = mean(x) - branching_true,
    sd_branching = sd(x),
    q05_branching = as.numeric(quantile(x, 0.05)),
    median_branching = as.numeric(quantile(x, 0.50)),
    q95_branching = as.numeric(quantile(x, 0.95))
  )
}

summary <- do.call(rbind, summary_list)
rownames(summary) <- NULL

print(summary)

write.csv(
  results,
  "results/sample_size_effect_replicates.csv",
  row.names = FALSE
)

write.csv(
  summary,
  "results/sample_size_effect_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Plots
# ------------------------------------------------------------

p_branching <- ggplot(
  results,
  aes(x = factor(T_end), y = branching_hat)
) +
  geom_boxplot(outlier.alpha = 0.35) +
  geom_hline(
    yintercept = branching_true,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(
    title = "Effect of observation time on branching recovery",
    subtitle = paste0(
      "Dashed line: true branching = ",
      round(branching_true, 3)
    ),
    x = "observation time T",
    y = "estimated branching ratio"
  ) +
  theme_minimal()

ggsave(
  "figures/sample_size_effect_branching_boxplots.pdf",
  p_branching,
  width = 7,
  height = 4
)

ggsave(
  "figures/sample_size_effect_branching_boxplots.png",
  p_branching,
  width = 7,
  height = 4,
  dpi = 300
)

p_events <- ggplot(
  summary,
  aes(x = T_end, y = mean_n_events)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    title = "Mean number of events versus observation time",
    x = "observation time T",
    y = "mean number of events"
  ) +
  theme_minimal()

ggsave(
  "figures/sample_size_effect_event_counts.pdf",
  p_events,
  width = 7,
  height = 4
)

ggsave(
  "figures/sample_size_effect_event_counts.png",
  p_events,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/sample_size_effect_replicates.csv\n")
cat("- results/sample_size_effect_summary.csv\n")
cat("- figures/sample_size_effect_branching_boxplots.pdf/png\n")
cat("- figures/sample_size_effect_event_counts.pdf/png\n")