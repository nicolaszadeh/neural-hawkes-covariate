# 13_false_positive_calibration.R
#
# False-positive calibration under H0.
#
# True model:
#   lambda(t) = exp(gamma0 + gamma1 X(t))
#
# There is no Hawkes excitation.
#
# Question:
# if we test at level 5%, how often do we falsely reject H0?

set.seed(13)

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

alpha_level <- 0.05

B_values <- c(20, 50, 100, 200)
n_rep <- 30

# ------------------------------------------------------------
# 2. Simulation under H0
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
# 4. Fitting and LR
# ------------------------------------------------------------

fit_covariate_poisson <- function(events, grid, X) {
  init <- c(log(max(length(events), 1) / max(grid)), 0)
  
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
    control = list(maxit = 8000)
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
    branching_hat = alpha_hat / beta_hat
  )
}

bootstrap_p_value <- function(events_obs, grid, X, B) {
  obs <- compute_LR(events_obs, grid, X)
  
  gamma0_null_hat <- obs$fit_null$par[1]
  gamma1_null_hat <- obs$fit_null$par[2]
  
  boot_LR <- numeric(B)
  
  for (b in seq_len(B)) {
    boot_events <- simulate_covariate_poisson(
      grid = grid,
      X = X,
      gamma0 = gamma0_null_hat,
      gamma1 = gamma1_null_hat
    )
    
    boot <- compute_LR(boot_events, grid, X)
    boot_LR[b] <- boot$LR
  }
  
  n_exceed <- sum(boot_LR >= obs$LR)
  
  list(
    observed_LR = obs$LR,
    branching_hat = obs$branching_hat,
    p_empirical = n_exceed / B,
    p_corrected = (1 + n_exceed) / (B + 1),
    max_boot_LR = max(boot_LR),
    mean_boot_LR = mean(boot_LR)
  )
}

# ------------------------------------------------------------
# 5. Calibration loop under H0
# ------------------------------------------------------------

results_list <- list()
counter <- 1

for (B in B_values) {
  for (rep in seq_len(n_rep)) {
    events <- simulate_covariate_poisson(
      grid = grid,
      X = X,
      gamma0 = gamma0_true,
      gamma1 = gamma1_true
    )
    
    test <- bootstrap_p_value(
      events_obs = events,
      grid = grid,
      X = X,
      B = B
    )
    
    reject <- test$p_corrected < alpha_level
    
    results_list[[counter]] <- data.frame(
      B = B,
      replicate = rep,
      n_events = length(events),
      observed_LR = test$observed_LR,
      branching_hat = test$branching_hat,
      p_empirical = test$p_empirical,
      p_corrected = test$p_corrected,
      reject = reject,
      max_boot_LR = test$max_boot_LR,
      mean_boot_LR = test$mean_boot_LR
    )
    
    cat(
      "B =", B,
      "- rep", rep, "of", n_rep,
      "- n =", length(events),
      "- LR =", round(test$observed_LR, 3),
      "- p =", round(test$p_corrected, 3),
      "- reject =", reject,
      "\n"
    )
    
    counter <- counter + 1
  }
}

results <- do.call(rbind, results_list)

# ------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------

summary_list <- list()

for (B in sort(unique(results$B))) {
  subset_B <- results[results$B == B, ]
  
  summary_list[[as.character(B)]] <- data.frame(
    B = B,
    mean_n_events = mean(subset_B$n_events),
    false_positive_rate = mean(subset_B$reject),
    mean_p_value = mean(subset_B$p_corrected),
    median_p_value = median(subset_B$p_corrected),
    mean_branching_hat = mean(subset_B$branching_hat),
    sd_branching_hat = sd(subset_B$branching_hat),
    mean_observed_LR = mean(subset_B$observed_LR)
  )
}

summary <- do.call(rbind, summary_list)
rownames(summary) <- NULL

print(summary)

write.csv(
  results,
  "results/false_positive_calibration_replicates.csv",
  row.names = FALSE
)

write.csv(
  summary,
  "results/false_positive_calibration_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Plots
# ------------------------------------------------------------

p_fpr <- ggplot(summary, aes(x = B, y = false_positive_rate)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = alpha_level,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  labs(
    title = "False-positive calibration under H0",
    subtitle = paste0(
      "Dashed line: nominal level ",
      alpha_level
    ),
    x = "number of bootstrap replicates B",
    y = "false-positive rate"
  ) +
  theme_minimal()

ggsave(
  "figures/false_positive_calibration_rate.pdf",
  p_fpr,
  width = 7,
  height = 4
)

ggsave(
  "figures/false_positive_calibration_rate.png",
  p_fpr,
  width = 7,
  height = 4,
  dpi = 300
)

p_pvalues <- ggplot(
  results,
  aes(x = p_corrected)
) +
  geom_histogram(bins = 20, boundary = 0) +
  facet_wrap(~ B) +
  labs(
    title = "Bootstrap p-values under H0",
    subtitle = "Under correct calibration, p-values should be roughly uniform",
    x = "corrected bootstrap p-value",
    y = "count"
  ) +
  theme_minimal()

ggsave(
  "figures/false_positive_calibration_pvalues.pdf",
  p_pvalues,
  width = 8,
  height = 5
)

ggsave(
  "figures/false_positive_calibration_pvalues.png",
  p_pvalues,
  width = 8,
  height = 5,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/false_positive_calibration_replicates.csv\n")
cat("- results/false_positive_calibration_summary.csv\n")
cat("- figures/false_positive_calibration_rate.pdf/png\n")
cat("- figures/false_positive_calibration_pvalues.pdf/png\n")