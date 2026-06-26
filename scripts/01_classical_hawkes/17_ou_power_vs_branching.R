# 17_ou_power_vs_branching.R
#
# OU power study:
# How often does the OU covariate-adjusted bootstrap LR test detect
# Hawkes excitation as the true branching ratio increases?

set.seed(17)

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

beta_true <- 1.2

branching_values <- c(
  0,
  0.05,
  0.10,
  0.20,
  0.30,
  0.40,
  0.50
)

kappa_true <- 0.4
sigma_true <- 0.8
X0 <- 0

n_rep <- 30
B_bootstrap <- 50
alpha_level <- 0.05

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

# ------------------------------------------------------------
# 3. Point-process simulation
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
# 5. Fitting and LR
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
    branching_hat = alpha_hat / beta_hat,
    alpha_hat = alpha_hat,
    beta_hat = beta_hat
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
    observed_branching_hat = obs$branching_hat,
    observed_alpha_hat = obs$alpha_hat,
    observed_beta_hat = obs$beta_hat,
    p_empirical = n_exceed / B,
    p_corrected = (1 + n_exceed) / (B + 1),
    max_boot_LR = max(boot_LR),
    mean_boot_LR = mean(boot_LR)
  )
}

# ------------------------------------------------------------
# 6. Power loop
# ------------------------------------------------------------

results_list <- list()
counter <- 1

for (branching_true in branching_values) {
  alpha_true <- beta_true * branching_true
  
  for (rep in seq_len(n_rep)) {
    X <- simulate_ou(
      grid = grid,
      kappa = kappa_true,
      sigma = sigma_true,
      X0 = X0
    )
    
    if (branching_true == 0) {
      events <- simulate_covariate_poisson(
        grid = grid,
        X = X,
        gamma0 = gamma0_true,
        gamma1 = gamma1_true
      )
    } else {
      events <- simulate_covariate_hawkes(
        grid = grid,
        X = X,
        gamma0 = gamma0_true,
        gamma1 = gamma1_true,
        alpha = alpha_true,
        beta = beta_true
      )
    }
    
    test <- bootstrap_p_value(
      events_obs = events,
      grid = grid,
      X = X,
      B = B_bootstrap
    )
    
    reject <- test$p_corrected < alpha_level
    
    results_list[[counter]] <- data.frame(
      branching_true = branching_true,
      alpha_true = alpha_true,
      beta_true = beta_true,
      replicate = rep,
      n_events = length(events),
      observed_LR = test$observed_LR,
      branching_hat = test$observed_branching_hat,
      alpha_hat = test$observed_alpha_hat,
      beta_hat = test$observed_beta_hat,
      p_empirical = test$p_empirical,
      p_corrected = test$p_corrected,
      reject = reject,
      max_boot_LR = test$max_boot_LR,
      mean_boot_LR = test$mean_boot_LR
    )
    
    cat(
      "branching =", branching_true,
      "- rep", rep, "of", n_rep,
      "- n =", length(events),
      "- LR =", round(test$observed_LR, 2),
      "- p =", round(test$p_corrected, 3),
      "- reject =", reject,
      "\n"
    )
    
    counter <- counter + 1
  }
}

results <- do.call(rbind, results_list)

# ------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------

summary_list <- list()

for (br in sort(unique(results$branching_true))) {
  subset_br <- results[results$branching_true == br, ]
  
  summary_list[[as.character(br)]] <- data.frame(
    branching_true = br,
    mean_n_events = mean(subset_br$n_events),
    rejection_rate = mean(subset_br$reject),
    mean_p_value = mean(subset_br$p_corrected),
    median_p_value = median(subset_br$p_corrected),
    mean_branching_hat = mean(subset_br$branching_hat),
    sd_branching_hat = sd(subset_br$branching_hat),
    mean_observed_LR = mean(subset_br$observed_LR)
  )
}

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

print(summary_df)

write.csv(
  results,
  "results/ou_power_vs_branching_replicates.csv",
  row.names = FALSE
)

write.csv(
  summary_df,
  "results/ou_power_vs_branching_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Plots
# ------------------------------------------------------------

p_power <- ggplot(
  summary_df,
  aes(x = branching_true, y = rejection_rate)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = alpha_level,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  labs(
    title = "Power of the OU covariate-adjusted Hawkes test",
    subtitle = paste0(
      "Dashed line: nominal level ",
      alpha_level
    ),
    x = "true branching ratio",
    y = "rejection rate"
  ) +
  theme_minimal()

ggsave(
  "figures/ou_power_vs_branching_rejection_rate.pdf",
  p_power,
  width = 7,
  height = 4
)

ggsave(
  "figures/ou_power_vs_branching_rejection_rate.png",
  p_power,
  width = 7,
  height = 4,
  dpi = 300
)

p_branching <- ggplot(
  results,
  aes(x = factor(branching_true), y = branching_hat)
) +
  geom_boxplot(outlier.alpha = 0.35) +
  labs(
    title = "OU estimated branching versus true branching",
    x = "true branching ratio",
    y = "estimated branching ratio"
  ) +
  theme_minimal()

ggsave(
  "figures/ou_power_vs_branching_estimates.pdf",
  p_branching,
  width = 7,
  height = 4
)

ggsave(
  "figures/ou_power_vs_branching_estimates.png",
  p_branching,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/ou_power_vs_branching_replicates.csv\n")
cat("- results/ou_power_vs_branching_summary.csv\n")
cat("- figures/ou_power_vs_branching_rejection_rate.pdf/png\n")
cat("- figures/ou_power_vs_branching_estimates.pdf/png\n")