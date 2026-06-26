# 15_ou_omitted_covariate_bias.R
#
# True model:
#   lambda(t) = exp(gamma0 + gamma1 X(t))
#
# where X(t) is an OU process.
#
# Fitted wrong model:
#   lambda(t) = mu + alpha R(t)
#
# Question:
# If we ignore X(t), does a Hawkes model invent
# self-excitation?

set.seed(15)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

source("R/load_all.R")

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

T_end <- 200
dt <- 0.01

gamma0_true <- -0.3

gamma1_values <- c(
  0,
  0.25,
  0.5,
  0.75,
  1.0,
  1.25,
  1.5,
  2.0
)

kappa_true <- 0.4
sigma_true <- 0.8

n_rep <- 30

# ------------------------------------------------------------
# OU simulation
# ------------------------------------------------------------

simulate_ou <- function(T_end, dt, kappa, sigma) {
  grid <- seq(0, T_end, by = dt)
  X <- numeric(length(grid))
  
  for (i in 2:length(grid)) {
    X[i] <- X[i - 1] -
      kappa * X[i - 1] * dt +
      sigma * sqrt(dt) * rnorm(1)
  }
  
  list(
    grid = grid,
    X = X
  )
}

# ------------------------------------------------------------
# Covariate Poisson simulation
# ------------------------------------------------------------

simulate_covariate_poisson <- function(
    grid,
    X,
    gamma0,
    gamma1
) {
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
    
    X_t <- approx(
      grid,
      X,
      xout = t,
      rule = 2
    )$y
    
    lambda_t <- exp(gamma0 + gamma1 * X_t)
    
    if (runif(1) <= lambda_t / lambda_max) {
      events <- c(events, t)
    }
  }
  
  events
}

# ------------------------------------------------------------
# Naive Hawkes likelihood
# ------------------------------------------------------------

hawkes_negloglik <- function(
    theta,
    events,
    T_end
) {
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
  lambda_events <- numeric(n)
  
  for (i in seq_len(n)) {
    if (i == 1) {
      R[i] <- 0
    } else {
      delta <- events[i] - events[i - 1]
      
      R[i] <- exp(-beta * delta) *
        (1 + R[i - 1])
    }
    
    lambda_events[i] <- mu + alpha * R[i]
    
    if (
      !is.finite(lambda_events[i]) ||
      lambda_events[i] <= 0
    ) {
      return(1e12)
    }
  }
  
  integral <- mu * T_end +
    sum(
      alpha / beta *
        (1 - exp(-beta * (T_end - events)))
    )
  
  -(
    sum(log(lambda_events)) -
      integral
  )
}

fit_hawkes <- function(events, T_end) {
  if (length(events) == 0) {
    return(NULL)
  }
  
  init <- c(
    log(length(events) / T_end),
    log(0.2),
    log(1)
  )
  
  optim(
    par = init,
    fn = hawkes_negloglik,
    events = events,
    T_end = T_end,
    method = "Nelder-Mead",
    control = list(maxit = 5000)
  )
}

# ------------------------------------------------------------
# Monte Carlo
# ------------------------------------------------------------

results_list <- list()
counter <- 1

for (gamma1 in gamma1_values) {
  for (rep in seq_len(n_rep)) {
    ou <- simulate_ou(
      T_end,
      dt,
      kappa_true,
      sigma_true
    )
    
    events <- simulate_covariate_poisson(
      ou$grid,
      ou$X,
      gamma0_true,
      gamma1
    )
    
    fit <- fit_hawkes(
      events,
      T_end
    )
    
    if (is.null(fit)) {
      next
    }
    
    mu_hat <- exp(fit$par[1])
    alpha_hat <- exp(fit$par[2])
    beta_hat <- exp(fit$par[3])
    branching_hat <- alpha_hat / beta_hat
    
    results_list[[counter]] <- data.frame(
      gamma1_true = gamma1,
      replicate = rep,
      n_events = length(events),
      mu_hat = mu_hat,
      alpha_hat = alpha_hat,
      beta_hat = beta_hat,
      branching_hat = branching_hat,
      neg_loglik = fit$value,
      convergence = fit$convergence
    )
    
    cat(
      "gamma1 =", gamma1,
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
# Summary
# ------------------------------------------------------------

results <- as.data.frame(results)
results$gamma1_true <- as.numeric(results$gamma1_true)
results$branching_hat <- as.numeric(results$branching_hat)

summary_list <- list()

for (g in sort(unique(results$gamma1_true))) {
  x <- results$branching_hat[results$gamma1_true == g]
  n <- results$n_events[results$gamma1_true == g]
  
  summary_list[[as.character(g)]] <- data.frame(
    gamma1_true = g,
    mean_n_events = mean(n),
    mean_branching = mean(x),
    sd_branching = sd(x),
    q05_branching = as.numeric(quantile(x, 0.05)),
    median_branching = as.numeric(quantile(x, 0.50)),
    q95_branching = as.numeric(quantile(x, 0.95))
  )
}

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

print(summary_df)

write.csv(
  results,
  "results/ou_omitted_covariate_bias_replicates.csv",
  row.names = FALSE
)

write.csv(
  summary_df,
  "results/ou_omitted_covariate_bias_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot(
  results,
  aes(
    x = gamma1_true,
    y = branching_hat
  )
) +
  geom_point(
    alpha = 0.35,
    position = position_jitter(
      width = 0.04,
      height = 0
    )
  ) +
  geom_line(
    data = summary_df,
    aes(
      x = gamma1_true,
      y = mean_branching
    ),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  geom_point(
    data = summary_df,
    aes(
      x = gamma1_true,
      y = mean_branching
    ),
    inherit.aes = FALSE,
    size = 2
  ) +
  labs(
    title = "OU omitted-covariate bias",
    subtitle =
      "Ignoring a stochastic covariate creates fake Hawkes branching",
    x = expression(gamma[1]),
    y = "estimated branching ratio"
  ) +
  theme_minimal()

ggsave(
  "figures/ou_omitted_covariate_bias.pdf",
  p,
  width = 7,
  height = 4
)

ggsave(
  "figures/ou_omitted_covariate_bias.png",
  p,
  width = 7,
  height = 4,
  dpi = 300
)

p_events <- ggplot(
  summary_df,
  aes(
    x = gamma1_true,
    y = mean_n_events
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    title = "Mean event count under OU covariate drive",
    x = expression(gamma[1]),
    y = "mean number of events"
  ) +
  theme_minimal()

ggsave(
  "figures/ou_omitted_covariate_event_counts.pdf",
  p_events,
  width = 7,
  height = 4
)

ggsave(
  "figures/ou_omitted_covariate_event_counts.png",
  p_events,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/ou_omitted_covariate_bias_replicates.csv\n")
cat("- results/ou_omitted_covariate_bias_summary.csv\n")
cat("- figures/ou_omitted_covariate_bias.pdf/png\n")
cat("- figures/ou_omitted_covariate_event_counts.pdf/png\n")