# 09_omitted_covariate_bias.R
#
# Omitted-covariate bias experiment.
#
# True model:
#   lambda(t) = exp(gamma0 + gamma1 X(t))
#
# Fitted wrong model:
#   lambda(t) = mu + alpha R(t)
#
# Question:
# if we ignore X(t), does the Hawkes model invent
# self-excitation?

set.seed(9)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

source("R/load_all.R")

# ------------------------------------------------------------
# 1. Grid and covariate
# ------------------------------------------------------------

T_end <- 200
dt <- 0.01
grid <- seq(0, T_end, by = dt)

X <- sin(2 * pi * grid / 40) + 0.5 * sin(2 * pi * grid / 13)

gamma0_true <- -0.3
gamma1_values <- c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)

n_rep <- 30

# ------------------------------------------------------------
# 2. Simulate covariate-only Poisson process
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
# 3. Naive Hawkes likelihood
# ------------------------------------------------------------

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
  lambda_events <- numeric(n)
  
  for (i in seq_len(n)) {
    if (i == 1) {
      R[i] <- 0
    } else {
      delta <- events[i] - events[i - 1]
      R[i] <- exp(-beta * delta) * (1 + R[i - 1])
    }
    
    lambda_events[i] <- mu + alpha * R[i]
    
    if (!is.finite(lambda_events[i]) || lambda_events[i] <= 0) {
      return(1e12)
    }
  }
  
  integral <- mu * T_end +
    sum(alpha / beta * (1 - exp(-beta * (T_end - events))))
  
  -(sum(log(lambda_events)) - integral)
}

fit_naive_hawkes <- function(events, T_end) {
  n_events <- length(events)
  
  if (n_events == 0) {
    return(NULL)
  }
  
  init <- c(
    log(n_events / T_end),
    log(0.2),
    log(1.0)
  )
  
  optim(
    par = init,
    fn = hawkes_naive_negloglik,
    events = events,
    T_end = T_end,
    method = "Nelder-Mead",
    control = list(maxit = 5000)
  )
}

# ------------------------------------------------------------
# 4. Monte Carlo loop
# ------------------------------------------------------------

results_list <- list()
counter <- 1

for (gamma1 in gamma1_values) {
  for (rep in seq_len(n_rep)) {
    events <- simulate_covariate_poisson(
      grid = grid,
      X = X,
      gamma0 = gamma0_true,
      gamma1 = gamma1
    )
    
    fit <- fit_naive_hawkes(events, T_end)
    
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
# 5. Summary
# ------------------------------------------------------------

results <- as.data.frame(results)
results$gamma1_true <- as.numeric(results$gamma1_true)
results$branching_hat <- as.numeric(results$branching_hat)

summary_list <- list()

for (g in sort(unique(results$gamma1_true))) {
  x <- results$branching_hat[results$gamma1_true == g]
  
  summary_list[[as.character(g)]] <- data.frame(
    gamma1_true = g,
    mean_branching = mean(x),
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
  "results/omitted_covariate_bias_replicates.csv",
  row.names = FALSE
)

write.csv(
  summary,
  "results/omitted_covariate_bias_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Plot
# ------------------------------------------------------------

p <- ggplot(results, aes(x = gamma1_true, y = branching_hat)) +
  geom_point(
    alpha = 0.35,
    position = position_jitter(width = 0.04, height = 0)
  ) +
  geom_line(
    data = summary,
    aes(x = gamma1_true, y = mean_branching),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  geom_point(
    data = summary,
    aes(x = gamma1_true, y = mean_branching),
    inherit.aes = FALSE,
    size = 2
  ) +
  labs(
    title = "Omitted-covariate bias",
    subtitle = paste(
      "True Hawkes excitation is zero;",
      "naive Hawkes fit invents branching when X(t) is ignored"
    ),
    x = expression("covariate strength " * gamma[1]),
    y = "estimated branching ratio"
  ) +
  theme_minimal()

ggsave(
  "figures/omitted_covariate_bias_branching.pdf",
  p,
  width = 7,
  height = 4
)

ggsave(
  "figures/omitted_covariate_bias_branching.png",
  p,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/omitted_covariate_bias_replicates.csv\n")
cat("- results/omitted_covariate_bias_summary.csv\n")
cat("- figures/omitted_covariate_bias_branching.pdf/png\n")