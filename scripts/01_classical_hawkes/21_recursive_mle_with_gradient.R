# 21_recursive_mle_with_gradient.R
#
# Refactored version.
# Compares Nelder-Mead likelihood-only optimization with
# BFGS using the recursive analytical score.

set.seed(21)

source("R/load_all.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

T_end <- 1000

mu_true <- 0.8
alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

events <- simulate_hawkes(
  T_end = T_end,
  mu = mu_true,
  alpha = alpha_true,
  beta = beta_true
)

cat("Number of events:", length(events), "\n")

theta_init <- log(c(
  length(events) / T_end,
  0.2,
  1.0
))

cat("Initial theta:", theta_init, "\n")

time_nm <- system.time({
  fit_nm <- optim(
    par = theta_init,
    fn = hawkes_negloglik_theta,
    events = events,
    T_end = T_end,
    method = "Nelder-Mead",
    control = list(maxit = 10000)
  )
})

time_bfgs <- system.time({
  fit_bfgs <- optim(
    par = theta_init,
    fn = hawkes_negloglik_theta,
    gr = hawkes_negscore_theta,
    events = events,
    T_end = T_end,
    method = "BFGS",
    control = list(maxit = 10000)
  )
})

extract_fit <- function(fit, method, elapsed_time) {
  mu_hat <- exp(fit$par[1])
  alpha_hat <- exp(fit$par[2])
  beta_hat <- exp(fit$par[3])
  branching_hat <- alpha_hat / beta_hat

  data.frame(
    method = method,
    neg_loglik = fit$value,
    convergence = fit$convergence,
    elapsed_time = as.numeric(elapsed_time["elapsed"]),
    mu_hat = mu_hat,
    alpha_hat = alpha_hat,
    beta_hat = beta_hat,
    branching_hat = branching_hat,
    log_mu_hat = fit$par[1],
    log_alpha_hat = fit$par[2],
    log_beta_hat = fit$par[3]
  )
}

fit_summary <- rbind(
  extract_fit(fit_nm, "Nelder-Mead", time_nm),
  extract_fit(fit_bfgs, "BFGS_gradient", time_bfgs)
)

truth <- data.frame(
  parameter = c("mu", "alpha", "beta", "branching"),
  truth = c(
    mu_true,
    alpha_true,
    beta_true,
    branching_true
  )
)

score_nm <- hawkes_score_theta(
  fit_nm$par,
  events,
  T_end
)

score_bfgs <- hawkes_score_theta(
  fit_bfgs$par,
  events,
  T_end
)

gradient_summary <- data.frame(
  method = c("Nelder-Mead", "BFGS_gradient"),
  max_abs_score = c(
    max(abs(score_nm)),
    max(abs(score_bfgs))
  ),
  score_log_mu = c(score_nm[1], score_bfgs[1]),
  score_log_alpha = c(score_nm[2], score_bfgs[2]),
  score_log_beta = c(score_nm[3], score_bfgs[3])
)

print(fit_summary)
print(truth)
print(gradient_summary)

write.csv(
  fit_summary,
  "results/recursive_mle_with_gradient_fit_summary.csv",
  row.names = FALSE
)

write.csv(
  truth,
  "results/recursive_mle_with_gradient_truth.csv",
  row.names = FALSE
)

write.csv(
  gradient_summary,
  "results/recursive_mle_with_gradient_score_summary.csv",
  row.names = FALSE
)

estimate_long <- rbind(
  data.frame(
    method = fit_summary$method,
    parameter = "mu",
    estimate = fit_summary$mu_hat,
    truth = mu_true
  ),
  data.frame(
    method = fit_summary$method,
    parameter = "alpha",
    estimate = fit_summary$alpha_hat,
    truth = alpha_true
  ),
  data.frame(
    method = fit_summary$method,
    parameter = "beta",
    estimate = fit_summary$beta_hat,
    truth = beta_true
  ),
  data.frame(
    method = fit_summary$method,
    parameter = "branching",
    estimate = fit_summary$branching_hat,
    truth = branching_true
  )
)

truth_long <- unique(
  estimate_long[, c("parameter", "truth")]
)

p_est <- ggplot(
  estimate_long,
  aes(
    x = parameter,
    y = estimate,
    fill = method
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_point(
    data = truth_long,
    aes(
      x = parameter,
      y = truth
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3
  ) +
  labs(
    title = "MLE estimates with and without analytical gradient",
    subtitle = "Crosses indicate true parameter values",
    x = "parameter",
    y = "estimate"
  ) +
  theme_minimal()

ggsave(
  "figures/recursive_mle_with_gradient_estimates.pdf",
  p_est,
  width = 7,
  height = 4
)

ggsave(
  "figures/recursive_mle_with_gradient_estimates.png",
  p_est,
  width = 7,
  height = 4,
  dpi = 300
)

p_time <- ggplot(
  fit_summary,
  aes(
    x = method,
    y = elapsed_time
  )
) +
  geom_col(width = 0.6) +
  labs(
    title = "Optimization time",
    x = "method",
    y = "elapsed time in seconds"
  ) +
  theme_minimal()

ggsave(
  "figures/recursive_mle_with_gradient_time.pdf",
  p_time,
  width = 6,
  height = 4
)

ggsave(
  "figures/recursive_mle_with_gradient_time.png",
  p_time,
  width = 6,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/recursive_mle_with_gradient_fit_summary.csv\n")
cat("- results/recursive_mle_with_gradient_truth.csv\n")
cat("- results/recursive_mle_with_gradient_score_summary.csv\n")
cat("- figures/recursive_mle_with_gradient_estimates.pdf/png\n")
cat("- figures/recursive_mle_with_gradient_time.pdf/png\n")
