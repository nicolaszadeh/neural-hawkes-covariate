set.seed(1)

# ============================================================
# 1. Parameters
# ============================================================

T_end <- 500
dt <- 0.01
time <- seq(0, T_end, by = dt)
n_steps <- length(time)

# Ornstein-Uhlenbeck covariate parameters
gamma_true <- 0.8
sigma_true <- 1.0

# Hawkes parameters
mu0_true <- -0.5
mu1_true <- 0.7
alpha_true <- 0.7
beta_true <- 2.5

# Stability check for the one-dimensional Hawkes kernel:
# alpha / beta should be < 1.
cat("alpha / beta =", alpha_true / beta_true, "\n")


# ============================================================
# 2. Simulate the OU covariate X(t)
# ============================================================

simulate_ou <- function(time, gamma, sigma, x0 = 0) {
  dt <- time[2] - time[1]
  n <- length(time)
  
  X <- numeric(n)
  X[1] <- x0
  
  for (k in 2:n) {
    X[k] <- X[k - 1] -
      gamma * X[k - 1] * dt +
      sigma * sqrt(dt) * rnorm(1)
  }
  
  X
}

X <- simulate_ou(time, gamma_true, sigma_true)


# ============================================================
# 3. Simulate the Hawkes process with covariate baseline
# ============================================================

simulate_hawkes_discrete <- function(
    X, dt, mu0, mu1, alpha, beta
) {
  n <- length(X)
  
  dN <- integer(n)
  lambda <- numeric(n)
  H <- 0
  
  for (k in 1:n) {
    baseline <- exp(mu0 + mu1 * X[k])
    lambda[k] <- baseline + H
    
    dN[k] <- rpois(1, lambda[k] * dt)
    
    H <- H * exp(-beta * dt) + alpha * dN[k]
  }
  
  list(dN = dN, lambda = lambda)
}

sim <- simulate_hawkes_discrete(
  X = X,
  dt = dt,
  mu0 = mu0_true,
  mu1 = mu1_true,
  alpha = alpha_true,
  beta = beta_true
)

dN <- sim$dN
lambda_true <- sim$lambda

event_times <- time[dN > 0]

cat("Number of events:", sum(dN), "\n")
cat("Max number of events in one bin:", max(dN), "\n")


# ============================================================
# 4. Log-likelihood
# ============================================================
#
# We use the binned point-process likelihood:
#
#   sum_k dN_k log(lambda_k) - lambda_k dt
#
# We omit constants independent of the parameters.

neg_loglik_full <- function(par, X, dN, dt) {
  mu0 <- par[1]
  mu1 <- par[2]
  alpha <- exp(par[3])
  beta <- exp(par[4])
  
  n <- length(X)
  H <- 0
  loglik <- 0
  
  for (k in 1:n) {
    lambda <- exp(mu0 + mu1 * X[k]) + H
    lambda <- max(lambda, 1e-12)
    
    loglik <- loglik +
      dN[k] * log(lambda) -
      lambda * dt
    
    H <- H * exp(-beta * dt) + alpha * dN[k]
  }
  
  -loglik
}


# ============================================================
# 5. Fit the full model
# ============================================================

start_full <- c(
  mu0 = log(sum(dN) / T_end + 1e-6),
  mu1 = 0,
  log_alpha = log(0.2),
  log_beta = log(1.0)
)

fit_full <- optim(
  par = start_full,
  fn = neg_loglik_full,
  X = X,
  dN = dN,
  dt = dt,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 1000)
)

est_full <- c(
  mu0 = fit_full$par[1],
  mu1 = fit_full$par[2],
  alpha = exp(fit_full$par[3]),
  beta = exp(fit_full$par[4])
)

cat("\nFull model estimates:\n")
print(est_full)

cat("\nTrue values:\n")
print(c(
  mu0 = mu0_true,
  mu1 = mu1_true,
  alpha = alpha_true,
  beta = beta_true
))

loglik_full <- -fit_full$value


# ============================================================
# 6. Null model: no covariate effect, mu1 = 0
# ============================================================

neg_loglik_no_covariate <- function(par, X, dN, dt) {
  mu0 <- par[1]
  mu1 <- 0
  alpha <- exp(par[2])
  beta <- exp(par[3])
  
  n <- length(X)
  H <- 0
  loglik <- 0
  
  for (k in 1:n) {
    lambda <- exp(mu0 + mu1 * X[k]) + H
    lambda <- max(lambda, 1e-12)
    
    loglik <- loglik +
      dN[k] * log(lambda) -
      lambda * dt
    
    H <- H * exp(-beta * dt) + alpha * dN[k]
  }
  
  -loglik
}

start_no_cov <- c(
  mu0 = log(sum(dN) / T_end + 1e-6),
  log_alpha = log(0.2),
  log_beta = log(1.0)
)

fit_no_cov <- optim(
  par = start_no_cov,
  fn = neg_loglik_no_covariate,
  X = X,
  dN = dN,
  dt = dt,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 1000)
)

loglik_no_cov <- -fit_no_cov$value

LR_cov <- 2 * (loglik_full - loglik_no_cov)
p_cov <- pchisq(LR_cov, df = 1, lower.tail = FALSE)

cat("\nTest of covariate relevance:\n")
cat("H0: mu1 = 0\n")
cat("LR statistic =", LR_cov, "\n")
cat("p-value =", p_cov, "\n")


# ============================================================
# 7. Null model: no Hawkes self-excitation, alpha = 0
# ============================================================

neg_loglik_no_hawkes <- function(par, X, dN, dt) {
  mu0 <- par[1]
  mu1 <- par[2]
  
  lambda <- exp(mu0 + mu1 * X)
  lambda <- pmax(lambda, 1e-12)
  
  loglik <- sum(dN * log(lambda) - lambda * dt)
  
  -loglik
}

start_no_hawkes <- c(
  mu0 = log(sum(dN) / T_end + 1e-6),
  mu1 = 0
)

fit_no_hawkes <- optim(
  par = start_no_hawkes,
  fn = neg_loglik_no_hawkes,
  X = X,
  dN = dN,
  dt = dt,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 1000)
)

loglik_no_hawkes <- -fit_no_hawkes$value

LR_hawkes <- 2 * (loglik_full - loglik_no_hawkes)
p_hawkes <- pchisq(LR_hawkes, df = 1, lower.tail = FALSE)

cat("\nTest of Hawkes self-excitation:\n")
cat("H0: alpha = 0\n")
cat("LR statistic =", LR_hawkes, "\n")
cat("p-value =", p_hawkes, "\n")
cat("Warning: this p-value is approximate because alpha = 0\n")
cat("is on the boundary of the parameter space.\n")


# ============================================================
# 8. Reconstruct fitted intensity
# ============================================================

compute_intensity <- function(X, dN, dt, mu0, mu1, alpha, beta) {
  n <- length(X)
  lambda <- numeric(n)
  H <- 0
  
  for (k in 1:n) {
    lambda[k] <- exp(mu0 + mu1 * X[k]) + H
    H <- H * exp(-beta * dt) + alpha * dN[k]
  }
  
  lambda
}

lambda_fit <- compute_intensity(
  X = X,
  dN = dN,
  dt = dt,
  mu0 = est_full["mu0"],
  mu1 = est_full["mu1"],
  alpha = est_full["alpha"],
  beta = est_full["beta"]
)


# ============================================================
# 9. Plots
# ============================================================

par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

plot(
  time, X,
  type = "l",
  xlab = "time",
  ylab = "X(t)",
  main = "OU covariate"
)

plot(
  time, lambda_true,
  type = "l",
  xlab = "time",
  ylab = "lambda(t)",
  main = "True and fitted intensity"
)

lines(time, lambda_fit, lty = 2)

legend(
  "topright",
  legend = c("true", "fitted"),
  lty = c(1, 2),
  bty = "n"
)

plot(
  time, dN,
  type = "h",
  xlab = "time",
  ylab = "dN",
  main = "Observed events"
)

par(mfrow = c(1, 1))