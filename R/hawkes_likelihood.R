# ============================================================
# Likelihood functions for Hawkes models
# ============================================================
#
# Contains both early binned likelihoods and later continuous-time
# likelihoods for exponential Hawkes models.

# ============================================================
# Binned full model:
# lambda(t) = exp(mu0 + mu1 X(t)) + Hawkes memory
# ============================================================

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
# Binned null model for covariate test: H0 mu1 = 0
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

# ============================================================
# Binned null model for Hawkes test: H0 alpha = 0
# ============================================================

neg_loglik_no_hawkes <- function(par, X, dN, dt) {
  mu0 <- par[1]
  mu1 <- par[2]

  lambda <- exp(mu0 + mu1 * X)
  lambda <- pmax(lambda, 1e-12)

  loglik <- sum(dN * log(lambda) - lambda * dt)

  -loglik
}

# ============================================================
# Continuous-time exponential Hawkes log-likelihood
# theta = (log_mu, log_alpha, log_beta)
# ============================================================

hawkes_loglik_theta <- function(theta, events, T_end) {
  mu <- exp(theta[1])
  alpha <- exp(theta[2])
  beta <- exp(theta[3])

  if (alpha >= beta) {
    return(-Inf)
  }

  n <- length(events)

  if (n == 0) {
    return(-mu * T_end)
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

    if (!is.finite(lambda_events[i]) ||
        lambda_events[i] <= 0) {
      return(-Inf)
    }
  }

  integral <- mu * T_end +
    sum(alpha / beta * (1 - exp(-beta * (T_end - events))))

  sum(log(lambda_events)) - integral
}

hawkes_negloglik_theta <- function(theta, events, T_end) {
  value <- hawkes_loglik_theta(theta, events, T_end)

  if (!is.finite(value)) {
    return(1e12)
  }

  -value
}

# Parameterized by natural parameters rather than log parameters.
hawkes_loglik <- function(events, T_end, mu, alpha, beta) {
  hawkes_loglik_theta(
    theta = log(c(mu, alpha, beta)),
    events = events,
    T_end = T_end
  )
}

# ============================================================
# Continuous-time covariate-Hawkes log-likelihood
# theta = (gamma0, gamma1, log_alpha, log_beta)
# ============================================================

covariate_hawkes_loglik <- function(theta, events, grid, X) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]
  alpha <- exp(theta[3])
  beta <- exp(theta[4])

  if (alpha >= beta) {
    return(-Inf)
  }

  T_end <- max(grid)
  dt <- grid[2] - grid[1]
  n <- length(events)

  baseline_grid <- exp(gamma0 + gamma1 * X)
  baseline_integral <- sum(baseline_grid) * dt

  if (n == 0) {
    return(-baseline_integral)
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

    if (!is.finite(lambda_events[i]) ||
        lambda_events[i] <= 0) {
      return(-Inf)
    }
  }

  hawkes_integral <- sum(
    alpha / beta *
      (1 - exp(-beta * (T_end - events)))
  )

  sum(log(lambda_events)) -
    baseline_integral -
    hawkes_integral
}

covariate_hawkes_negloglik <- function(theta, events, grid, X) {
  value <- covariate_hawkes_loglik(
    theta = theta,
    events = events,
    grid = grid,
    X = X
  )

  if (!is.finite(value)) {
    return(1e12)
  }

  -value
}

# ============================================================
# Continuous-time covariate-only Poisson log-likelihood
# theta = (gamma0, gamma1)
# ============================================================

covariate_poisson_loglik <- function(theta, events, grid, X) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]

  dt <- grid[2] - grid[1]
  lambda_grid <- exp(gamma0 + gamma1 * X)
  integral <- sum(lambda_grid) * dt

  if (length(events) == 0) {
    return(-integral)
  }

  X_events <- approx(grid, X, xout = events, rule = 2)$y
  sum(gamma0 + gamma1 * X_events) - integral
}

covariate_poisson_negloglik <- function(theta, events, grid, X) {
  value <- covariate_poisson_loglik(theta, events, grid, X)

  if (!is.finite(value)) {
    return(1e12)
  }

  -value
}
