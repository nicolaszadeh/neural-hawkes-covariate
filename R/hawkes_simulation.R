# ============================================================
# Simulation and intensity utilities for Hawkes models
# ============================================================
#
# This file contains both the original discrete-time helpers used
# in the early scripts and the continuous-time helpers introduced
# later in the project.

# ============================================================
# Ornstein-Uhlenbeck covariate
# ============================================================

simulate_ou <- function(
    time = NULL,
    grid = NULL,
    gamma = NULL,
    sigma,
    x0 = 0,
    stationary_start = FALSE,
    kappa = NULL,
    X0 = NULL
) {
  if (is.null(time)) {
    time <- grid
  }
  # Backward compatibility:
  # earlier scripts used simulate_ou(time, gamma, sigma, x0),
  # later scripts used simulate_ou(grid, kappa, sigma, X0).
  if (is.null(gamma)) {
    gamma <- kappa
  }
  if (!is.null(X0)) {
    x0 <- X0
  }

  dt <- time[2] - time[1]
  n <- length(time)

  if (gamma <= 0) {
    stop("gamma/kappa must be positive.")
  }

  X <- numeric(n)

  if (stationary_start) {
    X[1] <- rnorm(
      1,
      mean = 0,
      sd = sigma / sqrt(2 * gamma)
    )
  } else {
    X[1] <- x0
  }

  # Exact OU transition on a regular grid.
  a <- exp(-gamma * dt)
  innovation_sd <- sigma *
    sqrt((1 - exp(-2 * gamma * dt)) / (2 * gamma))

  for (k in 2:n) {
    X[k] <- a * X[k - 1] + innovation_sd * rnorm(1)
  }

  X
}

# ============================================================
# Continuous-time univariate exponential Hawkes simulation
# ============================================================

simulate_hawkes <- function(T_end, mu, alpha, beta) {
  if (mu <= 0 || alpha < 0 || beta <= 0) {
    stop("Require mu > 0, alpha >= 0, beta > 0.")
  }
  if (alpha >= beta) {
    warning("alpha >= beta: branching ratio >= 1; simulation may be unstable.")
  }

  t <- 0
  events <- numeric(0)
  R <- 0

  while (t < T_end) {
    lambda_upper <- mu + alpha * R

    if (!is.finite(lambda_upper) || lambda_upper <= 0) {
      stop("Invalid upper intensity.")
    }

    t_candidate <- t + rexp(1, rate = lambda_upper)

    if (t_candidate > T_end) {
      break
    }

    R_candidate <- R * exp(-beta * (t_candidate - t))
    lambda_candidate <- mu + alpha * R_candidate

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

# ============================================================
# Continuous-time covariate-only Poisson simulation
# ============================================================

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

# ============================================================
# Continuous-time covariate-Hawkes simulation
# ============================================================

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

# ============================================================
# Discrete-time Hawkes simulation
# ============================================================

simulate_hawkes_discrete <- function(
    X,
    dt,
    mu0,
    mu1,
    alpha,
    beta
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

  list(
    dN = dN,
    lambda = lambda
  )
}

# ============================================================
# Simulation under the no-Hawkes null model alpha = 0
# ============================================================

simulate_no_hawkes <- function(
    X,
    dt,
    mu0,
    mu1
) {
  lambda <- exp(mu0 + mu1 * X)
  dN <- rpois(length(X), lambda * dt)

  list(
    dN = dN,
    lambda = lambda
  )
}

# ============================================================
# Reconstruct intensity from binned event counts
# ============================================================

compute_intensity_from_binned_counts <- function(
    X,
    dN,
    dt,
    mu0,
    mu1,
    alpha,
    beta
) {
  n <- length(X)
  lambda <- numeric(n)
  H <- 0

  for (k in 1:n) {
    lambda[k] <- exp(mu0 + mu1 * X[k]) + H
    H <- H * exp(-beta * dt) + alpha * dN[k]
  }

  lambda
}

compute_intensity <- compute_intensity_from_binned_counts
