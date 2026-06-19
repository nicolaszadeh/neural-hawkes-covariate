# ============================================================
# Simulation and intensity utilities for Hawkes models
# ============================================================


# ============================================================
# Exact simulation of an Ornstein-Uhlenbeck covariate
# ============================================================

simulate_ou <- function(
    time,
    gamma,
    sigma,
    x0 = 0,
    stationary_start = FALSE
) {
  dt <- time[2] - time[1]
  n <- length(time)
  
  if (gamma <= 0) {
    stop("gamma must be positive.")
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
  
  a <- exp(-gamma * dt)
  
  innovation_sd <- sigma *
    sqrt((1 - exp(-2 * gamma * dt)) / (2 * gamma))
  
  for (k in 2:n) {
    X[k] <- a * X[k - 1] +
      innovation_sd * rnorm(1)
  }
  
  X
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


# ============================================================
# Backward-compatible shorter alias
# ============================================================

compute_intensity <- compute_intensity_from_binned_counts