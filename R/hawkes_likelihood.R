# ============================================================
# Likelihood functions for binned Hawkes models
# ============================================================


# ============================================================
# Full model:
#
# lambda(t)
# =
# exp(mu0 + mu1 X(t)) + Hawkes memory
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
# Null model for covariate test:
#
# H0: mu1 = 0
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
# Null model for Hawkes test:
#
# H0: alpha = 0
# ============================================================

neg_loglik_no_hawkes <- function(par, X, dN, dt) {
  mu0 <- par[1]
  mu1 <- par[2]
  
  lambda <- exp(mu0 + mu1 * X)
  lambda <- pmax(lambda, 1e-12)
  
  loglik <- sum(dN * log(lambda) - lambda * dt)
  
  -loglik
}