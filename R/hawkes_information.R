# ============================================================
# Information-matrix utilities
# ============================================================

observed_information_from_hessian <- function(
    f,
    theta_hat,
    T_end,
    eps = 1e-4
) {
  if (!exists("numerical_hessian")) {
    stop("Please source R/hawkes_gradient.R first.")
  }

  H_total <- numerical_hessian(
    f = f,
    theta = theta_hat,
    eps = eps
  )

  I_hat <- H_total / T_end
  I_inv_hat <- solve(I_hat)

  list(
    H_total = H_total,
    I_hat = I_hat,
    I_inv_hat = I_inv_hat,
    cov_total = solve(H_total)
  )
}

matrix_inv_sqrt <- function(M, ridge = 1e-10) {
  M_sym <- (M + t(M)) / 2
  eig <- eigen(M_sym)

  values <- pmax(eig$values, ridge)
  vectors <- eig$vectors

  out <- vectors %*%
    diag(1 / sqrt(values), nrow = length(values)) %*%
    t(vectors)

  rownames(out) <- rownames(M)
  colnames(out) <- colnames(M)

  out
}

information_summary <- function(theta_hat, H_total, I_inv_hat) {
  out <- data.frame(
    parameter = names(theta_hat),
    se_total = sqrt(diag(solve(H_total))),
    sigma_hat_asymptotic = sqrt(diag(I_inv_hat))
  )

  rownames(out) <- NULL
  out
}

hawkes_observed_information <- function(events, T_end, theta_hat) {
  f <- function(theta) {
    hawkes_negloglik_theta(
      theta = theta,
      events = events,
      T_end = T_end
    )
  }

  observed_information_from_hessian(
    f = f,
    theta_hat = theta_hat,
    T_end = T_end
  )
}

covariate_hawkes_observed_information <- function(
    events,
    grid,
    X,
    theta_hat
) {
  T_end <- max(grid)

  f <- function(theta) {
    covariate_hawkes_negloglik(
      theta = theta,
      events = events,
      grid = grid,
      X = X
    )
  }

  observed_information_from_hessian(
    f = f,
    theta_hat = theta_hat,
    T_end = T_end
  )
}
