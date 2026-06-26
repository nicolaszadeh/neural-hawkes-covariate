# ============================================================
# Gradients, recursive quantities, and numerical derivatives
# ============================================================

hawkes_score_theta <- function(theta, events, T_end) {
  mu <- exp(theta[1])
  alpha <- exp(theta[2])
  beta <- exp(theta[3])

  if (alpha >= beta) {
    return(rep(NA_real_, 3))
  }

  n <- length(events)

  if (n == 0) {
    return(c(
      log_mu = -mu * T_end,
      log_alpha = 0,
      log_beta = 0
    ))
  }

  R <- numeric(n)
  dR_dbeta <- numeric(n)

  grad_mu <- 0
  grad_alpha <- 0
  grad_beta <- 0

  for (i in seq_len(n)) {
    if (i == 1) {
      R[i] <- 0
      dR_dbeta[i] <- 0
    } else {
      delta <- events[i] - events[i - 1]
      decay <- exp(-beta * delta)

      R[i] <- decay * (1 + R[i - 1])

      dR_dbeta[i] <- decay * dR_dbeta[i - 1] -
        delta * decay * (1 + R[i - 1])
    }

    lambda_i <- mu + alpha * R[i]

    grad_mu <- grad_mu + 1 / lambda_i
    grad_alpha <- grad_alpha + R[i] / lambda_i
    grad_beta <- grad_beta +
      alpha * dR_dbeta[i] / lambda_i
  }

  tail_times <- T_end - events

  grad_mu <- grad_mu - T_end

  grad_alpha <- grad_alpha -
    sum((1 / beta) * (1 - exp(-beta * tail_times)))

  grad_beta <- grad_beta -
    sum(
      alpha *
        (
          -1 / beta^2 * (1 - exp(-beta * tail_times)) +
            (tail_times / beta) * exp(-beta * tail_times)
        )
    )

  c(
    log_mu = mu * grad_mu,
    log_alpha = alpha * grad_alpha,
    log_beta = beta * grad_beta
  )
}

hawkes_negscore_theta <- function(theta, events, T_end) {
  score <- hawkes_score_theta(theta, events, T_end)

  if (any(!is.finite(score))) {
    return(rep(1e6, 3))
  }

  -score
}

compute_hawkes_event_quantities <- function(events, mu, alpha, beta) {
  n <- length(events)

  R <- numeric(n)
  dR_dbeta <- numeric(n)
  lambda <- numeric(n)

  grad_lambda <- matrix(0, nrow = n, ncol = 3)
  colnames(grad_lambda) <- c("log_mu", "log_alpha", "log_beta")

  for (i in seq_len(n)) {
    if (i == 1) {
      R[i] <- 0
      dR_dbeta[i] <- 0
    } else {
      delta <- events[i] - events[i - 1]
      decay <- exp(-beta * delta)

      R[i] <- decay * (1 + R[i - 1])

      dR_dbeta[i] <- decay * dR_dbeta[i - 1] -
        delta * decay * (1 + R[i - 1])
    }

    lambda[i] <- mu + alpha * R[i]

    grad_lambda[i, "log_mu"] <- mu
    grad_lambda[i, "log_alpha"] <- alpha * R[i]
    grad_lambda[i, "log_beta"] <- beta * alpha * dR_dbeta[i]
  }

  list(
    R = R,
    dR_dbeta = dR_dbeta,
    lambda = lambda,
    grad_lambda = grad_lambda
  )
}

finite_difference_gradient <- function(f, theta, eps = 1e-5) {
  p <- length(theta)
  grad <- numeric(p)

  for (j in seq_len(p)) {
    e <- rep(0, p)
    e[j] <- eps

    grad[j] <- (
      f(theta + e) -
        f(theta - e)
    ) / (2 * eps)
  }

  grad
}

numerical_hessian <- function(f, theta, eps = 1e-4) {
  p <- length(theta)
  H <- matrix(0, nrow = p, ncol = p)

  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      ei <- rep(0, p)
      ej <- rep(0, p)

      ei[i] <- eps
      ej[j] <- eps

      H[i, j] <- (
        f(theta + ei + ej) -
          f(theta + ei - ej) -
          f(theta - ei + ej) +
          f(theta - ei - ej)
      ) / (4 * eps^2)
    }
  }

  if (!is.null(names(theta))) {
    rownames(H) <- names(theta)
    colnames(H) <- names(theta)
  }

  H
}
