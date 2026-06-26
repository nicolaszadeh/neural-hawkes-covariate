# ============================================================
# SPBDBS2025-style Algorithm 2 helpers
# ============================================================

hawkes_rho_hat <- function(events, T_end, mu, alpha, beta) {
  tail_times <- T_end - events

  integral_R <- sum(
    (1 / beta) *
      (1 - exp(-beta * tail_times))
  )

  integral_dR_dbeta <- sum(
    -(
      1 / beta^2 *
        (1 - exp(-beta * tail_times)) -
        (tail_times / beta) *
        exp(-beta * tail_times)
    )
  )

  c(
    log_mu = mu,
    log_alpha = alpha * integral_R / T_end,
    log_beta = beta * alpha * integral_dR_dbeta / T_end
  )
}

hawkes_I_hat_algorithm2 <- function(events, T_end, mu, alpha, beta) {
  eq <- compute_hawkes_event_quantities(
    events = events,
    mu = mu,
    alpha = alpha,
    beta = beta
  )

  grad_lambda <- eq$grad_lambda
  lambda_events <- eq$lambda

  I_hat <- matrix(0, nrow = 3, ncol = 3)

  for (i in seq_along(events)) {
    g <- matrix(grad_lambda[i, ], ncol = 1)
    I_hat <- I_hat + (g %*% t(g)) / (lambda_events[i]^2)
  }

  I_hat <- I_hat / T_end

  rownames(I_hat) <- colnames(grad_lambda)
  colnames(I_hat) <- colnames(grad_lambda)

  I_hat
}

hawkes_algorithm2_corrected_increments <- function(
    events,
    T_end,
    mu,
    alpha,
    beta,
    x_gaussian = rnorm(3)
) {
  Lambda_events <- hawkes_compensator_at_events(
    events = events,
    mu = mu,
    alpha = alpha,
    beta = beta
  )

  delta_Lambda <- diff(c(0, Lambda_events))
  delta_t <- diff(c(0, events))

  rho_hat <- hawkes_rho_hat(
    events = events,
    T_end = T_end,
    mu = mu,
    alpha = alpha,
    beta = beta
  )

  I_hat <- hawkes_I_hat_algorithm2(
    events = events,
    T_end = T_end,
    mu = mu,
    alpha = alpha,
    beta = beta
  )

  I_inv_sqrt <- matrix_inv_sqrt(I_hat)

  correction_scalar <- as.numeric(
    rho_hat %*% I_inv_sqrt %*% x_gaussian
  )

  corrected <- delta_Lambda +
    (delta_t / sqrt(T_end)) * correction_scalar

  list(
    usual = delta_Lambda,
    corrected = corrected,
    positive_corrected = corrected[corrected > 0],
    n_nonpositive = sum(corrected <= 0),
    rho_hat = rho_hat,
    I_hat = I_hat,
    I_inv_sqrt = I_inv_sqrt,
    correction_scalar = correction_scalar
  )
}

hawkes_algorithm2_test <- function(
    events,
    T_end,
    mu,
    alpha,
    beta,
    alpha_level = 0.05
) {
  inc <- hawkes_algorithm2_corrected_increments(
    events = events,
    T_end = T_end,
    mu = mu,
    alpha = alpha,
    beta = beta
  )

  N_events <- length(events)
  m_subsample <- floor(N_events^(2 / 3))

  subsample_indices <- sample(
    seq_along(inc$positive_corrected),
    size = min(m_subsample, length(inc$positive_corrected)),
    replace = FALSE
  )

  subsampled_corrected <- inc$positive_corrected[subsample_indices]

  ks_usual <- ks.test(inc$usual, "pexp", rate = 1)
  ks_corrected <- ks.test(inc$positive_corrected, "pexp", rate = 1)
  ks_subsampled <- ks.test(subsampled_corrected, "pexp", rate = 1)

  data.frame(
    n_events = N_events,
    m_subsample = m_subsample,
    n_nonpositive_corrected = inc$n_nonpositive,
    correction_scalar = inc$correction_scalar,
    usual_D = as.numeric(ks_usual$statistic),
    usual_p = ks_usual$p.value,
    corrected_D = as.numeric(ks_corrected$statistic),
    corrected_p = ks_corrected$p.value,
    subsampled_D = as.numeric(ks_subsampled$statistic),
    subsampled_p = ks_subsampled$p.value,
    usual_reject = ks_usual$p.value < alpha_level,
    corrected_reject = ks_corrected$p.value < alpha_level,
    subsampled_reject = ks_subsampled$p.value < alpha_level
  )
}
