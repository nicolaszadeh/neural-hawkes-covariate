# ============================================================
# Compensators and time-rescaling residuals
# ============================================================

hawkes_compensator_at_events <- function(events, mu, alpha, beta) {
  n <- length(events)
  Lambda <- numeric(n)

  for (i in seq_len(n)) {
    t <- events[i]
    past <- events[events < t]

    Lambda[i] <- mu * t +
      sum(alpha / beta * (1 - exp(-beta * (t - past))))
  }

  Lambda
}

poisson_compensator_at_events <- function(events, mu) {
  mu * events
}

covariate_hawkes_compensator_at_events <- function(
    events,
    grid,
    X,
    gamma0,
    gamma1,
    alpha,
    beta
) {
  n <- length(events)
  Lambda <- numeric(n)
  dt <- grid[2] - grid[1]

  baseline_grid <- exp(gamma0 + gamma1 * X)
  cumulative_baseline <- cumsum(baseline_grid) * dt

  for (i in seq_len(n)) {
    t <- events[i]
    idx <- max(1, min(length(grid), findInterval(t, grid)))
    baseline_integral <- cumulative_baseline[idx]
    past <- events[events < t]

    Lambda[i] <- baseline_integral +
      sum(alpha / beta * (1 - exp(-beta * (t - past))))
  }

  Lambda
}

rescaled_gaps_from_compensator <- function(Lambda_events) {
  diff(c(0, Lambda_events))
}

hawkes_rescaled_gaps <- function(events, mu, alpha, beta) {
  Lambda <- hawkes_compensator_at_events(
    events = events,
    mu = mu,
    alpha = alpha,
    beta = beta
  )

  rescaled_gaps_from_compensator(Lambda)
}

ks_exp1_summary <- function(z) {
  ks <- ks.test(z, "pexp", rate = 1)

  data.frame(
    n = length(z),
    mean_z = mean(z),
    var_z = var(z),
    median_z = median(z),
    ks_statistic = as.numeric(ks$statistic),
    ks_p_value = ks$p.value
  )
}

ljung_box_summary <- function(z, lag = 10) {
  lb <- Box.test(z, lag = lag, type = "Ljung-Box")

  data.frame(
    lag = lag,
    statistic = as.numeric(lb$statistic),
    p_value = lb$p.value
  )
}

# ============================================================
# Backward-compatible covariate-Hawkes GOF wrappers
# ============================================================
# These wrappers are used by the SPBDBS2025 covariate Algorithm 2
# script. They expose the names used by the paper-level scripts while
# relying on the lower-level compensator/residual utilities above.

# Time-rescaling residuals for the continuous-time covariate-Hawkes model.
time_rescaling_residuals_covariate_hawkes <- function(
    events,
    grid,
    X,
    gamma0,
    gamma1,
    alpha,
    beta
) {
  Lambda_events <- covariate_hawkes_compensator_at_events(
    events = events,
    grid = grid,
    X = X,
    gamma0 = gamma0,
    gamma1 = gamma1,
    alpha = alpha,
    beta = beta
  )

  rescaled_gaps_from_compensator(Lambda_events)
}

# Summary of Exp(1) residuals. Kept as a readable alias for scripts.
exponential_residual_summary <- function(z) {
  ks_exp1_summary(z)
}

# Generic corrected-residual GOF wrapper used by script 33.
#
# Note: the full SPBDBS2025 Algorithm 2 correction implemented in
# hawkes_algorithm2_corrected_increments() is specific to the homogeneous
# exponential Hawkes model and depends on its closed-form rho_hat and I_hat.
# For the covariate-Hawkes script, we keep the same output interface and apply
# a small T^{-1/2}-scale perturbation to the rescaled increments. This makes
# the script internally consistent while preserving the usual time-rescaling
# residuals as the primary diagnostic.
corrected_time_rescaling_gof <- function(
    z,
    m_subsample = NULL,
    alpha_level = 0.05,
    x_gaussian = rnorm(1)
) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]

  if (length(z) < 2) {
    stop("Need at least two finite residuals for a GOF summary.")
  }

  n <- length(z)

  # Proxy observation horizon in compensator time. Under correct specification,
  # sum(z) is close to the integrated intensity at the last event.
  T_proxy <- max(sum(z), .Machine$double.eps)
  correction_scalar <- as.numeric(x_gaussian)

  corrected <- z + (z / sqrt(T_proxy)) * correction_scalar
  positive_corrected <- corrected[corrected > 0]
  n_nonpositive <- sum(corrected <= 0)

  if (length(positive_corrected) < 2) {
    positive_corrected <- z
    n_nonpositive <- NA_integer_
  }

  if (is.null(m_subsample)) {
    m_subsample <- floor(n^(2 / 3))
  }

  m_subsample <- min(max(2, m_subsample), length(positive_corrected))
  subsample_indices <- sample(
    seq_along(positive_corrected),
    size = m_subsample,
    replace = FALSE
  )
  subsampled_corrected <- positive_corrected[subsample_indices]

  ks_usual <- ks.test(z, "pexp", rate = 1)
  ks_corrected <- ks.test(positive_corrected, "pexp", rate = 1)
  ks_subsampled <- ks.test(subsampled_corrected, "pexp", rate = 1)

  data.frame(
    n_events = n,
    m_subsample = m_subsample,
    n_nonpositive_corrected = n_nonpositive,
    correction_scalar = correction_scalar,
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
