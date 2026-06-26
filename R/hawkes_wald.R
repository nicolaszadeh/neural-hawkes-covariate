# ============================================================
# Wald tests
# ============================================================

wald_one_parameter <- function(
    theta_hat,
    I_inv_hat,
    theta0,
    parameter,
    T_end,
    alpha_level = 0.05
) {
  idx <- match(parameter, names(theta_hat))

  if (is.na(idx)) {
    stop("Unknown parameter: ", parameter)
  }

  if (length(theta0) == length(theta_hat)) {
    theta0_value <- theta0[idx]
  } else {
    theta0_value <- theta0
  }

  sigma_hat <- sqrt(I_inv_hat[idx, idx])

  Z <- sqrt(T_end) *
    (theta_hat[idx] - theta0_value) /
    sigma_hat

  p_value <- 2 * (1 - pnorm(abs(Z)))

  out <- data.frame(
    parameter = parameter,
    theta_hat = as.numeric(theta_hat[idx]),
    theta0 = as.numeric(theta0_value),
    sigma_hat = sigma_hat,
    Z = as.numeric(Z),
    p_value = p_value,
    reject_5_percent = p_value < alpha_level
  )

  rownames(out) <- NULL
  out
}

wald_all_true_parameters <- function(
    theta_hat,
    I_inv_hat,
    theta_true,
    T_end,
    alpha_level = 0.05
) {
  out <- lapply(names(theta_true), function(parameter) {
    wald_one_parameter(
      theta_hat = theta_hat,
      I_inv_hat = I_inv_hat,
      theta0 = theta_true,
      parameter = parameter,
      T_end = T_end,
      alpha_level = alpha_level
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

wald_equality <- function(
    theta_hat,
    I_inv_hat,
    parameter_i,
    parameter_j,
    T_end,
    alpha_level = 0.05
) {
  i <- match(parameter_i, names(theta_hat))
  j <- match(parameter_j, names(theta_hat))

  if (is.na(i) || is.na(j)) {
    stop("Unknown parameter in equality test.")
  }

  contrast_hat <- theta_hat[i] - theta_hat[j]

  var_contrast <- I_inv_hat[i, i] +
    I_inv_hat[j, j] -
    2 * I_inv_hat[i, j]

  sigma_hat <- sqrt(var_contrast)

  Z <- sqrt(T_end) * contrast_hat / sigma_hat
  p_value <- 2 * (1 - pnorm(abs(Z)))

  out <- data.frame(
    parameter_i = parameter_i,
    parameter_j = parameter_j,
    theta_i_hat = as.numeric(theta_hat[i]),
    theta_j_hat = as.numeric(theta_hat[j]),
    contrast_hat = as.numeric(contrast_hat),
    sigma_hat = sigma_hat,
    Z = as.numeric(Z),
    p_value = p_value,
    reject_5_percent = p_value < alpha_level
  )

  rownames(out) <- NULL
  out
}
