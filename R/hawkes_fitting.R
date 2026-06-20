# ============================================================
# Fitting and likelihood-ratio tests for binned Hawkes models
# ============================================================


# ============================================================
# Small utility: empirical event rate
# ============================================================

compute_event_rate <- function(dN, dt) {
  T_obs <- max((length(dN) - 1) * dt, dt)
  sum(dN) / T_obs
}


# ============================================================
# Small utility: remove accidental row names
# ============================================================

clean_output_table <- function(df) {
  rownames(df) <- NULL
  df
}


# ============================================================
# Check that likelihood functions have been loaded
# ============================================================

check_hawkes_likelihood_loaded <- function() {
  needed <- c(
    "neg_loglik_full",
    "neg_loglik_no_covariate",
    "neg_loglik_no_hawkes"
  )
  
  missing <- needed[!vapply(needed, exists, logical(1))]
  
  if (length(missing) > 0) {
    stop(
      "Missing likelihood functions: ",
      paste(missing, collapse = ", "),
      "\nPlease run source('R/hawkes_likelihood.R') first."
    )
  }
  
  invisible(TRUE)
}


# ============================================================
# Fit the full model
# ============================================================

fit_full_model <- function(
    X,
    dN,
    dt,
    hessian = FALSE,
    maxit = 1000
) {
  check_hawkes_likelihood_loaded()
  
  event_rate <- compute_event_rate(dN, dt)
  
  start_full <- c(
    mu0 = log(event_rate + 1e-6),
    mu1 = 0,
    log_alpha = log(0.2),
    log_beta = log(1.0)
  )
  
  fit <- optim(
    par = start_full,
    fn = neg_loglik_full,
    X = X,
    dN = dN,
    dt = dt,
    method = "BFGS",
    hessian = hessian,
    control = list(maxit = maxit)
  )
  
  est <- c(
    mu0 = unname(fit$par[1]),
    mu1 = unname(fit$par[2]),
    alpha = unname(exp(fit$par[3])),
    beta = unname(exp(fit$par[4]))
  )
  
  list(
    fit = fit,
    est = est,
    loglik = -fit$value
  )
}


# ============================================================
# Fit the no-covariate null model
# ============================================================

fit_no_covariate_model <- function(
    X,
    dN,
    dt,
    hessian = FALSE,
    maxit = 1000
) {
  check_hawkes_likelihood_loaded()
  
  event_rate <- compute_event_rate(dN, dt)
  
  start_no_cov <- c(
    mu0 = log(event_rate + 1e-6),
    log_alpha = log(0.2),
    log_beta = log(1.0)
  )
  
  fit <- optim(
    par = start_no_cov,
    fn = neg_loglik_no_covariate,
    X = X,
    dN = dN,
    dt = dt,
    method = "BFGS",
    hessian = hessian,
    control = list(maxit = maxit)
  )
  
  est <- c(
    mu0 = unname(fit$par[1]),
    mu1 = 0,
    alpha = unname(exp(fit$par[2])),
    beta = unname(exp(fit$par[3]))
  )
  
  list(
    fit = fit,
    est = est,
    loglik = -fit$value
  )
}


# ============================================================
# Fit the no-Hawkes null model
# ============================================================

fit_no_hawkes_model <- function(
    X,
    dN,
    dt,
    hessian = FALSE,
    maxit = 1000
) {
  check_hawkes_likelihood_loaded()
  
  event_rate <- compute_event_rate(dN, dt)
  
  start_no_hawkes <- c(
    mu0 = log(event_rate + 1e-6),
    mu1 = 0
  )
  
  fit <- optim(
    par = start_no_hawkes,
    fn = neg_loglik_no_hawkes,
    X = X,
    dN = dN,
    dt = dt,
    method = "BFGS",
    hessian = hessian,
    control = list(maxit = maxit)
  )
  
  est <- c(
    mu0 = unname(fit$par[1]),
    mu1 = unname(fit$par[2]),
    alpha = 0,
    beta = NA
  )
  
  list(
    fit = fit,
    est = est,
    loglik = -fit$value
  )
}


# ============================================================
# Fit all three models and compute both LR tests
# ============================================================

fit_all_models_and_tests <- function(
    X,
    dN,
    dt,
    hessian = FALSE,
    maxit = 1000
) {
  full <- fit_full_model(
    X = X,
    dN = dN,
    dt = dt,
    hessian = hessian,
    maxit = maxit
  )
  
  no_cov <- fit_no_covariate_model(
    X = X,
    dN = dN,
    dt = dt,
    hessian = hessian,
    maxit = maxit
  )
  
  no_hawkes <- fit_no_hawkes_model(
    X = X,
    dN = dN,
    dt = dt,
    hessian = hessian,
    maxit = maxit
  )
  
  LR_cov <- 2 * (full$loglik - no_cov$loglik)
  LR_cov <- max(LR_cov, 0)
  
  p_cov <- pchisq(
    LR_cov,
    df = 1,
    lower.tail = FALSE
  )
  
  LR_hawkes <- 2 * (full$loglik - no_hawkes$loglik)
  LR_hawkes <- max(LR_hawkes, 0)
  
  p_hawkes <- pchisq(
    LR_hawkes,
    df = 1,
    lower.tail = FALSE
  )
  
  out <- data.frame(
    n_events = sum(dN),
    max_dN = max(dN),
    
    mu0_hat = unname(full$est["mu0"]),
    mu1_hat = unname(full$est["mu1"]),
    alpha_hat = unname(full$est["alpha"]),
    beta_hat = unname(full$est["beta"]),
    
    loglik_full = full$loglik,
    loglik_no_cov = no_cov$loglik,
    loglik_no_hawkes = no_hawkes$loglik,
    
    LR_cov = LR_cov,
    p_cov = p_cov,
    
    LR_hawkes = LR_hawkes,
    p_hawkes = p_hawkes,
    
    conv_full = full$fit$convergence,
    conv_no_cov = no_cov$fit$convergence,
    conv_no_hawkes = no_hawkes$fit$convergence
  )
  
  clean_output_table(out)
}


# ============================================================
# Fit only the Hawkes LR test
# ============================================================

fit_hawkes_test <- function(
    X,
    dN,
    dt,
    hessian = FALSE,
    maxit = 1000
) {
  full <- fit_full_model(
    X = X,
    dN = dN,
    dt = dt,
    hessian = hessian,
    maxit = maxit
  )
  
  no_hawkes <- fit_no_hawkes_model(
    X = X,
    dN = dN,
    dt = dt,
    hessian = hessian,
    maxit = maxit
  )
  
  LR <- 2 * (full$loglik - no_hawkes$loglik)
  LR <- max(LR, 0)
  
  p_chisq <- pchisq(
    LR,
    df = 1,
    lower.tail = FALSE
  )
  
  out <- data.frame(
    n_events = sum(dN),
    
    mu0_null = unname(no_hawkes$est["mu0"]),
    mu1_null = unname(no_hawkes$est["mu1"]),
    
    mu0_full = unname(full$est["mu0"]),
    mu1_full = unname(full$est["mu1"]),
    alpha_full = unname(full$est["alpha"]),
    beta_full = unname(full$est["beta"]),
    
    loglik_null = no_hawkes$loglik,
    loglik_full = full$loglik,
    
    LR = LR,
    p_chisq = p_chisq,
    
    conv_null = no_hawkes$fit$convergence,
    conv_full = full$fit$convergence
  )
  
  clean_output_table(out)
}


# ============================================================
# Fit only the covariate LR test
# ============================================================

fit_covariate_test <- function(
    X,
    dN,
    dt,
    hessian = FALSE,
    maxit = 1000
) {
  full <- fit_full_model(
    X = X,
    dN = dN,
    dt = dt,
    hessian = hessian,
    maxit = maxit
  )
  
  no_cov <- fit_no_covariate_model(
    X = X,
    dN = dN,
    dt = dt,
    hessian = hessian,
    maxit = maxit
  )
  
  LR <- 2 * (full$loglik - no_cov$loglik)
  LR <- max(LR, 0)
  
  p_chisq <- pchisq(
    LR,
    df = 1,
    lower.tail = FALSE
  )
  
  out <- data.frame(
    n_events = sum(dN),
    
    mu0_null = unname(no_cov$est["mu0"]),
    alpha_null = unname(no_cov$est["alpha"]),
    beta_null = unname(no_cov$est["beta"]),
    
    mu0_full = unname(full$est["mu0"]),
    mu1_full = unname(full$est["mu1"]),
    alpha_full = unname(full$est["alpha"]),
    beta_full = unname(full$est["beta"]),
    
    loglik_null = no_cov$loglik,
    loglik_full = full$loglik,
    
    LR = LR,
    p_chisq = p_chisq,
    
    conv_null = no_cov$fit$convergence,
    conv_full = full$fit$convergence
  )
  
  clean_output_table(out)
}