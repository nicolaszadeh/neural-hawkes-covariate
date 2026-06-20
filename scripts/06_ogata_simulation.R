# 06_ogata_simulation.R
#
# One-dimensional exponential Hawkes process.
#
# Model:
#
#   lambda(t) = mu + sum_{t_k < t} alpha exp(-beta (t - t_k))
#
# Ogata thinning simulation + likelihood + MLE check.
#
# Stability condition:
#
#   alpha / beta < 1


simulate_hawkes_ogata_1d <- function(mu,
                                     alpha,
                                     beta,
                                     T,
                                     max_events = 1e6) {
  if (mu <= 0) {
    stop("mu must be positive.")
  }
  if (alpha < 0) {
    stop("alpha must be non-negative.")
  }
  if (beta <= 0) {
    stop("beta must be positive.")
  }
  if (T <= 0) {
    stop("T must be positive.")
  }
  
  t <- 0
  excitation <- 0
  events <- numeric(0)
  
  while (t < T && length(events) < max_events) {
    lambda_bar <- mu + excitation
    
    if (lambda_bar <= 0) {
      break
    }
    
    w <- rexp(1, rate = lambda_bar)
    t_candidate <- t + w
    
    if (t_candidate > T) {
      break
    }
    
    excitation_candidate <- excitation * exp(-beta * w)
    lambda_candidate <- mu + excitation_candidate
    
    accept_prob <- lambda_candidate / lambda_bar
    
    if (runif(1) <= accept_prob) {
      events <- c(events, t_candidate)
      excitation <- excitation_candidate + alpha
    } else {
      excitation <- excitation_candidate
    }
    
    t <- t_candidate
  }
  
  if (length(events) >= max_events) {
    warning("max_events reached. The process may be unstable.")
  }
  
  events
}


hawkes_intensity_at_events_1d <- function(events,
                                          mu,
                                          alpha,
                                          beta) {
  n <- length(events)
  
  if (n == 0) {
    return(numeric(0))
  }
  
  intensities <- numeric(n)
  excitation <- 0
  last_t <- 0
  
  for (k in seq_len(n)) {
    dt <- events[k] - last_t
    excitation <- excitation * exp(-beta * dt)
    
    intensities[k] <- mu + excitation
    
    excitation <- excitation + alpha
    last_t <- events[k]
  }
  
  intensities
}


hawkes_loglik_1d <- function(events,
                             T,
                             mu,
                             alpha,
                             beta) {
  if (mu <= 0 || alpha < 0 || beta <= 0 || T <= 0) {
    return(-Inf)
  }
  
  if (length(events) > 0) {
    if (any(events <= 0) || any(events > T)) {
      return(-Inf)
    }
    
    if (is.unsorted(events)) {
      events <- sort(events)
    }
  }
  
  intensities <- hawkes_intensity_at_events_1d(
    events = events,
    mu = mu,
    alpha = alpha,
    beta = beta
  )
  
  if (any(intensities <= 0)) {
    return(-Inf)
  }
  
  log_part <- sum(log(intensities))
  
  compensator <- mu * T
  
  if (length(events) > 0 && alpha > 0) {
    compensator <- compensator +
      sum(alpha / beta * (1 - exp(-beta * (T - events))))
  }
  
  log_part - compensator
}


theta_to_params_1d <- function(theta) {
  mu <- exp(theta[1])
  alpha <- exp(theta[2])
  beta <- alpha + exp(theta[3])
  
  c(mu = mu, alpha = alpha, beta = beta)
}


hawkes_negloglik_theta_1d <- function(theta,
                                      events,
                                      T) {
  pars <- theta_to_params_1d(theta)
  
  ll <- hawkes_loglik_1d(
    events = events,
    T = T,
    mu = pars["mu"],
    alpha = pars["alpha"],
    beta = pars["beta"]
  )
  
  if (!is.finite(ll)) {
    return(1e100)
  }
  
  -ll
}


fit_hawkes_mle_1d <- function(events,
                              T,
                              init_mu = NULL,
                              init_alpha = 0.05,
                              init_beta = 1.00) {
  if (length(events) == 0) {
    stop("Cannot fit Hawkes model with zero events.")
  }
  
  if (is.null(init_mu)) {
    init_mu <- max(length(events) / T, 1e-3)
  }
  
  init_alpha <- max(init_alpha, 1e-6)
  
  if (init_beta <= init_alpha) {
    init_beta <- init_alpha + 1
  }
  
  theta0 <- c(
    log(init_mu),
    log(init_alpha),
    log(init_beta - init_alpha)
  )
  
  opt <- optim(
    par = theta0,
    fn = hawkes_negloglik_theta_1d,
    events = events,
    T = T,
    method = "BFGS",
    control = list(maxit = 2000)
  )
  
  pars <- theta_to_params_1d(opt$par)
  
  list(
    par = pars,
    branching_ratio = pars["alpha"] / pars["beta"],
    loglik = -opt$value,
    convergence = opt$convergence,
    optim = opt
  )
}


summarize_fit_1d <- function(label,
                             true_mu,
                             true_alpha,
                             true_beta,
                             events,
                             fit) {
  true_eta <- true_alpha / true_beta
  est <- fit$par
  
  data.frame(
    case = label,
    n_events = length(events),
    
    true_mu = true_mu,
    est_mu = unname(est["mu"]),
    
    true_alpha = true_alpha,
    est_alpha = unname(est["alpha"]),
    
    true_beta = true_beta,
    est_beta = unname(est["beta"]),
    
    true_branching = true_eta,
    est_branching = unname(fit$branching_ratio),
    
    loglik = fit$loglik,
    convergence = fit$convergence
  )
}


run_ogata_demo <- function() {
  set.seed(123)
  
  T <- 500
  
  negative_true <- list(
    mu = 0.50,
    alpha = 0.00,
    beta = 1.50
  )
  
  positive_true <- list(
    mu = 0.50,
    alpha = 0.60,
    beta = 1.50
  )
  
  events_neg <- simulate_hawkes_ogata_1d(
    mu = negative_true$mu,
    alpha = negative_true$alpha,
    beta = negative_true$beta,
    T = T
  )
  
  events_pos <- simulate_hawkes_ogata_1d(
    mu = positive_true$mu,
    alpha = positive_true$alpha,
    beta = positive_true$beta,
    T = T
  )
  
  fit_neg <- fit_hawkes_mle_1d(
    events = events_neg,
    T = T
  )
  
  fit_pos <- fit_hawkes_mle_1d(
    events = events_pos,
    T = T
  )
  
  res_neg <- summarize_fit_1d(
    label = "negative_control_poisson",
    true_mu = negative_true$mu,
    true_alpha = negative_true$alpha,
    true_beta = negative_true$beta,
    events = events_neg,
    fit = fit_neg
  )
  
  res_pos <- summarize_fit_1d(
    label = "positive_control_hawkes",
    true_mu = positive_true$mu,
    true_alpha = positive_true$alpha,
    true_beta = positive_true$beta,
    events = events_pos,
    fit = fit_pos
  )
  
  results <- rbind(res_neg, res_pos)
  
  print(results)
  
  invisible(list(
    negative_events = events_neg,
    positive_events = events_pos,
    negative_fit = fit_neg,
    positive_fit = fit_pos,
    table = results
  ))
}


demo <- run_ogata_demo()