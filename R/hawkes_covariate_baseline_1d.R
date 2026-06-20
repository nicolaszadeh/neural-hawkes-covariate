# R/hawkes_covariate_baseline_1d.R
#
# One-dimensional exponential Hawkes process with covariate-driven
# baseline:
#
#   lambda(t) = exp(gamma0 + gamma1 X(t))
#             + sum_{t_k < t} alpha exp(-beta (t - t_k))
#
# Here X(t) is deterministic and observed.
#
# This is the next step toward diffusion-driven baselines.


source("R/hawkes_td_baseline_1d.R")


make_sinusoidal_covariate_1d <- function(period = 100,
                                         phase = 0) {
  if (period <= 0) {
    stop("period must be positive.")
  }
  
  x_fun <- function(t) {
    sin(2 * pi * t / period + phase)
  }
  
  list(
    x_fun = x_fun,
    x_lower = -1,
    x_upper = 1,
    period = period,
    phase = phase
  )
}


make_covariate_baseline_1d <- function(x_fun,
                                       gamma0,
                                       gamma1,
                                       x_lower = -1,
                                       x_upper = 1) {
  mu_fun <- function(t) {
    exp(gamma0 + gamma1 * x_fun(t))
  }
  
  if (gamma1 >= 0) {
    mu_upper <- exp(gamma0 + gamma1 * x_upper)
  } else {
    mu_upper <- exp(gamma0 + gamma1 * x_lower)
  }
  
  list(
    mu_fun = mu_fun,
    mu_upper = mu_upper,
    x_fun = x_fun,
    gamma0 = gamma0,
    gamma1 = gamma1,
    x_lower = x_lower,
    x_upper = x_upper
  )
}


covariate_baseline_integral_1d <- function(x_fun,
                                           T,
                                           gamma0,
                                           gamma1) {
  integrate(
    f = function(t) {
      exp(gamma0 + gamma1 * x_fun(t))
    },
    lower = 0,
    upper = T,
    subdivisions = 1000,
    rel.tol = 1e-8
  )$value
}


simulate_covariate_hawkes_1d <- function(x_fun,
                                         gamma0,
                                         gamma1,
                                         x_lower,
                                         x_upper,
                                         alpha,
                                         beta,
                                         T) {
  baseline <- make_covariate_baseline_1d(
    x_fun = x_fun,
    gamma0 = gamma0,
    gamma1 = gamma1,
    x_lower = x_lower,
    x_upper = x_upper
  )
  
  simulate_hawkes_ogata_td_baseline_1d(
    mu_fun = baseline$mu_fun,
    mu_upper = baseline$mu_upper,
    alpha = alpha,
    beta = beta,
    T = T
  )
}


covariate_poisson_loglik_1d <- function(events,
                                        T,
                                        x_fun,
                                        gamma0,
                                        gamma1,
                                        baseline_integral = NULL) {
  if (is.null(baseline_integral)) {
    baseline_integral <- covariate_baseline_integral_1d(
      x_fun = x_fun,
      T = T,
      gamma0 = gamma0,
      gamma1 = gamma1
    )
  }
  
  if (length(events) == 0) {
    return(-baseline_integral)
  }
  
  if (any(events <= 0) || any(events > T)) {
    return(-Inf)
  }
  
  log_part <- sum(gamma0 + gamma1 * x_fun(events))
  
  log_part - baseline_integral
}


covariate_poisson_negloglik_theta_1d <- function(theta,
                                                 events,
                                                 T,
                                                 x_fun) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]
  
  ll <- covariate_poisson_loglik_1d(
    events = events,
    T = T,
    x_fun = x_fun,
    gamma0 = gamma0,
    gamma1 = gamma1
  )
  
  if (!is.finite(ll)) {
    return(1e100)
  }
  
  -ll
}


fit_covariate_poisson_mle_1d <- function(events,
                                         T,
                                         x_fun,
                                         use_multistart = TRUE) {
  n <- length(events)
  
  init_gamma0 <- log(max(n / T, 1e-6))
  
  starts <- data.frame(
    gamma0 = init_gamma0,
    gamma1 = 0
  )
  
  if (use_multistart) {
    extra <- data.frame(
      gamma0 = c(init_gamma0, init_gamma0, init_gamma0),
      gamma1 = c(-1, 0.5, 1)
    )
    
    starts <- rbind(starts, extra)
  }
  
  best <- NULL
  
  for (s in seq_len(nrow(starts))) {
    theta0 <- c(
      starts$gamma0[s],
      starts$gamma1[s]
    )
    
    opt <- tryCatch({
      optim(
        par = theta0,
        fn = covariate_poisson_negloglik_theta_1d,
        events = events,
        T = T,
        x_fun = x_fun,
        method = "BFGS",
        control = list(maxit = 2000)
      )
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(opt) && is.finite(opt$value)) {
      if (is.null(best) || opt$value < best$value) {
        best <- opt
      }
    }
  }
  
  if (is.null(best)) {
    stop("Poisson covariate MLE failed.")
  }
  
  list(
    par = c(
      gamma0 = best$par[1],
      gamma1 = best$par[2]
    ),
    loglik = -best$value,
    convergence = best$convergence,
    optim = best
  )
}


theta_to_cov_hawkes_params_1d <- function(theta) {
  gamma0 <- theta[1]
  gamma1 <- theta[2]
  
  alpha <- exp(theta[3])
  beta <- alpha + exp(theta[4])
  
  c(
    gamma0 = gamma0,
    gamma1 = gamma1,
    alpha = alpha,
    beta = beta
  )
}


cov_hawkes_params_to_theta_1d <- function(gamma0,
                                          gamma1,
                                          alpha,
                                          beta) {
  alpha <- max(alpha, 1e-8)
  
  if (beta <= alpha) {
    beta <- alpha + 1
  }
  
  c(
    gamma0,
    gamma1,
    log(alpha),
    log(beta - alpha)
  )
}


covariate_hawkes_loglik_1d <- function(events,
                                       T,
                                       x_fun,
                                       gamma0,
                                       gamma1,
                                       alpha,
                                       beta,
                                       baseline_integral = NULL) {
  if (alpha < 0 || beta <= 0 || T <= 0) {
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
  
  n <- length(events)
  
  excitation <- 0
  last_t <- 0
  log_part <- 0
  
  if (n > 0) {
    for (k in seq_len(n)) {
      dt <- events[k] - last_t
      excitation <- excitation * exp(-beta * dt)
      
      mu_k <- exp(gamma0 + gamma1 * x_fun(events[k]))
      lambda_k <- mu_k + excitation
      
      if (lambda_k <= 0) {
        return(-Inf)
      }
      
      log_part <- log_part + log(lambda_k)
      
      excitation <- excitation + alpha
      last_t <- events[k]
    }
  }
  
  if (is.null(baseline_integral)) {
    baseline_integral <- covariate_baseline_integral_1d(
      x_fun = x_fun,
      T = T,
      gamma0 = gamma0,
      gamma1 = gamma1
    )
  }
  
  compensator <- baseline_integral
  
  if (n > 0 && alpha > 0) {
    compensator <- compensator +
      sum(alpha / beta * (1 - exp(-beta * (T - events))))
  }
  
  log_part - compensator
}


covariate_hawkes_negloglik_theta_1d <- function(theta,
                                                events,
                                                T,
                                                x_fun) {
  pars <- theta_to_cov_hawkes_params_1d(theta)
  
  ll <- covariate_hawkes_loglik_1d(
    events = events,
    T = T,
    x_fun = x_fun,
    gamma0 = pars["gamma0"],
    gamma1 = pars["gamma1"],
    alpha = pars["alpha"],
    beta = pars["beta"]
  )
  
  if (!is.finite(ll)) {
    return(1e100)
  }
  
  -ll
}


fit_covariate_hawkes_mle_1d <- function(events,
                                        T,
                                        x_fun,
                                        poisson_fit = NULL,
                                        use_multistart = TRUE) {
  if (length(events) == 0) {
    stop("Cannot fit Hawkes model with zero events.")
  }
  
  if (is.null(poisson_fit)) {
    poisson_fit <- fit_covariate_poisson_mle_1d(
      events = events,
      T = T,
      x_fun = x_fun
    )
  }
  
  gamma0_0 <- unname(poisson_fit$par["gamma0"])
  gamma1_0 <- unname(poisson_fit$par["gamma1"])
  
  starts <- data.frame(
    gamma0 = gamma0_0,
    gamma1 = gamma1_0,
    alpha = 0.05,
    beta = 1.00
  )
  
  if (use_multistart) {
    extra <- data.frame(
      gamma0 = c(gamma0_0, gamma0_0, gamma0_0, gamma0_0),
      gamma1 = c(gamma1_0, gamma1_0, gamma1_0, gamma1_0),
      alpha = c(0.01, 0.10, 0.30, 0.60),
      beta = c(1.00, 1.50, 2.00, 2.50)
    )
    
    starts <- rbind(starts, extra)
  }
  
  best <- NULL
  
  for (s in seq_len(nrow(starts))) {
    theta0 <- cov_hawkes_params_to_theta_1d(
      gamma0 = starts$gamma0[s],
      gamma1 = starts$gamma1[s],
      alpha = starts$alpha[s],
      beta = starts$beta[s]
    )
    
    opt <- tryCatch({
      optim(
        par = theta0,
        fn = covariate_hawkes_negloglik_theta_1d,
        events = events,
        T = T,
        x_fun = x_fun,
        method = "BFGS",
        control = list(maxit = 3000)
      )
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(opt) && is.finite(opt$value)) {
      if (is.null(best) || opt$value < best$value) {
        best <- opt
      }
    }
  }
  
  if (is.null(best)) {
    stop("Covariate Hawkes MLE failed.")
  }
  
  pars <- theta_to_cov_hawkes_params_1d(best$par)
  
  list(
    par = pars,
    branching_ratio = pars["alpha"] / pars["beta"],
    loglik = -best$value,
    convergence = best$convergence,
    optim = best
  )
}


compute_lrt_covariate_baseline_1d <- function(events,
                                              T,
                                              x_fun) {
  n <- length(events)
  
  if (n == 0) {
    return(list(
      n_events = 0,
      poisson_fit = NULL,
      hawkes_fit = NULL,
      poisson_loglik = NA,
      hawkes_loglik = NA,
      lrt_stat = NA,
      error = TRUE
    ))
  }
  
  out <- tryCatch({
    poisson_fit <- fit_covariate_poisson_mle_1d(
      events = events,
      T = T,
      x_fun = x_fun
    )
    
    hawkes_fit <- fit_covariate_hawkes_mle_1d(
      events = events,
      T = T,
      x_fun = x_fun,
      poisson_fit = poisson_fit
    )
    
    lrt <- 2 * (
      hawkes_fit$loglik - poisson_fit$loglik
    )
    
    list(
      n_events = n,
      poisson_fit = poisson_fit,
      hawkes_fit = hawkes_fit,
      poisson_loglik = poisson_fit$loglik,
      hawkes_loglik = hawkes_fit$loglik,
      lrt_stat = lrt,
      error = FALSE
    )
  }, error = function(e) {
    list(
      n_events = n,
      poisson_fit = NULL,
      hawkes_fit = NULL,
      poisson_loglik = NA,
      hawkes_loglik = NA,
      lrt_stat = NA,
      error = TRUE
    )
  })
  
  out
}


summarize_covariate_baseline_fit_1d <- function(label,
                                                true_gamma0,
                                                true_gamma1,
                                                true_alpha,
                                                true_beta,
                                                events,
                                                test) {
  if (is.null(test$poisson_fit)) {
    null_gamma0 <- NA
    null_gamma1 <- NA
  } else {
    null_gamma0 <- unname(test$poisson_fit$par["gamma0"])
    null_gamma1 <- unname(test$poisson_fit$par["gamma1"])
  }
  
  if (is.null(test$hawkes_fit)) {
    hawkes_gamma0 <- NA
    hawkes_gamma1 <- NA
    est_alpha <- NA
    est_beta <- NA
    est_branching <- NA
    convergence <- NA
  } else {
    hp <- test$hawkes_fit$par
    
    hawkes_gamma0 <- unname(hp["gamma0"])
    hawkes_gamma1 <- unname(hp["gamma1"])
    est_alpha <- unname(hp["alpha"])
    est_beta <- unname(hp["beta"])
    est_branching <- unname(test$hawkes_fit$branching_ratio)
    convergence <- test$hawkes_fit$convergence
  }
  
  data.frame(
    case = label,
    n_events = length(events),
    
    true_gamma0 = true_gamma0,
    poisson_gamma0 = null_gamma0,
    hawkes_gamma0 = hawkes_gamma0,
    
    true_gamma1 = true_gamma1,
    poisson_gamma1 = null_gamma1,
    hawkes_gamma1 = hawkes_gamma1,
    
    true_alpha = true_alpha,
    est_alpha = est_alpha,
    
    true_beta = true_beta,
    est_beta = est_beta,
    
    true_branching = true_alpha / true_beta,
    est_branching = est_branching,
    
    poisson_loglik = test$poisson_loglik,
    hawkes_loglik = test$hawkes_loglik,
    lrt_stat = test$lrt_stat,
    
    convergence = convergence,
    error = test$error
  )
}


plot_covariate_and_baseline_1d <- function(covariate,
                                           baseline,
                                           T,
                                           output_dir) {
  png(
    filename = file.path(
      output_dir,
      "covariate_and_baseline.png"
    ),
    width = 900,
    height = 700
  )
  
  old_par <- par(no.readonly = TRUE)
  
  par(mfrow = c(2, 1))
  
  curve(
    covariate$x_fun(x),
    from = 0,
    to = T,
    xlab = "t",
    ylab = "X(t)",
    main = "Observed covariate"
  )
  
  curve(
    baseline$mu_fun(x),
    from = 0,
    to = T,
    xlab = "t",
    ylab = expression(mu(t)),
    main = expression(mu(t) == exp(gamma[0] + gamma[1] * X(t)))
  )
  
  par(old_par)
  dev.off()
  
  invisible(NULL)
}


run_covariate_baseline_demo_1d <- function(T = 500,
                                           seed = 20260620) {
  set.seed(seed)
  
  output_dir <- "results/covariate_baseline"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  covariate <- make_sinusoidal_covariate_1d(
    period = 100,
    phase = 0
  )
  
  true_gamma0 <- log(0.50)
  true_gamma1 <- 0.80
  
  baseline <- make_covariate_baseline_1d(
    x_fun = covariate$x_fun,
    gamma0 = true_gamma0,
    gamma1 = true_gamma1,
    x_lower = covariate$x_lower,
    x_upper = covariate$x_upper
  )
  
  negative_true <- list(
    alpha = 0.00,
    beta = 1.50
  )
  
  positive_true <- list(
    alpha = 0.60,
    beta = 1.50
  )
  
  events_neg <- simulate_covariate_hawkes_1d(
    x_fun = covariate$x_fun,
    gamma0 = true_gamma0,
    gamma1 = true_gamma1,
    x_lower = covariate$x_lower,
    x_upper = covariate$x_upper,
    alpha = negative_true$alpha,
    beta = negative_true$beta,
    T = T
  )
  
  events_pos <- simulate_covariate_hawkes_1d(
    x_fun = covariate$x_fun,
    gamma0 = true_gamma0,
    gamma1 = true_gamma1,
    x_lower = covariate$x_lower,
    x_upper = covariate$x_upper,
    alpha = positive_true$alpha,
    beta = positive_true$beta,
    T = T
  )
  
  test_neg <- compute_lrt_covariate_baseline_1d(
    events = events_neg,
    T = T,
    x_fun = covariate$x_fun
  )
  
  test_pos <- compute_lrt_covariate_baseline_1d(
    events = events_pos,
    T = T,
    x_fun = covariate$x_fun
  )
  
  res_neg <- summarize_covariate_baseline_fit_1d(
    label = "covariate_negative_control",
    true_gamma0 = true_gamma0,
    true_gamma1 = true_gamma1,
    true_alpha = negative_true$alpha,
    true_beta = negative_true$beta,
    events = events_neg,
    test = test_neg
  )
  
  res_pos <- summarize_covariate_baseline_fit_1d(
    label = "covariate_positive_control",
    true_gamma0 = true_gamma0,
    true_gamma1 = true_gamma1,
    true_alpha = positive_true$alpha,
    true_beta = positive_true$beta,
    events = events_pos,
    test = test_pos
  )
  
  results <- rbind(res_neg, res_pos)
  
  write.csv(
    results,
    file = file.path(output_dir, "covariate_baseline_demo.csv"),
    row.names = FALSE
  )
  
  plot_covariate_and_baseline_1d(
    covariate = covariate,
    baseline = baseline,
    T = T,
    output_dir = output_dir
  )
  
  cat("\nCovariate baseline demo:\n")
  print(results)
  
  cat("\nFiles written to:\n")
  cat(file.path(output_dir, "covariate_baseline_demo.csv"), "\n")
  cat(file.path(output_dir, "covariate_and_baseline.png"), "\n")
  
  invisible(list(
    covariate = covariate,
    baseline = baseline,
    negative_events = events_neg,
    positive_events = events_pos,
    negative_test = test_neg,
    positive_test = test_pos,
    table = results
  ))
}
run_one_covariate_baseline_mc_case_1d <- function(rep_id,
                                                  label,
                                                  covariate,
                                                  true_gamma0,
                                                  true_gamma1,
                                                  true_alpha,
                                                  true_beta,
                                                  T) {
  out <- tryCatch({
    events <- simulate_covariate_hawkes_1d(
      x_fun = covariate$x_fun,
      gamma0 = true_gamma0,
      gamma1 = true_gamma1,
      x_lower = covariate$x_lower,
      x_upper = covariate$x_upper,
      alpha = true_alpha,
      beta = true_beta,
      T = T
    )
    
    test <- compute_lrt_covariate_baseline_1d(
      events = events,
      T = T,
      x_fun = covariate$x_fun
    )
    
    if (is.null(test$poisson_fit)) {
      poisson_gamma0 <- NA
      poisson_gamma1 <- NA
      poisson_convergence <- NA
    } else {
      poisson_gamma0 <- unname(test$poisson_fit$par["gamma0"])
      poisson_gamma1 <- unname(test$poisson_fit$par["gamma1"])
      poisson_convergence <- test$poisson_fit$convergence
    }
    
    if (is.null(test$hawkes_fit)) {
      hawkes_gamma0 <- NA
      hawkes_gamma1 <- NA
      est_alpha <- NA
      est_beta <- NA
      est_branching <- NA
      hawkes_convergence <- NA
    } else {
      hp <- test$hawkes_fit$par
      
      hawkes_gamma0 <- unname(hp["gamma0"])
      hawkes_gamma1 <- unname(hp["gamma1"])
      est_alpha <- unname(hp["alpha"])
      est_beta <- unname(hp["beta"])
      est_branching <- unname(test$hawkes_fit$branching_ratio)
      hawkes_convergence <- test$hawkes_fit$convergence
    }
    
    data.frame(
      rep = rep_id,
      case = label,
      n_events = length(events),
      
      true_gamma0 = true_gamma0,
      poisson_gamma0 = poisson_gamma0,
      hawkes_gamma0 = hawkes_gamma0,
      
      true_gamma1 = true_gamma1,
      poisson_gamma1 = poisson_gamma1,
      hawkes_gamma1 = hawkes_gamma1,
      
      true_alpha = true_alpha,
      est_alpha = est_alpha,
      
      true_beta = true_beta,
      est_beta = est_beta,
      
      true_branching = true_alpha / true_beta,
      est_branching = est_branching,
      
      poisson_loglik = test$poisson_loglik,
      hawkes_loglik = test$hawkes_loglik,
      lrt_stat = test$lrt_stat,
      
      poisson_convergence = poisson_convergence,
      hawkes_convergence = hawkes_convergence,
      
      error = test$error
    )
  }, error = function(e) {
    data.frame(
      rep = rep_id,
      case = label,
      n_events = NA,
      
      true_gamma0 = true_gamma0,
      poisson_gamma0 = NA,
      hawkes_gamma0 = NA,
      
      true_gamma1 = true_gamma1,
      poisson_gamma1 = NA,
      hawkes_gamma1 = NA,
      
      true_alpha = true_alpha,
      est_alpha = NA,
      
      true_beta = true_beta,
      est_beta = NA,
      
      true_branching = true_alpha / true_beta,
      est_branching = NA,
      
      poisson_loglik = NA,
      hawkes_loglik = NA,
      lrt_stat = NA,
      
      poisson_convergence = NA,
      hawkes_convergence = NA,
      
      error = TRUE
    )
  })
  
  out
}


summarize_covariate_mc_results_1d <- function(results) {
  cases <- split(results, results$case)
  
  summaries <- lapply(names(cases), function(case_name) {
    d <- cases[[case_name]]
    
    d_ok <- d[
      !d$error &
        d$poisson_convergence == 0 &
        d$hawkes_convergence == 0,
    ]
    
    data.frame(
      case = case_name,
      n_rep = nrow(d),
      n_ok = nrow(d_ok),
      
      mean_n_events = mean(d_ok$n_events),
      
      mean_poisson_gamma0 = mean(d_ok$poisson_gamma0),
      mean_hawkes_gamma0 = mean(d_ok$hawkes_gamma0),
      
      mean_poisson_gamma1 = mean(d_ok$poisson_gamma1),
      mean_hawkes_gamma1 = mean(d_ok$hawkes_gamma1),
      
      mean_est_alpha = mean(d_ok$est_alpha),
      mean_est_beta = mean(d_ok$est_beta),
      
      mean_est_branching = mean(d_ok$est_branching),
      sd_est_branching = sd(d_ok$est_branching),
      
      q05_branching = unname(quantile(
        d_ok$est_branching,
        probs = 0.05
      )),
      
      q50_branching = unname(quantile(
        d_ok$est_branching,
        probs = 0.50
      )),
      
      q95_branching = unname(quantile(
        d_ok$est_branching,
        probs = 0.95
      )),
      
      mean_lrt_stat = mean(d_ok$lrt_stat),
      
      q95_lrt_stat = unname(quantile(
        d_ok$lrt_stat,
        probs = 0.95
      )),
      
      q99_lrt_stat = unname(quantile(
        d_ok$lrt_stat,
        probs = 0.99
      ))
    )
  })
  
  do.call(rbind, summaries)
}


plot_covariate_mc_branching_1d <- function(results, output_dir) {
  ok <- results[
    !results$error &
      results$poisson_convergence == 0 &
      results$hawkes_convergence == 0,
  ]
  
  if (nrow(ok) == 0) {
    return(invisible(NULL))
  }
  
  png(
    filename = file.path(
      output_dir,
      "covariate_mc_branching_boxplot.png"
    ),
    width = 900,
    height = 700
  )
  
  boxplot(
    est_branching ~ case,
    data = ok,
    ylab = "Estimated branching ratio",
    xlab = "",
    main = "Covariate baseline: branching estimates"
  )
  
  abline(h = 0, lty = 2)
  abline(h = 0.4, lty = 2)
  
  dev.off()
  
  invisible(NULL)
}


plot_covariate_mc_lrt_1d <- function(results, output_dir) {
  ok <- results[
    !results$error &
      results$poisson_convergence == 0 &
      results$hawkes_convergence == 0,
  ]
  
  if (nrow(ok) == 0) {
    return(invisible(NULL))
  }
  
  png(
    filename = file.path(
      output_dir,
      "covariate_mc_lrt_boxplot.png"
    ),
    width = 900,
    height = 700
  )
  
  boxplot(
    lrt_stat ~ case,
    data = ok,
    ylab = "Likelihood-ratio statistic",
    xlab = "",
    main = "Covariate baseline: LRT"
  )
  
  abline(h = 0, lty = 2)
  
  dev.off()
  
  invisible(NULL)
}


plot_covariate_mc_gamma1_1d <- function(results, output_dir) {
  ok <- results[
    !results$error &
      results$poisson_convergence == 0 &
      results$hawkes_convergence == 0,
  ]
  
  if (nrow(ok) == 0) {
    return(invisible(NULL))
  }
  
  gamma_df <- rbind(
    data.frame(
      case = ok$case,
      model = "Poisson baseline",
      gamma1 = ok$poisson_gamma1
    ),
    data.frame(
      case = ok$case,
      model = "Hawkes baseline",
      gamma1 = ok$hawkes_gamma1
    )
  )
  
  png(
    filename = file.path(
      output_dir,
      "covariate_mc_gamma1_boxplot.png"
    ),
    width = 1000,
    height = 700
  )
  
  boxplot(
    gamma1 ~ case + model,
    data = gamma_df,
    ylab = expression(hat(gamma)[1]),
    xlab = "",
    main = expression("Covariate effect estimate " * gamma[1])
  )
  
  abline(h = 0.8, lty = 2)
  
  dev.off()
  
  invisible(NULL)
}


run_covariate_baseline_monte_carlo_1d <- function(B = 50,
                                                  T = 500,
                                                  seed = 20260620) {
  set.seed(seed)
  
  output_dir <- "results/covariate_baseline_mc"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  covariate <- make_sinusoidal_covariate_1d(
    period = 100,
    phase = 0
  )
  
  true_gamma0 <- log(0.50)
  true_gamma1 <- 0.80
  
  negative_true <- list(
    alpha = 0.00,
    beta = 1.50
  )
  
  positive_true <- list(
    alpha = 0.60,
    beta = 1.50
  )
  
  all_results <- list()
  counter <- 1
  
  for (r in seq_len(B)) {
    cat("rep", r, "of", B, "- covariate negative control\n")
    
    all_results[[counter]] <-
      run_one_covariate_baseline_mc_case_1d(
        rep_id = r,
        label = "covariate_negative_control",
        covariate = covariate,
        true_gamma0 = true_gamma0,
        true_gamma1 = true_gamma1,
        true_alpha = negative_true$alpha,
        true_beta = negative_true$beta,
        T = T
      )
    
    counter <- counter + 1
    
    cat("rep", r, "of", B, "- covariate positive control\n")
    
    all_results[[counter]] <-
      run_one_covariate_baseline_mc_case_1d(
        rep_id = r,
        label = "covariate_positive_control",
        covariate = covariate,
        true_gamma0 = true_gamma0,
        true_gamma1 = true_gamma1,
        true_alpha = positive_true$alpha,
        true_beta = positive_true$beta,
        T = T
      )
    
    counter <- counter + 1
  }
  
  results <- do.call(rbind, all_results)
  summary <- summarize_covariate_mc_results_1d(results)
  
  write.csv(
    results,
    file = file.path(output_dir, "covariate_mc_raw.csv"),
    row.names = FALSE
  )
  
  write.csv(
    summary,
    file = file.path(output_dir, "covariate_mc_summary.csv"),
    row.names = FALSE
  )
  
  plot_covariate_mc_branching_1d(
    results = results,
    output_dir = output_dir
  )
  
  plot_covariate_mc_lrt_1d(
    results = results,
    output_dir = output_dir
  )
  
  plot_covariate_mc_gamma1_1d(
    results = results,
    output_dir = output_dir
  )
  
  cat("\nCovariate baseline Monte Carlo summary:\n")
  print(summary)
  
  cat("\nFiles written to:\n")
  cat(file.path(output_dir, "covariate_mc_raw.csv"), "\n")
  cat(file.path(output_dir, "covariate_mc_summary.csv"), "\n")
  cat(file.path(output_dir, "covariate_mc_branching_boxplot.png"), "\n")
  cat(file.path(output_dir, "covariate_mc_lrt_boxplot.png"), "\n")
  cat(file.path(output_dir, "covariate_mc_gamma1_boxplot.png"), "\n")
  
  invisible(list(
    covariate = covariate,
    raw = results,
    summary = summary
  ))
}
bootstrap_lrt_pvalue_covariate_baseline_1d <- function(events,
                                                       T,
                                                       covariate,
                                                       B = 100,
                                                       seed = NULL,
                                                       label = "dataset") {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  observed <- compute_lrt_covariate_baseline_1d(
    events = events,
    T = T,
    x_fun = covariate$x_fun
  )
  
  if (observed$error || !is.finite(observed$lrt_stat)) {
    stop("Could not compute observed covariate-baseline LRT.")
  }
  
  null_gamma0 <- unname(observed$poisson_fit$par["gamma0"])
  null_gamma1 <- unname(observed$poisson_fit$par["gamma1"])
  
  boot <- data.frame(
    b = seq_len(B),
    n_events = NA_integer_,
    
    null_gamma0 = NA_real_,
    null_gamma1 = NA_real_,
    
    hawkes_gamma0 = NA_real_,
    hawkes_gamma1 = NA_real_,
    
    est_alpha = NA_real_,
    est_beta = NA_real_,
    est_branching = NA_real_,
    
    lrt_stat = NA_real_,
    
    poisson_convergence = NA_integer_,
    hawkes_convergence = NA_integer_,
    
    error = NA
  )
  
  for (b in seq_len(B)) {
    if (b %% 10 == 0) {
      cat(label, "- bootstrap", b, "of", B, "\n")
    }
    
    boot_events <- simulate_covariate_hawkes_1d(
      x_fun = covariate$x_fun,
      gamma0 = null_gamma0,
      gamma1 = null_gamma1,
      x_lower = covariate$x_lower,
      x_upper = covariate$x_upper,
      alpha = 0,
      beta = 1,
      T = T
    )
    
    boot_test <- compute_lrt_covariate_baseline_1d(
      events = boot_events,
      T = T,
      x_fun = covariate$x_fun
    )
    
    boot$n_events[b] <- boot_test$n_events
    boot$lrt_stat[b] <- boot_test$lrt_stat
    boot$error[b] <- boot_test$error
    
    if (!is.null(boot_test$poisson_fit)) {
      boot$null_gamma0[b] <-
        unname(boot_test$poisson_fit$par["gamma0"])
      boot$null_gamma1[b] <-
        unname(boot_test$poisson_fit$par["gamma1"])
      boot$poisson_convergence[b] <-
        boot_test$poisson_fit$convergence
    }
    
    if (!is.null(boot_test$hawkes_fit)) {
      hp <- boot_test$hawkes_fit$par
      
      boot$hawkes_gamma0[b] <- unname(hp["gamma0"])
      boot$hawkes_gamma1[b] <- unname(hp["gamma1"])
      boot$est_alpha[b] <- unname(hp["alpha"])
      boot$est_beta[b] <- unname(hp["beta"])
      boot$est_branching[b] <-
        unname(boot_test$hawkes_fit$branching_ratio)
      boot$hawkes_convergence[b] <-
        boot_test$hawkes_fit$convergence
    }
  }
  
  boot_ok <- boot[
    !boot$error &
      is.finite(boot$lrt_stat) &
      !is.na(boot$lrt_stat) &
      boot$poisson_convergence == 0 &
      boot$hawkes_convergence == 0,
  ]
  
  p_value <- (
    1 + sum(boot_ok$lrt_stat >= observed$lrt_stat)
  ) / (
    1 + nrow(boot_ok)
  )
  
  list(
    label = label,
    observed = observed,
    boot = boot,
    boot_ok = boot_ok,
    p_value = p_value
  )
}


summarize_covariate_bootstrap_test_1d <- function(test) {
  poisson_fit <- test$observed$poisson_fit
  hawkes_fit <- test$observed$hawkes_fit
  
  if (is.null(poisson_fit)) {
    poisson_gamma0 <- NA
    poisson_gamma1 <- NA
  } else {
    poisson_gamma0 <- unname(poisson_fit$par["gamma0"])
    poisson_gamma1 <- unname(poisson_fit$par["gamma1"])
  }
  
  if (is.null(hawkes_fit)) {
    hawkes_gamma0 <- NA
    hawkes_gamma1 <- NA
    est_alpha <- NA
    est_beta <- NA
    est_branching <- NA
  } else {
    hp <- hawkes_fit$par
    
    hawkes_gamma0 <- unname(hp["gamma0"])
    hawkes_gamma1 <- unname(hp["gamma1"])
    est_alpha <- unname(hp["alpha"])
    est_beta <- unname(hp["beta"])
    est_branching <- unname(hawkes_fit$branching_ratio)
  }
  
  data.frame(
    label = test$label,
    n_events = test$observed$n_events,
    
    poisson_gamma0 = poisson_gamma0,
    hawkes_gamma0 = hawkes_gamma0,
    
    poisson_gamma1 = poisson_gamma1,
    hawkes_gamma1 = hawkes_gamma1,
    
    est_alpha = est_alpha,
    est_beta = est_beta,
    est_branching = est_branching,
    
    poisson_loglik = test$observed$poisson_loglik,
    hawkes_loglik = test$observed$hawkes_loglik,
    observed_lrt = test$observed$lrt_stat,
    
    boot_mean_lrt = mean(test$boot_ok$lrt_stat),
    
    boot_q95_lrt = unname(quantile(
      test$boot_ok$lrt_stat,
      probs = 0.95
    )),
    
    boot_q99_lrt = unname(quantile(
      test$boot_ok$lrt_stat,
      probs = 0.99
    )),
    
    p_value = test$p_value,
    n_boot_ok = nrow(test$boot_ok)
  )
}


plot_covariate_bootstrap_lrt_1d <- function(test, output_dir) {
  filename <- paste0(
    "covariate_bootstrap_lrt_",
    gsub("[^A-Za-z0-9_]+", "_", test$label),
    ".png"
  )
  
  png(
    filename = file.path(output_dir, filename),
    width = 900,
    height = 700
  )
  
  hist(
    test$boot_ok$lrt_stat,
    breaks = 30,
    main = paste("Covariate baseline bootstrap LRT:", test$label),
    xlab = "LRT under fitted covariate-Poisson null"
  )
  
  abline(
    v = test$observed$lrt_stat,
    lwd = 3
  )
  
  legend(
    "topright",
    legend = paste(
      "observed LRT =",
      round(test$observed$lrt_stat, 3)
    ),
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
}


run_covariate_baseline_bootstrap_pvalue_demo_1d <-
  function(B = 100,
           T = 500,
           seed = 20260620) {
    set.seed(seed)
    
    output_dir <- "results/covariate_baseline_bootstrap"
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    covariate <- make_sinusoidal_covariate_1d(
      period = 100,
      phase = 0
    )
    
    true_gamma0 <- log(0.50)
    true_gamma1 <- 0.80
    
    negative_true <- list(
      alpha = 0.00,
      beta = 1.50
    )
    
    positive_true <- list(
      alpha = 0.60,
      beta = 1.50
    )
    
    events_neg <- simulate_covariate_hawkes_1d(
      x_fun = covariate$x_fun,
      gamma0 = true_gamma0,
      gamma1 = true_gamma1,
      x_lower = covariate$x_lower,
      x_upper = covariate$x_upper,
      alpha = negative_true$alpha,
      beta = negative_true$beta,
      T = T
    )
    
    events_pos <- simulate_covariate_hawkes_1d(
      x_fun = covariate$x_fun,
      gamma0 = true_gamma0,
      gamma1 = true_gamma1,
      x_lower = covariate$x_lower,
      x_upper = covariate$x_upper,
      alpha = positive_true$alpha,
      beta = positive_true$beta,
      T = T
    )
    
    cat("\nRunning covariate bootstrap test: negative control\n")
    
    test_neg <- bootstrap_lrt_pvalue_covariate_baseline_1d(
      events = events_neg,
      T = T,
      covariate = covariate,
      B = B,
      seed = seed + 1,
      label = "covariate_negative_control"
    )
    
    cat("\nRunning covariate bootstrap test: positive control\n")
    
    test_pos <- bootstrap_lrt_pvalue_covariate_baseline_1d(
      events = events_pos,
      T = T,
      covariate = covariate,
      B = B,
      seed = seed + 2,
      label = "covariate_positive_control"
    )
    
    summary <- rbind(
      summarize_covariate_bootstrap_test_1d(test_neg),
      summarize_covariate_bootstrap_test_1d(test_pos)
    )
    
    write.csv(
      summary,
      file = file.path(
        output_dir,
        "covariate_bootstrap_pvalue_summary.csv"
      ),
      row.names = FALSE
    )
    
    write.csv(
      test_neg$boot,
      file = file.path(
        output_dir,
        "covariate_bootstrap_raw_negative.csv"
      ),
      row.names = FALSE
    )
    
    write.csv(
      test_pos$boot,
      file = file.path(
        output_dir,
        "covariate_bootstrap_raw_positive.csv"
      ),
      row.names = FALSE
    )
    
    plot_covariate_bootstrap_lrt_1d(test_neg, output_dir)
    plot_covariate_bootstrap_lrt_1d(test_pos, output_dir)
    
    cat("\nCovariate baseline bootstrap p-value summary:\n")
    print(summary)
    
    cat("\nFiles written to:\n")
    cat(
      file.path(
        output_dir,
        "covariate_bootstrap_pvalue_summary.csv"
      ),
      "\n"
    )
    cat(
      file.path(
        output_dir,
        "covariate_bootstrap_raw_negative.csv"
      ),
      "\n"
    )
    cat(
      file.path(
        output_dir,
        "covariate_bootstrap_raw_positive.csv"
      ),
      "\n"
    )
    
    invisible(list(
      covariate = covariate,
      negative = test_neg,
      positive = test_pos,
      summary = summary
    ))
  }