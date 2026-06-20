# R/hawkes_td_baseline_1d.R
#
# One-dimensional exponential Hawkes process with known
# time-dependent deterministic baseline:
#
#   lambda(t) = mu(t)
#             + sum_{t_k < t} alpha exp(-beta (t - t_k))
#
# This is the first bridge toward diffusion-driven baselines.


source("R/hawkes_1d.R")


make_sinusoidal_baseline_1d <- function(base_mu = 0.50,
                                        amplitude = 0.50,
                                        period = 100,
                                        phase = 0) {
  if (base_mu <= 0) {
    stop("base_mu must be positive.")
  }
  if (amplitude < 0 || amplitude >= 1) {
    stop("amplitude must satisfy 0 <= amplitude < 1.")
  }
  if (period <= 0) {
    stop("period must be positive.")
  }
  
  mu_fun <- function(t) {
    base_mu * (
      1 + amplitude * sin(2 * pi * t / period + phase)
    )
  }
  
  list(
    mu_fun = mu_fun,
    mu_upper = base_mu * (1 + amplitude),
    base_mu = base_mu,
    amplitude = amplitude,
    period = period,
    phase = phase
  )
}


baseline_integral_1d <- function(mu_fun, T) {
  if (T <= 0) {
    stop("T must be positive.")
  }
  
  integrate(
    f = function(t) mu_fun(t),
    lower = 0,
    upper = T,
    subdivisions = 1000,
    rel.tol = 1e-8
  )$value
}


simulate_hawkes_ogata_td_baseline_1d <- function(mu_fun,
                                                 mu_upper,
                                                 alpha,
                                                 beta,
                                                 T,
                                                 max_events = 1e6) {
  if (mu_upper <= 0) {
    stop("mu_upper must be positive.")
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
    lambda_bar <- mu_upper + excitation
    
    if (lambda_bar <= 0) {
      break
    }
    
    w <- rexp(1, rate = lambda_bar)
    t_candidate <- t + w
    
    if (t_candidate > T) {
      break
    }
    
    excitation_candidate <- excitation * exp(-beta * w)
    
    lambda_candidate <- mu_fun(t_candidate) +
      excitation_candidate
    
    accept_prob <- lambda_candidate / lambda_bar
    
    if (accept_prob > 1 + 1e-10) {
      stop("Dominating intensity failed. Increase mu_upper.")
    }
    
    accept_prob <- min(accept_prob, 1)
    
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


hawkes_loglik_td_baseline_1d <- function(events,
                                         T,
                                         mu_fun,
                                         alpha,
                                         beta,
                                         mu_integral = NULL) {
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
      
      lambda_k <- mu_fun(events[k]) + excitation
      
      if (lambda_k <= 0) {
        return(-Inf)
      }
      
      log_part <- log_part + log(lambda_k)
      
      excitation <- excitation + alpha
      last_t <- events[k]
    }
  }
  
  if (is.null(mu_integral)) {
    mu_integral <- baseline_integral_1d(mu_fun, T)
  }
  
  compensator <- mu_integral
  
  if (n > 0 && alpha > 0) {
    compensator <- compensator +
      sum(alpha / beta * (1 - exp(-beta * (T - events))))
  }
  
  log_part - compensator
}


td_poisson_loglik_known_baseline_1d <- function(events,
                                                T,
                                                mu_fun,
                                                mu_integral = NULL) {
  if (is.null(mu_integral)) {
    mu_integral <- baseline_integral_1d(mu_fun, T)
  }
  
  if (length(events) == 0) {
    return(-mu_integral)
  }
  
  if (any(events <= 0) || any(events > T)) {
    return(-Inf)
  }
  
  values <- mu_fun(events)
  
  if (any(values <= 0)) {
    return(-Inf)
  }
  
  sum(log(values)) - mu_integral
}


theta_to_excitation_params_1d <- function(theta) {
  alpha <- exp(theta[1])
  beta <- alpha + exp(theta[2])
  
  c(alpha = alpha, beta = beta)
}


excitation_params_to_theta_1d <- function(alpha,
                                          beta) {
  alpha <- max(alpha, 1e-8)
  
  if (beta <= alpha) {
    beta <- alpha + 1
  }
  
  c(
    log(alpha),
    log(beta - alpha)
  )
}


hawkes_td_negloglik_theta_1d <- function(theta,
                                         events,
                                         T,
                                         mu_fun,
                                         mu_integral) {
  pars <- theta_to_excitation_params_1d(theta)
  
  ll <- hawkes_loglik_td_baseline_1d(
    events = events,
    T = T,
    mu_fun = mu_fun,
    alpha = pars["alpha"],
    beta = pars["beta"],
    mu_integral = mu_integral
  )
  
  if (!is.finite(ll)) {
    return(1e100)
  }
  
  -ll
}


fit_hawkes_td_baseline_mle_1d <- function(events,
                                          T,
                                          mu_fun,
                                          init_alpha = 0.05,
                                          init_beta = 1.00,
                                          use_multistart = TRUE) {
  if (length(events) == 0) {
    stop("Cannot fit Hawkes model with zero events.")
  }
  
  mu_integral <- baseline_integral_1d(mu_fun, T)
  
  starts <- data.frame(
    alpha = init_alpha,
    beta = init_beta
  )
  
  if (use_multistart) {
    extra <- data.frame(
      alpha = c(0.01, 0.10, 0.30, 0.60),
      beta = c(1.00, 1.50, 2.00, 2.50)
    )
    
    starts <- rbind(starts, extra)
  }
  
  best <- NULL
  
  for (s in seq_len(nrow(starts))) {
    theta0 <- excitation_params_to_theta_1d(
      alpha = starts$alpha[s],
      beta = starts$beta[s]
    )
    
    opt <- tryCatch({
      optim(
        par = theta0,
        fn = hawkes_td_negloglik_theta_1d,
        events = events,
        T = T,
        mu_fun = mu_fun,
        mu_integral = mu_integral,
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
    stop("All optimizations failed.")
  }
  
  pars <- theta_to_excitation_params_1d(best$par)
  
  list(
    par = pars,
    branching_ratio = pars["alpha"] / pars["beta"],
    loglik = -best$value,
    convergence = best$convergence,
    optim = best,
    mu_integral = mu_integral
  )
}


compute_lrt_td_baseline_1d <- function(events,
                                       T,
                                       mu_fun) {
  n <- length(events)
  mu_integral <- baseline_integral_1d(mu_fun, T)
  
  if (n == 0) {
    return(list(
      n_events = 0,
      hawkes_fit = NULL,
      poisson_loglik = -mu_integral,
      hawkes_loglik = -mu_integral,
      lrt_stat = 0,
      convergence = NA,
      error = FALSE
    ))
  }
  
  out <- tryCatch({
    fit <- fit_hawkes_td_baseline_mle_1d(
      events = events,
      T = T,
      mu_fun = mu_fun
    )
    
    ll_poisson <- td_poisson_loglik_known_baseline_1d(
      events = events,
      T = T,
      mu_fun = mu_fun,
      mu_integral = mu_integral
    )
    
    ll_hawkes <- fit$loglik
    lrt <- 2 * (ll_hawkes - ll_poisson)
    
    list(
      n_events = n,
      hawkes_fit = fit,
      poisson_loglik = ll_poisson,
      hawkes_loglik = ll_hawkes,
      lrt_stat = lrt,
      convergence = fit$convergence,
      error = FALSE
    )
  }, error = function(e) {
    list(
      n_events = n,
      hawkes_fit = NULL,
      poisson_loglik = NA,
      hawkes_loglik = NA,
      lrt_stat = NA,
      convergence = NA,
      error = TRUE
    )
  })
  
  out
}


summarize_td_fit_1d <- function(label,
                                true_alpha,
                                true_beta,
                                events,
                                test) {
  fit <- test$hawkes_fit
  
  if (is.null(fit)) {
    est_alpha <- NA
    est_beta <- NA
    est_branching <- NA
    loglik <- NA
    convergence <- NA
  } else {
    est_alpha <- unname(fit$par["alpha"])
    est_beta <- unname(fit$par["beta"])
    est_branching <- unname(fit$branching_ratio)
    loglik <- fit$loglik
    convergence <- fit$convergence
  }
  
  data.frame(
    case = label,
    n_events = length(events),
    
    true_alpha = true_alpha,
    est_alpha = est_alpha,
    
    true_beta = true_beta,
    est_beta = est_beta,
    
    true_branching = true_alpha / true_beta,
    est_branching = est_branching,
    
    poisson_loglik = test$poisson_loglik,
    hawkes_loglik = loglik,
    lrt_stat = test$lrt_stat,
    
    convergence = convergence,
    error = test$error
  )
}


plot_td_baseline_1d <- function(baseline,
                                T,
                                output_dir) {
  png(
    filename = file.path(output_dir, "td_baseline_mu_t.png"),
    width = 900,
    height = 700
  )
  
  curve(
    baseline$mu_fun(x),
    from = 0,
    to = T,
    xlab = "t",
    ylab = expression(mu(t)),
    main = "Known time-dependent baseline"
  )
  
  dev.off()
  
  invisible(NULL)
}


run_td_baseline_demo_1d <- function(T = 500,
                                    seed = 20260620) {
  set.seed(seed)
  
  output_dir <- "results/time_dependent_baseline"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  baseline <- make_sinusoidal_baseline_1d(
    base_mu = 0.50,
    amplitude = 0.50,
    period = 100
  )
  
  negative_true <- list(
    alpha = 0.00,
    beta = 1.50
  )
  
  positive_true <- list(
    alpha = 0.60,
    beta = 1.50
  )
  
  events_neg <- simulate_hawkes_ogata_td_baseline_1d(
    mu_fun = baseline$mu_fun,
    mu_upper = baseline$mu_upper,
    alpha = negative_true$alpha,
    beta = negative_true$beta,
    T = T
  )
  
  events_pos <- simulate_hawkes_ogata_td_baseline_1d(
    mu_fun = baseline$mu_fun,
    mu_upper = baseline$mu_upper,
    alpha = positive_true$alpha,
    beta = positive_true$beta,
    T = T
  )
  
  test_neg <- compute_lrt_td_baseline_1d(
    events = events_neg,
    T = T,
    mu_fun = baseline$mu_fun
  )
  
  test_pos <- compute_lrt_td_baseline_1d(
    events = events_pos,
    T = T,
    mu_fun = baseline$mu_fun
  )
  
  res_neg <- summarize_td_fit_1d(
    label = "td_baseline_negative_control",
    true_alpha = negative_true$alpha,
    true_beta = negative_true$beta,
    events = events_neg,
    test = test_neg
  )
  
  res_pos <- summarize_td_fit_1d(
    label = "td_baseline_positive_control",
    true_alpha = positive_true$alpha,
    true_beta = positive_true$beta,
    events = events_pos,
    test = test_pos
  )
  
  results <- rbind(res_neg, res_pos)
  
  write.csv(
    results,
    file = file.path(output_dir, "td_baseline_demo.csv"),
    row.names = FALSE
  )
  
  plot_td_baseline_1d(
    baseline = baseline,
    T = T,
    output_dir = output_dir
  )
  
  cat("\nTime-dependent baseline demo:\n")
  print(results)
  
  cat("\nFiles written to:\n")
  cat(file.path(output_dir, "td_baseline_demo.csv"), "\n")
  cat(file.path(output_dir, "td_baseline_mu_t.png"), "\n")
  
  invisible(list(
    baseline = baseline,
    negative_events = events_neg,
    positive_events = events_pos,
    negative_test = test_neg,
    positive_test = test_pos,
    table = results
  ))
}
run_one_td_baseline_mc_case_1d <- function(rep_id,
                                           label,
                                           baseline,
                                           true_alpha,
                                           true_beta,
                                           T) {
  out <- tryCatch({
    events <- simulate_hawkes_ogata_td_baseline_1d(
      mu_fun = baseline$mu_fun,
      mu_upper = baseline$mu_upper,
      alpha = true_alpha,
      beta = true_beta,
      T = T
    )
    
    test <- compute_lrt_td_baseline_1d(
      events = events,
      T = T,
      mu_fun = baseline$mu_fun
    )
    
    fit <- test$hawkes_fit
    
    if (is.null(fit)) {
      est_alpha <- NA
      est_beta <- NA
      est_branching <- NA
      convergence <- NA
    } else {
      est_alpha <- unname(fit$par["alpha"])
      est_beta <- unname(fit$par["beta"])
      est_branching <- unname(fit$branching_ratio)
      convergence <- fit$convergence
    }
    
    data.frame(
      rep = rep_id,
      case = label,
      n_events = length(events),
      
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
  }, error = function(e) {
    data.frame(
      rep = rep_id,
      case = label,
      n_events = NA,
      
      true_alpha = true_alpha,
      est_alpha = NA,
      
      true_beta = true_beta,
      est_beta = NA,
      
      true_branching = true_alpha / true_beta,
      est_branching = NA,
      
      poisson_loglik = NA,
      hawkes_loglik = NA,
      lrt_stat = NA,
      
      convergence = NA,
      error = TRUE
    )
  })
  
  out
}


summarize_td_mc_results_1d <- function(results) {
  cases <- split(results, results$case)
  
  summaries <- lapply(names(cases), function(case_name) {
    d <- cases[[case_name]]
    d_ok <- d[!d$error & d$convergence == 0, ]
    
    data.frame(
      case = case_name,
      n_rep = nrow(d),
      n_ok = nrow(d_ok),
      
      mean_n_events = mean(d_ok$n_events),
      
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


plot_td_mc_branching_1d <- function(results, output_dir) {
  ok <- results[!results$error & results$convergence == 0, ]
  
  if (nrow(ok) == 0) {
    return(invisible(NULL))
  }
  
  png(
    filename = file.path(
      output_dir,
      "td_baseline_mc_branching_boxplot.png"
    ),
    width = 900,
    height = 700
  )
  
  boxplot(
    est_branching ~ case,
    data = ok,
    ylab = "Estimated branching ratio",
    xlab = "",
    main = "Time-dependent baseline: branching estimates"
  )
  
  abline(h = 0, lty = 2)
  abline(h = 0.4, lty = 2)
  
  dev.off()
  
  invisible(NULL)
}


plot_td_mc_lrt_1d <- function(results, output_dir) {
  ok <- results[!results$error & results$convergence == 0, ]
  
  if (nrow(ok) == 0) {
    return(invisible(NULL))
  }
  
  png(
    filename = file.path(
      output_dir,
      "td_baseline_mc_lrt_boxplot.png"
    ),
    width = 900,
    height = 700
  )
  
  boxplot(
    lrt_stat ~ case,
    data = ok,
    ylab = "Likelihood-ratio statistic",
    xlab = "",
    main = "Time-dependent baseline: LRT"
  )
  
  abline(h = 0, lty = 2)
  
  dev.off()
  
  invisible(NULL)
}


run_td_baseline_monte_carlo_1d <- function(B = 50,
                                           T = 500,
                                           seed = 20260620) {
  set.seed(seed)
  
  output_dir <- "results/time_dependent_baseline_mc"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  baseline <- make_sinusoidal_baseline_1d(
    base_mu = 0.50,
    amplitude = 0.50,
    period = 100
  )
  
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
    cat("rep", r, "of", B, "- td negative control\n")
    
    all_results[[counter]] <- run_one_td_baseline_mc_case_1d(
      rep_id = r,
      label = "td_baseline_negative_control",
      baseline = baseline,
      true_alpha = negative_true$alpha,
      true_beta = negative_true$beta,
      T = T
    )
    
    counter <- counter + 1
    
    cat("rep", r, "of", B, "- td positive control\n")
    
    all_results[[counter]] <- run_one_td_baseline_mc_case_1d(
      rep_id = r,
      label = "td_baseline_positive_control",
      baseline = baseline,
      true_alpha = positive_true$alpha,
      true_beta = positive_true$beta,
      T = T
    )
    
    counter <- counter + 1
  }
  
  results <- do.call(rbind, all_results)
  summary <- summarize_td_mc_results_1d(results)
  
  write.csv(
    results,
    file = file.path(output_dir, "td_baseline_mc_raw.csv"),
    row.names = FALSE
  )
  
  write.csv(
    summary,
    file = file.path(output_dir, "td_baseline_mc_summary.csv"),
    row.names = FALSE
  )
  
  plot_td_mc_branching_1d(
    results = results,
    output_dir = output_dir
  )
  
  plot_td_mc_lrt_1d(
    results = results,
    output_dir = output_dir
  )
  
  cat("\nTime-dependent baseline Monte Carlo summary:\n")
  print(summary)
  
  cat("\nFiles written to:\n")
  cat(file.path(output_dir, "td_baseline_mc_raw.csv"), "\n")
  cat(file.path(output_dir, "td_baseline_mc_summary.csv"), "\n")
  cat(file.path(output_dir, "td_baseline_mc_branching_boxplot.png"), "\n")
  cat(file.path(output_dir, "td_baseline_mc_lrt_boxplot.png"), "\n")
  
  invisible(list(
    baseline = baseline,
    raw = results,
    summary = summary
  ))
}
bootstrap_lrt_pvalue_td_baseline_1d <- function(events,
                                                T,
                                                baseline,
                                                B = 200,
                                                seed = NULL,
                                                label = "dataset") {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  observed <- compute_lrt_td_baseline_1d(
    events = events,
    T = T,
    mu_fun = baseline$mu_fun
  )
  
  if (observed$error || !is.finite(observed$lrt_stat)) {
    stop("Could not compute observed LRT.")
  }
  
  boot <- data.frame(
    b = seq_len(B),
    n_events = NA_integer_,
    lrt_stat = NA_real_,
    convergence = NA_integer_,
    error = NA
  )
  
  for (b in seq_len(B)) {
    if (b %% 20 == 0) {
      cat(label, "- bootstrap", b, "of", B, "\n")
    }
    
    boot_events <- simulate_hawkes_ogata_td_baseline_1d(
      mu_fun = baseline$mu_fun,
      mu_upper = baseline$mu_upper,
      alpha = 0,
      beta = 1,
      T = T
    )
    
    boot_lrt <- compute_lrt_td_baseline_1d(
      events = boot_events,
      T = T,
      mu_fun = baseline$mu_fun
    )
    
    boot$n_events[b] <- boot_lrt$n_events
    boot$lrt_stat[b] <- boot_lrt$lrt_stat
    boot$convergence[b] <- boot_lrt$convergence
    boot$error[b] <- boot_lrt$error
  }
  
  boot_ok <- boot[
    !boot$error &
      is.finite(boot$lrt_stat) &
      !is.na(boot$lrt_stat),
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


summarize_td_bootstrap_test_1d <- function(test) {
  fit <- test$observed$hawkes_fit
  
  if (is.null(fit)) {
    est_alpha <- NA
    est_beta <- NA
    est_branching <- NA
  } else {
    est_alpha <- unname(fit$par["alpha"])
    est_beta <- unname(fit$par["beta"])
    est_branching <- unname(fit$branching_ratio)
  }
  
  data.frame(
    label = test$label,
    n_events = test$observed$n_events,
    
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


plot_td_bootstrap_lrt_1d <- function(test, output_dir) {
  filename <- paste0(
    "td_bootstrap_lrt_",
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
    main = paste("TD baseline bootstrap LRT:", test$label),
    xlab = "LRT under time-dependent Poisson null"
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


run_td_baseline_bootstrap_pvalue_demo_1d <- function(B = 200,
                                                     T = 500,
                                                     seed = 20260620) {
  set.seed(seed)
  
  output_dir <- "results/time_dependent_baseline_bootstrap"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  baseline <- make_sinusoidal_baseline_1d(
    base_mu = 0.50,
    amplitude = 0.50,
    period = 100
  )
  
  negative_true <- list(
    alpha = 0.00,
    beta = 1.50
  )
  
  positive_true <- list(
    alpha = 0.60,
    beta = 1.50
  )
  
  events_neg <- simulate_hawkes_ogata_td_baseline_1d(
    mu_fun = baseline$mu_fun,
    mu_upper = baseline$mu_upper,
    alpha = negative_true$alpha,
    beta = negative_true$beta,
    T = T
  )
  
  events_pos <- simulate_hawkes_ogata_td_baseline_1d(
    mu_fun = baseline$mu_fun,
    mu_upper = baseline$mu_upper,
    alpha = positive_true$alpha,
    beta = positive_true$beta,
    T = T
  )
  
  cat("\nRunning TD bootstrap test: negative control\n")
  
  test_neg <- bootstrap_lrt_pvalue_td_baseline_1d(
    events = events_neg,
    T = T,
    baseline = baseline,
    B = B,
    seed = seed + 1,
    label = "td_baseline_negative_control"
  )
  
  cat("\nRunning TD bootstrap test: positive control\n")
  
  test_pos <- bootstrap_lrt_pvalue_td_baseline_1d(
    events = events_pos,
    T = T,
    baseline = baseline,
    B = B,
    seed = seed + 2,
    label = "td_baseline_positive_control"
  )
  
  summary <- rbind(
    summarize_td_bootstrap_test_1d(test_neg),
    summarize_td_bootstrap_test_1d(test_pos)
  )
  
  write.csv(
    summary,
    file = file.path(output_dir, "td_bootstrap_pvalue_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    test_neg$boot,
    file = file.path(output_dir, "td_bootstrap_raw_negative.csv"),
    row.names = FALSE
  )
  
  write.csv(
    test_pos$boot,
    file = file.path(output_dir, "td_bootstrap_raw_positive.csv"),
    row.names = FALSE
  )
  
  plot_td_bootstrap_lrt_1d(test_neg, output_dir)
  plot_td_bootstrap_lrt_1d(test_pos, output_dir)
  
  cat("\nTime-dependent baseline bootstrap p-value summary:\n")
  print(summary)
  
  cat("\nFiles written to:\n")
  cat(file.path(output_dir, "td_bootstrap_pvalue_summary.csv"), "\n")
  cat(file.path(output_dir, "td_bootstrap_raw_negative.csv"), "\n")
  cat(file.path(output_dir, "td_bootstrap_raw_positive.csv"), "\n")
  
  invisible(list(
    baseline = baseline,
    negative = test_neg,
    positive = test_pos,
    summary = summary
  ))
}