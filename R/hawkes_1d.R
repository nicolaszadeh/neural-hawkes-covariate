# R/hawkes_1d.R
#
# Reusable functions for one-dimensional exponential Hawkes processes.
#
# Model:
#
#   lambda(t) = mu + sum_{t_k < t} alpha exp(-beta (t - t_k))
#
# Stability:
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


simulate_poisson_1d <- function(mu, T) {
  if (mu < 0) {
    stop("mu must be non-negative.")
  }
  if (T <= 0) {
    stop("T must be positive.")
  }
  
  n <- rpois(1, mu * T)
  
  if (n == 0) {
    return(numeric(0))
  }
  
  sort(runif(n, min = 0, max = T))
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


poisson_loglik_mle_1d <- function(events, T) {
  n <- length(events)
  
  if (n == 0) {
    return(0)
  }
  
  mu_hat <- n / T
  n * log(mu_hat) - mu_hat * T
}


theta_to_params_1d <- function(theta) {
  mu <- exp(theta[1])
  alpha <- exp(theta[2])
  beta <- alpha + exp(theta[3])
  
  c(mu = mu, alpha = alpha, beta = beta)
}


params_to_theta_1d <- function(mu,
                               alpha,
                               beta) {
  mu <- max(mu, 1e-8)
  alpha <- max(alpha, 1e-8)
  
  if (beta <= alpha) {
    beta <- alpha + 1
  }
  
  c(
    log(mu),
    log(alpha),
    log(beta - alpha)
  )
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
                              init_beta = 1.00,
                              use_multistart = TRUE) {
  if (length(events) == 0) {
    stop("Cannot fit Hawkes model with zero events.")
  }
  
  if (is.null(init_mu)) {
    init_mu <- max(length(events) / T, 1e-3)
  }
  
  starts <- data.frame(
    mu = init_mu,
    alpha = init_alpha,
    beta = init_beta
  )
  
  if (use_multistart) {
    extra <- data.frame(
      mu = c(init_mu, 0.8 * init_mu, 0.5 * init_mu),
      alpha = c(0.01, 0.10, 0.30),
      beta = c(1.00, 1.50, 2.00)
    )
    
    starts <- rbind(starts, extra)
  }
  
  best <- NULL
  
  for (s in seq_len(nrow(starts))) {
    theta0 <- params_to_theta_1d(
      mu = starts$mu[s],
      alpha = starts$alpha[s],
      beta = starts$beta[s]
    )
    
    opt <- tryCatch({
      optim(
        par = theta0,
        fn = hawkes_negloglik_theta_1d,
        events = events,
        T = T,
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
  
  pars <- theta_to_params_1d(best$par)
  
  list(
    par = pars,
    branching_ratio = pars["alpha"] / pars["beta"],
    loglik = -best$value,
    convergence = best$convergence,
    optim = best
  )
}


compute_lrt_1d <- function(events, T) {
  n <- length(events)
  
  if (n == 0) {
    return(list(
      n_events = 0,
      mu_null = 0,
      hawkes_fit = NULL,
      poisson_loglik = 0,
      hawkes_loglik = 0,
      lrt_stat = 0,
      convergence = NA,
      error = FALSE
    ))
  }
  
  out <- tryCatch({
    fit <- fit_hawkes_mle_1d(
      events = events,
      T = T
    )
    
    ll_poisson <- poisson_loglik_mle_1d(
      events = events,
      T = T
    )
    
    ll_hawkes <- fit$loglik
    lrt <- 2 * (ll_hawkes - ll_poisson)
    
    list(
      n_events = n,
      mu_null = n / T,
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
      mu_null = n / T,
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


run_ogata_demo_1d <- function(T = 500,
                              seed = 123) {
  set.seed(seed)
  
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


run_one_ogata_mc_case_1d <- function(rep_id,
                                     label,
                                     true_mu,
                                     true_alpha,
                                     true_beta,
                                     T) {
  out <- tryCatch({
    events <- simulate_hawkes_ogata_1d(
      mu = true_mu,
      alpha = true_alpha,
      beta = true_beta,
      T = T
    )
    
    fit <- fit_hawkes_mle_1d(
      events = events,
      T = T
    )
    
    ll_null <- poisson_loglik_mle_1d(
      events = events,
      T = T
    )
    
    ll_hawkes <- fit$loglik
    lrt_stat <- 2 * (ll_hawkes - ll_null)
    
    est <- fit$par
    
    data.frame(
      rep = rep_id,
      case = label,
      n_events = length(events),
      
      true_mu = true_mu,
      est_mu = unname(est["mu"]),
      
      true_alpha = true_alpha,
      est_alpha = unname(est["alpha"]),
      
      true_beta = true_beta,
      est_beta = unname(est["beta"]),
      
      true_branching = true_alpha / true_beta,
      est_branching = unname(fit$branching_ratio),
      
      poisson_loglik = ll_null,
      hawkes_loglik = ll_hawkes,
      lrt_stat = lrt_stat,
      
      convergence = fit$convergence,
      error = FALSE
    )
  }, error = function(e) {
    data.frame(
      rep = rep_id,
      case = label,
      n_events = NA,
      
      true_mu = true_mu,
      est_mu = NA,
      
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


summarize_mc_results_1d <- function(results) {
  cases <- split(results, results$case)
  
  summaries <- lapply(names(cases), function(case_name) {
    d <- cases[[case_name]]
    d_ok <- d[!d$error & d$convergence == 0, ]
    
    data.frame(
      case = case_name,
      n_rep = nrow(d),
      n_ok = nrow(d_ok),
      
      mean_n_events = mean(d_ok$n_events),
      
      mean_est_mu = mean(d_ok$est_mu),
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
      ))
    )
  })
  
  do.call(rbind, summaries)
}


plot_mc_branching_1d <- function(results, output_dir) {
  ok <- results[!results$error & results$convergence == 0, ]
  
  if (nrow(ok) == 0) {
    return(invisible(NULL))
  }
  
  png(
    filename = file.path(output_dir, "ogata_mc_branching_boxplot.png"),
    width = 900,
    height = 700
  )
  
  boxplot(
    est_branching ~ case,
    data = ok,
    ylab = "Estimated branching ratio",
    xlab = "",
    main = "Ogata Monte Carlo: estimated branching ratio"
  )
  
  abline(h = 0, lty = 2)
  abline(h = 0.4, lty = 2)
  
  dev.off()
  
  invisible(NULL)
}


run_ogata_monte_carlo_1d <- function(B = 50,
                                     T = 500,
                                     seed = 20260620) {
  set.seed(seed)
  
  output_dir <- "results/ogata_monte_carlo"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
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
  
  all_results <- list()
  counter <- 1
  
  for (r in seq_len(B)) {
    cat("rep", r, "of", B, "- negative control\n")
    
    all_results[[counter]] <- run_one_ogata_mc_case_1d(
      rep_id = r,
      label = "negative_control_poisson",
      true_mu = negative_true$mu,
      true_alpha = negative_true$alpha,
      true_beta = negative_true$beta,
      T = T
    )
    
    counter <- counter + 1
    
    cat("rep", r, "of", B, "- positive control\n")
    
    all_results[[counter]] <- run_one_ogata_mc_case_1d(
      rep_id = r,
      label = "positive_control_hawkes",
      true_mu = positive_true$mu,
      true_alpha = positive_true$alpha,
      true_beta = positive_true$beta,
      T = T
    )
    
    counter <- counter + 1
  }
  
  results <- do.call(rbind, all_results)
  summary <- summarize_mc_results_1d(results)
  
  write.csv(
    results,
    file = file.path(output_dir, "ogata_mc_raw.csv"),
    row.names = FALSE
  )
  
  write.csv(
    summary,
    file = file.path(output_dir, "ogata_mc_summary.csv"),
    row.names = FALSE
  )
  
  plot_mc_branching_1d(
    results = results,
    output_dir = output_dir
  )
  
  cat("\nMonte Carlo summary:\n")
  print(summary)
  
  cat("\nFiles written to:\n")
  cat(file.path(output_dir, "ogata_mc_raw.csv"), "\n")
  cat(file.path(output_dir, "ogata_mc_summary.csv"), "\n")
  cat(file.path(output_dir, "ogata_mc_branching_boxplot.png"), "\n")
  
  invisible(list(
    raw = results,
    summary = summary
  ))
}


bootstrap_lrt_pvalue_1d <- function(events,
                                    T,
                                    B = 200,
                                    seed = NULL,
                                    label = "dataset") {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  observed <- compute_lrt_1d(
    events = events,
    T = T
  )
  
  if (observed$error || !is.finite(observed$lrt_stat)) {
    stop("Could not compute observed LRT.")
  }
  
  mu_null <- observed$mu_null
  
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
    
    boot_events <- simulate_poisson_1d(
      mu = mu_null,
      T = T
    )
    
    boot_lrt <- compute_lrt_1d(
      events = boot_events,
      T = T
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


summarize_bootstrap_test_1d <- function(test) {
  fit <- test$observed$hawkes_fit
  
  if (is.null(fit)) {
    est_mu <- NA
    est_alpha <- NA
    est_beta <- NA
    est_branching <- NA
  } else {
    est_mu <- unname(fit$par["mu"])
    est_alpha <- unname(fit$par["alpha"])
    est_beta <- unname(fit$par["beta"])
    est_branching <- unname(fit$branching_ratio)
  }
  
  data.frame(
    label = test$label,
    n_events = test$observed$n_events,
    
    mu_null = test$observed$mu_null,
    
    est_mu = est_mu,
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


plot_bootstrap_lrt_1d <- function(test, output_dir) {
  filename <- paste0(
    "bootstrap_lrt_",
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
    main = paste("Bootstrap null LRT:", test$label),
    xlab = "LRT under fitted Poisson null"
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


run_bootstrap_pvalue_demo_1d <- function(B = 200,
                                         T = 500,
                                         seed = 20260620) {
  set.seed(seed)
  
  output_dir <- "results/ogata_bootstrap_pvalue"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
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
  
  cat("\nRunning bootstrap test: negative control\n")
  
  test_neg <- bootstrap_lrt_pvalue_1d(
    events = events_neg,
    T = T,
    B = B,
    seed = seed + 1,
    label = "negative_control_poisson"
  )
  
  cat("\nRunning bootstrap test: positive control\n")
  
  test_pos <- bootstrap_lrt_pvalue_1d(
    events = events_pos,
    T = T,
    B = B,
    seed = seed + 2,
    label = "positive_control_hawkes"
  )
  
  summary <- rbind(
    summarize_bootstrap_test_1d(test_neg),
    summarize_bootstrap_test_1d(test_pos)
  )
  
  write.csv(
    summary,
    file = file.path(output_dir, "bootstrap_pvalue_summary.csv"),
    row.names = FALSE
  )
  
  write.csv(
    test_neg$boot,
    file = file.path(output_dir, "bootstrap_raw_negative.csv"),
    row.names = FALSE
  )
  
  write.csv(
    test_pos$boot,
    file = file.path(output_dir, "bootstrap_raw_positive.csv"),
    row.names = FALSE
  )
  
  plot_bootstrap_lrt_1d(test_neg, output_dir)
  plot_bootstrap_lrt_1d(test_pos, output_dir)
  
  cat("\nBootstrap p-value summary:\n")
  print(summary)
  
  cat("\nFiles written to:\n")
  cat(file.path(output_dir, "bootstrap_pvalue_summary.csv"), "\n")
  cat(file.path(output_dir, "bootstrap_raw_negative.csv"), "\n")
  cat(file.path(output_dir, "bootstrap_raw_positive.csv"), "\n")
  
  invisible(list(
    negative = test_neg,
    positive = test_pos,
    summary = summary
  ))
}