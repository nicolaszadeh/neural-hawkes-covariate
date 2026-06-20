# R/hawkes_diffusion_baseline_1d.R
#
# Diffusion-driven baseline prototype.
#
# We simulate an Ornstein-Uhlenbeck covariate X_t:
#
#   dX_t = kappa (theta - X_t) dt + sigma dW_t
#
# and use the baseline
#
#   mu(t) = exp(gamma0 + gamma1 X_t).
#
# Conditional on the observed path X_t, the inference is the same as
# the covariate-baseline Hawkes inference from script 12.


source("R/hawkes_covariate_baseline_1d.R")

# Override the generic covariate-baseline integral when the
# covariate is stored on a grid, as for simulated diffusion paths.
#
# The generic version uses integrate(), which is fine for smooth
# deterministic covariates but fragile for interpolated stochastic paths.

covariate_baseline_integral_1d <- function(x_fun,
                                           T,
                                           gamma0,
                                           gamma1) {
  grid_times <- attr(x_fun, "grid_times")
  grid_values <- attr(x_fun, "grid_values")
  
  if (!is.null(grid_times) && !is.null(grid_values)) {
    keep <- grid_times >= 0 & grid_times <= T
    
    tt <- grid_times[keep]
    xx <- grid_values[keep]
    
    if (length(tt) == 0 || tt[1] > 0) {
      tt <- c(0, tt)
      xx <- c(as.numeric(x_fun(0)), xx)
    }
    
    if (tail(tt, 1) < T) {
      tt <- c(tt, T)
      xx <- c(xx, as.numeric(x_fun(T)))
    }
    
    yy <- exp(gamma0 + gamma1 * xx)
    
    return(sum(
      diff(tt) * (head(yy, -1) + tail(yy, -1)) / 2
    ))
  }
  
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

simulate_ou_path_1d <- function(T = 500,
                                dt = 0.10,
                                kappa = 0.20,
                                theta = 0.00,
                                sigma = 0.60,
                                x0 = 0.00) {
  if (T <= 0) {
    stop("T must be positive.")
  }
  if (dt <= 0) {
    stop("dt must be positive.")
  }
  if (kappa <= 0) {
    stop("kappa must be positive.")
  }
  if (sigma < 0) {
    stop("sigma must be non-negative.")
  }
  
  times <- seq(0, T, by = dt)
  n <- length(times)
  
  x <- numeric(n)
  x[1] <- x0
  
  for (k in seq_len(n - 1)) {
    dW <- sqrt(dt) * rnorm(1)
    
    x[k + 1] <- x[k] +
      kappa * (theta - x[k]) * dt +
      sigma * dW
  }
  
  x_fun <- approxfun(
    x = times,
    y = x,
    method = "linear",
    rule = 2
  )
  
  attr(x_fun, "grid_times") <- times
  attr(x_fun, "grid_values") <- x
  
  list(
    times = times,
    values = x,
    x_fun = x_fun,
    x_lower = min(x),
    x_upper = max(x),
    T = T,
    dt = dt,
    kappa = kappa,
    theta = theta,
    sigma = sigma,
    x0 = x0
  )
}


make_diffusion_baseline_from_path_1d <- function(path,
                                                 gamma0,
                                                 gamma1) {
  make_covariate_baseline_1d(
    x_fun = path$x_fun,
    gamma0 = gamma0,
    gamma1 = gamma1,
    x_lower = path$x_lower,
    x_upper = path$x_upper
  )
}


simulate_diffusion_baseline_hawkes_1d <- function(path,
                                                  gamma0,
                                                  gamma1,
                                                  alpha,
                                                  beta,
                                                  T) {
  simulate_covariate_hawkes_1d(
    x_fun = path$x_fun,
    gamma0 = gamma0,
    gamma1 = gamma1,
    x_lower = path$x_lower,
    x_upper = path$x_upper,
    alpha = alpha,
    beta = beta,
    T = T
  )
}


plot_diffusion_path_and_baseline_1d <- function(path,
                                                baseline,
                                                T,
                                                output_dir) {
  png(
    filename = file.path(
      output_dir,
      "diffusion_path_and_baseline.png"
    ),
    width = 900,
    height = 800
  )
  
  old_par <- par(no.readonly = TRUE)
  
  par(mfrow = c(2, 1))
  
  plot(
    path$times,
    path$values,
    type = "l",
    xlab = "t",
    ylab = expression(X[t]),
    main = "Simulated diffusion covariate"
  )
  
  curve(
    baseline$mu_fun(x),
    from = 0,
    to = T,
    xlab = "t",
    ylab = expression(mu(t)),
    main = expression(mu(t) == exp(gamma[0] + gamma[1] * X[t]))
  )
  
  par(old_par)
  dev.off()
  
  invisible(NULL)
}


run_diffusion_baseline_demo_1d <- function(T = 500,
                                           dt = 0.10,
                                           seed = 20260620) {
  set.seed(seed)
  
  output_dir <- "results/diffusion_baseline"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  path <- simulate_ou_path_1d(
    T = T,
    dt = dt,
    kappa = 0.20,
    theta = 0.00,
    sigma = 0.60,
    x0 = 0.00
  )
  
  true_gamma0 <- log(0.40)
  true_gamma1 <- 0.50
  
  baseline <- make_diffusion_baseline_from_path_1d(
    path = path,
    gamma0 = true_gamma0,
    gamma1 = true_gamma1
  )
  
  negative_true <- list(
    alpha = 0.00,
    beta = 1.50
  )
  
  positive_true <- list(
    alpha = 0.60,
    beta = 1.50
  )
  
  events_neg <- simulate_diffusion_baseline_hawkes_1d(
    path = path,
    gamma0 = true_gamma0,
    gamma1 = true_gamma1,
    alpha = negative_true$alpha,
    beta = negative_true$beta,
    T = T
  )
  
  events_pos <- simulate_diffusion_baseline_hawkes_1d(
    path = path,
    gamma0 = true_gamma0,
    gamma1 = true_gamma1,
    alpha = positive_true$alpha,
    beta = positive_true$beta,
    T = T
  )
  
  test_neg <- compute_lrt_covariate_baseline_1d(
    events = events_neg,
    T = T,
    x_fun = path$x_fun
  )
  
  test_pos <- compute_lrt_covariate_baseline_1d(
    events = events_pos,
    T = T,
    x_fun = path$x_fun
  )
  
  res_neg <- summarize_covariate_baseline_fit_1d(
    label = "diffusion_negative_control",
    true_gamma0 = true_gamma0,
    true_gamma1 = true_gamma1,
    true_alpha = negative_true$alpha,
    true_beta = negative_true$beta,
    events = events_neg,
    test = test_neg
  )
  
  res_pos <- summarize_covariate_baseline_fit_1d(
    label = "diffusion_positive_control",
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
    file = file.path(output_dir, "diffusion_baseline_demo.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      t = path$times,
      X = path$values
    ),
    file = file.path(output_dir, "diffusion_path.csv"),
    row.names = FALSE
  )
  
  plot_diffusion_path_and_baseline_1d(
    path = path,
    baseline = baseline,
    T = T,
    output_dir = output_dir
  )
  
  cat("\nDiffusion baseline demo:\n")
  print(results)
  
  cat("\nDiffusion path range:\n")
  print(c(
    min_X = path$x_lower,
    max_X = path$x_upper
  ))
  
  cat("\nFiles written to:\n")
  cat(file.path(output_dir, "diffusion_baseline_demo.csv"), "\n")
  cat(file.path(output_dir, "diffusion_path.csv"), "\n")
  cat(file.path(output_dir, "diffusion_path_and_baseline.png"), "\n")
  
  invisible(list(
    path = path,
    baseline = baseline,
    negative_events = events_neg,
    positive_events = events_pos,
    negative_test = test_neg,
    positive_test = test_pos,
    table = results
  ))
}
run_one_diffusion_baseline_mc_case_1d <- function(rep_id,
                                                  label,
                                                  path,
                                                  true_gamma0,
                                                  true_gamma1,
                                                  true_alpha,
                                                  true_beta,
                                                  T) {
  out <- tryCatch({
    events <- simulate_diffusion_baseline_hawkes_1d(
      path = path,
      gamma0 = true_gamma0,
      gamma1 = true_gamma1,
      alpha = true_alpha,
      beta = true_beta,
      T = T
    )
    
    test <- compute_lrt_covariate_baseline_1d(
      events = events,
      T = T,
      x_fun = path$x_fun
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
      
      path_min = path$x_lower,
      path_max = path$x_upper,
      
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
      
      path_min = path$x_lower,
      path_max = path$x_upper,
      
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


summarize_diffusion_mc_results_1d <- function(results) {
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
      
      mean_path_min = mean(d_ok$path_min),
      mean_path_max = mean(d_ok$path_max),
      
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


plot_diffusion_mc_branching_1d <- function(results, output_dir) {
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
      "diffusion_mc_branching_boxplot.png"
    ),
    width = 900,
    height = 700
  )
  
  boxplot(
    est_branching ~ case,
    data = ok,
    ylab = "Estimated branching ratio",
    xlab = "",
    main = "Diffusion baseline: branching estimates"
  )
  
  abline(h = 0, lty = 2)
  abline(h = 0.4, lty = 2)
  
  dev.off()
  
  invisible(NULL)
}


plot_diffusion_mc_lrt_1d <- function(results, output_dir) {
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
      "diffusion_mc_lrt_boxplot.png"
    ),
    width = 900,
    height = 700
  )
  
  boxplot(
    lrt_stat ~ case,
    data = ok,
    ylab = "Likelihood-ratio statistic",
    xlab = "",
    main = "Diffusion baseline: LRT"
  )
  
  abline(h = 0, lty = 2)
  
  dev.off()
  
  invisible(NULL)
}


plot_diffusion_mc_gamma1_1d <- function(results, output_dir) {
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
      "diffusion_mc_gamma1_boxplot.png"
    ),
    width = 1000,
    height = 700
  )
  
  boxplot(
    gamma1 ~ case + model,
    data = gamma_df,
    ylab = expression(hat(gamma)[1]),
    xlab = "",
    main = expression("Diffusion baseline: " * gamma[1])
  )
  
  abline(h = 0.5, lty = 2)
  
  dev.off()
  
  invisible(NULL)
}


run_diffusion_baseline_monte_carlo_1d <- function(B = 50,
                                                  T = 500,
                                                  dt = 0.10,
                                                  seed = 20260620) {
  set.seed(seed)
  
  output_dir <- "results/diffusion_baseline_mc"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  true_gamma0 <- log(0.40)
  true_gamma1 <- 0.50
  
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
    cat("rep", r, "of", B, "- simulating diffusion path\n")
    
    path <- simulate_ou_path_1d(
      T = T,
      dt = dt,
      kappa = 0.20,
      theta = 0.00,
      sigma = 0.60,
      x0 = 0.00
    )
    
    cat("rep", r, "of", B, "- diffusion negative control\n")
    
    all_results[[counter]] <-
      run_one_diffusion_baseline_mc_case_1d(
        rep_id = r,
        label = "diffusion_negative_control",
        path = path,
        true_gamma0 = true_gamma0,
        true_gamma1 = true_gamma1,
        true_alpha = negative_true$alpha,
        true_beta = negative_true$beta,
        T = T
      )
    
    counter <- counter + 1
    
    cat("rep", r, "of", B, "- diffusion positive control\n")
    
    all_results[[counter]] <-
      run_one_diffusion_baseline_mc_case_1d(
        rep_id = r,
        label = "diffusion_positive_control",
        path = path,
        true_gamma0 = true_gamma0,
        true_gamma1 = true_gamma1,
        true_alpha = positive_true$alpha,
        true_beta = positive_true$beta,
        T = T
      )
    
    counter <- counter + 1
  }
  
  results <- do.call(rbind, all_results)
  summary <- summarize_diffusion_mc_results_1d(results)
  
  write.csv(
    results,
    file = file.path(output_dir, "diffusion_mc_raw.csv"),
    row.names = FALSE
  )
  
  write.csv(
    summary,
    file = file.path(output_dir, "diffusion_mc_summary.csv"),
    row.names = FALSE
  )
  
  plot_diffusion_mc_branching_1d(
    results = results,
    output_dir = output_dir
  )
  
  plot_diffusion_mc_lrt_1d(
    results = results,
    output_dir = output_dir
  )
  
  plot_diffusion_mc_gamma1_1d(
    results = results,
    output_dir = output_dir
  )
  
  cat("\nDiffusion baseline Monte Carlo summary:\n")
  print(summary)
  
  cat("\nFiles written to:\n")
  cat(file.path(output_dir, "diffusion_mc_raw.csv"), "\n")
  cat(file.path(output_dir, "diffusion_mc_summary.csv"), "\n")
  cat(file.path(output_dir, "diffusion_mc_branching_boxplot.png"), "\n")
  cat(file.path(output_dir, "diffusion_mc_lrt_boxplot.png"), "\n")
  cat(file.path(output_dir, "diffusion_mc_gamma1_boxplot.png"), "\n")
  
  invisible(list(
    raw = results,
    summary = summary
  ))
}