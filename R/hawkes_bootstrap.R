# ============================================================
# Bootstrap utilities for Hawkes likelihood-ratio tests
# ============================================================


# ============================================================
# Check that required functions have been loaded
# ============================================================

check_hawkes_bootstrap_dependencies <- function() {
  needed <- c(
    "simulate_no_hawkes",
    "fit_hawkes_test"
  )

  missing <- needed[!vapply(needed, exists, logical(1))]

  if (length(missing) > 0) {
    stop(
      "Missing functions: ",
      paste(missing, collapse = ", "),
      "\nPlease source the simulation and fitting files first."
    )
  }

  invisible(TRUE)
}


# ============================================================
# Parametric bootstrap for H0: alpha = 0
# ============================================================

run_hawkes_parametric_bootstrap <- function(
    X,
    dt,
    mu0_null_hat,
    mu1_null_hat,
    LR_obs,
    B_boot = 100,
    verbose = TRUE
) {
  check_hawkes_bootstrap_dependencies()

  boot_results <- vector("list", B_boot)

  for (b in 1:B_boot) {
    if (verbose) {
      cat("Bootstrap repetition", b, "of", B_boot, "\n")
    }

    sim_boot <- simulate_no_hawkes(
      X = X,
      dt = dt,
      mu0 = mu0_null_hat,
      mu1 = mu1_null_hat
    )

    boot_fit <- tryCatch(
      {
        fit_hawkes_test(
          X = X,
          dN = sim_boot$dN,
          dt = dt,
          hessian = FALSE
        )
      },
      error = function(e) {
        cat("Bootstrap fit error:", conditionMessage(e), "\n")

        data.frame(
          n_events = NA,

          mu0_null = NA,
          mu1_null = NA,

          mu0_full = NA,
          mu1_full = NA,
          alpha_full = NA,
          beta_full = NA,

          loglik_null = NA,
          loglik_full = NA,

          LR = NA,
          p_chisq = NA,

          conv_null = NA,
          conv_full = NA
        )
      }
    )

    boot_fit$boot_id <- b
    boot_results[[b]] <- boot_fit
  }

  boot_results <- do.call(rbind, boot_results)

  boot_ok <- complete.cases(boot_results$LR)
  boot_LR <- boot_results$LR[boot_ok]

  p_boot <- (1 + sum(boot_LR >= LR_obs)) /
    (1 + length(boot_LR))

  crit_boot_95 <- unname(quantile(
    boot_LR,
    probs = 0.95
  ))

  crit_chisq_95 <- qchisq(
    0.95,
    df = 1
  )

  list(
    boot_results = boot_results,
    boot_LR = boot_LR,
    p_boot = p_boot,
    n_success = length(boot_LR),
    crit_boot_95 = crit_boot_95,
    crit_chisq_95 = crit_chisq_95
  )
}