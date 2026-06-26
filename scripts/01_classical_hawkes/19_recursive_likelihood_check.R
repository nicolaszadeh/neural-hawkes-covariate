# 19_recursive_likelihood_check.R
#
# Refactored version.
# Checks that naive and recursive Hawkes likelihoods agree,
# and compares their computational cost.

set.seed(19)

source("R/load_all.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(ggplot2)

T_end <- 500

mu_true <- 0.8
alpha_true <- 0.35
beta_true <- 1.2
branching_true <- alpha_true / beta_true

events <- simulate_hawkes(
  T_end = T_end,
  mu = mu_true,
  alpha = alpha_true,
  beta = beta_true
)

cat("Number of events:", length(events), "\n")

hawkes_loglik_naive <- function(events, T_end, mu, alpha, beta) {
  n <- length(events)

  if (n == 0) {
    return(-mu * T_end)
  }

  lambda_events <- numeric(n)

  for (i in seq_len(n)) {
    if (i == 1) {
      R_i <- 0
    } else {
      past_events <- events[seq_len(i - 1)]
      R_i <- sum(exp(-beta * (events[i] - past_events)))
    }

    lambda_events[i] <- mu + alpha * R_i

    if (!is.finite(lambda_events[i]) ||
        lambda_events[i] <= 0) {
      return(-Inf)
    }
  }

  integral <- mu * T_end +
    sum(alpha / beta * (1 - exp(-beta * (T_end - events))))

  sum(log(lambda_events)) - integral
}

loglik_naive <- hawkes_loglik_naive(
  events = events,
  T_end = T_end,
  mu = mu_true,
  alpha = alpha_true,
  beta = beta_true
)

loglik_recursive <- hawkes_loglik(
  events = events,
  T_end = T_end,
  mu = mu_true,
  alpha = alpha_true,
  beta = beta_true
)

difference <- abs(loglik_naive - loglik_recursive)

cat("Naive log-likelihood:", loglik_naive, "\n")
cat("Recursive log-likelihood:", loglik_recursive, "\n")
cat("Absolute difference:", difference, "\n")

n_repeat <- 200

event_sizes <- c(50, 100, 200, 400, 600)
event_sizes <- event_sizes[event_sizes <= length(events)]

timing_rows <- list()

for (m in event_sizes) {
  sub_events <- events[seq_len(m)]
  sub_T <- max(sub_events)

  ll_naive_m <- hawkes_loglik_naive(
    sub_events,
    sub_T,
    mu_true,
    alpha_true,
    beta_true
  )

  ll_recursive_m <- hawkes_loglik(
    sub_events,
    sub_T,
    mu_true,
    alpha_true,
    beta_true
  )

  naive_time <- system.time({
    for (k in seq_len(n_repeat)) {
      hawkes_loglik_naive(
        sub_events,
        sub_T,
        mu_true,
        alpha_true,
        beta_true
      )
    }
  })["elapsed"]

  recursive_time <- system.time({
    for (k in seq_len(n_repeat)) {
      hawkes_loglik(
        sub_events,
        sub_T,
        mu_true,
        alpha_true,
        beta_true
      )
    }
  })["elapsed"]

  naive_time <- as.numeric(naive_time)
  recursive_time <- as.numeric(recursive_time)

  timing_rows[[length(timing_rows) + 1]] <- data.frame(
    n_events = m,
    n_repeat = n_repeat,
    loglik_naive = ll_naive_m,
    loglik_recursive = ll_recursive_m,
    abs_difference = abs(ll_naive_m - ll_recursive_m),
    naive_time_total = naive_time,
    recursive_time_total = recursive_time,
    naive_time_per_eval = naive_time / n_repeat,
    recursive_time_per_eval = recursive_time / n_repeat,
    speedup = naive_time / max(recursive_time, .Machine$double.eps)
  )

  cat(
    "n =", m,
    "- diff =", signif(abs(ll_naive_m - ll_recursive_m), 3),
    "- naive total =", round(naive_time, 4),
    "s - recursive total =", round(recursive_time, 4),
    "s - speedup =",
    round(naive_time / max(recursive_time, .Machine$double.eps), 2),
    "\n"
  )
}

timing_df <- do.call(rbind, timing_rows)

summary_df <- data.frame(
  T_end = T_end,
  n_events = length(events),
  mu_true = mu_true,
  alpha_true = alpha_true,
  beta_true = beta_true,
  branching_true = branching_true,
  loglik_naive = loglik_naive,
  loglik_recursive = loglik_recursive,
  abs_difference = difference
)

print(summary_df)
print(timing_df)

write.csv(
  summary_df,
  "results/recursive_likelihood_check_summary.csv",
  row.names = FALSE
)

write.csv(
  timing_df,
  "results/recursive_likelihood_timing.csv",
  row.names = FALSE
)

timing_long <- rbind(
  data.frame(
    n_events = timing_df$n_events,
    method = "naive",
    elapsed_time = timing_df$naive_time_per_eval
  ),
  data.frame(
    n_events = timing_df$n_events,
    method = "recursive",
    elapsed_time = timing_df$recursive_time_per_eval
  )
)

p_time <- ggplot(
  timing_long,
  aes(
    x = n_events,
    y = elapsed_time,
    linetype = method
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    title = "Naive versus recursive Hawkes likelihood",
    subtitle = "Recursive evaluation avoids summing over the full past",
    x = "number of events",
    y = "elapsed time per evaluation in seconds"
  ) +
  theme_minimal()

ggsave(
  "figures/recursive_likelihood_timing.pdf",
  p_time,
  width = 7,
  height = 4
)

ggsave(
  "figures/recursive_likelihood_timing.png",
  p_time,
  width = 7,
  height = 4,
  dpi = 300
)

p_diff <- ggplot(
  timing_df,
  aes(
    x = n_events,
    y = abs_difference
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    title = "Naive and recursive likelihood agree",
    x = "number of events",
    y = "absolute difference"
  ) +
  theme_minimal()

ggsave(
  "figures/recursive_likelihood_difference.pdf",
  p_diff,
  width = 7,
  height = 4
)

ggsave(
  "figures/recursive_likelihood_difference.png",
  p_diff,
  width = 7,
  height = 4,
  dpi = 300
)

cat("\nSaved:\n")
cat("- results/recursive_likelihood_check_summary.csv\n")
cat("- results/recursive_likelihood_timing.csv\n")
cat("- figures/recursive_likelihood_timing.pdf/png\n")
cat("- figures/recursive_likelihood_difference.pdf/png\n")
