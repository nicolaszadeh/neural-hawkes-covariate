# 14_covariate_baseline_bootstrap_pvalue.R

source("R/hawkes_covariate_baseline_1d.R")

cov_boot <- run_covariate_baseline_bootstrap_pvalue_demo_1d(
  B = 100
)