# ============================================================
# Paper/code checklist for SPBDBS2025 algorithms
# ============================================================

# This script records the correspondence between paper items,
# executable scripts, and expected result files.

setwd(normalizePath("~/GitHub/biomedical-hawkes-covariates"))

dir.create("results", showWarnings = FALSE)

checklist <- data.frame(
  paper_item = c(
    "Algorithm 1: Wald test for one parameter",
    "Algorithm 2: corrected goodness-of-fit test",
    "Algorithm 3: equality Wald test",
    "Algorithm 1 calibration",
    "Algorithm 2 GOF calibration",
    "Algorithm 3 calibration",
    "Covariate-Hawkes model example",
    "Covariate Wald calibration",
    "Covariate Algorithm 2 GOF",
    "Project summary"
  ),
  script = c(
    "26_SPBDBS2025_algorithm1_wald_test.R",
    "27_SPBDBS2025_algorithm2_corrected_gof.R",
    "28_SPBDBS2025_algorithm3_equality_wald_test.R",
    "29_SPBDBS2025_wald_calibration.R",
    "30_SPBDBS2025_algorithm2_gof_calibration.R",
    "29_SPBDBS2025_wald_calibration.R",
    "31_SPBDBS2025_covariate_model_example.R",
    "32_SPBDBS2025_covariate_wald_calibration.R",
    "33_SPBDBS2025_covariate_algorithm2_gof.R",
    "99_collect_SPBDBS2025_results.R"
  ),
  result_files = c(
    "SPBDBS2025_algorithm1_wald_results.csv",
    "SPBDBS2025_algorithm2_test_summary.csv",
    "SPBDBS2025_algorithm3_equality_wald_results.csv",
    "SPBDBS2025_wald_calibration_algorithm1_summary.csv",
    "SPBDBS2025_algorithm2_gof_calibration_summary.csv",
    "SPBDBS2025_wald_calibration_algorithm3_summary.csv",
    "SPBDBS2025_covariate_model_fit_summary.csv",
    "SPBDBS2025_covariate_wald_calibration_summary.csv",
    "SPBDBS2025_covariate_algorithm2_gof_summary.csv",
    "SPBDBS2025_project_summary.csv"
  ),
  status = c(
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented_with_caveat",
    "implemented"
  ),
  comment = c(
    "Check notation and scaling of sigma_hat.",
    "Check corrected/subsampled statistic against paper.",
    "Equality contrast implemented; min p-value can be 0 numerically.",
    "Looks calibrated around 5 percent.",
    "Subsampled version appears best calibrated.",
    "Looks calibrated around 5 percent.",
    "Good recovery at T = 1000.",
    "Mean rejection close to 5 percent.",
    "Runs end-to-end. Current corrected covariate GOF is an interface-level extension of residual diagnostics, not a separate theoretical SPBDBS covariate Algorithm 2.",
    "Includes scripts 26--33."
  ),
  stringsAsFactors = FALSE
)

checklist$script_exists <- file.exists(file.path("scripts/03_SPBDBS2025_algorithms", checklist$script))
checklist$result_exists <- file.exists(file.path("results", checklist$result_files))

print(checklist)

write.csv(
  checklist,
  "results/SPBDBS2025_paper_code_checklist.csv",
  row.names = FALSE
)

cat("\nSaved:\n")
cat("- results/SPBDBS2025_paper_code_checklist.csv\n")
