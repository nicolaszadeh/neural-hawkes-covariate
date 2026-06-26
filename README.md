# Biomedical Hawkes Covariates

This repository contains R implementations and numerical experiments for Hawkes processes with biomedical covariates.

The project is motivated by statistical inference for event data in biomedical and neuroscience-inspired settings, where event intensities may depend both on self-excitation and on external covariate signals.

The repository includes:

- simulation tools for Hawkes and covariate-Hawkes processes;
- likelihood-based estimation routines;
- likelihood-ratio and Wald-type tests;
- goodness-of-fit diagnostics based on time-rescaling;
- Monte Carlo calibration experiments;
- reproducible scripts for the SPBDBS2025 algorithmic study.

---

## Repository structure

```text
R/
```

Core R functions for simulation, likelihoods, fitting, residual diagnostics, Wald tests, bootstrap routines, and plotting.

```text
scripts/
```

Executable scripts grouped by numerical experiment or project stage.

```text
scripts/03_SPBDBS2025_algorithms/
```

Main reproducible scripts for the SPBDBS2025 algorithmic experiments.

```text
results/
```

CSV summaries produced by the scripts.

```text
figures/
```

PDF and PNG diagnostic plots produced by the scripts.

---

## Loading the project

From the repository root, load all project functions with

```r
source("R/load_all.R")
```

Most scripts assume that the working directory is the repository root.

---

## Main SPBDBS2025 experiments

The main algorithm scripts are located in

```text
scripts/03_SPBDBS2025_algorithms/
```

They cover:

| Script | Purpose |
|---|---|
| `26_SPBDBS2025_algorithm1_wald_test.R` | Single-parameter Wald test. |
| `27_SPBDBS2025_algorithm2_corrected_gof.R` | Corrected goodness-of-fit test. |
| `28_SPBDBS2025_algorithm3_equality_wald_test.R` | Equality Wald test. |
| `29_SPBDBS2025_wald_calibration.R` | Monte Carlo calibration of Algorithms 1 and 3. |
| `30_SPBDBS2025_algorithm2_gof_calibration.R` | Monte Carlo calibration of Algorithm 2. |
| `31_SPBDBS2025_covariate_model_example.R` | Covariate-Hawkes simulation and estimation example. |
| `32_SPBDBS2025_covariate_wald_calibration.R` | Covariate-Hawkes Wald calibration. |
| `33_SPBDBS2025_covariate_algorithm2_gof.R` | Covariate-Hawkes residual goodness-of-fit diagnostics. |
| `99_collect_SPBDBS2025_results.R` | Consolidated numerical summary. |

A more detailed description is provided in

```text
README_SPBDBS2025.md
```

---

## Running the SPBDBS2025 suite

From the repository root:

```r
source("R/load_all.R")

scripts <- file.path(
  "scripts/03_SPBDBS2025_algorithms",
  c(
    "26_SPBDBS2025_algorithm1_wald_test.R",
    "27_SPBDBS2025_algorithm2_corrected_gof.R",
    "28_SPBDBS2025_algorithm3_equality_wald_test.R",
    "29_SPBDBS2025_wald_calibration.R",
    "30_SPBDBS2025_algorithm2_gof_calibration.R",
    "31_SPBDBS2025_covariate_model_example.R",
    "32_SPBDBS2025_covariate_wald_calibration.R",
    "33_SPBDBS2025_covariate_algorithm2_gof.R",
    "99_collect_SPBDBS2025_results.R"
  )
)

for (s in scripts) {
  cat("\nRunning:", s, "\n")
  source(s)
}
```

The checklist can be regenerated with

```r
source("scripts/03_SPBDBS2025_algorithms/00_SPBDBS2025_paper_code_checklist.R")
```

and the global numerical summary with

```r
source("scripts/03_SPBDBS2025_algorithms/99_collect_SPBDBS2025_results.R")
```

---

## Main outputs

The scripts write numerical summaries to

```text
results/
```

including

```text
results/SPBDBS2025_paper_code_checklist.csv
results/SPBDBS2025_project_summary.csv
```

Diagnostic plots are written to

```text
figures/
```

usually in both PDF and PNG formats.

---

## Current numerical status

The current SPBDBS2025 scripts run end-to-end.

The numerical experiments show that:

- the Wald tests are empirically calibrated near the nominal level;
- the equality Wald test behaves as expected in the tested regimes;
- the subsampled corrected goodness-of-fit diagnostic gives the best empirical calibration for Algorithm 2;
- the covariate-Hawkes maximum likelihood estimator recovers the simulated parameters well;
- the covariate-Hawkes residual diagnostics are consistent with the exponential time-rescaling principle.

The covariate goodness-of-fit script is currently marked with a caveat: the residual diagnostics run correctly and are useful empirically, but the fully justified SPBDBS Algorithm 2 correction has only been implemented for the non-covariate Hawkes model.

---

## Requirements

The code is written in R.

The scripts mainly use standard R functionality together with common plotting and data-handling packages. Project functions are loaded through

```r
source("R/load_all.R")
```

---

## Status

This repository is a research and portfolio project. It is intended to provide reproducible code for Hawkes-process inference experiments and to document the development of statistical tools for event data with biomedical covariates.
