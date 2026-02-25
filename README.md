# miceFast <a href='https://github.com/polkas/miceFast'><img src='man/figures/miceFast_logo.png' align="right" width="200" /></a>

**Author**: Maciej Nasinski  

[**Check the miceFast website for more details**](https://polkas.github.io/miceFast/index.html)

[![R build status](https://github.com/polkas/miceFast/workflows/R-CMD-check/badge.svg)](https://github.com/polkas/miceFast/actions)
[![CRAN](http://www.r-pkg.org/badges/version/miceFast)](https://cran.r-project.org/package=miceFast)
[![codecov](https://codecov.io/gh/Polkas/miceFast/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Polkas/miceFast)
[![Dependencies](https://tinyverse.netlify.app/badge/miceFast)](https://cran.r-project.org/package=miceFast)

## Overview

**miceFast** provides fast methods for imputing missing data, leveraging an object-oriented programming paradigm and optimized linear algebra routines.  
The package includes convenient helper functions compatible with **data.table**, **dplyr**, and other popular R packages.

Major speed improvements occur when:  
- Using a **grouping variable**, where the data is automatically sorted by group, significantly reducing computation time.
- Performing **multiple imputations**, by evaluating the underlying quantitative model only once for multiple draws.
- Running **Predictive Mean Matching (PMM)**, thanks to presorting and binary search.

For performance details, see `performance_validity.R` in the `extdata` folder.

**Vignettes:**

- [Introduction and Advanced Usage](https://polkas.github.io/miceFast/articles/miceFast-intro.html)
- [Missing Data Mechanisms and Multiple Imputation](https://polkas.github.io/miceFast/articles/missing-data-and-imputation.html)

## Multiple Imputation and mice

[mice](https://cran.r-project.org/package=mice) implements the full MI pipeline (impute, analyze, pool). **miceFast** focuses on the computationally expensive part: fitting the imputation models. Two usage modes:

1. **MI with Rubin's rules** — call `fill_NA()` with a stochastic model (`lm_bayes`, `lm_noise`, or `lda` with a random `ridge`) in a loop to create *m* completed datasets, then `pool()` the fitted models.

2. **Single-dataset averaging** — `fill_NA_N()` returns the mean of *k* draws per missing value. Handy for exploration, but not for Rubin's rules (between-imputation variance is lost).

See the [MI vignette](https://polkas.github.io/miceFast/articles/missing-data-and-imputation.html) for worked examples.

## Installation

You can install **miceFast** from CRAN:
```r
install.packages("miceFast")
```
Or install the development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("polkas/miceFast")
```

## Quick Example

```r
library(miceFast)
library(dplyr)

set.seed(1234)
data(air_miss)

# Visualize the NA structure
upset_NA(air_miss, 6)

# Naive imputation (quick, but biased — see ?naive_fill_NA)
naive_fill_NA(air_miss)

# Model-based single imputation with fill_NA
air_miss %>%
  mutate(Ozone_imp = fill_NA(
    x = ., model = "lm_bayes",
    posit_y = "Ozone", posit_x = c("Solar.R", "Wind", "Temp")
  ))

# Proper MI: impute m times, fit models, pool with Rubin's rules
completed <- lapply(1:5, function(i) {
  air_miss %>%
    mutate(Ozone_imp = fill_NA(
      x = ., model = "lm_bayes",
      posit_y = "Ozone", posit_x = c("Solar.R", "Wind", "Temp")
    ))
})
fits <- lapply(completed, function(d) lm(Ozone_imp ~ Wind + Temp, data = d))
pool(fits)
```

See the [Introduction vignette](https://polkas.github.io/miceFast/articles/miceFast-intro.html) for grouped imputation, data.table syntax, the OOP interface, and more.

---

## Key Features

- **Object-Oriented Interface** via `miceFast` objects (Rcpp modules).
- **Convenient Helpers**:  
  - `fill_NA()`: Single imputation (`lda`, `lm_pred`, `lm_bayes`, `lm_noise`).  
  - `fill_NA_N()`: Multiple imputations (`pmm`, `lm_bayes`, `lm_noise`).  
  - `pool()`: Pool multiply imputed results using Rubin's rules.  
  - `VIF()`: Variance Inflation Factor calculations.  
  - `naive_fill_NA()`: Automatic naive imputations.  
  - `compare_imp()`: Compare original vs. imputed values.  
  - `upset_NA()`: Visualize NA structure using [UpSetR](https://cran.r-project.org/package=UpSetR).

**Quick Reference Table**:

| Function        | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `new(miceFast)` | Creates an OOP instance with numerous imputation methods (see the vignette). |
| `fill_NA()`     | Single imputation: `lda`, `lm_pred`, `lm_bayes`, `lm_noise`.                   |
| `fill_NA_N()`   | Averaged multiple imputations (mean of N draws): `pmm`, `lm_bayes`, `lm_noise`. |
| `pool()`        | Pools estimates from *m* imputed datasets using Rubin's rules. Works with any model that has `coef()` and `vcov()`. |
| `VIF()`         | Computes Variance Inflation Factors.                                         |
| `naive_fill_NA()` | Performs automatic, naive imputations.                                     |
| `compare_imp()` | Compares imputations vs. original data.                                      |
| `upset_NA()`    | Visualizes NA structure using an UpSet plot.                                 |

---

## Performance Highlights

Median timings on 100k rows, 10 variables, 100 groups (R 4.4.3, macOS M3 Pro, [optimized BLAS/LAPACK](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Which-BLAS-is-used-and-how-can-it-be-changed_003f)):

Imputation quality (SSE) is comparable to mice across all models.

![](man/figures/g_summary.png)

Full benchmark script: `inst/extdata/performance_validity.R`.
