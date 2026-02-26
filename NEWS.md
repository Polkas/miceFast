# miceFast 0.9.1

## Bug fixes

* PMM returned predicted values instead of observed values (C++): The `pmm` model returned predicted $\hat{y}$ for missing rows instead of the nearest observed $y$ values. Now it follows Little and Rubin (2002).
* PMM with character/factor variables (R): `fill_NA_N()` with `model = "pmm"` and a character dependent variable failed because it attempted `as.numeric()` on non-numeric strings, producing all NAs.
* Character dependent variable with lm models: `fill_NA()` and `fill_NA_N()` with `model = "lm_pred"`, `"lm_bayes"`, or `"lm_noise"` silently returned all NAs when the dependent variable was character with non-numeric labels (e.g., `"apple"`, `"banana"`).

## Documentation

* README: added sequential-chain MI examples (dplyr and data.table) showing how to impute multiple variables and pool with Rubin's rules.
* Introduction vignette: added full imputation workflow with sequential ordering (impute variables whose predictors are complete first), FCS (chained equations) section with data.table example, and PMM note for the OOP interface.
* MI vignette: expanded Rubin's rules derivations, added PMM MI example using the OOP interface, expanded "Important caveat" section with OOP and data.table FCS code snippets for non-monotone patterns.
* Documented PMM as a proper MI method throughout vignettes and README.
* Improved prose throughout vignettes and README.

## Tests

* Added 20 PMM-specific tests (`test-pmm.R`): observed-value returns, factor/character support, weighted PMM, grouped data.table, reproducibility, stochasticity.
* Added 31 FCS tests (`test-fcs.R`): data.table, data.frame, and OOP FCS helpers; joint-missingness handling; MI+pool workflow; comparison with `mice` (pooled estimates and imputed means).
* Added tests for character dependent variables with non-numeric labels across all models and data types.
* Test suite expanded from 243 to 311 tests.

# miceFast 0.9.0

Kota Hattori, thank you for your feedback and for motivating me for this deep update.

## New features

* `pool()` function for combining results from multiply imputed datasets (Rubin's rules, Barnard-Rubin df adjustment). Works with `lm`, `glm`, and other models that support `coef()` and `vcov()`. Validated against `mice`.
* `print` and `summary` methods for pooled results.

## Bug fixes

* fixed residual variance estimator in `lm_noise` and `lm_bayes` stochastic models: divisor changed from `n-p-1` to `n-p`, where `p` already counts the intercept column supplied by the user. The previous formula over-corrected by one degree of freedom.

## Documentation and internals

* new vignette on missing data mechanisms (MCAR/MAR/MNAR) and MI workflows.
* refactored introduction vignette with `pool()` examples.
* improved README with MI section and benchmark table.
* test suite for `pool()`, including comparison against `mice::pool()`.
* new weighted regression validation test against `lm.wfit()`.
* refactored C++ source code for clarity.
* fixed typos in error messages and documentation.
* regenerated performance benchmarks on R 4.4.3, macOS M3 Pro.

# miceFast 0.8.5

* cran related update, `OMP_THREAD_LIMIT`.

# miceFast 0.8.4

* fixed CRAN Notes.
* style the cpp code.
* VIF() should be more stable.

# miceFast 0.8.2

* simplified `naive_fill_NA`, It is a regular sampling imputation now.
* Fixed `dontrun` examples.
* replace `ggplot2::aes_string` with  `ggplot2::aes`, as the former is depreciated.
* regenerate performance benchmarks on R 4.2.1.
* styler over the code.
* improve documentation.

# miceFast 0.8.1

* `tinyverse` world, less dependencies.
* fixed imputations for character variables under linear models.
* speed up the `pmm` model.
* more tests, higher `covr`.
* rerun performance tests.

# miceFast 0.7.1

* update URL inside README.

# miceFast 0.7.0

* improve coverage.
* use drop = FALSE when subsetting the data.frame
* healthy DESCRIPTION file, fix spaces.
* more input validation.

# miceFast 0.6.8

* update broken vignette links

# miceFast 0.6.6

* solve broken UpSetR::upset reference links

# miceFast 0.6.5

* upset_NA based on UpSetR::upset plot function
* compare_imp plot function
* new logo
* remove times argument

# miceFast 0.6.2

* R CRAN r-oldrel-windows-ix86+x86_64 problems

# miceFast 0.6.1

* lifecycle problems

# miceFast 0.6.0

* fill_NA_N has a new model which is pmm - predictive mean matching
* fast PMM - presorting and binary search
* naive_fill_NA - auto function for data.frames - bayes mean and lda
* ridge argument for lm models - adding small disturbance to diag of X'X
* lm_bayes provide more disturbance
* new tests
* codecov

# miceFast 0.5.1

* remove old urls form vignettes

# miceFast 0.5.0

* providing a more comfortable environment for data.table/dplyr users
* expand vignette and documentation
* updated performance benchmarks
* fix a glitch - e.g. lack of correct warning for a lda model with zero variance variables

# miceFast 0.2.1-3

* data.table problem - jump to R 3.5.0
* valgrind -  a lot of optimizations - problem with arma::exp and arma::randn
* optimize a lot of code
* methods/functions  resistant to glitches

# miceFast 0.2.0

* fix imputations with a grouping variable - error if there is precisly one NA at any group
* add data.table to benchmarks - model with a grouping variable
* add R functions (`fill_NA_N`,`fill_NA`,`VIF`) which could be used by a data.table user

# miceFast 0.1.0

* add `impute_N` method - optimized multiple imputations
* add `vif` method -  Variance inflation factors

# miceFast 0.0.3

* vignette,readme,description,todo

# miceFast 0.0.2

* adjust to solaris
* reference - set a grouping variable by a reference but as a numeric vector - integer vector do not work (randomly lost pointer)
