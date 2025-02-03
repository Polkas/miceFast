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
