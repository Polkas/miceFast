# 0.6.2

* R CRAN r-oldrel-windows-ix86+x86_64 problems

# 0.6.1

* lifecycle problems

# 0.6.0

* fill_NA_N has a new model which is pmm - predictive mean matching
* fast PMM - presorting and binary search
* naive_fill_NA - auto function for data.frames - bayes mean and lda
* ridge argument for lm models - adding small disturbance to diag of X'X
* lm_bayes provide more disturbance
* new tests
* codecov

# 0.5.1

* remove old urls form vignettes

# 0.5.0

* providing a more comfortable environment for data.table/dplyr users
* expand vignette and documentation
* updated performance benchmarks
* fix a glitch - e.g. lack of correct warning for a lda model with zero variance variables

# 0.2.1-3

* data.table problem - jump to R 3.5.0
* valgrind -  a lot of optimizations - problem with arma::exp and arma::randn
* optimize a lot of code
* methods/functions  resistant to glitches

## 0.2.0

* fix imputations with a grouping variable - error if there is precisly one NA at any group
* add data.table to benchmarks - model with a grouping variable
* add R functions (`fill_NA_N`,`fill_NA`,`VIF`) which could be used by a data.table user

### 0.1.0

* add `impute_N` method - optimized multiple imputations
* add `vif` method -  Variance inflation factors

### 0.0.3

* vignette,readme,description,todo

### 0.0.2

* adjust to solaris
* reference - set a grouping variable by a reference but as a numeric vector - integer vector do not work (randomly lost pointer)
