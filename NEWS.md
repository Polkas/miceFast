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
