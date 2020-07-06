# miceFast <a href='https://github.com/polkas/miceFast'><img src='miceFast_logo.png' align="right" height="200" /></a>

Maciej Nasinski  
GitHub:  https://github.com/polkas/miceFast

[**Check the R CRAN for more details**](https://CRAN.R-project.org/package=miceFast)

[![R build status](https://github.com/polkas/miceFast/workflows/R-CMD-check/badge.svg)](https://github.com/polkas/miceFast/actions)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/miceFast?color=brightgreen)](http://www.r-pkg.org/pkg/miceFast)
[![CRAN](http://www.r-pkg.org/badges/version/miceFast)](https://cran.r-project.org/package=miceFast)
[![codecov](https://codecov.io/gh/Polkas/miceFast/branch/master/graph/badge.svg)](https://codecov.io/gh/Polkas/miceFast)

Fast imputations under the object-oriented programming paradigm. There was used quantitative models with a closed-form solution. Thus package is based on linear algebra operations. The biggest improvement in time performance could be achieve for a calculation where a grouping variable have to be used. A single evaluation of a quantitative model for the multiple imputations is another major enhancement. Moreover there are offered a few functions built to work with popular R packages (data.table/dplyr).
A new major improvement is the one of the fastest predictive mean matching in the R world, based on pre-sorting and binary search not knn algorithms or O(N^2) loops.


Performance benchmarks (check performance_validity.R file at extdata).

[Advanced Usage - Vignette](https://cran.r-project.org/web/packages/miceFast/vignettes/miceFast-intro.html)

## Installation

```r
install.packages('miceFast')
```

or

```r
# install.packages("devtools")
devtools::install_github("polkas/miceFast")
```

## Quick Implementation

```r
library(miceFast)
library(dplyr)

set.seed(1234)
data(air_miss)

naive_fill_NA(air_miss)

#Check vigniette for an advance usage
#there is required a thorough examination

#Other packages - popular simple solutions
#Hmisc
data.frame(Map(function(x) Hmisc::impute(x,'random'), air_miss))


#mice
mice::complete(mice::mice(air_miss, printFlag = F))

```
