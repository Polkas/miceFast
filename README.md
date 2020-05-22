# miceFast

Maciej Nasinski  
GitHub:  https://github.com/polkas/miceFast

[**Check the R CRAN for more details**](https://CRAN.R-project.org/package=miceFast)

[![Build Status](https://travis-ci.org/Polkas/miceFast.svg?branch=master)](https://travis-ci.org/Polkas/miceFast) 
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/miceFast?color=brightgreen)](http://www.r-pkg.org/pkg/miceFast)
[![CRAN](http://www.r-pkg.org/badges/version/miceFast)](https://cran.r-project.org/package=miceFast)

Fast imputations under the object-oriented programming paradigm.  
Moreover there are offered a few functions built to work with popular R packages (data.table/dplyr).

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
set.seed(1234)
data(air_miss)

naive_fill_NA(air_miss)

#Check vigniette for an advance usage
#there is required a thorough examination

#Other packages - popular simple solutions
#Hmisc
data.frame(Map(function(x) Hmisc::impute(x,'random'),air_miss))
#mice
mice::complete(mice::mice(air_miss,printFlag = F))

```

