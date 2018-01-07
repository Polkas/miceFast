# miceFast

Maciej Nasinski  
GitHub:  https://github.com/polkas/miceFast

Travis badge - click on the image:

[![Build Status](https://travis-ci.org/Polkas/miceFast.svg?branch=master)](https://travis-ci.org/Polkas/miceFast) 

Object-oriented R programming by [Rcpp Module](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-modules.pdf)

R package was built under Rcpp and Armadillo for the purpose of fast multiple imputations.
There was used quantitative models with closed-form solution. Thus package is based on linear algebra operations.
The biggest improvement in time performance could be achieved for a calculation where a grouping variable have to be used (x40).
Another performance boost could be achived for Linear Discriminant Analysis (x5).

Implemented classes:

- `miceFast` (methods:`impute`,`get_model`,`get_models`,`is_vars_updated`,`get_index_NA_R`,`get_index_full_R`)
- `corrData` (methods:`fill`)

The first module offers capabilities of multiple imputations models with a closed-form solution. The main upgrade is possibility of including a grouping and/or weighting (only for linear models) variable.
The second module was made for purpose of presenting the miceFast usage and performance. It provides functionality of generating correlated data with a discrete, binomial or continuous dependent variable and continous independent variables.

Performance benchmarks (check performance_validity.R file at extdata).

miceFast module usage:

```r
library(miceFast)
library(mice)

model = new(miceFast,as.matrix(nhanes))

cbind(nhanes[,1],
      model$impute("lm_pred",2,1,TRUE)$imputations,
      model$impute("lda",3,c(1,2),TRUE)$imputations,
      model$impute("lm_bayes",4,c(1,2,3),FALSE)$imputations)

```

For more complex examples check the vigniette
