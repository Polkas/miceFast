# miceFast

Maciej Nasinski  
GitHub:  https://github.com/polkas/miceFast

Travis badge - click on the image:

[![Build Status](https://travis-ci.org/Polkas/miceFast.svg?branch=master)](https://travis-ci.org/Polkas/miceFast) 

Object-oriented R programming by [Rcpp Module](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-modules.pdf)

R package was built under Rcpp and Armadillo for the purpose of fast multiple imputations in the face of **Big Data World** (but still fitting into the RAM).
There was used quantitative models with closed-form solution. Thus package is based on linear algebra operations.
The biggest improvement in time performance could be achieved for a calculation where a grouping variable have to be used (around x50 depending on data dimensions and number of groups and even more than x1000).
Another performance boost could be achieved for Linear Discriminant Analysis (x10).

Implemented classes:

- `miceFast` (methods:`set_data`,`set_w`,`set_g`,`impute`,`impute_force`,`which_updated`,`get_model`,`get_models`,...)
- `corrData` (methods:`fill`)

The first module offers capabilities of multiple imputations models with a closed-form solution. The main upgrade is possibility of including a grouping and/or weighting (only for linear models) variable.
The second module was made for purpose of presenting the miceFast usage and performance. It provides functionality of generating correlated data with a discrete, binomial or continuous dependent variable and continuous independent variables.

Performance benchmarks (check performance_validity.R file at extdata).

miceFast module usage:

```r
library(miceFast)
#install.packages("mice")

data = cbind(as.matrix(mice::nhanes),1)

model = new(miceFast)

model$set_data(data) #providing data by a reference

# Setting imputations directly is preferable when working with "Big Data"

data[,2] = model$impute_force("lm_pred",2,5)$imputations #Simple mean - permanent imputation at the object
data[,3] = model$impute_force("lda",3,c(1,2))$imputations # permanent imputation at the object
data[,4] = rowMeans(sapply(1:10,function(x) model$impute("lm_bayes",4,c(1,2,3))$imputations))

# Be careful with impute_force because of the permanent update at the object - check the documentation
# That is why impute_force could be used only ones for a certain column
# check which variables was updated - only inside the object
model$which_updated()

data

rm(model)

###################################
###Model with additional parameters
###################################

data = cbind(as.matrix(airquality[,-5]),1) # adding a intercept
weights = rgamma(nrow(data),3,3) # positive numeric values
groups = airquality[,5] # vector of integers

model = new(miceFast)

model$set_data(data) # providing data by a reference
model$set_w(weights)
model$set_g(groups)

#impute adapt to provided parmaters like w or g
data[,1] = model$impute_force("lm_pred",1,c(6))$imputations #Simple mean - permanent imputation at the object
data[,2] = rowMeans(sapply(1:10,function(x) model$impute("lm_bayes",2,c(1,3,4,5,6))$imputations))
data = cbind(data,groups)

data

rm(model)

```

For more complex examples check the vignette
