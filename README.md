# miceFast

Maciej Nasinski  
GitHub:  https://github.com/polkas/miceFast

Travis badge - click on the image:

[![Build Status](https://travis-ci.org/Polkas/miceFast.svg?branch=master)](https://travis-ci.org/Polkas/miceFast) 
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/miceFast?color=brightgreen)](http://www.r-pkg.org/pkg/miceFast)
[![CRAN](http://www.r-pkg.org/badges/version/miceFast)](https://cran.r-project.org/package=miceFast)

Object-oriented R programming by [Rcpp Module](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-modules.pdf)

R package was built under Rcpp and Armadillo for the purpose of fast imputations in the face of **Big Data World** (but still fitting into the RAM).
There was used quantitative models with closed-form solution. Thus package is based on linear algebra operations.
The biggest improvement in time performance could be achieved for a calculation where a grouping variable have to be used (around x50 depending on data dimensions and number of groups and even more than x1000).
Another performance boost could be achieved for Linear Discriminant Analysis (x10).

Implemented classes:

- `miceFast` (methods:`set_data`,`get_data`,`set_w`,`get_w`,`set_g`,`get_g`,`impute`,`update_var`,`which_updated`,`get_model`,
                      `get_models`,`sort_byg`,`is_sorted_byg`,`get_index`,`vifs`,...)
- `corrData` (methods:`fill`)

The first module offers capabilities of imputations models with a closed-form solution. The main upgrade is possibility of including a grouping and/or weighting (only for linear models) variable.
The second module was made for purpose of presenting the miceFast usage and performance. It provides functionality of generating correlated data with a discrete, binomial or continuous dependent variable and continuous independent variables.

Performance benchmarks (check performance_validity.R file at extdata).

## Installation

```r
# install.packages("devtools")
devtools::install_github("polkas/miceFast")
```

## miceFast module usage:

Remeber that a matrix could be build only under a one data type so factor variables have to be melted
use `model.matrix` to get numeric matrix from `data.frame`

```r
#install.packages("mice")
data = cbind(as.matrix(mice::nhanes),intercept=1,index=1:nrow(mice::nhanes))
model = new(miceFast)
model$set_data(data) #providing data by a reference

model$update_var(2,model$impute("lm_pred",2,5)$imputations)
#OR not recommended
#data[,2] = model$impute("lm_pred",2,5)$imputations
#model$set_data(data) #Updating the object

model$update_var(3,model$impute("lda",3,c(1,2))$imputations) 
model$update_var(4,rowMeans(sapply(1:10,function(x) 
  model$impute("lm_bayes",4,c(1,2,3))$imputations))
  )

#When working with 'Big Data'
#it is recommended to occasionally manually invoke a garbage collector `gc()`

# Be careful with `update_var` because of the permanent update at the object and data
# That is why `update_var` could be used only ones for a certain column
# check which variables was updated - inside the object
model$which_updated()
head(model$get_data(),4)
head(data,4)
head(mice::nhanes,4)
rm(model)

########################################################################
###Model with additional parameters: - sorted by the grouping variable
########################################################################

data = cbind(as.matrix(airquality[,-5]),intercept=1,index=1:nrow(airquality))
weights = rgamma(nrow(data),3,3) # a numeric vector - positive values
groups = as.numeric(airquality[,5]) # a numeric vector not integers - positive values - sorted increasingly

model = new(miceFast)
model$set_data(data) # providing by a reference
model$set_w(weights) # providing by a reference
model$set_g(groups)  # providing by a reference

#impute adapt to provided parmaters like w or g
#Simple mean - permanent imputation at the object and data
model$update_var(1,model$impute("lm_pred",1,c(6))$imputations)

model$update_var(2,rowMeans(sapply(1:10,function(x) 
  model$impute("lm_bayes",2,c(1,3,4,5,6))$imputations))
  )
#Printing data and retrieving an old order
head(cbind(model$get_data(),model$get_g(),model$get_w())[order(model$get_index()),],4)

head(airquality,4)

head(cbind(model$get_data(),model$get_g(),model$get_w()),4)

head(cbind(data,groups,weights),4)

rm(model)

############################################################################
###Model with additional parameters:** - data not sorted by the grouping variable
############################################################################

data = cbind(as.matrix(airquality[,-5]),intercept=1,index=1:nrow(airquality))
weights = rgamma(nrow(data),3,3) # a numeric vector - positive values
#groups = as.numeric(airquality[,5]) # a numeric vector not integers - positive values
groups = as.numeric(sample(1:3,nrow(data),replace=T)) # a numeric vector not integers - positive values

model = new(miceFast)
model$set_data(data) # providing by a reference
model$set_w(weights) # providing by a reference
model$set_g(groups)  # providing by a reference
#impute adapt to provided parmaters like w or g
#Warning - if data is not sorted increasingly by the g then it would be done automatically 
#during a first imputation
#Simple mean - permanent imputation at the object and data
model$update_var(1,model$impute("lm_pred",1,c(6))$imputations)

model$update_var(2,rowMeans(sapply(1:10,function(x) 
  model$impute("lm_bayes",2,c(1,3,4,5,6))$imputations))
  )
#Printing data and retrieving an old order
head(cbind(model$get_data(),model$get_g(),model$get_w())[order(model$get_index()),],4)

head(airquality,4)

head(cbind(model$get_data(),model$get_g(),model$get_w()),4) #is ordered by g

head(cbind(data,groups,weights),4) #is sorted by g cause we provide data by a reference

rm(model)

```

## Tips

**matrix from data.frame**

Remeber that a matrix could be build only under a one data type so factor variables have to be melted
use `model.matrix` to get numeric matrix from `data.frame`

```r
str(mtcars)

mtcars$cyl= factor(mtcars$cyl)
mtcars$gear= factor(mtcars$gear)
mtcars_mat = model.matrix(~.,mtcars)

str(mtcars_mat)
```

**VIF**

Variance inflation factors (VIF) measure how much the variance of the estimated regression coefficients are inflated. It helps to indentify when the predictor variables are linearly related. You have to decide which variable should be delete from regression independent variables.

```r
airquality2 = airquality
airquality2$Temp2 = airquality2$Temp**2
#install.packages("car")
car::vif(lm(Ozone ~ ., data=airquality2))

airquality2_mat = as.matrix(airquality2)
model = new(miceFast)
model$set_data(airquality2_mat)

as.vector(model$vifs(1,c(2,3,4,5,6,7)))
```
