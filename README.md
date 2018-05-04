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
Summing up, miceFast offer a relevant reduction of a calculations time for:

- Linear Discriminant Analysis **(x10)**
- where a grouping variable have to be used **(around x50 depending on data dimensions and number of groups and even more than x1000 although compared to data.table only a few times faster or even the same)** because data is sorted by grouping variable*
- multiple imputations is faster around **x(number of multiple imputatition)** because the core of a model is evaluated only ones.

Environment: MRO 3.4.4 Intel MKL - i7 6700HQ and 24GB DDR4 2133. MRO (Microsoft R Open) provide to R a sophisticated library for linear algebra operations so remember about that when reading a performance comparison. 

Implemented classes:

- `miceFast` (methods:`set_data`,`get_data`,`set_w`,`get_w`,`set_g`,`get_g`,`impute`,`impute_N`,`update_var`,`which_updated`,`get_model`,
`get_models`,`sort_byg`,`is_sorted_byg`,`get_index`,`vifs`,...)
- `corrData` (methods:`fill`)

The first module offers capabilities of imputations models with a closed-form solution. The main upgrade is possibility of including a grouping and/or weighting (only for linear models) variable.
The second module was made for purpose of presenting the miceFast usage and performance. It provides functionality of generating correlated data with a discrete, binomial or continuous dependent variable and continuous independent variables.

Performance benchmarks (check performance_validity.R file at extdata).

Moreover there are offered a few functions (`fill_NA`, `fill_NA_N` and `VIF`) built to work with the popular R packages such as 'data.table'.

## Installation

```r
# install.packages("devtools")
devtools::install_github("polkas/miceFast")
```

## Introduction for data.table users - using additional functions from miceFast:

Usage of `fill_NA` and `fill_NA_N`  functions from miceFast

```r
library(miceFast)
library(data.table)
library(magrittr)
```

```r
data = cbind(as.matrix(mice::nhanes),intercept=1,index=1:nrow(mice::nhanes))
data_DT = data.table(data)

# simple mean imputation - intercept at position 5
data_DT[,bmi_imp:=fill_NA(x=as.matrix(.SD),
                         model="lm_bayes",
                         posit_y=2,
                         posit_x=5)] %>% 
# there is a new variable at position 7 - bmi_imp
  .[,hyp_imp:=fill_NA(x=as.matrix(.SD),
                     model="lda",
                     posit_y=3,
                     posit_x=c(1,7)),] %>% 
  .[,chl_imp:=fill_NA_N(x=as.matrix(.SD),
                       model="lm_noise",
                       posit_y=4,
                       posit_x=c(1,7,8),
                       times=10),]

head(data_DT,3)
```

**Model with additional parameters:** - data with the grouping variable

```r
data = cbind(as.matrix(airquality[,-5]),intercept=1,index=1:nrow(airquality),
             # a numeric vector - positive values 
             weights = round(rgamma(nrow(airquality),3,3),1),
             # as.numeric is needed only for miceFast - see on next pages
             groups = airquality[,5])
data_DT = data.table(data)

# simple mean imputation - intercept at position 6
data_DT[,Ozone_imp:=fill_NA(x=as.matrix(.SD), 
                           model="lm_pred",
                           posit_y=1,
                           posit_x=c(6),w=.SD[['weights']]),by=.(groups)] %>% 
# avg of 10 multiple imputations - last posit_x equal to 9 not 10 
# because the groups variable is not included in .SD
  .[,Solar_R_imp:=fill_NA_N(as.matrix(.SD),
                           model="lm_bayes",
                           posit_y=2,
                           posit_x=c(3,4,5,6,9),w=.SD[['weights']],times=10),by=.(groups)]
head(data_DT,10)
```

## miceFast module usage:

```r
library(miceFast)
library(data.table)
library(magrittr)
```

Remember that a matrix could be build only under a one data type so factor variables have to be melted
use `model.matrix` to get numeric matrix from `data.frame` - see Tips in this document

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

#Old slow syntax
#model$update_var(4,rowMeans(sapply(1:10,function(x) 
#  model$impute("lm_bayes",4,c(1,2,3))$imputations))
#  )

model$update_var(4,model$impute_N("lm_bayes",4,c(1,2,3),10)$imputations)

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

model$update_var(2,model$impute_N("lm_bayes",2,c(1,3,4,5,6),10)$imputations)

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

model$update_var(2,model$impute_N("lm_bayes",2,c(1,3,4,5,6),10)$imputations)

#Printing data and retrieving an old order
head(cbind(model$get_data(),model$get_g(),model$get_w())[order(model$get_index()),],4)

head(airquality,4)

head(cbind(model$get_data(),model$get_g(),model$get_w()),4) #is ordered by g

head(cbind(data,groups,weights),4) #is sorted by g cause we provide data by a reference

rm(model)

```

## Tips

**matrix from data.frame**

Remember that a matrix could be build only under a one data type so factor/character variables have to be melted
use `model.matrix` to get numeric matrix from `data.frame`

```r
str(mtcars)

mtcars$cyl= factor(mtcars$cyl)
mtcars$gear= factor(mtcars$gear)
mtcars_mat = model.matrix.lm(~.,mtcars,na.action="na.pass")

str(mtcars_mat)
```

**VIF**

Variance inflation factors (VIF) measure how much the variance of the estimated regression coefficients are inflated. It helps to identify when the predictor variables are linearly related. You have to decide which variable should be delete from regression independent variables.

```r
airquality2 = airquality
airquality2$Temp2 = airquality2$Temp**2
#install.packages("car")
#car::vif(lm(Ozone ~ ., data=airquality2))

airquality2_mat = as.matrix(airquality2)
model = new(miceFast)
model$set_data(airquality2_mat)

as.vector(model$vifs(1,c(2,3,4,5,6,7)))
```
