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

- Linear Discriminant Analysis **(x5)**
- where a grouping variable have to be used **(around x10 depending on data dimensions and number of groups and even more than x100 although compared to data.table only a few times faster or even the same)** because data is sorted by grouping variable
- multiple imputations is faster around **x(number of multiple imputations)** because the core of a model is evaluated only ones.

Environment: MRO 3.5.3 Intel MKL - i7 6700HQ and 24GB DDR4 2133. MRO (Microsoft R Open) provide to R a sophisticated library for linear algebra operations so remember about that when reading a performance comparison - Date 2019-08-07. 

Example:

![Performance Summary](./inst/extdata/images/g_summary.png)

Implemented classes:

- `miceFast` (methods:`set_data`,`get_data`,`set_w`,`get_w`,`set_g`,`get_g`,`impute`,`impute_N`,`update_var`,`which_updated`,`get_model`,
`get_models`,`sort_byg`,`is_sorted_byg`,`get_index`,`vifs`,...)
- `corrData` (methods:`fill`)

The first module offers capabilities of imputations models with a closed-form solution. The main upgrade is possibility of including a grouping and/or weighting (only for linear models) variable.
The second module was made for purpose of presenting the miceFast usage and performance. It provides functionality of generating correlated data with a discrete, binomial or continuous dependent variable and continuous independent variables.

Performance benchmarks (check performance_validity.R file at extdata).

Moreover there are offered a few functions (`fill_NA`, `fill_NA_N` and `VIF`) built to work with the popular R packages such as 'data.table'.
This functions should be resistant to glitches from an user activity perspective and a data structure.

## Installation

```r
# install.packages("devtools")
devtools::install_github("polkas/miceFast")
```

## Introduction for data.table/dplyr users - using additional functions from miceFast:

```
#install.packages('pacman')
pacman::p_load(miceFast,data.table,magrittr,mice,car,dplyr)
```

Usage of `fill_NA` and `fill_NA_N` functions from miceFast - this functions should be resistant to glitches from an user activity perspective and a data structure.

### Intro:data.table

#### Working with names

```{r,echo=TRUE}
data = cbind(as.matrix(mice::nhanes),intercept=1,index=1:nrow(mice::nhanes))
data = do.call(rbind,replicate(10,data,simplify = F))
data_df = as.data.frame(data)
data_DT = data.table(data)

data_DT[,bmi_imp:=fill_NA(x=.SD,
                         model="lm_pred",
                         posit_y='bmi',
                         posit_x='intercept')] %>% 
  .[,hyp_imp:=fill_NA(x=.SD,
                     model="lda",
                     posit_y='hyp',
                     posit_x=c('age','bmi_imp')),] %>% 
  .[,chl_imp:=fill_NA_N(x=.SD,
                       model="lm_noise",
                       posit_y='chl',
                       posit_x=c('age','bmi_imp','hyp_imp'),
                       times=10),]

head(data_DT,2)
```

#### Working with positions

```{r,echo=TRUE}
data_DT[,bmi_imp:=fill_NA(x=.SD,
                         model="lm_pred",
                         posit_y=2,
                         posit_x=5)] %>% 
# there is a new variable at position 7 - bmi_imp
  .[,hyp_imp:=fill_NA(x=.SD,
                     model="lda",
                     posit_y=3,
                     posit_x=c(1,7)),] %>% 
  .[,chl_imp:=fill_NA_N(x=.SD,
                       model="lm_noise",
                       posit_y=4,
                       posit_x=c(1,7,8),
                       times=10),]

head(data_DT,2)
```

### Intro:dplyr

#### Working with names

```{r,echo=TRUE}
data_df = data_df %>% mutate(bmi_imp=fill_NA(x=.,
                         model="lm_pred",
                         posit_y='bmi',
                         posit_x='intercept')) %>% 
  mutate(hyp_imp=fill_NA(x=.,
                     model="lda",
                     posit_y='hyp',
                     posit_x=c('age','bmi_imp'))) %>% 
  mutate(chl_imp=fill_NA_N(x=.,
                       model="lm_noise",
                       posit_y='chl',
                       posit_x=c('age','bmi_imp','hyp_imp'),
                       times=10))

head(data_df,2)
```


#### Working with positions

```{r,echo=TRUE}
data_df = data_df %>% mutate(bmi_imp=fill_NA(x=.,
                         model="lm_pred",
                         posit_y=2,
                         posit_x=5)) %>% 
# there is a new variable at position 7 - bmi_imp
  mutate(hyp_imp=fill_NA(x=.,
                     model="lda",
                     posit_y=3,
                     posit_x=c(1,7))) %>% 
  mutate(chl_imp=fill_NA_N(x=.,
                       model="lm_noise",
                       posit_y=4,
                       posit_x=c(1,7,8),
                       times=10))

head(data_df,2)
```

**Model with additional parameters:** - data with the grouping/weighting variable - data.table recommended

```{r,echo=TRUE}
data = cbind(as.matrix(airquality[,-5]),Intercept=1,index=1:nrow(airquality),
             # a numeric vector - positive values 
             weights = round(rgamma(nrow(airquality),3,3),1),
             groups = airquality[,5])

data = do.call(rbind,replicate(10,data,simplify = F))
data_DT = data.table(data)
```

#### Working with names

```{r,echo=TRUE}
data_DT[,Ozone_imp:=fill_NA(x=.SD, 
                           model="lm_pred",
                           posit_y='Ozone',
                           posit_x='Intercept',w=.SD[['weights']]),by=.(groups)] %>% 
  .[,Solar_R_imp:=fill_NA_N(.SD,
                           model="lm_bayes",
                           posit_y='Solar.R',
                           posit_x=c('Wind','Temp','Day','Intercept','Ozone_imp'),
                           w=.SD[['weights']],
                           times=10),by=.(groups)]

data_DT[which(is.na(data_DT[,1]))[1],]
```

#### Working with positions

```{r,echo=TRUE}
# simple mean imputation - intercept at position 6
data_DT[,Ozone_imp:=fill_NA(x=.SD, 
                           model="lm_pred",
                           posit_y=1,
                           posit_x=c(6),w=.SD[['weights']]),by=.(groups)] %>% 
# avg of 10 multiple imputations - last posit_x equal to 9 not 10 
# because the groups variable is not included in .SD
  .[,Solar_R_imp:=fill_NA_N(.SD,
                           model="lm_bayes",
                           posit_y=2,
                           posit_x=c(3,4,5,6,9),
                           w=.SD[['weights']],
                           times=10),by=.(groups)]

data_DT[which(is.na(data_DT[,1]))[1],]

```


**VIF**

Variance inflation factors (VIF) measure how much the variance of the estimated regression coefficients are inflated. It helps to identify when the predictor variables are linearly related. You have to decide which variable should be delete from regression independent variables.

```r
airquality2 = airquality
airquality2$Temp2 = airquality2$Temp**2
airquality2$Month = factor(airquality2$Month)

car::vif(lm(Ozone ~ ., data=airquality2))
```

```{echo=TRUE}
data_DT = data.table(airquality2)
data_DT[,.(vifs=VIF(x=.SD,
                posit_y='Ozone',
                posit_x=c('Solar.R','Wind','Temp','Month','Day','Temp2'),correct=FALSE))][['vifs']]

data_DT[,.(vifs=VIF(x=.SD,
                posit_y=1,
                posit_x=c(2,3,4,5,6,7),correct = TRUE))][['vifs']]
```
