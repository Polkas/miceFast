# miceFast

Maciej Nasinski  
GitHub:  https://github.com/polkas/miceFast

[**Check the vignette for more details**](https://CRAN.R-project.org/package=miceFast)

Travis badge - click on the image:

[![Build Status](https://travis-ci.org/Polkas/miceFast.svg?branch=master)](https://travis-ci.org/Polkas/miceFast) 
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/miceFast?color=brightgreen)](http://www.r-pkg.org/pkg/miceFast)
[![CRAN](http://www.r-pkg.org/badges/version/miceFast)](https://cran.r-project.org/package=miceFast)

Fast imputations under the object-oriented programming paradigm.
Moreover there are offered a few functions built to work with popular R packages (data.table/dplyr).

Summing up, miceFast offer a relevant reduction of a calculations time for:

- Linear Discriminant Analysis **(x5)**
- where a grouping variable have to be used **(around x10 depending on data dimensions and number of groups and even more than x100 although compared to data.table only a few times faster or even the same)** because data is sorted by grouping variable
- multiple imputations is faster around **x(number of multiple imputations)** because the core of a model is evaluated only ones
- Variance inflation factors (VIF) **(x5)** because the unnecessary linear regression is not evaluated - we need only inverse of X'X

Example (2019-08-07):

![Performance Summary](./inst/extdata/images/g_summary.png)

Performance benchmarks (check performance_validity.R file at extdata).


## Installation

```r
# install.packages("devtools")
devtools::install_github("polkas/miceFast")
```

## Introduction for data.table users - using additional functions from miceFast:

```r
#install.packages('pacman')
pacman::p_load(miceFast,data.table,magrittr,mice,car,dplyr,ggplot2)
```

Usage of `fill_NA` and `fill_NA_N` functions from miceFast - this functions should be resistant to glitches from an user activity perspective and a data structure.

### Fast Intro

```r
data = cbind(as.matrix(airquality[,-5]),Intercept=1,index=1:nrow(airquality),
             # a numeric vector - positive values 
             weights = round(rgamma(nrow(airquality),3,3),1),
             groups = airquality[,5])
        
data_DT = data.table(data)
data_DT$groups = factor(data_DT$groups)
#additional vars
# Ozone has log-normal distiribution so we should regress a log variable and next take exp of it
data_DT$Ozone_log = log(data_DT$Ozone) 
data_DT$x_character = as.character(cut(data_DT$Solar.R,seq(0,350,70)))
data_DT$y_chac = as.character(cut(data_DT$Ozone,c(0,20,40,60,80,100,120,140,160)))

data_DT[,.(VIF(.SD,posit_y='Ozone',posit_x=c('Solar.R','Wind','Temp','x_character','Day','index','weights','groups')))]

data_DT[,Solar_R_imp:=fill_NA_N(.SD,
                           model="lm_bayes",
                           posit_y='Solar.R',
                           posit_x=c('Wind','Temp','Intercept'),
                           w=.SD[['weights']],
                           times=100),by=.(groups)] %>%
  .[,x_character_imp:=fill_NA(.SD,
                           model="lda",
                           posit_y='x_character',
                           posit_x=c('Wind','Temp','groups'))] %>%
  .[,Ozone_imp1:=exp(fill_NA_N(x=.SD, 
                           model="lm_noise",
                           posit_y='Ozone_log',
                           posit_x=c('Intercept','x_character_imp','Wind','Temp','Solar_R_imp'),
                           w=.SD[['weights']],
                           times=10))] %>% 
  .[,Ozone_imp2:=exp(fill_NA_N(x=.SD, 
                           model="lm_bayes",
                           posit_y='Ozone_log',
                           posit_x=c('Intercept','x_character_imp','Wind','Temp'),
                           w=.SD[['weights']],
                           times=10))]  %>% 
  .[,Ozone_imp3:=exp(fill_NA(x=.SD, 
                           model="lm_pred",
                           posit_y='Ozone_log',
                           posit_x=c('Intercept','x_character_imp','Wind','Temp'),
                           w=.SD[['weights']])),.(groups)]%>%
  # average of a few methods
  .[,Ozone_imp_mix := apply(.SD,1,mean),.SDcols=Ozone_imp1:Ozone_imp3] %>% 
  #Protecting against collinearity or low number of observations - across small groups
  .[,y_chac_imp:=tryCatch(fill_NA(x=.SD, 
                                   model="lda",
                                   posit_y='y_chac',
                                   posit_x=c('Intercept','Month','Day','Temp','x_character_imp','Solar_R_imp'),
                                   w=.SD[['weights']]),
                          error=function(e) .SD[['y_chac']]),.(groups)] 
                           
                           

data_DT[which(is.na(data_DT[,1]))[1:5],]


data_DT$Ozone_NA = ifelse(is.na(data_DT$Ozone),'missing','complete')

data_DT[,c('Ozone','Ozone_imp_mix','Ozone_NA')] %>% 
melt(id=c('Ozone_NA'),measure=c('Ozone','Ozone_imp_mix')) %>% 
.[!(!(Ozone_NA=='missing') & (variable=='Ozone_imp_mix')),] %>% 
ggplot2::ggplot(.,ggplot2::aes(x=value,fill=variable)) + 
ggplot2::geom_histogram(alpha=0.3) + ggplot2::facet_wrap(Ozone_NA ~.,scales='free')

```

Be carful when using a data.table grouping option because of lack of protection against collinearity or low number of observations. 
There could be used a tryCatch(fill_NA(...),error=function(e) return(...))
