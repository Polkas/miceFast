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

## Introduction for data.table/dplyr users - using additional functions from miceFast:

```r
#install.packages('pacman')
pacman::p_load(miceFast,data.table,magrittr,car,dplyr,ggplot2)
```

Usage of `fill_NA` and `fill_NA_N` functions from miceFast - this functions should be resistant to glitches from an user activity perspective and a data structures.

### Data 

```r
# airquality dataset with additional variables
data(air_miss)
```

### Intro: data.table

```r
# VIF - values bigger than 10 (around) suggest that there might be a collinearity problem.
# VIF is high for Solar.R and x_character which is obvious - x_character is a factor version of numeric Solar.R
air_miss[,.(VIF(.SD,posit_y='Ozone',
                   posit_x=c('Solar.R',
                             'Wind',
                             'Temp',
                             'x_character',
                             'Day',
                             'weights',
                             'groups')))] 
                             
# IMPUTATIONS
# Imputations with a grouping option (models are separately assessed for each group) 
# taking into account provided weights
air_miss[,Solar_R_imp := fill_NA_N(x=.SD,
                                   model="lm_bayes",
                                   posit_y='Solar.R',
                                   posit_x=c('Wind','Temp','Intercept'),
                                   w=.SD[['weights']],
                                   times=100),by=.(groups)] %>%
# Imputations - discrete variable
  .[,x_character_imp := fill_NA(x=.SD,
                                model="lda",
                                posit_y='x_character',
                                posit_x=c('Wind','Temp','groups'))] %>%
# logreg was used because almost log-normal distribution of Ozone
# imputations around mean
.[,Ozone_imp1 := fill_NA(x=.SD, 
                         model="lm_bayes",
                         posit_y='Ozone',
                         posit_x=c('Intercept'),
                         logreg=TRUE)] %>% 
# imputations using positions - Intercept, Temp
.[,Ozone_imp2 := fill_NA(x=.SD, 
                         model="lm_bayes",
                         posit_y=1,
                         posit_x=c(4,6)),
                         logreg=TRUE] %>% 
# model with a factor independent variable 
# multiple imputations (average of x30 imputations) 
# with a factor independent variable, weights and logreg options
  .[,Ozone_imp3 := fill_NA_N(x=.SD, 
                             model="lm_noise",
                             posit_y='Ozone',
                             posit_x=c('Intercept','x_character_imp','Wind','Temp'),
                             w=.SD[['weights']],
                             logreg=TRUE,
                             times=30)] %>% 
  .[,Ozone_imp4 := fill_NA_N(x=.SD, 
                             model="lm_bayes",
                             posit_y='Ozone',
                             posit_x=c('Intercept','x_character_imp','Wind','Temp'),
                             w=.SD[['weights']],
                             logreg=TRUE,
                             times=30)] %>% 
  .[,Ozone_imp5 := fill_NA(x=.SD, 
                           model="lm_pred",
                           posit_y='Ozone',
                           posit_x=c('Intercept','x_character_imp','Wind','Temp'),
                           w=.SD[['weights']],
                           logreg=TRUE),.(groups)] %>%
                           
# Average of a few methods
# First imputations was skipped, lack of logreg option brings a bad prediction model
  .[,Ozone_imp_mix := apply(.SD,1,mean),.SDcols=Ozone_imp3:Ozone_imp5] %>% 
  
# Protecting against collinearity or low number of observations - across small groups
# Be carful when using a data.table grouping option 
# because of lack of protection against collinearity or low number of observations. 
# There could be used a tryCatch(fill_NA(...),error=function(e) return(...))

  .[,Ozone_chac_imp := tryCatch(fill_NA(x=.SD, 
                                 model="lda",
                                 posit_y='Ozone_chac',
                                 posit_x=c('Intercept','Month','Day','Temp','x_character_imp'),
                                 w=.SD[['weights']]),
                                error=function(e) .SD[['Ozone_chac']]),.(groups)] 
                           
# Sample of results                           
air_miss[which(is.na(air_miss[,1]))[1:5],]
```

### Intro: dplyr

```r
# VIF - values bigger than 10 (around) suggest that there might be a collinearity problem.
# VIF is high for Solar.R and x_character which is obvious - x_character is a factor version of numeric Solar.R
air_miss %>% do(vifs=VIF(.,posit_y='Ozone',
                           posit_x=c('Solar.R',
                                     'Wind',
                                     'Temp',
                                     'x_character',
                                     'Day',
                                     'weights',
                                     'groups'))) %>% unlist()
                             
# IMPUTATIONS
air_miss = air_miss %>% 
# Imputations with a grouping option (models are separately assessed for each group) 
# taking into account provided weights
group_by(groups) %>% 
do(mutate(.,Solar_R_imp = fill_NA(x=.,
                                  model="lm_pred",
                                  posit_y='Solar.R',
                                  posit_x=c('Wind','Temp','Intercept'),
                                  w=.[['weights']]))) %>%
ungroup() %>%
# Imputations - discrete variable
mutate(x_character_imp = fill_NA(x=.,
                                 model="lda",
                                 posit_y='x_character',
                                 posit_x=c('Wind','Temp'))) %>%
# logreg was used because almost log-normal distribution of Ozone
# imputations around mean
mutate(Ozone_imp1 = fill_NA(x=., 
                            model="lm_bayes",
                            posit_y='Ozone',
                            posit_x=c('Intercept'),
                            logreg=TRUE)) %>% 
# imputations using positions - Intercept, Temp
mutate(Ozone_imp2 = fill_NA(x=., 
                            model="lm_bayes",
                            posit_y=1,
                            posit_x=c(4,6),
                            logreg=TRUE)) %>% 
# multiple imputations (average of x30 imputations) 
# with a factor independent variable, weights and logreg options
mutate(Ozone_imp3 = fill_NA_N(x=., 
                              model="lm_noise",
                              posit_y='Ozone',
                              posit_x=c('Intercept','x_character_imp','Wind','Temp'),
                              w=.[['weights']],
                              logreg=TRUE,
                              times=30)) %>% 
mutate(Ozone_imp4 = fill_NA_N(x=., 
                              model="lm_bayes",
                              posit_y='Ozone',
                              posit_x=c('Intercept','x_character_imp','Wind','Temp'),
                              w=.[['weights']],
                              logreg=TRUE,
                              times=30)) %>% 
group_by(groups) %>%
do(mutate(.,Ozone_imp5 = fill_NA(x=., 
                                 model="lm_pred",
                                 posit_y='Ozone',
                                 posit_x=c('Intercept','x_character_imp','Wind','Temp'),
                                 w=.[['weights']],
                                 logreg=TRUE))) %>%
ungroup() %>%
# Average of a few methods
mutate(Ozone_imp_mix = rowMeans(select(.,starts_with("Ozone_imp")))) %>%

# Protecting against collinearity or low number of observations - across small groups
# Be carful when using a data.table grouping option 
# because of lack of protection against collinearity or low number of observations. 
# There could be used a tryCatch(fill_NA(...),error=function(e) return(...))
group_by(groups) %>%
do(mutate(.,Ozone_chac_imp = tryCatch(fill_NA(x=., 
                                       model="lda",
                                       posit_y='Ozone_chac',
                                       posit_x=c('Intercept','Month','Day','Temp','x_character_imp'),
                                       w=.[['weights']]),
                                      error=function(e) .[['Ozone_chac']]))) %>%
ungroup()

# Sample of results                           
air_miss[which(is.na(air_miss[,1]))[1:5],]
```

### results - visualization

```r
# Distribution of imputations vs Distribution of initial data
air_miss$Ozone_NA = ifelse(is.na(air_miss$Ozone),'imputations','complete')
air_miss = as.data.table(air_miss)

air_miss[,c('Ozone','Ozone_imp_mix','Ozone_NA')] %>% 
melt(id=c('Ozone_NA'),measure=c('Ozone','Ozone_imp_mix')) %>% 
.[(((Ozone_NA=='imputations') & (variable=='Ozone_imp_mix'))|((Ozone_NA=='complete') & (variable=='Ozone'))),] %>% 
ggplot2::ggplot(.,ggplot2::aes(x=value,fill=variable)) + 
ggplot2::geom_density() + ggplot2::facet_wrap(Ozone_NA ~.,scales='free')
```
