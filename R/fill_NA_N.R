
#' \code{fill_NA_N} function for the multiple imputations purpose.
#'
#' @description Multiple imputations to fill the missing data.
#' Non missing independent variables are used to approximate a missing observations for a dependent variable.
#' Quantitative models were built under Rcpp packages and the C++ library Armadillo.
#'
#' @param x a numeric matrix or data.frame/data.table (factor/character/numeric) - variables
#' @param model a character - posibble options ("lm_bayes","lm_noise")
#' @param posit_y an integer/character - a position/name of dependent variable
#' @param posit_x an integer/character vector - positions/names of independent variables
#' @param w  a numeric vector - a weighting variable - only positive values
#' @param times an integer - a number of multiple imputations - default 10
#'
#' @return load imputations in a numeric/factor (depend on an input type) vector format
#'
#' @note
#' The lda model is assessed only if there are more than 15 complete observations
#' and for the lms models if number of variables is smaller than number of observations.
#'
#' @seealso \code{\link{fill_NA}} \code{\link{VIF}}
#'
#' @examples
#' \dontrun{
#' library(miceFast)
#' library(data.table)
#' library(magrittr)
#' library(mice)
#' library(dplyr)
#'
#' ### Intro:data.table
#' #### Working with names
#' data = cbind(as.matrix(mice::nhanes),intercept=1,index=1:nrow(mice::nhanes))
#' data = do.call(rbind,replicate(10,data,simplify = F))
#' data_df = as.data.frame(data)
#' data_DT = data.table(data)
#'
#' data_DT[,bmi_imp:=fill_NA(x=.SD,
#'                          model="lm_pred",
#'                          posit_y='bmi',
#'                          posit_x='intercept')] %>%
#'  .[,hyp_imp:=fill_NA(x=.SD,
#'                      model="lda",
#'                      posit_y='hyp',
#'                      posit_x=c('age','bmi_imp')),] %>%
#'  .[,chl_imp:=fill_NA_N(x=.SD,
#'                        model="lm_noise",
#'                        posit_y='chl',
#'                        posit_x=c('age','bmi_imp','hyp_imp'),
#'                        times=10),]
#'
#' head(data_DT,2)
#'
#' #### Working with positions
#'
#' data_DT[,bmi_imp:=fill_NA(x=.SD,
#'                          model="lm_pred",
#'                          posit_y=2,
#'                          posit_x=5)] %>%
#'  # there is a new variable at position 7 - bmi_imp
#'  .[,hyp_imp:=fill_NA(x=.SD,
#'                      model="lda",
#'                      posit_y=3,
#'                      posit_x=c(1,7)),] %>%
#'  .[,chl_imp:=fill_NA_N(x=.SD,
#'                        model="lm_noise",
#'                        posit_y=4,
#'                        posit_x=c(1,7,8),
#'                        times=10),]
#'
#' head(data_DT,2)
#'
#' ### Intro:dplyr
#'
#' #### Working with names
#'
#' data_df = data_df %>% mutate(bmi_imp=fill_NA(x=.,
#'                                             model="lm_pred",
#'                                             posit_y='bmi',
#'                                             posit_x='intercept')) %>%
#'  mutate(hyp_imp=fill_NA(x=.,
#'                         model="lda",
#'                         posit_y='hyp',
#'                         posit_x=c('age','bmi_imp'))) %>%
#'  mutate(chl_imp=fill_NA_N(x=.,
#'                           model="lm_noise",
#'                           posit_y='chl',
#'                           posit_x=c('age','bmi_imp','hyp_imp'),
#'                           times=10))
#'
#' head(data_df,2)
#'
#'
#'
#' #### Working with positions
#'
#' data_df = data_df %>% mutate(bmi_imp=fill_NA(x=.,
#'                                             model="lm_pred",
#'                                             posit_y=2,
#'                                             posit_x=5)) %>%
#'  #there is a new variable at position 7 - bmi_imp
#'  mutate(hyp_imp=fill_NA(x=.,
#'                         model="lda",
#'                         posit_y=3,
#'                         posit_x=c(1,7))) %>%
#'  mutate(chl_imp=fill_NA_N(x=.,
#'                           model="lm_noise",
#'                           posit_y=4,
#'                           posit_x=c(1,7,8),
#'                           times=10))
#'
#' head(data_df,2)
#'
#' ### Model with additional parameters: - data with the grouping/weighting variable
#' ### data.table recommended
#'
#' data = cbind(as.matrix(airquality[,-5]),Intercept=1,index=1:nrow(airquality),
#'             # a numeric vector - positive values
#'             weights = round(rgamma(nrow(airquality),3,3),1),
#'             groups = airquality[,5])
#'
#' data = do.call(rbind,replicate(10,data,simplify = F))
#'
#' data_DT = data.table(data)
#'
#'
#' #### Working with names
#'
#' data_DT[,Ozone_imp:=fill_NA(x=.SD,
#'                            model="lm_pred",
#'                            posit_y='Ozone',
#'                            posit_x='Intercept',w=.SD[['weights']]),by=.(groups)] %>%
#'  .[,Solar_R_imp:=fill_NA_N(.SD,
#'                            model="lm_bayes",
#'                            posit_y='Solar.R',
#'                            posit_x=c('Wind','Temp','Day','Intercept','Ozone_imp'),
#'                            w=.SD[['weights']],
#'                            times=10),by=.(groups)]
#'
#' data_DT[which(is.na(data_DT[,1]))[1],]
#'
#'
#' #### Working with positions
#'
#' # simple mean imputation - intercept at position 6
#' data_DT[,Ozone_imp:=fill_NA(x=.SD,
#'                            model="lm_pred",
#'                            posit_y=1,
#'                            posit_x=c(6),w=.SD[['weights']]),by=.(groups)] %>%
#'  # avg of 10 multiple imputations - last posit_x equal to 9 not 10
#'  # because the groups variable is not included in .SD
#'  .[,Solar_R_imp:=fill_NA_N(.SD,
#'                            model="lm_bayes",
#'                            posit_y=2,
#'                            posit_x=c(3,4,5,6,9),
#'                            w=.SD[['weights']],times=10),by=.(groups)]
#'
#' data_DT[which(is.na(data_DT[,1]))[1],]
#' }
#'
#' @name fill_NA_N
#'
#' @export

fill_NA_N <- function(x, model, posit_y, posit_x, w = 0L,times=1) {

  UseMethod('fill_NA_N')

}

#' @rdname fill_NA_N

fill_NA_N.data.frame <- function(x, model, posit_y, posit_x, w = 0L,times=1) {

  cols = colnames(x)
  # posit as character vector
  if(is.character(posit_x)) {
    posit_x = pmatch(posit_x,cols)
    posit_x = posit_x[!is.na(posit_x)]
    if(length(posit_x)==0) stop('posit_x is empty')
  }

  if(is.character(posit_y)){
    posit_y = pmatch(posit_y,cols)
    if(length(posit_y)==0) stop('posit_y is empty')
  }

  #Intercept


  if(inherits(x,'data.frame')){

    x_small = subset(x,select=posit_x)

    types = lapply(x_small,class)

    x_ncols = length(posit_x)

    contains_intercept = any(unlist(lapply(x_small,function(i) all(duplicated(i)[-1L]))))

    p_x_factor_character = which(unlist(lapply(types,function(i) !all(is.na(match(c('factor','character'),i))))))

    len_p_x_factor_character = length(p_x_factor_character)

    yy = x[[posit_y]]

    xx = vector('list',2)

    if(len_p_x_factor_character>0){

      x_fc = subset(x,select=posit_x[p_x_factor_character])
      x_fc = model.matrix.lm(~.,x_fc,na.action="na.pass")

      if(contains_intercept) x_fc = x_fc[,-1]

      xx[[1]] = x_fc

    }

    if(x_ncols>len_p_x_factor_character){
      x_ni = as.matrix(subset(x,select=setdiff(posit_x,posit_x[p_x_factor_character])))
      xx[[2]] = x_ni
    }

    xx = do.call(cbind,xx[!is.null(xx)])

    if(is.factor(yy)){

      l=levels(yy)
      yy = as.numeric(yy)

      f = round(fill_NA_N_(cbind(yy,xx), model, 1, 2:(ncol(xx)+1), w,times))

      f[f<=0] = 1
      f[f>length(l)] = length(l)
      ff = factor(l[f])

    } else if(is.character(yy)){

      yy = factor(yy)
      l=levels(yy)
      yy = as.numeric(yy)
      yy = yy

      f = round(fill_NA_N_(cbind(yy,xx), model, 1, 2:(ncol(xx)+1), w,times))
      f[f<=0] = 1
      f[f>length(l)] = length(l)
      ff = factor(l[f])

    } else {
      ff = fill_NA_N_(cbind(yy,xx), model, 1, 2:(ncol(xx)+1), w,times)
    }
  } else {stop("wrong data type - it should be data.frame")}

  attr(ff,'dim') = attributes(ff)$dim[1]

  return(ff)


}

#' @rdname fill_NA_N

fill_NA_N.matrix <- function(x, model, posit_y, posit_x, w = 0L,times=1) {

  cols = colnames(x)
  # posit as character vector
  if(is.character(posit_x)) {
    posit_x = pmatch(posit_x,cols)
    posit_x = posit_x[!is.na(posit_x)]
    if(length(posit_x)==0) stop('posit_x is empty')
  }

  if(is.character(posit_y)){
    posit_y = pmatch(posit_y,cols)
    if(length(posit_y)==0) stop('posit_y is empty')
  }

  if(is.matrix(x)){

  ff = fill_NA_N_(x, model, posit_y, posit_x, w,times)

  } else {stop("wrong data type - it should be data.frame or matrix")}

  attr(ff,'dim') = attributes(ff)$dim[1]

  return(ff)

}
