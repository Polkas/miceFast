#' \code{VIF} function for assessing VIF.
#'
#' @description VIF measure how much the variance of the estimated regression coefficients are inflated.
#' It helps to identify when the predictor variables are linearly related.
#' You have to decide which variable should be delete. Values higher than 10 signal a potential collinearity problem.
#'
#' @param x a numeric matrix or data.frame/data.table (factor/character/numeric) - variables
#' @param posit_y an integer/character - a position/name of dependent variable
#' @param posit_x an integer/character vector - positions/names of independent variables
#' @param correct a boolean - basic or corrected - Default: FALSE
#'
#' @note vif_corrected = vif_basic^{(1/(2*df))}
#'
#' @return load a numeric vector with VIF for all variables provided by posit_x
#'
#' @seealso \code{\link{fill_NA}} \code{\link{fill_NA_N}}
#'
#' @examples
#' \dontrun{
#' library(miceFast)
#' library(data.table)
#'
#' airquality2 = airquality
#' airquality2$Temp2 = airquality2$Temp**2
#' #install.packages("car")
#' #car::vif(lm(Ozone ~ ., data=airquality2))
#'
#' data_DT = data.table(airquality2)
#' data_DT[,.(vifs=VIF(x=.SD,
#'                    posit_y='Ozone',
#'                    posit_x=c('Solar.R','Wind','Temp','Month','Day','Temp2'),
#'                    correct=FALSE))][['vifs']]
#'
#' data_DT[,.(vifs=VIF(x=.SD,
#'                    posit_y=1,
#'                    posit_x=c(2,3,4,5,6,7),
#'                    correct=TRUE))][['vifs']]
#'
#' }
#'
#' @name VIF
#'
#' @export

VIF = function(x, posit_y, posit_x, correct = FALSE){
UseMethod('VIF')
}


#' @describeIn  VIF

VIF.data.frame <- function(x, posit_y, posit_x, correct = FALSE ) {


if(inherits(x,'data.frame')){

  is_DT = inherits(x,'data.table')

  if(posit_y %in% posit_x){stop("the same variable is dependent and indepentent");}

  #contains_intercept = any(unlist(lapply(x_small,function(i) all(duplicated(i)[-1L]))))
  #if(contains_intercept){stop("Do not include an intercept");}

  cols = colnames(x)
  # posit as character vector
  if(is.character(posit_x)) {
    posit_x = pmatch(posit_x,cols)
    posit_x = posit_x[!is.na(posit_x)]
    if(length(posit_x)==0) stop('posit_x is empty')
  }

  if(length(posit_x)<2) stop("at least two independent variables should be provided")


  if(is.character(posit_y)){
    posit_y = pmatch(posit_y,cols)
    if(length(posit_y)==0) stop('posit_y is empty')
  }



      x_small = if(is_DT) x[,c(posit_y,posit_x),with=FALSE] else x[,c(posit_y,posit_x)]

      x_small[[1]] = as.numeric( x_small[[1]] )

      xx = model.matrix.lm(~.,x_small,na.action="na.pass")

      aa = attributes(xx)$assign

      ll = 3:ncol(xx)

      VIF_(xx, 2 , ll,aa[ll],correct)

}
}

#' @describeIn  VIF

VIF.matrix <- function(x, posit_y, posit_x, correct = FALSE ) {

if(inherits(x,'matrix')){

if(posit_y %in% posit_x){stop("the same variable is dependent and indepentent");}

#contains_intercept = any(unlist(apply(x[,posit_x],2,function(i) all(duplicated(i)[-1L]))))
#if(contains_intercept){stop("Do not include an intercept");}

cols = colnames(x)
# posit as character vector
if(is.character(posit_x)) {
  posit_x = match(posit_x,cols)
  posit_x = posit_x[!is.na(posit_x)]
  if(length(posit_x)==0) stop('posit_x is empty')
}

if(length(posit_x)<2) stop("at least two independent variables should be provided")

if(is.character(posit_y)){
  posit_y = match(posit_y,cols)
  if(length(posit_y)==0) stop('posit_y is empty')
}

x_small = x[,c(posit_y,posit_x)]

    ncol_x = ncol(x_small)

    VIF_(x_small, 1 , 2:ncol_x,2:ncol_x,correct)

} else {stop("wrong data type - it should be data.frame or matrix")}


}



