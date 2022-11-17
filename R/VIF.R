#' \code{VIF} function for assessing VIF.
#'
#' @description VIF measure how much the variance of the estimated regression coefficients are inflated.
#' It helps to identify when the predictor variables are linearly related.
#' You have to decide which variable should be delete.
#' Usually values higher than 10 (around), mean a collinearity problem.
#'
#' @param x a numeric matrix or data.frame/data.table (factor/character/numeric) - variables
#' @param posit_y an integer/character - a position/name of dependent variable.
#' This variable is taken into account only for getting complete cases.
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
#' airquality2 <- airquality
#' airquality2$Temp2 <- airquality2$Temp**2
#' airquality2$Month <- factor(airquality2$Month)
#' data_DT <- data.table(airquality2)
#' data_DT[, .(vifs = VIF(
#'   x = .SD,
#'   posit_y = "Ozone",
#'   posit_x = c("Solar.R", "Wind", "Temp", "Month", "Day", "Temp2"),
#'   correct = FALSE
#' ))][["vifs.V1"]]
#'
#' data_DT[, .(vifs = VIF(
#'   x = .SD,
#'   posit_y = 1,
#'   posit_x = c(2, 3, 4, 5, 6, 7),
#'   correct = TRUE
#' ))][["vifs.V1"]]
#' }
#'
#' @name VIF
#'
#' @export

VIF <- function(x, posit_y, posit_x, correct = FALSE) {
  if (inherits(x, "data.frame") || inherits(x, "matrix") || inherits(x, "data.table")) {
    stopifnot(is.logical(correct))

    if (posit_y %in% posit_x) {
      stop("the same variable is dependent and indepentent")
    }

    cols <- colnames(x)
    # posit as character vector
    if (is.character(posit_x)) {
      posit_x <- match(posit_x, cols)
      posit_x <- posit_x[!is.na(posit_x)]
      if (length(posit_x) == 0) stop("posit_x is empty")
    } else {
      stopifnot(posit_x %in% seq_len(dim(x)[2]))
    }

    if (length(posit_x) < 2) stop("at least two independent variables should be provided")

    # contains_intercept <- any(unlist(apply(x[, posit_x], 2, function(i) all(duplicated(i)[-1L]))))
    # if (contains_intercept) stop("Do not include an intercept")

    if (is.character(posit_y)) {
      posit_y <- match(posit_y, cols)
      posit_y <- posit_y[!is.na(posit_y)]
      if (length(posit_y) == 0) stop("posit_y is empty")
    } else {
      stopifnot(posit_y %in% seq_len(dim(x)[2]))
    }

    UseMethod("VIF", x)
  } else {
    stop("wrong data type - it should be data.frame, matrix or data.table")
  }
}


#' @describeIn  VIF

VIF.data.frame <- function(x, posit_y, posit_x, correct = FALSE) {
  x_small <- x[, c(posit_y, posit_x)]
  x_small[[1]] <- as.numeric(x_small[[1]])
  xx <- model.matrix.lm(~.-1, x_small, na.action = "na.pass")
  aa <- attributes(xx)$assign
  ll <- 2:ncol(xx)
  VIF_(xx, 1, ll, aa[ll], correct)
}

#' @describeIn  VIF

VIF.data.table <- function(x, posit_y, posit_x, correct = FALSE) {
  x_small <- x[, c(posit_y, posit_x), with = FALSE]
  x_small[[1]] <- as.numeric(x_small[[1]])
  xx <- model.matrix.lm(~., x_small, na.action = "na.pass")
  aa <- attributes(xx)$assign
  ll <- 3:ncol(xx)
  VIF_(xx, 2, ll, aa[ll], correct)
}

#' @describeIn  VIF

VIF.matrix <- function(x, posit_y, posit_x, correct = FALSE) {
  x_small <- x[, c(posit_y, posit_x)]
  ncol_x <- ncol(x_small)
  VIF_(x_small, 1, 2:ncol_x, 2:ncol_x, correct)
}
