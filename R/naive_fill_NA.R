#' \code{naive_fill_NA} function for the simple and automatic imputation
#'
#' @description
#' Automatically fill the missing data with a simple imputation method,
#' impute with sampling the non missing values.
#' It is recommended to use this function for each categorical variable separately.
#'
#' @param x a numeric matrix or data.frame/data.table (factor/character/numeric/logical variables)
#' @return object with a similar structure to the input but without missing values.
#'
#' @note this is a very simple and fast solution but not recommended,
#' for more complex solutions please check the vignette.
#'
#' @seealso \code{\link{fill_NA}} \code{\link{fill_NA_N}}  \code{\link{VIF}}
#'
#' @examples
#' \dontrun{
#' library(miceFast)
#' data(air_miss)
#' naive_fill_NA(air_miss)
#' # Could be useful to run it separately for each group level
#' do.call(rbind, Map(naive_fill_NA, split(air_miss, air_miss$groups)))
#' }
#'
#' @name naive_fill_NA
#'
#' @export
naive_fill_NA <- function(x) {
  if (inherits(x, "data.frame") || inherits(x, "matrix") || inherits(x, "data.table") ) {
    UseMethod("naive_fill_NA", x)
  } else {
    stop("wrong data type - it should be data.frame, matrix or data.table")
  }
}

#' @describeIn naive_fill_NA S3 method for data.frame
naive_fill_NA.data.frame <- function(x) {
  any_na <- vapply(
    X = x,
    FUN = function(c) anyNA(c),
    FUN.VALUE = logical(1)
  )
  if (any(any_na)) {
    cols <- which(any_na)
    for (icol in cols) {
      vec <- x[[icol]]
      is_na <- is.na(vec)
      sum_na <- sum(is_na)
      sum_nonna <- length(vec) - sum_na
      if (any(is_na) && (sum_nonna >= 1)) {
        x[is_na, icol] <- get_replacements(vec, !is_na, sum_nonna, sum_na)
      }
    }
  }
  return(x)
}

#' @describeIn naive_fill_NA S3 method for data.table
naive_fill_NA.data.table <- function(x) {
  any_na <- vapply(
    X = x,
    FUN = function(c) anyNA(c),
    FUN.VALUE = logical(1)
  )
  if (any(any_na)) {
    cols <- which(any_na)
    for (icol in cols) {
      vec <- x[[icol]]
      is_na <- is.na(vec)
      sum_na <- sum(is_na)
      sum_nonna <- length(vec) - sum_na
      if (any(is_na) && (sum_nonna >= 1)) {
        data.table::set(
          x = x,
          i = as.integer(which(is_na)),
          j = icol,
          value = get_replacements(vec, !is_na, sum_nonna, sum_na)
        )
      }
    }
  }
  return(x)
}

#' @describeIn naive_fill_NA S3 method for matrix
naive_fill_NA.matrix <- function(x) {
  any_na <- vapply(
    X = seq_len(ncol(x)),
    FUN = function(c) anyNA(x[, c]),
    FUN.VALUE = logical(1)
  )
  if (any(any_na)) {
    cols <- which(any_na)
    for (icol in cols) {
      vec <- x[, icol]
      is_na <- is.na(vec)
      sum_na <- sum(is_na)
      sum_nonna <- length(vec) - sum_na
      if (any(is_na) && (sum_nonna >= 1)) {
        x[is_na, icol] <- get_replacements(vec, !is_na, sum_nonna, sum_na)
      }
    }
  }
  return(x)
}


#' @keywords internal
get_replacements <- function(
    vec,
    is_not_na = !is.na(vec),
    sum_nonna = sum(is_not_na),
    sum_na = length(vec) - sum_nonna
  ) {
  if (sum_nonna == 1) {
    rep(x = vec[is_not_na], times = sum_na)
  } else {
    base::sample(
      x = vec[is_not_na],
      size = sum_na,
      replace = TRUE
    )
  }
}
