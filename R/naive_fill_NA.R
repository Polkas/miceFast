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
  if (inherits(x, "data.frame") || inherits(x, "matrix")) {
    UseMethod("naive_fill_NA", x)
  } else {
    stop("wrong data type - it should be data.frame, matrix or data.table")
  }
}

#' @describeIn naive_fill_NA S3 method for data.frame
naive_fill_NA.data.frame <- function(x) {
  naive_fill_NA.matrix(x)
}

#' @describeIn naive_fill_NA S3 method for data.table
naive_fill_NA.data.table <- function(x) {
  any_na <- vapply(x, function(c) anyNA(c), logical(1))
  if (any(any_na)) {
    cols <- which(any_na)
    for (icol in cols) {
      vec <- as.vector(x[[icol]])
      is_na <- is.na(vec)
      sum_na <- sum(is_na)
      sum_nonna <- sum(!is_na)
      if (any(is_na) && (sum_nonna > 1)) {
        replacements <- base::sample(
          vec[!is_na],
          sum_na,
          replace = TRUE
        )
        data.table::set(x, as.integer(which(is_na)), icol, replacements)
      }
    }
  }
  return(x)
}

#' @describeIn naive_fill_NA S3 method for matrix
naive_fill_NA.matrix <- function(x) {
  any_na <- vapply(seq_len(ncol(x)), function(c) anyNA(x[, c]), logical(1))
  if (any(any_na)) {
    cols <- which(any_na)
    for (icol in cols) {
      vec <- as.vector(x[, icol])
      is_na <- is.na(vec)
      sum_na <- sum(is_na)
      sum_nonna <- sum(!is_na)
      if (any(is_na) && (sum_nonna > 1)) {
        replacements <- base::sample(
          vec[!is_na],
          sum_na,
          replace = TRUE
        )
        x[is_na, icol] <- replacements
      }
    }
  }
  return(x)
}


