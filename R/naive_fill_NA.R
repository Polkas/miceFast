#' \code{naive_fill_NA} function for the imputations purpose.
#'
#' @description
#' Regular imputations to fill the missing data.
#' Non missing independent variables are used to approximate a missing observations for a dependent variable.
#' For numeric columns with any missing data a simple bayesian mean will be used.
#' Next all numeric variables will be utilized to impute character/factoer variables by Linear Discriminant Analysis.
#'
#' @param x a numeric matrix or data.frame/data.table (factor/character/numeric/logical)  - variables
#'
#' @return load imputations in a similar to the input type
#'
#' @note this is a very simple and fast solution but not recommended, for more complex solutions check the vignette
#'
#' @seealso \code{\link{fill_NA}} \code{\link{fill_NA_N}}  \code{\link{VIF}}
#'
#' @examples
#' \dontrun{
#' library(naive_fill_NAast)
#' data(air_miss)
#' naive_fill_NA(air_miss)
#' }
#'
#' @name naive_fill_NA
#'
#' @export
naive_fill_NA <- function(x) {
  if (inherits(x, "data.frame") || inherits(x, "matrix")) {
    UseMethod("naive_fill_NA")
  } else {
    stop("wrong data type - it should be data.frame, matrix or data.table")
  }
}

# check which have less than 15 levels
# matrix option

#' @describeIn naive_fill_NA S3 method for data.frame
naive_fill_NA.data.frame <- function(x) {
  types <- vapply(x, class, FUN.VALUE = character(1))

  types_n <- types %in% c("numeric", "integer")

  types_fcl <- types %in% c("factor", "character", "logical")

  na_col_p <- vapply(x, function(i) mean(is.na(i)), FUN.VALUE = numeric(1))

  col_names <- colnames(x)

  mm <- match(c("weights", "wei", "weis", "w"), tolower(colnames(x)))

  ww <- if (all(is.na(mm))) vector() else x[[na.omit(mm)[1]]]

  ww <- if (any(is.na(ww))) vector() else ww

  nn <- col_names[types_n & (na_col_p > 0)]

  for (posit_y in nn) {

    # full_n <- (types_n) & (na_col_p_new == 0)

    yy <- x[[posit_y]]

    # p_posit_y <- match(posit_y, col_names)

    yy_class <- class(yy)

    all_pos_y <- FALSE

    if (is.numeric(yy)) {
      all_pos_y <- !any(yy < 0, na.rm = TRUE)
    }

    mean_over_median <- mean(yy, na.rm = T) > median(yy, na.rm = T)

    log_transform <- all_pos_y && mean_over_median

    pp <- as.matrix(rep(1, length(yy)))

    # num_add <- setdiff(which(full_n), p_posit_y)

    # xx <- as.matrix(if (is_DT) x[, num_add, with = FALSE] else x[, num_add])

    if (log_transform) {
      yy <- log(yy + 1e-8)
      data_temp <- cbind(yy, 1)
      ff <- exp(fill_NA_(data_temp, "lm_bayes", 1, 2, ww))
    } else {
      data_temp <- cbind(yy, 1)
      ff <- fill_NA_(data_temp, "lm_bayes", 1, 2, ww)
    }

    x[[posit_y]] <- as.vector(ff)
  }


  if (any(types_n)) {
    for (posit_y in col_names[types_fcl & (na_col_p > 0)]) {
      yy <- x[[posit_y]]

      n_unique <- length(unique(yy))

      if ((n_unique < 2) || (n_unique > 15)) next

      yy_class <- class(yy)

      is_yy_factor <- is.factor(yy)

      yy <- if (is_yy_factor) yy else factor(yy)
      l <- levels(yy)
      yy <- as.numeric(yy)
      pp <- prcomp(Filter(is.numeric, x))$x[, 1]
      f <- round(fill_NA_(cbind(yy, pp), "lda", 1, 2, vector()))
      f[f <= 0] <- 1
      f[f > length(l)] <- length(l)
      ff <- if (is_yy_factor) factor(l[f]) else l[f]

      class(ff) <- yy_class

      attr(ff, "dim") <- attributes(ff)$dim[1]

      x[[posit_y]] <- ff
    }
  }

  return(x)
}

#' @describeIn naive_fill_NA S3 method for matrix

naive_fill_NA.matrix <- function(x) {
  stopifnot(is.numeric(x))

  na_col_p <- vapply(x, function(i) mean(is.na(i)), FUN.VALUE = numeric(1))

  col_names <- colnames(x)

  mm <- match(c("weights", "wei", "weis", "w"), tolower(colnames(x)))

  ww <- if (all(is.na(mm))) vector() else x[[na.omit(mm)[1]]]

  ww <- if (any(is.na(ww))) vector() else ww

  nn <- col_names[(na_col_p > 0)]

  for (posit_y in nn) {
    yy <- x[[posit_y]]

    all_pos_y <- !any(yy < 0, na.rm = TRUE)

    mean_over_median <- mean(yy, na.rm = T) > median(yy, na.rm = T)

    log_transform <- all_pos_y && mean_over_median

    if (log_transform) {
      yy <- log(yy + 1e-8)
      data_temp <- cbind(yy, 1)
      ff <- exp(fill_NA_(data_temp, "lm_bayes", 1, 2, ww))
    } else {
      data_temp <- cbind(yy, 1)
      ff <- fill_NA_(data_temp, "lm_bayes", 1, 2, ww)
    }

    x[[posit_y]] <- as.vector(ff)
  }

  return(x)
}
