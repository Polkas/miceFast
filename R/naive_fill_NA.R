#' \code{naive_fill_NA} function for the imputations purpose.
#'
#' @description Regular imputations to fill the missing data.
#' Non missing independent variables are used to approximate a missing observations for a dependent variable.
#' Quantitative models were built under Rcpp packages and the C++ library Armadillo.
#'
#' @param x a numeric matrix or data.frame/data.table (factor/character/numeric/logical)  - variables
#'
#' @return load imputations in a similar to the input type
#'
#' @note
#' There is assumed that users add the intercept by their own.
#' The miceFast module provides the most efficient environment, the second recommended option is to use data.table and the numeric matrix data type.
#' The lda model is assessed only if there are more than 15 complete observations
#' and for the lms models if number of independent variables is smaller than number of observations.
#'
#' @seealso \code{\link{fill_NA}} \code{\link{fill_NA_N}}  \code{\link{VIF}}
#'
#' @examples
#' \dontrun{
#' library(miceFast)
#' data(air_miss)
#' naive_fill_NA(air_miss)
#' }
#' 
#' @name naive_fill_NA
#'
#' @export
naive_fill_NA <- function(x) {
  UseMethod("naive_fill_NA")
}

# check which have less than 15 levels
# matrix option

#' @describeIn naive_fill_NA
naive_fill_NA.data.frame <- function(x) {
  if (inherits(x, "data.frame")) {
    is_DT <- inherits(x, "data.table")

    types <- vapply(x, class, FUN.VALUE = character(1))

    types_n <- types %in% c("numeric", "integer")

    types_fcl <- types %in% c("factor", "character", "logical")

    na_col_p <- vapply(x, function(i) mean(is.na(i)), FUN.VALUE = numeric(1))

    col_names <- colnames(x)

    models <- list(character = "lda", factor = "lda", numeric = "lm_bayes", integer = "lm_bayes")

    mm <- match(c("weights", "wei", "weis", "w"), tolower(colnames(x)))

    ww <- if (all(is.na(mm))) vector() else x[[na.omit(mm)[1]]]

    nn <- col_names[types_n & (na_col_p > 0)]

    for (posit_y in nn) {
      na_col_p_new <- vapply(x, function(i) mean(is.na(i)), FUN.VALUE = numeric(1))

      full_n <- (types_n) & (na_col_p_new == 0)

      yy <- x[[posit_y]]

      p_posit_y <- match(posit_y, col_names)

      yy_class <- class(yy)

      model_y <- models[[yy_class]]

      all_pos_y <- FALSE

      if (yy_class %in% c("numeric", "integer")) {
        all_pos_y <- !any(yy < 0, na.rm = TRUE)
      }

      mean_over_median <- mean(yy, na.rm = T) > median(yy, na.rm = T)

      log_transform <- all_pos_y && mean_over_median

      pp <- as.matrix(rep(1, length(yy)))

      num_add <- setdiff(which(full_n), p_posit_y)

      xx <- as.matrix(if (is_DT) x[, num_add, with = FALSE] else x[, num_add])

      len_yy <- length(yy)

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

        yy <- if (yy_class == "factor") yy else factor(yy)
        l <- levels(yy)
        yy <- as.numeric(yy)
        pp <- prcomp(Filter(is.numeric, x))$x[, 1]
        f <- round(fill_NA_(cbind(yy, pp), "lda", 1, 2, vector()))
        f[f <= 0] <- 1
        f[f > length(l)] <- length(l)
        ff <- if (yy_class == "factor") factor(l[f]) else l[f]

        class(ff) <- yy_class

        x[[posit_y]] <- as.vector(ff)
      }
    }

    return(x)
  }
}

#' @describeIn naive_fill_NA

naive_fill_NA.matrix <- function(x) {
  if (inherits(x, "matrix")) {
    stopifnot(typeof(x) %in% c("numeric", "integer", "double"))

    na_col_p <- vapply(x, function(i) mean(is.na(i)), FUN.VALUE = numeric(1))

    col_names <- colnames(x)

    mm <- match(c("weights", "wei", "weis", "w"), tolower(colnames(x)))

    ww <- if (all(is.na(mm))) vector() else x[[na.omit(mm)[1]]]

    nn <- col_names[(na_col_p > 0)]

    for (posit_y in nn) {
      yy <- x[[posit_y]]

      p_posit_y <- match(posit_y, col_names)

      yy_class <- class(yy)

      model_y <- models[[yy_class]]

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
}
