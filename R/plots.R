#' Comparing imputations and original data distributions
#' @description ggplot2 visualization to support which imputation method choose
#' @param df data.frame with origin variable and the new one with imputations
#' @param origin character value - the name of origin variable with values before any imputations
#' @param target character vector - names of variables with applied imputations
#' @return ggplot2 object
#' @examples
#' data(air_miss)
#' air_miss$Ozone_imp <- fill_NA(
#'     x = air_miss,
#'     model = "lm_bayes",
#'     posit_y = 1,
#'     posit_x = c(4, 6),
#'     logreg = TRUE
#' )
#' air_miss$Ozone_imp2 <- fill_NA_N(
#'     x = air_miss,
#'     model = "pmm",
#'     posit_y = 1,
#'     posit_x = c(4, 6),
#'     logreg = TRUE
#' )
#'
#' compare_imp(air_miss, origin = "Ozone", "Ozone_imp")
#' compare_imp(air_miss, origin = "Ozone", c("Ozone_imp", 'Ozone_imp2'))
#'
compare_imp = function(df, origin, target) {
  assert_that(inherits(df, "data.frame"))
  assert_that(inherits(origin, "character"))
  assert_that(inherits(target, "character"))

  data = as.data.frame(df)
  data$origin_NA <- ifelse(is.na(data[[origin]]) , "missing", "complete")
  data_long <- tidyr::pivot_longer(data[,c(origin, "origin_NA", target)], !all_of("origin_NA"))
  data_final <- data_long[(((data_long$origin_NA == "missing") & (data_long$name %in% target)) | ((data_long$origin_NA == "complete") & (data_long$name == origin))), ]
  ggplot2::ggplot(data_final, ggplot2::aes_string(x = "value", fill = "name", group = "name")) +
  ggplot2::geom_density(alpha = 0.4)
}

#' upset plot for NA values
#' @description wrapper around UpSetR::upset for vizualization of NA values
#' @description Visualization of set intersections using novel UpSet matrix design.
#' @param ... all arguments accepted by UpSetR::upset where the first one is expected to be a data.
#' @details Visualization of set data in the layout described by Lex and Gehlenborg in \url{https://www.nature.com/articles/nmeth.3033}.
#' UpSet also allows for visualization of queries on intersections and elements, along with custom queries queries implemented using
#' Hadley Wickham's apply function. To further analyze the data contained in the intersections, the user may select additional attribute plots
#' to be displayed alongside the UpSet plot. The user also has the the ability to pass their own plots into the function to further analyze
#' data belonging to queries of interest. Most aspects of the UpSet plot are customizable, allowing the user to select the plot that best suits their style.
#' Depending on how the features are selected, UpSet can display between 25-65 sets and between 40-100 intersections.
#' @note Data set must be formatted as described on the original UpSet github page: \url{https://github.com/VCG/upset/wiki}.
#' @references Lex et al. (2014). UpSet: Visualization of Intersecting Sets
#' IEEE Transactions on Visualization and Computer Graphics (Proceedings of InfoVis 2014), vol 20, pp. 1983-1992, (2014).
#' @references Lex and Gehlenborg (2014). Points of view: Sets and intersections. Nature Methods 11, 779 (2014). \url{https://www.nature.com/articles/nmeth.3033}
#' @examples
#' upset_NA(airquality)
#' upset_NA(air_miss, 6)
#'
upset_NA <- function(...) {
  args <- list(...)
  assert_that(inherits(args[[1]], "data.frame"))
  args[[1]] <- data.frame(Map(function(x) as.integer(is.na(x)), args[[1]]))
  do.call(UpSetR::upset, args)
}
