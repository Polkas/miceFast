testthat::test_that("compare_imp", {
  data(air_miss)
  air_miss$Ozone_imp <- fill_NA(
    x = air_miss,
    model = "lm_bayes",
    posit_y = 1,
    posit_x = c(4, 6),
    logreg = TRUE
  )

  testthat::expect_silent(compare_imp(air_miss, origin = "Ozone", "Ozone_imp"))

  testthat::expect_error(compare_imp(air_miss, origin = "WRONG", "Ozone_imp"))
  testthat::expect_error(compare_imp(2, origin = "WRONG", "Ozone_imp"))
  testthat::expect_error(compare_imp(air_miss, origin = 2, "Ozone_imp"))
})

testthat::test_that("upset_NA", {
  testthat::expect_silent(upset_NA(airquality))
})

testthat::test_that("upset_NA error", {
  testthat::expect_error(upset_NA(2))
})
