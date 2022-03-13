test_that("naive_fill_NA matrix all NA filled", {
  data(air_miss)

  testthat::expect_true(sum(is.na(naive_fill_NA(as.matrix(airquality)))) == 0)
})

test_that("naive_fill_NA model.matrix.lm", {
  data(air_miss)

  air_miss_mat <- model.matrix.lm(~., data = air_miss)

  testthat::expect_true(sum(is.na(naive_fill_NA(air_miss_mat))) == 0)
})

test_that("naive_fill_NA matrix persist class", {
  data(air_miss)

  testthat::expect_true(is.numeric(as.matrix(airquality)))
})

test_that("naive_fill_NA data.frame all NA filled", {
  data(air_miss)

  testthat::expect_true(sum(is.na(naive_fill_NA(air_miss))) == 0)
})

test_that("naive_fill_NA data.frame persist classes", {
  data(air_miss)

  testthat::expect_identical(vapply(naive_fill_NA(air_miss), class, character(1)), vapply(air_miss, class, character(1)))
})
