
context("miceFast-auto_imputes_funcs")


test_that("auto_imputes", {
  data(air_miss)

  air_miss_mat <- model.matrix.lm(~., data = air_miss)

  test <- sum(is.na(naive_fill_NA(air_miss))) == 0

  test_1 <- identical(vapply(naive_fill_NA(air_miss), class, character(1)), vapply(air_miss, class, character(1)))

  test_2 <- sum(is.na(naive_fill_NA(air_miss_mat))) == 0

  test_3 <- sum(is.na(naive_fill_NA(as.data.frame(air_miss)))) == 0


  test_all <- expect_true(all(c(test, test_1, test_2, test_3)))
})
