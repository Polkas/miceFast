test_that("naive_fill_NA on empty", {
  testthat::expect_identical(naive_fill_NA(data.frame()), data.frame())
  testthat::expect_identical(naive_fill_NA(data.table::data.table()), data.table::data.table())
  testthat::expect_identical(naive_fill_NA(matrix()), matrix())
})

test_that("naive_fill_NA on one col", {
  data(air_miss)
  testthat::expect_false(anyNA(naive_fill_NA(air_miss[, 1, drop = FALSE])))
  testthat::expect_false(anyNA(naive_fill_NA(data.table::as.data.table(air_miss)[, 1])))
  testthat::expect_false(anyNA(naive_fill_NA(as.matrix(air_miss[, 1, drop = FALSE]))))
})

test_that("naive_fill_NA all NA filled", {
  data(air_miss)
  testthat::expect_true(sum(is.na(naive_fill_NA(air_miss))) == 0)
  testthat::expect_true(sum(is.na(naive_fill_NA(as.matrix(air_miss)))) == 0)
  testthat::expect_true(sum(is.na(naive_fill_NA(data.table::as.data.table(air_miss)))) == 0)
})

test_that("naive_fill_NA data.frame/data.table persist classes", {
  data(air_miss)
  cols_types <- vapply(air_miss, class, character(1))
  testthat::expect_identical(vapply(naive_fill_NA(air_miss), class, character(1)), cols_types)
  testthat::expect_identical(vapply(naive_fill_NA(data.table::as.data.table(air_miss)), class, character(1)), cols_types)
})

test_that("naive_fill_NA fill in loop", {
  testthat::expect_false(
    anyNA(do.call(rbind, Map(naive_fill_NA, split(air_miss, air_miss$groups))))
  )
})
