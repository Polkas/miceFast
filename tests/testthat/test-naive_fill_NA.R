test_that("get_replacements", {
  testthat::expect_identical(get_replacements(c(1,2,3)), numeric(0))
  testthat::expect_identical(get_replacements(1:3), integer(0))
  set.seed(1234)
  testthat::expect_identical(get_replacements(c(1:3, NA)), 2L)
  set.seed(1234)
  testthat::expect_identical(get_replacements(c(1:3, NA), c(TRUE, TRUE, TRUE, FALSE)), 2L)
  set.seed(1234)
  testthat::expect_identical(get_replacements(c(1:3, NA), c(TRUE, TRUE, TRUE, FALSE), 3), 2L)
  set.seed(1234)
  testthat::expect_identical(get_replacements(c(1:3, NA), c(TRUE, TRUE, TRUE, FALSE), 3, 1), 2L)
})

test_that("naive_fill_NA on empty", {
  testthat::expect_identical(naive_fill_NA(data.frame()), data.frame())
  testthat::expect_identical(naive_fill_NA(data.table::data.table()), data.table::data.table())
  testthat::expect_identical(naive_fill_NA(matrix()), matrix())
})

test_that("naive_fill_NA with one non missing", {
  testthat::expect_identical(naive_fill_NA(data.frame(a = c(1, NA))), data.frame(a = c(1, 1)))
  testthat::expect_identical(naive_fill_NA(data.table::data.table(a = c(1, NA))), data.table::data.table(a = c(1, 1)))
  testthat::expect_identical(naive_fill_NA(matrix(c(1, NA), ncol = 1)), matrix(c(1, 1), ncol = 1))

  testthat::expect_identical(naive_fill_NA(data.frame(a = c(1, NA, NA))), data.frame(a = c(1, 1, 1)))
  testthat::expect_identical(naive_fill_NA(data.table::data.table(a = c(1, NA, NA))), data.table::data.table(a = c(1, 1, 1)))
  testthat::expect_identical(naive_fill_NA(matrix(c(1, NA, NA), ncol = 1)), matrix(c(1, 1, 1), ncol = 1))
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
