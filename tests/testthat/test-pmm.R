# ===========================================================================
# PMM fix: imputed values must be actual observed y values
# ===========================================================================

testthat::test_that("PMM returns observed y values (numeric, fill_NA_N)", {
    set.seed(42)
    data(air_miss)

    result <- air_miss %>%
        mutate(
            Ozone_pmm = fill_NA_N(
                x = .,
                model = "pmm",
                posit_y = "Ozone",
                posit_x = c("Solar.R", "Wind", "Temp"),
                k = 5
            )
        )

    observed_vals <- na.omit(air_miss$Ozone)
    na_idx <- which(is.na(air_miss$Ozone))
    # rows where all predictors are also observed (can actually be imputed)
    can_impute <- na_idx[complete.cases(air_miss[
        na_idx,
        c("Solar.R", "Wind", "Temp")
    ])]

    imputed_vals <- result$Ozone_pmm[can_impute]
    testthat::expect_true(
        all(imputed_vals %in% observed_vals),
        info = "PMM must return actual observed y values, not predicted values"
    )
})

testthat::test_that("PMM returns observed y values (numeric, OOP interface)", {
    set.seed(42)
    data(air_miss)

    dat <- as.matrix(air_miss[, c("Ozone", "Wind", "Temp")])
    obj <- new(miceFast)
    obj$set_data(dat)

    res <- obj$impute_N("pmm", 1, c(2, 3), 5)
    imputed_vals <- res$imputations

    observed_vals <- na.omit(dat[, 1])
    testthat::expect_true(
        all(imputed_vals %in% observed_vals),
        info = "OOP PMM must return actual observed y values"
    )
})

testthat::test_that("PMM returns observed y values (numeric, data.table)", {
    set.seed(42)
    data(air_miss)
    setDT(air_miss)

    air_miss[,
        Ozone_pmm := fill_NA_N(
            x = .SD,
            model = "pmm",
            posit_y = "Ozone",
            posit_x = c("Solar.R", "Wind", "Temp"),
            k = 5
        )
    ]

    observed_vals <- na.omit(air_miss$Ozone)
    na_idx <- which(is.na(air_miss$Ozone))
    can_impute <- na_idx[complete.cases(air_miss[
        na_idx,
        c("Solar.R", "Wind", "Temp")
    ])]

    imputed_vals <- air_miss$Ozone_pmm[can_impute]
    testthat::expect_true(all(imputed_vals %in% observed_vals))
})

testthat::test_that("PMM with k=1 returns nearest observed value", {
    set.seed(42)
    data(air_miss)

    result <- air_miss %>%
        mutate(
            Ozone_pmm = fill_NA_N(
                x = .,
                model = "pmm",
                posit_y = "Ozone",
                posit_x = c("Solar.R", "Wind", "Temp"),
                k = 1
            )
        )

    observed_vals <- na.omit(air_miss$Ozone)
    na_idx <- which(is.na(air_miss$Ozone))
    can_impute <- na_idx[complete.cases(air_miss[
        na_idx,
        c("Solar.R", "Wind", "Temp")
    ])]

    imputed_vals <- result$Ozone_pmm[can_impute]
    testthat::expect_true(all(imputed_vals %in% observed_vals))
})

testthat::test_that("PMM respects observed range", {
    set.seed(42)
    data(air_miss)

    result <- air_miss %>%
        mutate(
            Ozone_pmm = fill_NA_N(
                x = .,
                model = "pmm",
                posit_y = "Ozone",
                posit_x = c("Solar.R", "Wind", "Temp"),
                k = 10
            )
        )

    obs_range <- range(air_miss$Ozone, na.rm = TRUE)
    na_idx <- which(is.na(air_miss$Ozone))
    can_impute <- na_idx[complete.cases(air_miss[
        na_idx,
        c("Solar.R", "Wind", "Temp")
    ])]

    imputed_vals <- result$Ozone_pmm[can_impute]
    testthat::expect_true(all(imputed_vals >= obs_range[1]))
    testthat::expect_true(all(imputed_vals <= obs_range[2]))
})

testthat::test_that("PMM weighted returns observed y values", {
    set.seed(42)
    data(air_miss)

    result <- air_miss %>%
        mutate(
            Ozone_pmm = fill_NA_N(
                x = .,
                model = "pmm",
                posit_y = "Ozone",
                posit_x = c("Wind", "Temp"),
                w = weights,
                k = 5
            )
        )

    observed_vals <- na.omit(air_miss$Ozone)
    na_idx <- which(is.na(air_miss$Ozone))
    can_impute <- na_idx[complete.cases(air_miss[na_idx, c("Wind", "Temp")])]

    imputed_vals <- result$Ozone_pmm[can_impute]
    testthat::expect_true(
        all(imputed_vals %in% observed_vals),
        info = "Weighted PMM must also return observed y values"
    )
})

# ===========================================================================
# PMM with categorical variables (factor and character)
# ===========================================================================

testthat::test_that("PMM works with factor y variable", {
    set.seed(42)
    data(air_miss)

    # Ozone_chac is character but let's test a proper factor
    air_test <- air_miss
    air_test$oz_factor <- factor(ifelse(air_miss$Ozone > 30, "high", "low"))
    # Inject some NAs
    air_test$oz_factor[sample(which(!is.na(air_test$oz_factor)), 20)] <- NA

    result <- air_test %>%
        mutate(
            oz_imp = fill_NA_N(
                x = .,
                model = "pmm",
                posit_y = "oz_factor",
                posit_x = c("Wind", "Temp"),
                k = 5
            )
        )

    valid_levels <- levels(air_test$oz_factor)
    # All imputed values should be valid levels
    na_idx <- which(is.na(air_test$oz_factor))
    testthat::expect_true(all(result$oz_imp[na_idx] %in% valid_levels))
    testthat::expect_equal(sum(is.na(result$oz_imp)), 0)
})

testthat::test_that("PMM works with character y variable (non-numeric labels)", {
    set.seed(42)
    data(air_miss)

    # x_character has non-numeric labels like "(0,70]", "(70,140]"
    testthat::expect_true(is.character(air_miss$x_character))
    n_na_before <- sum(is.na(air_miss$x_character))
    testthat::expect_true(n_na_before > 0)

    result <- air_miss %>%
        mutate(
            x_char_imp = fill_NA_N(
                x = .,
                model = "pmm",
                posit_y = "x_character",
                posit_x = c("Wind", "Temp"),
                k = 5
            )
        )

    observed_levels <- unique(na.omit(air_miss$x_character))

    testthat::expect_true(is.character(result$x_char_imp))
    # All imputed values should be valid observed levels
    na_idx <- which(is.na(air_miss$x_character))
    can_impute <- na_idx[complete.cases(air_miss[na_idx, c("Wind", "Temp")])]
    testthat::expect_true(
        all(result$x_char_imp[can_impute] %in% observed_levels),
        info = "PMM with character y must return valid observed categories"
    )
    testthat::expect_equal(sum(is.na(result$x_char_imp)), 0)
})

testthat::test_that("PMM with character y on data.table", {
    set.seed(42)
    data(air_miss)
    setDT(air_miss)

    n_na_before <- sum(is.na(air_miss$x_character))

    air_miss[,
        x_char_imp := fill_NA_N(
            x = .SD,
            model = "pmm",
            posit_y = "x_character",
            posit_x = c("Wind", "Temp"),
            k = 5
        )
    ]

    observed_levels <- unique(na.omit(air_miss$x_character))
    testthat::expect_true(is.character(air_miss$x_char_imp))
    testthat::expect_equal(sum(is.na(air_miss$x_char_imp)), 0)
    testthat::expect_true(all(air_miss$x_char_imp %in% observed_levels))
})

# ===========================================================================
# PMM reproducibility
# ===========================================================================

testthat::test_that("PMM is reproducible with set.seed", {
    data(air_miss)

    set.seed(999)
    r1 <- air_miss %>%
        mutate(
            Ozone_pmm = fill_NA_N(
                x = .,
                model = "pmm",
                posit_y = "Ozone",
                posit_x = c("Solar.R", "Wind", "Temp"),
                k = 5
            )
        )

    set.seed(999)
    r2 <- air_miss %>%
        mutate(
            Ozone_pmm = fill_NA_N(
                x = .,
                model = "pmm",
                posit_y = "Ozone",
                posit_x = c("Solar.R", "Wind", "Temp"),
                k = 5
            )
        )

    testthat::expect_identical(r1$Ozone_pmm, r2$Ozone_pmm)
})

# ===========================================================================
# PMM with grouped imputation
# ===========================================================================

testthat::test_that("PMM with grouped data.table returns observed values", {
    set.seed(42)
    data(air_miss)
    setDT(air_miss)

    air_miss[,
        Ozone_pmm := fill_NA_N(
            x = .SD,
            model = "pmm",
            posit_y = "Ozone",
            posit_x = c("Wind", "Temp", "Intercept"),
            k = 5
        ),
        by = .(groups)
    ]

    observed_vals <- na.omit(air_miss$Ozone)
    na_idx <- which(is.na(air_miss$Ozone))
    can_impute <- na_idx[complete.cases(air_miss[na_idx, c("Wind", "Temp")])]

    imputed_vals <- air_miss$Ozone_pmm[can_impute]
    testthat::expect_true(all(imputed_vals %in% observed_vals))
})

# ===========================================================================
# PMM on matrix input
# ===========================================================================

testthat::test_that("PMM returns observed values with matrix input", {
    set.seed(42)
    data(air_miss)

    mat <- as.matrix(air_miss[, c("Ozone", "Solar.R", "Wind", "Temp")])
    result <- fill_NA_N(
        mat,
        model = "pmm",
        posit_y = 1,
        posit_x = c(3, 4),
        k = 5
    )

    observed_vals <- na.omit(mat[, 1])
    na_idx <- which(is.na(mat[, 1]))
    can_impute <- na_idx[complete.cases(mat[na_idx, c(3, 4)])]

    testthat::expect_true(all(result[can_impute] %in% observed_vals))
})
