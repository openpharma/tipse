library(testthat)
library(dplyr)

test_that("average_km throws error for invalid km_data input", {
  expect_error(
    average_km(NULL, arm = "A"),
    "'km_data' must be a non-empty list"
  )

  expect_error(
    average_km(list(), arm = "A"),
    "'km_data' must be a non-empty list"
  )
})

test_that("average_km throws error when required columns are missing", {
  bad_df <- data.frame(AVAL = 1:3, EVENT = c(1, 0, 1))
  expect_error(
    average_km(list(bad_df), arm = "A"),
    "Each km_data element must contain columns: AVAL, EVENT, TRT01P, iter"
  )
})

test_that("average_km returns expected structure", {
  # minimal reproducible KM dataset for two imputations
  df1 <- data.frame(
    AVAL = c(2, 4, 6),
    EVENT = c(1, 1, 1),
    TRT01P = "A",
    iter = 1,
    tipping_param = 10
  )

  df2 <- data.frame(
    AVAL = c(2, 4, 6), # same time points as df1
    EVENT = c(1, 1, 0),
    TRT01P = "A",
    iter = 2,
    tipping_param = 10
  )

  res <- average_km(km_data = list(df1, df2), arm = "A")

  expect_s3_class(res, "data.frame")
  expect_true(all(c(
    "time", "parameter", "survival_comb",
    "survival_lcl_comb", "survival_ucl_comb", "stderr"
  ) %in% names(res)))
})

test_that("average_km correctly filters by arm and iter", {
  df1 <- data.frame(
    AVAL = c(5, 8, 10),
    EVENT = c(1, 0, 0),
    TRT01P = c("A", "A", "B"), # arm B appears but should not be pooled
    iter = 1,
    tipping = 20
  )

  df2 <- data.frame(
    AVAL = c(5, 8),
    EVENT = c(1, 0),
    TRT01P = c("A", "A"),
    iter = 2,
    tipping = 20
  )

  res <- average_km(km_data = list(df1, df2), arm = "A")

  expect_true(all(res$parameter == 20))
})

test_that("average_km produces valid survival probabilities", {
  df1 <- data.frame(
    AVAL = c(2, 4, 6),
    EVENT = c(1, 1, 0),
    TRT01P = c("A", "A", "A"),
    iter = 1,
    tipping = 30
  )

  df2 <- data.frame(
    AVAL = c(2, 4, 6),
    EVENT = c(1, 0, 1),
    TRT01P = c("A", "A", "A"),
    iter = 2,
    tipping = 30
  )

  res <- average_km(km_data = list(df1, df2), arm = "A")

  expect_true(all(res$survival_comb >= 0 & res$survival_comb <= 1))
  expect_true(all(res$survival_lcl_comb >= 0))
  expect_true(all(res$survival_ucl_comb <= 1))
})

test_that("average_km throws an error when required columns are missing", {
  bad_df <- data.frame(
    AVAL = c(5, 6, 7),
    EVENT = c(1, 0, 1),
    # TRT01P missing
    iter = c(1, 1, 1),
    tipping = 50
  )

  expect_error(
    average_km(km_data = list(bad_df), arm = "A"),
    "must contain columns"
  )
})


test_that("average_km back transformation gives valid survival probabilities", {
  df1 <- data.frame(
    AVAL = c(2, 4, 6),
    EVENT = c(1, 1, 0),
    TRT01P = c("A", "A", "A"),
    iter = 1,
    tipping = 30
  )

  df2 <- data.frame(
    AVAL = c(2, 4, 6),
    EVENT = c(1, 0, 1),
    TRT01P = c("A", "A", "A"),
    iter = 2,
    tipping = 30
  )

  res <- average_km(km_data = list(df1, df2), arm = "A")

  expect_true(all(res$survival_comb > 0 & res$survival_comb < 1))
  expect_true(all(res$survival_lcl_comb > 0 & res$survival_lcl_comb < 1))
  expect_true(all(res$survival_ucl_comb > 0 & res$survival_ucl_comb < 1))

  # Lower ≤ estimate ≤ upper
  expect_true(all(res$survival_lcl_comb <= res$survival_comb))
  expect_true(all(res$survival_comb <= res$survival_ucl_comb))
})
