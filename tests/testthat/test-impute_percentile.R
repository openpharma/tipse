library(testthat)
library(dplyr)
library(purrr)

test_that("impute_percentile returns a list of data.frames with correct structure", {
  res <- impute_percentile(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    percentile = 20, J = 3, seed = 123
  )

  expect_type(res, "list")
  expect_length(res, 3) # J = 3 imputations

  for (df in res) {
    expect_s3_class(df, "data.frame")
    expect_true(all(c("TRT01P", "AVAL", "EVENT", "CNSRRS", "AVALo", "EVENTo", "iter", "percentile") %in% names(df)))
  }
})

test_that("impute_percentile imputes values correctly (AVAL changes for imputed rows)", {
  res <- impute_percentile(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    percentile = 50, J = 2, seed = 123
  )

  for (df in res) {
    expect_true(any(df$AVAL[df$impute] != df$AVALo[df$impute]))
  }
})

test_that("impute_percentile throws error for invalid impute arm", {
  expect_error(
    impute_percentile(codebreak200,
      reason = "Early dropout", impute = "invalid_arm",
      percentile = 20, J = 2
    ),
    "Argument 'impute' must be one of the arms"
  )
})

test_that("impute_percentile warns when multiple arms passed for impute", {
  expect_error(
    impute_percentile(codebreak200,
      reason = "Early dropout", impute = c("control", "sotorasib"),
      percentile = 20, J = 2
    ),
    "Imputation only in one arm is allowed"
  )
})

test_that("impute_percentile throws error for invalid percentile", {
  expect_error(
    impute_percentile(codebreak200,
      reason = "Early dropout", impute = "sotorasib",
      percentile = 0, J = 2
    ),
    "Argument 'percentile' must be greater than 0 and less or equal to 100"
  )
  expect_error(
    impute_percentile(codebreak200,
      reason = "Early dropout", impute = "sotorasib",
      percentile = 150, J = 2
    ),
    "Argument 'percentile' must be greater than 0 and less or equal to 100"
  )
})

test_that("impute_percentile handles case with no matching reason (empty ids)", {
  expect_error(
    impute_percentile(codebreak200,
      reason = "Nonexistent reason", impute = "sotorasib",
      percentile = 20, J = 2
    ),
    "No subject identified meeting imputation criteria."
  )
})

test_that("impute_percentile respects seed for reproducibility", {
  res1 <- impute_percentile(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    percentile = 20, J = 2, seed = 123
  )
  res2 <- impute_percentile(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    percentile = 20, J = 2, seed = 123
  )

  expect_equal(res1, res2)
})

test_that("impute_percentile computes EVENT correctly after imputation", {
  res <- impute_percentile(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    percentile = 20, J = 2, seed = 123
  )

  for (df in res) {
    expect_true(all(df$EVENT %in% c(0, 1)))
  }
})
