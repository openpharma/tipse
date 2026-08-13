library(testthat)
library(dplyr)

test_that("impute_landmark returns a list of data.frames with correct structure", {
  res <- impute_landmark(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    npts = 2, J = 3, seed = 123
  )

  expect_type(res, "list")
  expect_length(res, 3) # J = 3 imputations

  for (df in res) {
    expect_s3_class(df, "data.frame")
    expect_true(all(c("TRT01P", "AVAL", "EVENT", "CNSRRS", "AVALo", "EVENTo", "iter", "npts") %in% names(df)))
  }
})

test_that("impute_landmark imputes values correctly for treatment arm (EVENT=1 for imputed rows)", {
  res <- impute_landmark(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    npts = 2, J = 2, seed = 123
  )

  for (df in res) {
    imputed_rows <- df %>% filter(impute)
    expect_true(any(imputed_rows$EVENT == 1)) # EVENT should be 1 for imputed rows in trt arm
  }
})

test_that("impute_landmark imputes values correctly for control arm (AVAL extended to longest_control)", {
  res <- impute_landmark(codebreak200,
    reason = "Early dropout", impute = "docetaxel",
    npts = 2, J = 2, seed = 123
  )

  longest_control <- max(codebreak200$AVAL)
  for (df in res) {
    imputed_rows <- df %>% filter(impute)
    expect_true(any(imputed_rows$AVAL == longest_control)) # AVAL extended for control arm
  }
})

test_that("impute_landmark throws error for invalid impute arm", {
  expect_error(
    impute_landmark(codebreak200,
      reason = "Early dropout", impute = "invalid_arm",
      npts = 2, J = 2, seed = 123
    ),
    "Argument 'impute' must be one of the arms"
  )
})

test_that("impute_landmark warns when multiple arms passed for impute", {
  expect_error(
    impute_landmark(codebreak200,
      reason = "Early dropout", impute = c("docetaxel", "sotorasib"),
      npts = 2, J = 2, seed = 123
    ),
    "Imputation only in one arm is allowed"
  )
})

test_that("impute_landmark throws error for invalid npts", {
  expect_error(
    impute_landmark(codebreak200,
      reason = "Early dropout", impute = "sotorasib",
      npts = -1, J = 2, seed = 123
    ),
    "Argument 'npts' must be greater than 0"
  )
  expect_error(
    impute_landmark(codebreak200,
      reason = "Early dropout", impute = "sotorasib",
      npts = 100, J = 2, seed = 123
    ),
    "Argument 'npts' must be greater than 0"
  )
})

test_that("impute_landmark handles case with no matching reason (empty ids)", {
  expect_error(
    impute_landmark(codebreak200,
      reason = "Nonexistent reason", impute = "sotorasib",
      npts = 0, J = 2, seed = 123
    ),
    "No subject identified meeting imputation criteria."
  )
})

test_that("impute_landmark respects seed for reproducibility", {
  res1 <- impute_landmark(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    npts = 2, J = 2, seed = 123
  )
  res2 <- impute_landmark(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    npts = 2, J = 2, seed = 123
  )

  expect_equal(res1, res2)
})
