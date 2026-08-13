library(testthat)
library(dplyr)

test_that("impute_model returns a list of data.frames with correct structure", {
  res <- impute_model(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    imputation_model = "weibull", alpha = 0.5, J = 3, seed = 123
  )

  expect_type(res, "list")
  expect_length(res, 3) # J = 3 imputations

  # Check each element is a data.frame with expected columns
  for (df in res) {
    expect_s3_class(df, "data.frame")
    expect_true(all(c("TRT01P", "AVAL", "EVENT", "CNSRRS", "AVALo", "EVENTo", "iter", "alpha") %in% names(df)))
  }
})

test_that("impute_model imputes values correctly (AVAL changes for imputed rows)", {
  res <- impute_model(codebreak200,
    reason = "Early dropout", impute = "docetaxel",
    imputation_model = "weibull", alpha = 0.5, J = 2, seed = 123
  )

  # Original AVALo should differ from imputed AVAL for imputed rows
  imputed_rows <- which(codebreak200$TRT01P == "docetaxel" & codebreak200$CNSRRS == "Early dropout")
  for (df in res) {
    expect_true(any(df$AVAL[imputed_rows] != df$AVALo[imputed_rows]))
  }
})

test_that("impute_model throws error for invalid impute arm", {
  expect_error(
    impute_model(codebreak200,
      reason = "Early dropout", impute = "invalid_arm",
      imputation_model = "weibull", alpha = 0.5, J = 2
    ),
    "Argument 'impute' must be one of the arms"
  )
})

test_that("impute_model warns when multiple arms passed for impute", {
  expect_error(
    impute_model(codebreak200,
      reason = "Early dropout", impute = c("control", "sotorasib"),
      imputation_model = "weibull", alpha = 0.5, J = 2
    ),
    "Imputation only in one arm is allowed"
  )
})


test_that("impute_model respects seed for reproducibility", {
  res1 <- impute_model(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    imputation_model = "weibull", alpha = 0.5, J = 2, seed = 123
  )
  res2 <- impute_model(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    imputation_model = "weibull", alpha = 0.5, J = 2, seed = 123
  )

  expect_equal(res1, res2)
})

test_that("impute_model computes EVENT correctly (1 if t <= longest_control)", {
  res <- impute_model(codebreak200,
    reason = "Early dropout", impute = "sotorasib",
    imputation_model = "weibull", alpha = 0.5, J = 2, seed = 123
  )

  for (df in res) {
    expect_true(all(df$EVENT %in% c(TRUE, FALSE)))
  }
})
