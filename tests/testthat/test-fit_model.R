library(testthat)
library(dplyr)
library(survival)

test_that("fit_model returns expected columns", {
  res <- fit_model(codebreak200, reason = "Early dropout", impute = "sotorasib", imputation_model = "weibull")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("TRT01P", "AVAL", "EVENT", "CNSRRS", "AVALo", "EVENTo", "impute", "a", "b", "cdf") %in% names(res)))
})

test_that("fit_model flags imputed subjects correctly", {
  res <- fit_model(codebreak200, reason = "Early dropout", impute = "sotorasib", imputation_model = "weibull")
  flagged_ids <- which(codebreak200$TRT01P == "sotorasib" & codebreak200$CNSRRS == "Early dropout")
  expect_true(all(res$impute[flagged_ids]))
  expect_false(any(res$impute[!flagged_ids]))
})

test_that("fit_model computes correct cdf for imputed rows", {
  res <- fit_model(codebreak200, reason = "Early dropout", impute = "sotorasib", imputation_model = "weibull")
  imputed_rows <- which(res$impute)
  for (i in imputed_rows) {
    expected_cdf <- 1 - exp(-(res$AVAL[i] / res$b[i])^res$a[i])
    expect_equal(res$cdf[i], expected_cdf)
  }
})

test_that("fit_model throws error for invalid impute arm", {
  expect_error(
    fit_model(codebreak200, reason = "Early dropout", impute = "invalid_arm", imputation_model = "weibull"),
    "Argument 'impute' must be one of the arms"
  )
})

test_that("fit_model warns when multiple arms passed for impute", {
  expect_error(
    fit_model(codebreak200, reason = "Early dropout", impute = c("control", "sotorasib"), imputation_model = "weibull"),
    "Imputation only in one arm is allowed"
  )
})

test_that("fit_model works with exponential distribution", {
  res <- fit_model(codebreak200, reason = "Early dropout", impute = "sotorasib", imputation_model = "exponential")
  expect_true(all(res$a == 1)) # For exponential, shape parameter should be 1
  expect_true(all(res$b > 0))
})

test_that("fit_model handles case with no matching reason (empty ids)", {
  expect_error(
    fit_model(codebreak200, reason = "Nonexistent reason", impute = "sotorasib", imputation_model = "weibull"),
    "No subject identified meeting imputation criteria."
  )
})

test_that("fit_model computes positive parameters", {
  res <- fit_model(codebreak200, reason = "Early dropout", impute = "sotorasib", imputation_model = "weibull")
  expect_true(all(res$a > 0))
  expect_true(all(res$b > 0))
})
