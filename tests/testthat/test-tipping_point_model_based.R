library(testthat)
library(survival)
library(dplyr)

# Cox model for testing
cox_fit <- coxph(Surv(AVAL, EVENT) ~ TRT01P, ties = "exact", data = codebreak200)

test_that("tipping_point_model_based returns correct structure for model-based imputations", {
  res <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "docetaxel",
    imputation_model = "weibull",
    J = 2,
    tipping_range = seq(0.1, 1, by = 0.1),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )

  expect_s3_class(res, "tipse")
  expect_true(all(c(
    "original_data", "original_HR", "reason_to_impute", "arm_to_impute",
    "method_to_impute", "imputation_results", "imputation_data"
  ) %in% names(res)))

  # Check imputation_results structure
  expect_true(all(c("HR", "HR_upperCI", "HR_lowerCI", "parameter_docetaxel", "tipping_point") %in% names(res$imputation_results)))

  # imputation_data should be a named list
  expect_true(is.list(res$imputation_data))
  expect_equal(length(res$imputation_data), length(seq(0.1, 1, by = 0.1)))
})

test_that("tipping_point_model_based throws error for invalid cox_fit", {
  expect_error(
    tipping_point_model_based(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "docetaxel",
      imputation_model = "weibull",
      J = 2,
      tipping_range = seq(0.1, 1, by = 0.1),
      cox_fit = "not_a_coxph"
    ),
    "Argument 'cox_fit' must be a valid cox model object"
  )
})

test_that("tipping_point_model_based throws error for invalid impute arm", {
  expect_error(
    tipping_point_model_based(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "invalid_arm",
      imputation_model = "weibull",
      J = 2,
      tipping_range = seq(0.1, 1, by = 0.1),
      cox_fit = cox_fit
    ),
    "Argument 'impute' must be one or both of the arms"
  )
})

test_that("tipping_point_model_based throws error when reason is empty", {
  expect_error(
    tipping_point_model_based(
      dat = codebreak200,
      reason = character(0),
      impute = "sotorasib",
      imputation_model = "weibull",
      J = 2,
      tipping_range = seq(0.1, 1, by = 0.1),
      cox_fit = cox_fit
    ),
    "Argument 'reason' must specify at least one censoring reason"
  )
})

test_that("tipping_point_model_based validates distribution argument", {
  expect_error(
    tipping_point_model_based(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "sotorasib",
      imputation_model = "invalid_dist",
      J = 2,
      tipping_range = seq(0.1, 1, by = 0.1),
      cox_fit = cox_fit
    ),
    "'arg' should be one of \"weibull\", \"exponential\""
  )
})

test_that("tipping_point_model_based sets tipping_point flag correctly", {
  res <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "docetaxel",
    imputation_model = "weibull",
    J = 2,
    tipping_range = seq(0.1, 1, by = 0.1),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )

  expect_true("tipping_point" %in% names(res$imputation_results))
  expect_true(any(res$imputation_results$tipping_point %in% c(TRUE, FALSE)))
})

test_that("tipping_point_model_based respects seed for reproducibility", {
  res1 <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "sotorasib",
    imputation_model = "weibull",
    J = 2,
    tipping_range = seq(0.1, 1, by = 0.1),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )

  res2 <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "sotorasib",
    imputation_model = "weibull",
    J = 2,
    tipping_range = seq(0.1, 1, by = 0.1),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )

  expect_equal(res1$imputation_results, res2$imputation_results)
})

test_that("tipping_point_model_based detects method correctly based on tipping_range", {
  res_deflation <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "docetaxel",
    imputation_model = "weibull",
    J = 2,
    tipping_range = seq(0.05, 1, by = 0.05),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )
  expect_equal(res_deflation$method_to_impute, "hazard multiplication")

  res_inflation <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "sotorasib",
    imputation_model = "weibull",
    J = 2,
    tipping_range = c(1.1, 1.2, 1.5),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )
  expect_equal(res_inflation$method_to_impute, "hazard multiplication")
})

test_that("tipping_point_model_based two-arm: returns correct structure", {
  res <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = c("docetaxel", "sotorasib"),
    imputation_model = "weibull",
    J = 2,
    tipping_range = c(0.7, 0.9),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )

  expect_s3_class(res, "tipse")
  expect_equal(res$arm_to_impute, c("docetaxel", "sotorasib"))
  expect_equal(res$method_to_impute, "hazard multiplication")

  # Two-arm results should have arm-named parameter columns, not generic ones
  expect_true("parameter_docetaxel" %in% names(res$imputation_results))
  expect_true("parameter_sotorasib" %in% names(res$imputation_results))
  expect_false("parameter" %in% names(res$imputation_results))
  expect_false("parameter_1" %in% names(res$imputation_results))

  # Grid: 2 x 2 = 4 rows
  expect_equal(nrow(res$imputation_results), 4L)
  expect_equal(length(res$imputation_data), 4L)
})

test_that("tipping_point_model_based two-arm: per-arm list tipping_range produces Cartesian grid", {
  res <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = c("docetaxel", "sotorasib"),
    imputation_model = "weibull",
    J = 2,
    tipping_range = list(c(0.6, 0.7, 0.8), c(0.8, 0.9)),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )

  # 3 x 2 = 6 rows
  expect_equal(nrow(res$imputation_results), 6L)
  expect_equal(length(unique(res$imputation_results$parameter_docetaxel)), 3L)
  expect_equal(length(unique(res$imputation_results$parameter_sotorasib)), 2L)
})

test_that("tipping_point_model_based two-arm: mixed direction (deflation + inflation) runs without error", {
  res <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = c("docetaxel", "sotorasib"),
    imputation_model = "weibull",
    J = 2,
    tipping_range = list(c(0.7, 0.9), c(1.1, 1.3)),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )

  expect_s3_class(res, "tipse")
  expect_equal(nrow(res$imputation_results), 4L)
  expect_true(all(res$imputation_results$parameter_docetaxel <= 1))
  expect_true(all(res$imputation_results$parameter_sotorasib >= 1))
})

test_that("tipping_point_model_based two-arm: km_data_list names use pipe-separated arm parameters", {
  res <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = c("docetaxel", "sotorasib"),
    imputation_model = "weibull",
    J = 2,
    tipping_range = c(0.8, 0.9),
    cox_fit = cox_fit,
    verbose = FALSE,
    seed = 123
  )

  expect_true(all(grepl("|", names(res$imputation_data), fixed = TRUE)))
})

test_that("tipping_point_model_based two-arm: errors on duplicate arms", {
  expect_error(
    tipping_point_model_based(
      dat = codebreak200,
      reason = "Early dropout",
      impute = c("docetaxel", "docetaxel"),
      imputation_model = "weibull",
      J = 2,
      tipping_range = c(0.7, 0.9),
      cox_fit = cox_fit
    ),
    "must not contain the same arm twice"
  )
})

test_that("tipping_point_model_based two-arm: errors when list tipping_range has wrong length", {
  expect_error(
    tipping_point_model_based(
      dat = codebreak200,
      reason = "Early dropout",
      impute = c("docetaxel", "sotorasib"),
      imputation_model = "weibull",
      J = 2,
      tipping_range = list(c(0.7, 0.9)),
      cox_fit = cox_fit
    ),
    "exactly two elements"
  )
})

test_that("tipping_point_model_based two-arm: list tipping_range rejected for single arm", {
  expect_error(
    tipping_point_model_based(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "docetaxel",
      imputation_model = "weibull",
      J = 2,
      tipping_range = list(c(0.7, 0.9), c(0.8, 1.0)),
      cox_fit = cox_fit
    ),
    "only allowed when two arms"
  )
})

test_that("tipping_point_model_based two-arm: reproducible with seed", {
  args <- list(
    dat = codebreak200, reason = "Early dropout",
    impute = c("docetaxel", "sotorasib"), imputation_model = "weibull",
    J = 2, tipping_range = c(0.7, 0.9), cox_fit = cox_fit,
    verbose = FALSE, seed = 42
  )
  res1 <- do.call(tipping_point_model_based, args)
  res2 <- do.call(tipping_point_model_based, args)
  expect_equal(res1$imputation_results, res2$imputation_results)
})
